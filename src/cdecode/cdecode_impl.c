/*
 * cdecode_impl.c
 *
 * High-performance numeric-to-string decoding for Stata
 * Part of the ctools suite
 *
 * Algorithm:
 * 1. Parse value-label mappings passed from Stata
 * 2. Build hash table mapping integer values to string labels
 * 3. Bulk load numeric source variable via ctools_data_load()
 * 4. For each observation:
 *    a. Look up corresponding label in hash table (from loaded data)
 *    b. Write label string to destination variable (must be sequential)
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_timer.h"
#include "ctools_config.h"
#include "ctools_arena.h"
#include "ctools_hash.h"
#include "cdecode_impl.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CDECODE_HASH_INIT_SIZE  256
#define CDECODE_HASH_LOAD_FACTOR 0.75
#define CDECODE_MAX_LABELS      65536
#define CDECODE_MAX_LABEL_LEN   2045

/* ============================================================================
 * Hash Table for Integer -> String Mapping
 * Uses unified ctools_int_hash_table from ctools_hash.h
 * ============================================================================ */

/* Compatibility typedefs and macros for cdecode */
typedef ctools_int_hash_table cdecode_hash_table;

#define cdecode_hash_init   ctools_int_hash_init
#define cdecode_hash_free   ctools_int_hash_free
#define cdecode_hash_insert ctools_int_hash_insert
#define cdecode_hash_lookup ctools_int_hash_lookup

/* ============================================================================
 * Label Parsing
 * ============================================================================ */

/*
 * Unescape a label string.
 * Converts \| to |, \\ to \, and \" to "
 */
static void unescape_label(const char *label, char *unescaped, size_t max_len)
{
    size_t j = 0;
    for (size_t i = 0; label[i] && j < max_len; i++) {
        if (label[i] == '\\' && label[i+1]) {
            i++;
            if (label[i] == '|') {
                unescaped[j++] = '|';
            } else if (label[i] == '\\') {
                unescaped[j++] = '\\';
            } else if (label[i] == '"') {
                unescaped[j++] = '"';
            } else {
                unescaped[j++] = label[i];
            }
        } else {
            unescaped[j++] = label[i];
        }
    }
    unescaped[j] = '\0';
}

/*
 * Parse value-label pairs from a file.
 * File format: one "value|label" per line (labels are escaped)
 * Returns 0 on success, -1 on error.
 */
static int parse_labels_from_file(const char *filepath, cdecode_hash_table *ht, int *max_len)
{
    if (!filepath || *filepath == '\0') {
        return 0;  /* No file is OK */
    }

    FILE *fp = fopen(filepath, "r");
    if (!fp) {
        return -1;
    }

    *max_len = 0;
    char line[CDECODE_MAX_LABEL_LEN + 32];

    while (fgets(line, sizeof(line), fp)) {
        /* Remove trailing newline */
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) {
            line[--len] = '\0';
        }

        if (len == 0) continue;

        /* Find the pipe separator between value and label */
        char *pipe = strchr(line, '|');
        if (!pipe) continue;

        *pipe = '\0';
        int value;
        if (!ctools_safe_atoi(line, &value)) {
            fclose(fp);
            return -1;  /* Invalid value in labels file */
        }
        char *label = pipe + 1;

        /* Unescape the label */
        char unescaped[CDECODE_MAX_LABEL_LEN + 1];
        unescape_label(label, unescaped, CDECODE_MAX_LABEL_LEN);

        /* Track max length */
        int lbl_len = (int)strlen(unescaped);
        if (lbl_len > *max_len) {
            *max_len = lbl_len;
        }

        /* Insert into hash table */
        if (cdecode_hash_insert(ht, value, unescaped) != 0) {
            fclose(fp);
            return -1;
        }
    }

    fclose(fp);
    return 0;
}

/* ============================================================================
 * Main Implementation
 * ============================================================================ */

ST_retcode cdecode_main(const char *args)
{
    double t_start, t_parse, t_decode, t_total;
    int maxlen = 0;
    char labelsfile[1024] = "";
    size_t nobs;
    ST_retcode rc = 0;

    /*
     * IMPORTANT: In Stata's Plugin Interface, SF_var_is_string() and other
     * SF_* functions access variables by their 1-based index within the
     * varlist passed to the plugin call, NOT by their position in the full
     * dataset.
     *
     * cdecode is called as: plugin call ctools_plugin src_var dst_var, args
     * So the source variable is always index 1, and destination is index 2.
     */
    const int src_idx = 1;  /* First variable passed to plugin = source */
    const int dst_idx = 2;  /* Second variable passed to plugin = destination */

    t_start = ctools_timer_seconds();

    /* Parse arguments */
    if (args == NULL || *args == '\0') {
        SF_error("cdecode: no arguments specified\n");
        return 198;
    }

    /* Make a copy for parsing */
    size_t args_len = strlen(args);
    char *args_copy = malloc(args_len + 1);
    if (!args_copy) {
        SF_error("cdecode: memory allocation failed\n");
        return 920;
    }
    memcpy(args_copy, args, args_len + 1);  /* Include null terminator */

    /* Find labelsfile= parameter */
    char *labelsfile_ptr = strstr(args_copy, "labelsfile=");
    if (labelsfile_ptr) {
        char *filepath_start = labelsfile_ptr + 11;
        /* Copy until end of string or whitespace */
        size_t i = 0;
        while (filepath_start[i] && filepath_start[i] != ' ' && filepath_start[i] != '\t' && i < sizeof(labelsfile) - 1) {
            labelsfile[i] = filepath_start[i];
            i++;
        }
        labelsfile[i] = '\0';
        *labelsfile_ptr = '\0'; /* Terminate args before labelsfile= */
    }

    /* Parse other arguments */
    char *token = strtok(args_copy, " \t");
    while (token != NULL) {
        if (strncmp(token, "maxlen=", 7) == 0) {
            if (!ctools_safe_atoi(token + 7, &maxlen)) {
                free(args_copy);
                SF_error("cdecode: invalid maxlen value\n");
                return 198;
            }
        }
        token = strtok(NULL, " \t");
    }

    /* Source should be numeric (check index 1 in plugin's varlist) */
    if (SF_var_is_string(src_idx)) {
        free(args_copy);
        SF_error("cdecode: source variable must be numeric\n");
        return 198;
    }

    t_parse = ctools_timer_seconds();

    /* Build hash table from labels file */
    cdecode_hash_table ht;
    if (cdecode_hash_init(&ht, CDECODE_HASH_INIT_SIZE) != 0) {
        free(args_copy);
        SF_error("cdecode: memory allocation failed for hash table\n");
        return 920;
    }

    int detected_maxlen = 0;
    if (labelsfile[0] != '\0') {
        if (parse_labels_from_file(labelsfile, &ht, &detected_maxlen) != 0) {
            cdecode_hash_free(&ht);
            free(args_copy);
            SF_error("cdecode: failed to parse labels file\n");
            return 920;
        }
    }

    free(args_copy);
    args_copy = NULL;

    /* Use detected maxlen if not specified */
    if (maxlen == 0) {
        maxlen = detected_maxlen > 0 ? detected_maxlen : 1;
    }

    /* ========================================================================
     * Bulk load source numeric variable with if/in filtering
     * Uses ctools_data_load() to load only filtered observations.
     * ======================================================================== */
    int var_indices[1] = { src_idx };
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    stata_retcode load_rc = ctools_data_load(&filtered, var_indices, 1, 0, 0, 0);

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        cdecode_hash_free(&ht);
        SF_error("cdecode: failed to load source data\n");
        return 920;
    }

    nobs = filtered.data.nobs;
    perm_idx_t *obs_map = filtered.obs_map;

    if (nobs == 0) {
        ctools_filtered_data_free(&filtered);
        cdecode_hash_free(&ht);
        SF_scal_save("_cdecode_n_decoded", 0);
        SF_scal_save("_cdecode_n_missing", 0);
        SF_scal_save("_cdecode_n_unlabeled", 0);
        return 0;
    }

    double *src_values = filtered.data.vars[0].data.dbl;

    /* Counters */
    int n_decoded = 0;
    int n_missing = 0;
    int n_unlabeled = 0;

    /* Process filtered observations - lookups from loaded data, string writes use obs_map */
    for (size_t i = 0; i < nobs; i++) {
        ST_int obs = (ST_int)obs_map[i];  /* Get Stata observation from obs_map */

        double dval = src_values[i];

        if (SF_is_missing(dval)) {
            /* Missing value -> empty string */
            SF_sstore(dst_idx, obs, "");
            n_missing++;
            continue;
        }

        /* Check if value is close to an integer (like native decode) */
        double rounded = round(dval);
        if (fabs(dval - rounded) > 1e-9) {
            /* Non-integer value - no label match (matches native decode behavior) */
            SF_sstore(dst_idx, obs, "");
            n_unlabeled++;
            continue;
        }
        int ival = (int)rounded;

        /* Look up label */
        const char *label = cdecode_hash_lookup(&ht, ival);

        if (label) {
            SF_sstore(dst_idx, obs, (char *)label);
            n_decoded++;
        } else {
            /* No label found - store empty string */
            SF_sstore(dst_idx, obs, "");
            n_unlabeled++;
        }
    }

    /* Free loaded data */
    ctools_filtered_data_free(&filtered);

    t_decode = ctools_timer_seconds();
    t_total = t_decode - t_start;

    /* Save count before freeing */
    size_t n_labels = ht.count;
    cdecode_hash_free(&ht);

    /* Store results */
    SF_scal_save("_cdecode_n_decoded", (double)n_decoded);
    SF_scal_save("_cdecode_n_missing", (double)n_missing);
    SF_scal_save("_cdecode_n_unlabeled", (double)n_unlabeled);
    SF_scal_save("_cdecode_n_labels", (double)n_labels);
    SF_scal_save("_cdecode_max_len", (double)detected_maxlen);
    SF_scal_save("_cdecode_time_parse", t_parse - t_start);
    SF_scal_save("_cdecode_time_decode", t_decode - t_parse);
    SF_scal_save("_cdecode_time_total", t_total);

#ifdef _OPENMP
    SF_scal_save("_cdecode_openmp_enabled", 1.0);
    SF_scal_save("_cdecode_threads_max", (double)omp_get_max_threads());
#else
    SF_scal_save("_cdecode_openmp_enabled", 0.0);
    SF_scal_save("_cdecode_threads_max", 1.0);
#endif

    return rc;
}

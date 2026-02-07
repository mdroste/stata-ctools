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
#include <math.h>

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_runtime.h"
#include "ctools_config.h"
#include "ctools_hash.h"
#include "ctools_parse.h"
#include "cdecode_impl.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CDECODE_HASH_INIT_SIZE  256

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

    ctools_parse_string_option(args, "labelsfile", labelsfile, sizeof(labelsfile));
    ctools_parse_int_option(args, "maxlen", &maxlen);

    /* Source should be numeric (check index 1 in plugin's varlist) */
    if (SF_var_is_string(src_idx)) {
        SF_error("cdecode: source variable must be numeric\n");
        return 198;
    }

    t_parse = ctools_timer_seconds();

    /* Build hash table from labels file */
    ctools_int_hash_table ht;
    if (ctools_int_hash_init(&ht, CDECODE_HASH_INIT_SIZE) != 0) {
        SF_error("cdecode: memory allocation failed for hash table\n");
        return 920;
    }

    int detected_maxlen = 0;
    if (labelsfile[0] != '\0') {
        if (ctools_label_parse_file(labelsfile, &ht, &detected_maxlen) != 0) {
            ctools_int_hash_free(&ht);
            SF_error("cdecode: failed to parse labels file\n");
            return 920;
        }
    }

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
        ctools_int_hash_free(&ht);
        SF_error("cdecode: failed to load source data\n");
        return 920;
    }

    nobs = filtered.data.nobs;
    perm_idx_t *obs_map = filtered.obs_map;

    if (nobs == 0) {
        ctools_filtered_data_free(&filtered);
        ctools_int_hash_free(&ht);
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

    /* Decode and write directly to Stata — no intermediate array needed */
    for (size_t i = 0; i < nobs; i++) {
        double dval = src_values[i];
        const char *label = "";

        if (SF_is_missing(dval)) {
            n_missing++;
        } else {
            /* Check if value is close to an integer (like native decode) */
            double rounded = round(dval);
            if (fabs(dval - rounded) > 1e-9) {
                n_unlabeled++;
            } else {
                const char *found = ctools_int_hash_lookup(&ht, (int)rounded);
                if (found) {
                    label = found;
                    n_decoded++;
                } else {
                    n_unlabeled++;
                }
            }
        }

        SF_sstore(dst_idx, (ST_int)obs_map[i], (char *)label);
    }

    /* Free loaded data */
    ctools_filtered_data_free(&filtered);

    t_decode = ctools_timer_seconds();
    t_total = t_decode - t_start;

    /* Save count before freeing */
    size_t n_labels = ht.count;
    ctools_int_hash_free(&ht);

    /* Store results */
    SF_scal_save("_cdecode_n_decoded", (double)n_decoded);
    SF_scal_save("_cdecode_n_missing", (double)n_missing);
    SF_scal_save("_cdecode_n_unlabeled", (double)n_unlabeled);
    SF_scal_save("_cdecode_n_labels", (double)n_labels);
    SF_scal_save("_cdecode_max_len", (double)detected_maxlen);
    SF_scal_save("_cdecode_time_parse", t_parse - t_start);
    SF_scal_save("_cdecode_time_decode", t_decode - t_parse);
    SF_scal_save("_cdecode_time_total", t_total);

    CTOOLS_SAVE_THREAD_INFO("_cdecode");

    return rc;
}

/*
 * cencode_impl.c
 *
 * High-performance string encoding for Stata
 * Part of the ctools suite
 *
 * Algorithm:
 * 1. Bulk load source string variable via ctools_data_load_selective()
 * 2. Build hash table of unique values from loaded strings
 * 3. Sort unique strings alphabetically, assign codes 1, 2, 3, ...
 * 4. Map original strings to sorted codes, write numeric output
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_timer.h"
#include "ctools_config.h"
#include "ctools_arena.h"
#include "ctools_hash.h"
#include "cencode_impl.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CENCODE_HASH_INIT_SIZE  1024
#define CENCODE_HASH_LOAD_FACTOR 0.75
#define CENCODE_MAX_LABELS      65536
#define CENCODE_STR_BUF_SIZE    2048
#define CENCODE_EXISTING_BUF    32000

/* ============================================================================
 * Hash Table for String -> Integer Mapping
 * Uses unified ctools_str_hash_table from ctools_hash.h
 * ============================================================================ */

/* Compatibility typedefs and macros for cencode */
typedef ctools_str_hash_table cencode_hash_table;

#define cencode_hash_string      ctools_str_hash_compute
#define cencode_hash_init        ctools_str_hash_init
#define cencode_hash_free        ctools_str_hash_free
#define cencode_hash_insert      ctools_str_hash_insert
#define cencode_hash_insert_value ctools_str_hash_insert_value
#define cencode_hash_lookup      ctools_str_hash_lookup

/* Parse existing labels from format: "value1|string1||value2|string2||..." */
static int cencode_parse_existing_labels(cencode_hash_table *ht, const char *existing_str)
{
    if (!existing_str || *existing_str == '\0') return 0;

    size_t len = strlen(existing_str);
    char *buf = malloc(len + 1);
    if (!buf) return -1;
    memcpy(buf, existing_str, len + 1);  /* Include null terminator */

    char *pos = buf;
    while (*pos) {
        /* Parse value */
        char *pipe = strchr(pos, '|');
        if (!pipe) break;
        *pipe = '\0';
        int value;
        if (!ctools_safe_atoi(pos, &value)) {
            free(buf);
            return -1;  /* Invalid value in existing labels */
        }
        pos = pipe + 1;

        /* Find end of string (|| or end of buffer) */
        char *end = strstr(pos, "||");
        if (end) {
            *end = '\0';
        }

        /* Unescape the string */
        char *dst = pos;
        char *src = pos;
        while (*src) {
            if (*src == '\\' && (*(src+1) == '|' || *(src+1) == '\\')) {
                src++;
            }
            *dst++ = *src++;
        }
        *dst = '\0';

        /* Insert into hash table */
        if (cencode_hash_insert_value(ht, pos, value) < 0) {
            free(buf);
            return -1;
        }

        if (end) {
            pos = end + 2;  /* Skip || */
        } else {
            break;
        }
    }

    free(buf);
    return 0;
}

/* ============================================================================
 * Sorted String List
 * ============================================================================ */

typedef struct {
    char *str;
    int original_code;
    int final_code;
} cencode_string_entry;

static int cencode_string_compare(const void *a, const void *b)
{
    const cencode_string_entry *ea = (const cencode_string_entry *)a;
    const cencode_string_entry *eb = (const cencode_string_entry *)b;
    return strcmp(ea->str, eb->str);
}

/* ============================================================================
 * Main Implementation
 * ============================================================================ */

ST_retcode cencode_main(const char *args)
{
    double t_start, t_parse, t_collect, t_sort, t_encode, t_labels, t_total;
    int var_idx = 0, gen_idx = 0;
    char label_name[256] = "";
    int noextend = 0;
    char *existing_labels_str = NULL;
    size_t nobs;
    ST_retcode rc = 0;

    /*
     * cencode is called with ALL dataset variables in the varlist:
     *   plugin call ctools_plugin allvars, "cencode var_idx gen_idx ..."
     *
     * This makes varlist indices equal to dataset column indices.
     * var_idx = source string variable's position in dataset
     * gen_idx = destination numeric variable's position in dataset
     */

    t_start = ctools_timer_seconds();

    /* Parse arguments */
    if (args == NULL || *args == '\0') {
        SF_error("cencode: no arguments specified\n");
        return 198;
    }

    /* Make a larger copy for existing= which can be very long */
    size_t args_len = strlen(args);
    char *args_copy = malloc(args_len + 1);
    if (!args_copy) {
        SF_error("cencode: memory allocation failed\n");
        return 920;
    }
    memcpy(args_copy, args, args_len + 1);  /* Include null terminator */

    /* Parse arguments: first two numbers are var_idx and gen_idx */
    char *token = strtok(args_copy, " ");
    if (token && (token[0] >= '0' && token[0] <= '9')) {
        if (!ctools_safe_atoi(token, &var_idx)) {
            free(args_copy);
            SF_error("cencode: invalid var_idx value\n");
            return 198;
        }
        token = strtok(NULL, " ");
    }
    if (token && (token[0] >= '0' && token[0] <= '9')) {
        if (!ctools_safe_atoi(token, &gen_idx)) {
            free(args_copy);
            SF_error("cencode: invalid gen_idx value\n");
            return 198;
        }
        token = strtok(NULL, " ");
    }

    /* Parse remaining optional arguments */
    while (token != NULL) {
        if (strncmp(token, "label=", 6) == 0) {
            strncpy(label_name, token + 6, sizeof(label_name) - 1);
            label_name[sizeof(label_name) - 1] = '\0';
        } else if (strcmp(token, "noextend") == 0) {
            noextend = 1;
        } else if (strncmp(token, "existing=", 9) == 0) {
            /* existing= contains label mappings: value1|string1||value2|string2||... */
            existing_labels_str = token + 9;
        }
        token = strtok(NULL, " ");
    }

    /* Build existing labels hash table if noextend with existing labels */
    cencode_hash_table existing_ht;
    int have_existing = 0;
    if (noextend && existing_labels_str && *existing_labels_str) {
        if (cencode_hash_init(&existing_ht, CENCODE_HASH_INIT_SIZE) != 0) {
            free(args_copy);
            SF_error("cencode: memory allocation failed for existing labels\n");
            return 920;
        }
        if (cencode_parse_existing_labels(&existing_ht, existing_labels_str) != 0) {
            cencode_hash_free(&existing_ht);
            free(args_copy);
            SF_error("cencode: failed to parse existing labels\n");
            return 920;
        }
        have_existing = 1;
    }

    if (var_idx < 1 || gen_idx < 1) {
        if (have_existing) cencode_hash_free(&existing_ht);
        free(args_copy);
        SF_error("cencode: invalid variable indices\n");
        return 198;
    }

    /* Check source is a string variable */
    if (!SF_var_is_string(var_idx)) {
        if (have_existing) cencode_hash_free(&existing_ht);
        free(args_copy);
        SF_error("cencode: source variable must be a string variable\n");
        return 198;
    }

    ST_int obs1 = SF_in1();
    ST_int obs2 = SF_in2();
    nobs = (size_t)(obs2 - obs1 + 1);

    if (nobs == 0) {
        /* Empty dataset - match Stata's encode behavior (return success) */
        if (have_existing) cencode_hash_free(&existing_ht);
        free(args_copy);
        SF_scal_save("_cencode_n_unique", 0);
        SF_scal_save("_cencode_n_chunks", 1);
        SF_macro_save("_cencode_labels_0", "");
        SF_scal_save("_cencode_time_parse", ctools_timer_seconds() - t_start);
        SF_scal_save("_cencode_time_load", 0);
        SF_scal_save("_cencode_time_collect", 0);
        SF_scal_save("_cencode_time_sort", 0);
        SF_scal_save("_cencode_time_encode", 0);
        SF_scal_save("_cencode_time_labels", 0);
        SF_scal_save("_cencode_time_total", ctools_timer_seconds() - t_start);
        return 0;
    }

    t_parse = ctools_timer_seconds();
    double time_parse = t_parse - t_start;

    /* ========================================================================
     * Bulk load source string variable using ctools_data_load_selective()
     * ======================================================================== */
    int src_var_indices[1] = { var_idx };
    stata_data src_data;
    stata_data_init(&src_data);
    stata_retcode load_rc = ctools_data_load_selective(&src_data, src_var_indices, 1, obs1, obs2);

    if (load_rc != STATA_OK) {
        if (have_existing) cencode_hash_free(&existing_ht);
        free(args_copy);
        SF_error("cencode: failed to load source string data\n");
        return 920;
    }

    char **src_strings = src_data.vars[0].data.str;

    /* Build if-condition mask (sequential — SPI calls) */
    unsigned char *if_mask = malloc(nobs);
    if (!if_mask) {
        stata_data_free(&src_data);
        if (have_existing) cencode_hash_free(&existing_ht);
        free(args_copy);
        SF_error("cencode: memory allocation failed\n");
        return 920;
    }
    for (size_t i = 0; i < nobs; i++) {
        if_mask[i] = SF_ifobs((ST_int)(obs1 + (ST_int)i)) ? 1 : 0;
    }

    double time_load = ctools_timer_seconds() - t_parse;

    /* ========================================================================
     * Process loaded strings: collect unique values AND store codes
     * ======================================================================== */

    cencode_hash_table ht;
    int need_ht = !have_existing;  /* Only need hash table if not using existing labels */
    if (need_ht) {
        if (cencode_hash_init(&ht, CENCODE_HASH_INIT_SIZE) != 0) {
            free(if_mask);
            stata_data_free(&src_data);
            if (have_existing) cencode_hash_free(&existing_ht);
            free(args_copy);
            SF_error("cencode: memory allocation failed for hash table\n");
            return 920;
        }
    }

    /* Allocate array to store codes for each observation */
    int *obs_codes = malloc(nobs * sizeof(int));
    if (!obs_codes) {
        if (need_ht) cencode_hash_free(&ht);
        free(if_mask);
        stata_data_free(&src_data);
        if (have_existing) cencode_hash_free(&existing_ht);
        free(args_copy);
        SF_error("cencode: memory allocation failed for observation codes\n");
        return 920;
    }
    memset(obs_codes, 0, nobs * sizeof(int));  /* 0 = missing/empty */

    int collect_error = 0;

    for (size_t i = 0; i < nobs && !collect_error; i++) {
        /* Skip observations that don't meet the if condition */
        if (!if_mask[i]) {
            obs_codes[i] = 0;  /* not in sample */
            continue;
        }

        char *str_ptr = src_strings[i];

        /* Empty or NULL string -> missing */
        if (!str_ptr || str_ptr[0] == '\0') {
            obs_codes[i] = 0;  /* missing */
            continue;
        }

        if (have_existing) {
            /* Use existing label mapping - look up string in existing_ht */
            int code = cencode_hash_lookup(&existing_ht, str_ptr);
            obs_codes[i] = code;  /* 0 if not found (will become missing) */
        } else {
            /* Normal path - collect unique strings */
            uint32_t hash = cencode_hash_string(str_ptr);
            int code = cencode_hash_insert(&ht, str_ptr, hash);
            if (code < 0) {
                collect_error = 1;
            } else {
                obs_codes[i] = code;  /* Store the original code */
            }
        }
    }

    /* Free if_mask - no longer needed */
    free(if_mask);

    if (collect_error) {
        free(obs_codes);
        if (need_ht) cencode_hash_free(&ht);
        stata_data_free(&src_data);
        if (have_existing) cencode_hash_free(&existing_ht);
        free(args_copy);
        SF_error("cencode: memory allocation failed during collection\n");
        return 920;
    }

    t_collect = ctools_timer_seconds();
    double time_collect = t_collect - t_parse - time_load;

    /* ========================================================================
     * Handle existing labels path (noextend with pre-existing label)
     * ======================================================================== */

    if (have_existing) {
        /* Write encoded values directly - obs_codes already contain final codes */
        for (size_t i = 0; i < nobs; i++) {
            ST_int obs = obs1 + (ST_int)i;
            int code = obs_codes[i];

            if (code > 0) {
                SF_vstore(gen_idx, obs, (double)code);
            } else {
                SF_vstore(gen_idx, obs, SV_missval);
            }
        }

        free(obs_codes);
        stata_data_free(&src_data);
        cencode_hash_free(&existing_ht);
        free(args_copy);

        t_total = ctools_timer_seconds();

        /* No new labels to define - using existing */
        SF_scal_save("_cencode_n_unique", 0);  /* Signal no new labels needed */
        SF_scal_save("_cencode_n_chunks", 1);
        SF_macro_save("_cencode_labels_0", "");

        SF_scal_save("_cencode_time_parse", time_parse);
        SF_scal_save("_cencode_time_load", time_load);
        SF_scal_save("_cencode_time_collect", time_collect);
        SF_scal_save("_cencode_time_sort", 0);
        SF_scal_save("_cencode_time_encode", t_total - t_collect);
        SF_scal_save("_cencode_time_labels", 0);
        SF_scal_save("_cencode_time_total", t_total - t_start);

#ifdef _OPENMP
        SF_scal_save("_cencode_openmp_enabled", 1.0);
        SF_scal_save("_cencode_threads_max", (double)omp_get_max_threads());
#else
        SF_scal_save("_cencode_openmp_enabled", 0.0);
        SF_scal_save("_cencode_threads_max", 1.0);
#endif

        return 0;
    }

    /* ========================================================================
     * Normal path - sort and create new labels
     * ======================================================================== */

    size_t n_unique = ht.count;

    if (n_unique == 0) {
        for (size_t i = 0; i < nobs; i++) {
            SF_vstore(gen_idx, obs1 + (ST_int)i, SV_missval);
        }
        SF_scal_save("_cencode_n_unique", 0);
        SF_scal_save("_cencode_n_chunks", 1);
        SF_macro_save("_cencode_labels_0", "");

        free(obs_codes);
        stata_data_free(&src_data);
        cencode_hash_free(&ht);
        free(args_copy);

        t_total = ctools_timer_seconds();
        SF_scal_save("_cencode_time_parse", time_parse);
        SF_scal_save("_cencode_time_load", time_load);
        SF_scal_save("_cencode_time_collect", time_collect);
        SF_scal_save("_cencode_time_sort", 0);
        SF_scal_save("_cencode_time_encode", 0);
        SF_scal_save("_cencode_time_labels", 0);
        SF_scal_save("_cencode_time_total", t_total - t_start);
        return 0;
    }

    if (n_unique > CENCODE_MAX_LABELS) {
        char msg[256];
        snprintf(msg, sizeof(msg), "cencode: too many unique values (%zu, max %d)\n",
                 n_unique, CENCODE_MAX_LABELS);
        SF_error(msg);
        free(obs_codes);
        stata_data_free(&src_data);
        cencode_hash_free(&ht);
        free(args_copy);
        return 198;
    }

    /* ========================================================================
     * Sort unique strings and build code mapping
     * ======================================================================== */

    cencode_string_entry *sorted_strings = malloc(n_unique * sizeof(cencode_string_entry));
    if (!sorted_strings) {
        free(obs_codes);
        stata_data_free(&src_data);
        cencode_hash_free(&ht);
        free(args_copy);
        SF_error("cencode: memory allocation failed\n");
        return 920;
    }

    size_t idx = 0;
    for (size_t i = 0; i < ht.capacity && idx < n_unique; i++) {
        if (ht.entries[i].key) {
            sorted_strings[idx].str = ht.entries[i].key;
            sorted_strings[idx].original_code = ht.entries[i].value;
            idx++;
        }
    }

    qsort(sorted_strings, n_unique, sizeof(cencode_string_entry), cencode_string_compare);

    /* Check for overflow in (n_unique + 1) */
    if (n_unique >= SIZE_MAX) {
        free(sorted_strings);
        free(obs_codes);
        stata_data_free(&src_data);
        cencode_hash_free(&ht);
        free(args_copy);
        SF_error("cencode: too many unique values\n");
        return 920;
    }

    /* Use overflow-safe allocation */
    int *code_map = ctools_safe_malloc2(n_unique + 1, sizeof(int));
    if (!code_map) {
        free(sorted_strings);
        free(obs_codes);
        stata_data_free(&src_data);
        cencode_hash_free(&ht);
        free(args_copy);
        SF_error("cencode: memory allocation failed\n");
        return 920;
    }
    memset(code_map, 0, (n_unique + 1) * sizeof(int));

    for (size_t i = 0; i < n_unique; i++) {
        sorted_strings[i].final_code = (int)(i + 1);
        code_map[sorted_strings[i].original_code] = sorted_strings[i].final_code;
    }

    t_sort = ctools_timer_seconds();
    double time_sort = t_sort - t_collect;

    /* ========================================================================
     * Pass 2: Apply code_map to stored codes, write encoded values
     * No need to re-read strings - we stored the original codes in Pass 1
     * ======================================================================== */

    for (size_t i = 0; i < nobs; i++) {
        ST_int obs = obs1 + (ST_int)i;
        int orig_code = obs_codes[i];

        if (orig_code > 0) {
            SF_vstore(gen_idx, obs, (double)code_map[orig_code]);
        } else {
            SF_vstore(gen_idx, obs, SV_missval);
        }
    }

    free(obs_codes);

    t_encode = ctools_timer_seconds();
    double time_encode = t_encode - t_sort;

    /* ========================================================================
     * Store label mapping in Stata macros
     * ======================================================================== */

    SF_scal_save("_cencode_n_unique", (double)n_unique);

    size_t chunk = 0;
    size_t labels_in_chunk = 0;
    char macro_name[64];
    char *macro_buf = malloc(32000);
    size_t macro_pos = 0;

    if (!macro_buf) {
        free(code_map);
        free(sorted_strings);
        free(obs_codes);
        stata_data_free(&src_data);
        cencode_hash_free(&ht);
        free(args_copy);
        SF_error("cencode: memory allocation failed\n");
        return 920;
    }

    macro_buf[0] = '\0';

    for (size_t i = 0; i < n_unique; i++) {
        int code = sorted_strings[i].final_code;
        const char *str = sorted_strings[i].str;

        char escaped[2048];
        size_t esc_pos = 0;
        for (const char *p = str; *p && esc_pos < sizeof(escaped) - 2; p++) {
            if (*p == '"' || *p == '\\' || *p == '|') {
                escaped[esc_pos++] = '\\';
            }
            escaped[esc_pos++] = *p;
        }
        escaped[esc_pos] = '\0';

        char entry[2200];
        int entry_len = snprintf(entry, sizeof(entry), "%s%d|%s",
                                 (labels_in_chunk > 0) ? "||" : "",
                                 code, escaped);

        if (macro_pos + entry_len >= 31000) {
            snprintf(macro_name, sizeof(macro_name), "_cencode_labels_%zu", chunk);
            SF_macro_save(macro_name, macro_buf);

            chunk++;
            labels_in_chunk = 0;
            macro_pos = 0;
            macro_buf[0] = '\0';

            entry_len = snprintf(entry, sizeof(entry), "%d|%s", code, escaped);
        }

        memcpy(macro_buf + macro_pos, entry, entry_len + 1);
        macro_pos += entry_len;
        labels_in_chunk++;
    }

    if (labels_in_chunk > 0) {
        snprintf(macro_name, sizeof(macro_name), "_cencode_labels_%zu", chunk);
        SF_macro_save(macro_name, macro_buf);
    }

    SF_scal_save("_cencode_n_chunks", (double)(chunk + 1));

    if (label_name[0] != '\0') {
        SF_macro_save("_cencode_label_name", label_name);
    }

    free(macro_buf);

    t_labels = ctools_timer_seconds();
    double time_labels = t_labels - t_encode;

    t_total = ctools_timer_seconds();
    double time_total = t_total - t_start;

    free(code_map);
    free(sorted_strings);
    stata_data_free(&src_data);
    cencode_hash_free(&ht);
    free(args_copy);

    /* Store timing */
    SF_scal_save("_cencode_time_parse", time_parse);
    SF_scal_save("_cencode_time_load", time_load);
    SF_scal_save("_cencode_time_collect", time_collect);
    SF_scal_save("_cencode_time_sort", time_sort);
    SF_scal_save("_cencode_time_encode", time_encode);
    SF_scal_save("_cencode_time_labels", time_labels);
    SF_scal_save("_cencode_time_total", time_total);

#ifdef _OPENMP
    SF_scal_save("_cencode_openmp_enabled", 1.0);
    SF_scal_save("_cencode_threads_max", (double)omp_get_max_threads());
#else
    SF_scal_save("_cencode_openmp_enabled", 0.0);
    SF_scal_save("_cencode_threads_max", 1.0);
#endif

    return rc;
}

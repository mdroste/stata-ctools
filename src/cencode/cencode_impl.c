/*
 * cencode_impl.c
 *
 * High-performance string encoding for Stata
 * Part of the ctools suite
 *
 * Algorithm:
 * 1. Bulk load source string variable via ctools_data_load()
 * 2. Build hash table of unique values from loaded strings
 * 3. Sort unique strings alphabetically, assign codes 1, 2, 3, ...
 * 4. Map original strings to sorted codes, write numeric output
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_runtime.h"
#include "ctools_config.h"
#include "ctools_hash.h"
#include "ctools_parse.h"
#include "cencode_impl.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CENCODE_HASH_INIT_SIZE  1024

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
    size_t nobs;

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

    /* First two tokens are positional: var_idx gen_idx */
    const char *cursor = args;
    if (ctools_parse_next_int(&cursor, &var_idx) != 0) {
        SF_error("cencode: invalid var_idx value\n");
        return 198;
    }
    if (ctools_parse_next_int(&cursor, &gen_idx) != 0) {
        SF_error("cencode: invalid gen_idx value\n");
        return 198;
    }

    /* Named options */
    ctools_parse_string_option(args, "label", label_name, sizeof(label_name));
    noextend = ctools_parse_bool_option(args, "noextend");

    char existfile[1024] = "";
    ctools_parse_string_option(args, "existfile", existfile, sizeof(existfile));

    char labelfile[1024] = "";
    ctools_parse_string_option(args, "labelfile", labelfile, sizeof(labelfile));

    /* Build existing labels hash table from label save file (noextend path) */
    ctools_str_hash_table existing_ht;
    int have_existing = 0;
    if (noextend && existfile[0] != '\0') {
        if (ctools_str_hash_init(&existing_ht, CENCODE_HASH_INIT_SIZE) != 0) {
            SF_error("cencode: memory allocation failed for existing labels\n");
            return 920;
        }
        if (ctools_label_parse_stata_file_str(existfile, &existing_ht) != 0) {
            ctools_str_hash_free(&existing_ht);
            SF_error("cencode: failed to parse existing labels\n");
            return 920;
        }
        have_existing = 1;
    }

    if (var_idx < 1 || gen_idx < 1) {
        if (have_existing) ctools_str_hash_free(&existing_ht);
        SF_error("cencode: invalid variable indices\n");
        return 198;
    }

    /* Check source is a string variable */
    if (!SF_var_is_string(var_idx)) {
        if (have_existing) ctools_str_hash_free(&existing_ht);
        SF_error("cencode: source variable must be a string variable\n");
        return 198;
    }

    t_parse = ctools_timer_seconds();
    double time_parse = t_parse - t_start;

    /* ========================================================================
     * Bulk load source string variable with if/in filtering
     * Uses ctools_data_load() to load only filtered observations,
     * eliminating the need for a separate if_mask array.
     * ======================================================================== */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    int load_var_idx = var_idx;
    stata_retcode load_rc = ctools_data_load(&filtered, &load_var_idx, 1, 0, 0, 0);

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        if (have_existing) ctools_str_hash_free(&existing_ht);
        SF_error("cencode: failed to load source string data\n");
        return 920;
    }

    nobs = filtered.data.nobs;  /* Number of filtered observations */
    perm_idx_t *obs_map = filtered.obs_map;  /* Maps filtered index -> Stata obs */

    if (nobs == 0) {
        /* No observations match if/in - match Stata's encode behavior (return success) */
        ctools_filtered_data_free(&filtered);
        if (have_existing) ctools_str_hash_free(&existing_ht);
        SF_scal_save("_cencode_n_unique", 0);
        /* (no label file or macros needed for this early-return path) */
        SF_scal_save("_cencode_time_parse", time_parse);
        SF_scal_save("_cencode_time_load", ctools_timer_seconds() - t_parse);
        SF_scal_save("_cencode_time_collect", 0);
        SF_scal_save("_cencode_time_sort", 0);
        SF_scal_save("_cencode_time_encode", 0);
        SF_scal_save("_cencode_time_labels", 0);
        SF_scal_save("_cencode_time_total", ctools_timer_seconds() - t_start);
        return 0;
    }

    char **src_strings = filtered.data.vars[0].data.str;

    double time_load = ctools_timer_seconds() - t_parse;

    /* ========================================================================
     * Process loaded strings: collect unique values AND store codes
     * ======================================================================== */

    ctools_str_hash_table ht;
    int need_ht = !have_existing;  /* Only need hash table if not using existing labels */
    if (need_ht) {
        if (ctools_str_hash_init(&ht, CENCODE_HASH_INIT_SIZE) != 0) {
            ctools_filtered_data_free(&filtered);
            if (have_existing) ctools_str_hash_free(&existing_ht);
            SF_error("cencode: memory allocation failed for hash table\n");
            return 920;
        }
    }

    /* Allocate array to store codes for each filtered observation */
    int *obs_codes = malloc(nobs * sizeof(int));
    if (!obs_codes) {
        if (need_ht) ctools_str_hash_free(&ht);
        ctools_filtered_data_free(&filtered);
        if (have_existing) ctools_str_hash_free(&existing_ht);
        SF_error("cencode: memory allocation failed for observation codes\n");
        return 920;
    }
    int collect_error = 0;

    /* Process filtered observations - no if_mask check needed since data is pre-filtered */
    for (size_t i = 0; i < nobs && !collect_error; i++) {
        char *str_ptr = src_strings[i];

        /* Empty or NULL string -> missing */
        if (!str_ptr || str_ptr[0] == '\0') {
            obs_codes[i] = 0;  /* missing */
            continue;
        }

        if (have_existing) {
            /* Use existing label mapping - look up string in existing_ht */
            int code = ctools_str_hash_lookup(&existing_ht, str_ptr);
            obs_codes[i] = code;  /* 0 if not found (will become missing) */
        } else {
            /* Normal path - collect unique strings */
            uint32_t hash = ctools_str_hash_compute(str_ptr);
            int code = ctools_str_hash_insert(&ht, str_ptr, hash);
            if (code < 0) {
                collect_error = 1;
            } else {
                obs_codes[i] = code;  /* Store the original code */
            }
        }
    }

    if (collect_error) {
        free(obs_codes);
        if (need_ht) ctools_str_hash_free(&ht);
        ctools_filtered_data_free(&filtered);
        if (have_existing) ctools_str_hash_free(&existing_ht);
        SF_error("cencode: memory allocation failed during collection\n");
        return 920;
    }

    t_collect = ctools_timer_seconds();
    double time_collect = t_collect - t_parse - time_load;

    /* ========================================================================
     * Handle existing labels path (noextend with pre-existing label)
     * ======================================================================== */

    if (have_existing) {
        /* Build double array and use row-parallel store */
        double *enc_values = malloc(nobs * sizeof(double));
        if (!enc_values) {
            free(obs_codes);
            ctools_filtered_data_free(&filtered);
            ctools_str_hash_free(&existing_ht);
            SF_error("cencode: memory allocation failed\n");
            return 920;
        }
        for (size_t i = 0; i < nobs; i++) {
            int code = obs_codes[i];
            enc_values[i] = (code > 0) ? (double)code : SV_missval;
        }
        ctools_store_filtered_rowpar(enc_values, nobs, gen_idx, obs_map);
        free(enc_values);

        free(obs_codes);
        ctools_filtered_data_free(&filtered);
        ctools_str_hash_free(&existing_ht);

        t_total = ctools_timer_seconds();

        /* No new labels to define - using existing */
        SF_scal_save("_cencode_n_unique", 0);  /* Signal no new labels needed */
        /* (no label file or macros needed for this early-return path) */

        SF_scal_save("_cencode_time_parse", time_parse);
        SF_scal_save("_cencode_time_load", time_load);
        SF_scal_save("_cencode_time_collect", time_collect);
        SF_scal_save("_cencode_time_sort", 0);
        SF_scal_save("_cencode_time_encode", t_total - t_collect);
        SF_scal_save("_cencode_time_labels", 0);
        SF_scal_save("_cencode_time_total", t_total - t_start);

        CTOOLS_SAVE_THREAD_INFO("_cencode");

        return 0;
    }

    /* ========================================================================
     * Normal path - sort and create new labels
     * ======================================================================== */

    size_t n_unique = ht.count;

    if (n_unique == 0) {
        /* All values were empty/missing - write all missing */
        double *miss_values = malloc(nobs * sizeof(double));
        if (miss_values) {
            for (size_t i = 0; i < nobs; i++) miss_values[i] = SV_missval;
            ctools_store_filtered_rowpar(miss_values, nobs, gen_idx, obs_map);
            free(miss_values);
        } else {
            for (size_t i = 0; i < nobs; i++)
                SF_vstore(gen_idx, (ST_int)obs_map[i], SV_missval);
        }
        SF_scal_save("_cencode_n_unique", 0);
        /* (no label file or macros needed for this early-return path) */

        free(obs_codes);
        ctools_filtered_data_free(&filtered);
        ctools_str_hash_free(&ht);

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

    if (n_unique > CTOOLS_MAX_LABELS) {
        char msg[256];
        snprintf(msg, sizeof(msg), "cencode: too many unique values (%zu, max %d)\n",
                 n_unique, CTOOLS_MAX_LABELS);
        SF_error(msg);
        free(obs_codes);
        ctools_filtered_data_free(&filtered);
        ctools_str_hash_free(&ht);
        return 198;
    }

    /* ========================================================================
     * Sort unique strings and build code mapping
     * ======================================================================== */

    cencode_string_entry *sorted_strings = malloc(n_unique * sizeof(cencode_string_entry));
    if (!sorted_strings) {
        free(obs_codes);
        ctools_filtered_data_free(&filtered);
        ctools_str_hash_free(&ht);
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
        ctools_filtered_data_free(&filtered);
        ctools_str_hash_free(&ht);
        SF_error("cencode: too many unique values\n");
        return 920;
    }

    /* Use overflow-safe allocation */
    int *code_map = ctools_safe_malloc2(n_unique + 1, sizeof(int));
    if (!code_map) {
        free(sorted_strings);
        free(obs_codes);
        ctools_filtered_data_free(&filtered);
        ctools_str_hash_free(&ht);
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
     * via ctools_store_filtered
     * ======================================================================== */

    /* Build double array with mapped codes and use row-parallel store */
    {
        double *enc_values = malloc(nobs * sizeof(double));
        if (!enc_values) {
            free(obs_codes);
            free(code_map);
            free(sorted_strings);
            ctools_filtered_data_free(&filtered);
            ctools_str_hash_free(&ht);
            SF_error("cencode: memory allocation failed\n");
            return 920;
        }
        for (size_t i = 0; i < nobs; i++) {
            int orig_code = obs_codes[i];
            enc_values[i] = (orig_code > 0) ? (double)code_map[orig_code] : SV_missval;
        }
        ctools_store_filtered_rowpar(enc_values, nobs, gen_idx, obs_map);
        free(enc_values);
    }

    free(obs_codes);

    t_encode = ctools_timer_seconds();
    double time_encode = t_encode - t_sort;

    /* ========================================================================
     * Write label definitions to .do file for the .ado to run
     * ======================================================================== */

    SF_scal_save("_cencode_n_unique", (double)n_unique);

    if (labelfile[0] != '\0') {
        /* Build arrays for file writing */
        const char **label_strings = malloc(n_unique * sizeof(const char *));
        int *label_codes = malloc(n_unique * sizeof(int));
        if (!label_strings || !label_codes) {
            free(label_strings);
            free(label_codes);
            free(code_map);
            free(sorted_strings);
            ctools_filtered_data_free(&filtered);
            ctools_str_hash_free(&ht);
            SF_error("cencode: memory allocation failed\n");
            return 920;
        }

        for (size_t i = 0; i < n_unique; i++) {
            label_strings[i] = sorted_strings[i].str;
            label_codes[i] = sorted_strings[i].final_code;
        }

        if (ctools_label_write_stata_file(label_strings, label_codes,
                                           n_unique, label_name,
                                           labelfile) != 0) {
            free(label_strings);
            free(label_codes);
            free(code_map);
            free(sorted_strings);
            ctools_filtered_data_free(&filtered);
            ctools_str_hash_free(&ht);
            SF_error("cencode: failed to write label file\n");
            return 920;
        }

        free(label_strings);
        free(label_codes);
    }

    t_labels = ctools_timer_seconds();
    double time_labels = t_labels - t_encode;

    t_total = ctools_timer_seconds();
    double time_total = t_total - t_start;

    free(code_map);
    free(sorted_strings);
    ctools_filtered_data_free(&filtered);
    ctools_str_hash_free(&ht);

    /* Store timing */
    SF_scal_save("_cencode_time_parse", time_parse);
    SF_scal_save("_cencode_time_load", time_load);
    SF_scal_save("_cencode_time_collect", time_collect);
    SF_scal_save("_cencode_time_sort", time_sort);
    SF_scal_save("_cencode_time_encode", time_encode);
    SF_scal_save("_cencode_time_labels", time_labels);
    SF_scal_save("_cencode_time_total", time_total);

    CTOOLS_SAVE_THREAD_INFO("_cencode");

    return 0;
}

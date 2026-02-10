/*
 * cdecode_impl.c
 *
 * High-performance numeric-to-string decoding for Stata
 * Part of the ctools suite
 *
 * Algorithm:
 * 1. Parse Stata `label save` format file directly into hash table
 * 2. Build flat lookup array for dense label values (fast path)
 * 3. Bulk load numeric source variable via ctools_data_load()
 * 4. For each observation, look up label and write to destination
 *    (skipping SF_sstore for missing/unlabeled values since dest is "")
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

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

/* Safe range for double-to-int conversion without undefined behavior */
#define CDECODE_INT_MIN_DBL (-2147483648.0)
#define CDECODE_INT_MAX_DBL  (2147483647.0)

/* ============================================================================
 * Scan-only entry point: returns max label length without decoding
 * ============================================================================ */

ST_retcode cdecode_scan_main(const char *args)
{
    char labelsfile[1024] = "";
    ctools_parse_string_option(args, "labelsfile", labelsfile, sizeof(labelsfile));

    if (labelsfile[0] == '\0') {
        SF_scal_save("_cdecode_maxlen", 1.0);
        return 0;
    }

    int maxlen = 0;
    if (ctools_label_scan_stata_file(labelsfile, &maxlen) != 0) {
        SF_error("cdecode_scan: failed to parse label file\n");
        return 920;
    }

    if (maxlen < 1) maxlen = 1;
    SF_scal_save("_cdecode_maxlen", (double)maxlen);
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
     * cdecode is called as: plugin call ctools_plugin src_var dst_var, args
     * So the source variable is always index 1, and destination is index 2.
     */
    const int src_idx = 1;
    const int dst_idx = 2;

    t_start = ctools_timer_seconds();

    /* Parse arguments */
    if (args == NULL || *args == '\0') {
        SF_error("cdecode: no arguments specified\n");
        return 198;
    }

    ctools_parse_string_option(args, "labelsfile", labelsfile, sizeof(labelsfile));
    ctools_parse_int_option(args, "maxlen", &maxlen);

    /* Source should be numeric */
    if (SF_var_is_string(src_idx)) {
        SF_error("cdecode: source variable must be numeric\n");
        return 198;
    }

    t_parse = ctools_timer_seconds();

    /* Build hash table directly from Stata label save file */
    ctools_int_hash_table ht;
    if (ctools_int_hash_init(&ht, CDECODE_HASH_INIT_SIZE) != 0) {
        SF_error("cdecode: memory allocation failed for hash table\n");
        return 920;
    }

    int detected_maxlen = 0;
    if (labelsfile[0] != '\0') {
        if (ctools_label_parse_stata_file_int(labelsfile, &ht, &detected_maxlen) != 0) {
            ctools_int_hash_free(&ht);
            SF_error("cdecode: failed to parse labels file\n");
            return 920;
        }
    }

    if (maxlen == 0) {
        maxlen = detected_maxlen > 0 ? detected_maxlen : 1;
    }

    /* ========================================================================
     * Build flat array for dense label values (fast path)
     *
     * When label values are dense (e.g., 1..K), a flat array gives O(1)
     * lookup without hash computation. Falls back to hash table for sparse
     * values (e.g., {1, 100, 10000}).
     * ======================================================================== */

    int use_flat = 0;
    const char **flat_labels = NULL;
    int flat_min_key = 0;
    int flat_max_key = 0;

    if (ht.count > 0) {
        int min_key = INT_MAX, max_key = INT_MIN;
        for (size_t i = 0; i < ht.capacity; i++) {
            if (ht.entries[i].occupied) {
                if (ht.entries[i].key < min_key) min_key = ht.entries[i].key;
                if (ht.entries[i].key > max_key) max_key = ht.entries[i].key;
            }
        }

        size_t range = (size_t)(max_key - min_key) + 1;
        /* Use flat array if range is reasonable: at most 4x count, under 1M */
        if (range <= 4 * ht.count && range <= 1000000) {
            flat_labels = calloc(range, sizeof(const char *));
            if (flat_labels) {
                for (size_t i = 0; i < ht.capacity; i++) {
                    if (ht.entries[i].occupied) {
                        flat_labels[ht.entries[i].key - min_key] = ht.entries[i].value;
                    }
                }
                flat_min_key = min_key;
                flat_max_key = max_key;
                use_flat = 1;
            }
        }
    }

    /* ========================================================================
     * Bulk load source numeric variable with if/in filtering
     * ======================================================================== */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    stata_retcode load_rc = ctools_data_load_single_var_rowpar(&filtered, src_idx, 0, 0, 0);

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        free(flat_labels);
        ctools_int_hash_free(&ht);
        SF_error("cdecode: failed to load source data\n");
        return 920;
    }

    nobs = filtered.data.nobs;
    perm_idx_t *obs_map = filtered.obs_map;

    if (nobs == 0) {
        ctools_filtered_data_free(&filtered);
        free(flat_labels);
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

    /* ========================================================================
     * Decode loop
     *
     * Two variants: flat array (no hashing) vs hash table (with cache).
     * Both skip SF_sstore for missing/unlabeled values since the
     * destination variable was initialized to "" by the .ado.
     *
     * Uses fast integer check: (int)dval == dval instead of round()+fabs().
     * Guarded by range check to avoid undefined behavior on overflow.
     * ======================================================================== */

    if (use_flat) {
        /* Fast path: direct array lookup, no hashing */
        for (size_t i = 0; i < nobs; i++) {
            double dval = src_values[i];

            if (SF_is_missing(dval)) {
                n_missing++;
                continue;
            }

            if (dval >= CDECODE_INT_MIN_DBL && dval <= CDECODE_INT_MAX_DBL) {
                int key = (int)dval;
                if (dval == (double)key &&
                    key >= flat_min_key && key <= flat_max_key) {
                    const char *found = flat_labels[key - flat_min_key];
                    if (found) {
                        SF_sstore(dst_idx, (ST_int)obs_map[i], (char *)found);
                        n_decoded++;
                        continue;
                    }
                }
            }
            n_unlabeled++;
        }
    } else {
        /* Hash table path with last-value cache */
        int cache_key = INT_MIN;
        const char *cache_label = NULL;

        for (size_t i = 0; i < nobs; i++) {
            double dval = src_values[i];

            if (SF_is_missing(dval)) {
                n_missing++;
                continue;
            }

            if (dval >= CDECODE_INT_MIN_DBL && dval <= CDECODE_INT_MAX_DBL) {
                int key = (int)dval;
                if (dval == (double)key) {
                    const char *label;
                    if (key == cache_key) {
                        label = cache_label;
                    } else {
                        label = ctools_int_hash_lookup(&ht, key);
                        cache_key = key;
                        cache_label = label;
                    }
                    if (label) {
                        SF_sstore(dst_idx, (ST_int)obs_map[i], (char *)label);
                        n_decoded++;
                        continue;
                    }
                }
            }
            n_unlabeled++;
        }
    }

    /* Free loaded data */
    ctools_filtered_data_free(&filtered);
    free(flat_labels);

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

/*
    cmerge_impl.c
    Optimized C-Accelerated Merge for Stata

    High-performance merge with MINIMAL data transfer:
    - Only loads key variables + keepusing variables into C
    - Master non-key variables are STREAMED (read-permute-write) without loading
    - Memory usage: O(nobs * (nkeys + nkeepusing)) instead of O(nobs * all_vars)

    Architecture (Two-Phase Plugin Call):
    Phase 1 - "load_using":
        - Load ONLY key + keepusing variables from using dataset
        - Sort using data on key variables
        - Store in static cache for phase 2

    Phase 2 - "execute":
        - Load ONLY master key variables + _orig_row index
        - Perform sorted merge join to produce output row mapping
        - Write _merge and keepusing variables directly
        - STREAM master non-key variables (parallel gather-scatter)

    Merge Types: 1:1, m:1, 1:m, m:m - all supported with streaming
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <pthread.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_timer.h"
#include "ctools_threads.h"
#include "cmerge_impl.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CMERGE_MAX_KEYVARS 32
#define CMERGE_MAX_VARS 1024

/* Merge type codes */
typedef enum {
    MERGE_1_1 = 0,
    MERGE_M_1 = 1,
    MERGE_1_M = 2,
    MERGE_M_M = 3
} cmerge_type_t;

/* Merge result codes (matches Stata's _merge values) */
typedef enum {
    MERGE_RESULT_MASTER_ONLY = 1,
    MERGE_RESULT_USING_ONLY = 2,
    MERGE_RESULT_BOTH = 3
} cmerge_result_t;

#define cmerge_get_time_ms() ctools_timer_ms()

/* ============================================================================
 * Output row specification - minimal per-row data
 * ============================================================================ */

typedef struct {
    int64_t master_sorted_row;  /* Row in SORTED master (-1 if using-only) */
    int64_t using_sorted_row;   /* Row in SORTED using (-1 if master-only) */
    int8_t merge_result;        /* 1, 2, or 3 */
} cmerge_output_spec_t;

/* ============================================================================
 * Using data cache (persists between plugin calls)
 * ============================================================================ */

typedef struct {
    stata_data keys;            /* Key variables only */
    stata_data keepusing;       /* Keepusing variables only */
    size_t nobs;
    int nkeys;
    int n_keepusing;
    int loaded;
} cmerge_using_cache_t;

static cmerge_using_cache_t g_using_cache = {0};

/* ============================================================================
 * Key comparison functions
 * ============================================================================ */

static int cmerge_compare_keys(stata_data *data_a, size_t row_a,
                                stata_data *data_b, size_t row_b,
                                int nkeys)
{
    for (int k = 0; k < nkeys; k++) {
        stata_variable *var_a = &data_a->vars[k];
        stata_variable *var_b = &data_b->vars[k];

        if (var_a->type == STATA_TYPE_DOUBLE) {
            double val_a = var_a->data.dbl[row_a];
            double val_b = var_b->data.dbl[row_b];

            int miss_a = SF_is_missing(val_a);
            int miss_b = SF_is_missing(val_b);

            if (miss_a && miss_b) continue;
            if (miss_a) return 1;
            if (miss_b) return -1;

            if (val_a < val_b) return -1;
            if (val_a > val_b) return 1;
        } else {
            char *str_a = var_a->data.str[row_a];
            char *str_b = var_b->data.str[row_b];

            if (str_a == NULL) str_a = "";
            if (str_b == NULL) str_b = "";

            int cmp = strcmp(str_a, str_b);
            if (cmp != 0) return (cmp < 0) ? -1 : 1;
        }
    }
    return 0;
}

static int cmerge_compare_keys_same(stata_data *data, size_t row_a, size_t row_b, int nkeys)
{
    return cmerge_compare_keys(data, row_a, data, row_b, nkeys);
}

/* ============================================================================
 * Sorted merge join - produces output specifications
 * ============================================================================ */

static int64_t cmerge_sorted_join(
    stata_data *master_keys, stata_data *using_keys,
    int nkeys, cmerge_type_t merge_type,
    cmerge_output_spec_t **output_specs_out)
{
    size_t m_nobs = master_keys->nobs;
    size_t u_nobs = using_keys->nobs;

    /* Initial allocation estimate */
    size_t capacity;
    switch (merge_type) {
        case MERGE_1_1:
            capacity = (m_nobs > u_nobs ? m_nobs : u_nobs) + (m_nobs < u_nobs ? m_nobs : u_nobs) / 4;
            break;
        case MERGE_M_1:
            capacity = m_nobs + u_nobs / 4;
            break;
        case MERGE_1_M:
            capacity = u_nobs + m_nobs / 4;
            break;
        case MERGE_M_M:
        default:
            capacity = m_nobs + u_nobs;
            break;
    }
    if (capacity < 1024) capacity = 1024;

    cmerge_output_spec_t *specs = malloc(capacity * sizeof(cmerge_output_spec_t));
    if (!specs) return -1;

    size_t count = 0;
    size_t m_idx = 0;
    size_t u_idx = 0;

    while (m_idx < m_nobs || u_idx < u_nobs) {
        /* Grow if needed */
        if (count >= capacity) {
            capacity *= 2;
            cmerge_output_spec_t *new_specs = realloc(specs, capacity * sizeof(cmerge_output_spec_t));
            if (!new_specs) { free(specs); return -1; }
            specs = new_specs;
        }

        if (m_idx >= m_nobs) {
            /* Master exhausted - remaining using rows */
            specs[count].master_sorted_row = -1;
            specs[count].using_sorted_row = (int64_t)u_idx;
            specs[count].merge_result = MERGE_RESULT_USING_ONLY;
            count++;
            u_idx++;
        }
        else if (u_idx >= u_nobs) {
            /* Using exhausted - remaining master rows */
            specs[count].master_sorted_row = (int64_t)m_idx;
            specs[count].using_sorted_row = -1;
            specs[count].merge_result = MERGE_RESULT_MASTER_ONLY;
            count++;
            m_idx++;
        }
        else {
            int cmp = cmerge_compare_keys(master_keys, m_idx, using_keys, u_idx, nkeys);

            if (cmp < 0) {
                /* Master only */
                specs[count].master_sorted_row = (int64_t)m_idx;
                specs[count].using_sorted_row = -1;
                specs[count].merge_result = MERGE_RESULT_MASTER_ONLY;
                count++;
                m_idx++;
            }
            else if (cmp > 0) {
                /* Using only */
                specs[count].master_sorted_row = -1;
                specs[count].using_sorted_row = (int64_t)u_idx;
                specs[count].merge_result = MERGE_RESULT_USING_ONLY;
                count++;
                u_idx++;
            }
            else {
                /* Keys match - count groups */
                size_t m_start = m_idx;
                size_t u_start = u_idx;

                size_t m_count = 1;
                while (m_idx + m_count < m_nobs &&
                       cmerge_compare_keys_same(master_keys, m_start, m_idx + m_count, nkeys) == 0) {
                    m_count++;
                }

                size_t u_count = 1;
                while (u_idx + u_count < u_nobs &&
                       cmerge_compare_keys_same(using_keys, u_start, u_idx + u_count, nkeys) == 0) {
                    u_count++;
                }

                /* Generate output based on merge type */
                switch (merge_type) {
                    case MERGE_1_1:
                        if (count >= capacity) {
                            capacity *= 2;
                            specs = realloc(specs, capacity * sizeof(cmerge_output_spec_t));
                            if (!specs) return -1;
                        }
                        specs[count].master_sorted_row = (int64_t)m_start;
                        specs[count].using_sorted_row = (int64_t)u_start;
                        specs[count].merge_result = MERGE_RESULT_BOTH;
                        count++;
                        break;

                    case MERGE_M_1:
                        /* Multiple master rows map to one using row */
                        for (size_t i = 0; i < m_count; i++) {
                            if (count >= capacity) {
                                capacity *= 2;
                                specs = realloc(specs, capacity * sizeof(cmerge_output_spec_t));
                                if (!specs) return -1;
                            }
                            specs[count].master_sorted_row = (int64_t)(m_start + i);
                            specs[count].using_sorted_row = (int64_t)u_start;
                            specs[count].merge_result = MERGE_RESULT_BOTH;
                            count++;
                        }
                        break;

                    case MERGE_1_M:
                        /* One master row maps to multiple using rows */
                        for (size_t i = 0; i < u_count; i++) {
                            if (count >= capacity) {
                                capacity *= 2;
                                specs = realloc(specs, capacity * sizeof(cmerge_output_spec_t));
                                if (!specs) return -1;
                            }
                            specs[count].master_sorted_row = (int64_t)m_start;
                            specs[count].using_sorted_row = (int64_t)(u_start + i);
                            specs[count].merge_result = MERGE_RESULT_BOTH;
                            count++;
                        }
                        break;

                    case MERGE_M_M:
                        /* Sequential pairing within groups */
                        {
                            size_t pairs = (m_count < u_count) ? m_count : u_count;
                            for (size_t i = 0; i < pairs; i++) {
                                if (count >= capacity) {
                                    capacity *= 2;
                                    specs = realloc(specs, capacity * sizeof(cmerge_output_spec_t));
                                    if (!specs) return -1;
                                }
                                specs[count].master_sorted_row = (int64_t)(m_start + i);
                                specs[count].using_sorted_row = (int64_t)(u_start + i);
                                specs[count].merge_result = MERGE_RESULT_BOTH;
                                count++;
                            }
                            /* Excess master rows */
                            for (size_t i = pairs; i < m_count; i++) {
                                if (count >= capacity) {
                                    capacity *= 2;
                                    specs = realloc(specs, capacity * sizeof(cmerge_output_spec_t));
                                    if (!specs) return -1;
                                }
                                specs[count].master_sorted_row = (int64_t)(m_start + i);
                                specs[count].using_sorted_row = -1;
                                specs[count].merge_result = MERGE_RESULT_MASTER_ONLY;
                                count++;
                            }
                            /* Excess using rows */
                            for (size_t i = pairs; i < u_count; i++) {
                                if (count >= capacity) {
                                    capacity *= 2;
                                    specs = realloc(specs, capacity * sizeof(cmerge_output_spec_t));
                                    if (!specs) return -1;
                                }
                                specs[count].master_sorted_row = -1;
                                specs[count].using_sorted_row = (int64_t)(u_start + i);
                                specs[count].merge_result = MERGE_RESULT_USING_ONLY;
                                count++;
                            }
                        }
                        break;
                }

                m_idx = m_start + m_count;
                u_idx = u_start + u_count;
            }
        }
    }

    *output_specs_out = specs;
    return (int64_t)count;
}

/* ============================================================================
 * Streaming thread for master non-key variables
 * ============================================================================ */

typedef struct {
    int stata_var_idx;              /* 1-based Stata variable index */
    int64_t *master_orig_rows;      /* master_orig_rows[i] = original row for output i */
    size_t output_nobs;
    int is_key;                     /* Is this a key variable? */
    int key_idx;                    /* If key, which key (0-based in cache) */
    cmerge_output_spec_t *specs;    /* For using key fallback */
    int success;
} stream_var_args_t;

static void *stream_master_var_thread(void *arg)
{
    stream_var_args_t *a = (stream_var_args_t *)arg;
    size_t output_nobs = a->output_nobs;
    ST_int var_idx = (ST_int)a->stata_var_idx;
    int64_t *src_rows = a->master_orig_rows;
    int is_string = SF_var_is_string(var_idx);

    a->success = 0;

    if (is_string) {
        char **buf = malloc(output_nobs * sizeof(char *));
        if (!buf) return (void *)1;

        char strbuf[2048];

        /* GATHER from original master positions */
        for (size_t i = 0; i < output_nobs; i++) {
            if (src_rows[i] >= 0) {
                SF_sdata(var_idx, (ST_int)(src_rows[i] + 1), strbuf);
                buf[i] = strdup(strbuf);
            } else if (a->is_key && g_using_cache.loaded) {
                /* Using-only row with key fallback */
                int64_t using_row = a->specs[i].using_sorted_row;
                if (using_row >= 0 && a->key_idx >= 0) {
                    char *s = g_using_cache.keys.vars[a->key_idx].data.str[using_row];
                    buf[i] = strdup(s ? s : "");
                } else {
                    buf[i] = strdup("");
                }
            } else {
                buf[i] = strdup("");
            }
            if (!buf[i]) {
                for (size_t k = 0; k < i; k++) free(buf[k]);
                free(buf);
                return (void *)1;
            }
        }

        /* SCATTER to output positions */
        for (size_t i = 0; i < output_nobs; i++) {
            SF_sstore(var_idx, (ST_int)(i + 1), buf[i]);
            free(buf[i]);
        }
        free(buf);

    } else {
        double *buf = malloc(output_nobs * sizeof(double));
        if (!buf) return (void *)1;

        /* GATHER from original master positions */
        for (size_t i = 0; i < output_nobs; i++) {
            if (src_rows[i] >= 0) {
                SF_vdata(var_idx, (ST_int)(src_rows[i] + 1), &buf[i]);
            } else if (a->is_key && g_using_cache.loaded) {
                /* Using-only row with key fallback */
                int64_t using_row = a->specs[i].using_sorted_row;
                if (using_row >= 0 && a->key_idx >= 0) {
                    buf[i] = g_using_cache.keys.vars[a->key_idx].data.dbl[using_row];
                } else {
                    buf[i] = SV_missval;
                }
            } else {
                buf[i] = SV_missval;
            }
        }

        /* SCATTER to output positions (unrolled) */
        size_t i_end = output_nobs - (output_nobs % 8);
        for (size_t i = 0; i < i_end; i += 8) {
            SF_vstore(var_idx, (ST_int)(i + 1), buf[i]);
            SF_vstore(var_idx, (ST_int)(i + 2), buf[i + 1]);
            SF_vstore(var_idx, (ST_int)(i + 3), buf[i + 2]);
            SF_vstore(var_idx, (ST_int)(i + 4), buf[i + 3]);
            SF_vstore(var_idx, (ST_int)(i + 5), buf[i + 4]);
            SF_vstore(var_idx, (ST_int)(i + 6), buf[i + 5]);
            SF_vstore(var_idx, (ST_int)(i + 7), buf[i + 6]);
            SF_vstore(var_idx, (ST_int)(i + 8), buf[i + 7]);
        }
        for (size_t i = i_end; i < output_nobs; i++) {
            SF_vstore(var_idx, (ST_int)(i + 1), buf[i]);
        }

        free(buf);
    }

    a->success = 1;
    return NULL;
}

/* ============================================================================
 * Phase 1: Load using dataset (keys + keepusing only)
 * ============================================================================ */

static ST_retcode cmerge_load_using(const char *args)
{
    char msg[256];
    ST_retcode rc;

    /* Parse: "load_using <nkeys> <key_idx1>...<key_idxN> n_keepusing <count>
              keepusing_indices <idx1>...<idxM> [verbose]" */
    char args_copy[4096];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    int nkeys = 0;
    int key_indices[CMERGE_MAX_KEYVARS];
    int n_keepusing = 0;
    int keepusing_indices[CMERGE_MAX_VARS];
    int verbose = 0;

    char *token = strtok(args_copy, " ");
    int arg_idx = 0;
    int in_keepusing_indices = 0;
    int keepusing_idx = 0;

    while (token != NULL) {
        if (arg_idx == 0) {
            nkeys = atoi(token);
            if (nkeys > CMERGE_MAX_KEYVARS) nkeys = CMERGE_MAX_KEYVARS;
        }
        else if (arg_idx <= nkeys) {
            key_indices[arg_idx - 1] = atoi(token);
        }
        else if (strcmp(token, "n_keepusing") == 0) {
            token = strtok(NULL, " ");
            if (token) n_keepusing = atoi(token);
        }
        else if (strcmp(token, "keepusing_indices") == 0) {
            in_keepusing_indices = 1;
            keepusing_idx = 0;
        }
        else if (in_keepusing_indices && keepusing_idx < n_keepusing) {
            keepusing_indices[keepusing_idx++] = atoi(token);
        }
        else if (strcmp(token, "verbose") == 0) {
            verbose = 1;
        }
        arg_idx++;
        token = strtok(NULL, " ");
    }

    /* Free any previous cache */
    if (g_using_cache.loaded) {
        stata_data_free(&g_using_cache.keys);
        stata_data_free(&g_using_cache.keepusing);
        g_using_cache.loaded = 0;
    }

    double t_start = cmerge_get_time_ms();

    /* Load ONLY key variables */
    if (verbose) {
        SF_display("cmerge: Loading using keys...\n");
    }

    rc = ctools_data_load_selective(&g_using_cache.keys, key_indices, nkeys, 0, 0);
    if (rc != STATA_OK) {
        SF_error("cmerge: Failed to load using keys\n");
        return 459;
    }

    /* Load ONLY keepusing variables */
    if (n_keepusing > 0) {
        if (verbose) {
            SF_display("cmerge: Loading keepusing variables...\n");
        }
        rc = ctools_data_load_selective(&g_using_cache.keepusing, keepusing_indices, n_keepusing, 0, 0);
        if (rc != STATA_OK) {
            stata_data_free(&g_using_cache.keys);
            SF_error("cmerge: Failed to load keepusing variables\n");
            return 459;
        }
    } else {
        stata_data_init(&g_using_cache.keepusing);
    }

    g_using_cache.nobs = g_using_cache.keys.nobs;
    g_using_cache.nkeys = nkeys;
    g_using_cache.n_keepusing = n_keepusing;

    double t_load = cmerge_get_time_ms() - t_start;

    /* Sort using data on keys */
    if (verbose) {
        SF_display("cmerge: Sorting using data...\n");
    }
    double t_sort_start = cmerge_get_time_ms();

    int *sort_vars = malloc(nkeys * sizeof(int));
    for (int i = 0; i < nkeys; i++) {
        sort_vars[i] = i + 1;  /* 1-based within keys structure */
    }

    /* Allocate permutation array to capture ordering before it's applied */
    size_t nobs = g_using_cache.nobs;
    size_t *perm = malloc(nobs * sizeof(size_t));
    if (!perm) {
        free(sort_vars);
        stata_data_free(&g_using_cache.keys);
        stata_data_free(&g_using_cache.keepusing);
        SF_error("cmerge: Memory allocation failed for permutation\n");
        return 920;
    }

    /* Use _with_perm version to get permutation BEFORE it's reset to identity */
    rc = ctools_sort_radix_lsd_with_perm(&g_using_cache.keys, sort_vars, nkeys, perm);
    if (rc != STATA_OK) {
        free(perm);
        free(sort_vars);
        stata_data_free(&g_using_cache.keys);
        stata_data_free(&g_using_cache.keepusing);
        SF_error("cmerge: Failed to sort using data\n");
        return 459;
    }

    /* Apply same permutation to keepusing
       perm[sorted_idx] = original_idx, so:
       new_data[sorted_idx] = old_data[perm[sorted_idx]] */
    if (n_keepusing > 0) {
        for (int v = 0; v < n_keepusing; v++) {
            stata_variable *var = &g_using_cache.keepusing.vars[v];
            if (var->type == STATA_TYPE_DOUBLE) {
                double *new_data = malloc(nobs * sizeof(double));
                for (size_t i = 0; i < nobs; i++) {
                    new_data[i] = var->data.dbl[perm[i]];
                }
                free(var->data.dbl);
                var->data.dbl = new_data;
            } else {
                char **new_data = malloc(nobs * sizeof(char *));
                for (size_t i = 0; i < nobs; i++) {
                    new_data[i] = var->data.str[perm[i]];
                }
                free(var->data.str);
                var->data.str = new_data;
            }
        }
    }

    free(perm);
    free(sort_vars);

    double t_sort = cmerge_get_time_ms() - t_sort_start;
    double t_total = cmerge_get_time_ms() - t_start;

    g_using_cache.loaded = 1;

    if (verbose) {
        snprintf(msg, sizeof(msg),
                 "  Loaded %zu obs (%d keys, %d keepusing) in %.1f ms, sorted in %.1f ms\n",
                 g_using_cache.nobs, nkeys, n_keepusing, t_load, t_sort);
        SF_display(msg);
    }

    SF_scal_save("_cmerge_using_nobs", (double)g_using_cache.nobs);

    return 0;
}

/* ============================================================================
 * Phase 2: Execute merge with streaming
 * ============================================================================ */

static ST_retcode cmerge_execute(const char *args)
{
    char msg[512];
    ST_retcode rc;
    double t_start, t_total_start;

    t_total_start = cmerge_get_time_ms();

    if (!g_using_cache.loaded) {
        SF_error("cmerge: Using data not loaded. Call load_using first.\n");
        return 459;
    }

    /* Parse arguments */
    char args_copy[8192];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    cmerge_type_t merge_type = MERGE_1_1;
    int nkeys = 0;
    int master_key_indices[CMERGE_MAX_KEYVARS];
    int orig_row_idx = 0;
    size_t master_nobs = 0;
    size_t master_nvars = 0;
    int n_keepusing = 0;
    int keepusing_placeholder_indices[CMERGE_MAX_VARS];
    int merge_var_idx = 0;
    int preserve_order = 0;
    int verbose = 0;

    int arg_idx = 0;
    int in_keepusing_placeholders = 0;
    int keepusing_idx = 0;

    char *token = strtok(args_copy, " ");
    while (token != NULL) {
        if (arg_idx == 0) {
            merge_type = (cmerge_type_t)atoi(token);
        }
        else if (arg_idx == 1) {
            nkeys = atoi(token);
            if (nkeys > CMERGE_MAX_KEYVARS) nkeys = CMERGE_MAX_KEYVARS;
        }
        else if (arg_idx >= 2 && arg_idx < 2 + nkeys) {
            master_key_indices[arg_idx - 2] = atoi(token);
        }
        else if (strcmp(token, "orig_row_idx") == 0) {
            token = strtok(NULL, " ");
            if (token) orig_row_idx = atoi(token);
        }
        else if (strcmp(token, "master_nobs") == 0) {
            token = strtok(NULL, " ");
            if (token) master_nobs = (size_t)atol(token);
        }
        else if (strcmp(token, "master_nvars") == 0) {
            token = strtok(NULL, " ");
            if (token) master_nvars = (size_t)atol(token);
        }
        else if (strcmp(token, "n_keepusing") == 0) {
            token = strtok(NULL, " ");
            if (token) n_keepusing = atoi(token);
        }
        else if (strcmp(token, "keepusing_placeholders") == 0) {
            in_keepusing_placeholders = 1;
            keepusing_idx = 0;
        }
        else if (in_keepusing_placeholders && keepusing_idx < n_keepusing) {
            keepusing_placeholder_indices[keepusing_idx++] = atoi(token);
        }
        else if (strcmp(token, "merge_var_idx") == 0) {
            token = strtok(NULL, " ");
            if (token) merge_var_idx = atoi(token);
            in_keepusing_placeholders = 0;  /* End of keepusing list */
        }
        else if (strcmp(token, "preserve_order") == 0) {
            token = strtok(NULL, " ");
            if (token) preserve_order = atoi(token);
        }
        else if (strcmp(token, "verbose") == 0) {
            verbose = 1;
        }
        arg_idx++;
        token = strtok(NULL, " ");
    }

    /* ===================================================================
     * Step 1: Load ONLY master keys + _orig_row
     * =================================================================== */

    if (verbose) SF_display("cmerge: Loading master keys only...\n");
    t_start = cmerge_get_time_ms();

    int n_to_load = nkeys + 1;  /* keys + _orig_row */
    int *load_indices = malloc(n_to_load * sizeof(int));
    for (int i = 0; i < nkeys; i++) {
        load_indices[i] = master_key_indices[i];
    }
    load_indices[nkeys] = orig_row_idx;

    stata_data master_minimal;
    rc = ctools_data_load_selective(&master_minimal, load_indices, n_to_load, 1, master_nobs);
    free(load_indices);

    if (rc != STATA_OK) {
        SF_error("cmerge: Failed to load master keys\n");
        return rc;
    }

    double t_load = cmerge_get_time_ms() - t_start;
    if (verbose) {
        snprintf(msg, sizeof(msg), "  Loaded %zu obs, %d vars in %.1f ms\n",
                 master_minimal.nobs, n_to_load, t_load);
        SF_display(msg);
    }

    /* ===================================================================
     * Step 2: Sort master keys
     * =================================================================== */

    if (verbose) SF_display("cmerge: Sorting master keys...\n");
    t_start = cmerge_get_time_ms();

    int *sort_vars = malloc(nkeys * sizeof(int));
    for (int i = 0; i < nkeys; i++) {
        sort_vars[i] = i + 1;  /* 1-based within master_minimal */
    }

    rc = ctools_sort_radix_lsd(&master_minimal, sort_vars, nkeys);
    free(sort_vars);

    if (rc != STATA_OK) {
        stata_data_free(&master_minimal);
        return rc;
    }

    double t_sort = cmerge_get_time_ms() - t_start;
    if (verbose) {
        snprintf(msg, sizeof(msg), "  Sorted in %.1f ms\n", t_sort);
        SF_display(msg);
    }

    /* ===================================================================
     * Step 3: Perform sorted merge join
     * =================================================================== */

    if (verbose) SF_display("cmerge: Performing merge join...\n");
    t_start = cmerge_get_time_ms();

    cmerge_output_spec_t *output_specs = NULL;
    int64_t output_nobs_signed = cmerge_sorted_join(
        &master_minimal, &g_using_cache.keys,
        nkeys, merge_type, &output_specs);

    if (output_nobs_signed < 0) {
        stata_data_free(&master_minimal);
        SF_error("cmerge: Merge join failed\n");
        return 920;
    }

    size_t output_nobs = (size_t)output_nobs_signed;

    double t_merge = cmerge_get_time_ms() - t_start;
    if (verbose) {
        snprintf(msg, sizeof(msg), "  Merge produced %zu rows in %.1f ms\n",
                 output_nobs, t_merge);
        SF_display(msg);
    }

    /* ===================================================================
     * Step 4: Build master_orig_row mapping
     * =================================================================== */

    int64_t *master_orig_rows = malloc(output_nobs * sizeof(int64_t));
    if (!master_orig_rows) {
        free(output_specs);
        stata_data_free(&master_minimal);
        return 920;
    }

    /* _orig_row is the last variable in master_minimal */
    double *orig_row_data = master_minimal.vars[nkeys].data.dbl;

    for (size_t i = 0; i < output_nobs; i++) {
        if (output_specs[i].master_sorted_row >= 0) {
            size_t sorted_idx = (size_t)output_specs[i].master_sorted_row;
            master_orig_rows[i] = (int64_t)orig_row_data[sorted_idx] - 1;  /* 0-based */
        } else {
            master_orig_rows[i] = -1;  /* Using-only */
        }
    }

    /* ===================================================================
     * Step 5: Apply preserve_order if requested
     * =================================================================== */

    if (preserve_order && output_nobs > 0) {
        if (verbose) SF_display("cmerge: Restoring original order...\n");

        /* Compute order keys */
        size_t *order_keys = malloc(output_nobs * sizeof(size_t));
        size_t *indices = malloc(output_nobs * sizeof(size_t));
        for (size_t i = 0; i < output_nobs; i++) {
            indices[i] = i;
            if (master_orig_rows[i] >= 0) {
                order_keys[i] = (size_t)master_orig_rows[i];
            } else {
                order_keys[i] = master_nobs + (size_t)output_specs[i].using_sorted_row;
            }
        }

        /* Simple insertion sort for order restoration (stable) */
        for (size_t i = 1; i < output_nobs; i++) {
            size_t key = order_keys[i];
            size_t idx = indices[i];
            cmerge_output_spec_t spec = output_specs[i];
            int64_t orig_row = master_orig_rows[i];

            size_t j = i;
            while (j > 0 && order_keys[j - 1] > key) {
                order_keys[j] = order_keys[j - 1];
                indices[j] = indices[j - 1];
                output_specs[j] = output_specs[j - 1];
                master_orig_rows[j] = master_orig_rows[j - 1];
                j--;
            }
            order_keys[j] = key;
            indices[j] = idx;
            output_specs[j] = spec;
            master_orig_rows[j] = orig_row;
        }

        free(order_keys);
        free(indices);
    }

    /* Count merge results */
    size_t n_master_only = 0, n_using_only = 0, n_matched = 0;
    for (size_t i = 0; i < output_nobs; i++) {
        switch (output_specs[i].merge_result) {
            case MERGE_RESULT_MASTER_ONLY: n_master_only++; break;
            case MERGE_RESULT_USING_ONLY: n_using_only++; break;
            case MERGE_RESULT_BOTH: n_matched++; break;
        }
    }

    /* ===================================================================
     * Step 6: Write output to Stata
     * =================================================================== */

    if (verbose) SF_display("cmerge: Writing output...\n");
    t_start = cmerge_get_time_ms();

    /* 6a: Write _merge variable */
    if (merge_var_idx > 0) {
        for (size_t i = 0; i < output_nobs; i++) {
            SF_vstore(merge_var_idx, (ST_int)(i + 1),
                      (double)output_specs[i].merge_result);
        }
    }

    /* 6b: Write keepusing variables from cache
       For shared variables (dest_idx <= master_nvars): only write for using-only rows
       For new variables (dest_idx > master_nvars): write for all rows with using data */
    for (int kv = 0; kv < n_keepusing; kv++) {
        stata_variable *src = &g_using_cache.keepusing.vars[kv];
        ST_int dest_idx = (ST_int)keepusing_placeholder_indices[kv];
        int is_shared = (dest_idx <= (ST_int)master_nvars);

        if (src->type == STATA_TYPE_DOUBLE) {
            for (size_t i = 0; i < output_nobs; i++) {
                int64_t using_row = output_specs[i].using_sorted_row;
                int64_t master_row = output_specs[i].master_sorted_row;

                /* For shared vars: only write for using-only rows (no master data)
                   For new vars: write for all rows with using data */
                if (is_shared && master_row >= 0) {
                    /* Shared var + has master data: skip (preserve master value) */
                    continue;
                }

                double val = (using_row >= 0) ? src->data.dbl[using_row] : SV_missval;
                SF_vstore(dest_idx, (ST_int)(i + 1), val);
            }
        } else {
            for (size_t i = 0; i < output_nobs; i++) {
                int64_t using_row = output_specs[i].using_sorted_row;
                int64_t master_row = output_specs[i].master_sorted_row;

                if (is_shared && master_row >= 0) {
                    continue;
                }

                const char *val = (using_row >= 0 && src->data.str[using_row]) ?
                                   src->data.str[using_row] : "";
                SF_sstore(dest_idx, (ST_int)(i + 1), (char *)val);
            }
        }
    }

    double t_write_keepusing = cmerge_get_time_ms() - t_start;

    /* 6c: Stream master variables (parallel) */
    if (verbose) SF_display("cmerge: Streaming master variables...\n");
    t_start = cmerge_get_time_ms();

    /* Determine which variables to stream:
       - All master vars from 1 to master_nvars
       - Exclude _orig_row (which is at position orig_row_idx)
       - Exclude shared variables (handled by keepusing write) */
    size_t n_stream_vars = 0;
    stream_var_args_t *stream_args = malloc(master_nvars * sizeof(stream_var_args_t));

    for (size_t v = 1; v <= master_nvars; v++) {
        if ((int)v == orig_row_idx) continue;  /* Skip _orig_row */

        /* Check if this variable is a shared keepusing variable
           (will be handled by the keepusing write, not streaming) */
        int is_shared_keepusing = 0;
        for (int kv = 0; kv < n_keepusing; kv++) {
            if (keepusing_placeholder_indices[kv] == (int)v) {
                is_shared_keepusing = 1;
                break;
            }
        }
        if (is_shared_keepusing) continue;  /* Skip - handled by keepusing */

        stream_args[n_stream_vars].stata_var_idx = (int)v;
        stream_args[n_stream_vars].master_orig_rows = master_orig_rows;
        stream_args[n_stream_vars].output_nobs = output_nobs;
        stream_args[n_stream_vars].specs = output_specs;
        stream_args[n_stream_vars].success = 0;

        /* Check if this is a key variable */
        stream_args[n_stream_vars].is_key = 0;
        stream_args[n_stream_vars].key_idx = -1;
        for (int k = 0; k < nkeys; k++) {
            if (master_key_indices[k] == (int)v) {
                stream_args[n_stream_vars].is_key = 1;
                stream_args[n_stream_vars].key_idx = k;
                break;
            }
        }

        n_stream_vars++;
    }

    /* Execute streaming in parallel */
    if (n_stream_vars >= 2) {
        ctools_thread_pool pool;
        if (ctools_pool_init(&pool, n_stream_vars, stream_args, sizeof(stream_var_args_t)) != 0) {
            free(stream_args);
            free(output_specs);
            free(master_orig_rows);
            stata_data_free(&master_minimal);
            return 920;
        }
        ctools_pool_run(&pool, stream_master_var_thread);
        ctools_pool_free(&pool);
    } else {
        for (size_t v = 0; v < n_stream_vars; v++) {
            stream_master_var_thread(&stream_args[v]);
        }
    }

    free(stream_args);

    double t_stream = cmerge_get_time_ms() - t_start;

    if (verbose) {
        snprintf(msg, sizeof(msg),
                 "  Wrote keepusing in %.1f ms, streamed %zu master vars in %.1f ms\n",
                 t_write_keepusing, n_stream_vars, t_stream);
        SF_display(msg);
    }

    /* ===================================================================
     * Cleanup and return
     * =================================================================== */

    free(output_specs);
    free(master_orig_rows);
    stata_data_free(&master_minimal);

    /* Free using cache */
    stata_data_free(&g_using_cache.keys);
    stata_data_free(&g_using_cache.keepusing);
    g_using_cache.loaded = 0;

    /* Save results */
    SF_scal_save("_cmerge_N", (double)output_nobs);
    SF_scal_save("_cmerge_N_1", (double)n_master_only);
    SF_scal_save("_cmerge_N_2", (double)n_using_only);
    SF_scal_save("_cmerge_N_3", (double)n_matched);

    double t_total = cmerge_get_time_ms() - t_total_start;
    SF_scal_save("_cmerge_time", t_total / 1000.0);

    if (verbose) {
        snprintf(msg, sizeof(msg), "\nTotal C time: %.1f ms\n", t_total);
        SF_display(msg);
        snprintf(msg, sizeof(msg), "  Load: %.1f ms, Sort: %.1f ms, Merge: %.1f ms\n",
                 t_load, t_sort, t_merge);
        SF_display(msg);
    }

    return STATA_OK;
}

/* ============================================================================
 * Clear cached using data
 * ============================================================================ */

static ST_retcode cmerge_clear(void)
{
    if (g_using_cache.loaded) {
        stata_data_free(&g_using_cache.keys);
        stata_data_free(&g_using_cache.keepusing);
        g_using_cache.loaded = 0;
    }
    return 0;
}

/* ============================================================================
 * Plugin entry point
 * ============================================================================ */

ST_retcode cmerge_main(const char *args)
{
    if (args == NULL || strlen(args) == 0) {
        SF_error("cmerge: no arguments specified\n");
        return 198;
    }

    char args_copy[8192];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    char *phase = strtok(args_copy, " ");
    const char *rest = args + strlen(phase);
    while (*rest == ' ') rest++;

    if (strcmp(phase, "load_using") == 0) {
        return cmerge_load_using(rest);
    }
    else if (strcmp(phase, "execute") == 0) {
        return cmerge_execute(rest);
    }
    else if (strcmp(phase, "clear") == 0) {
        return cmerge_clear();
    }
    else {
        char msg[256];
        snprintf(msg, sizeof(msg), "cmerge: unknown phase '%s'\n", phase);
        SF_error(msg);
        return 198;
    }
}

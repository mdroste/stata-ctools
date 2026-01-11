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

/* Order pair for preserve_order sorting */
typedef struct {
    size_t order_key;
    size_t orig_idx;
} cmerge_order_pair_t;

/* Comparison function for qsort - makes sort stable by comparing orig_idx on ties */
static int cmerge_compare_order_pairs(const void *a, const void *b) {
    const cmerge_order_pair_t *pa = (const cmerge_order_pair_t *)a;
    const cmerge_order_pair_t *pb = (const cmerge_order_pair_t *)b;
    if (pa->order_key < pb->order_key) return -1;
    if (pa->order_key > pb->order_key) return 1;
    /* Stable sort: preserve original order for equal keys */
    if (pa->orig_idx < pb->orig_idx) return -1;
    if (pa->orig_idx > pb->orig_idx) return 1;
    return 0;
}

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
 * Key comparison functions - with specialized numeric fast path
 * ============================================================================ */

/* Check if all keys are numeric (for fast path selection) */
static int cmerge_all_keys_numeric(stata_data *data, int nkeys)
{
    for (int k = 0; k < nkeys; k++) {
        if (data->vars[k].type != STATA_TYPE_DOUBLE) {
            return 0;
        }
    }
    return 1;
}

/* Fast path: compare all-numeric keys without type checking in loop */
static inline int cmerge_compare_keys_numeric(stata_data *data_a, size_t row_a,
                                               stata_data *data_b, size_t row_b,
                                               int nkeys)
{
    for (int k = 0; k < nkeys; k++) {
        double val_a = data_a->vars[k].data.dbl[row_a];
        double val_b = data_b->vars[k].data.dbl[row_b];

        /* Branchless missing value handling for common case (no missing) */
        int miss_a = SF_is_missing(val_a);
        int miss_b = SF_is_missing(val_b);

        if (miss_a | miss_b) {
            if (miss_a && miss_b) continue;
            return miss_a ? 1 : -1;
        }

        /* Use subtraction for branchless comparison when possible */
        if (val_a < val_b) return -1;
        if (val_a > val_b) return 1;
    }
    return 0;
}

/* Fast path: compare same-dataset numeric keys */
static inline int cmerge_compare_keys_numeric_same(stata_data *data, size_t row_a,
                                                    size_t row_b, int nkeys)
{
    for (int k = 0; k < nkeys; k++) {
        double *col = data->vars[k].data.dbl;
        double val_a = col[row_a];
        double val_b = col[row_b];

        int miss_a = SF_is_missing(val_a);
        int miss_b = SF_is_missing(val_b);

        if (miss_a | miss_b) {
            if (miss_a && miss_b) continue;
            return miss_a ? 1 : -1;
        }

        if (val_a < val_b) return -1;
        if (val_a > val_b) return 1;
    }
    return 0;
}

/* General path: handles mixed string/numeric keys */
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
 * Hybrid linear/binary search for match group boundaries (optimization)
 * ============================================================================ */

/* Threshold for switching from linear to binary search.
   For small groups, linear scan is faster due to branch prediction and cache locality. */
#define GROUP_SEARCH_LINEAR_THRESHOLD 8

/* Find the end of a matching group using hybrid linear/binary search.
   Returns the first index where the key differs from key at start_idx.
   Assumes data is sorted on keys.

   OPTIMIZATION: Uses linear scan for first few elements (cache-friendly),
   then switches to binary search if group is larger. */
static size_t cmerge_find_group_end_hybrid(stata_data *data, size_t start_idx,
                                            size_t max_idx, int nkeys,
                                            int all_numeric)
{
    if (start_idx >= max_idx) return start_idx + 1;

    /* Phase 1: Linear scan for first GROUP_SEARCH_LINEAR_THRESHOLD elements */
    size_t linear_end = start_idx + GROUP_SEARCH_LINEAR_THRESHOLD;
    if (linear_end > max_idx) linear_end = max_idx;

    size_t pos = start_idx + 1;
    while (pos <= linear_end) {
        int cmp;
        if (all_numeric) {
            cmp = cmerge_compare_keys_numeric_same(data, start_idx, pos, nkeys);
        } else {
            cmp = cmerge_compare_keys_same(data, start_idx, pos, nkeys);
        }
        if (cmp != 0) {
            return pos;  /* Found boundary within linear scan */
        }
        pos++;
    }

    /* If we reached here and pos > max_idx, entire remaining data matches */
    if (pos > max_idx) return max_idx + 1;

    /* Phase 2: Binary search for larger groups */
    size_t lo = pos;
    size_t hi = max_idx;

    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int cmp;

        if (all_numeric) {
            cmp = cmerge_compare_keys_numeric_same(data, start_idx, mid, nkeys);
        } else {
            cmp = cmerge_compare_keys_same(data, start_idx, mid, nkeys);
        }

        if (cmp == 0) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    return lo;
}

static int64_t cmerge_sorted_join(
    stata_data *master_keys, stata_data *using_keys,
    int nkeys, cmerge_type_t merge_type,
    cmerge_output_spec_t **output_specs_out)
{
    size_t m_nobs = master_keys->nobs;
    size_t u_nobs = using_keys->nobs;

    /* Determine if we can use the numeric fast path */
    int master_numeric = cmerge_all_keys_numeric(master_keys, nkeys);
    int using_numeric = cmerge_all_keys_numeric(using_keys, nkeys);
    int all_numeric = master_numeric && using_numeric;

    /* OPTIMIZATION 1: Better initial capacity estimation
       Based on merge type and dataset sizes to minimize reallocations */
    size_t capacity;
    switch (merge_type) {
        case MERGE_1_1:
            /* Max possible: all unique = m + u, typical: max(m,u) */
            capacity = m_nobs + u_nobs;
            break;
        case MERGE_M_1:
            /* Output = master rows (each gets one using match or none) */
            capacity = m_nobs + u_nobs;
            break;
        case MERGE_1_M:
            /* Output can expand: each master can match multiple using */
            capacity = m_nobs + u_nobs + u_nobs / 2;
            break;
        case MERGE_M_M:
        default:
            /* Worst case: cross product, typical: sum + some expansion */
            capacity = m_nobs + u_nobs + (m_nobs < u_nobs ? m_nobs : u_nobs);
            break;
    }
    /* Round up to power of 2 for better realloc behavior */
    size_t pow2 = 1024;
    while (pow2 < capacity) pow2 *= 2;
    capacity = pow2;

    cmerge_output_spec_t *specs = malloc(capacity * sizeof(cmerge_output_spec_t));
    if (!specs) return -1;

    size_t count = 0;
    size_t m_idx = 0;
    size_t u_idx = 0;

    /* Macro to ensure capacity - reduces code duplication */
    #define ENSURE_CAPACITY(needed) \
        do { \
            if (count + (needed) > capacity) { \
                while (capacity < count + (needed)) capacity *= 2; \
                cmerge_output_spec_t *new_specs = realloc(specs, capacity * sizeof(cmerge_output_spec_t)); \
                if (!new_specs) { free(specs); return -1; } \
                specs = new_specs; \
            } \
        } while (0)

    while (m_idx < m_nobs || u_idx < u_nobs) {
        if (m_idx >= m_nobs) {
            /* Master exhausted - bulk add remaining using rows */
            size_t remaining = u_nobs - u_idx;
            ENSURE_CAPACITY(remaining);
            for (size_t i = 0; i < remaining; i++) {
                specs[count].master_sorted_row = -1;
                specs[count].using_sorted_row = (int64_t)(u_idx + i);
                specs[count].merge_result = MERGE_RESULT_USING_ONLY;
                count++;
            }
            u_idx = u_nobs;
        }
        else if (u_idx >= u_nobs) {
            /* Using exhausted - bulk add remaining master rows */
            size_t remaining = m_nobs - m_idx;
            ENSURE_CAPACITY(remaining);
            for (size_t i = 0; i < remaining; i++) {
                specs[count].master_sorted_row = (int64_t)(m_idx + i);
                specs[count].using_sorted_row = -1;
                specs[count].merge_result = MERGE_RESULT_MASTER_ONLY;
                count++;
            }
            m_idx = m_nobs;
        }
        else {
            /* OPTIMIZATION 2: Use specialized numeric comparator when possible */
            int cmp;
            if (all_numeric) {
                cmp = cmerge_compare_keys_numeric(master_keys, m_idx, using_keys, u_idx, nkeys);
            } else {
                cmp = cmerge_compare_keys(master_keys, m_idx, using_keys, u_idx, nkeys);
            }

            if (cmp < 0) {
                /* Master only */
                ENSURE_CAPACITY(1);
                specs[count].master_sorted_row = (int64_t)m_idx;
                specs[count].using_sorted_row = -1;
                specs[count].merge_result = MERGE_RESULT_MASTER_ONLY;
                count++;
                m_idx++;
            }
            else if (cmp > 0) {
                /* Using only */
                ENSURE_CAPACITY(1);
                specs[count].master_sorted_row = -1;
                specs[count].using_sorted_row = (int64_t)u_idx;
                specs[count].merge_result = MERGE_RESULT_USING_ONLY;
                count++;
                u_idx++;
            }
            else {
                /* Keys match - find group boundaries */
                size_t m_start = m_idx;
                size_t u_start = u_idx;

                /* OPTIMIZATION 3: Use binary search for group boundaries
                   instead of linear scan when groups might be large */
                size_t m_end = cmerge_find_group_end_hybrid(master_keys, m_start,
                                                            m_nobs, nkeys, master_numeric);
                size_t u_end = cmerge_find_group_end_hybrid(using_keys, u_start,
                                                            u_nobs, nkeys, using_numeric);

                size_t m_count = m_end - m_start;
                size_t u_count = u_end - u_start;

                /* Pre-calculate output size for this group */
                size_t group_output;
                switch (merge_type) {
                    case MERGE_1_1:
                        group_output = 1;
                        break;
                    case MERGE_M_1:
                        group_output = m_count;
                        break;
                    case MERGE_1_M:
                        group_output = u_count;
                        break;
                    case MERGE_M_M:
                    default:
                        group_output = (m_count > u_count ? m_count : u_count);
                        break;
                }
                ENSURE_CAPACITY(group_output);

                /* Generate output based on merge type */
                switch (merge_type) {
                    case MERGE_1_1:
                        specs[count].master_sorted_row = (int64_t)m_start;
                        specs[count].using_sorted_row = (int64_t)u_start;
                        specs[count].merge_result = MERGE_RESULT_BOTH;
                        count++;
                        break;

                    case MERGE_M_1:
                        /* Multiple master rows map to one using row */
                        for (size_t i = 0; i < m_count; i++) {
                            specs[count].master_sorted_row = (int64_t)(m_start + i);
                            specs[count].using_sorted_row = (int64_t)u_start;
                            specs[count].merge_result = MERGE_RESULT_BOTH;
                            count++;
                        }
                        break;

                    case MERGE_1_M:
                        /* One master row maps to multiple using rows */
                        for (size_t i = 0; i < u_count; i++) {
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
                                specs[count].master_sorted_row = (int64_t)(m_start + i);
                                specs[count].using_sorted_row = (int64_t)(u_start + i);
                                specs[count].merge_result = MERGE_RESULT_BOTH;
                                count++;
                            }
                            /* Excess master rows */
                            for (size_t i = pairs; i < m_count; i++) {
                                specs[count].master_sorted_row = (int64_t)(m_start + i);
                                specs[count].using_sorted_row = -1;
                                specs[count].merge_result = MERGE_RESULT_MASTER_ONLY;
                                count++;
                            }
                            /* Excess using rows */
                            for (size_t i = pairs; i < u_count; i++) {
                                specs[count].master_sorted_row = -1;
                                specs[count].using_sorted_row = (int64_t)(u_start + i);
                                specs[count].merge_result = MERGE_RESULT_USING_ONLY;
                                count++;
                            }
                        }
                        break;
                }

                m_idx = m_end;
                u_idx = u_end;
            }
        }
    }

    #undef ENSURE_CAPACITY

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

/* Prefetch distance for software prefetching - 8 cache lines ahead */
#define STREAM_PREFETCH_DISTANCE 64

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
            int64_t row = src_rows[i];
            if (row >= 0) {
                SF_sdata(var_idx, (ST_int)(row + 1), strbuf);
                buf[i] = strdup(strbuf);
            } else if (a->is_key && g_using_cache.loaded) {
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
        size_t str_end8 = output_nobs - (output_nobs % 8);

        /* SCATTER to output positions - 8x unrolled */
        for (size_t i = 0; i < str_end8; i += 8) {
            SF_sstore(var_idx, (ST_int)(i + 1), buf[i]);
            SF_sstore(var_idx, (ST_int)(i + 2), buf[i + 1]);
            SF_sstore(var_idx, (ST_int)(i + 3), buf[i + 2]);
            SF_sstore(var_idx, (ST_int)(i + 4), buf[i + 3]);
            SF_sstore(var_idx, (ST_int)(i + 5), buf[i + 4]);
            SF_sstore(var_idx, (ST_int)(i + 6), buf[i + 5]);
            SF_sstore(var_idx, (ST_int)(i + 7), buf[i + 6]);
            SF_sstore(var_idx, (ST_int)(i + 8), buf[i + 7]);
            free(buf[i]); free(buf[i+1]); free(buf[i+2]); free(buf[i+3]);
            free(buf[i+4]); free(buf[i+5]); free(buf[i+6]); free(buf[i+7]);
        }
        for (size_t i = str_end8; i < output_nobs; i++) {
            SF_sstore(var_idx, (ST_int)(i + 1), buf[i]);
            free(buf[i]);
        }
        free(buf);

    } else {
        double *buf = malloc(output_nobs * sizeof(double));
        if (!buf) return (void *)1;

        /* GATHER from original master positions - 8x unrolled with prefetching */
        size_t i_end8 = output_nobs - (output_nobs % 8);
        for (size_t i = 0; i < i_end8; i += 8) {
            /* Prefetch future source row indices */
            if (i + STREAM_PREFETCH_DISTANCE < output_nobs) {
                __builtin_prefetch(&src_rows[i + STREAM_PREFETCH_DISTANCE], 0, 0);
            }

            /* Read 8 source row indices at once for better cache usage */
            int64_t r0 = src_rows[i], r1 = src_rows[i+1];
            int64_t r2 = src_rows[i+2], r3 = src_rows[i+3];
            int64_t r4 = src_rows[i+4], r5 = src_rows[i+5];
            int64_t r6 = src_rows[i+6], r7 = src_rows[i+7];

            /* Handle each row - branch-reduced pattern */
            if (r0 >= 0) SF_vdata(var_idx, (ST_int)(r0 + 1), &buf[i]);
            else if (a->is_key && g_using_cache.loaded && a->specs[i].using_sorted_row >= 0 && a->key_idx >= 0)
                buf[i] = g_using_cache.keys.vars[a->key_idx].data.dbl[a->specs[i].using_sorted_row];
            else buf[i] = SV_missval;

            if (r1 >= 0) SF_vdata(var_idx, (ST_int)(r1 + 1), &buf[i+1]);
            else if (a->is_key && g_using_cache.loaded && a->specs[i+1].using_sorted_row >= 0 && a->key_idx >= 0)
                buf[i+1] = g_using_cache.keys.vars[a->key_idx].data.dbl[a->specs[i+1].using_sorted_row];
            else buf[i+1] = SV_missval;

            if (r2 >= 0) SF_vdata(var_idx, (ST_int)(r2 + 1), &buf[i+2]);
            else if (a->is_key && g_using_cache.loaded && a->specs[i+2].using_sorted_row >= 0 && a->key_idx >= 0)
                buf[i+2] = g_using_cache.keys.vars[a->key_idx].data.dbl[a->specs[i+2].using_sorted_row];
            else buf[i+2] = SV_missval;

            if (r3 >= 0) SF_vdata(var_idx, (ST_int)(r3 + 1), &buf[i+3]);
            else if (a->is_key && g_using_cache.loaded && a->specs[i+3].using_sorted_row >= 0 && a->key_idx >= 0)
                buf[i+3] = g_using_cache.keys.vars[a->key_idx].data.dbl[a->specs[i+3].using_sorted_row];
            else buf[i+3] = SV_missval;

            if (r4 >= 0) SF_vdata(var_idx, (ST_int)(r4 + 1), &buf[i+4]);
            else if (a->is_key && g_using_cache.loaded && a->specs[i+4].using_sorted_row >= 0 && a->key_idx >= 0)
                buf[i+4] = g_using_cache.keys.vars[a->key_idx].data.dbl[a->specs[i+4].using_sorted_row];
            else buf[i+4] = SV_missval;

            if (r5 >= 0) SF_vdata(var_idx, (ST_int)(r5 + 1), &buf[i+5]);
            else if (a->is_key && g_using_cache.loaded && a->specs[i+5].using_sorted_row >= 0 && a->key_idx >= 0)
                buf[i+5] = g_using_cache.keys.vars[a->key_idx].data.dbl[a->specs[i+5].using_sorted_row];
            else buf[i+5] = SV_missval;

            if (r6 >= 0) SF_vdata(var_idx, (ST_int)(r6 + 1), &buf[i+6]);
            else if (a->is_key && g_using_cache.loaded && a->specs[i+6].using_sorted_row >= 0 && a->key_idx >= 0)
                buf[i+6] = g_using_cache.keys.vars[a->key_idx].data.dbl[a->specs[i+6].using_sorted_row];
            else buf[i+6] = SV_missval;

            if (r7 >= 0) SF_vdata(var_idx, (ST_int)(r7 + 1), &buf[i+7]);
            else if (a->is_key && g_using_cache.loaded && a->specs[i+7].using_sorted_row >= 0 && a->key_idx >= 0)
                buf[i+7] = g_using_cache.keys.vars[a->key_idx].data.dbl[a->specs[i+7].using_sorted_row];
            else buf[i+7] = SV_missval;
        }
        /* Handle remainder */
        for (size_t i = i_end8; i < output_nobs; i++) {
            if (src_rows[i] >= 0) {
                SF_vdata(var_idx, (ST_int)(src_rows[i] + 1), &buf[i]);
            } else if (a->is_key && g_using_cache.loaded) {
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

        /* SCATTER with 16-way unrolling for better cache utilization */
        size_t i_end16 = output_nobs - (output_nobs % 16);
        for (size_t i = 0; i < i_end16; i += 16) {
            SF_vstore(var_idx, (ST_int)(i + 1), buf[i]);
            SF_vstore(var_idx, (ST_int)(i + 2), buf[i + 1]);
            SF_vstore(var_idx, (ST_int)(i + 3), buf[i + 2]);
            SF_vstore(var_idx, (ST_int)(i + 4), buf[i + 3]);
            SF_vstore(var_idx, (ST_int)(i + 5), buf[i + 4]);
            SF_vstore(var_idx, (ST_int)(i + 6), buf[i + 5]);
            SF_vstore(var_idx, (ST_int)(i + 7), buf[i + 6]);
            SF_vstore(var_idx, (ST_int)(i + 8), buf[i + 7]);
            SF_vstore(var_idx, (ST_int)(i + 9), buf[i + 8]);
            SF_vstore(var_idx, (ST_int)(i + 10), buf[i + 9]);
            SF_vstore(var_idx, (ST_int)(i + 11), buf[i + 10]);
            SF_vstore(var_idx, (ST_int)(i + 12), buf[i + 11]);
            SF_vstore(var_idx, (ST_int)(i + 13), buf[i + 12]);
            SF_vstore(var_idx, (ST_int)(i + 14), buf[i + 13]);
            SF_vstore(var_idx, (ST_int)(i + 15), buf[i + 14]);
            SF_vstore(var_idx, (ST_int)(i + 16), buf[i + 15]);
        }
        for (size_t i = i_end16; i < output_nobs; i++) {
            SF_vstore(var_idx, (ST_int)(i + 1), buf[i]);
        }

        free(buf);
    }

    a->success = 1;
    return NULL;
}

/* ============================================================================
 * OPTIMIZATION 5: Parallel thread for keepusing variable writes
 * ============================================================================ */

typedef struct {
    int keepusing_idx;              /* Index in keepusing array */
    ST_int dest_idx;                /* Destination Stata variable index */
    int is_shared;                  /* Is this a shared variable? */
    size_t master_nvars;            /* For determining shared status */
    cmerge_output_spec_t *specs;    /* Output specifications */
    size_t output_nobs;             /* Number of output observations */
    int success;                    /* Thread result */
} keepusing_write_args_t;

static void *write_keepusing_var_thread(void *arg)
{
    keepusing_write_args_t *a = (keepusing_write_args_t *)arg;
    size_t output_nobs = a->output_nobs;
    ST_int dest_idx = a->dest_idx;
    int is_shared = a->is_shared;
    cmerge_output_spec_t *specs = a->specs;
    stata_variable *src = &g_using_cache.keepusing.vars[a->keepusing_idx];

    a->success = 0;

    if (src->type == STATA_TYPE_DOUBLE) {
        /* Numeric keepusing variable */
        for (size_t i = 0; i < output_nobs; i++) {
            int64_t using_row = specs[i].using_sorted_row;
            int64_t master_row = specs[i].master_sorted_row;

            /* For shared vars: only write for using-only rows */
            if (is_shared && master_row >= 0) {
                continue;
            }

            double val = (using_row >= 0) ? src->data.dbl[using_row] : SV_missval;
            SF_vstore(dest_idx, (ST_int)(i + 1), val);
        }
    } else {
        /* String keepusing variable */
        for (size_t i = 0; i < output_nobs; i++) {
            int64_t using_row = specs[i].using_sorted_row;
            int64_t master_row = specs[i].master_sorted_row;

            if (is_shared && master_row >= 0) {
                continue;
            }

            const char *val = (using_row >= 0 && src->data.str[using_row]) ?
                               src->data.str[using_row] : "";
            SF_sstore(dest_idx, (ST_int)(i + 1), (char *)val);
        }
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
     * Uses qsort O(n log n) - stable and well-tested
     * =================================================================== */

    if (preserve_order && output_nobs > 0) {
        if (verbose) SF_display("cmerge: Restoring original order...\n");
        double t_reorder = cmerge_get_time_ms();

        /* Create (order_key, original_index) pairs for sorting */
        cmerge_order_pair_t *pairs = malloc(output_nobs * sizeof(cmerge_order_pair_t));
        if (!pairs) {
            free(output_specs);
            free(master_orig_rows);
            stata_data_free(&master_minimal);
            return 920;
        }

        for (size_t i = 0; i < output_nobs; i++) {
            pairs[i].orig_idx = i;
            if (master_orig_rows[i] >= 0) {
                pairs[i].order_key = (size_t)master_orig_rows[i];
            } else {
                pairs[i].order_key = master_nobs + (size_t)output_specs[i].using_sorted_row;
            }
        }

        qsort(pairs, output_nobs, sizeof(cmerge_order_pair_t), cmerge_compare_order_pairs);

        /* Apply permutation to output_specs and master_orig_rows */
        cmerge_output_spec_t *new_specs = malloc(output_nobs * sizeof(cmerge_output_spec_t));
        int64_t *new_orig_rows = malloc(output_nobs * sizeof(int64_t));
        if (!new_specs || !new_orig_rows) {
            free(pairs);
            if (new_specs) free(new_specs);
            if (new_orig_rows) free(new_orig_rows);
            free(output_specs);
            free(master_orig_rows);
            stata_data_free(&master_minimal);
            return 920;
        }

        for (size_t i = 0; i < output_nobs; i++) {
            size_t src_idx = pairs[i].orig_idx;
            new_specs[i] = output_specs[src_idx];
            new_orig_rows[i] = master_orig_rows[src_idx];
        }

        /* Swap in the new sorted arrays */
        free(output_specs);
        free(master_orig_rows);
        output_specs = new_specs;
        master_orig_rows = new_orig_rows;

        free(pairs);

        if (verbose) {
            snprintf(msg, sizeof(msg), "  Reordered in %.1f ms\n",
                     cmerge_get_time_ms() - t_reorder);
            SF_display(msg);
        }
    }

    /* Count merge results - unrolled for better performance */
    size_t n_master_only = 0, n_using_only = 0, n_matched = 0;
    size_t cnt_end8 = output_nobs - (output_nobs % 8);

    /* Main loop - 8x unrolling */
    for (size_t i = 0; i < cnt_end8; i += 8) {
        /* Branchless counting using comparison */
        int8_t m0 = output_specs[i].merge_result;
        int8_t m1 = output_specs[i + 1].merge_result;
        int8_t m2 = output_specs[i + 2].merge_result;
        int8_t m3 = output_specs[i + 3].merge_result;
        int8_t m4 = output_specs[i + 4].merge_result;
        int8_t m5 = output_specs[i + 5].merge_result;
        int8_t m6 = output_specs[i + 6].merge_result;
        int8_t m7 = output_specs[i + 7].merge_result;

        n_master_only += (m0 == 1) + (m1 == 1) + (m2 == 1) + (m3 == 1) +
                         (m4 == 1) + (m5 == 1) + (m6 == 1) + (m7 == 1);
        n_using_only += (m0 == 2) + (m1 == 2) + (m2 == 2) + (m3 == 2) +
                        (m4 == 2) + (m5 == 2) + (m6 == 2) + (m7 == 2);
        n_matched += (m0 == 3) + (m1 == 3) + (m2 == 3) + (m3 == 3) +
                     (m4 == 3) + (m5 == 3) + (m6 == 3) + (m7 == 3);
    }
    /* Handle remainder */
    for (size_t i = cnt_end8; i < output_nobs; i++) {
        int8_t m = output_specs[i].merge_result;
        n_master_only += (m == 1);
        n_using_only += (m == 2);
        n_matched += (m == 3);
    }

    /* ===================================================================
     * Step 6: Write output to Stata
     * =================================================================== */

    if (verbose) SF_display("cmerge: Writing output...\n");
    t_start = cmerge_get_time_ms();

    /* 6a: Write _merge variable - 16x unrolled for better cache performance */
    if (merge_var_idx > 0) {
        ST_int mvar = (ST_int)merge_var_idx;
        size_t merge_end16 = output_nobs - (output_nobs % 16);

        /* Main loop - 16x unrolling (128 bytes = 2 cache lines of merge results) */
        for (size_t i = 0; i < merge_end16; i += 16) {
            SF_vstore(mvar, (ST_int)(i + 1),  (double)output_specs[i].merge_result);
            SF_vstore(mvar, (ST_int)(i + 2),  (double)output_specs[i + 1].merge_result);
            SF_vstore(mvar, (ST_int)(i + 3),  (double)output_specs[i + 2].merge_result);
            SF_vstore(mvar, (ST_int)(i + 4),  (double)output_specs[i + 3].merge_result);
            SF_vstore(mvar, (ST_int)(i + 5),  (double)output_specs[i + 4].merge_result);
            SF_vstore(mvar, (ST_int)(i + 6),  (double)output_specs[i + 5].merge_result);
            SF_vstore(mvar, (ST_int)(i + 7),  (double)output_specs[i + 6].merge_result);
            SF_vstore(mvar, (ST_int)(i + 8),  (double)output_specs[i + 7].merge_result);
            SF_vstore(mvar, (ST_int)(i + 9),  (double)output_specs[i + 8].merge_result);
            SF_vstore(mvar, (ST_int)(i + 10), (double)output_specs[i + 9].merge_result);
            SF_vstore(mvar, (ST_int)(i + 11), (double)output_specs[i + 10].merge_result);
            SF_vstore(mvar, (ST_int)(i + 12), (double)output_specs[i + 11].merge_result);
            SF_vstore(mvar, (ST_int)(i + 13), (double)output_specs[i + 12].merge_result);
            SF_vstore(mvar, (ST_int)(i + 14), (double)output_specs[i + 13].merge_result);
            SF_vstore(mvar, (ST_int)(i + 15), (double)output_specs[i + 14].merge_result);
            SF_vstore(mvar, (ST_int)(i + 16), (double)output_specs[i + 15].merge_result);
        }
        /* Handle remainder */
        for (size_t i = merge_end16; i < output_nobs; i++) {
            SF_vstore(mvar, (ST_int)(i + 1), (double)output_specs[i].merge_result);
        }
    }

    /* 6b: Write keepusing variables from cache (PARALLELIZED)
       For shared variables (dest_idx <= master_nvars): only write for using-only rows
       For new variables (dest_idx > master_nvars): write for all rows with using data */
    if (n_keepusing > 0) {
        keepusing_write_args_t *ku_args = malloc(n_keepusing * sizeof(keepusing_write_args_t));
        if (!ku_args) {
            free(output_specs);
            free(master_orig_rows);
            return 920;
        }

        /* Set up thread arguments */
        for (int kv = 0; kv < n_keepusing; kv++) {
            ku_args[kv].keepusing_idx = kv;
            ku_args[kv].dest_idx = (ST_int)keepusing_placeholder_indices[kv];
            ku_args[kv].is_shared = (ku_args[kv].dest_idx <= (ST_int)master_nvars);
            ku_args[kv].master_nvars = master_nvars;
            ku_args[kv].specs = output_specs;
            ku_args[kv].output_nobs = output_nobs;
            ku_args[kv].success = 0;
        }

        /* Execute in parallel if multiple keepusing vars, otherwise sequential */
        if (n_keepusing >= 2) {
            ctools_thread_pool ku_pool;
            if (ctools_pool_init(&ku_pool, n_keepusing, ku_args, sizeof(keepusing_write_args_t)) != 0) {
                /* Fall back to sequential */
                for (int kv = 0; kv < n_keepusing; kv++) {
                    write_keepusing_var_thread(&ku_args[kv]);
                }
            } else {
                ctools_pool_launch(&ku_pool, write_keepusing_var_thread);
                ctools_pool_join(&ku_pool);
                ctools_pool_free(&ku_pool);
            }
        } else {
            /* Single keepusing var - run directly */
            write_keepusing_var_thread(&ku_args[0]);
        }

        free(ku_args);
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

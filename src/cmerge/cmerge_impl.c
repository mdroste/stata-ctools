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
#include "ctools_spi.h"
#include "ctools_config.h"
#include "cmerge_impl.h"

/* ============================================================================
 * Memory Allocation Utilities
 * ============================================================================ */

/* Allocate memory aligned to cache line boundary for optimal memory access */
static inline void *cmerge_aligned_alloc(size_t size)
{
    void *ptr = NULL;
#if defined(_WIN32)
    ptr = _aligned_malloc(size, CACHE_LINE_SIZE);
#else
    if (posix_memalign(&ptr, CACHE_LINE_SIZE, size) != 0) {
        return NULL;
    }
#endif
    return ptr;
}

/* Free cache-line aligned memory */
static inline void cmerge_aligned_free(void *ptr)
{
#if defined(_WIN32)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

/* ============================================================================
 * String Arena Allocator for Merge Operations
 *
 * Allocates strings from a contiguous block to reduce per-string malloc
 * overhead and enable O(1) bulk free instead of O(n) individual frees.
 * ============================================================================ */

typedef struct {
    char *base;         /* Base pointer to arena memory */
    size_t capacity;    /* Total arena capacity in bytes */
    size_t used;        /* Currently used bytes */
} cmerge_string_arena;

/* Create a string arena with given capacity */
static cmerge_string_arena *cmerge_arena_create(size_t capacity)
{
    cmerge_string_arena *arena = (cmerge_string_arena *)malloc(sizeof(cmerge_string_arena));
    if (arena == NULL) return NULL;

    arena->base = (char *)malloc(capacity);
    if (arena->base == NULL) {
        free(arena);
        return NULL;
    }

    arena->capacity = capacity;
    arena->used = 0;
    return arena;
}

/* Allocate a string copy from the arena. Falls back to strdup if full. */
static char *cmerge_arena_strdup(cmerge_string_arena *arena, const char *s)
{
    size_t len = strlen(s) + 1;

    if (arena != NULL && arena->used + len <= arena->capacity) {
        char *ptr = arena->base + arena->used;
        memcpy(ptr, s, len);
        arena->used += len;
        return ptr;
    }

    /* Fallback to regular strdup */
    return strdup(s);
}

/* Check if a pointer was allocated from the arena */
static inline int cmerge_arena_owns(cmerge_string_arena *arena, const char *ptr)
{
    if (arena == NULL || ptr == NULL) return 0;
    return (ptr >= arena->base && ptr < arena->base + arena->capacity);
}

/* Free the arena. Note: strings from fallback strdup must be freed separately */
static void cmerge_arena_free(cmerge_string_arena *arena)
{
    if (arena != NULL) {
        free(arena->base);
        free(arena);
    }
}

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

/*
 * Parallel LSD Radix Sort for order pairs.
 *
 * O(n) stable sort optimized for bounded integer keys with parallel execution.
 * Uses 4 passes of 8-bit radix sort (sufficient for Stata's 2^31 obs limit).
 * Parallel histogram + parallel scatter for high throughput.
 * Falls back to qsort for very small arrays where parallel overhead dominates.
 */

/* Thread arguments for parallel histogram */
typedef struct {
    cmerge_order_pair_t *pairs;
    size_t start;
    size_t end;
    int shift;
    size_t *local_counts;  /* [256] */
} cmerge_hist_args_t;

/* Thread arguments for parallel scatter */
typedef struct {
    cmerge_order_pair_t *src;
    cmerge_order_pair_t *dst;
    size_t start;
    size_t end;
    int shift;
    size_t *offsets;  /* Starting offset for each bucket for this thread */
} cmerge_scatter_args_t;

static void *cmerge_histogram_thread(void *arg)
{
    cmerge_hist_args_t *args = (cmerge_hist_args_t *)arg;
    cmerge_order_pair_t *pairs = args->pairs;
    int shift = args->shift;
    size_t *counts = args->local_counts;

    memset(counts, 0, 256 * sizeof(size_t));

    for (size_t i = args->start; i < args->end; i++) {
        uint8_t digit = (pairs[i].order_key >> shift) & 0xFF;
        counts[digit]++;
    }

    return NULL;
}

static void *cmerge_scatter_thread(void *arg)
{
    cmerge_scatter_args_t *args = (cmerge_scatter_args_t *)arg;
    cmerge_order_pair_t *src = args->src;
    cmerge_order_pair_t *dst = args->dst;
    int shift = args->shift;
    size_t *offsets = args->offsets;

    for (size_t i = args->start; i < args->end; i++) {
        uint8_t digit = (src[i].order_key >> shift) & 0xFF;
        dst[offsets[digit]++] = src[i];
    }

    return NULL;
}

static void cmerge_radix_sort_order_pairs(cmerge_order_pair_t *pairs, size_t n)
{
    /* For small arrays, qsort is faster due to parallel overhead */
    if (n < 10000) {
        qsort(pairs, n, sizeof(cmerge_order_pair_t), cmerge_compare_order_pairs);
        return;
    }

    /* Determine thread count */
    int num_threads = NUM_THREADS;
    if (num_threads > CTOOLS_IO_MAX_THREADS) num_threads = CTOOLS_IO_MAX_THREADS;
    if ((size_t)num_threads > n / 1000) num_threads = (int)(n / 1000);
    if (num_threads < 1) num_threads = 1;

    /* Allocate auxiliary buffer and thread resources */
    cmerge_order_pair_t *aux = (cmerge_order_pair_t *)malloc(n * sizeof(cmerge_order_pair_t));
    pthread_t *threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
    cmerge_hist_args_t *hist_args = (cmerge_hist_args_t *)malloc(num_threads * sizeof(cmerge_hist_args_t));
    cmerge_scatter_args_t *scatter_args = (cmerge_scatter_args_t *)malloc(num_threads * sizeof(cmerge_scatter_args_t));
    size_t *all_counts = (size_t *)malloc(num_threads * 256 * sizeof(size_t));
    size_t *all_offsets = (size_t *)malloc(num_threads * 256 * sizeof(size_t));

    if (!aux || !threads || !hist_args || !scatter_args || !all_counts || !all_offsets) {
        /* Fallback to qsort on allocation failure */
        free(aux);
        free(threads);
        free(hist_args);
        free(scatter_args);
        free(all_counts);
        free(all_offsets);
        qsort(pairs, n, sizeof(cmerge_order_pair_t), cmerge_compare_order_pairs);
        return;
    }

    /* Calculate chunk boundaries */
    size_t chunk_size = n / num_threads;
    size_t remainder = n % num_threads;

    cmerge_order_pair_t *src = pairs;
    cmerge_order_pair_t *dst = aux;

    /* 4 passes of 8-bit radix sort on order_key (32 bits sufficient for Stata) */
    for (int pass = 0; pass < 4; pass++) {
        int shift = pass * 8;

        /* Phase 1: Parallel histogram */
        size_t offset = 0;
        for (int t = 0; t < num_threads; t++) {
            hist_args[t].pairs = src;
            hist_args[t].start = offset;
            hist_args[t].end = offset + chunk_size + (t < (int)remainder ? 1 : 0);
            hist_args[t].shift = shift;
            hist_args[t].local_counts = all_counts + t * 256;
            offset = hist_args[t].end;

            pthread_create(&threads[t], NULL, cmerge_histogram_thread, &hist_args[t]);
        }

        for (int t = 0; t < num_threads; t++) {
            pthread_join(threads[t], NULL);
        }

        /* Compute global counts and per-thread offsets */
        size_t global_counts[256];
        memset(global_counts, 0, sizeof(global_counts));

        for (int t = 0; t < num_threads; t++) {
            size_t *local = all_counts + t * 256;
            for (int b = 0; b < 256; b++) {
                global_counts[b] += local[b];
            }
        }

        /* Convert global counts to prefix sums */
        size_t total = 0;
        size_t global_offsets[256];
        for (int b = 0; b < 256; b++) {
            global_offsets[b] = total;
            total += global_counts[b];
        }

        /* Compute per-thread starting offsets for each bucket */
        for (int b = 0; b < 256; b++) {
            size_t bucket_offset = global_offsets[b];
            for (int t = 0; t < num_threads; t++) {
                all_offsets[t * 256 + b] = bucket_offset;
                bucket_offset += all_counts[t * 256 + b];
            }
        }

        /* Phase 2: Parallel scatter */
        offset = 0;
        for (int t = 0; t < num_threads; t++) {
            scatter_args[t].src = src;
            scatter_args[t].dst = dst;
            scatter_args[t].start = offset;
            scatter_args[t].end = offset + chunk_size + (t < (int)remainder ? 1 : 0);
            scatter_args[t].shift = shift;
            scatter_args[t].offsets = all_offsets + t * 256;
            offset = scatter_args[t].end;

            pthread_create(&threads[t], NULL, cmerge_scatter_thread, &scatter_args[t]);
        }

        for (int t = 0; t < num_threads; t++) {
            pthread_join(threads[t], NULL);
        }

        /* Swap buffers */
        cmerge_order_pair_t *tmp = src;
        src = dst;
        dst = tmp;
    }

    /* If result is in aux, copy back to pairs */
    if (src != pairs) {
        memcpy(pairs, src, n * sizeof(cmerge_order_pair_t));
    }

    free(aux);
    free(threads);
    free(hist_args);
    free(scatter_args);
    free(all_counts);
    free(all_offsets);
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
    int merge_by_n;             /* 1 if merging by observation number (_n) */
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

   @param data       Sorted data to search
   @param start_idx  Starting index of the group
   @param nobs       Total number of observations (exclusive upper bound)
   @param nkeys      Number of key variables
   @param all_numeric True if all keys are numeric (for fast path)

   @return Index of first non-matching element, or nobs if all remaining match

   OPTIMIZATION: Uses linear scan for first few elements (cache-friendly),
   then switches to binary search if group is larger. */
static size_t cmerge_find_group_end_hybrid(stata_data *data, size_t start_idx,
                                            size_t nobs, int nkeys,
                                            int all_numeric)
{
    /* If start_idx is at or past the last valid index, return nobs */
    if (start_idx + 1 >= nobs) return nobs;

    /* Phase 1: Linear scan for first GROUP_SEARCH_LINEAR_THRESHOLD elements
       Note: nobs is exclusive upper bound, so valid indices are 0..nobs-1 */
    size_t linear_end = start_idx + 1 + GROUP_SEARCH_LINEAR_THRESHOLD;
    if (linear_end > nobs) linear_end = nobs;

    size_t pos = start_idx + 1;
    while (pos < linear_end) {  /* Use < not <= since nobs is exclusive */
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

    /* If we've scanned all remaining data, return nobs */
    if (pos >= nobs) return nobs;

    /* Phase 2: Binary search for larger groups
       Binary search for the first index where key differs from start_idx's key.
       lo = first unchecked index, hi = nobs (exclusive upper bound) */
    size_t lo = pos;
    size_t hi = nobs;

    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int cmp;

        if (all_numeric) {
            cmp = cmerge_compare_keys_numeric_same(data, start_idx, mid, nkeys);
        } else {
            cmp = cmerge_compare_keys_same(data, start_idx, mid, nkeys);
        }

        if (cmp == 0) {
            lo = mid + 1;  /* Key matches, search right half */
        } else {
            hi = mid;      /* Key differs, search left half */
        }
    }

    return lo;  /* lo is now the first index where key differs, or nobs if all match */
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
            /* Sequential pairing with expansion: output count = max of group counts */
            capacity = m_nobs + u_nobs;
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
                        group_output = (m_count > u_count) ? m_count : u_count;  /* Max for sequential pairing */
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
                        /* Sequential pairing: pair row 1 with row 1, row 2 with row 2, etc.
                         * If one side has more rows, excess rows get the last row of the other side.
                         * This matches Stata's m:m merge behavior. */
                        {
                            size_t max_count = (m_count > u_count) ? m_count : u_count;
                            for (size_t i = 0; i < max_count; i++) {
                                /* For master: use row i if available, else use last row */
                                size_t mi = (i < m_count) ? i : (m_count - 1);
                                /* For using: use row i if available, else use last row */
                                size_t ui = (i < u_count) ? i : (u_count - 1);

                                specs[count].master_sorted_row = (int64_t)(m_start + mi);
                                specs[count].using_sorted_row = (int64_t)(u_start + ui);
                                specs[count].merge_result = MERGE_RESULT_BOTH;
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

/* NOTE: This function is currently unused but retained for potential future use */
__attribute__((unused))
static void *stream_master_var_thread(void *arg)
{
    stream_var_args_t *a = (stream_var_args_t *)arg;
    size_t output_nobs = a->output_nobs;
    ST_int var_idx = (ST_int)a->stata_var_idx;
    int64_t *src_rows = a->master_orig_rows;
    int is_string = SF_var_is_string(var_idx);

    a->success = 0;

    if (is_string) {
        /* STRING VARIABLE STREAMING - Sequential I/O with C buffer
         *
         * CRITICAL: Stata's SPI (SF_sdata/SF_sstore) is NOT thread-safe.
         * OpenMP parallel for on SPI calls causes data corruption.
         * We must follow the same pattern as numerics:
         * 1. Read ALL source strings sequentially into C buffer
         * 2. Permute in C memory
         * 3. Write back to Stata sequentially
         */

        /* Find max source row to determine source buffer size */
        int64_t max_src_row = -1;
        for (size_t i = 0; i < output_nobs; i++) {
            if (src_rows[i] > max_src_row) max_src_row = src_rows[i];
        }
        size_t src_nobs = (max_src_row >= 0) ? (size_t)(max_src_row + 1) : 0;

        /* Allocate source and output string buffers */
        char **src_buf = NULL;
        char **out_buf = (char **)calloc(output_nobs, sizeof(char *));
        if (!out_buf) return (void *)1;

        if (src_nobs > 0) {
            src_buf = (char **)calloc(src_nobs, sizeof(char *));
            if (!src_buf) {
                free(out_buf);
                return (void *)1;
            }

            /* STEP 1: Sequential read from Stata into source buffer */
            char strbuf[2049];
            for (size_t i = 0; i < src_nobs; i++) {
                SF_sdata(var_idx, (ST_int)(i + 1), strbuf);
                src_buf[i] = strdup(strbuf);
                if (!src_buf[i]) {
                    /* Cleanup on allocation failure */
                    for (size_t j = 0; j < i; j++) free(src_buf[j]);
                    free(src_buf);
                    free(out_buf);
                    return (void *)1;
                }
            }
        }

        /* STEP 2: Permute in C memory */
        char **using_key_data = (a->is_key && g_using_cache.loaded && a->key_idx >= 0) ?
                                 g_using_cache.keys.vars[a->key_idx].data.str : NULL;

        for (size_t i = 0; i < output_nobs; i++) {
            int64_t row = src_rows[i];
            if (row >= 0 && src_buf && src_buf[row]) {
                out_buf[i] = strdup(src_buf[row]);  /* Copy - multiple outputs may reference same source row */
            } else if (using_key_data != NULL) {
                int64_t using_row = a->specs[i].using_sorted_row;
                if (using_row >= 0 && using_key_data[using_row]) {
                    out_buf[i] = strdup(using_key_data[using_row]);
                } else {
                    out_buf[i] = strdup("");
                }
            } else {
                out_buf[i] = strdup("");
            }
        }

        /* Free remaining src_buf entries (unused rows) */
        if (src_buf) {
            for (size_t i = 0; i < src_nobs; i++) {
                if (src_buf[i]) free(src_buf[i]);
            }
            free(src_buf);
        }

        /* STEP 3: Sequential write to Stata */
        for (size_t i = 0; i < output_nobs; i++) {
            SF_sstore(var_idx, (ST_int)(i + 1), out_buf[i] ? out_buf[i] : "");
        }

        /* Free output buffer */
        for (size_t i = 0; i < output_nobs; i++) {
            if (out_buf[i]) free(out_buf[i]);
        }
        free(out_buf);

    } else {
        /* NUMERIC VARIABLE STREAMING - Optimized Sequential I/O
         *
         * Key insight: Stata's SPI (SF_vdata/SF_vstore) is inherently sequential I/O.
         * Row-level OpenMP adds overhead without benefit. Instead, we use:
         * 1. 16x batched SPI calls to reduce function call overhead
         * 2. Double buffering for software pipelining (overlap load/store)
         * 3. Prefetching for improved memory latency hiding
         * 4. Variable-level parallelism via pthreads (handled by caller)
         *
         * This matches the efficient pattern in ctools_data_io.c.
         */

        /* Find max source row to determine source buffer size */
        int64_t max_src_row = -1;
        for (size_t i = 0; i < output_nobs; i++) {
            if (src_rows[i] > max_src_row) max_src_row = src_rows[i];
        }

        size_t src_nobs = (max_src_row >= 0) ? (size_t)(max_src_row + 1) : 0;

        /* Allocate aligned source and output buffers */
        double *src_buf = NULL;
        double *out_buf = (double *)cmerge_aligned_alloc(output_nobs * sizeof(double));
        if (!out_buf) return (void *)1;

        if (src_nobs > 0) {
            src_buf = (double *)cmerge_aligned_alloc(src_nobs * sizeof(double));
            if (!src_buf) {
                cmerge_aligned_free(out_buf);
                return (void *)1;
            }

            /* STEP 1: Sequential batched read from Stata into source buffer
             * Uses 16x unrolling with software pipelining for maximum throughput */
            double v0[16], v1[16];
            size_t i_end_16 = src_nobs - (src_nobs % 16);
            size_t prefetch_end = (src_nobs > 32) ? (i_end_16 - 32) : 0;

            size_t i = 0;
            if (prefetch_end > 0) {
                /* Prime the pipeline */
                SF_VDATA_BATCH16(var_idx, 1, v0);

                for (; i < prefetch_end; i += 16) {
                    /* Prefetch destination ahead */
                    CTOOLS_PREFETCH(&src_buf[i + 32]);
                    CTOOLS_PREFETCH(&src_buf[i + 40]);

                    /* Load next batch while storing current */
                    SF_VDATA_BATCH16(var_idx, i + 17, v1);

                    /* Store to C buffer */
                    memcpy(&src_buf[i], v0, sizeof(v0));

                    /* Swap buffers */
                    memcpy(v0, v1, sizeof(v0));
                }
            }

            /* Remaining 16-element blocks */
            for (; i < i_end_16; i += 16) {
                SF_VDATA_BATCH16(var_idx, i + 1, v0);
                memcpy(&src_buf[i], v0, sizeof(v0));
            }

            /* Handle remainder */
            for (; i < src_nobs; i++) {
                SF_vdata(var_idx, (ST_int)(i + 1), &src_buf[i]);
            }
        }

        /* STEP 2: Permute in C memory (fast - pure memory operations) */
        if (!a->is_key) {
            /* FAST PATH: Non-key variable - simple gather */
            for (size_t i = 0; i < output_nobs; i++) {
                int64_t row = src_rows[i];
                out_buf[i] = (row >= 0) ? src_buf[row] : SV_missval;
            }
        } else {
            /* SLOW PATH: Key variable with using-cache fallback */
            double *using_key_data = (g_using_cache.loaded && a->key_idx >= 0) ?
                                      g_using_cache.keys.vars[a->key_idx].data.dbl : NULL;

            for (size_t i = 0; i < output_nobs; i++) {
                int64_t row = src_rows[i];
                if (row >= 0) {
                    out_buf[i] = src_buf[row];
                } else if (using_key_data != NULL) {
                    int64_t using_row = a->specs[i].using_sorted_row;
                    out_buf[i] = (using_row >= 0) ? using_key_data[using_row] : SV_missval;
                } else {
                    out_buf[i] = SV_missval;
                }
            }
        }

        /* STEP 3: Sequential batched write to Stata
         * Uses 16x unrolling with software pipelining */
        {
            double v0[16], v1[16];
            size_t i_end_16 = output_nobs - (output_nobs % 16);
            size_t prefetch_end = (output_nobs > 32) ? (i_end_16 - 32) : 0;

            size_t i = 0;
            if (prefetch_end > 0) {
                /* Prime the pipeline */
                memcpy(v0, &out_buf[0], sizeof(v0));

                for (; i < prefetch_end; i += 16) {
                    /* Prefetch source ahead */
                    CTOOLS_PREFETCH(&out_buf[i + 32]);
                    CTOOLS_PREFETCH(&out_buf[i + 40]);

                    /* Load next batch */
                    memcpy(v1, &out_buf[i + 16], sizeof(v1));

                    /* Write current batch to Stata */
                    SF_VSTORE_BATCH16(var_idx, i + 1, v0);

                    /* Swap buffers */
                    memcpy(v0, v1, sizeof(v0));
                }
            }

            /* Remaining 16-element blocks */
            for (; i < i_end_16; i += 16) {
                memcpy(v0, &out_buf[i], sizeof(v0));
                SF_VSTORE_BATCH16(var_idx, i + 1, v0);
            }

            /* Handle remainder */
            for (; i < output_nobs; i++) {
                SF_vstore(var_idx, (ST_int)(i + 1), out_buf[i]);
            }
        }

        if (src_buf) cmerge_aligned_free(src_buf);
        cmerge_aligned_free(out_buf);
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
    int update_mode;                /* Update missing master values */
    int replace_mode;               /* Replace all master values */
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
    int update_mode = a->update_mode;
    int replace_mode = a->replace_mode;
    cmerge_output_spec_t *specs = a->specs;
    stata_variable *src = &g_using_cache.keepusing.vars[a->keepusing_idx];

    a->success = 0;

    if (src->type == STATA_TYPE_DOUBLE) {
        /* Numeric keepusing variable */
        for (size_t i = 0; i < output_nobs; i++) {
            int64_t using_row = specs[i].using_sorted_row;
            int8_t merge_result = specs[i].merge_result;

            /* Determine if we should write using value */
            int should_write = 0;

            if (merge_result == MERGE_RESULT_USING_ONLY) {
                /* Always write for using-only rows */
                should_write = 1;
            }
            else if (merge_result == MERGE_RESULT_BOTH) {
                /* Matched row - depends on shared status and update/replace */
                if (!is_shared) {
                    /* Non-shared var: always write using value for matched rows */
                    should_write = 1;
                }
                else if (replace_mode) {
                    /* Shared var with replace: always overwrite */
                    should_write = 1;
                }
                else if (update_mode) {
                    /* Shared var with update: only if master value is missing */
                    double current_val;
                    SF_vdata(dest_idx, (ST_int)(i + 1), &current_val);
                    should_write = SF_is_missing(current_val);
                }
                /* else: shared var without update/replace - don't overwrite master */
            }

            if (should_write && using_row >= 0) {
                SF_vstore(dest_idx, (ST_int)(i + 1), src->data.dbl[using_row]);
            }
        }
    } else {
        /* String keepusing variable */
        char str_buf[2049];
        for (size_t i = 0; i < output_nobs; i++) {
            int64_t using_row = specs[i].using_sorted_row;
            int8_t merge_result = specs[i].merge_result;

            /* Determine if we should write using value */
            int should_write = 0;

            if (merge_result == MERGE_RESULT_USING_ONLY) {
                /* Always write for using-only rows */
                should_write = 1;
            }
            else if (merge_result == MERGE_RESULT_BOTH) {
                /* Matched row - depends on shared status and update/replace */
                if (!is_shared) {
                    /* Non-shared var: always write using value for matched rows */
                    should_write = 1;
                }
                else if (replace_mode) {
                    /* Shared var with replace: always overwrite */
                    should_write = 1;
                }
                else if (update_mode) {
                    /* Shared var with update: only if master value is missing (empty string) */
                    SF_sdata(dest_idx, (ST_int)(i + 1), str_buf);
                    should_write = (str_buf[0] == '\0');
                }
                /* else: shared var without update/replace - don't overwrite master */
            }

            if (should_write && using_row >= 0) {
                const char *val = src->data.str[using_row] ? src->data.str[using_row] : "";
                SF_sstore(dest_idx, (ST_int)(i + 1), (char *)val);
            }
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
    int sorted = 0;
    int merge_by_n = 0;

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
        else if (strcmp(token, "sorted") == 0) {
            sorted = 1;
        }
        else if (strcmp(token, "merge_by_n") == 0) {
            merge_by_n = 1;
        }
        arg_idx++;
        token = strtok(NULL, " ");
    }

    (void)verbose;  /* Timing output moved to ado file */


    /* Free any previous cache */
    if (g_using_cache.loaded) {
        stata_data_free(&g_using_cache.keys);
        stata_data_free(&g_using_cache.keepusing);
        g_using_cache.loaded = 0;
    }

    double t_phase1_start = cmerge_get_time_ms();

    double t_load_keys = 0;
    double t_load_keepusing = 0;

    /* For merge_by_n, we don't need key variables - just keepusing */
    if (merge_by_n) {
        /* Initialize empty keys structure */
        stata_data_init(&g_using_cache.keys);
        g_using_cache.keys.nobs = SF_nobs();  /* Get nobs from Stata */

        /* Load keepusing variables */
        if (n_keepusing > 0) {
            double t_load_keepusing_start = cmerge_get_time_ms();
            rc = ctools_data_load_selective(&g_using_cache.keepusing, keepusing_indices, n_keepusing, 0, 0);
            if (rc != STATA_OK) {
                SF_error("cmerge: Failed to load keepusing variables\n");
                return 459;
            }
            if (g_using_cache.keepusing.vars == NULL) {
                SF_error("cmerge: FATAL - keepusing.vars is NULL after load\n");
                return 920;
            }
            t_load_keepusing = cmerge_get_time_ms() - t_load_keepusing_start;
            g_using_cache.keys.nobs = g_using_cache.keepusing.nobs;
        } else {
            stata_data_init(&g_using_cache.keepusing);
        }
    } else {
        /* Standard merge: Load ONLY key variables */
        double t_load_keys_start = cmerge_get_time_ms();
        rc = ctools_data_load_selective(&g_using_cache.keys, key_indices, nkeys, 0, 0);
        if (rc != STATA_OK) {
            SF_error("cmerge: Failed to load using keys\n");
            return 459;
        }
        /* Critical null check - prevents crash if load failed silently */
        if (g_using_cache.keys.vars == NULL) {
            SF_error("cmerge: FATAL - keys.vars is NULL after load\n");
            return 920;
        }
        t_load_keys = cmerge_get_time_ms() - t_load_keys_start;

        /* Load ONLY keepusing variables */
        if (n_keepusing > 0) {
            double t_load_keepusing_start = cmerge_get_time_ms();
            rc = ctools_data_load_selective(&g_using_cache.keepusing, keepusing_indices, n_keepusing, 0, 0);
            if (rc != STATA_OK) {
                stata_data_free(&g_using_cache.keys);
                SF_error("cmerge: Failed to load keepusing variables\n");
                return 459;
            }
            /* Critical null check */
            if (g_using_cache.keepusing.vars == NULL) {
                stata_data_free(&g_using_cache.keys);
                SF_error("cmerge: FATAL - keepusing.vars is NULL after load\n");
                return 920;
            }
            t_load_keepusing = cmerge_get_time_ms() - t_load_keepusing_start;
        } else {
            stata_data_init(&g_using_cache.keepusing);
        }
    }

    g_using_cache.nobs = g_using_cache.keys.nobs;
    g_using_cache.nkeys = nkeys;
    g_using_cache.n_keepusing = n_keepusing;
    g_using_cache.merge_by_n = merge_by_n;

    double t_sort = 0;
    double t_apply_perm = 0;

    /* Sort using data on keys (skip if sorted option specified or merge_by_n) */
    if (!sorted && !merge_by_n) {
        double t_sort_start = cmerge_get_time_ms();

        int *sort_vars = malloc(nkeys * sizeof(int));
        if (!sort_vars) {
            stata_data_free(&g_using_cache.keys);
            stata_data_free(&g_using_cache.keepusing);
            SF_error("cmerge: Memory allocation failed for sort_vars\n");
            return 920;
        }
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
        rc = ctools_sort_ips4o_with_perm(&g_using_cache.keys, sort_vars, nkeys, perm);
        if (rc != STATA_OK) {
            free(perm);
            free(sort_vars);
            stata_data_free(&g_using_cache.keys);
            stata_data_free(&g_using_cache.keepusing);
            SF_error("cmerge: Failed to sort using data\n");
            return 459;
        }

        t_sort = cmerge_get_time_ms() - t_sort_start;

        /* Apply same permutation to keepusing */
        double t_apply_perm_start = cmerge_get_time_ms();
        if (n_keepusing > 0 && g_using_cache.keepusing.vars != NULL) {
            for (int v = 0; v < n_keepusing; v++) {
                stata_variable *var = &g_using_cache.keepusing.vars[v];
                if (var->type == STATA_TYPE_DOUBLE) {
                    double *new_data = malloc(nobs * sizeof(double));
                    if (!new_data) {
                        free(perm);
                        free(sort_vars);
                        stata_data_free(&g_using_cache.keys);
                        stata_data_free(&g_using_cache.keepusing);
                        SF_error("cmerge: Memory allocation failed in permutation\n");
                        return 920;
                    }
                    for (size_t i = 0; i < nobs; i++) {
                        /* Bounds check on permutation index */
                        if (perm[i] >= nobs) {
                            free(new_data);
                            free(perm);
                            free(sort_vars);
                            stata_data_free(&g_using_cache.keys);
                            stata_data_free(&g_using_cache.keepusing);
                            SF_error("cmerge: Invalid permutation index\n");
                            return 459;
                        }
                        new_data[i] = var->data.dbl[perm[i]];
                    }
                    free(var->data.dbl);
                    var->data.dbl = new_data;
                } else {
                    char **new_data = malloc(nobs * sizeof(char *));
                    if (!new_data) {
                        free(perm);
                        free(sort_vars);
                        stata_data_free(&g_using_cache.keys);
                        stata_data_free(&g_using_cache.keepusing);
                        SF_error("cmerge: Memory allocation failed in permutation\n");
                        return 920;
                    }
                    for (size_t i = 0; i < nobs; i++) {
                        /* Bounds check on permutation index */
                        if (perm[i] >= nobs) {
                            free(new_data);
                            free(perm);
                            free(sort_vars);
                            stata_data_free(&g_using_cache.keys);
                            stata_data_free(&g_using_cache.keepusing);
                            SF_error("cmerge: Invalid permutation index\n");
                            return 459;
                        }
                        new_data[i] = var->data.str[perm[i]];
                    }
                    free(var->data.str);
                    var->data.str = new_data;
                }
            }
        }
        t_apply_perm = cmerge_get_time_ms() - t_apply_perm_start;

        free(perm);
        free(sort_vars);
    }

    g_using_cache.loaded = 1;

    double t_phase1_total = cmerge_get_time_ms() - t_phase1_start;

    /* Save Phase 1 timing to Stata scalars (in seconds) */
    SF_scal_save("_cmerge_using_nobs", (double)g_using_cache.nobs);
    SF_scal_save("_cmerge_p1_load_keys", t_load_keys / 1000.0);
    SF_scal_save("_cmerge_p1_load_keepusing", t_load_keepusing / 1000.0);
    SF_scal_save("_cmerge_p1_sort", t_sort / 1000.0);
    SF_scal_save("_cmerge_p1_apply_perm", t_apply_perm / 1000.0);
    SF_scal_save("_cmerge_p1_total", t_phase1_total / 1000.0);

    return 0;
}

/* ============================================================================
 * Phase 2: Execute merge with full load/permute/store
 *
 * This implementation loads ALL master variables into C memory, applies the
 * merge permutation in C, then stores everything back using parallel I/O.
 * This is more memory-intensive than streaming but enables parallel I/O.
 * ============================================================================ */

static ST_retcode cmerge_execute(const char *args)
{
    ST_retcode rc;
    double t_start, t_total_start;

    t_total_start = cmerge_get_time_ms();

    if (!g_using_cache.loaded) {
        SF_error("cmerge: Using data not loaded. Call load_using first.\n");
        return 459;
    }

    /* Critical null check - prevents crash if cache was corrupted
       Note: For merge_by_n mode, keys.vars is intentionally NULL */
    if (g_using_cache.keys.vars == NULL && !g_using_cache.merge_by_n) {
        SF_error("cmerge: FATAL - g_using_cache.keys.vars is NULL!\n");
        return 459;
    }

    /* Parse arguments - use static buffer to avoid stack issues */
    static char args_copy[8192];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    cmerge_type_t merge_type = MERGE_1_1;
    int nkeys = 0;
    static int master_key_indices[CMERGE_MAX_KEYVARS];
    int orig_row_idx = 0;
    size_t master_nobs = 0;
    size_t master_nvars = 0;
    int n_keepusing = 0;
    static int keepusing_placeholder_indices[CMERGE_MAX_VARS];
    int merge_var_idx = 0;
    int preserve_order = 0;
    int verbose = 0;
    int sorted = 0;
    int update_mode = 0;
    int replace_mode = 0;
    static int shared_flags[CMERGE_MAX_VARS];

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
            in_keepusing_placeholders = 0;
        }
        else if (strcmp(token, "preserve_order") == 0) {
            token = strtok(NULL, " ");
            if (token) preserve_order = atoi(token);
        }
        else if (strcmp(token, "verbose") == 0) {
            verbose = 1;
        }
        else if (strcmp(token, "sorted") == 0) {
            sorted = 1;
        }
        else if (strcmp(token, "update") == 0) {
            update_mode = 1;
        }
        else if (strcmp(token, "replace") == 0) {
            replace_mode = 1;
        }
        else if (strcmp(token, "shared_flags") == 0) {
            for (int i = 0; i < n_keepusing; i++) {
                token = strtok(NULL, " ");
                shared_flags[i] = token ? atoi(token) : 0;
            }
        }
        else if (strcmp(token, "merge_by_n") == 0) {
            /* Parsed but value comes from cache */
        }
        arg_idx++;
        token = strtok(NULL, " ");
    }

    /* Get merge_by_n from cache (set during load_using phase) */
    int merge_by_n = g_using_cache.merge_by_n;

    (void)verbose;  /* Timing output moved to ado file */

    /* ===================================================================
     * Step 1: Load ALL master variables (parallel I/O)
     * =================================================================== */

    t_start = cmerge_get_time_ms();

    /* Build array of all variable indices to load */
    int *all_var_indices = malloc(master_nvars * sizeof(int));
    if (!all_var_indices) {
        SF_error("cmerge: Failed to allocate var indices\n");
        return 920;
    }
    for (size_t v = 0; v < master_nvars; v++) {
        all_var_indices[v] = (int)(v + 1);
    }

    /* Load all master variables using parallel I/O */
    stata_data master_data;
    rc = ctools_data_load_selective(&master_data, all_var_indices, master_nvars, 1, master_nobs);
    if (rc != STATA_OK) {
        free(all_var_indices);
        SF_error("cmerge: Failed to load master data\n");
        return rc;
    }

    /* Critical null check - prevents crash if load failed silently */
    if (master_data.vars == NULL) {
        free(all_var_indices);
        SF_error("cmerge: FATAL - master_data.vars is NULL after load!\n");
        return 920;
    }
    double t_load = cmerge_get_time_ms() - t_start;

    /* ===================================================================
     * Step 2: Sort master data on keys (skip if sorted option specified or merge_by_n)
     * =================================================================== */

    double t_sort = 0;
    size_t *sort_perm = NULL;

    if (!sorted && !merge_by_n) {
        t_start = cmerge_get_time_ms();

        /* Build sort variable indices (within master_data, 1-based) */
        int *sort_vars = malloc(nkeys * sizeof(int));
        if (!sort_vars) {
            stata_data_free(&master_data);
            free(all_var_indices);
            SF_error("cmerge: Failed to allocate sort_vars\n");
            return 920;
        }

        for (int i = 0; i < nkeys; i++) {
            sort_vars[i] = master_key_indices[i];  /* Already 1-based Stata index */
        }

        /* Allocate permutation array to capture sort order */
        sort_perm = malloc(master_nobs * sizeof(size_t));
        if (!sort_perm) {
            free(sort_vars);
            stata_data_free(&master_data);
            free(all_var_indices);
            return 920;
        }

        rc = ctools_sort_ips4o_with_perm(&master_data, sort_vars, nkeys, sort_perm);
        free(sort_vars);

        if (rc != STATA_OK) {
            free(sort_perm);
            stata_data_free(&master_data);
            free(all_var_indices);
            return rc;
        }

        t_sort = cmerge_get_time_ms() - t_start;
    }

    /* ===================================================================
     * Step 3: Perform merge join
     * For merge_by_n: simple row-by-row matching
     * For standard merge: sorted join using keys
     * =================================================================== */

    t_start = cmerge_get_time_ms();

    cmerge_output_spec_t *output_specs = NULL;
    size_t output_nobs = 0;

    if (merge_by_n) {
        /* Merge by observation number (_n): simple row-by-row matching */
        size_t m_nobs = master_nobs;
        size_t u_nobs = g_using_cache.nobs;
        size_t max_nobs = (m_nobs > u_nobs) ? m_nobs : u_nobs;

        output_specs = malloc(max_nobs * sizeof(cmerge_output_spec_t));
        if (!output_specs) {
            if (sort_perm) free(sort_perm);
            stata_data_free(&master_data);
            free(all_var_indices);
            return 920;
        }

        /* Build output specs: row i matches row i */
        for (size_t i = 0; i < max_nobs; i++) {
            if (i < m_nobs && i < u_nobs) {
                /* Both master and using have this row - matched */
                output_specs[i].master_sorted_row = (int64_t)i;
                output_specs[i].using_sorted_row = (int64_t)i;
                output_specs[i].merge_result = MERGE_RESULT_BOTH;
            } else if (i < m_nobs) {
                /* Only master has this row */
                output_specs[i].master_sorted_row = (int64_t)i;
                output_specs[i].using_sorted_row = -1;
                output_specs[i].merge_result = MERGE_RESULT_MASTER_ONLY;
            } else {
                /* Only using has this row */
                output_specs[i].master_sorted_row = -1;
                output_specs[i].using_sorted_row = (int64_t)i;
                output_specs[i].merge_result = MERGE_RESULT_USING_ONLY;
            }
        }
        output_nobs = max_nobs;
    } else {
        /* Standard merge: create keys view and perform sorted join */
        stata_data master_keys_view;
        stata_data_init(&master_keys_view);
        master_keys_view.nobs = master_data.nobs;
        master_keys_view.nvars = nkeys;
        master_keys_view.vars = malloc(nkeys * sizeof(stata_variable));
        if (!master_keys_view.vars) {
            if (sort_perm) free(sort_perm);
            stata_data_free(&master_data);
            free(all_var_indices);
            return 920;
        }

        /* Point to the key variables in master_data (0-indexed in vars array) */
        for (int k = 0; k < nkeys; k++) {
            int var_idx = master_key_indices[k] - 1;  /* Convert to 0-based */

            /* Critical bounds check - prevents crash from bad indices */
            if (var_idx < 0 || (size_t)var_idx >= master_data.nvars) {
                SF_error("cmerge: key variable index out of bounds\n");
                free(master_keys_view.vars);
                if (sort_perm) free(sort_perm);
                stata_data_free(&master_data);
                free(all_var_indices);
                return 459;
            }

            master_keys_view.vars[k] = master_data.vars[var_idx];
        }

        int64_t output_nobs_signed = cmerge_sorted_join(
            &master_keys_view, &g_using_cache.keys,
            nkeys, merge_type, &output_specs);

        /* Free the view (but not the underlying data which belongs to master_data) */
        free(master_keys_view.vars);

        if (output_nobs_signed < 0) {
            if (sort_perm) free(sort_perm);
            stata_data_free(&master_data);
            free(all_var_indices);
            SF_error("cmerge: Merge join failed\n");
            return 920;
        }

        output_nobs = (size_t)output_nobs_signed;
    }

    double t_merge = cmerge_get_time_ms() - t_start;

    /* ===================================================================
     * Step 4: Build master_orig_row mapping using _orig_row variable
     * =================================================================== */

    int64_t *master_orig_rows = malloc(output_nobs * sizeof(int64_t));
    if (!master_orig_rows) {
        free(output_specs);
        if (sort_perm) free(sort_perm);
        stata_data_free(&master_data);
        free(all_var_indices);
        return 920;
    }

    /* _orig_row is at index orig_row_idx-1 in master_data.vars */
    /* Bounds check for orig_row_idx */
    if (orig_row_idx <= 0 || (size_t)(orig_row_idx - 1) >= master_data.nvars) {
        free(output_specs);
        if (sort_perm) free(sort_perm);
        stata_data_free(&master_data);
        free(all_var_indices);
        free(master_orig_rows);
        SF_error("cmerge: orig_row_idx out of bounds\n");
        return 459;
    }

    double *orig_row_data = master_data.vars[orig_row_idx - 1].data.dbl;

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

    double t_reorder = 0;
    if (preserve_order && output_nobs > 0) {
        double t_reorder_start = cmerge_get_time_ms();

        cmerge_order_pair_t *pairs = malloc(output_nobs * sizeof(cmerge_order_pair_t));
        if (!pairs) {
            free(output_specs);
            free(master_orig_rows);
            if (sort_perm) free(sort_perm);
            stata_data_free(&master_data);
            free(all_var_indices);
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

        /* Use radix sort for O(n) performance instead of O(n log n) qsort */
        cmerge_radix_sort_order_pairs(pairs, output_nobs);

        cmerge_output_spec_t *new_specs = malloc(output_nobs * sizeof(cmerge_output_spec_t));
        int64_t *new_orig_rows = malloc(output_nobs * sizeof(int64_t));
        if (!new_specs || !new_orig_rows) {
            free(pairs);
            if (new_specs) free(new_specs);
            if (new_orig_rows) free(new_orig_rows);
            free(output_specs);
            free(master_orig_rows);
            if (sort_perm) free(sort_perm);
            stata_data_free(&master_data);
            free(all_var_indices);
            return 920;
        }

        for (size_t i = 0; i < output_nobs; i++) {
            size_t src_idx = pairs[i].orig_idx;
            new_specs[i] = output_specs[src_idx];
            new_orig_rows[i] = master_orig_rows[src_idx];
        }

        free(output_specs);
        free(master_orig_rows);
        output_specs = new_specs;
        master_orig_rows = new_orig_rows;
        free(pairs);

        t_reorder = cmerge_get_time_ms() - t_reorder_start;
    }

    /* Count merge results */
    size_t n_master_only = 0, n_using_only = 0, n_matched = 0;
    for (size_t i = 0; i < output_nobs; i++) {
        int8_t m = output_specs[i].merge_result;
        n_master_only += (m == 1);
        n_using_only += (m == 2);
        n_matched += (m == 3);
    }

    /* ===================================================================
     * Step 6: Create output data structure with permuted data
     * =================================================================== */

    t_start = cmerge_get_time_ms();

    /* Determine which variables to include in output (exclude _orig_row and shared keepusing) */
    size_t n_output_vars = 0;
    int *output_var_indices = malloc(master_nvars * sizeof(int));
    int *output_var_stata_idx = malloc(master_nvars * sizeof(int));
    int *output_var_is_key = malloc(master_nvars * sizeof(int));
    int *output_var_key_idx = malloc(master_nvars * sizeof(int));

    if (!output_var_indices || !output_var_stata_idx || !output_var_is_key || !output_var_key_idx) {
        if (output_var_indices) free(output_var_indices);
        if (output_var_stata_idx) free(output_var_stata_idx);
        if (output_var_is_key) free(output_var_is_key);
        if (output_var_key_idx) free(output_var_key_idx);
        free(output_specs);
        free(master_orig_rows);
        if (sort_perm) free(sort_perm);
        stata_data_free(&master_data);
        free(all_var_indices);
        return 920;
    }

    for (size_t v = 0; v < master_nvars; v++) {
        int stata_idx = (int)(v + 1);
        if (stata_idx == orig_row_idx) continue;  /* Skip _orig_row */

        /* Check if shared keepusing (handled separately) */
        int is_shared_keepusing = 0;
        for (int kv = 0; kv < n_keepusing; kv++) {
            if (keepusing_placeholder_indices[kv] == stata_idx) {
                is_shared_keepusing = 1;
                break;
            }
        }
        if (is_shared_keepusing) continue;

        output_var_indices[n_output_vars] = (int)v;  /* Index in master_data.vars */
        output_var_stata_idx[n_output_vars] = stata_idx;

        /* Check if key variable */
        output_var_is_key[n_output_vars] = 0;
        output_var_key_idx[n_output_vars] = -1;
        for (int k = 0; k < nkeys; k++) {
            if (master_key_indices[k] == stata_idx) {
                output_var_is_key[n_output_vars] = 1;
                output_var_key_idx[n_output_vars] = k;
                break;
            }
        }

        n_output_vars++;
    }

    /* Create output data structure */
    stata_data output_data;
    stata_data_init(&output_data);
    output_data.nobs = output_nobs;
    output_data.nvars = n_output_vars;
    output_data.vars = (stata_variable *)cmerge_aligned_alloc(n_output_vars * sizeof(stata_variable));
    if (!output_data.vars) {
        free(output_var_indices);
        free(output_var_stata_idx);
        free(output_var_is_key);
        free(output_var_key_idx);
        free(output_specs);
        free(master_orig_rows);
        if (sort_perm) free(sort_perm);
        stata_data_free(&master_data);
        free(all_var_indices);
        return 920;
    }
    memset(output_data.vars, 0, n_output_vars * sizeof(stata_variable));

    /* Count string variables and estimate arena size */
    size_t n_string_vars = 0;
    for (size_t vi = 0; vi < n_output_vars; vi++) {
        int src_var_idx = output_var_indices[vi];
        if (master_data.vars[src_var_idx].type == STATA_TYPE_STRING) {
            n_string_vars++;
        }
    }

    /* Create string arena for merge output: estimate 64 bytes per string cell.
     * Arena enables O(1) bulk free instead of O(n*m) individual frees. */
    cmerge_string_arena *str_arena = NULL;
    if (n_string_vars > 0) {
        size_t arena_size = n_string_vars * output_nobs * 64;
        str_arena = cmerge_arena_create(arena_size);
        /* If arena creation fails, we fall back to strdup (no error) */
    }

    /* Apply permutation to each variable (can be parallelized with OpenMP) */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t vi = 0; vi < n_output_vars; vi++) {
        int src_var_idx = output_var_indices[vi];
        stata_variable *src_var = &master_data.vars[src_var_idx];
        stata_variable *dst_var = &output_data.vars[vi];

        dst_var->nobs = output_nobs;
        dst_var->type = src_var->type;
        dst_var->_arena = NULL;

        if (src_var->type == STATA_TYPE_DOUBLE) {
            dst_var->data.dbl = (double *)cmerge_aligned_alloc(output_nobs * sizeof(double));
            if (!dst_var->data.dbl) continue;  /* Error handled below */

            int is_key = output_var_is_key[vi];
            int key_idx = output_var_key_idx[vi];
            double *using_key_data = (is_key && g_using_cache.loaded && key_idx >= 0) ?
                                      g_using_cache.keys.vars[key_idx].data.dbl : NULL;

            /* Use master_sorted_row to index into sorted master_data */
            for (size_t i = 0; i < output_nobs; i++) {
                int64_t sorted_row = output_specs[i].master_sorted_row;
                if (sorted_row >= 0) {
                    dst_var->data.dbl[i] = src_var->data.dbl[sorted_row];
                } else if (using_key_data != NULL) {
                    int64_t using_row = output_specs[i].using_sorted_row;
                    dst_var->data.dbl[i] = (using_row >= 0) ? using_key_data[using_row] : SV_missval;
                } else {
                    dst_var->data.dbl[i] = SV_missval;
                }
            }
        } else {
            /* String variable - use arena allocator to reduce malloc overhead */
            dst_var->data.str = (char **)calloc(output_nobs, sizeof(char *));
            if (!dst_var->data.str) continue;

            /* Mark this variable as using the shared arena */
            dst_var->_arena = str_arena;

            int is_key = output_var_is_key[vi];
            int key_idx = output_var_key_idx[vi];
            char **using_key_data = (is_key && g_using_cache.loaded && key_idx >= 0) ?
                                     g_using_cache.keys.vars[key_idx].data.str : NULL;

            /* Use master_sorted_row to index into sorted master_data */
            for (size_t i = 0; i < output_nobs; i++) {
                int64_t sorted_row = output_specs[i].master_sorted_row;
                if (sorted_row >= 0 && src_var->data.str[sorted_row]) {
                    dst_var->data.str[i] = cmerge_arena_strdup(str_arena, src_var->data.str[sorted_row]);
                } else if (using_key_data != NULL) {
                    int64_t using_row = output_specs[i].using_sorted_row;
                    if (using_row >= 0 && using_key_data[using_row]) {
                        dst_var->data.str[i] = cmerge_arena_strdup(str_arena, using_key_data[using_row]);
                    } else {
                        dst_var->data.str[i] = cmerge_arena_strdup(str_arena, "");
                    }
                } else {
                    dst_var->data.str[i] = cmerge_arena_strdup(str_arena, "");
                }
            }
        }
    }

    double t_permute = cmerge_get_time_ms() - t_start;

    /* ===================================================================
     * Step 7: Store output data using parallel I/O
     * =================================================================== */

    t_start = cmerge_get_time_ms();

    rc = ctools_data_store_selective(&output_data, output_var_stata_idx, n_output_vars, 1);
    if (rc != STATA_OK) {
        /* Cleanup on error - free only non-arena strings */
        for (size_t vi = 0; vi < n_output_vars; vi++) {
            if (output_data.vars[vi].type == STATA_TYPE_STRING && output_data.vars[vi].data.str) {
                /* Only free strings not owned by arena */
                for (size_t i = 0; i < output_nobs; i++) {
                    if (output_data.vars[vi].data.str[i] &&
                        !cmerge_arena_owns(str_arena, output_data.vars[vi].data.str[i])) {
                        free(output_data.vars[vi].data.str[i]);
                    }
                }
                free(output_data.vars[vi].data.str);
            } else if (output_data.vars[vi].data.dbl) {
                cmerge_aligned_free(output_data.vars[vi].data.dbl);
            }
        }
        cmerge_aligned_free(output_data.vars);
        cmerge_arena_free(str_arena);  /* Free arena in one operation */
        free(output_var_indices);
        free(output_var_stata_idx);
        free(output_var_is_key);
        free(output_var_key_idx);
        free(output_specs);
        free(master_orig_rows);
        if (sort_perm) free(sort_perm);
        stata_data_free(&master_data);
        free(all_var_indices);
        SF_error("cmerge: Failed to store output data\n");
        return rc;
    }

    double t_store = cmerge_get_time_ms() - t_start;

    /* ===================================================================
     * Step 8: Write _merge variable
     * =================================================================== */

    t_start = cmerge_get_time_ms();

    if (merge_var_idx > 0) {
        ST_int mvar = (ST_int)merge_var_idx;
        for (size_t i = 0; i < output_nobs; i++) {
            SF_vstore(mvar, (ST_int)(i + 1), (double)output_specs[i].merge_result);
        }
    }

    /* ===================================================================
     * Step 9: Write keepusing variables from cache (sequential)
     * =================================================================== */

    if (n_keepusing > 0) {
        for (int kv = 0; kv < n_keepusing; kv++) {
            keepusing_write_args_t ku_args;
            ku_args.keepusing_idx = kv;
            ku_args.dest_idx = (ST_int)keepusing_placeholder_indices[kv];
            ku_args.is_shared = shared_flags[kv];
            ku_args.update_mode = update_mode;
            ku_args.replace_mode = replace_mode;
            ku_args.specs = output_specs;
            ku_args.output_nobs = output_nobs;
            ku_args.success = 0;
            write_keepusing_var_thread(&ku_args);
        }
    }

    double t_write_meta = cmerge_get_time_ms() - t_start;

    /* ===================================================================
     * Cleanup
     * =================================================================== */

    double t_cleanup_start = cmerge_get_time_ms();

    /* Free output data - use arena for O(1) string cleanup */
    for (size_t vi = 0; vi < n_output_vars; vi++) {
        if (output_data.vars[vi].type == STATA_TYPE_STRING && output_data.vars[vi].data.str) {
            /* Only free strings that fell back to strdup (not owned by arena) */
            for (size_t i = 0; i < output_nobs; i++) {
                if (output_data.vars[vi].data.str[i] &&
                    !cmerge_arena_owns(str_arena, output_data.vars[vi].data.str[i])) {
                    free(output_data.vars[vi].data.str[i]);
                }
            }
            free(output_data.vars[vi].data.str);
        } else if (output_data.vars[vi].data.dbl) {
            cmerge_aligned_free(output_data.vars[vi].data.dbl);
        }
    }
    cmerge_aligned_free(output_data.vars);
    cmerge_arena_free(str_arena);  /* O(1) cleanup of all arena strings */

    free(output_var_indices);
    free(output_var_stata_idx);
    free(output_var_is_key);
    free(output_var_key_idx);
    free(output_specs);
    free(master_orig_rows);
    if (sort_perm) free(sort_perm);
    stata_data_free(&master_data);
    free(all_var_indices);

    /* Free using cache */
    stata_data_free(&g_using_cache.keys);
    stata_data_free(&g_using_cache.keepusing);
    g_using_cache.loaded = 0;

    double t_cleanup = cmerge_get_time_ms() - t_cleanup_start;

    /* Calculate total time */
    double t_total = cmerge_get_time_ms() - t_total_start;

    /* Save merge results to Stata scalars */
    SF_scal_save("_cmerge_N", (double)output_nobs);
    SF_scal_save("_cmerge_N_1", (double)n_master_only);
    SF_scal_save("_cmerge_N_2", (double)n_using_only);
    SF_scal_save("_cmerge_N_3", (double)n_matched);

    /* Save Phase 2 timing to Stata scalars (in seconds) */
    SF_scal_save("_cmerge_p2_load_master", t_load / 1000.0);
    SF_scal_save("_cmerge_p2_sort_master", t_sort / 1000.0);
    SF_scal_save("_cmerge_p2_merge_join", t_merge / 1000.0);
    SF_scal_save("_cmerge_p2_reorder", t_reorder / 1000.0);
    SF_scal_save("_cmerge_p2_permute", t_permute / 1000.0);
    SF_scal_save("_cmerge_p2_store", t_store / 1000.0);
    SF_scal_save("_cmerge_p2_write_meta", t_write_meta / 1000.0);
    SF_scal_save("_cmerge_p2_cleanup", t_cleanup / 1000.0);
    SF_scal_save("_cmerge_p2_total", t_total / 1000.0);
    SF_scal_save("_cmerge_n_output_vars", (double)n_output_vars);

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

    /* Use static buffer to avoid potential stack issues */
    static char args_copy[8192];
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

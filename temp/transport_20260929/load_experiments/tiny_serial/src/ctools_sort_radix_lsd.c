/*
    ctools_sort_radix_lsd.c
    Parallel LSD Radix Sort Module

    Sorts data loaded from Stata using LSD (Least Significant Digit) radix sort.
    This module is designed to be reusable for any sorting operation on data
    already loaded into stata_data structures.

    Algorithm:
    - LSD radix sort: processes bytes from least to most significant
    - 8-bit radix (256 buckets per pass)
    - 8 passes for 64-bit doubles
    - Multi-key stable sort: keys processed from last to first

    Numeric Sorting:
    - IEEE 754 doubles converted to sortable uint64 via bit manipulation
    - Negative numbers: flip all bits (makes them sort before positives)
    - Positive numbers: flip sign bit only
    - Stata missing values: sort to end, in Stata's . < .a < ... < .z order

    String Sorting:
    - Byte-by-byte LSD sort from rightmost character
    - Shorter strings padded with zeros (sort before longer strings)
    - Pre-cached string lengths to avoid repeated strlen() calls

    Parallelization:
    - Parallel histogram: each thread counts its chunk, then combine
    - Parallel scatter: each thread writes to pre-computed offsets
    - Threshold: MIN_OBS_PER_THREAD * 2 observations to enable parallel
    - Reuses thread allocations across all 8 radix passes

    Optimizations:
    - Pointer swapping instead of memcpy between passes
    - Early exit for uniform byte distributions (no-op passes)
    - Aligned memory allocation for better cache/SIMD performance
    - Pre-cached string lengths for string sorting

    Post-Sort Permutation:
    - After sorting, applies permutation to all variables in parallel
    - One thread per variable for cache-friendly access pattern
    - Data becomes physically sorted (sequential access for store phase)
*/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_threads.h"

/*
    NOTE: Aligned memory allocation is now provided by ctools_config.h
    Use ctools_aligned_alloc() and ctools_aligned_free() for cross-platform support.
    On Windows, _aligned_malloc requires _aligned_free - using regular free() causes heap corruption.

    NOTE: Double-to-sortable conversion is now provided by ctools_types.h
    Use ctools_double_to_sortable(d) for sorting doubles.
*/

/* ============================================================================
   Parallel radix sort structures and functions

   Permutation arrays use perm_idx_t (uint32_t) instead of size_t for:
   - 50% memory reduction (4 bytes vs 8 bytes per index)
   - Better cache utilization during sorting
   - Sufficient for Stata's 2^31 observation limit
   ============================================================================ */

/*
    Reusable allocation structure for parallel radix sort.
    Allocated once and reused across all 8 passes to avoid repeated malloc/free.
*/
typedef struct {
    int num_threads;
    size_t **all_local_counts;   /* [num_threads][RADIX_SIZE] */
    size_t *global_counts;       /* [RADIX_SIZE] */
    size_t *global_offsets;      /* [RADIX_SIZE] */
    size_t **thread_offsets;     /* [num_threads][RADIX_SIZE] */
} radix_sort_context_t;

/*
    Allocate reusable context for parallel radix sort.
    Returns NULL on allocation failure.
*/
static radix_sort_context_t *radix_context_alloc(int num_threads)
{
    radix_sort_context_t *ctx;
    int t;

    ctx = (radix_sort_context_t *)calloc(1, sizeof(radix_sort_context_t));
    if (ctx == NULL) return NULL;

    ctx->num_threads = num_threads;
    ctx->all_local_counts = (size_t **)calloc(num_threads, sizeof(size_t *));
    ctx->global_counts = (size_t *)calloc(RADIX_SIZE, sizeof(size_t));
    ctx->global_offsets = (size_t *)malloc(RADIX_SIZE * sizeof(size_t));
    ctx->thread_offsets = (size_t **)calloc(num_threads, sizeof(size_t *));

    if (!ctx->all_local_counts || !ctx->global_counts ||
        !ctx->global_offsets || !ctx->thread_offsets) {
        goto fail;
    }

    /* Initialize pointer arrays to NULL for safe cleanup */
    for (t = 0; t < num_threads; t++) {
        ctx->all_local_counts[t] = NULL;
        ctx->thread_offsets[t] = NULL;
    }

    /* Allocate per-thread arrays */
    for (t = 0; t < num_threads; t++) {
        ctx->all_local_counts[t] = (size_t *)calloc(RADIX_SIZE, sizeof(size_t));
        ctx->thread_offsets[t] = (size_t *)malloc(RADIX_SIZE * sizeof(size_t));
        if (!ctx->all_local_counts[t] || !ctx->thread_offsets[t]) {
            goto fail;
        }
    }

    return ctx;

fail:
    if (ctx) {
        if (ctx->all_local_counts) {
            for (t = 0; t < num_threads; t++) {
                free(ctx->all_local_counts[t]);
            }
        }
        if (ctx->thread_offsets) {
            for (t = 0; t < num_threads; t++) {
                free(ctx->thread_offsets[t]);
            }
        }
        free(ctx->all_local_counts);
        free(ctx->global_counts);
        free(ctx->global_offsets);
        free(ctx->thread_offsets);
        free(ctx);
    }
    return NULL;
}

/*
    Free radix sort context.
*/
static void radix_context_free(radix_sort_context_t *ctx)
{
    int t;
    if (ctx == NULL) return;

    for (t = 0; t < ctx->num_threads; t++) {
        free(ctx->all_local_counts[t]);
        free(ctx->thread_offsets[t]);
    }
    free(ctx->all_local_counts);
    free(ctx->global_counts);
    free(ctx->global_offsets);
    free(ctx->thread_offsets);
    free(ctx);
}

/*
    Parallel radix sort pass for numeric data.
    Uses reusable context to avoid per-pass allocations.
    Returns 1 if pass was skipped (uniform distribution), 0 otherwise.
*/
static int radix_sort_pass_numeric_parallel(perm_idx_t *order,
                                            perm_idx_t *temp_order,
                                            uint64_t *keys,
                                            size_t nobs,
                                            int byte_pos,
                                            radix_sort_context_t *ctx)
{
    size_t chunk_size;
    int shift = byte_pos * RADIX_BITS;
    int t, b;
    int num_threads = ctx->num_threads;

    /* Phase 1: Parallel histogram computation (OpenMP) */
    chunk_size = (nobs + num_threads - 1) / num_threads;

    #pragma omp parallel for num_threads(num_threads) schedule(static)
    for (int tid = 0; tid < num_threads; tid++) {
        size_t my_start = (size_t)tid * chunk_size;
        size_t my_end = (size_t)(tid + 1) * chunk_size;
        if (my_end > nobs) my_end = nobs;
        if (my_start >= nobs) my_start = my_end = nobs;

        size_t *my_counts = ctx->all_local_counts[tid];
        memset(my_counts, 0, RADIX_SIZE * sizeof(size_t));

        /* 4x unrolled histogram with prefetching */
        size_t i = my_start;
        for (; i + 4 <= my_end; i += 4) {
            if (i + PREFETCH_DISTANCE < my_end) {
                CTOOLS_PREFETCH(&keys[order[i + PREFETCH_DISTANCE]]);
            }
            uint8_t b0 = (keys[order[i+0]] >> shift) & RADIX_MASK;
            uint8_t b1 = (keys[order[i+1]] >> shift) & RADIX_MASK;
            uint8_t b2 = (keys[order[i+2]] >> shift) & RADIX_MASK;
            uint8_t b3 = (keys[order[i+3]] >> shift) & RADIX_MASK;
            my_counts[b0]++;
            my_counts[b1]++;
            my_counts[b2]++;
            my_counts[b3]++;
        }
        for (; i < my_end; i++) {
            my_counts[(keys[order[i]] >> shift) & RADIX_MASK]++;
        }
    }

    /* Combine local histograms into global counts */
    /* Optimization: iterate bucket-first for better cache locality */
    /* Memory layout is all_local_counts[thread][bucket], so b-outer is cache-friendly */
    memset(ctx->global_counts, 0, RADIX_SIZE * sizeof(size_t));
    for (b = 0; b < RADIX_SIZE; b++) {
        for (t = 0; t < num_threads; t++) {
            ctx->global_counts[b] += ctx->all_local_counts[t][b];
        }
    }

    /* Optimization 3: Check for uniform distribution (all in one bucket) */
    {
        int non_empty_buckets = 0;
        for (b = 0; b < RADIX_SIZE; b++) {
            if (ctx->global_counts[b] > 0) {
                non_empty_buckets++;
                if (non_empty_buckets > 1) break;
            }
        }
        if (non_empty_buckets <= 1) {
            /* All elements in same bucket - this pass is a no-op */
            return 1;
        }
    }

    /* Compute global prefix sums (bucket starting positions) */
    ctx->global_offsets[0] = 0;
    for (b = 1; b < RADIX_SIZE; b++) {
        ctx->global_offsets[b] = ctx->global_offsets[b - 1] + ctx->global_counts[b - 1];
    }

    /* Compute per-thread offsets within each bucket */
    for (b = 0; b < RADIX_SIZE; b++) {
        size_t offset = ctx->global_offsets[b];
        for (t = 0; t < num_threads; t++) {
            ctx->thread_offsets[t][b] = offset;
            offset += ctx->all_local_counts[t][b];
        }
    }

    /* Phase 2: Parallel scatter (OpenMP) */
    #pragma omp parallel for num_threads(num_threads) schedule(static)
    for (int tid = 0; tid < num_threads; tid++) {
        size_t my_start = (size_t)tid * chunk_size;
        size_t my_end = (size_t)(tid + 1) * chunk_size;
        if (my_end > nobs) my_end = nobs;
        if (my_start >= nobs) my_start = my_end = nobs;

        size_t local_offsets[RADIX_SIZE];
        memcpy(local_offsets, ctx->thread_offsets[tid], RADIX_SIZE * sizeof(size_t));

        for (size_t i = my_start; i < my_end; i++) {
            uint8_t byte_val = (keys[order[i]] >> shift) & RADIX_MASK;
            temp_order[local_offsets[byte_val]++] = order[i];
        }
    }

    return 0;  /* Pass was not skipped */
}

/*
    Sequential radix sort pass (for small datasets or fallback).
    Returns 1 if pass was skipped (uniform distribution), 0 otherwise.
*/
static int radix_sort_pass_numeric(perm_idx_t *order,
                                   perm_idx_t *temp_order,
                                   uint64_t *keys,
                                   size_t nobs,
                                   int byte_pos)
{
    size_t counts[RADIX_SIZE] = {0};
    size_t offsets[RADIX_SIZE];
    size_t i;
    uint8_t byte_val;
    int shift = byte_pos * RADIX_BITS;
    int non_empty_buckets = 0;

    /* Count occurrences of each byte value */
    for (i = 0; i < nobs; i++) {
        byte_val = (keys[order[i]] >> shift) & RADIX_MASK;
        counts[byte_val]++;
    }

    /* Optimization 3: Check for uniform distribution */
    for (i = 0; i < RADIX_SIZE; i++) {
        if (counts[i] > 0) {
            non_empty_buckets++;
            if (non_empty_buckets > 1) break;
        }
    }
    if (non_empty_buckets <= 1) {
        return 1;  /* Skip this pass */
    }

    /* Compute prefix sums (starting positions for each bucket) */
    offsets[0] = 0;
    for (i = 1; i < RADIX_SIZE; i++) {
        offsets[i] = offsets[i - 1] + counts[i - 1];
    }

    /* Place elements in sorted order */
    for (i = 0; i < nobs; i++) {
        byte_val = (keys[order[i]] >> shift) & RADIX_MASK;
        temp_order[offsets[byte_val]++] = order[i];
    }

    return 0;  /* Pass was not skipped */
}

/* ============================================================================
   String sorting with pre-cached lengths and parallel histogram/scatter
   ============================================================================ */

/*
    String sort context for reusable allocations.
*/
typedef struct {
    int num_threads;
    size_t **all_local_counts;  /* [num_threads][RADIX_SIZE+1] */
    size_t *global_counts;      /* [RADIX_SIZE+1] */
    size_t *global_offsets;     /* [RADIX_SIZE+1] */
    size_t **thread_offsets;    /* [num_threads][RADIX_SIZE+1] */
} string_sort_context_t;

static string_sort_context_t *string_context_alloc(int num_threads)
{
    string_sort_context_t *ctx;
    int t;

    ctx = (string_sort_context_t *)calloc(1, sizeof(string_sort_context_t));
    if (ctx == NULL) return NULL;

    ctx->num_threads = num_threads;
    ctx->all_local_counts = (size_t **)calloc(num_threads, sizeof(size_t *));
    ctx->global_counts = (size_t *)calloc(RADIX_SIZE + 1, sizeof(size_t));
    ctx->global_offsets = (size_t *)malloc((RADIX_SIZE + 1) * sizeof(size_t));
    ctx->thread_offsets = (size_t **)calloc(num_threads, sizeof(size_t *));

    if (!ctx->all_local_counts || !ctx->global_counts ||
        !ctx->global_offsets || !ctx->thread_offsets) {
        goto fail;
    }

    for (t = 0; t < num_threads; t++) {
        ctx->all_local_counts[t] = NULL;
        ctx->thread_offsets[t] = NULL;
    }

    for (t = 0; t < num_threads; t++) {
        ctx->all_local_counts[t] = (size_t *)calloc(RADIX_SIZE + 1, sizeof(size_t));
        ctx->thread_offsets[t] = (size_t *)malloc((RADIX_SIZE + 1) * sizeof(size_t));
        if (!ctx->all_local_counts[t] || !ctx->thread_offsets[t]) {
            goto fail;
        }
    }

    return ctx;

fail:
    if (ctx) {
        if (ctx->all_local_counts) {
            for (t = 0; t < num_threads; t++) {
                free(ctx->all_local_counts[t]);
            }
        }
        if (ctx->thread_offsets) {
            for (t = 0; t < num_threads; t++) {
                free(ctx->thread_offsets[t]);
            }
        }
        free(ctx->all_local_counts);
        free(ctx->global_counts);
        free(ctx->global_offsets);
        free(ctx->thread_offsets);
        free(ctx);
    }
    return NULL;
}

static void string_context_free(string_sort_context_t *ctx)
{
    int t;
    if (ctx == NULL) return;

    for (t = 0; t < ctx->num_threads; t++) {
        free(ctx->all_local_counts[t]);
        free(ctx->thread_offsets[t]);
    }
    free(ctx->all_local_counts);
    free(ctx->global_counts);
    free(ctx->global_offsets);
    free(ctx->thread_offsets);
    free(ctx);
}

/*
    Parallel string radix sort pass.
    Returns 1 if skipped (uniform), 0 otherwise.
*/
static int radix_sort_pass_string_parallel(perm_idx_t *order,
                                           perm_idx_t *temp_order,
                                           char **strings,
                                           size_t *str_lengths,
                                           size_t nobs,
                                           size_t char_pos,
                                           string_sort_context_t *ctx)
{
    size_t chunk_size;
    int t, b;
    int num_threads = ctx->num_threads;

    chunk_size = (nobs + num_threads - 1) / num_threads;

    /* Phase 1: Parallel histogram (OpenMP) */
    #pragma omp parallel for num_threads(num_threads) schedule(static)
    for (int tid = 0; tid < num_threads; tid++) {
        size_t my_start = (size_t)tid * chunk_size;
        size_t my_end = (size_t)(tid + 1) * chunk_size;
        if (my_end > nobs) my_end = nobs;
        if (my_start >= nobs) my_start = my_end = nobs;

        size_t *my_counts = ctx->all_local_counts[tid];
        memset(my_counts, 0, (RADIX_SIZE + 1) * sizeof(size_t));

        for (size_t i = my_start; i < my_end; i++) {
            perm_idx_t idx = order[i];
            size_t len = str_lengths[idx];
            unsigned char byte_val;
            if (char_pos < len && strings[idx] != NULL) {
                byte_val = (unsigned char)strings[idx][char_pos];
            } else {
                byte_val = 0;
            }
            my_counts[byte_val]++;
        }
    }

    /* Combine histograms */
    /* Optimization: iterate bucket-first for better cache locality */
    memset(ctx->global_counts, 0, (RADIX_SIZE + 1) * sizeof(size_t));
    for (b = 0; b <= RADIX_SIZE; b++) {
        for (t = 0; t < num_threads; t++) {
            ctx->global_counts[b] += ctx->all_local_counts[t][b];
        }
    }

    /* Check for uniform distribution */
    {
        int non_empty = 0;
        for (b = 0; b <= RADIX_SIZE; b++) {
            if (ctx->global_counts[b] > 0) {
                non_empty++;
                if (non_empty > 1) break;
            }
        }
        if (non_empty <= 1) {
            return 1;
        }
    }

    /* Compute prefix sums */
    ctx->global_offsets[0] = 0;
    for (b = 1; b <= RADIX_SIZE; b++) {
        ctx->global_offsets[b] = ctx->global_offsets[b - 1] + ctx->global_counts[b - 1];
    }

    /* Compute per-thread offsets */
    for (b = 0; b <= RADIX_SIZE; b++) {
        size_t offset = ctx->global_offsets[b];
        for (t = 0; t < num_threads; t++) {
            ctx->thread_offsets[t][b] = offset;
            offset += ctx->all_local_counts[t][b];
        }
    }

    /* Phase 2: Parallel scatter (OpenMP) */
    #pragma omp parallel for num_threads(num_threads) schedule(static)
    for (int tid = 0; tid < num_threads; tid++) {
        size_t my_start = (size_t)tid * chunk_size;
        size_t my_end = (size_t)(tid + 1) * chunk_size;
        if (my_end > nobs) my_end = nobs;
        if (my_start >= nobs) my_start = my_end = nobs;

        size_t local_offsets[RADIX_SIZE + 1];
        memcpy(local_offsets, ctx->thread_offsets[tid], (RADIX_SIZE + 1) * sizeof(size_t));

        for (size_t i = my_start; i < my_end; i++) {
            perm_idx_t idx = order[i];
            size_t len = str_lengths[idx];
            unsigned char byte_val;
            if (char_pos < len && strings[idx] != NULL) {
                byte_val = (unsigned char)strings[idx][char_pos];
            } else {
                byte_val = 0;
            }
            temp_order[local_offsets[byte_val]++] = idx;
        }
    }

    return 0;
}

/*
    Sequential string radix sort pass with pre-cached lengths.
    Returns 1 if skipped, 0 otherwise.
*/
static int radix_sort_pass_string(perm_idx_t *order,
                                  perm_idx_t *temp_order,
                                  char **strings,
                                  size_t *str_lengths,
                                  size_t nobs,
                                  size_t char_pos)
{
    size_t counts[RADIX_SIZE + 1] = {0};
    size_t offsets[RADIX_SIZE + 1];
    size_t i, idx, len;
    unsigned char byte_val;
    int non_empty = 0;

    /* Count occurrences using pre-cached lengths */
    for (i = 0; i < nobs; i++) {
        idx = order[i];
        len = str_lengths[idx];
        if (char_pos < len && strings[idx] != NULL) {
            byte_val = (unsigned char)strings[idx][char_pos];
        } else {
            byte_val = 0;
        }
        counts[byte_val]++;
    }

    /* Check for uniform distribution */
    for (i = 0; i <= RADIX_SIZE; i++) {
        if (counts[i] > 0) {
            non_empty++;
            if (non_empty > 1) break;
        }
    }
    if (non_empty <= 1) {
        return 1;
    }

    /* Compute prefix sums */
    offsets[0] = 0;
    for (i = 1; i <= RADIX_SIZE; i++) {
        offsets[i] = offsets[i - 1] + counts[i - 1];
    }

    /* Scatter */
    for (i = 0; i < nobs; i++) {
        idx = order[i];
        len = str_lengths[idx];
        if (char_pos < len && strings[idx] != NULL) {
            byte_val = (unsigned char)strings[idx][char_pos];
        } else {
            byte_val = 0;
        }
        temp_order[offsets[byte_val]++] = idx;
    }

    return 0;
}

/*
    Sort by a single numeric variable using radix sort.
    Uses parallel implementation for large datasets.
    Optimization 1: Swaps pointers instead of memcpy.
    Optimization 3: Skips uniform passes.
    Optimization 4: Reuses allocations across passes.
    Optimization 6: Uses aligned memory for keys.
*/
static stata_retcode sort_by_numeric_var(stata_data *data, int var_idx)
{
    perm_idx_t *order_a;
    perm_idx_t *order_b;
    perm_idx_t *current_order;
    perm_idx_t *temp_order;
    uint64_t *keys;
    size_t i;
    int byte_pos;
    double *dbl_data;
    int use_parallel;
    int num_threads;
    radix_sort_context_t *ctx = NULL;
    int swapped = 0;

    /* Allocate temporary order array (uses perm_idx_t for 50% memory savings) */
    order_a = data->sort_order;
    order_b = (perm_idx_t *)ctools_aligned_alloc(CACHE_LINE_SIZE, data->nobs * sizeof(perm_idx_t));

    /* Allocate aligned keys array */
    keys = (uint64_t *)ctools_aligned_alloc(CACHE_LINE_SIZE, data->nobs * sizeof(uint64_t));

    if (order_b == NULL || keys == NULL) {
        ctools_aligned_free(order_b);
        ctools_aligned_free(keys);
        return STATA_ERR_MEMORY;
    }

    /* Convert doubles to sortable uint64 keys (parallel for large datasets) */
    dbl_data = data->vars[var_idx].data.dbl;
    #pragma omp parallel for schedule(static) if(data->nobs >= MIN_OBS_PER_THREAD * 2)
    for (i = 0; i < data->nobs; i++) {
        keys[i] = ctools_double_to_sortable(dbl_data[i]);
    }

    /* Decide whether to use parallel sort */
    use_parallel = (data->nobs >= MIN_OBS_PER_THREAD * 2);
    num_threads = ctools_get_openmp_threads();
    if (data->nobs < (size_t)MIN_OBS_PER_THREAD * (size_t)num_threads) {
        num_threads = (int)(data->nobs / MIN_OBS_PER_THREAD);
        if (num_threads < 2) {
            use_parallel = 0;
        }
    }

    /* Allocate reusable context for parallel sort */
    if (use_parallel) {
        ctx = radix_context_alloc(num_threads);
        if (ctx == NULL) {
            ctools_aligned_free(order_b);
            ctools_aligned_free(keys);
            return STATA_ERR_MEMORY;
        }
    }

    /* Start with order_a as current, order_b as temp */
    current_order = order_a;
    temp_order = order_b;

    /* Perform radix sort passes (LSD - from least significant to most) */
    for (byte_pos = 0; byte_pos < 8; byte_pos++) {
        int skipped;

        if (use_parallel) {
            skipped = radix_sort_pass_numeric_parallel(current_order, temp_order, keys,
                                                        data->nobs, byte_pos, ctx);
        } else {
            skipped = radix_sort_pass_numeric(current_order, temp_order, keys,
                                              data->nobs, byte_pos);
        }

        /* Optimization 1: Swap pointers instead of memcpy */
        if (!skipped) {
            perm_idx_t *tmp = current_order;
            current_order = temp_order;
            temp_order = tmp;
            swapped = !swapped;
        }
    }

    /* Ensure final result is in data->sort_order */
    if (current_order != data->sort_order) {
        memcpy(data->sort_order, current_order, data->nobs * sizeof(perm_idx_t));
    }

    radix_context_free(ctx);
    ctools_aligned_free(order_b);
    ctools_aligned_free(keys);
    return STATA_OK;
}

/*
    Sort by a single string variable using radix sort.
    Optimization 2: Pre-caches string lengths.
    Optimization 5: Uses parallel sort for large datasets.
*/
static stata_retcode sort_by_string_var(stata_data *data, int var_idx)
{
    perm_idx_t *order_a;
    perm_idx_t *order_b;
    perm_idx_t *current_order;
    perm_idx_t *temp_order;
    size_t *str_lengths;
    size_t max_len = 0;
    size_t i, len;
    int char_pos;
    char **str_data;
    int use_parallel;
    int num_threads;
    string_sort_context_t *ctx = NULL;
    int swapped = 0;

    str_data = data->vars[var_idx].data.str;

    /* Safety check for NULL string data array */
    if (str_data == NULL) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Optimization 2: Pre-cache all string lengths */
    str_lengths = (size_t *)malloc(data->nobs * sizeof(size_t));
    if (str_lengths == NULL) {
        return STATA_ERR_MEMORY;
    }

    for (i = 0; i < data->nobs; i++) {
        /* Handle NULL string pointers gracefully */
        if (str_data[i] == NULL) {
            len = 0;
        } else {
            len = strlen(str_data[i]);
        }
        str_lengths[i] = len;
        if (len > max_len) {
            max_len = len;
        }
    }

    if (max_len == 0) {
        free(str_lengths);
        return STATA_OK;
    }

    /* Allocate temporary order array - aligned for cache efficiency (uses perm_idx_t) */
    order_a = data->sort_order;
    order_b = (perm_idx_t *)ctools_aligned_alloc(CACHE_LINE_SIZE, data->nobs * sizeof(perm_idx_t));
    if (order_b == NULL) {
        free(str_lengths);
        return STATA_ERR_MEMORY;
    }

    /* Decide on parallelization */
    use_parallel = (data->nobs >= MIN_OBS_PER_THREAD * 2);
    num_threads = ctools_get_openmp_threads();
    if (data->nobs < (size_t)MIN_OBS_PER_THREAD * (size_t)num_threads) {
        num_threads = (int)(data->nobs / MIN_OBS_PER_THREAD);
        if (num_threads < 2) {
            use_parallel = 0;
        }
    }

    if (use_parallel) {
        ctx = string_context_alloc(num_threads);
        if (ctx == NULL) {
            ctools_aligned_free(order_b);
            free(str_lengths);
            return STATA_ERR_MEMORY;
        }
    }

    current_order = order_a;
    temp_order = order_b;

    /* Perform radix sort passes (LSD - from rightmost character to leftmost) */
    for (char_pos = (int)max_len - 1; char_pos >= 0; char_pos--) {
        int skipped;

        if (use_parallel) {
            skipped = radix_sort_pass_string_parallel(current_order, temp_order, str_data,
                                                       str_lengths, data->nobs, (size_t)char_pos, ctx);
        } else {
            skipped = radix_sort_pass_string(current_order, temp_order, str_data,
                                             str_lengths, data->nobs, (size_t)char_pos);
        }

        if (!skipped) {
            perm_idx_t *tmp = current_order;
            current_order = temp_order;
            temp_order = tmp;
            swapped = !swapped;
        }
    }

    /* Ensure final result is in data->sort_order */
    if (current_order != data->sort_order) {
        memcpy(data->sort_order, current_order, data->nobs * sizeof(perm_idx_t));
    }

    string_context_free(ctx);
    ctools_aligned_free(order_b);
    free(str_lengths);
    return STATA_OK;
}

/*
    ctools_sort_radix_lsd_order_only - LSD radix sort without applying permutation

    Computes sort_order but does NOT apply the permutation to data.
    After this call, data->sort_order contains the permutation but data is unchanged.
    Call ctools_apply_permutation() separately to apply the permutation.

    @param data       Dataset to sort
    @param sort_vars  Array of 1-based variable indices to sort by
    @param nsort      Number of sort variables

    @return STATA_OK on success, or error code
*/
stata_retcode ctools_sort_radix_lsd_order_only(stata_data *data, int *sort_vars, size_t nsort)
{
    int k;
    int var_idx;
    stata_retcode rc;

    if (data == NULL || sort_vars == NULL || data->nobs == 0 || nsort == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    /*
        For stable LSD radix sort with multiple keys:
        Sort from the LAST (least significant) key to the FIRST (most significant).
    */
    for (k = (int)nsort - 1; k >= 0; k--) {
        var_idx = sort_vars[k] - 1;  /* Convert to 0-based index */

        if (var_idx < 0 || var_idx >= (int)data->nvars) {
            return STATA_ERR_INVALID_INPUT;
        }

        if (data->vars[var_idx].type == STATA_TYPE_DOUBLE) {
            rc = sort_by_numeric_var(data, var_idx);
        } else {
            rc = sort_by_string_var(data, var_idx);
        }

        if (rc != STATA_OK) {
            return rc;
        }
    }

    /* Do NOT apply permutation - caller will do it separately */
    return STATA_OK;
}

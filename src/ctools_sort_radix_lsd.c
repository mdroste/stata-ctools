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
    - Stata missing values: sort to end (mapped to UINT64_MAX)

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
#include <pthread.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"

/*
    Aligned memory allocation wrapper.
    Falls back to regular malloc if posix_memalign unavailable.
*/
static void *aligned_alloc_wrapper(size_t alignment, size_t size)
{
    void *ptr = NULL;
#if defined(__APPLE__) || defined(__linux__)
    if (posix_memalign(&ptr, alignment, size) != 0) {
        return NULL;
    }
    return ptr;
#else
    /* Fallback for systems without posix_memalign */
    return malloc(size);
#endif
}

/*
    Convert IEEE 754 double to uint64 that sorts correctly with unsigned comparison.

    Problem: Raw IEEE 754 bit patterns don't sort correctly as unsigned integers:
    - Negative numbers have sign bit set (MSB=1), so they appear "larger"
    - Negative magnitudes are inverted (e.g., -2.0 > -1.0 in raw bits)

    Solution: Transform bits so unsigned comparison produces correct numeric order.
    - Positive: flip sign bit only (0x8... → 0x0...)
    - Negative: flip ALL bits (~bits)
    - Missing: map to UINT64_MAX (sorts to end)

    @param d  Input double value
    @return   Sortable uint64 representation
*/
static uint64_t double_to_sortable_uint64(double d)
{
    uint64_t bits;
    memcpy(&bits, &d, sizeof(bits));

    /* Check if the value is a Stata missing value */
    if (SF_is_missing(d)) {
        /* Sort missing values to the end (largest possible value) */
        return UINT64_MAX;
    }

    /* If negative (sign bit set), flip all bits */
    /* If positive, flip only the sign bit */
    if (bits & ((uint64_t)1 << 63)) {
        bits = ~bits;
    } else {
        bits ^= ((uint64_t)1 << 63);
    }

    return bits;
}

/* ============================================================================
   Parallel radix sort structures and functions
   ============================================================================ */

/* Thread argument structure for parallel histogram */
typedef struct {
    size_t *order;           /* Input order array */
    uint64_t *keys;          /* Key array */
    size_t start;            /* Start index for this thread */
    size_t end;              /* End index (exclusive) for this thread */
    int shift;               /* Bit shift for current byte */
    size_t *local_counts;    /* Thread-local histogram (RADIX_SIZE elements) */
} histogram_args_t;

/* Thread argument structure for parallel scatter */
typedef struct {
    size_t *order;           /* Input order array */
    size_t *temp_order;      /* Output order array */
    uint64_t *keys;          /* Key array */
    size_t start;            /* Start index for this thread */
    size_t end;              /* End index (exclusive) for this thread */
    int shift;               /* Bit shift for current byte */
    size_t *global_offsets;  /* Global starting offsets for each bucket */
    size_t *local_offsets;   /* This thread's starting offset within each bucket */
} scatter_args_t;

/*
    Reusable allocation structure for parallel radix sort.
    Allocated once and reused across all 8 passes to avoid repeated malloc/free.
*/
typedef struct {
    int num_threads;
    pthread_t *threads;
    histogram_args_t *hist_args;
    scatter_args_t *scatter_args;
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
    ctx->threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
    ctx->hist_args = (histogram_args_t *)malloc(num_threads * sizeof(histogram_args_t));
    ctx->scatter_args = (scatter_args_t *)malloc(num_threads * sizeof(scatter_args_t));
    ctx->all_local_counts = (size_t **)malloc(num_threads * sizeof(size_t *));
    ctx->global_counts = (size_t *)calloc(RADIX_SIZE, sizeof(size_t));
    ctx->global_offsets = (size_t *)malloc(RADIX_SIZE * sizeof(size_t));
    ctx->thread_offsets = (size_t **)malloc(num_threads * sizeof(size_t *));

    if (!ctx->threads || !ctx->hist_args || !ctx->scatter_args ||
        !ctx->all_local_counts || !ctx->global_counts ||
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
        free(ctx->threads);
        free(ctx->hist_args);
        free(ctx->scatter_args);
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
    free(ctx->threads);
    free(ctx->hist_args);
    free(ctx->scatter_args);
    free(ctx->all_local_counts);
    free(ctx->global_counts);
    free(ctx->global_offsets);
    free(ctx->thread_offsets);
    free(ctx);
}

/* Thread function: compute local histogram */
static void *histogram_thread(void *arg)
{
    histogram_args_t *args = (histogram_args_t *)arg;
    size_t i;
    uint8_t byte_val;
    size_t *order = args->order;
    uint64_t *keys = args->keys;
    int shift = args->shift;

    /* Zero out local counts */
    memset(args->local_counts, 0, RADIX_SIZE * sizeof(size_t));

    /* Count occurrences in this thread's range */
    for (i = args->start; i < args->end; i++) {
        byte_val = (keys[order[i]] >> shift) & RADIX_MASK;
        args->local_counts[byte_val]++;
    }

    return NULL;
}

/* Thread function: scatter elements to sorted positions */
static void *scatter_thread(void *arg)
{
    scatter_args_t *args = (scatter_args_t *)arg;
    size_t i;
    uint8_t byte_val;
    size_t local_offsets[RADIX_SIZE];  /* Stack allocation for speed */
    size_t *order = args->order;
    uint64_t *keys = args->keys;
    size_t *temp_order = args->temp_order;
    int shift = args->shift;

    /* Copy offsets to local array */
    memcpy(local_offsets, args->local_offsets, RADIX_SIZE * sizeof(size_t));

    /* Scatter elements from this thread's range */
    for (i = args->start; i < args->end; i++) {
        byte_val = (keys[order[i]] >> shift) & RADIX_MASK;
        temp_order[local_offsets[byte_val]++] = order[i];
    }

    return NULL;
}

/*
    Parallel radix sort pass for numeric data.
    Uses reusable context to avoid per-pass allocations.
    Returns 1 if pass was skipped (uniform distribution), 0 otherwise.
*/
static int radix_sort_pass_numeric_parallel(size_t *order,
                                            size_t *temp_order,
                                            uint64_t *keys,
                                            size_t nobs,
                                            int byte_pos,
                                            radix_sort_context_t *ctx)
{
    size_t chunk_size, start, end;
    int shift = byte_pos * RADIX_BITS;
    int t, b;
    int num_threads = ctx->num_threads;

    /* Phase 1: Parallel histogram computation */
    chunk_size = (nobs + num_threads - 1) / num_threads;
    for (t = 0; t < num_threads; t++) {
        start = (size_t)t * chunk_size;
        end = (size_t)(t + 1) * chunk_size;
        if (end > nobs) end = nobs;
        if (start >= nobs) start = end = nobs;

        ctx->hist_args[t].order = order;
        ctx->hist_args[t].keys = keys;
        ctx->hist_args[t].start = start;
        ctx->hist_args[t].end = end;
        ctx->hist_args[t].shift = shift;
        ctx->hist_args[t].local_counts = ctx->all_local_counts[t];

        pthread_create(&ctx->threads[t], NULL, histogram_thread, &ctx->hist_args[t]);
    }

    /* Wait for histogram threads */
    for (t = 0; t < num_threads; t++) {
        pthread_join(ctx->threads[t], NULL);
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

    /* Phase 2: Parallel scatter */
    for (t = 0; t < num_threads; t++) {
        start = (size_t)t * chunk_size;
        end = (size_t)(t + 1) * chunk_size;
        if (end > nobs) end = nobs;
        if (start >= nobs) start = end = nobs;

        ctx->scatter_args[t].order = order;
        ctx->scatter_args[t].temp_order = temp_order;
        ctx->scatter_args[t].keys = keys;
        ctx->scatter_args[t].start = start;
        ctx->scatter_args[t].end = end;
        ctx->scatter_args[t].shift = shift;
        ctx->scatter_args[t].global_offsets = ctx->global_offsets;
        ctx->scatter_args[t].local_offsets = ctx->thread_offsets[t];

        pthread_create(&ctx->threads[t], NULL, scatter_thread, &ctx->scatter_args[t]);
    }

    /* Wait for scatter threads */
    for (t = 0; t < num_threads; t++) {
        pthread_join(ctx->threads[t], NULL);
    }

    return 0;  /* Pass was not skipped */
}

/*
    Sequential radix sort pass (for small datasets or fallback).
    Returns 1 if pass was skipped (uniform distribution), 0 otherwise.
*/
static int radix_sort_pass_numeric(size_t *order,
                                   size_t *temp_order,
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

/* Thread argument for parallel string histogram */
typedef struct {
    size_t *order;
    char **strings;
    size_t *str_lengths;     /* Pre-cached string lengths */
    size_t start;
    size_t end;
    size_t char_pos;
    size_t *local_counts;
} string_histogram_args_t;

/* Thread argument for parallel string scatter */
typedef struct {
    size_t *order;
    size_t *temp_order;
    char **strings;
    size_t *str_lengths;
    size_t start;
    size_t end;
    size_t char_pos;
    size_t *local_offsets;
} string_scatter_args_t;

/* Thread function: compute string histogram */
static void *string_histogram_thread(void *arg)
{
    string_histogram_args_t *args = (string_histogram_args_t *)arg;
    size_t i;
    unsigned char byte_val;
    size_t idx, len;

    memset(args->local_counts, 0, (RADIX_SIZE + 1) * sizeof(size_t));

    for (i = args->start; i < args->end; i++) {
        idx = args->order[i];
        len = args->str_lengths[idx];
        if (args->char_pos < len) {
            byte_val = (unsigned char)args->strings[idx][args->char_pos];
        } else {
            byte_val = 0;
        }
        args->local_counts[byte_val]++;
    }

    return NULL;
}

/* Thread function: scatter strings */
static void *string_scatter_thread(void *arg)
{
    string_scatter_args_t *args = (string_scatter_args_t *)arg;
    size_t i;
    unsigned char byte_val;
    size_t idx, len;
    size_t local_offsets[RADIX_SIZE + 1];

    memcpy(local_offsets, args->local_offsets, (RADIX_SIZE + 1) * sizeof(size_t));

    for (i = args->start; i < args->end; i++) {
        idx = args->order[i];
        len = args->str_lengths[idx];
        if (args->char_pos < len) {
            byte_val = (unsigned char)args->strings[idx][args->char_pos];
        } else {
            byte_val = 0;
        }
        args->temp_order[local_offsets[byte_val]++] = idx;
    }

    return NULL;
}

/*
    String sort context for reusable allocations.
*/
typedef struct {
    int num_threads;
    pthread_t *threads;
    string_histogram_args_t *hist_args;
    string_scatter_args_t *scatter_args;
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
    ctx->threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
    ctx->hist_args = (string_histogram_args_t *)malloc(num_threads * sizeof(string_histogram_args_t));
    ctx->scatter_args = (string_scatter_args_t *)malloc(num_threads * sizeof(string_scatter_args_t));
    ctx->all_local_counts = (size_t **)malloc(num_threads * sizeof(size_t *));
    ctx->global_counts = (size_t *)calloc(RADIX_SIZE + 1, sizeof(size_t));
    ctx->global_offsets = (size_t *)malloc((RADIX_SIZE + 1) * sizeof(size_t));
    ctx->thread_offsets = (size_t **)malloc(num_threads * sizeof(size_t *));

    if (!ctx->threads || !ctx->hist_args || !ctx->scatter_args ||
        !ctx->all_local_counts || !ctx->global_counts ||
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
        free(ctx->threads);
        free(ctx->hist_args);
        free(ctx->scatter_args);
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
    free(ctx->threads);
    free(ctx->hist_args);
    free(ctx->scatter_args);
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
static int radix_sort_pass_string_parallel(size_t *order,
                                           size_t *temp_order,
                                           char **strings,
                                           size_t *str_lengths,
                                           size_t nobs,
                                           size_t char_pos,
                                           string_sort_context_t *ctx)
{
    size_t chunk_size, start, end;
    int t, b;
    int num_threads = ctx->num_threads;

    chunk_size = (nobs + num_threads - 1) / num_threads;

    /* Phase 1: Parallel histogram */
    for (t = 0; t < num_threads; t++) {
        start = (size_t)t * chunk_size;
        end = (size_t)(t + 1) * chunk_size;
        if (end > nobs) end = nobs;
        if (start >= nobs) start = end = nobs;

        ctx->hist_args[t].order = order;
        ctx->hist_args[t].strings = strings;
        ctx->hist_args[t].str_lengths = str_lengths;
        ctx->hist_args[t].start = start;
        ctx->hist_args[t].end = end;
        ctx->hist_args[t].char_pos = char_pos;
        ctx->hist_args[t].local_counts = ctx->all_local_counts[t];

        pthread_create(&ctx->threads[t], NULL, string_histogram_thread, &ctx->hist_args[t]);
    }

    for (t = 0; t < num_threads; t++) {
        pthread_join(ctx->threads[t], NULL);
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

    /* Phase 2: Parallel scatter */
    for (t = 0; t < num_threads; t++) {
        start = (size_t)t * chunk_size;
        end = (size_t)(t + 1) * chunk_size;
        if (end > nobs) end = nobs;
        if (start >= nobs) start = end = nobs;

        ctx->scatter_args[t].order = order;
        ctx->scatter_args[t].temp_order = temp_order;
        ctx->scatter_args[t].strings = strings;
        ctx->scatter_args[t].str_lengths = str_lengths;
        ctx->scatter_args[t].start = start;
        ctx->scatter_args[t].end = end;
        ctx->scatter_args[t].char_pos = char_pos;
        ctx->scatter_args[t].local_offsets = ctx->thread_offsets[t];

        pthread_create(&ctx->threads[t], NULL, string_scatter_thread, &ctx->scatter_args[t]);
    }

    for (t = 0; t < num_threads; t++) {
        pthread_join(ctx->threads[t], NULL);
    }

    return 0;
}

/*
    Sequential string radix sort pass with pre-cached lengths.
    Returns 1 if skipped, 0 otherwise.
*/
static int radix_sort_pass_string(size_t *order,
                                  size_t *temp_order,
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
        if (char_pos < len) {
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
        if (char_pos < len) {
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
    size_t *order_a;
    size_t *order_b;
    size_t *current_order;
    size_t *temp_order;
    uint64_t *keys;
    size_t i;
    int byte_pos;
    double *dbl_data;
    int use_parallel;
    int num_threads;
    radix_sort_context_t *ctx = NULL;
    int swapped = 0;

    /* Allocate temporary order array */
    order_a = data->sort_order;
    order_b = (size_t *)aligned_alloc_wrapper(CACHE_LINE_SIZE, data->nobs * sizeof(size_t));

    /* Allocate aligned keys array */
    keys = (uint64_t *)aligned_alloc_wrapper(CACHE_LINE_SIZE, data->nobs * sizeof(uint64_t));

    if (order_b == NULL || keys == NULL) {
        free(order_b);
        free(keys);
        return STATA_ERR_MEMORY;
    }

    /* Convert doubles to sortable uint64 keys */
    dbl_data = data->vars[var_idx].data.dbl;
    for (i = 0; i < data->nobs; i++) {
        keys[i] = double_to_sortable_uint64(dbl_data[i]);
    }

    /* Decide whether to use parallel sort */
    use_parallel = (data->nobs >= MIN_OBS_PER_THREAD * 2);
    num_threads = NUM_THREADS;
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
            free(order_b);
            free(keys);
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
            size_t *tmp = current_order;
            current_order = temp_order;
            temp_order = tmp;
            swapped = !swapped;
        }
    }

    /* Ensure final result is in data->sort_order */
    if (current_order != data->sort_order) {
        memcpy(data->sort_order, current_order, data->nobs * sizeof(size_t));
    }

    radix_context_free(ctx);
    free(order_b);
    free(keys);
    return STATA_OK;
}

/*
    Sort by a single string variable using radix sort.
    Optimization 2: Pre-caches string lengths.
    Optimization 5: Uses parallel sort for large datasets.
*/
static stata_retcode sort_by_string_var(stata_data *data, int var_idx)
{
    size_t *order_a;
    size_t *order_b;
    size_t *current_order;
    size_t *temp_order;
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

    /* Optimization 2: Pre-cache all string lengths */
    str_lengths = (size_t *)malloc(data->nobs * sizeof(size_t));
    if (str_lengths == NULL) {
        return STATA_ERR_MEMORY;
    }

    for (i = 0; i < data->nobs; i++) {
        len = strlen(str_data[i]);
        str_lengths[i] = len;
        if (len > max_len) {
            max_len = len;
        }
    }

    if (max_len == 0) {
        free(str_lengths);
        return STATA_OK;
    }

    /* Allocate temporary order array - aligned for cache efficiency */
    order_a = data->sort_order;
    order_b = (size_t *)aligned_alloc_wrapper(CACHE_LINE_SIZE, data->nobs * sizeof(size_t));
    if (order_b == NULL) {
        free(str_lengths);
        return STATA_ERR_MEMORY;
    }

    /* Decide on parallelization */
    use_parallel = (data->nobs >= MIN_OBS_PER_THREAD * 2);
    num_threads = NUM_THREADS;
    if (data->nobs < (size_t)MIN_OBS_PER_THREAD * (size_t)num_threads) {
        num_threads = (int)(data->nobs / MIN_OBS_PER_THREAD);
        if (num_threads < 2) {
            use_parallel = 0;
        }
    }

    if (use_parallel) {
        ctx = string_context_alloc(num_threads);
        if (ctx == NULL) {
            free(order_b);
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
            size_t *tmp = current_order;
            current_order = temp_order;
            temp_order = tmp;
            swapped = !swapped;
        }
    }

    /* Ensure final result is in data->sort_order */
    if (current_order != data->sort_order) {
        memcpy(data->sort_order, current_order, data->nobs * sizeof(size_t));
    }

    string_context_free(ctx);
    free(order_b);
    free(str_lengths);
    return STATA_OK;
}

/*
    Thread argument structure for parallel permutation application.
    Each thread receives one of these to permute a single variable.
*/
typedef struct {
    stata_variable *var;     /* [in/out] Variable to permute in place */
    size_t *sort_order;      /* [in]     Permutation: new[i] = old[sort_order[i]] */
    size_t nobs;             /* [in]     Number of observations */
    int success;             /* [out]    1 on success, 0 on failure */
} apply_permute_args_t;

/*
    Thread function: Apply permutation to a single variable's data.

    Allocates new array, copies data in permuted order, replaces original.
    For strings, only pointer array is reallocated (strings themselves stay).

    @param arg  Pointer to apply_permute_args_t with input/output parameters
    @return     NULL (success/failure indicated in args->success)
*/
static void *apply_permute_thread(void *arg)
{
    apply_permute_args_t *args = (apply_permute_args_t *)arg;
    size_t i;
    size_t nobs = args->nobs;
    size_t *perm = args->sort_order;

    if (args->var->type == STATA_TYPE_DOUBLE) {
        double *old_data = args->var->data.dbl;
        double *new_data = (double *)aligned_alloc_wrapper(CACHE_LINE_SIZE, nobs * sizeof(double));
        if (new_data == NULL) {
            args->success = 0;
            return NULL;
        }

        /* Apply permutation: new_data[i] = old_data[perm[i]] */
        for (i = 0; i < nobs; i++) {
            new_data[i] = old_data[perm[i]];
        }

        /* Replace old data with new */
        free(old_data);
        args->var->data.dbl = new_data;
    } else {
        /* String variable */
        char **old_data = args->var->data.str;
        char **new_data = (char **)malloc(nobs * sizeof(char *));
        if (new_data == NULL) {
            args->success = 0;
            return NULL;
        }

        /* Apply permutation: new_data[i] = old_data[perm[i]] */
        for (i = 0; i < nobs; i++) {
            new_data[i] = old_data[perm[i]];
        }

        /* Replace old data with new (don't free strings, just the array) */
        free(old_data);
        args->var->data.str = new_data;
    }

    args->success = 1;
    return NULL;
}

/*
    Apply sort_order permutation to all variables in parallel.

    After this function returns, data->vars[j].data arrays contain
    physically sorted data, enabling sequential writes in the store phase.

    @param data  [in/out] stata_data with computed sort_order:
                 - Input: data->sort_order contains permutation
                 - Output: data->vars[j].data physically reordered
                 - Output: data->sort_order reset to identity [0,1,2,...]

    @return STATA_OK on success, STATA_ERR_MEMORY on allocation failure
*/
static stata_retcode apply_permutation_to_all_vars(stata_data *data)
{
    size_t j;
    size_t nvars = data->nvars;
    pthread_t *threads;
    apply_permute_args_t *args;
    int all_success = 1;

    if (nvars == 0) {
        return STATA_OK;
    }

    threads = (pthread_t *)malloc(nvars * sizeof(pthread_t));
    args = (apply_permute_args_t *)malloc(nvars * sizeof(apply_permute_args_t));

    if (threads == NULL || args == NULL) {
        free(threads);
        free(args);
        return STATA_ERR_MEMORY;
    }

    /* Launch threads to apply permutation to each variable */
    for (j = 0; j < nvars; j++) {
        args[j].var = &data->vars[j];
        args[j].sort_order = data->sort_order;
        args[j].nobs = data->nobs;
        args[j].success = 0;

        pthread_create(&threads[j], NULL, apply_permute_thread, &args[j]);
    }

    /* Wait for all threads */
    for (j = 0; j < nvars; j++) {
        pthread_join(threads[j], NULL);
        if (!args[j].success) {
            all_success = 0;
        }
    }

    free(threads);
    free(args);

    /* Reset sort_order to identity since data is now physically sorted */
    for (j = 0; j < data->nobs; j++) {
        data->sort_order[j] = j;
    }

    return all_success ? STATA_OK : STATA_ERR_MEMORY;
}

/*
    Sort data using LSD radix sort and apply permutation to all variables.

    This is the main entry point for the sort module. It performs a stable
    multi-key sort using LSD radix sort, then physically reorders all variable
    data so the store phase can write sequentially.

    @param data       [in/out] stata_data structure:
                      - Input: data->vars contains unsorted data
                      - Input: data->sort_order is identity permutation
                      - Output: data->vars contains physically sorted data
                      - Output: data->sort_order reset to identity
    @param sort_vars  [in] Array of 1-based variable indices specifying sort keys
    @param nsort      [in] Number of sort key variables

    @return STATA_OK on success, or:
            STATA_ERR_INVALID_INPUT if data/sort_vars is NULL, nobs=0, or nsort=0
            STATA_ERR_MEMORY on allocation failure

    Multi-Key Sort Order:
    - Keys processed from last (least significant) to first (most significant)
    - This ensures stable sort: ties in key[0] broken by key[1], etc.
    - Example: sort by (name, age) → first sort by age, then by name
*/
stata_retcode ctools_sort_radix_lsd(stata_data *data, int *sort_vars, size_t nsort)
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
        This ensures that ties in more significant keys are broken by less significant keys.
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

    /* Apply the permutation to all variables so data is physically sorted */
    rc = apply_permutation_to_all_vars(data);
    if (rc != STATA_OK) {
        return rc;
    }

    return STATA_OK;
}

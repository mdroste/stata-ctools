/*
    ctools_sort_radix_msd.c
    Parallel MSD (Most Significant Digit) Radix Sort Module

    MSD radix sort is optimized for variable-length strings where it can
    short-circuit on unique prefixes. Unlike LSD which must process all
    characters, MSD stops as soon as elements are distinguishable.

    Key advantages over LSD for strings:
    - Early termination: stops when all elements in a bucket are identical
    - Better cache locality: strings with common prefixes stay together
    - O(N*k) where k is the average distinguishing prefix length (not max length)

    Algorithm:
    - Process characters from most significant (leftmost) to least
    - Recursively sort each non-empty bucket
    - Use insertion sort for small buckets (better cache behavior)

    Parallelization Strategy:
    - Top-level histogram and scatter are parallelized with pthreads
    - Recursive calls on independent buckets run in parallel via OpenMP tasks
    - Falls back to sequential insertion sort for small subarrays
    - Reusable context allocation to avoid per-pass malloc overhead

    Numeric Sorting:
    - Uses same IEEE 754 bit manipulation as LSD
    - Processes bytes from MSB (byte 7) to LSB (byte 0)
    - Single pass when data has low entropy in high bytes

    Optimizations:
    - Pointer swapping instead of memcpy between passes
    - Parallel histogram and scatter at top level
    - Reusable allocation context across all passes
    - Cache-aligned memory allocation
*/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_threads.h"

/* Threshold for switching to insertion sort */
#define MSD_INSERTION_THRESHOLD 64

/* Minimum bucket size for parallel recursion */
#define MSD_PARALLEL_THRESHOLD 10000

/* Minimum observations for parallel histogram/scatter */
#define MSD_PARALLEL_HIST_THRESHOLD (MIN_OBS_PER_THREAD * 2)

/* Maximum recursion depth before switching to a non-recursive fallback
   (strings: merge sort; numeric keys have 8 bytes and never get this deep).
   Each recursive call uses ~4KB stack space, limit to prevent overflow */
#define MSD_MAX_RECURSION_DEPTH 32

/* ============================================================================
   Utility Functions
   ============================================================================ */

/*
    NOTE: Aligned memory allocation is now provided by ctools_config.h
    Use ctools_aligned_alloc() and ctools_aligned_free() for cross-platform support.

    NOTE: Double-to-sortable conversion is now provided by ctools_types.h
    Use ctools_double_to_sortable(d) for sorting doubles.
*/

/* ============================================================================
   Reusable Context for Parallel MSD Sort
   ============================================================================ */

/* Thread argument for parallel histogram (uses perm_idx_t for 50% memory savings) */
typedef struct {
    perm_idx_t *order;
    uint64_t *keys;
    size_t start;
    size_t end;
    int shift;
    size_t *local_counts;
} msd_histogram_args_t;

/* Thread argument for parallel scatter (uses perm_idx_t for 50% memory savings) */
typedef struct {
    perm_idx_t *order;
    perm_idx_t *temp_order;
    uint64_t *keys;
    size_t start;
    size_t end;
    int shift;
    size_t *local_offsets;
} msd_scatter_args_t;

/*
    Reusable allocation context for parallel MSD sort.
    Allocated once per variable and reused across all byte passes.
*/
typedef struct {
    int num_threads;
    pthread_t *threads;
    msd_histogram_args_t *hist_args;
    msd_scatter_args_t *scatter_args;
    size_t **all_local_counts;   /* [num_threads][RADIX_SIZE] */
    size_t *global_counts;       /* [RADIX_SIZE] */
    size_t *global_offsets;      /* [RADIX_SIZE] */
    size_t **thread_offsets;     /* [num_threads][RADIX_SIZE] */
} msd_sort_context_t;

/*
    Allocate reusable context for parallel MSD sort.
*/
static msd_sort_context_t *msd_context_alloc(int num_threads)
{
    msd_sort_context_t *ctx;
    int t;

    ctx = (msd_sort_context_t *)calloc(1, sizeof(msd_sort_context_t));
    if (ctx == NULL) return NULL;

    ctx->num_threads = num_threads;
    ctx->threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
    ctx->hist_args = (msd_histogram_args_t *)malloc(num_threads * sizeof(msd_histogram_args_t));
    ctx->scatter_args = (msd_scatter_args_t *)malloc(num_threads * sizeof(msd_scatter_args_t));
    ctx->all_local_counts = (size_t **)calloc(num_threads, sizeof(size_t *));
    ctx->global_counts = (size_t *)calloc(RADIX_SIZE, sizeof(size_t));
    ctx->global_offsets = (size_t *)malloc(RADIX_SIZE * sizeof(size_t));
    ctx->thread_offsets = (size_t **)calloc(num_threads, sizeof(size_t *));

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
    Free MSD sort context.
*/
static void msd_context_free(msd_sort_context_t *ctx)
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

/* Thread function: compute local histogram for MSD */
static void *msd_histogram_thread(void *arg)
{
    msd_histogram_args_t *args = (msd_histogram_args_t *)arg;
    size_t i;
    uint8_t byte_val;
    perm_idx_t *order = args->order;
    uint64_t *keys = args->keys;
    int shift = args->shift;

    memset(args->local_counts, 0, RADIX_SIZE * sizeof(size_t));

    for (i = args->start; i < args->end; i++) {
        byte_val = (keys[order[i]] >> shift) & RADIX_MASK;
        args->local_counts[byte_val]++;
    }

    return NULL;
}

/* Thread function: scatter elements to sorted positions */
static void *msd_scatter_thread(void *arg)
{
    msd_scatter_args_t *args = (msd_scatter_args_t *)arg;
    size_t i;
    uint8_t byte_val;
    size_t local_offsets[RADIX_SIZE];
    perm_idx_t *order = args->order;
    uint64_t *keys = args->keys;
    perm_idx_t *temp_order = args->temp_order;
    int shift = args->shift;

    memcpy(local_offsets, args->local_offsets, RADIX_SIZE * sizeof(size_t));

    for (i = args->start; i < args->end; i++) {
        byte_val = (keys[order[i]] >> shift) & RADIX_MASK;
        temp_order[local_offsets[byte_val]++] = order[i];
    }

    return NULL;
}

/* ============================================================================
   Insertion Sort (for small subarrays)
   ============================================================================ */

/*
    Insertion sort for numeric data using precomputed keys.
    Used when subarray size falls below MSD_INSERTION_THRESHOLD.
*/
static void insertion_sort_numeric(perm_idx_t *order, uint64_t *keys,
                                   size_t start, size_t end)
{
    size_t i, j;
    perm_idx_t temp_idx;
    uint64_t temp_key;

    for (i = start + 1; i < end; i++) {
        temp_idx = order[i];
        temp_key = keys[temp_idx];
        j = i;

        while (j > start && keys[order[j - 1]] > temp_key) {
            order[j] = order[j - 1];
            j--;
        }
        order[j] = temp_idx;
    }
}

/*
    Compare two strings from char_pos onwards (earlier bytes are shared).
    A string that ends first sorts first, as in the radix passes.
*/
static inline int compare_string_from(char **strings, size_t *str_lengths,
                                      perm_idx_t a, perm_idx_t b, size_t char_pos)
{
    size_t len1 = str_lengths[a] > char_pos ? str_lengths[a] - char_pos : 0;
    size_t len2 = str_lengths[b] > char_pos ? str_lengths[b] - char_pos : 0;
    size_t min_len = len1 < len2 ? len1 : len2;
    int cmp = min_len ? memcmp(strings[a] + char_pos, strings[b] + char_pos, min_len) : 0;

    if (cmp == 0) {
        cmp = (len1 > len2) - (len1 < len2);
    }
    return cmp;
}

/*
    Insertion sort for string data.
    Compares strings starting from a given character position.
*/
static void insertion_sort_string(perm_idx_t *order, char **strings,
                                  size_t *str_lengths, size_t start,
                                  size_t end, size_t char_pos)
{
    size_t i, j;
    perm_idx_t temp_idx;

    for (i = start + 1; i < end; i++) {
        temp_idx = order[i];
        j = i;

        while (j > start &&
               compare_string_from(strings, str_lengths, order[j - 1], temp_idx, char_pos) > 0) {
            order[j] = order[j - 1];
            j--;
        }
        order[j] = temp_idx;
    }
}

/*
    First position at or after char_pos where the strings in order[start, end)
    stop sharing a byte (two differ or one ends). All strings must be at least
    char_pos bytes long. Lets MSD jump over a long common prefix at once
    instead of making one radix pass per shared byte.
*/
static size_t common_prefix_end(perm_idx_t *order, char **strings, size_t *str_lengths,
                                size_t start, size_t end, size_t char_pos)
{
    const char *first = strings[order[start]];
    size_t lcp_end = str_lengths[order[start]];
    size_t i;

    for (i = start + 1; i < end && lcp_end > char_pos; i++) {
        const char *s = strings[order[i]];
        size_t limit = str_lengths[order[i]] < lcp_end ? str_lengths[order[i]] : lcp_end;
        size_t p = char_pos;

        while (p < limit && s[p] == first[p]) p++;
        lcp_end = p;
    }
    return lcp_end > char_pos ? lcp_end : char_pos;
}

/*
    Stable bottom-up merge sort for string data, comparing from char_pos.
    O(n log n) fallback once MSD recursion gets too deep. Sorts order[start, end)
    using temp[start, end) as scratch; returns the array holding the result.
*/
static perm_idx_t *merge_sort_string(perm_idx_t *order, perm_idx_t *temp,
                                     char **strings, size_t *str_lengths,
                                     size_t start, size_t end, size_t char_pos)
{
    perm_idx_t *src = order;
    perm_idx_t *dst = temp;
    size_t width, lo;

    for (lo = start; lo < end; lo += MSD_INSERTION_THRESHOLD) {
        size_t hi = (end - lo > MSD_INSERTION_THRESHOLD) ? lo + MSD_INSERTION_THRESHOLD : end;
        insertion_sort_string(src, strings, str_lengths, lo, hi, char_pos);
    }

    for (width = MSD_INSERTION_THRESHOLD; width < end - start; width *= 2) {
        for (lo = start; lo < end; lo += 2 * width) {
            size_t mid = (end - lo > width) ? lo + width : end;
            size_t hi = (end - mid > width) ? mid + width : end;
            size_t i = lo, j = mid, k = lo;

            /* Take from the right run only when strictly smaller (stable) */
            while (i < mid && j < hi) {
                if (compare_string_from(strings, str_lengths, src[j], src[i], char_pos) < 0) {
                    dst[k++] = src[j++];
                } else {
                    dst[k++] = src[i++];
                }
            }
            while (i < mid) dst[k++] = src[i++];
            while (j < hi) dst[k++] = src[j++];
        }
        perm_idx_t *t = src; src = dst; dst = t;
    }

    return src;
}

/* ============================================================================
   MSD Radix Sort for Numeric Data
   ============================================================================ */

/*
    Parallel MSD radix sort pass for numeric data.
    Uses context for thread-local allocations.
    Returns 1 if pass was skipped (uniform distribution), 0 otherwise.
*/
static int msd_sort_pass_numeric_parallel(perm_idx_t *order, perm_idx_t *temp_order,
                                          uint64_t *keys, size_t start,
                                          size_t end, int byte_pos,
                                          msd_sort_context_t *ctx,
                                          size_t *bucket_starts_out,
                                          size_t *counts_out)
{
    size_t n = end - start;
    size_t chunk_size;
    int shift = byte_pos * RADIX_BITS;
    int t, b;
    int num_threads = ctx->num_threads;

    /* Phase 1: Parallel histogram */
    chunk_size = (n + num_threads - 1) / num_threads;
    int threads_created = 0;
    for (t = 0; t < num_threads; t++) {
        size_t t_start = start + (size_t)t * chunk_size;
        size_t t_end = start + (size_t)(t + 1) * chunk_size;
        if (t_end > end) t_end = end;
        if (t_start >= end) t_start = t_end = end;

        ctx->hist_args[t].order = order;
        ctx->hist_args[t].keys = keys;
        ctx->hist_args[t].start = t_start;
        ctx->hist_args[t].end = t_end;
        ctx->hist_args[t].shift = shift;
        ctx->hist_args[t].local_counts = ctx->all_local_counts[t];

        if (pthread_create(&ctx->threads[threads_created], NULL, msd_histogram_thread, &ctx->hist_args[t]) == 0) {
            threads_created++;
        } else {
            /* Run in current thread as fallback */
            msd_histogram_thread(&ctx->hist_args[t]);
        }
    }

    for (t = 0; t < threads_created; t++) {
        pthread_join(ctx->threads[t], NULL);
    }

    /* Combine histograms - iterate bucket-first for cache locality */
    memset(ctx->global_counts, 0, RADIX_SIZE * sizeof(size_t));
    for (b = 0; b < RADIX_SIZE; b++) {
        for (t = 0; t < num_threads; t++) {
            ctx->global_counts[b] += ctx->all_local_counts[t][b];
        }
    }

    /* Check for uniform distribution */
    {
        int non_empty = 0;
        for (b = 0; b < RADIX_SIZE; b++) {
            if (ctx->global_counts[b] > 0) {
                non_empty++;
                if (non_empty > 1) break;
            }
        }
        if (non_empty <= 1) {
            return 1;  /* Skip this pass */
        }
    }

    /* Compute prefix sums */
    ctx->global_offsets[0] = start;
    for (b = 1; b < RADIX_SIZE; b++) {
        ctx->global_offsets[b] = ctx->global_offsets[b - 1] + ctx->global_counts[b - 1];
    }

    /* Save bucket starts for recursive calls */
    if (bucket_starts_out) {
        for (b = 0; b < RADIX_SIZE; b++) {
            bucket_starts_out[b] = ctx->global_offsets[b];
        }
    }
    if (counts_out) {
        memcpy(counts_out, ctx->global_counts, RADIX_SIZE * sizeof(size_t));
    }

    /* Compute per-thread offsets */
    for (b = 0; b < RADIX_SIZE; b++) {
        size_t offset = ctx->global_offsets[b];
        for (t = 0; t < num_threads; t++) {
            ctx->thread_offsets[t][b] = offset;
            offset += ctx->all_local_counts[t][b];
        }
    }

    /* Phase 2: Parallel scatter */
    threads_created = 0;
    for (t = 0; t < num_threads; t++) {
        size_t t_start = start + (size_t)t * chunk_size;
        size_t t_end = start + (size_t)(t + 1) * chunk_size;
        if (t_end > end) t_end = end;
        if (t_start >= end) t_start = t_end = end;

        ctx->scatter_args[t].order = order;
        ctx->scatter_args[t].temp_order = temp_order;
        ctx->scatter_args[t].keys = keys;
        ctx->scatter_args[t].start = t_start;
        ctx->scatter_args[t].end = t_end;
        ctx->scatter_args[t].shift = shift;
        ctx->scatter_args[t].local_offsets = ctx->thread_offsets[t];

        if (pthread_create(&ctx->threads[threads_created], NULL, msd_scatter_thread, &ctx->scatter_args[t]) == 0) {
            threads_created++;
        } else {
            /* Run in current thread as fallback */
            msd_scatter_thread(&ctx->scatter_args[t]);
        }
    }

    for (t = 0; t < threads_created; t++) {
        pthread_join(ctx->threads[t], NULL);
    }

    return 0;
}

/*
    Sequential MSD radix sort pass (for small subarrays or recursive calls).
    Returns 1 if skipped (uniform), 0 otherwise.
*/
static int msd_sort_pass_numeric_seq(perm_idx_t *order, perm_idx_t *temp_order,
                                     uint64_t *keys, size_t start, size_t end,
                                     int byte_pos, size_t *bucket_starts_out,
                                     size_t *counts_out)
{
    size_t counts[RADIX_SIZE] = {0};
    size_t offsets[RADIX_SIZE];
    size_t i;
    uint8_t byte_val;
    int shift = byte_pos * RADIX_BITS;
    int non_empty = 0;
    int b;

    for (i = start; i < end; i++) {
        byte_val = (keys[order[i]] >> shift) & RADIX_MASK;
        counts[byte_val]++;
    }

    for (b = 0; b < RADIX_SIZE; b++) {
        if (counts[b] > 0) {
            non_empty++;
            if (non_empty > 1) break;
        }
    }
    if (non_empty <= 1) {
        return 1;
    }

    offsets[0] = start;
    for (b = 1; b < RADIX_SIZE; b++) {
        offsets[b] = offsets[b - 1] + counts[b - 1];
    }

    if (bucket_starts_out) {
        for (b = 0; b < RADIX_SIZE; b++) {
            bucket_starts_out[b] = (b == 0) ? start : offsets[b];
        }
    }
    if (counts_out) {
        memcpy(counts_out, counts, RADIX_SIZE * sizeof(size_t));
    }

    for (i = start; i < end; i++) {
        byte_val = (keys[order[i]] >> shift) & RADIX_MASK;
        temp_order[offsets[byte_val]++] = order[i];
    }

    return 0;
}

/*
    Recursive MSD radix sort for numeric data (for subarrays).
    Uses pointer swapping to avoid memcpy.
    order_is_primary indicates whether order currently holds the data.
    depth tracks recursion depth to prevent stack overflow.
*/
static void msd_sort_numeric_recursive_impl(perm_idx_t *order_a, perm_idx_t *order_b,
                                            uint64_t *keys, size_t start,
                                            size_t end, int byte_pos,
                                            int a_is_current, int depth)
{
    size_t counts[RADIX_SIZE];
    size_t bucket_starts[RADIX_SIZE];
    size_t n = end - start;
    perm_idx_t *current_order = a_is_current ? order_a : order_b;
    perm_idx_t *temp_order = a_is_current ? order_b : order_a;
    int b;
    int skipped;

    /* Base case: use insertion sort for small arrays */
    if (n <= MSD_INSERTION_THRESHOLD) {
        insertion_sort_numeric(current_order, keys, start, end);
        /* If result is in wrong buffer, copy it */
        if (!a_is_current) {
            memcpy(&order_a[start], &order_b[start], n * sizeof(perm_idx_t));
        }
        return;
    }

    /* Base case: processed all bytes */
    if (byte_pos < 0) {
        if (!a_is_current) {
            memcpy(&order_a[start], &order_b[start], n * sizeof(perm_idx_t));
        }
        return;
    }

    /* Fallback: if recursion too deep, use insertion sort to prevent stack overflow */
    if (depth >= MSD_MAX_RECURSION_DEPTH) {
        insertion_sort_numeric(current_order, keys, start, end);
        if (!a_is_current) {
            memcpy(&order_a[start], &order_b[start], n * sizeof(perm_idx_t));
        }
        return;
    }

    skipped = msd_sort_pass_numeric_seq(current_order, temp_order, keys,
                                        start, end, byte_pos,
                                        bucket_starts, counts);

    if (skipped) {
        /* Skip this byte, continue with same buffer */
        msd_sort_numeric_recursive_impl(order_a, order_b, keys, start, end,
                                        byte_pos - 1, a_is_current, depth + 1);
        return;
    }

    /* Result is now in temp_order, so flip the flag */
    int next_is_a = !a_is_current;

    /* Recursively sort each bucket */
    #ifdef _OPENMP
    #pragma omp parallel if(n >= MSD_PARALLEL_THRESHOLD)
    {
        #pragma omp single nowait
        {
    #endif
            for (b = 0; b < RADIX_SIZE; b++) {
                if (counts[b] > 1) {
                    size_t bucket_start = bucket_starts[b];
                    size_t bucket_count = counts[b];
                    size_t bucket_end = bucket_start + bucket_count;
    #ifdef _OPENMP
                    if (bucket_count >= MSD_PARALLEL_THRESHOLD) {
                        #pragma omp task firstprivate(bucket_start, bucket_end, next_is_a, byte_pos, depth)
                        msd_sort_numeric_recursive_impl(order_a, order_b, keys,
                                                        bucket_start, bucket_end,
                                                        byte_pos - 1, next_is_a, depth + 1);
                    } else {
                        msd_sort_numeric_recursive_impl(order_a, order_b, keys,
                                                        bucket_start, bucket_end,
                                                        byte_pos - 1, next_is_a, depth + 1);
                    }
    #else
                    msd_sort_numeric_recursive_impl(order_a, order_b, keys,
                                                    bucket_start, bucket_end,
                                                    byte_pos - 1, next_is_a, depth + 1);
    #endif
                } else if (counts[b] == 1 && !next_is_a) {
                    /* Single element in wrong buffer - copy it */
                    order_a[bucket_starts[b]] = order_b[bucket_starts[b]];
                }
            }
    #ifdef _OPENMP
        }
        #pragma omp taskwait
    }
    #endif
}

/* Wrapper for backward compatibility */
static void msd_sort_numeric_recursive(perm_idx_t *order_a, perm_idx_t *order_b,
                                       uint64_t *keys, size_t start,
                                       size_t end, int byte_pos,
                                       int a_is_current)
{
    msd_sort_numeric_recursive_impl(order_a, order_b, keys, start, end,
                                    byte_pos, a_is_current, 0);
}

/*
    Sort by a single numeric variable using MSD radix sort.
    Uses parallel histogram/scatter at top level, then recursive for buckets.
*/
static stata_retcode msd_sort_by_numeric_var(stata_data *data, int var_idx)
{
    perm_idx_t *order_a;
    perm_idx_t *order_b;
    uint64_t *keys;
    size_t i;
    double *dbl_data;
    int use_parallel;
    int num_threads;
    msd_sort_context_t *ctx = NULL;
    size_t bucket_starts[RADIX_SIZE];
    size_t counts[RADIX_SIZE];
    int skipped;
    int b;

    order_a = data->sort_order;
    order_b = (perm_idx_t *)ctools_aligned_alloc(CACHE_LINE_SIZE,
                                           data->nobs * sizeof(perm_idx_t));
    keys = (uint64_t *)ctools_aligned_alloc(CACHE_LINE_SIZE,
                                          data->nobs * sizeof(uint64_t));

    if (order_b == NULL || keys == NULL) {
        ctools_aligned_free(order_b);
        ctools_aligned_free(keys);
        return STATA_ERR_MEMORY;
    }

    /* Convert doubles to sortable uint64 keys */
    dbl_data = data->vars[var_idx].data.dbl;
    #ifdef _OPENMP
    #pragma omp parallel for if(data->nobs >= MSD_PARALLEL_HIST_THRESHOLD)
    #endif
    for (i = 0; i < data->nobs; i++) {
        keys[i] = ctools_double_to_sortable(dbl_data[i]);
    }

    /* Decide on parallel sort */
    use_parallel = (data->nobs >= MSD_PARALLEL_HIST_THRESHOLD);
    num_threads = ctools_get_openmp_threads();
    if (data->nobs < (size_t)MIN_OBS_PER_THREAD * (size_t)num_threads) {
        num_threads = (int)(data->nobs / MIN_OBS_PER_THREAD);
        if (num_threads < 2) {
            use_parallel = 0;
        }
    }

    if (use_parallel) {
        ctx = msd_context_alloc(num_threads);
        if (ctx == NULL) {
            ctools_aligned_free(order_b);
            ctools_aligned_free(keys);
            return STATA_ERR_MEMORY;
        }

        /* Top-level parallel sort pass */
        skipped = msd_sort_pass_numeric_parallel(order_a, order_b, keys,
                                                  0, data->nobs, 7, ctx,
                                                  bucket_starts, counts);
        msd_context_free(ctx);

        if (skipped) {
            /* All in one bucket at byte 7, recurse sequentially */
            msd_sort_numeric_recursive(order_a, order_b, keys, 0, data->nobs, 6, 1);
        } else {
            /* Result is in order_b, recurse on each bucket */
            #ifdef _OPENMP
            #pragma omp parallel
            {
                #pragma omp single nowait
                {
            #endif
                    for (b = 0; b < RADIX_SIZE; b++) {
                        if (counts[b] > 1) {
                            size_t bucket_start = bucket_starts[b];
                            size_t bucket_count = counts[b];
                            size_t bucket_end = bucket_start + bucket_count;
            #ifdef _OPENMP
                            if (bucket_count >= MSD_PARALLEL_THRESHOLD) {
                                #pragma omp task firstprivate(bucket_start, bucket_end)
                                msd_sort_numeric_recursive_impl(order_a, order_b, keys,
                                                                bucket_start, bucket_end,
                                                                6, 0, 1);  /* 0 = data is in order_b, depth=1 */
                            } else {
                                msd_sort_numeric_recursive_impl(order_a, order_b, keys,
                                                                bucket_start, bucket_end,
                                                                6, 0, 1);
                            }
            #else
                            msd_sort_numeric_recursive_impl(order_a, order_b, keys,
                                                            bucket_start, bucket_end,
                                                            6, 0, 1);
            #endif
                        } else if (counts[b] == 1) {
                            /* Single element - copy from order_b to order_a */
                            order_a[bucket_starts[b]] = order_b[bucket_starts[b]];
                        }
                    }
            #ifdef _OPENMP
                }
                #pragma omp taskwait
            }
            #endif
        }
    } else {
        /* Sequential sort for small datasets */
        msd_sort_numeric_recursive(order_a, order_b, keys, 0, data->nobs, 7, 1);
    }

    ctools_aligned_free(order_b);
    ctools_aligned_free(keys);
    return STATA_OK;
}

/* ============================================================================
   MSD Radix Sort for String Data
   ============================================================================ */

/* String bucket count (0 for end-of-string, 1-256 for characters) */
#define STRING_BUCKET_SIZE 257

/* Keeps the string pass out of line: inlined into the recursion, its two
   257-entry arrays would double every level's stack frame (and slow the sort) */
#if defined(__GNUC__) || defined(__clang__)
    #define MSD_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
    #define MSD_NOINLINE __declspec(noinline)
#else
    #define MSD_NOINLINE
#endif

/*
    Sequential MSD string sort pass.
    Returns 1 if skipped, 0 otherwise.
*/
static MSD_NOINLINE int msd_sort_pass_string_seq(perm_idx_t *order, perm_idx_t *temp_order,
                                                 char **strings, size_t *str_lengths,
                                                 size_t start, size_t end, size_t char_pos,
                                                 size_t *bucket_starts_out, size_t *counts_out)
{
    size_t counts[STRING_BUCKET_SIZE] = {0};
    size_t offsets[STRING_BUCKET_SIZE];
    size_t i;
    perm_idx_t idx;
    size_t len;
    unsigned int bucket;
    int non_empty = 0;
    int b;

    for (i = start; i < end; i++) {
        idx = order[i];
        len = str_lengths[idx];
        bucket = (char_pos >= len) ? 0 : ((unsigned char)strings[idx][char_pos] + 1);
        counts[bucket]++;
    }

    for (b = 0; b < STRING_BUCKET_SIZE; b++) {
        if (counts[b] > 0) {
            non_empty++;
            if (non_empty > 1) break;
        }
    }
    if (non_empty <= 1) {
        /* Check if all strings ended */
        if (counts[0] == (end - start)) {
            return 2;  /* All done, no more recursion needed */
        }
        return 1;  /* Skip this pass but continue */
    }

    offsets[0] = start;
    for (b = 1; b < STRING_BUCKET_SIZE; b++) {
        offsets[b] = offsets[b - 1] + counts[b - 1];
    }

    if (bucket_starts_out) {
        for (b = 0; b < STRING_BUCKET_SIZE; b++) {
            bucket_starts_out[b] = (b == 0) ? start : offsets[b];
        }
    }
    if (counts_out) {
        memcpy(counts_out, counts, STRING_BUCKET_SIZE * sizeof(size_t));
    }

    for (i = start; i < end; i++) {
        idx = order[i];
        len = str_lengths[idx];
        bucket = (char_pos >= len) ? 0 : ((unsigned char)strings[idx][char_pos] + 1);
        temp_order[offsets[bucket]++] = idx;
    }

    return 0;
}

/*
    Recursive MSD radix sort for string data with pointer swapping.
    a_is_current indicates whether order_a currently holds the data.
    depth tracks recursion depth to prevent stack overflow.
*/
static void msd_sort_string_recursive_impl(perm_idx_t *order_a, perm_idx_t *order_b,
                                           char **strings, size_t *str_lengths,
                                           size_t start, size_t end, size_t char_pos,
                                           int a_is_current, int depth)
{
    size_t counts[STRING_BUCKET_SIZE];
    size_t bucket_starts[STRING_BUCKET_SIZE];
    size_t n = end - start;
    perm_idx_t *current_order = a_is_current ? order_a : order_b;
    perm_idx_t *temp_order = a_is_current ? order_b : order_a;
    int b;
    int result;

    /* Base case: use insertion sort for small arrays */
    if (n <= MSD_INSERTION_THRESHOLD) {
        insertion_sort_string(current_order, strings, str_lengths, start, end, char_pos);
        if (!a_is_current) {
            memcpy(&order_a[start], &order_b[start], n * sizeof(perm_idx_t));
        }
        return;
    }

    /* Fallback: if recursion too deep, finish with an O(n log n) stable merge
       sort to prevent stack overflow */
    if (depth >= MSD_MAX_RECURSION_DEPTH) {
        perm_idx_t *sorted = merge_sort_string(current_order, temp_order, strings,
                                               str_lengths, start, end, char_pos);
        if (sorted != order_a) {
            memcpy(&order_a[start], &sorted[start], n * sizeof(perm_idx_t));
        }
        return;
    }

    /* When every string shares the byte, jump to the end of the common prefix
       instead of recursing: long shared prefixes neither recurse nor count
       toward the depth limit */
    for (;;) {
        result = msd_sort_pass_string_seq(current_order, temp_order, strings, str_lengths,
                                          start, end, char_pos, bucket_starts, counts);
        if (result != 1) break;
        char_pos = common_prefix_end(current_order, strings, str_lengths,
                                     start, end, char_pos + 1);
    }

    if (result == 2) {
        /* All strings ended */
        if (!a_is_current) {
            memcpy(&order_a[start], &order_b[start], n * sizeof(perm_idx_t));
        }
        return;
    }

    /* Result is in temp_order, so flip the flag */
    int next_is_a = !a_is_current;

    /* Recursively sort each bucket (skip bucket 0 - strings that ended) */
    #ifdef _OPENMP
    #pragma omp parallel if(n >= MSD_PARALLEL_THRESHOLD)
    {
        #pragma omp single nowait
        {
    #endif
            for (b = 1; b < STRING_BUCKET_SIZE; b++) {
                if (counts[b] > 1) {
                    size_t bucket_start = bucket_starts[b];
                    size_t bucket_count = counts[b];
                    size_t bucket_end = bucket_start + bucket_count;
    #ifdef _OPENMP
                    if (bucket_count >= MSD_PARALLEL_THRESHOLD) {
                        #pragma omp task firstprivate(bucket_start, bucket_end, char_pos, next_is_a, depth)
                        msd_sort_string_recursive_impl(order_a, order_b, strings, str_lengths,
                                                       bucket_start, bucket_end,
                                                       char_pos + 1, next_is_a, depth + 1);
                    } else {
                        msd_sort_string_recursive_impl(order_a, order_b, strings, str_lengths,
                                                       bucket_start, bucket_end,
                                                       char_pos + 1, next_is_a, depth + 1);
                    }
    #else
                    msd_sort_string_recursive_impl(order_a, order_b, strings, str_lengths,
                                                   bucket_start, bucket_end,
                                                   char_pos + 1, next_is_a, depth + 1);
    #endif
                } else if (counts[b] == 1 && !next_is_a) {
                    /* Single element in wrong buffer */
                    order_a[bucket_starts[b]] = order_b[bucket_starts[b]];
                }
            }
            /* Handle bucket 0 (ended strings) */
            if (counts[0] > 0 && !next_is_a) {
                for (size_t j = bucket_starts[0]; j < bucket_starts[0] + counts[0]; j++) {
                    order_a[j] = order_b[j];
                }
            }
    #ifdef _OPENMP
        }
        #pragma omp taskwait
    }
    #endif
}

/* Wrapper for backward compatibility */
static void msd_sort_string_recursive(perm_idx_t *order_a, perm_idx_t *order_b,
                                      char **strings, size_t *str_lengths,
                                      size_t start, size_t end, size_t char_pos,
                                      int a_is_current)
{
    msd_sort_string_recursive_impl(order_a, order_b, strings, str_lengths,
                                   start, end, char_pos, a_is_current, 0);
}

/*
    Sort by a single string variable using MSD radix sort.
*/
static stata_retcode msd_sort_by_string_var(stata_data *data, int var_idx)
{
    perm_idx_t *order_a;
    perm_idx_t *order_b;
    size_t *str_lengths;
    size_t i;
    char **str_data;

    str_data = data->vars[var_idx].data.str;
    order_a = data->sort_order;

    /* Pre-cache string lengths with parallel computation */
    str_lengths = (size_t *)malloc(data->nobs * sizeof(size_t));
    if (str_lengths == NULL) {
        return STATA_ERR_MEMORY;
    }

    #ifdef _OPENMP
    #pragma omp parallel for if(data->nobs >= MSD_PARALLEL_HIST_THRESHOLD)
    #endif
    for (i = 0; i < data->nobs; i++) {
        str_lengths[i] = (str_data[i] != NULL) ? strlen(str_data[i]) : 0;
    }

    /* Allocate temporary order array (uses perm_idx_t for 50% memory savings) */
    order_b = (perm_idx_t *)ctools_aligned_alloc(CACHE_LINE_SIZE,
                                           data->nobs * sizeof(perm_idx_t));
    if (order_b == NULL) {
        free(str_lengths);
        return STATA_ERR_MEMORY;
    }

    /* Sort starting from position 0 */
    msd_sort_string_recursive(order_a, order_b, str_data, str_lengths,
                              0, data->nobs, 0, 1);  /* 1 = data starts in order_a */

    ctools_aligned_free(order_b);
    free(str_lengths);
    return STATA_OK;
}

/* ============================================================================
   Public API
   ============================================================================ */

/*
    ctools_sort_radix_msd_order_only - MSD radix sort without applying permutation

    Computes sort_order but does NOT apply the permutation to data.
    After this call, data->sort_order contains the permutation but data is unchanged.
    Call ctools_apply_permutation() separately to apply the permutation.

    @param data       Dataset to sort
    @param sort_vars  Array of 1-based variable indices to sort by
    @param nsort      Number of sort variables

    @return STATA_OK on success, or error code
*/
stata_retcode ctools_sort_radix_msd_order_only(stata_data *data, int *sort_vars, size_t nsort)
{
    int k;
    int var_idx;
    stata_retcode rc;

    if (data == NULL || sort_vars == NULL || data->nobs == 0 || nsort == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    for (k = (int)nsort - 1; k >= 0; k--) {
        var_idx = sort_vars[k] - 1;

        if (var_idx < 0 || var_idx >= (int)data->nvars) {
            return STATA_ERR_INVALID_INPUT;
        }

        if (data->vars[var_idx].type == STATA_TYPE_DOUBLE) {
            rc = msd_sort_by_numeric_var(data, var_idx);
        } else {
            rc = msd_sort_by_string_var(data, var_idx);
        }

        if (rc != STATA_OK) {
            return rc;
        }
    }

    /* Do NOT apply permutation - caller will do it separately */
    return STATA_OK;
}

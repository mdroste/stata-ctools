/*
    ctools_sort_merge.c
    Optimized Parallel Merge Sort Module

    A stable, predictable O(n log n) parallel sorting algorithm.

    Algorithm Overview:
    Phase 1 - Parallel Block Sort:
        - Divide data into p blocks of n/p elements
        - Each thread sorts its block independently using radix sort
        - All blocks can be sorted in parallel

    Phase 2 - Parallel K-way Merge:
        - Iteratively merge pairs of blocks
        - Use parallel merge via divide-and-conquer on output positions
        - Each merge level can process pairs in parallel

    Optimizations in this version:
    - OpenMP for efficient parallelism (replaces pthreads)
    - Prefetching in hot loops
    - Loop unrolling for better ILP
    - Cache-aligned allocations
    - Parallel pairwise merge using OpenMP tasks
    - Per-thread temp buffers for block sorting
*/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"

/* Minimum block size for parallel merge sort */
#define MERGE_MIN_BLOCK_SIZE 10000

/* Threshold for switching to sequential merge */
#define MERGE_SEQ_THRESHOLD 2048

/* Threshold for insertion sort in block sort */
#define MERGE_INSERTION_THRESHOLD 32

/* Maximum threads */
#define MERGE_MAX_THREADS 16

/* Use PREFETCH_DISTANCE from ctools_config.h */
#define MERGE_PREFETCH_DISTANCE PREFETCH_DISTANCE

/* Use centralized macros from ctools_config.h */
#define MERGE_RESTRICT CTOOLS_RESTRICT
#define MERGE_PREFETCH(addr) CTOOLS_PREFETCH(addr)
#define MERGE_PREFETCH_W(addr) CTOOLS_PREFETCH_W(addr)

/* ============================================================================
   Utility Functions
   ============================================================================ */

static void *merge_aligned_alloc(size_t alignment, size_t size)
{
    void *ptr = NULL;
#if defined(__APPLE__) || defined(__linux__)
    if (posix_memalign(&ptr, alignment, size) != 0) {
        return NULL;
    }
    return ptr;
#else
    return malloc(size);
#endif
}

static inline uint64_t merge_double_to_sortable(double d)
{
    uint64_t bits;
    memcpy(&bits, &d, sizeof(bits));

    if (SF_is_missing(d)) {
        return UINT64_MAX;
    }

    if (bits & ((uint64_t)1 << 63)) {
        bits = ~bits;
    } else {
        bits ^= ((uint64_t)1 << 63);
    }

    return bits;
}

/* ============================================================================
   Block Radix Sort (for Phase 1)
   ============================================================================ */

static void block_insertion_sort_numeric(size_t * MERGE_RESTRICT order,
                                         const uint64_t * MERGE_RESTRICT keys,
                                         size_t start, size_t len)
{
    for (size_t i = 1; i < len; i++) {
        size_t j = i;
        size_t tmp = order[start + j];
        uint64_t tmp_key = keys[tmp];
        while (j > 0 && keys[order[start + j - 1]] > tmp_key) {
            order[start + j] = order[start + j - 1];
            j--;
        }
        order[start + j] = tmp;
    }
}

static void block_radix_sort_numeric(size_t * MERGE_RESTRICT order,
                                     const uint64_t * MERGE_RESTRICT keys,
                                     size_t start, size_t len,
                                     size_t * MERGE_RESTRICT temp)
{
    size_t counts[256];
    size_t offsets[256];
    size_t *src = order + start;
    size_t *dst = temp;

    if (len < 2) return;

    if (len < MERGE_INSERTION_THRESHOLD) {
        block_insertion_sort_numeric(order, keys, start, len);
        return;
    }

    for (int byte_pos = 0; byte_pos < 8; byte_pos++) {
        int shift = byte_pos * 8;

        memset(counts, 0, sizeof(counts));

        /* Counting with prefetch and unrolling */
        size_t i = 0;
        for (; i + 4 <= len; i += 4) {
            if (i + MERGE_PREFETCH_DISTANCE < len) {
                MERGE_PREFETCH(&keys[src[i + MERGE_PREFETCH_DISTANCE]]);
            }
            counts[(keys[src[i]] >> shift) & 0xFF]++;
            counts[(keys[src[i+1]] >> shift) & 0xFF]++;
            counts[(keys[src[i+2]] >> shift) & 0xFF]++;
            counts[(keys[src[i+3]] >> shift) & 0xFF]++;
        }
        for (; i < len; i++) {
            counts[(keys[src[i]] >> shift) & 0xFF]++;
        }

        /* Check for uniform distribution */
        int non_empty = 0;
        for (size_t b = 0; b < 256; b++) {
            if (counts[b] > 0) {
                non_empty++;
                if (non_empty > 1) break;
            }
        }
        if (non_empty <= 1) continue;

        /* Prefix sum */
        offsets[0] = 0;
        for (size_t b = 1; b < 256; b++) {
            offsets[b] = offsets[b-1] + counts[b-1];
        }

        /* Scatter with prefetch */
        for (i = 0; i < len; i++) {
            if (i + MERGE_PREFETCH_DISTANCE < len) {
                MERGE_PREFETCH(&keys[src[i + MERGE_PREFETCH_DISTANCE]]);
            }
            uint8_t byte_val = (keys[src[i]] >> shift) & 0xFF;
            dst[offsets[byte_val]++] = src[i];
        }

        /* Swap */
        size_t *tmp = src;
        src = dst;
        dst = tmp;
    }

    /* Copy back if needed */
    if (src != order + start) {
        memcpy(order + start, src, len * sizeof(size_t));
    }
}

static void block_insertion_sort_string(size_t * MERGE_RESTRICT order,
                                        char * const * MERGE_RESTRICT strings,
                                        size_t start, size_t len)
{
    for (size_t i = 1; i < len; i++) {
        size_t j = i;
        size_t tmp = order[start + j];
        const char *tmp_str = strings[tmp];
        while (j > 0 && strcmp(strings[order[start + j - 1]], tmp_str) > 0) {
            order[start + j] = order[start + j - 1];
            j--;
        }
        order[start + j] = tmp;
    }
}

static void block_radix_sort_string(size_t * MERGE_RESTRICT order,
                                    char * const * MERGE_RESTRICT strings,
                                    size_t start, size_t len,
                                    size_t * MERGE_RESTRICT temp,
                                    size_t max_len)
{
    size_t counts[257];
    size_t offsets[257];

    if (len < 2) return;

    if (len < MERGE_INSERTION_THRESHOLD) {
        block_insertion_sort_string(order, strings, start, len);
        return;
    }

    /* MSD radix sort */
    for (size_t char_pos = 0; char_pos < max_len; char_pos++) {
        memset(counts, 0, sizeof(counts));

        for (size_t i = 0; i < len; i++) {
            const char *s = strings[order[start + i]];
            size_t slen = strlen(s);
            int c = (char_pos < slen) ? (unsigned char)s[char_pos] + 1 : 0;
            counts[c]++;
        }

        /* Check for uniform */
        int non_empty = 0;
        for (size_t b = 0; b < 257; b++) {
            if (counts[b] > 0) {
                non_empty++;
                if (non_empty > 1) break;
            }
        }
        if (non_empty <= 1) continue;

        /* Prefix sum */
        offsets[0] = 0;
        for (size_t b = 1; b < 257; b++) {
            offsets[b] = offsets[b-1] + counts[b-1];
        }

        /* Scatter */
        for (size_t i = 0; i < len; i++) {
            const char *s = strings[order[start + i]];
            size_t slen = strlen(s);
            int c = (char_pos < slen) ? (unsigned char)s[char_pos] + 1 : 0;
            temp[offsets[c]++] = order[start + i];
        }

        /* Copy back */
        memcpy(order + start, temp, len * sizeof(size_t));
    }
}

/* ============================================================================
   Binary Search for Merge
   ============================================================================ */

static inline size_t binary_search_numeric(const size_t * MERGE_RESTRICT order,
                                           const uint64_t * MERGE_RESTRICT keys,
                                           size_t start, size_t end,
                                           uint64_t target)
{
    while (start < end) {
        size_t mid = start + (end - start) / 2;
        if (keys[order[mid]] < target) {
            start = mid + 1;
        } else {
            end = mid;
        }
    }
    return start;
}

static inline size_t binary_search_string(const size_t * MERGE_RESTRICT order,
                                          char * const * MERGE_RESTRICT strings,
                                          size_t start, size_t end,
                                          const char *target)
{
    while (start < end) {
        size_t mid = start + (end - start) / 2;
        if (strcmp(strings[order[mid]], target) < 0) {
            start = mid + 1;
        } else {
            end = mid;
        }
    }
    return start;
}

/* ============================================================================
   Sequential Merge
   ============================================================================ */

static void seq_merge_numeric(const size_t * MERGE_RESTRICT order,
                              size_t * MERGE_RESTRICT output,
                              const uint64_t * MERGE_RESTRICT keys,
                              size_t start1, size_t len1,
                              size_t start2, size_t len2,
                              size_t out_start)
{
    size_t i = 0, j = 0, k = out_start;
    size_t end1 = start1 + len1;
    size_t end2 = start2 + len2;

    /* Main merge loop with prefetching */
    while (start1 + i < end1 && start2 + j < end2) {
        if (i + MERGE_PREFETCH_DISTANCE < len1) {
            MERGE_PREFETCH(&keys[order[start1 + i + MERGE_PREFETCH_DISTANCE]]);
        }
        if (j + MERGE_PREFETCH_DISTANCE < len2) {
            MERGE_PREFETCH(&keys[order[start2 + j + MERGE_PREFETCH_DISTANCE]]);
        }

        if (keys[order[start1 + i]] <= keys[order[start2 + j]]) {
            output[k++] = order[start1 + i++];
        } else {
            output[k++] = order[start2 + j++];
        }
    }

    /* Copy remaining elements */
    while (start1 + i < end1) {
        output[k++] = order[start1 + i++];
    }
    while (start2 + j < end2) {
        output[k++] = order[start2 + j++];
    }
}

static void seq_merge_string(const size_t * MERGE_RESTRICT order,
                             size_t * MERGE_RESTRICT output,
                             char * const * MERGE_RESTRICT strings,
                             size_t start1, size_t len1,
                             size_t start2, size_t len2,
                             size_t out_start)
{
    size_t i = 0, j = 0, k = out_start;
    size_t end1 = start1 + len1;
    size_t end2 = start2 + len2;

    while (start1 + i < end1 && start2 + j < end2) {
        if (strcmp(strings[order[start1 + i]], strings[order[start2 + j]]) <= 0) {
            output[k++] = order[start1 + i++];
        } else {
            output[k++] = order[start2 + j++];
        }
    }

    while (start1 + i < end1) {
        output[k++] = order[start1 + i++];
    }
    while (start2 + j < end2) {
        output[k++] = order[start2 + j++];
    }
}

/* ============================================================================
   Parallel Merge using Divide and Conquer
   ============================================================================ */

static void parallel_merge_numeric(const size_t * MERGE_RESTRICT order,
                                   size_t * MERGE_RESTRICT output,
                                   const uint64_t * MERGE_RESTRICT keys,
                                   size_t start1, size_t len1,
                                   size_t start2, size_t len2,
                                   size_t out_start, int depth)
{
    /* Base case: sequential merge */
    if (len1 + len2 < MERGE_SEQ_THRESHOLD || depth > 3) {
        seq_merge_numeric(order, output, keys, start1, len1, start2, len2, out_start);
        return;
    }

    /* Divide: split the larger array at midpoint */
    if (len1 >= len2) {
        size_t mid1 = len1 / 2;
        uint64_t pivot = keys[order[start1 + mid1]];
        size_t mid2 = binary_search_numeric(order, keys, start2, start2 + len2, pivot) - start2;

        /* Spawn parallel tasks */
        #pragma omp task if(depth < 2)
        parallel_merge_numeric(order, output, keys,
                               start1, mid1,
                               start2, mid2,
                               out_start, depth + 1);

        #pragma omp task if(depth < 2)
        parallel_merge_numeric(order, output, keys,
                               start1 + mid1, len1 - mid1,
                               start2 + mid2, len2 - mid2,
                               out_start + mid1 + mid2, depth + 1);

        #pragma omp taskwait
    } else {
        size_t mid2 = len2 / 2;
        uint64_t pivot = keys[order[start2 + mid2]];
        size_t mid1 = binary_search_numeric(order, keys, start1, start1 + len1, pivot) - start1;

        #pragma omp task if(depth < 2)
        parallel_merge_numeric(order, output, keys,
                               start1, mid1,
                               start2, mid2,
                               out_start, depth + 1);

        /* Copy pivot */
        output[out_start + mid1 + mid2] = order[start2 + mid2];

        #pragma omp task if(depth < 2)
        parallel_merge_numeric(order, output, keys,
                               start1 + mid1, len1 - mid1,
                               start2 + mid2 + 1, len2 - mid2 - 1,
                               out_start + mid1 + mid2 + 1, depth + 1);

        #pragma omp taskwait
    }
}

static void parallel_merge_string(const size_t * MERGE_RESTRICT order,
                                  size_t * MERGE_RESTRICT output,
                                  char * const * MERGE_RESTRICT strings,
                                  size_t start1, size_t len1,
                                  size_t start2, size_t len2,
                                  size_t out_start, int depth)
{
    /* Base case: sequential merge */
    if (len1 + len2 < MERGE_SEQ_THRESHOLD || depth > 3) {
        seq_merge_string(order, output, strings, start1, len1, start2, len2, out_start);
        return;
    }

    if (len1 >= len2) {
        size_t mid1 = len1 / 2;
        const char *pivot = strings[order[start1 + mid1]];
        size_t mid2 = binary_search_string(order, strings, start2, start2 + len2, pivot) - start2;

        #pragma omp task if(depth < 2)
        parallel_merge_string(order, output, strings,
                              start1, mid1,
                              start2, mid2,
                              out_start, depth + 1);

        #pragma omp task if(depth < 2)
        parallel_merge_string(order, output, strings,
                              start1 + mid1, len1 - mid1,
                              start2 + mid2, len2 - mid2,
                              out_start + mid1 + mid2, depth + 1);

        #pragma omp taskwait
    } else {
        size_t mid2 = len2 / 2;
        const char *pivot = strings[order[start2 + mid2]];
        size_t mid1 = binary_search_string(order, strings, start1, start1 + len1, pivot) - start1;

        #pragma omp task if(depth < 2)
        parallel_merge_string(order, output, strings,
                              start1, mid1,
                              start2, mid2,
                              out_start, depth + 1);

        output[out_start + mid1 + mid2] = order[start2 + mid2];

        #pragma omp task if(depth < 2)
        parallel_merge_string(order, output, strings,
                              start1 + mid1, len1 - mid1,
                              start2 + mid2 + 1, len2 - mid2 - 1,
                              out_start + mid1 + mid2 + 1, depth + 1);

        #pragma omp taskwait
    }
}

/* ============================================================================
   Main Parallel Merge Sort Implementation
   ============================================================================ */

static stata_retcode parallel_merge_sort_numeric(size_t * MERGE_RESTRICT order,
                                                 const uint64_t * MERGE_RESTRICT keys,
                                                 size_t nobs, int num_threads)
{
    stata_retcode rc = STATA_OK;

    if (num_threads > MERGE_MAX_THREADS) num_threads = MERGE_MAX_THREADS;

    size_t block_size = (nobs + num_threads - 1) / num_threads;

    /* Allocate buffers */
    size_t *temp = (size_t *)merge_aligned_alloc(64, nobs * sizeof(size_t));
    size_t *merge_temp = (size_t *)merge_aligned_alloc(64, nobs * sizeof(size_t));
    size_t *block_starts = (size_t *)malloc(num_threads * sizeof(size_t));
    size_t *block_lens = (size_t *)malloc(num_threads * sizeof(size_t));
    size_t **thread_temps = (size_t **)malloc(num_threads * sizeof(size_t *));

    if (!temp || !merge_temp || !block_starts || !block_lens || !thread_temps) {
        rc = STATA_ERR_MEMORY;
        goto cleanup;
    }

    /* Pre-allocate per-thread temp buffers */
    for (int t = 0; t < num_threads; t++) {
        thread_temps[t] = NULL;
    }

    for (int t = 0; t < num_threads; t++) {
        size_t start = (size_t)t * block_size;
        size_t end = (size_t)(t + 1) * block_size;
        if (end > nobs) end = nobs;
        if (start >= nobs) start = end = nobs;

        block_starts[t] = start;
        block_lens[t] = end - start;

        if (block_lens[t] > 0) {
            thread_temps[t] = (size_t *)malloc(block_lens[t] * sizeof(size_t));
            if (!thread_temps[t]) {
                rc = STATA_ERR_MEMORY;
                goto cleanup;
            }
        }
    }

    /* Phase 1: Parallel block sort */
    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (int t = 0; t < num_threads; t++) {
        if (block_lens[t] > 0) {
            block_radix_sort_numeric(order, keys, block_starts[t],
                                     block_lens[t], thread_temps[t]);
        }
    }

    /* Phase 2: Iterative pairwise merge */
    size_t *current = order;
    size_t *next = merge_temp;

    int num_blocks = num_threads;
    while (num_blocks > 1) {
        int new_num_blocks = 0;

        /* Process pairs in parallel */
        #pragma omp parallel num_threads(num_threads)
        {
            #pragma omp single
            {
                for (int i = 0; i < num_blocks; i += 2) {
                    if (i + 1 < num_blocks) {
                        /* Merge blocks i and i+1 */
                        size_t out_start = block_starts[i];

                        #pragma omp task
                        {
                            parallel_merge_numeric(current, next, keys,
                                                   block_starts[i], block_lens[i],
                                                   block_starts[i + 1], block_lens[i + 1],
                                                   out_start, 0);
                        }

                        /* Update block info for next level (done after taskwait) */
                    } else {
                        /* Odd block - just copy */
                        #pragma omp task
                        {
                            memcpy(next + block_starts[i], current + block_starts[i],
                                   block_lens[i] * sizeof(size_t));
                        }
                    }
                }

                #pragma omp taskwait

                /* Update block metadata after all merges complete */
                new_num_blocks = 0;
                for (int i = 0; i < num_blocks; i += 2) {
                    if (i + 1 < num_blocks) {
                        block_starts[new_num_blocks] = block_starts[i];
                        block_lens[new_num_blocks] = block_lens[i] + block_lens[i + 1];
                    } else {
                        block_starts[new_num_blocks] = block_starts[i];
                        block_lens[new_num_blocks] = block_lens[i];
                    }
                    new_num_blocks++;
                }
            }
        }

        /* Swap buffers */
        size_t *tmp = current;
        current = next;
        next = tmp;

        num_blocks = new_num_blocks;
    }

    /* Copy result to order if needed */
    if (current != order) {
        memcpy(order, current, nobs * sizeof(size_t));
    }

cleanup:
    free(temp);
    free(merge_temp);
    free(block_starts);
    free(block_lens);
    if (thread_temps) {
        for (int t = 0; t < num_threads; t++) {
            free(thread_temps[t]);
        }
        free(thread_temps);
    }

    return rc;
}

static stata_retcode parallel_merge_sort_string(size_t * MERGE_RESTRICT order,
                                                char * const * MERGE_RESTRICT strings,
                                                size_t nobs, int num_threads,
                                                size_t max_len)
{
    stata_retcode rc = STATA_OK;

    if (num_threads > MERGE_MAX_THREADS) num_threads = MERGE_MAX_THREADS;

    size_t block_size = (nobs + num_threads - 1) / num_threads;

    /* Allocate buffers */
    size_t *temp = (size_t *)merge_aligned_alloc(64, nobs * sizeof(size_t));
    size_t *merge_temp = (size_t *)merge_aligned_alloc(64, nobs * sizeof(size_t));
    size_t *block_starts = (size_t *)malloc(num_threads * sizeof(size_t));
    size_t *block_lens = (size_t *)malloc(num_threads * sizeof(size_t));
    size_t **thread_temps = (size_t **)malloc(num_threads * sizeof(size_t *));

    if (!temp || !merge_temp || !block_starts || !block_lens || !thread_temps) {
        rc = STATA_ERR_MEMORY;
        goto cleanup;
    }

    for (int t = 0; t < num_threads; t++) {
        thread_temps[t] = NULL;
    }

    for (int t = 0; t < num_threads; t++) {
        size_t start = (size_t)t * block_size;
        size_t end = (size_t)(t + 1) * block_size;
        if (end > nobs) end = nobs;
        if (start >= nobs) start = end = nobs;

        block_starts[t] = start;
        block_lens[t] = end - start;

        if (block_lens[t] > 0) {
            thread_temps[t] = (size_t *)malloc(block_lens[t] * sizeof(size_t));
            if (!thread_temps[t]) {
                rc = STATA_ERR_MEMORY;
                goto cleanup;
            }
        }
    }

    /* Phase 1: Parallel block sort */
    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (int t = 0; t < num_threads; t++) {
        if (block_lens[t] > 0) {
            block_radix_sort_string(order, strings, block_starts[t],
                                    block_lens[t], thread_temps[t], max_len);
        }
    }

    /* Phase 2: Iterative pairwise merge */
    size_t *current = order;
    size_t *next = merge_temp;

    int num_blocks = num_threads;
    while (num_blocks > 1) {
        int new_num_blocks = 0;

        #pragma omp parallel num_threads(num_threads)
        {
            #pragma omp single
            {
                for (int i = 0; i < num_blocks; i += 2) {
                    if (i + 1 < num_blocks) {
                        size_t out_start = block_starts[i];

                        #pragma omp task
                        {
                            parallel_merge_string(current, next, strings,
                                                  block_starts[i], block_lens[i],
                                                  block_starts[i + 1], block_lens[i + 1],
                                                  out_start, 0);
                        }
                    } else {
                        #pragma omp task
                        {
                            memcpy(next + block_starts[i], current + block_starts[i],
                                   block_lens[i] * sizeof(size_t));
                        }
                    }
                }

                #pragma omp taskwait

                new_num_blocks = 0;
                for (int i = 0; i < num_blocks; i += 2) {
                    if (i + 1 < num_blocks) {
                        block_starts[new_num_blocks] = block_starts[i];
                        block_lens[new_num_blocks] = block_lens[i] + block_lens[i + 1];
                    } else {
                        block_starts[new_num_blocks] = block_starts[i];
                        block_lens[new_num_blocks] = block_lens[i];
                    }
                    new_num_blocks++;
                }
            }
        }

        size_t *tmp = current;
        current = next;
        next = tmp;

        num_blocks = new_num_blocks;
    }

    if (current != order) {
        memcpy(order, current, nobs * sizeof(size_t));
    }

cleanup:
    free(temp);
    free(merge_temp);
    free(block_starts);
    free(block_lens);
    if (thread_temps) {
        for (int t = 0; t < num_threads; t++) {
            free(thread_temps[t]);
        }
        free(thread_temps);
    }

    return rc;
}

/* ============================================================================
   Variable-Level Sort Functions
   ============================================================================ */

static stata_retcode merge_sort_by_numeric_var(stata_data *data, int var_idx)
{
    uint64_t *keys;
    double *dbl_data;
    int num_threads;
    stata_retcode rc;

    keys = (uint64_t *)merge_aligned_alloc(64, data->nobs * sizeof(uint64_t));
    if (keys == NULL) {
        return STATA_ERR_MEMORY;
    }

    /* Convert doubles to sortable keys in parallel */
    dbl_data = data->vars[var_idx].data.dbl;
    #pragma omp parallel for
    for (size_t i = 0; i < data->nobs; i++) {
        keys[i] = merge_double_to_sortable(dbl_data[i]);
    }

    /* Determine thread count */
    num_threads = omp_get_max_threads();
    if (num_threads > MERGE_MAX_THREADS) num_threads = MERGE_MAX_THREADS;

    if (data->nobs < MERGE_MIN_BLOCK_SIZE * 2) {
        num_threads = 1;
    } else if (data->nobs < (size_t)MERGE_MIN_BLOCK_SIZE * (size_t)num_threads) {
        num_threads = (int)(data->nobs / MERGE_MIN_BLOCK_SIZE);
        if (num_threads < 2) num_threads = 2;
    }

    if (num_threads < 2) {
        /* Sequential radix sort */
        size_t *temp = (size_t *)malloc(data->nobs * sizeof(size_t));
        if (temp == NULL) {
            free(keys);
            return STATA_ERR_MEMORY;
        }
        block_radix_sort_numeric(data->sort_order, keys, 0, data->nobs, temp);
        free(temp);
        rc = STATA_OK;
    } else {
        rc = parallel_merge_sort_numeric(data->sort_order, keys, data->nobs, num_threads);
    }

    free(keys);
    return rc;
}

static stata_retcode merge_sort_by_string_var(stata_data *data, int var_idx)
{
    char **strings = data->vars[var_idx].data.str;
    size_t max_len = 0;
    int num_threads;
    stata_retcode rc;

    /* Find max string length */
    for (size_t i = 0; i < data->nobs; i++) {
        size_t len = strlen(strings[i]);
        if (len > max_len) max_len = len;
    }

    /* Determine thread count */
    num_threads = omp_get_max_threads();
    if (num_threads > MERGE_MAX_THREADS) num_threads = MERGE_MAX_THREADS;

    if (data->nobs < MERGE_MIN_BLOCK_SIZE * 2) {
        num_threads = 1;
    } else if (data->nobs < (size_t)MERGE_MIN_BLOCK_SIZE * (size_t)num_threads) {
        num_threads = (int)(data->nobs / MERGE_MIN_BLOCK_SIZE);
        if (num_threads < 2) num_threads = 2;
    }

    if (num_threads < 2) {
        size_t *temp = (size_t *)malloc(data->nobs * sizeof(size_t));
        if (temp == NULL) {
            return STATA_ERR_MEMORY;
        }
        block_radix_sort_string(data->sort_order, strings, 0, data->nobs, temp, max_len);
        free(temp);
        rc = STATA_OK;
    } else {
        rc = parallel_merge_sort_string(data->sort_order, strings, data->nobs,
                                        num_threads, max_len);
    }

    return rc;
}

/* ============================================================================
   Permutation Application (Parallel)
   ============================================================================ */

static stata_retcode merge_apply_permutation(stata_data *data)
{
    size_t nvars = data->nvars;
    size_t nobs = data->nobs;
    int all_success = 1;

    if (nvars == 0) {
        return STATA_OK;
    }

    /* Apply permutation to each variable in parallel */
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t j = 0; j < nvars; j++) {
        stata_variable *var = &data->vars[j];
        size_t *perm = data->sort_order;

        if (var->type == STATA_TYPE_DOUBLE) {
            double *old_data = var->data.dbl;
            double *new_data = (double *)merge_aligned_alloc(64, nobs * sizeof(double));
            if (new_data == NULL) {
                #pragma omp atomic write
                all_success = 0;
            } else {
                for (size_t i = 0; i < nobs; i++) {
                    new_data[i] = old_data[perm[i]];
                }
                free(old_data);
                var->data.dbl = new_data;
            }
        } else {
            char **old_data = var->data.str;
            char **new_data = (char **)malloc(nobs * sizeof(char *));
            if (new_data == NULL) {
                #pragma omp atomic write
                all_success = 0;
            } else {
                for (size_t i = 0; i < nobs; i++) {
                    new_data[i] = old_data[perm[i]];
                }
                free(old_data);
                var->data.str = new_data;
            }
        }
    }

    /* Reset sort order */
    for (size_t j = 0; j < nobs; j++) {
        data->sort_order[j] = j;
    }

    return all_success ? STATA_OK : STATA_ERR_MEMORY;
}

/* ============================================================================
   Public API
   ============================================================================ */

stata_retcode ctools_sort_merge(stata_data *data, int *sort_vars, size_t nsort)
{
    int k;
    int var_idx;
    stata_retcode rc;

    if (data == NULL || sort_vars == NULL || data->nobs == 0 || nsort == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Sort from last key to first for stable multi-key sort */
    for (k = (int)nsort - 1; k >= 0; k--) {
        var_idx = sort_vars[k] - 1;

        if (var_idx < 0 || var_idx >= (int)data->nvars) {
            return STATA_ERR_INVALID_INPUT;
        }

        if (data->vars[var_idx].type == STATA_TYPE_DOUBLE) {
            rc = merge_sort_by_numeric_var(data, var_idx);
        } else {
            rc = merge_sort_by_string_var(data, var_idx);
        }

        if (rc != STATA_OK) {
            return rc;
        }
    }

    rc = merge_apply_permutation(data);
    if (rc != STATA_OK) {
        return rc;
    }

    return STATA_OK;
}

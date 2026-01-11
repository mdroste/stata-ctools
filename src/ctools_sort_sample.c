/*
    ctools_sort_sample.c
    Optimized Parallel Sample Sort Module

    Sample sort is the gold standard for parallel sorting. It achieves near-linear
    speedup by ensuring each of p processors handles ~n/p elements with minimal
    synchronization.

    Algorithm Overview:
    1. Sample Phase: Each thread samples s elements from its n/p chunk
    2. Splitter Selection: Combine samples, sort, pick p-1 splitters (quantiles)
    3. Classification: Each thread partitions its chunk into p buckets using splitters
    4. Bucket Size Computation: Global reduction to compute bucket offsets
    5. Data Exchange: Each thread scatters elements to correct global positions
    6. Local Sort: Each thread sorts its assigned bucket independently

    Optimizations in this version:
    - OpenMP for efficient parallelism (replaces pthreads)
    - Per-thread temp buffers to avoid race conditions
    - Software prefetching in hot loops
    - Loop unrolling for better ILP
    - Cache-aligned allocations to avoid false sharing
    - Branchless bucket finding where possible
*/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"

/* Oversampling factor: sample this many elements per splitter */
#define SAMPLE_OVERSAMPLE 16

/* Minimum observations to use sample sort (otherwise use radix) */
#define SAMPLE_SORT_THRESHOLD 50000

/* Prefetch distance for memory access optimization */
#define SAMPLE_PREFETCH_DISTANCE 16

/* Maximum threads to use */
#define SAMPLE_MAX_THREADS 16

/* Restrict keyword for compiler optimization */
#ifdef __GNUC__
#define SAMPLE_RESTRICT __restrict__
#else
#define SAMPLE_RESTRICT
#endif

/* Prefetch macro */
#if defined(__GNUC__) || defined(__clang__)
#define SAMPLE_PREFETCH(addr) __builtin_prefetch(addr, 0, 1)
#define SAMPLE_PREFETCH_W(addr) __builtin_prefetch(addr, 1, 1)
#else
#define SAMPLE_PREFETCH(addr) ((void)0)
#define SAMPLE_PREFETCH_W(addr) ((void)0)
#endif

/* ============================================================================
   Utility Functions
   ============================================================================ */

static void *sample_aligned_alloc(size_t alignment, size_t size)
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

/*
    Convert IEEE 754 double to sortable uint64.
*/
static inline uint64_t sample_double_to_sortable(double d)
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
   Optimized Radix Sort for Samples (replaces qsort)
   ============================================================================ */

static void radix_sort_samples(uint64_t * SAMPLE_RESTRICT samples, size_t n,
                               uint64_t * SAMPLE_RESTRICT temp)
{
    size_t counts[256];
    size_t offsets[256];
    uint64_t *src = samples;
    uint64_t *dst = temp;

    if (n < 32) {
        /* Insertion sort for tiny arrays */
        for (size_t i = 1; i < n; i++) {
            uint64_t key = src[i];
            size_t j = i;
            while (j > 0 && src[j-1] > key) {
                src[j] = src[j-1];
                j--;
            }
            src[j] = key;
        }
        return;
    }

    for (int byte_pos = 0; byte_pos < 8; byte_pos++) {
        int shift = byte_pos * 8;

        memset(counts, 0, sizeof(counts));

        /* Counting pass with prefetching and unrolling */
        size_t i = 0;
        for (; i + 4 <= n; i += 4) {
            if (i + SAMPLE_PREFETCH_DISTANCE < n) {
                SAMPLE_PREFETCH(&src[i + SAMPLE_PREFETCH_DISTANCE]);
            }
            counts[(src[i] >> shift) & 0xFF]++;
            counts[(src[i+1] >> shift) & 0xFF]++;
            counts[(src[i+2] >> shift) & 0xFF]++;
            counts[(src[i+3] >> shift) & 0xFF]++;
        }
        for (; i < n; i++) {
            counts[(src[i] >> shift) & 0xFF]++;
        }

        /* Check for uniform distribution (skip this byte) */
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

        /* Scatter pass */
        for (i = 0; i < n; i++) {
            uint8_t byte_val = (src[i] >> shift) & 0xFF;
            dst[offsets[byte_val]++] = src[i];
        }

        /* Swap buffers */
        uint64_t *tmp = src;
        src = dst;
        dst = tmp;
    }

    /* Copy back if needed */
    if (src != samples) {
        memcpy(samples, src, n * sizeof(uint64_t));
    }
}

/* ============================================================================
   Bucket Finding with Binary Search (optimized)
   ============================================================================ */

static inline int find_bucket_numeric(uint64_t key,
                                      const uint64_t * SAMPLE_RESTRICT splitters,
                                      int num_splitters)
{
    int lo = 0, hi = num_splitters;

    /* Unrolled binary search for common cases */
    while (hi - lo > 4) {
        int mid = lo + ((hi - lo) >> 1);
        if (key <= splitters[mid]) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }

    /* Linear scan for remaining */
    while (lo < hi && key > splitters[lo]) {
        lo++;
    }

    return lo;
}

static inline int find_bucket_string(const char *key,
                                     char * const * SAMPLE_RESTRICT splitters,
                                     int num_splitters)
{
    int lo = 0, hi = num_splitters;

    while (lo < hi) {
        int mid = lo + ((hi - lo) >> 1);
        if (strcmp(key, splitters[mid]) <= 0) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    return lo;
}

/* ============================================================================
   Local Bucket Sorting
   ============================================================================ */

/* Insertion sort for small buckets */
static void insertion_sort_numeric(size_t * SAMPLE_RESTRICT order,
                                   const uint64_t * SAMPLE_RESTRICT keys,
                                   size_t start, size_t len)
{
    for (size_t i = 1; i < len; i++) {
        size_t temp = order[start + i];
        uint64_t temp_key = keys[temp];
        size_t j = i;

        while (j > 0 && keys[order[start + j - 1]] > temp_key) {
            order[start + j] = order[start + j - 1];
            j--;
        }
        order[start + j] = temp;
    }
}

static void insertion_sort_string(size_t * SAMPLE_RESTRICT order,
                                  char * const * SAMPLE_RESTRICT strings,
                                  size_t start, size_t len)
{
    for (size_t i = 1; i < len; i++) {
        size_t temp = order[start + i];
        const char *temp_str = strings[temp];
        size_t j = i;

        while (j > 0 && strcmp(strings[order[start + j - 1]], temp_str) > 0) {
            order[start + j] = order[start + j - 1];
            j--;
        }
        order[start + j] = temp;
    }
}

/* LSD radix sort for bucket (numeric) */
static void radix_sort_bucket_numeric(size_t * SAMPLE_RESTRICT order,
                                      size_t * SAMPLE_RESTRICT temp,
                                      const uint64_t * SAMPLE_RESTRICT keys,
                                      size_t start, size_t len)
{
    size_t counts[256];
    size_t offsets[256];
    size_t *src = order + start;
    size_t *dst = temp;

    if (len < 32) {
        insertion_sort_numeric(order, keys, start, len);
        return;
    }

    for (int byte_pos = 0; byte_pos < 8; byte_pos++) {
        int shift = byte_pos * 8;

        memset(counts, 0, sizeof(counts));

        /* Counting with prefetch and unrolling */
        size_t i = 0;
        for (; i + 4 <= len; i += 4) {
            if (i + SAMPLE_PREFETCH_DISTANCE < len) {
                SAMPLE_PREFETCH(&keys[src[i + SAMPLE_PREFETCH_DISTANCE]]);
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
            if (i + SAMPLE_PREFETCH_DISTANCE < len) {
                SAMPLE_PREFETCH(&keys[src[i + SAMPLE_PREFETCH_DISTANCE]]);
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

/* MSD radix sort for bucket (string) */
static void radix_sort_bucket_string(size_t * SAMPLE_RESTRICT order,
                                     size_t * SAMPLE_RESTRICT temp,
                                     char * const * SAMPLE_RESTRICT strings,
                                     size_t start, size_t len,
                                     size_t char_pos, size_t max_len)
{
    size_t counts[257];
    size_t offsets[257];

    if (len < 32) {
        insertion_sort_string(order, strings, start, len);
        return;
    }

    if (char_pos >= max_len) return;

    /* Count by current character */
    memset(counts, 0, sizeof(counts));
    for (size_t i = 0; i < len; i++) {
        const char *s = strings[order[start + i]];
        size_t slen = strlen(s);
        int c = (char_pos < slen) ? (unsigned char)s[char_pos] + 1 : 0;
        counts[c]++;
    }

    /* Check for uniform */
    int non_empty = 0;
    for (size_t i = 0; i < 257; i++) {
        if (counts[i] > 0) {
            non_empty++;
            if (non_empty > 1) break;
        }
    }

    if (non_empty <= 1) {
        if (counts[0] < len) {
            radix_sort_bucket_string(order, temp, strings, start, len,
                                     char_pos + 1, max_len);
        }
        return;
    }

    /* Prefix sum */
    offsets[0] = 0;
    for (size_t i = 1; i < 257; i++) {
        offsets[i] = offsets[i-1] + counts[i-1];
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

    /* Recompute offsets for recursion */
    offsets[0] = 0;
    for (size_t i = 1; i < 257; i++) {
        offsets[i] = offsets[i-1] + counts[i-1];
    }

    /* Recurse on buckets */
    for (int i = 1; i < 257; i++) {
        size_t bucket_size = counts[i];
        if (bucket_size > 1) {
            radix_sort_bucket_string(order, temp, strings, start + offsets[i],
                                     bucket_size, char_pos + 1, max_len);
        }
    }
}

/* ============================================================================
   Main Sample Sort Implementation (Numeric)
   ============================================================================ */

static stata_retcode sample_sort_numeric_impl(size_t * SAMPLE_RESTRICT order,
                                              uint64_t * SAMPLE_RESTRICT keys,
                                              size_t nobs, int num_threads)
{
    stata_retcode rc = STATA_OK;

    if (num_threads > SAMPLE_MAX_THREADS) num_threads = SAMPLE_MAX_THREADS;

    int num_splitters = num_threads - 1;
    size_t chunk_size = (nobs + num_threads - 1) / num_threads;
    size_t samples_per_thread = SAMPLE_OVERSAMPLE * num_threads;
    if (samples_per_thread > chunk_size) samples_per_thread = chunk_size;
    size_t total_samples = samples_per_thread * num_threads;

    /* Allocate all buffers upfront */
    uint64_t *all_samples = (uint64_t *)malloc(total_samples * sizeof(uint64_t));
    uint64_t *sample_temp = (uint64_t *)malloc(total_samples * sizeof(uint64_t));
    uint64_t *splitters = (uint64_t *)malloc(num_splitters * sizeof(uint64_t));
    size_t *temp_order = (size_t *)sample_aligned_alloc(64, nobs * sizeof(size_t));

    /* Per-thread bucket counts: use cache-line padding to avoid false sharing */
    size_t padded_buckets = ((num_threads + 7) / 8) * 8;  /* Round up to 8 */
    size_t *bucket_counts = (size_t *)sample_aligned_alloc(64,
                                num_threads * padded_buckets * sizeof(size_t));
    size_t *bucket_sizes = (size_t *)calloc(num_threads, sizeof(size_t));
    size_t *bucket_offsets = (size_t *)malloc(num_threads * sizeof(size_t));

    /* Per-thread temp buffers for bucket sorting (FIX: avoids race condition) */
    size_t **thread_temps = (size_t **)malloc(num_threads * sizeof(size_t *));

    if (!all_samples || !sample_temp || !splitters || !temp_order ||
        !bucket_counts || !bucket_sizes || !bucket_offsets || !thread_temps) {
        rc = STATA_ERR_MEMORY;
        goto cleanup;
    }

    /* Initialize bucket counts to zero */
    memset(bucket_counts, 0, num_threads * padded_buckets * sizeof(size_t));

    /* Pre-allocate per-thread temp buffers */
    for (int t = 0; t < num_threads; t++) {
        thread_temps[t] = NULL;
    }

    /* Phase 1: Parallel sampling */
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        size_t start = (size_t)tid * chunk_size;
        size_t end = (size_t)(tid + 1) * chunk_size;
        if (end > nobs) end = nobs;

        size_t local_chunk = end - start;
        size_t step = local_chunk / samples_per_thread;
        if (step == 0) step = 1;

        uint64_t *my_samples = all_samples + tid * samples_per_thread;

        for (size_t i = 0; i < samples_per_thread && start + i * step < end; i++) {
            size_t idx = order[start + i * step];
            my_samples[i] = keys[idx];
        }
    }

    /* Phase 2: Sort samples and select splitters (sequential - small data) */
    radix_sort_samples(all_samples, total_samples, sample_temp);

    size_t step = total_samples / num_threads;
    for (int i = 0; i < num_splitters; i++) {
        splitters[i] = all_samples[(i + 1) * step];
    }

    /* Phase 3: Parallel classification - count elements per bucket per thread */
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        size_t start = (size_t)tid * chunk_size;
        size_t end = (size_t)(tid + 1) * chunk_size;
        if (end > nobs) end = nobs;

        size_t *my_counts = bucket_counts + tid * padded_buckets;

        /* Unrolled classification loop with prefetching */
        size_t i = start;
        for (; i + 4 <= end; i += 4) {
            if (i + SAMPLE_PREFETCH_DISTANCE < end) {
                SAMPLE_PREFETCH(&keys[order[i + SAMPLE_PREFETCH_DISTANCE]]);
            }

            int b0 = find_bucket_numeric(keys[order[i]], splitters, num_splitters);
            int b1 = find_bucket_numeric(keys[order[i+1]], splitters, num_splitters);
            int b2 = find_bucket_numeric(keys[order[i+2]], splitters, num_splitters);
            int b3 = find_bucket_numeric(keys[order[i+3]], splitters, num_splitters);

            my_counts[b0]++;
            my_counts[b1]++;
            my_counts[b2]++;
            my_counts[b3]++;
        }
        for (; i < end; i++) {
            int bucket = find_bucket_numeric(keys[order[i]], splitters, num_splitters);
            my_counts[bucket]++;
        }
    }

    /* Phase 4: Compute global bucket sizes and offsets */
    for (int b = 0; b < num_threads; b++) {
        bucket_sizes[b] = 0;
        for (int t = 0; t < num_threads; t++) {
            bucket_sizes[b] += bucket_counts[t * padded_buckets + b];
        }
    }

    bucket_offsets[0] = 0;
    for (int b = 1; b < num_threads; b++) {
        bucket_offsets[b] = bucket_offsets[b-1] + bucket_sizes[b-1];
    }

    /* Compute per-thread write offsets within each bucket */
    size_t *thread_bucket_offsets = (size_t *)malloc(num_threads * num_threads * sizeof(size_t));
    if (!thread_bucket_offsets) {
        rc = STATA_ERR_MEMORY;
        goto cleanup;
    }

    for (int b = 0; b < num_threads; b++) {
        size_t offset = bucket_offsets[b];
        for (int t = 0; t < num_threads; t++) {
            thread_bucket_offsets[t * num_threads + b] = offset;
            offset += bucket_counts[t * padded_buckets + b];
        }
    }

    /* Phase 5: Parallel scatter */
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        size_t start = (size_t)tid * chunk_size;
        size_t end = (size_t)(tid + 1) * chunk_size;
        if (end > nobs) end = nobs;

        /* Local copy of offsets to avoid race conditions */
        size_t local_offsets[SAMPLE_MAX_THREADS];
        for (int b = 0; b < num_threads; b++) {
            local_offsets[b] = thread_bucket_offsets[tid * num_threads + b];
        }

        /* Scatter with prefetching */
        for (size_t i = start; i < end; i++) {
            if (i + SAMPLE_PREFETCH_DISTANCE < end) {
                SAMPLE_PREFETCH(&keys[order[i + SAMPLE_PREFETCH_DISTANCE]]);
            }

            size_t idx = order[i];
            int bucket = find_bucket_numeric(keys[idx], splitters, num_splitters);
            temp_order[local_offsets[bucket]++] = idx;
        }
    }

    /* Copy scattered order back */
    memcpy(order, temp_order, nobs * sizeof(size_t));

    /* Phase 6: Allocate per-thread temp buffers for bucket sort */
    size_t max_bucket_size = 0;
    for (int b = 0; b < num_threads; b++) {
        if (bucket_sizes[b] > max_bucket_size) {
            max_bucket_size = bucket_sizes[b];
        }
    }

    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        thread_temps[tid] = (size_t *)malloc(max_bucket_size * sizeof(size_t));
    }

    /* Check allocation success */
    for (int t = 0; t < num_threads; t++) {
        if (thread_temps[t] == NULL && bucket_sizes[t] > 0) {
            rc = STATA_ERR_MEMORY;
            goto cleanup;
        }
    }

    /* Phase 6: Parallel local bucket sort */
    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (int b = 0; b < num_threads; b++) {
        if (bucket_sizes[b] > 0) {
            int tid = omp_get_thread_num();
            radix_sort_bucket_numeric(order, thread_temps[tid], keys,
                                      bucket_offsets[b], bucket_sizes[b]);
        }
    }

    free(thread_bucket_offsets);

cleanup:
    free(all_samples);
    free(sample_temp);
    free(splitters);
    free(temp_order);
    free(bucket_counts);
    free(bucket_sizes);
    free(bucket_offsets);
    if (thread_temps) {
        for (int t = 0; t < num_threads; t++) {
            free(thread_temps[t]);
        }
        free(thread_temps);
    }

    return rc;
}

/* ============================================================================
   Main Sample Sort Implementation (String)
   ============================================================================ */

static stata_retcode sample_sort_string_impl(size_t * SAMPLE_RESTRICT order,
                                             char * const * SAMPLE_RESTRICT strings,
                                             size_t nobs, int num_threads,
                                             size_t max_len)
{
    stata_retcode rc = STATA_OK;

    if (num_threads > SAMPLE_MAX_THREADS) num_threads = SAMPLE_MAX_THREADS;

    int num_splitters = num_threads - 1;
    size_t chunk_size = (nobs + num_threads - 1) / num_threads;
    size_t samples_per_thread = SAMPLE_OVERSAMPLE * num_threads;
    if (samples_per_thread > chunk_size) samples_per_thread = chunk_size;
    size_t total_samples = samples_per_thread * num_threads;

    /* Allocate buffers */
    char **all_samples = (char **)malloc(total_samples * sizeof(char *));
    char **splitters = (char **)malloc(num_splitters * sizeof(char *));
    size_t *temp_order = (size_t *)sample_aligned_alloc(64, nobs * sizeof(size_t));

    size_t padded_buckets = ((num_threads + 7) / 8) * 8;
    size_t *bucket_counts = (size_t *)sample_aligned_alloc(64,
                                num_threads * padded_buckets * sizeof(size_t));
    size_t *bucket_sizes = (size_t *)calloc(num_threads, sizeof(size_t));
    size_t *bucket_offsets = (size_t *)malloc(num_threads * sizeof(size_t));
    size_t **thread_temps = (size_t **)malloc(num_threads * sizeof(size_t *));

    if (!all_samples || !splitters || !temp_order ||
        !bucket_counts || !bucket_sizes || !bucket_offsets || !thread_temps) {
        rc = STATA_ERR_MEMORY;
        goto cleanup;
    }

    memset(bucket_counts, 0, num_threads * padded_buckets * sizeof(size_t));

    for (int t = 0; t < num_threads; t++) {
        thread_temps[t] = NULL;
    }

    /* Phase 1: Parallel sampling */
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        size_t start = (size_t)tid * chunk_size;
        size_t end = (size_t)(tid + 1) * chunk_size;
        if (end > nobs) end = nobs;

        size_t local_chunk = end - start;
        size_t step = local_chunk / samples_per_thread;
        if (step == 0) step = 1;

        char **my_samples = all_samples + tid * samples_per_thread;

        for (size_t i = 0; i < samples_per_thread && start + i * step < end; i++) {
            size_t idx = order[start + i * step];
            my_samples[i] = (char *)strings[idx];
        }
    }

    /* Phase 2: Sort samples (qsort for strings - harder to radix sort) */
    /* Simple insertion sort for samples */
    for (size_t i = 1; i < total_samples; i++) {
        char *key = all_samples[i];
        size_t j = i;
        while (j > 0 && strcmp(all_samples[j-1], key) > 0) {
            all_samples[j] = all_samples[j-1];
            j--;
        }
        all_samples[j] = key;
    }

    size_t step = total_samples / num_threads;
    for (int i = 0; i < num_splitters; i++) {
        splitters[i] = all_samples[(i + 1) * step];
    }

    /* Phase 3: Parallel classification */
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        size_t start = (size_t)tid * chunk_size;
        size_t end = (size_t)(tid + 1) * chunk_size;
        if (end > nobs) end = nobs;

        size_t *my_counts = bucket_counts + tid * padded_buckets;

        for (size_t i = start; i < end; i++) {
            int bucket = find_bucket_string(strings[order[i]], splitters, num_splitters);
            my_counts[bucket]++;
        }
    }

    /* Phase 4: Compute global offsets */
    for (int b = 0; b < num_threads; b++) {
        bucket_sizes[b] = 0;
        for (int t = 0; t < num_threads; t++) {
            bucket_sizes[b] += bucket_counts[t * padded_buckets + b];
        }
    }

    bucket_offsets[0] = 0;
    for (int b = 1; b < num_threads; b++) {
        bucket_offsets[b] = bucket_offsets[b-1] + bucket_sizes[b-1];
    }

    size_t *thread_bucket_offsets = (size_t *)malloc(num_threads * num_threads * sizeof(size_t));
    if (!thread_bucket_offsets) {
        rc = STATA_ERR_MEMORY;
        goto cleanup;
    }

    for (int b = 0; b < num_threads; b++) {
        size_t offset = bucket_offsets[b];
        for (int t = 0; t < num_threads; t++) {
            thread_bucket_offsets[t * num_threads + b] = offset;
            offset += bucket_counts[t * padded_buckets + b];
        }
    }

    /* Phase 5: Parallel scatter */
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        size_t start = (size_t)tid * chunk_size;
        size_t end = (size_t)(tid + 1) * chunk_size;
        if (end > nobs) end = nobs;

        size_t local_offsets[SAMPLE_MAX_THREADS];
        for (int b = 0; b < num_threads; b++) {
            local_offsets[b] = thread_bucket_offsets[tid * num_threads + b];
        }

        for (size_t i = start; i < end; i++) {
            size_t idx = order[i];
            int bucket = find_bucket_string(strings[idx], splitters, num_splitters);
            temp_order[local_offsets[bucket]++] = idx;
        }
    }

    memcpy(order, temp_order, nobs * sizeof(size_t));

    /* Phase 6: Allocate per-thread temps */
    size_t max_bucket_size = 0;
    for (int b = 0; b < num_threads; b++) {
        if (bucket_sizes[b] > max_bucket_size) {
            max_bucket_size = bucket_sizes[b];
        }
    }

    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        thread_temps[tid] = (size_t *)malloc(max_bucket_size * sizeof(size_t));
    }

    for (int t = 0; t < num_threads; t++) {
        if (thread_temps[t] == NULL && bucket_sizes[t] > 0) {
            rc = STATA_ERR_MEMORY;
            goto cleanup;
        }
    }

    /* Phase 6: Parallel local bucket sort */
    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (int b = 0; b < num_threads; b++) {
        if (bucket_sizes[b] > 0) {
            int tid = omp_get_thread_num();
            radix_sort_bucket_string(order, thread_temps[tid], strings,
                                     bucket_offsets[b], bucket_sizes[b], 0, max_len);
        }
    }

    free(thread_bucket_offsets);

cleanup:
    free(all_samples);
    free(splitters);
    free(temp_order);
    free(bucket_counts);
    free(bucket_sizes);
    free(bucket_offsets);
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

static stata_retcode sample_sort_by_numeric_var(stata_data *data, int var_idx)
{
    uint64_t *keys;
    double *dbl_data;
    int num_threads;
    stata_retcode rc;

    keys = (uint64_t *)sample_aligned_alloc(64, data->nobs * sizeof(uint64_t));
    if (keys == NULL) {
        return STATA_ERR_MEMORY;
    }

    /* Convert doubles to sortable keys in parallel */
    dbl_data = data->vars[var_idx].data.dbl;
    #pragma omp parallel for
    for (size_t i = 0; i < data->nobs; i++) {
        keys[i] = sample_double_to_sortable(dbl_data[i]);
    }

    /* Determine thread count */
    num_threads = omp_get_max_threads();
    if (num_threads > SAMPLE_MAX_THREADS) num_threads = SAMPLE_MAX_THREADS;

    if (data->nobs < (size_t)SAMPLE_SORT_THRESHOLD) {
        num_threads = 1;
    } else if (data->nobs < (size_t)MIN_OBS_PER_THREAD * (size_t)num_threads) {
        num_threads = (int)(data->nobs / MIN_OBS_PER_THREAD);
        if (num_threads < 2) num_threads = 2;
    }

    if (num_threads < 2) {
        /* Sequential radix sort for small data */
        size_t *temp = (size_t *)malloc(data->nobs * sizeof(size_t));
        if (temp == NULL) {
            free(keys);
            return STATA_ERR_MEMORY;
        }
        radix_sort_bucket_numeric(data->sort_order, temp, keys, 0, data->nobs);
        free(temp);
        rc = STATA_OK;
    } else {
        rc = sample_sort_numeric_impl(data->sort_order, keys, data->nobs, num_threads);
    }

    free(keys);
    return rc;
}

static stata_retcode sample_sort_by_string_var(stata_data *data, int var_idx)
{
    int num_threads;
    stata_retcode rc;
    size_t max_len = 0;

    /* Find max string length */
    char **strings = data->vars[var_idx].data.str;
    for (size_t i = 0; i < data->nobs; i++) {
        size_t len = strlen(strings[i]);
        if (len > max_len) max_len = len;
    }

    /* Determine thread count */
    num_threads = omp_get_max_threads();
    if (num_threads > SAMPLE_MAX_THREADS) num_threads = SAMPLE_MAX_THREADS;

    if (data->nobs < (size_t)SAMPLE_SORT_THRESHOLD) {
        num_threads = 1;
    } else if (data->nobs < (size_t)MIN_OBS_PER_THREAD * (size_t)num_threads) {
        num_threads = (int)(data->nobs / MIN_OBS_PER_THREAD);
        if (num_threads < 2) num_threads = 2;
    }

    if (num_threads < 2) {
        size_t *temp = (size_t *)malloc(data->nobs * sizeof(size_t));
        if (temp == NULL) {
            return STATA_ERR_MEMORY;
        }
        radix_sort_bucket_string(data->sort_order, temp, strings, 0, data->nobs, 0, max_len);
        free(temp);
        rc = STATA_OK;
    } else {
        rc = sample_sort_string_impl(data->sort_order, strings, data->nobs,
                                     num_threads, max_len);
    }

    return rc;
}

/* ============================================================================
   Permutation Application (Parallel)
   ============================================================================ */

static stata_retcode sample_apply_permutation(stata_data *data)
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
            double *new_data = (double *)sample_aligned_alloc(64, nobs * sizeof(double));
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

stata_retcode ctools_sort_sample(stata_data *data, int *sort_vars, size_t nsort)
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
            rc = sample_sort_by_numeric_var(data, var_idx);
        } else {
            rc = sample_sort_by_string_var(data, var_idx);
        }

        if (rc != STATA_OK) {
            return rc;
        }
    }

    rc = sample_apply_permutation(data);
    if (rc != STATA_OK) {
        return rc;
    }

    return STATA_OK;
}

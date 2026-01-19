/*
    ctools_sort_ips4o.c
    In-place Parallel Super Scalar Samplesort (IPS⁴o) Module - Optimized v2

    Parallelization strategy:
    - Parallel key conversion with SIMD-friendly loops
    - Parallel partition (counting + scatter)
    - Parallel bucket sorting using OpenMP parallel for
    - No nested parallelism, no recursive task spawning
    - Pre-allocated buffers to avoid malloc in parallel regions

    Reference:
    Axtmann, Witt, Ferizovic, Sanders. "In-place Parallel Super Scalar
    Samplesort (IPS⁴o)", ESA 2017.
*/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* ============================================================================
   Configuration Constants
   ============================================================================ */

/* Tuned constants for Apple M-series and modern x86 CPUs */
#define IPS4O_BLOCK_SIZE 256
#define IPS4O_LOG_BUCKETS 8
#define IPS4O_MAX_BUCKETS (1 << IPS4O_LOG_BUCKETS)
#define IPS4O_OVERSAMPLE_FACTOR 16
#define IPS4O_BASE_CASE_SIZE 16          /* Smaller base case for faster recursion */
#define IPS4O_RADIX_THRESHOLD 1024       /* Lower threshold to use radix sort earlier */
#define IPS4O_PARALLEL_THRESHOLD 30000   /* Lower threshold for earlier parallelism */
#define IPS4O_PREFETCH_DISTANCE 16       /* Larger prefetch distance */

/* ============================================================================
   Compiler Hints - Use centralized macros from ctools_config.h
   ============================================================================ */

#define IPS4O_LIKELY(x)   CTOOLS_LIKELY(x)
#define IPS4O_UNLIKELY(x) CTOOLS_UNLIKELY(x)
#define IPS4O_PREFETCH(addr) CTOOLS_PREFETCH(addr)
#define IPS4O_RESTRICT CTOOLS_RESTRICT
#define IPS4O_INLINE CTOOLS_INLINE

/* ============================================================================
   Utility Functions
   ============================================================================ */

/*
    NOTE: Aligned memory allocation is now provided by ctools_config.h
    Use ctools_aligned_alloc() and ctools_aligned_free() for cross-platform support.

    NOTE: Double-to-sortable conversion is now provided by ctools_types.h
    Use ctools_double_to_sortable(d, SF_is_missing(d)) for sorting doubles.
*/

/* ============================================================================
   Branchless Classifier
   ============================================================================ */

/* Unrolled branchless tree classification */
static IPS4O_INLINE int ips4o_classify(uint64_t key, const uint64_t * IPS4O_RESTRICT tree,
                                       int log_buckets, int num_buckets)
{
    int pos = 1;

    /* Fully unrolled for common cases */
    switch (log_buckets) {
        case 8: pos = 2 * pos + (key > tree[pos]); /* fallthrough */
        case 7: pos = 2 * pos + (key > tree[pos]); /* fallthrough */
        case 6: pos = 2 * pos + (key > tree[pos]); /* fallthrough */
        case 5: pos = 2 * pos + (key > tree[pos]); /* fallthrough */
        case 4: pos = 2 * pos + (key > tree[pos]); /* fallthrough */
        case 3: pos = 2 * pos + (key > tree[pos]); /* fallthrough */
        case 2: pos = 2 * pos + (key > tree[pos]); /* fallthrough */
        case 1: pos = 2 * pos + (key > tree[pos]); break;
        default:
            for (int level = 0; level < log_buckets; level++) {
                pos = 2 * pos + (key > tree[pos]);
            }
    }

    return pos - num_buckets;
}

static IPS4O_INLINE int ips4o_classify_string(const char *key, char ** IPS4O_RESTRICT tree,
                                              int log_buckets, int num_buckets)
{
    int pos = 1;
    for (int level = 0; level < log_buckets; level++) {
        pos = 2 * pos + (strcmp(key, tree[pos]) > 0);
    }
    return pos - num_buckets;
}

static int compare_uint64(const void *a, const void *b)
{
    uint64_t va = *(const uint64_t *)a;
    uint64_t vb = *(const uint64_t *)b;
    return (va > vb) - (va < vb);
}

static int compare_string_ptr(const void *a, const void *b)
{
    return strcmp(*(const char **)a, *(const char **)b);
}

/* Build binary search tree from sorted splitters */
static void build_tree(uint64_t *tree, const uint64_t *splitters, int node, int left, int right)
{
    if (left > right) return;
    int mid = (left + right) / 2;
    tree[node] = splitters[mid];
    build_tree(tree, splitters, 2 * node, left, mid - 1);
    build_tree(tree, splitters, 2 * node + 1, mid + 1, right);
}

static void build_tree_string(char **tree, char **splitters, int node, int left, int right)
{
    if (left > right) return;
    int mid = (left + right) / 2;
    tree[node] = splitters[mid];
    build_tree_string(tree, splitters, 2 * node, left, mid - 1);
    build_tree_string(tree, splitters, 2 * node + 1, mid + 1, right);
}

/* ============================================================================
   Base Case Sorts
   ============================================================================ */

static void insertion_sort_numeric(perm_idx_t * IPS4O_RESTRICT order,
                                   const uint64_t * IPS4O_RESTRICT keys,
                                   size_t start, size_t end)
{
    for (size_t i = start + 1; i < end; i++) {
        perm_idx_t temp = order[i];
        uint64_t temp_key = keys[temp];
        size_t j = i;
        while (j > start && keys[order[j - 1]] > temp_key) {
            order[j] = order[j - 1];
            j--;
        }
        order[j] = temp;
    }
}

static void insertion_sort_string(perm_idx_t * IPS4O_RESTRICT order,
                                  char ** IPS4O_RESTRICT strings,
                                  size_t start, size_t end)
{
    for (size_t i = start + 1; i < end; i++) {
        perm_idx_t temp = order[i];
        const char *temp_str = strings[temp];
        size_t j = i;
        while (j > start && strcmp(strings[order[j - 1]], temp_str) > 0) {
            order[j] = order[j - 1];
            j--;
        }
        order[j] = temp;
    }
}

/* LSD Radix sort for numeric data - optimized with prefetching and unrolling */
static void radix_sort_numeric(perm_idx_t * IPS4O_RESTRICT order,
                               perm_idx_t * IPS4O_RESTRICT temp,
                               const uint64_t * IPS4O_RESTRICT keys,
                               size_t start, size_t len)
{
    if (len < IPS4O_BASE_CASE_SIZE) {
        insertion_sort_numeric(order, keys, start, start + len);
        return;
    }

    size_t counts[256], offsets[256];
    perm_idx_t *src = order + start;
    perm_idx_t *dst = temp;

    for (int byte = 0; byte < 8; byte++) {
        int shift = byte * 8;

        /* Counting pass with prefetching */
        memset(counts, 0, sizeof(counts));

        size_t i = 0;
        /* Unrolled counting loop */
        for (; i + 4 <= len; i += 4) {
            IPS4O_PREFETCH(&keys[src[i + IPS4O_PREFETCH_DISTANCE]]);
            counts[(keys[src[i]] >> shift) & 0xFF]++;
            counts[(keys[src[i + 1]] >> shift) & 0xFF]++;
            counts[(keys[src[i + 2]] >> shift) & 0xFF]++;
            counts[(keys[src[i + 3]] >> shift) & 0xFF]++;
        }
        for (; i < len; i++) {
            counts[(keys[src[i]] >> shift) & 0xFF]++;
        }

        /* Skip if all in one bucket */
        size_t max_count = 0;
        for (int j = 0; j < 256; j++) {
            if (counts[j] > max_count) max_count = counts[j];
        }
        if (max_count == len) continue;

        /* Prefix sum */
        offsets[0] = 0;
        for (int j = 1; j < 256; j++) {
            offsets[j] = offsets[j - 1] + counts[j - 1];
        }

        /* Scatter pass with prefetching */
        for (i = 0; i < len; i++) {
            if (IPS4O_LIKELY(i + IPS4O_PREFETCH_DISTANCE < len)) {
                IPS4O_PREFETCH(&keys[src[i + IPS4O_PREFETCH_DISTANCE]]);
            }
            uint8_t b = (keys[src[i]] >> shift) & 0xFF;
            dst[offsets[b]++] = src[i];
        }

        perm_idx_t *t = src; src = dst; dst = t;
    }

    if (src != order + start) {
        memcpy(order + start, src, len * sizeof(perm_idx_t));
    }
}

/* MSD Radix sort for strings */
static void radix_sort_string(perm_idx_t * IPS4O_RESTRICT order,
                              perm_idx_t * IPS4O_RESTRICT temp,
                              char ** IPS4O_RESTRICT strings,
                              size_t start, size_t len,
                              size_t char_pos, size_t max_len)
{
    if (len < IPS4O_BASE_CASE_SIZE) {
        insertion_sort_string(order, strings, start, start + len);
        return;
    }
    if (char_pos >= max_len) return;

    size_t counts[257], offsets[257];
    memset(counts, 0, sizeof(counts));

    for (size_t i = 0; i < len; i++) {
        const char *s = strings[order[start + i]];
        size_t slen = strlen(s);
        int c = (char_pos < slen) ? (unsigned char)s[char_pos] + 1 : 0;
        counts[c]++;
    }

    size_t max_count = 0;
    for (int j = 0; j < 257; j++) {
        if (counts[j] > max_count) max_count = counts[j];
    }
    if (max_count == len) {
        if (counts[0] < len) {
            radix_sort_string(order, temp, strings, start, len, char_pos + 1, max_len);
        }
        return;
    }

    offsets[0] = 0;
    for (int j = 1; j < 257; j++) {
        offsets[j] = offsets[j - 1] + counts[j - 1];
    }

    for (size_t i = 0; i < len; i++) {
        const char *s = strings[order[start + i]];
        size_t slen = strlen(s);
        int c = (char_pos < slen) ? (unsigned char)s[char_pos] + 1 : 0;
        temp[offsets[c]++] = order[start + i];
    }

    memcpy(order + start, temp, len * sizeof(perm_idx_t));

    /* Recompute offsets for recursion */
    offsets[0] = 0;
    for (int j = 1; j < 257; j++) {
        offsets[j] = offsets[j - 1] + counts[j - 1];
    }

    for (int j = 1; j < 257; j++) {
        if (counts[j] > 1) {
            radix_sort_string(order, temp, strings, start + offsets[j],
                              counts[j], char_pos + 1, max_len);
        }
    }
}

/* ============================================================================
   Sequential IPS4o (for small data or within parallel buckets)
   ============================================================================ */

static void ips4o_sequential_numeric(perm_idx_t * IPS4O_RESTRICT order,
                                     perm_idx_t * IPS4O_RESTRICT temp,
                                     uint64_t * IPS4O_RESTRICT keys,
                                     size_t start, size_t len, int depth)
{
    if (len <= IPS4O_RADIX_THRESHOLD || depth > 20) {
        radix_sort_numeric(order, temp, keys, start, len);
        return;
    }

    /* Determine bucket count */
    int log_buckets = IPS4O_LOG_BUCKETS;
    while (log_buckets > 2 && (len / (1 << log_buckets)) < IPS4O_BLOCK_SIZE) {
        log_buckets--;
    }
    int num_buckets = 1 << log_buckets;

    /* Sample and select splitters */
    size_t num_samples = IPS4O_OVERSAMPLE_FACTOR * num_buckets;
    if (num_samples > len) num_samples = len;

    uint64_t *samples = (uint64_t *)malloc(num_samples * sizeof(uint64_t));
    if (!samples) { radix_sort_numeric(order, temp, keys, start, len); return; }

    size_t step = len / num_samples;
    if (step == 0) step = 1;
    for (size_t i = 0; i < num_samples; i++) {
        samples[i] = keys[order[start + i * step]];
    }
    qsort(samples, num_samples, sizeof(uint64_t), compare_uint64);

    uint64_t *splitters = (uint64_t *)malloc((num_buckets - 1) * sizeof(uint64_t));
    if (!splitters) { free(samples); radix_sort_numeric(order, temp, keys, start, len); return; }

    size_t splitter_step = num_samples / num_buckets;
    for (int i = 0; i < num_buckets - 1; i++) {
        splitters[i] = samples[(i + 1) * splitter_step];
    }
    free(samples);

    /* Build classifier tree */
    uint64_t *tree = (uint64_t *)ctools_aligned_alloc(64, 2 * num_buckets * sizeof(uint64_t));
    if (!tree) { free(splitters); radix_sort_numeric(order, temp, keys, start, len); return; }
    memset(tree, 0, 2 * num_buckets * sizeof(uint64_t));
    build_tree(tree, splitters, 1, 0, num_buckets - 2);
    free(splitters);

    /* Count and partition */
    size_t *bucket_sizes = (size_t *)calloc(num_buckets, sizeof(size_t));
    size_t *bucket_offsets = (size_t *)malloc((num_buckets + 1) * sizeof(size_t));
    if (!bucket_sizes || !bucket_offsets) {
        ctools_aligned_free(tree); free(bucket_sizes); free(bucket_offsets);
        radix_sort_numeric(order, temp, keys, start, len);
        return;
    }

    /* Count with prefetching */
    for (size_t i = 0; i < len; i++) {
        if (i + IPS4O_PREFETCH_DISTANCE < len) {
            IPS4O_PREFETCH(&keys[order[start + i + IPS4O_PREFETCH_DISTANCE]]);
        }
        bucket_sizes[ips4o_classify(keys[order[start + i]], tree, log_buckets, num_buckets)]++;
    }

    bucket_offsets[0] = 0;
    for (int b = 1; b <= num_buckets; b++) {
        bucket_offsets[b] = bucket_offsets[b - 1] + bucket_sizes[b - 1];
    }

    /* Scatter */
    size_t *write_pos = (size_t *)malloc(num_buckets * sizeof(size_t));
    if (!write_pos) {
        ctools_aligned_free(tree); free(bucket_sizes); free(bucket_offsets);
        radix_sort_numeric(order, temp, keys, start, len);
        return;
    }
    memcpy(write_pos, bucket_offsets, num_buckets * sizeof(size_t));

    for (size_t i = 0; i < len; i++) {
        int b = ips4o_classify(keys[order[start + i]], tree, log_buckets, num_buckets);
        temp[write_pos[b]++] = order[start + i];
    }
    memcpy(order + start, temp, len * sizeof(perm_idx_t));

    ctools_aligned_free(tree);
    free(write_pos);

    /* Recurse on buckets */
    size_t offset = start;
    for (int b = 0; b < num_buckets; b++) {
        if (bucket_sizes[b] > 1) {
            ips4o_sequential_numeric(order, temp, keys, offset, bucket_sizes[b], depth + 1);
        }
        offset += bucket_sizes[b];
    }

    free(bucket_sizes);
    free(bucket_offsets);
}

static void ips4o_sequential_string(perm_idx_t * IPS4O_RESTRICT order,
                                    perm_idx_t * IPS4O_RESTRICT temp,
                                    char ** IPS4O_RESTRICT strings,
                                    size_t start, size_t len,
                                    size_t max_len, int depth)
{
    if (len <= IPS4O_RADIX_THRESHOLD || depth > 20) {
        radix_sort_string(order, temp, strings, start, len, 0, max_len);
        return;
    }

    int log_buckets = IPS4O_LOG_BUCKETS;
    while (log_buckets > 2 && (len / (1 << log_buckets)) < IPS4O_BLOCK_SIZE) {
        log_buckets--;
    }
    int num_buckets = 1 << log_buckets;

    size_t num_samples = IPS4O_OVERSAMPLE_FACTOR * num_buckets;
    if (num_samples > len) num_samples = len;

    char **samples = (char **)malloc(num_samples * sizeof(char *));
    if (!samples) { radix_sort_string(order, temp, strings, start, len, 0, max_len); return; }

    size_t step = len / num_samples;
    if (step == 0) step = 1;
    for (size_t i = 0; i < num_samples; i++) {
        samples[i] = strings[order[start + i * step]];
    }
    qsort(samples, num_samples, sizeof(char *), compare_string_ptr);

    char **splitters = (char **)malloc((num_buckets - 1) * sizeof(char *));
    if (!splitters) { free(samples); radix_sort_string(order, temp, strings, start, len, 0, max_len); return; }

    size_t splitter_step = num_samples / num_buckets;
    for (int i = 0; i < num_buckets - 1; i++) {
        splitters[i] = samples[(i + 1) * splitter_step];
    }
    free(samples);

    char **tree = (char **)calloc(2 * num_buckets, sizeof(char *));
    if (!tree) { free(splitters); radix_sort_string(order, temp, strings, start, len, 0, max_len); return; }
    build_tree_string(tree, splitters, 1, 0, num_buckets - 2);
    free(splitters);

    size_t *bucket_sizes = (size_t *)calloc(num_buckets, sizeof(size_t));
    size_t *bucket_offsets = (size_t *)malloc((num_buckets + 1) * sizeof(size_t));
    if (!bucket_sizes || !bucket_offsets) {
        free(tree); free(bucket_sizes); free(bucket_offsets);
        radix_sort_string(order, temp, strings, start, len, 0, max_len);
        return;
    }

    for (size_t i = 0; i < len; i++) {
        bucket_sizes[ips4o_classify_string(strings[order[start + i]], tree, log_buckets, num_buckets)]++;
    }

    bucket_offsets[0] = 0;
    for (int b = 1; b <= num_buckets; b++) {
        bucket_offsets[b] = bucket_offsets[b - 1] + bucket_sizes[b - 1];
    }

    size_t *write_pos = (size_t *)malloc(num_buckets * sizeof(size_t));
    if (!write_pos) {
        free(tree); free(bucket_sizes); free(bucket_offsets);
        radix_sort_string(order, temp, strings, start, len, 0, max_len);
        return;
    }
    memcpy(write_pos, bucket_offsets, num_buckets * sizeof(size_t));

    for (size_t i = 0; i < len; i++) {
        int b = ips4o_classify_string(strings[order[start + i]], tree, log_buckets, num_buckets);
        temp[write_pos[b]++] = order[start + i];
    }
    memcpy(order + start, temp, len * sizeof(perm_idx_t));

    free(tree);
    free(write_pos);

    size_t offset = start;
    for (int b = 0; b < num_buckets; b++) {
        if (bucket_sizes[b] > 1) {
            ips4o_sequential_string(order, temp, strings, offset, bucket_sizes[b], max_len, depth + 1);
        }
        offset += bucket_sizes[b];
    }

    free(bucket_sizes);
    free(bucket_offsets);
}

/* ============================================================================
   Parallel IPS4o - Single-level parallelism, no nesting
   ============================================================================ */

#ifdef _OPENMP

/* Structure to hold per-bucket work for parallel sorting */
typedef struct {
    size_t start;
    size_t len;
} bucket_work_t;

static stata_retcode ips4o_parallel_numeric(perm_idx_t * IPS4O_RESTRICT order,
                                            uint64_t * IPS4O_RESTRICT keys,
                                            size_t nobs)
{
    int num_threads = omp_get_max_threads();
    if (num_threads > 16) num_threads = 16;  /* Cap threads */

    /* Allocate global temp buffer */
    perm_idx_t *temp = (perm_idx_t *)ctools_aligned_alloc(64, nobs * sizeof(perm_idx_t));
    if (!temp) return STATA_ERR_MEMORY;

    /* Use thread count as bucket count for top-level partition */
    int log_buckets = 0;
    while ((1 << log_buckets) < num_threads) log_buckets++;
    if (log_buckets > IPS4O_LOG_BUCKETS) log_buckets = IPS4O_LOG_BUCKETS;
    int num_buckets = 1 << log_buckets;

    /* Sample */
    size_t num_samples = IPS4O_OVERSAMPLE_FACTOR * num_buckets;
    if (num_samples > nobs) num_samples = nobs;

    uint64_t *samples = (uint64_t *)malloc(num_samples * sizeof(uint64_t));
    if (!samples) { free(temp); return STATA_ERR_MEMORY; }

    size_t step = nobs / num_samples;
    if (step == 0) step = 1;

    /* Parallel sampling */
    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < num_samples; i++) {
        samples[i] = keys[order[i * step]];
    }

    qsort(samples, num_samples, sizeof(uint64_t), compare_uint64);

    /* Select splitters */
    uint64_t *splitters = (uint64_t *)malloc((num_buckets - 1) * sizeof(uint64_t));
    if (!splitters) { free(samples); free(temp); return STATA_ERR_MEMORY; }

    size_t splitter_step = num_samples / num_buckets;
    for (int i = 0; i < num_buckets - 1; i++) {
        splitters[i] = samples[(i + 1) * splitter_step];
    }
    free(samples);

    /* Build tree */
    uint64_t *tree = (uint64_t *)ctools_aligned_alloc(64, 2 * num_buckets * sizeof(uint64_t));
    if (!tree) { free(splitters); free(temp); return STATA_ERR_MEMORY; }
    memset(tree, 0, 2 * num_buckets * sizeof(uint64_t));
    build_tree(tree, splitters, 1, 0, num_buckets - 2);
    free(splitters);

    /* Allocate per-thread count arrays */
    size_t **thread_counts = (size_t **)malloc(num_threads * sizeof(size_t *));
    if (!thread_counts) { free(tree); free(temp); return STATA_ERR_MEMORY; }

    for (int t = 0; t < num_threads; t++) {
        thread_counts[t] = (size_t *)calloc(num_buckets, sizeof(size_t));
        if (!thread_counts[t]) {
            for (int j = 0; j < t; j++) free(thread_counts[j]);
            free(thread_counts); free(tree); free(temp);
            return STATA_ERR_MEMORY;
        }
    }

    /* Parallel counting */
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        size_t chunk = (nobs + num_threads - 1) / num_threads;
        size_t start = tid * chunk;
        size_t end = start + chunk;
        if (end > nobs) end = nobs;

        size_t *local_counts = thread_counts[tid];

        for (size_t i = start; i < end; i++) {
            if (i + IPS4O_PREFETCH_DISTANCE < end) {
                IPS4O_PREFETCH(&keys[order[i + IPS4O_PREFETCH_DISTANCE]]);
            }
            local_counts[ips4o_classify(keys[order[i]], tree, log_buckets, num_buckets)]++;
        }
    }

    /* Merge counts */
    size_t *bucket_sizes = (size_t *)calloc(num_buckets, sizeof(size_t));
    if (!bucket_sizes) {
        for (int t = 0; t < num_threads; t++) free(thread_counts[t]);
        free(thread_counts); free(tree); free(temp);
        return STATA_ERR_MEMORY;
    }

    for (int b = 0; b < num_buckets; b++) {
        for (int t = 0; t < num_threads; t++) {
            bucket_sizes[b] += thread_counts[t][b];
        }
    }

    /* Compute offsets */
    size_t *bucket_offsets = (size_t *)malloc((num_buckets + 1) * sizeof(size_t));
    if (!bucket_offsets) {
        free(bucket_sizes);
        for (int t = 0; t < num_threads; t++) free(thread_counts[t]);
        free(thread_counts); free(tree); free(temp);
        return STATA_ERR_MEMORY;
    }

    bucket_offsets[0] = 0;
    for (int b = 1; b <= num_buckets; b++) {
        bucket_offsets[b] = bucket_offsets[b - 1] + bucket_sizes[b - 1];
    }

    /* Compute per-thread write positions */
    size_t **thread_offsets = (size_t **)malloc(num_threads * sizeof(size_t *));
    if (!thread_offsets) {
        free(bucket_offsets); free(bucket_sizes);
        for (int t = 0; t < num_threads; t++) free(thread_counts[t]);
        free(thread_counts); free(tree); free(temp);
        return STATA_ERR_MEMORY;
    }

    for (int t = 0; t < num_threads; t++) {
        thread_offsets[t] = (size_t *)malloc(num_buckets * sizeof(size_t));
        if (!thread_offsets[t]) {
            for (int j = 0; j < t; j++) free(thread_offsets[j]);
            free(thread_offsets); free(bucket_offsets); free(bucket_sizes);
            for (int j = 0; j < num_threads; j++) free(thread_counts[j]);
            free(thread_counts); free(tree); free(temp);
            return STATA_ERR_MEMORY;
        }
    }

    for (int b = 0; b < num_buckets; b++) {
        size_t offset = bucket_offsets[b];
        for (int t = 0; t < num_threads; t++) {
            thread_offsets[t][b] = offset;
            offset += thread_counts[t][b];
        }
    }

    /* Parallel scatter */
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        size_t chunk = (nobs + num_threads - 1) / num_threads;
        size_t start = tid * chunk;
        size_t end = start + chunk;
        if (end > nobs) end = nobs;

        /* Thread-local copy of write positions */
        size_t local_offsets[IPS4O_MAX_BUCKETS];
        memcpy(local_offsets, thread_offsets[tid], num_buckets * sizeof(size_t));

        for (size_t i = start; i < end; i++) {
            if (i + IPS4O_PREFETCH_DISTANCE < end) {
                IPS4O_PREFETCH(&keys[order[i + IPS4O_PREFETCH_DISTANCE]]);
            }
            int b = ips4o_classify(keys[order[i]], tree, log_buckets, num_buckets);
            temp[local_offsets[b]++] = order[i];
        }
    }

    /* Copy back */
    memcpy(order, temp, nobs * sizeof(perm_idx_t));

    /* Cleanup partition resources */
    free(tree);
    for (int t = 0; t < num_threads; t++) {
        free(thread_counts[t]);
        free(thread_offsets[t]);
    }
    free(thread_counts);
    free(thread_offsets);

    /* Pre-allocate per-thread temp buffers for bucket sorting */
    size_t max_bucket = 0;
    for (int b = 0; b < num_buckets; b++) {
        if (bucket_sizes[b] > max_bucket) max_bucket = bucket_sizes[b];
    }

    perm_idx_t **bucket_temps = (perm_idx_t **)malloc(num_threads * sizeof(perm_idx_t *));
    if (!bucket_temps) {
        free(bucket_sizes); free(bucket_offsets); free(temp);
        return STATA_ERR_MEMORY;
    }

    for (int t = 0; t < num_threads; t++) {
        bucket_temps[t] = (perm_idx_t *)malloc(max_bucket * sizeof(perm_idx_t));
        if (!bucket_temps[t]) {
            for (int j = 0; j < t; j++) free(bucket_temps[j]);
            free(bucket_temps); free(bucket_sizes); free(bucket_offsets); free(temp);
            return STATA_ERR_MEMORY;
        }
    }

    /* Parallel bucket sorting using dynamic scheduling */
    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (int b = 0; b < num_buckets; b++) {
        if (bucket_sizes[b] > 1) {
            int tid = omp_get_thread_num();
            ips4o_sequential_numeric(order, bucket_temps[tid], keys,
                                     bucket_offsets[b], bucket_sizes[b], 0);
        }
    }

    /* Cleanup */
    for (int t = 0; t < num_threads; t++) {
        free(bucket_temps[t]);
    }
    free(bucket_temps);
    free(bucket_sizes);
    free(bucket_offsets);
    free(temp);

    return STATA_OK;
}

static stata_retcode ips4o_parallel_string(perm_idx_t * IPS4O_RESTRICT order,
                                           char ** IPS4O_RESTRICT strings,
                                           size_t nobs, size_t max_len)
{
    int num_threads = omp_get_max_threads();
    if (num_threads > 16) num_threads = 16;

    perm_idx_t *temp = (perm_idx_t *)ctools_aligned_alloc(64, nobs * sizeof(perm_idx_t));
    if (!temp) return STATA_ERR_MEMORY;

    int log_buckets = 0;
    while ((1 << log_buckets) < num_threads) log_buckets++;
    if (log_buckets > IPS4O_LOG_BUCKETS) log_buckets = IPS4O_LOG_BUCKETS;
    int num_buckets = 1 << log_buckets;

    size_t num_samples = IPS4O_OVERSAMPLE_FACTOR * num_buckets;
    if (num_samples > nobs) num_samples = nobs;

    char **samples = (char **)malloc(num_samples * sizeof(char *));
    if (!samples) { free(temp); return STATA_ERR_MEMORY; }

    size_t step = nobs / num_samples;
    if (step == 0) step = 1;
    for (size_t i = 0; i < num_samples; i++) {
        samples[i] = strings[order[i * step]];
    }
    qsort(samples, num_samples, sizeof(char *), compare_string_ptr);

    char **splitters = (char **)malloc((num_buckets - 1) * sizeof(char *));
    if (!splitters) { free(samples); free(temp); return STATA_ERR_MEMORY; }

    size_t splitter_step = num_samples / num_buckets;
    for (int i = 0; i < num_buckets - 1; i++) {
        splitters[i] = samples[(i + 1) * splitter_step];
    }
    free(samples);

    char **tree = (char **)calloc(2 * num_buckets, sizeof(char *));
    if (!tree) { free(splitters); free(temp); return STATA_ERR_MEMORY; }
    build_tree_string(tree, splitters, 1, 0, num_buckets - 2);
    free(splitters);

    /* Serial partition for strings (strcmp is expensive, parallelism less beneficial) */
    size_t *bucket_sizes = (size_t *)calloc(num_buckets, sizeof(size_t));
    size_t *bucket_offsets = (size_t *)malloc((num_buckets + 1) * sizeof(size_t));
    if (!bucket_sizes || !bucket_offsets) {
        free(tree); free(bucket_sizes); free(bucket_offsets); free(temp);
        return STATA_ERR_MEMORY;
    }

    for (size_t i = 0; i < nobs; i++) {
        bucket_sizes[ips4o_classify_string(strings[order[i]], tree, log_buckets, num_buckets)]++;
    }

    bucket_offsets[0] = 0;
    for (int b = 1; b <= num_buckets; b++) {
        bucket_offsets[b] = bucket_offsets[b - 1] + bucket_sizes[b - 1];
    }

    size_t *write_pos = (size_t *)malloc(num_buckets * sizeof(size_t));
    if (!write_pos) {
        free(tree); free(bucket_sizes); free(bucket_offsets); free(temp);
        return STATA_ERR_MEMORY;
    }
    memcpy(write_pos, bucket_offsets, num_buckets * sizeof(size_t));

    for (size_t i = 0; i < nobs; i++) {
        int b = ips4o_classify_string(strings[order[i]], tree, log_buckets, num_buckets);
        temp[write_pos[b]++] = order[i];
    }
    memcpy(order, temp, nobs * sizeof(perm_idx_t));

    free(tree);
    free(write_pos);

    /* Pre-allocate per-thread temp buffers */
    size_t max_bucket = 0;
    for (int b = 0; b < num_buckets; b++) {
        if (bucket_sizes[b] > max_bucket) max_bucket = bucket_sizes[b];
    }

    perm_idx_t **bucket_temps = (perm_idx_t **)malloc(num_threads * sizeof(perm_idx_t *));
    if (!bucket_temps) {
        free(bucket_sizes); free(bucket_offsets); free(temp);
        return STATA_ERR_MEMORY;
    }

    for (int t = 0; t < num_threads; t++) {
        bucket_temps[t] = (perm_idx_t *)malloc(max_bucket * sizeof(perm_idx_t));
        if (!bucket_temps[t]) {
            for (int j = 0; j < t; j++) free(bucket_temps[j]);
            free(bucket_temps); free(bucket_sizes); free(bucket_offsets); free(temp);
            return STATA_ERR_MEMORY;
        }
    }

    /* Parallel bucket sorting */
    #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (int b = 0; b < num_buckets; b++) {
        if (bucket_sizes[b] > 1) {
            int tid = omp_get_thread_num();
            ips4o_sequential_string(order, bucket_temps[tid], strings,
                                    bucket_offsets[b], bucket_sizes[b], max_len, 0);
        }
    }

    for (int t = 0; t < num_threads; t++) {
        free(bucket_temps[t]);
    }
    free(bucket_temps);
    free(bucket_sizes);
    free(bucket_offsets);
    free(temp);

    return STATA_OK;
}

#endif /* _OPENMP */

/* ============================================================================
   Variable-Level Sort Functions
   ============================================================================ */

static stata_retcode ips4o_sort_by_numeric_var(stata_data *data, int var_idx)
{
    double *dbl_data = data->vars[var_idx].data.dbl;
    size_t nobs = data->nobs;

    /* Allocate keys */
    uint64_t *keys = (uint64_t *)ctools_aligned_alloc(64, nobs * sizeof(uint64_t));
    if (!keys) return STATA_ERR_MEMORY;

    /* Parallel key conversion */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (size_t i = 0; i < nobs; i++) {
        keys[i] = ctools_double_to_sortable(dbl_data[i], SF_is_missing(dbl_data[i]));
    }

    stata_retcode rc;

#ifdef _OPENMP
    if (nobs >= IPS4O_PARALLEL_THRESHOLD) {
        rc = ips4o_parallel_numeric(data->sort_order, keys, nobs);
    } else {
        perm_idx_t *temp = (perm_idx_t *)malloc(nobs * sizeof(perm_idx_t));
        if (!temp) { ctools_aligned_free(keys); return STATA_ERR_MEMORY; }
        ips4o_sequential_numeric(data->sort_order, temp, keys, 0, nobs, 0);
        free(temp);
        rc = STATA_OK;
    }
#else
    perm_idx_t *temp = (perm_idx_t *)malloc(nobs * sizeof(perm_idx_t));
    if (!temp) { ctools_aligned_free(keys); return STATA_ERR_MEMORY; }
    ips4o_sequential_numeric(data->sort_order, temp, keys, 0, nobs, 0);
    free(temp);
    rc = STATA_OK;
#endif

    ctools_aligned_free(keys);
    return rc;
}

static stata_retcode ips4o_sort_by_string_var(stata_data *data, int var_idx)
{
    char **strings = data->vars[var_idx].data.str;
    size_t nobs = data->nobs;

    /* Calculate max string length */
    size_t max_len = 0;
    #ifdef _OPENMP
    #pragma omp parallel
    {
        size_t local_max = 0;
        #pragma omp for nowait
        for (size_t i = 0; i < nobs; i++) {
            size_t len = strlen(strings[i]);
            if (len > local_max) local_max = len;
        }
        #pragma omp critical
        {
            if (local_max > max_len) max_len = local_max;
        }
    }
    #else
    for (size_t i = 0; i < nobs; i++) {
        size_t len = strlen(strings[i]);
        if (len > max_len) max_len = len;
    }
    #endif

    stata_retcode rc;

#ifdef _OPENMP
    if (nobs >= IPS4O_PARALLEL_THRESHOLD) {
        rc = ips4o_parallel_string(data->sort_order, strings, nobs, max_len);
    } else {
        perm_idx_t *temp = (perm_idx_t *)malloc(nobs * sizeof(perm_idx_t));
        if (!temp) return STATA_ERR_MEMORY;
        ips4o_sequential_string(data->sort_order, temp, strings, 0, nobs, max_len, 0);
        free(temp);
        rc = STATA_OK;
    }
#else
    perm_idx_t *temp = (perm_idx_t *)malloc(nobs * sizeof(perm_idx_t));
    if (!temp) return STATA_ERR_MEMORY;
    ips4o_sequential_string(data->sort_order, temp, strings, 0, nobs, max_len, 0);
    free(temp);
    rc = STATA_OK;
#endif

    return rc;
}

/* ============================================================================
   Permutation Application - Parallel
   ============================================================================ */

static stata_retcode ips4o_apply_permutation(stata_data *data)
{
    size_t nvars = data->nvars;
    size_t nobs = data->nobs;
    perm_idx_t *perm = data->sort_order;

    if (nvars == 0 || nobs == 0) return STATA_OK;

#ifdef _OPENMP
    int success = 1;

    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t j = 0; j < nvars; j++) {
        if (!success) continue;

        stata_variable *var = &data->vars[j];

        if (var->type == STATA_TYPE_DOUBLE) {
            double *old_data = var->data.dbl;
            double *new_data = (double *)ctools_aligned_alloc(64, nobs * sizeof(double));
            if (!new_data) {
                #pragma omp atomic write
                success = 0;
                continue;
            }

            for (size_t i = 0; i < nobs; i++) {
                new_data[i] = old_data[perm[i]];
            }

            ctools_aligned_free(old_data);
            var->data.dbl = new_data;
        } else {
            char **old_data = var->data.str;
            char **new_data = (char **)ctools_aligned_alloc(CACHE_LINE_SIZE, nobs * sizeof(char *));
            if (!new_data) {
                #pragma omp atomic write
                success = 0;
                continue;
            }

            for (size_t i = 0; i < nobs; i++) {
                new_data[i] = old_data[perm[i]];
            }

            ctools_aligned_free(old_data);
            var->data.str = new_data;
        }
    }

    if (!success) return STATA_ERR_MEMORY;
#else
    for (size_t j = 0; j < nvars; j++) {
        stata_variable *var = &data->vars[j];

        if (var->type == STATA_TYPE_DOUBLE) {
            double *old_data = var->data.dbl;
            double *new_data = (double *)ctools_aligned_alloc(64, nobs * sizeof(double));
            if (!new_data) return STATA_ERR_MEMORY;

            for (size_t i = 0; i < nobs; i++) {
                new_data[i] = old_data[perm[i]];
            }

            ctools_aligned_free(old_data);
            var->data.dbl = new_data;
        } else {
            char **old_data = var->data.str;
            char **new_data = (char **)ctools_aligned_alloc(CACHE_LINE_SIZE, nobs * sizeof(char *));
            if (!new_data) return STATA_ERR_MEMORY;

            for (size_t i = 0; i < nobs; i++) {
                new_data[i] = old_data[perm[i]];
            }

            ctools_aligned_free(old_data);
            var->data.str = new_data;
        }
    }
#endif

    /* Reset sort_order to identity */
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_t j = 0; j < nobs; j++) {
        perm[j] = (perm_idx_t)j;
    }

    return STATA_OK;
}

/* ============================================================================
   Public API
   ============================================================================ */

stata_retcode ctools_sort_ips4o(stata_data *data, int *sort_vars, size_t nsort)
{
    if (!data || !sort_vars || data->nobs == 0 || nsort == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Sort from last key to first for stable multi-key sort */
    for (int k = (int)nsort - 1; k >= 0; k--) {
        int var_idx = sort_vars[k] - 1;

        if (var_idx < 0 || var_idx >= (int)data->nvars) {
            return STATA_ERR_INVALID_INPUT;
        }

        stata_retcode rc;
        if (data->vars[var_idx].type == STATA_TYPE_DOUBLE) {
            rc = ips4o_sort_by_numeric_var(data, var_idx);
        } else {
            rc = ips4o_sort_by_string_var(data, var_idx);
        }

        if (rc != STATA_OK) return rc;
    }

    return ips4o_apply_permutation(data);
}

/*
    ctools_sort_ips4o_with_perm - IPS4o sort with permutation output

    Same as ctools_sort_ips4o but captures the permutation mapping before
    applying it. This is needed by cmerge to apply the same permutation to
    keepusing variables that are stored separately from the key variables.

    @param data       Dataset to sort (modified in place)
    @param sort_vars  Array of 1-based variable indices to sort by
    @param nsort      Number of sort variables
    @param perm_out   Output array for permutation (must be pre-allocated with nobs elements)
                      Maps sorted index -> original index

    @return STATA_OK on success, or error code
*/
stata_retcode ctools_sort_ips4o_with_perm(stata_data *data, int *sort_vars,
                                           size_t nsort, size_t *perm_out)
{
    if (!data || !sort_vars || data->nobs == 0 || nsort == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Sort from last key to first for stable multi-key sort */
    for (int k = (int)nsort - 1; k >= 0; k--) {
        int var_idx = sort_vars[k] - 1;

        if (var_idx < 0 || var_idx >= (int)data->nvars) {
            return STATA_ERR_INVALID_INPUT;
        }

        stata_retcode rc;
        if (data->vars[var_idx].type == STATA_TYPE_DOUBLE) {
            rc = ips4o_sort_by_numeric_var(data, var_idx);
        } else {
            rc = ips4o_sort_by_string_var(data, var_idx);
        }

        if (rc != STATA_OK) return rc;
    }

    /* If caller wants the permutation, copy it BEFORE apply_permutation resets it */
    if (perm_out != NULL) {
        size_t nobs = data->nobs;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t i = 0; i < nobs; i++) {
            perm_out[i] = data->sort_order[i];
        }
    }

    /* Apply the permutation to all variables (this resets sort_order to identity) */
    return ips4o_apply_permutation(data);
}

/*
    ctools_sort_ips4o_order_only - IPS4o sort without applying permutation

    Computes sort_order but does NOT apply the permutation to data.
    After this call, data->sort_order contains the permutation but data is unchanged.
    Call ctools_apply_permutation() separately to apply the permutation.

    This is useful when you need to time sort vs permutation separately,
    or when you want to inspect the permutation before applying it.

    @param data       Dataset to sort
    @param sort_vars  Array of 1-based variable indices to sort by
    @param nsort      Number of sort variables

    @return STATA_OK on success, or error code
*/
stata_retcode ctools_sort_ips4o_order_only(stata_data *data, int *sort_vars, size_t nsort)
{
    if (!data || !sort_vars || data->nobs == 0 || nsort == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Sort from last key to first for stable multi-key sort */
    for (int k = (int)nsort - 1; k >= 0; k--) {
        int var_idx = sort_vars[k] - 1;

        if (var_idx < 0 || var_idx >= (int)data->nvars) {
            return STATA_ERR_INVALID_INPUT;
        }

        stata_retcode rc;
        if (data->vars[var_idx].type == STATA_TYPE_DOUBLE) {
            rc = ips4o_sort_by_numeric_var(data, var_idx);
        } else {
            rc = ips4o_sort_by_string_var(data, var_idx);
        }

        if (rc != STATA_OK) return rc;
    }

    /* Do NOT apply permutation - caller will do it separately */
    return STATA_OK;
}

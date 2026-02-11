/*
    crangestat_impl.c
    High-performance range statistics for Stata

    Performance optimizations:
    - Parallel IPS4o sort (via ctools infrastructure)
    - Binary search for O(log n) window finding
    - Run-length optimization: same-key observations share results
    - Single-pass statistics without data copying
    - SIMD-accelerated window statistics
    - Small group batching to reduce parallel overhead
    - OpenMP parallelization across observations/runs

    Author: ctools project
    License: MIT
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_simd.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "crangestat_impl.h"

#define CRANGESTAT_MODULE "crangestat"

/* Thresholds for optimizations */
#define CRANGESTAT_RUN_OPT_THRESHOLD 1000    /* Min group size for run optimization */
#define CRANGESTAT_BATCH_THRESHOLD 100       /* Max group size for batching */
#define CRANGESTAT_PARALLEL_THRESHOLD 10000  /* Min work for parallelization */
#define CRANGESTAT_PARALLEL_GROUPS_MIN 8     /* Min groups for cross-group parallelism */
#define CRANGESTAT_SMALL_GROUP_THRESHOLD 500 /* Max avg group size for cross-group parallelism */
#define CRANGESTAT_PREFIX_SUM_THRESHOLD 64   /* Min window size to use prefix sums */
#define CRANGESTAT_SPARSE_TABLE_LOG 20       /* Max log2(n) for sparse table */

/* Statistics identifiers */
typedef enum {
    STAT_COUNT = 0,
    STAT_MEAN,
    STAT_SUM,
    STAT_MIN,
    STAT_MAX,
    STAT_SD,
    STAT_VARIANCE,
    STAT_MEDIAN,
    STAT_IQR,
    STAT_FIRST,
    STAT_LAST,
    STAT_FIRSTNM,
    STAT_LASTNM,
    STAT_P1, STAT_P5, STAT_P10, STAT_P25, STAT_P75, STAT_P90, STAT_P95, STAT_P99,
    STAT_SKEWNESS,
    STAT_KURTOSIS,
    STAT_COUNT_ALL  /* Total number of statistic types */
} stat_type;

/* Configuration */
typedef struct {
    int exclude_self;
    int verbose;
    double interval_low;   /* Offset for low bound (can be . for -inf) */
    double interval_high;  /* Offset for high bound (can be . for +inf) */
    int low_is_missing;    /* 1 if low bound is . (meaning -infinity) */
    int high_is_missing;   /* 1 if high bound is . (meaning +infinity) */
} crangestat_config;

/* Statistic specification */
typedef struct {
    stat_type type;
    int source_var_idx;   /* Index into loaded source variables */
    int result_var_idx;   /* 1-based Stata variable index for result */
} stat_spec;

/* ===========================================================================
   Binary search for window boundaries
   =========================================================================== */

/* Find first index where key[i] >= target */
static size_t lower_bound(const double *key, size_t start, size_t end, double target)
{
    while (start < end) {
        size_t mid = start + (end - start) / 2;
        if (key[mid] < target) {
            start = mid + 1;
        } else {
            end = mid;
        }
    }
    return start;
}

/* Find first index where key[i] > target */
static size_t upper_bound(const double *key, size_t start, size_t end, double target)
{
    while (start < end) {
        size_t mid = start + (end - start) / 2;
        if (key[mid] <= target) {
            start = mid + 1;
        } else {
            end = mid;
        }
    }
    return start;
}

/* ===========================================================================
   Quickselect for median/percentiles
   =========================================================================== */

static void swap_double(double *a, double *b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}

static size_t partition_double(double *arr, size_t left, size_t right)
{
    /* Median of three pivot selection */
    size_t mid = left + (right - left) / 2;
    if (arr[mid] < arr[left]) swap_double(&arr[left], &arr[mid]);
    if (arr[right] < arr[left]) swap_double(&arr[left], &arr[right]);
    if (arr[mid] < arr[right]) swap_double(&arr[mid], &arr[right]);
    double pivot = arr[right];

    size_t i = left;
    for (size_t j = left; j < right; j++) {
        if (arr[j] <= pivot) {
            swap_double(&arr[i], &arr[j]);
            i++;
        }
    }
    swap_double(&arr[i], &arr[right]);
    return i;
}

static double quickselect(double *arr, size_t n, size_t k)
{
    if (n == 0) return SV_missval;
    if (n == 1) return arr[0];
    if (k >= n) k = n - 1;

    size_t left = 0;
    size_t right = n - 1;

    while (left < right) {
        size_t pivot_idx = partition_double(arr, left, right);
        if (k == pivot_idx) {
            return arr[k];
        } else if (k < pivot_idx) {
            right = pivot_idx - 1;
        } else {
            left = pivot_idx + 1;
        }
    }
    return arr[left];
}

static double compute_percentile(double *work, size_t n, double pctl)
{
    if (n == 0) return SV_missval;
    if (n == 1) return work[0];

    double pos = (pctl / 100.0) * (double)(n - 1);
    size_t k = (size_t)floor(pos);
    if (k >= n) k = n - 1;

    double val = quickselect(work, n, k);

    double frac = pos - (double)k;
    if (frac > 0.0 && k + 1 < n) {
        /* Find min of remaining elements > val */
        double next_val = DBL_MAX;
        for (size_t i = k + 1; i < n; i++) {
            if (work[i] < next_val) next_val = work[i];
        }
        if (next_val < DBL_MAX) {
            val = val + frac * (next_val - val);
        }
    }
    return val;
}

/* ===========================================================================
   Group detection
   =========================================================================== */

typedef struct {
    size_t start;
    size_t count;
} group_info;

static void detect_groups(double **by_data, size_t nobs, size_t nby,
                          group_info *groups, size_t *ngroups)
{
    if (nobs == 0) {
        *ngroups = 0;
        return;
    }

    const double miss = SV_missval;
    size_t g = 0;
    groups[0].start = 0;
    groups[0].count = 1;

    for (size_t i = 1; i < nobs; i++) {
        int same_group = 1;

        for (size_t b = 0; b < nby && same_group; b++) {
            double prev = by_data[b][i - 1];
            double curr = by_data[b][i];

            int prev_miss = (prev >= miss);
            int curr_miss = (curr >= miss);

            if (prev_miss && curr_miss) continue;
            if (prev_miss || curr_miss || prev != curr) {
                same_group = 0;
            }
        }

        if (same_group) {
            groups[g].count++;
        } else {
            g++;
            groups[g].start = i;
            groups[g].count = 1;
        }
    }

    *ngroups = g + 1;
}

/* ===========================================================================
   Prefix Sum Arrays for O(1) Range Queries

   For sum, count, mean, sd, variance, we precompute prefix sums:
   - prefix_sum[i] = sum of data[0..i-1] (exclusive end)
   - prefix_sum2[i] = sum of data[0..i-1]^2
   - prefix_count[i] = count of non-missing in [0..i-1]

   Range query: sum[a..b) = prefix[b] - prefix[a]
   With excludeself: subtract self contribution
   =========================================================================== */

typedef struct {
    double *prefix_sum;    /* Prefix sums: prefix_sum[i] = sum of [0,i) */
    double *prefix_sum2;   /* Prefix sum of squares */
    size_t *prefix_count;  /* Prefix count of non-missing */
    size_t nobs;
} prefix_arrays;

/* Build prefix arrays for a source variable within a group */
static void build_prefix_arrays(const double *data, size_t start, size_t count,
                                prefix_arrays *pa)
{
    const double miss = SV_missval;

    /* prefix[0] = 0 (empty prefix) */
    pa->prefix_sum[0] = 0.0;
    pa->prefix_sum2[0] = 0.0;
    pa->prefix_count[0] = 0;

    /*
     * OPTIMIZATION: Parallel prefix sum using work-efficient parallel scan
     *
     * For large groups, we parallelize the prefix sum computation:
     * 1. Divide into chunks, compute local sums per chunk (parallel)
     * 2. Compute prefix of chunk totals (sequential, but small)
     * 3. Add chunk prefix to each element (parallel)
     *
     * For small groups, sequential is faster due to less overhead.
     */
    #ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    #else
    int nthreads = 1;
    #endif

    size_t parallel_threshold = 10000;

    if (count >= parallel_threshold && nthreads > 1) {
        /* Parallel prefix sum */
        size_t chunk_size = (count + nthreads - 1) / nthreads;

        /* Per-chunk totals — stack-allocated (nthreads is bounded by hardware) */
        double chunk_sum[64];
        double chunk_sum2[64];
        size_t chunk_count[64];

        if (nthreads > 64) {
            goto sequential_prefix;
        }

        /* Phase 1: Compute local prefix sums within each chunk (parallel) */
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            size_t my_start = tid * chunk_size;
            size_t my_end = (my_start + chunk_size < count) ? my_start + chunk_size : count;

            double local_sum = 0.0, local_sum2 = 0.0;
            size_t local_count = 0;

            for (size_t i = my_start; i < my_end; i++) {
                double val = data[start + i];
                size_t idx = i + 1;

                if (val < miss) {
                    local_sum += val;
                    local_sum2 += val * val;
                    local_count++;
                }

                pa->prefix_sum[idx] = local_sum;
                pa->prefix_sum2[idx] = local_sum2;
                pa->prefix_count[idx] = local_count;
            }

            chunk_sum[tid] = local_sum;
            chunk_sum2[tid] = local_sum2;
            chunk_count[tid] = local_count;
        }

        /* Phase 2: Compute prefix of chunk totals (sequential) */
        double prefix_chunk_sum = 0.0, prefix_chunk_sum2 = 0.0;
        size_t prefix_chunk_count = 0;

        for (int t = 0; t < nthreads; t++) {
            double old_sum = chunk_sum[t];
            double old_sum2 = chunk_sum2[t];
            size_t old_count = chunk_count[t];

            chunk_sum[t] = prefix_chunk_sum;
            chunk_sum2[t] = prefix_chunk_sum2;
            chunk_count[t] = prefix_chunk_count;

            prefix_chunk_sum += old_sum;
            prefix_chunk_sum2 += old_sum2;
            prefix_chunk_count += old_count;
        }

        /* Phase 3: Add chunk prefix to each element (parallel) */
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            if (tid > 0) {  /* Thread 0's chunk already correct */
                size_t my_start = tid * chunk_size;
                size_t my_end = (my_start + chunk_size < count) ? my_start + chunk_size : count;

                double add_sum = chunk_sum[tid];
                double add_sum2 = chunk_sum2[tid];
                size_t add_count = chunk_count[tid];

                for (size_t i = my_start; i < my_end; i++) {
                    size_t idx = i + 1;
                    pa->prefix_sum[idx] += add_sum;
                    pa->prefix_sum2[idx] += add_sum2;
                    pa->prefix_count[idx] += add_count;
                }
            }
        }

        /* chunk_sum, chunk_sum2, chunk_count are stack-allocated — no free needed */
    } else {
sequential_prefix:;
        /* Sequential prefix sum for small groups */
        double sum = 0.0, sum2 = 0.0;
        size_t cnt = 0;

        for (size_t i = 0; i < count; i++) {
            double val = data[start + i];
            size_t idx = i + 1;

            if (val < miss) {
                sum += val;
                sum2 += val * val;
                cnt++;
            }

            pa->prefix_sum[idx] = sum;
            pa->prefix_sum2[idx] = sum2;
            pa->prefix_count[idx] = cnt;
        }
    }
}

/* Query sum for range [win_start, win_end) within group starting at g_start */
static inline void query_prefix_sums(const prefix_arrays *pa, size_t g_start,
                                     size_t win_start, size_t win_end,
                                     double *out_sum, double *out_sum2, size_t *out_count)
{
    size_t a = win_start - g_start;  /* Convert to group-relative indices */
    size_t b = win_end - g_start;

    *out_sum = pa->prefix_sum[b] - pa->prefix_sum[a];
    *out_sum2 = pa->prefix_sum2[b] - pa->prefix_sum2[a];
    *out_count = pa->prefix_count[b] - pa->prefix_count[a];
}

/* ===========================================================================
   Sparse Table for O(1) Range Min/Max Queries

   sparse_min[i][j] = min of range [i, i + 2^j)
   sparse_max[i][j] = max of range [i, i + 2^j)

   Query [l, r]: k = log2(r - l), answer = min(sparse[l][k], sparse[r - 2^k][k])
   =========================================================================== */

typedef struct {
    double **sparse_min;   /* sparse_min[i][j] = min of [i, i + 2^j) */
    double **sparse_max;   /* sparse_max[i][j] = max of [i, i + 2^j) */
    int *log_table;        /* Precomputed log2 values */
    size_t nobs;
    int max_log;
} sparse_table;

/* Build sparse table for range min/max queries */
static int build_sparse_table(const double *data, size_t start, size_t count,
                              sparse_table *st, const double miss)
{
    if (count == 0) return 0;

    /* Compute max_log = floor(log2(count)) */
    int max_log = 0;
    size_t temp = count;
    while (temp > 1) {
        temp >>= 1;
        max_log++;
    }
    if (max_log >= CRANGESTAT_SPARSE_TABLE_LOG) {
        max_log = CRANGESTAT_SPARSE_TABLE_LOG - 1;
    }
    st->max_log = max_log;
    st->nobs = count;

    /* Allocate log table */
    st->log_table = (int *)malloc((count + 1) * sizeof(int));
    if (!st->log_table) return -1;

    st->log_table[0] = 0;
    st->log_table[1] = 0;
    for (size_t i = 2; i <= count; i++) {
        st->log_table[i] = st->log_table[i / 2] + 1;
    }

    /* Allocate pointer arrays for sparse tables */
    int nlevels = max_log + 1;
    st->sparse_min = (double **)malloc(nlevels * sizeof(double *));
    st->sparse_max = (double **)malloc(nlevels * sizeof(double *));
    if (!st->sparse_min || !st->sparse_max) {
        free(st->log_table);
        if (st->sparse_min) free(st->sparse_min);
        if (st->sparse_max) free(st->sparse_max);
        return -1;
    }

    /* Single contiguous allocation for all sparse table data:
     * 2 * nlevels * count doubles (min and max for each level) */
    double *block = (double *)malloc(2 * (size_t)nlevels * count * sizeof(double));
    if (!block) {
        free(st->sparse_min);
        free(st->sparse_max);
        free(st->log_table);
        return -1;
    }

    for (int j = 0; j < nlevels; j++) {
        st->sparse_min[j] = block + (size_t)j * count;
        st->sparse_max[j] = block + (size_t)(nlevels + j) * count;
    }

    /* Initialize level 0: individual elements */
    for (size_t i = 0; i < count; i++) {
        double val = data[start + i];
        if (val >= miss) {
            /* Missing values: use sentinel values */
            st->sparse_min[0][i] = DBL_MAX;
            st->sparse_max[0][i] = -DBL_MAX;
        } else {
            st->sparse_min[0][i] = val;
            st->sparse_max[0][i] = val;
        }
    }

    /* Build higher levels */
    for (int j = 1; j <= max_log; j++) {
        size_t range_size = (size_t)1 << j;
        size_t half_range = (size_t)1 << (j - 1);
        size_t limit = (range_size <= count) ? count - range_size + 1 : 0;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) if(count > 10000)
        #endif
        for (size_t i = 0; i < limit; i++) {
            double min1 = st->sparse_min[j-1][i];
            double min2 = st->sparse_min[j-1][i + half_range];
            double max1 = st->sparse_max[j-1][i];
            double max2 = st->sparse_max[j-1][i + half_range];

            st->sparse_min[j][i] = (min1 < min2) ? min1 : min2;
            st->sparse_max[j][i] = (max1 > max2) ? max1 : max2;
        }
    }

    return 0;
}

/* Query range min for [l, r) - l and r are group-relative */
static inline double query_range_min(const sparse_table *st, size_t l, size_t r)
{
    if (l >= r) return DBL_MAX;
    size_t len = r - l;
    int k = st->log_table[len];
    size_t range = (size_t)1 << k;

    double min1 = st->sparse_min[k][l];
    double min2 = st->sparse_min[k][r - range];

    return (min1 < min2) ? min1 : min2;
}

/* Query range max for [l, r) - l and r are group-relative */
static inline double query_range_max(const sparse_table *st, size_t l, size_t r)
{
    if (l >= r) return -DBL_MAX;
    size_t len = r - l;
    int k = st->log_table[len];
    size_t range = (size_t)1 << k;

    double max1 = st->sparse_max[k][l];
    double max2 = st->sparse_max[k][r - range];

    return (max1 > max2) ? max1 : max2;
}

/* Free sparse table */
static void free_sparse_table(sparse_table *st)
{
    if (st->sparse_min) {
        /* All level data is in a single contiguous block starting at sparse_min[0] */
        free(st->sparse_min[0]);
        free(st->sparse_min);
    }
    if (st->sparse_max) {
        /* sparse_max pointers point into the same block; just free the pointer array */
        free(st->sparse_max);
    }
    if (st->log_table) free(st->log_table);
}

/* Allocate sparse table buffers without filling data.
 * max_count is the maximum group size; actual groups may be smaller.
 * Returns 0 on success, -1 on allocation failure. */
static int alloc_sparse_table(sparse_table *st, size_t max_count)
{
    if (max_count == 0) return 0;

    /* Compute max_log for the largest possible group */
    int max_log = 0;
    size_t temp = max_count;
    while (temp > 1) {
        temp >>= 1;
        max_log++;
    }
    if (max_log >= CRANGESTAT_SPARSE_TABLE_LOG)
        max_log = CRANGESTAT_SPARSE_TABLE_LOG - 1;

    st->max_log = max_log;
    st->nobs = max_count;

    /* Allocate log table for [0, max_count] */
    st->log_table = (int *)malloc((max_count + 1) * sizeof(int));
    if (!st->log_table) return -1;

    /* Pre-fill log table */
    st->log_table[0] = 0;
    if (max_count >= 1) st->log_table[1] = 0;
    for (size_t i = 2; i <= max_count; i++) {
        st->log_table[i] = st->log_table[i / 2] + 1;
    }

    /* Allocate pointer arrays */
    int nlevels = max_log + 1;
    st->sparse_min = (double **)malloc(nlevels * sizeof(double *));
    st->sparse_max = (double **)malloc(nlevels * sizeof(double *));
    if (!st->sparse_min || !st->sparse_max) {
        free(st->log_table);
        if (st->sparse_min) free(st->sparse_min);
        if (st->sparse_max) free(st->sparse_max);
        st->log_table = NULL;
        st->sparse_min = NULL;
        st->sparse_max = NULL;
        return -1;
    }

    /* Single contiguous allocation: stride = max_count per level */
    double *block = (double *)malloc(2 * (size_t)nlevels * max_count * sizeof(double));
    if (!block) {
        free(st->sparse_min);
        free(st->sparse_max);
        free(st->log_table);
        st->sparse_min = NULL;
        st->sparse_max = NULL;
        st->log_table = NULL;
        return -1;
    }

    for (int j = 0; j < nlevels; j++) {
        st->sparse_min[j] = block + (size_t)j * max_count;
        st->sparse_max[j] = block + (size_t)(nlevels + j) * max_count;
    }

    return 0;
}

/* Rebuild sparse table data for a specific group using pre-allocated buffers.
 * st must have been allocated with alloc_sparse_table(st, max_count) where
 * max_count >= count. Pointer stride remains max_count (the allocation capacity).
 * Only count entries per level are meaningful; queries must use indices < count. */
static void rebuild_sparse_table(sparse_table *st, const double *data,
                                 size_t start, size_t count, double miss)
{
    if (count == 0) {
        st->nobs = 0;
        return;
    }

    /* Recompute max_log for this group's count */
    int max_log = 0;
    size_t temp = count;
    while (temp > 1) {
        temp >>= 1;
        max_log++;
    }
    if (max_log >= CRANGESTAT_SPARSE_TABLE_LOG)
        max_log = CRANGESTAT_SPARSE_TABLE_LOG - 1;

    st->max_log = max_log;
    st->nobs = count;

    /* Rebuild log table for [0, count] */
    st->log_table[0] = 0;
    if (count >= 1) st->log_table[1] = 0;
    for (size_t i = 2; i <= count; i++) {
        st->log_table[i] = st->log_table[i / 2] + 1;
    }

    /* Initialize level 0 */
    for (size_t i = 0; i < count; i++) {
        double val = data[start + i];
        if (val >= miss) {
            st->sparse_min[0][i] = DBL_MAX;
            st->sparse_max[0][i] = -DBL_MAX;
        } else {
            st->sparse_min[0][i] = val;
            st->sparse_max[0][i] = val;
        }
    }

    /* Build higher levels */
    for (int j = 1; j <= max_log; j++) {
        size_t range_size = (size_t)1 << j;
        size_t half_range = (size_t)1 << (j - 1);
        size_t limit = (range_size <= count) ? count - range_size + 1 : 0;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) if(count > 10000)
        #endif
        for (size_t i = 0; i < limit; i++) {
            double min1 = st->sparse_min[j-1][i];
            double min2 = st->sparse_min[j-1][i + half_range];
            double max1 = st->sparse_max[j-1][i];
            double max2 = st->sparse_max[j-1][i + half_range];

            st->sparse_min[j][i] = (min1 < min2) ? min1 : min2;
            st->sparse_max[j][i] = (max1 > max2) ? max1 : max2;
        }
    }
}

/* ===========================================================================
   Statistics computation for a window
   =========================================================================== */

/* Check if stat type needs work array (percentiles, etc.) */
static int stat_needs_work_array(stat_type type)
{
    switch (type) {
        case STAT_MEDIAN:
        case STAT_IQR:
        case STAT_P1:
        case STAT_P5:
        case STAT_P10:
        case STAT_P25:
        case STAT_P75:
        case STAT_P90:
        case STAT_P95:
        case STAT_P99:
            return 1;
        default:
            return 0;
    }
}

/* Check if stat can use O(1) prefix sum queries */
static int stat_uses_prefix_sums(stat_type type)
{
    switch (type) {
        case STAT_COUNT:
        case STAT_SUM:
        case STAT_MEAN:
        case STAT_SD:
        case STAT_VARIANCE:
            return 1;
        default:
            return 0;
    }
}

/* Check if stat can use O(1) sparse table queries */
static int stat_uses_sparse_table(stat_type type)
{
    switch (type) {
        case STAT_MIN:
        case STAT_MAX:
            return 1;
        default:
            return 0;
    }
}

/*
    compute_stat_with_prefix - O(1) stat computation using prefix sums

    For count, sum, mean, sd, variance, this is O(1) instead of O(window_size)!
    This is a massive speedup for large windows.
*/
static double compute_stat_with_prefix(stat_type type, const prefix_arrays *pa,
                                       const double *data, size_t g_start,
                                       size_t win_start, size_t win_end,
                                       size_t self_idx, int exclude_self)
{
    double sum, sum2;
    size_t n;

    /* Query prefix sums - O(1) */
    query_prefix_sums(pa, g_start, win_start, win_end, &sum, &sum2, &n);

    /* Handle excludeself by subtracting self contribution */
    if (exclude_self && self_idx >= win_start && self_idx < win_end) {
        double self_val = data[self_idx];
        if (self_val < SV_missval) {
            sum -= self_val;
            sum2 -= self_val * self_val;
            n -= 1;
        }
    }

    if (n == 0) {
        /* For count, 0 valid observations = count of 0, not missing */
        if (type == STAT_COUNT) return 0.0;
        return SV_missval;
    }

    switch (type) {
        case STAT_COUNT:
            return (double)n;

        case STAT_SUM:
            return sum;

        case STAT_MEAN:
            return sum / (double)n;

        case STAT_SD:
        case STAT_VARIANCE: {
            if (n < 2) return SV_missval;
            double mean = sum / (double)n;
            double var = (sum2 - sum * mean) / (double)(n - 1);
            if (var < 0.0) var = 0.0;  /* Numerical stability */
            return (type == STAT_SD) ? sqrt(var) : var;
        }

        default:
            return SV_missval;
    }
}

/*
    compute_minmax_with_sparse - O(1) min/max using sparse table

    For min/max without excludeself, this is O(1) instead of O(window_size)!
    With excludeself, we check if result equals self and rescan if needed.
*/
static double compute_minmax_with_sparse(stat_type type, const sparse_table *st,
                                         const double *data, size_t g_start,
                                         size_t win_start, size_t win_end,
                                         size_t self_idx, int exclude_self)
{
    if (win_start >= win_end) return SV_missval;

    /* Convert to group-relative indices */
    size_t l = win_start - g_start;
    size_t r = win_end - g_start;

    double result;
    if (type == STAT_MIN) {
        result = query_range_min(st, l, r);
        if (result >= DBL_MAX - 1.0) return SV_missval;  /* All missing */
    } else {
        result = query_range_max(st, l, r);
        if (result <= -DBL_MAX + 1.0) return SV_missval;  /* All missing */
    }

    /* Handle excludeself */
    if (exclude_self && self_idx >= win_start && self_idx < win_end) {
        double self_val = data[self_idx];
        if (self_val < SV_missval) {
            if ((type == STAT_MIN && self_val == result) ||
                (type == STAT_MAX && self_val == result)) {
                /*
                    Self equals the min/max - need to find second-best.
                    Query left and right of self separately.
                */
                size_t self_rel = self_idx - g_start;
                double left_result, right_result;

                if (type == STAT_MIN) {
                    left_result = (self_rel > l) ? query_range_min(st, l, self_rel) : DBL_MAX;
                    right_result = (self_rel + 1 < r) ? query_range_min(st, self_rel + 1, r) : DBL_MAX;
                    result = (left_result < right_result) ? left_result : right_result;
                    if (result >= DBL_MAX - 1.0) return SV_missval;
                } else {
                    left_result = (self_rel > l) ? query_range_max(st, l, self_rel) : -DBL_MAX;
                    right_result = (self_rel + 1 < r) ? query_range_max(st, self_rel + 1, r) : -DBL_MAX;
                    result = (left_result > right_result) ? left_result : right_result;
                    if (result <= -DBL_MAX + 1.0) return SV_missval;
                }
            }
        }
    }

    return result;
}

/* SIMD-accelerated sum for window statistics */
static inline void simd_sum_minmax(const double *data, size_t start, size_t end,
                                   double *out_sum, double *out_sum2,
                                   double *out_min, double *out_max, size_t *out_n)
{
    const double miss = SV_missval;
    double sum = 0.0, sum2 = 0.0;
    double minval = DBL_MAX, maxval = -DBL_MAX;
    size_t n = 0;
    size_t i = start;

    #if CTOOLS_HAS_AVX2
    /* AVX2: process 4 doubles at a time */
    if (end - start >= 8) {
        __m256d vsum = _mm256_setzero_pd();
        __m256d vsum2 = _mm256_setzero_pd();
        __m256d vmin = _mm256_set1_pd(DBL_MAX);
        __m256d vmax = _mm256_set1_pd(-DBL_MAX);
        __m256d vmiss = _mm256_set1_pd(miss);

        for (; i + 4 <= end; i += 4) {
            __m256d v = _mm256_loadu_pd(&data[i]);
            __m256d mask = ctools_valid_mask_avx2(v, vmiss);

            /* Masked operations */
            __m256d v_masked = _mm256_and_pd(v, mask);
            vsum = _mm256_add_pd(vsum, v_masked);
            vsum2 = _mm256_add_pd(vsum2, _mm256_mul_pd(v_masked, v_masked));

            /* Min/max with mask */
            __m256d vmin_cand = _mm256_blendv_pd(vmin, v, mask);
            __m256d vmax_cand = _mm256_blendv_pd(vmax, v, mask);
            vmin = _mm256_min_pd(vmin, vmin_cand);
            vmax = _mm256_max_pd(vmax, vmax_cand);

            /* Count non-missing */
            n += (size_t)ctools_count_valid_avx2(mask);
        }

        /* Horizontal reduction */
        sum = ctools_hsum_avx2(vsum);
        sum2 = ctools_hsum_avx2(vsum2);
        minval = ctools_hmin_avx2(vmin);
        maxval = ctools_hmax_avx2(vmax);
    }
    #elif CTOOLS_HAS_NEON
    /* NEON: process 2 doubles at a time */
    if (end - start >= 4) {
        float64x2_t vsum = vdupq_n_f64(0.0);
        float64x2_t vsum2 = vdupq_n_f64(0.0);
        float64x2_t vmin = vdupq_n_f64(DBL_MAX);
        float64x2_t vmax = vdupq_n_f64(-DBL_MAX);
        float64x2_t vmiss = vdupq_n_f64(miss);

        for (; i + 2 <= end; i += 2) {
            float64x2_t v = vld1q_f64(&data[i]);
            uint64x2_t mask = ctools_valid_mask_neon(v, vmiss);

            /* Count non-missing */
            n += ctools_count_valid_neon(mask);

            /* Masked operations using bit select */
            float64x2_t zero = vdupq_n_f64(0.0);
            float64x2_t v_masked = vbslq_f64(mask, v, zero);
            vsum = vaddq_f64(vsum, v_masked);
            vsum2 = vaddq_f64(vsum2, vmulq_f64(v_masked, v_masked));

            /* Min/max */
            float64x2_t vmin_cand = vbslq_f64(mask, v, vmin);
            float64x2_t vmax_cand = vbslq_f64(mask, v, vmax);
            vmin = vminq_f64(vmin, vmin_cand);
            vmax = vmaxq_f64(vmax, vmax_cand);
        }

        /* Horizontal reduction */
        sum = ctools_hsum_neon(vsum);
        sum2 = ctools_hsum_neon(vsum2);
        minval = ctools_hmin_neon(vmin);
        maxval = ctools_hmax_neon(vmax);
    }
    #endif

    /* Scalar tail */
    for (; i < end; i++) {
        double val = data[i];
        if (val < miss) {
            sum += val;
            sum2 += val * val;
            if (val < minval) minval = val;
            if (val > maxval) maxval = val;
            n++;
        }
    }

    *out_sum = sum;
    *out_sum2 = sum2;
    *out_min = minval;
    *out_max = maxval;
    *out_n = n;
}

/* Compute simple statistics directly without copying to work array */
static double compute_window_stat_simple(stat_type type, const double *data,
                                         size_t win_start, size_t win_end,
                                         size_t self_idx, int exclude_self)
{
    const double miss = SV_missval;
    size_t n = 0;
    double sum = 0.0, sum2 = 0.0;
    double minval = DBL_MAX, maxval = -DBL_MAX;
    double first_val = miss, last_val = miss;

    /* Check if self is in window and non-missing (for excludeself optimization) */
    int self_in_window = 0;
    int self_is_missing = 1;
    double self_val = 0.0;
    if (exclude_self && self_idx >= win_start && self_idx < win_end) {
        self_in_window = 1;
        self_val = data[self_idx];
        self_is_missing = (self_val >= miss);
    }

    /* Use SIMD for large windows */
    if ((win_end - win_start) >= 8) {
        simd_sum_minmax(data, win_start, win_end, &sum, &sum2, &minval, &maxval, &n);

        /* SIMD-then-subtract optimization for excludeself */
        if (exclude_self && self_in_window && !self_is_missing) {
            /* Subtract self's contribution from aggregates */
            n -= 1;
            sum -= self_val;
            sum2 -= self_val * self_val;

            /* For min/max: if self equals min/max, need to rescan */
            /* Only rescan if we actually need min or max */
            if (type == STAT_MIN && self_val == minval) {
                minval = DBL_MAX;
                for (size_t i = win_start; i < win_end; i++) {
                    if (i == self_idx) continue;
                    double val = data[i];
                    if (val < miss && val < minval) minval = val;
                }
            }
            if (type == STAT_MAX && self_val == maxval) {
                maxval = -DBL_MAX;
                for (size_t i = win_start; i < win_end; i++) {
                    if (i == self_idx) continue;
                    double val = data[i];
                    if (val < miss && val > maxval) maxval = val;
                }
            }
        }

        /* Find first and last non-missing (still need sequential scan) */
        if (type == STAT_FIRST || type == STAT_FIRSTNM ||
            type == STAT_LAST || type == STAT_LASTNM) {
            for (size_t i = win_start; i < win_end; i++) {
                if (exclude_self && i == self_idx) continue;
                if (data[i] < miss) {
                    first_val = data[i];
                    break;
                }
            }
            for (size_t i = win_end; i > win_start; i--) {
                if (exclude_self && (i - 1) == self_idx) continue;
                if (data[i - 1] < miss) {
                    last_val = data[i - 1];
                    break;
                }
            }
        }
    } else {
        /* Scalar path for small windows */
        for (size_t i = win_start; i < win_end; i++) {
            if (exclude_self && i == self_idx) continue;
            double val = data[i];
            if (val < miss) {
                if (n == 0) first_val = val;
                last_val = val;
                sum += val;
                sum2 += val * val;
                if (val < minval) minval = val;
                if (val > maxval) maxval = val;
                n++;
            }
        }
    }

    if (n == 0) return SV_missval;

    switch (type) {
        case STAT_COUNT:
            return (double)n;

        case STAT_SUM:
            return sum;

        case STAT_MEAN:
            return sum / (double)n;

        case STAT_MIN:
            return minval;

        case STAT_MAX:
            return maxval;

        case STAT_SD:
        case STAT_VARIANCE: {
            if (n < 2) return SV_missval;
            double mean = sum / (double)n;
            double var = (sum2 - sum * mean) / (double)(n - 1);
            if (var < 0.0) var = 0.0;
            return (type == STAT_SD) ? sqrt(var) : var;
        }

        case STAT_FIRST:
        case STAT_FIRSTNM:
            return first_val;

        case STAT_LAST:
        case STAT_LASTNM:
            return last_val;

        case STAT_SKEWNESS: {
            if (n < 3) return SV_missval;
            double mean = sum / (double)n;
            /* Need second pass for higher moments */
            double m2 = 0.0, m3 = 0.0;
            for (size_t i = win_start; i < win_end; i++) {
                if (exclude_self && i == self_idx) continue;
                double val = data[i];
                if (val < miss) {
                    double d = val - mean;
                    m2 += d * d;
                    m3 += d * d * d;
                }
            }
            m2 /= (double)n;
            m3 /= (double)n;
            if (m2 <= 0.0) return SV_missval;
            double sd = sqrt(m2);
            return m3 / (sd * sd * sd);
        }

        case STAT_KURTOSIS: {
            if (n < 4) return SV_missval;
            double mean = sum / (double)n;
            /* Need second pass for higher moments */
            double m2 = 0.0, m4 = 0.0;
            for (size_t i = win_start; i < win_end; i++) {
                if (exclude_self && i == self_idx) continue;
                double val = data[i];
                if (val < miss) {
                    double d = val - mean;
                    double d2 = d * d;
                    m2 += d2;
                    m4 += d2 * d2;
                }
            }
            m2 /= (double)n;
            m4 /= (double)n;
            if (m2 <= 0.0) return SV_missval;
            return m4 / (m2 * m2);
        }

        default:
            return SV_missval;
    }
}

/* Compute statistics that need work array (percentiles, etc.) */
static double compute_window_stat_complex(stat_type type, const double *data,
                                          size_t win_start, size_t win_end,
                                          size_t self_idx, int exclude_self,
                                          double *work)
{
    const double miss = SV_missval;

    /* Extract non-missing values into work array */
    size_t n = 0;
    for (size_t i = win_start; i < win_end; i++) {
        if (exclude_self && i == self_idx) continue;
        double val = data[i];
        if (val < miss) {
            work[n++] = val;
        }
    }

    if (n == 0) return SV_missval;

    switch (type) {
        case STAT_MEDIAN:
            return compute_percentile(work, n, 50.0);

        case STAT_IQR: {
            double p25 = compute_percentile(work, n, 25.0);
            double p75 = compute_percentile(work, n, 75.0);
            if (p25 >= miss || p75 >= miss) return SV_missval;
            return p75 - p25;
        }

        case STAT_P1:  return compute_percentile(work, n, 1.0);
        case STAT_P5:  return compute_percentile(work, n, 5.0);
        case STAT_P10: return compute_percentile(work, n, 10.0);
        case STAT_P25: return compute_percentile(work, n, 25.0);
        case STAT_P75: return compute_percentile(work, n, 75.0);
        case STAT_P90: return compute_percentile(work, n, 90.0);
        case STAT_P95: return compute_percentile(work, n, 95.0);
        case STAT_P99: return compute_percentile(work, n, 99.0);

        default:
            return SV_missval;
    }
}

/* ===========================================================================
   Argument parsing
   =========================================================================== */

static int parse_option(const char *args, const char *name)
{
    if (!args || !name) return 0;
    const char *p = strstr(args, name);
    if (!p) return 0;
    if (p > args && *(p-1) != ' ' && *(p-1) != '\t') return 0;
    char after = *(p + strlen(name));
    if (after != '\0' && after != ' ' && after != '\t' && after != '=') return 0;
    return 1;
}

static double parse_double_option(const char *args, const char *key, double def)
{
    if (!args || !key) return def;
    char pat[64];
    snprintf(pat, sizeof(pat), "%s=", key);
    const char *p = strstr(args, pat);
    if (!p) return def;

    /* Use strtod directly since value may be followed by other options */
    char *endptr;
    double result = strtod(p + strlen(pat), &endptr);
    if (endptr == p + strlen(pat)) {
        return def;  /* No conversion performed */
    }
    return result;
}

/* ===========================================================================
   Main Entry Point
   =========================================================================== */

ST_retcode crangestat_main(const char *args)
{
    double t_start, t_load, t_sort, t_groups, t_compute, t_store;
    int *source_indices = NULL;
    int *result_indices = NULL;
    int *by_indices = NULL;
    int *load_indices = NULL;
    double **source_data = NULL;
    double **result_data = NULL;
    double **by_data = NULL;
    double *key_data = NULL;
    double **work_arrays = NULL;
    stat_spec *specs = NULL;
    group_info *groups = NULL;
    size_t *inv_perm = NULL;
    perm_idx_t *fwd_perm = NULL;
    ctools_filtered_data filtered;
    perm_idx_t *obs_map = NULL;
    int num_threads = 1;

    t_start = ctools_timer_seconds();

    if (!args || strlen(args) == 0) {
        ctools_error(CRANGESTAT_MODULE, "no arguments specified");
        return 198;
    }

    /* Parse options */
    crangestat_config config;
    config.exclude_self = parse_option(args, "excludeself");
    config.verbose = parse_option(args, "verbose");
    config.interval_low = parse_double_option(args, "low", 0.0);
    config.interval_high = parse_double_option(args, "high", 0.0);
    config.low_is_missing = parse_option(args, "low_miss");
    config.high_is_missing = parse_option(args, "high_miss");

    /*
       Argument format:
       "nstats nsource nby key_idx source_indices... result_indices... by_indices... stat_types... [options]"

       stat_types are space-separated integers corresponding to stat_type enum
    */

    const char *p = args;
    while (*p == ' ' || *p == '\t') p++;

    char *end;
    size_t nstats = (size_t)strtol(p, &end, 10);
    if (end == p || nstats == 0) {
        ctools_error(CRANGESTAT_MODULE, "invalid nstats");
        return 198;
    }
    p = end;

    while (*p == ' ' || *p == '\t') p++;
    size_t nsource = (size_t)strtol(p, &end, 10);
    if (end == p || nsource == 0) {
        ctools_error(CRANGESTAT_MODULE, "invalid nsource");
        return 198;
    }
    p = end;

    while (*p == ' ' || *p == '\t') p++;
    size_t nby = (size_t)strtol(p, &end, 10);
    if (end == p) {
        ctools_error(CRANGESTAT_MODULE, "invalid nby");
        return 198;
    }
    p = end;

    while (*p == ' ' || *p == '\t') p++;
    int key_idx = (int)strtol(p, &end, 10);
    if (end == p || key_idx < 1) {
        ctools_error(CRANGESTAT_MODULE, "invalid key_idx");
        return 198;
    }
    p = end;

    /* Allocate arrays */
    source_indices = (int *)malloc(nsource * sizeof(int));
    result_indices = (int *)malloc(nstats * sizeof(int));
    specs = (stat_spec *)malloc(nstats * sizeof(stat_spec));

    if (!source_indices || !result_indices || !specs) {
        if (source_indices) free(source_indices);
        if (result_indices) free(result_indices);
        if (specs) free(specs);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    /* Parse source variable indices */
    for (size_t i = 0; i < nsource; i++) {
        while (*p == ' ' || *p == '\t') p++;
        source_indices[i] = (int)strtol(p, &end, 10);
        if (end == p) {
            free(source_indices);
            free(result_indices);
            free(specs);
            ctools_error(CRANGESTAT_MODULE, "failed to parse source indices");
            return 198;
        }
        p = end;
    }

    /* Parse result variable indices */
    for (size_t i = 0; i < nstats; i++) {
        while (*p == ' ' || *p == '\t') p++;
        result_indices[i] = (int)strtol(p, &end, 10);
        if (end == p) {
            free(source_indices);
            free(result_indices);
            free(specs);
            ctools_error(CRANGESTAT_MODULE, "failed to parse result indices");
            return 198;
        }
        p = end;
    }

    /* Parse by-variable indices if present */
    if (nby > 0) {
        by_indices = (int *)malloc(nby * sizeof(int));
        if (!by_indices) {
            free(source_indices);
            free(result_indices);
            free(specs);
            ctools_error_alloc(CRANGESTAT_MODULE);
            return 920;
        }

        for (size_t i = 0; i < nby; i++) {
            while (*p == ' ' || *p == '\t') p++;
            by_indices[i] = (int)strtol(p, &end, 10);
            if (end == p) {
                free(source_indices);
                free(result_indices);
                free(specs);
                free(by_indices);
                ctools_error(CRANGESTAT_MODULE, "failed to parse by indices");
                return 198;
            }
            p = end;
        }
    }

    /* Parse stat types (as a group) */
    for (size_t i = 0; i < nstats; i++) {
        while (*p == ' ' || *p == '\t') p++;
        int stat_int = (int)strtol(p, &end, 10);
        if (end == p) {
            free(source_indices);
            free(result_indices);
            free(specs);
            if (by_indices) free(by_indices);
            ctools_error(CRANGESTAT_MODULE, "failed to parse stat types");
            return 198;
        }
        p = end;
        specs[i].type = (stat_type)stat_int;
        specs[i].result_var_idx = result_indices[i];
    }

    /* Parse source var mappings (as a group) */
    for (size_t i = 0; i < nstats; i++) {
        while (*p == ' ' || *p == '\t') p++;
        int src_idx = (int)strtol(p, &end, 10);
        if (end == p) {
            free(source_indices);
            free(result_indices);
            free(specs);
            if (by_indices) free(by_indices);
            ctools_error(CRANGESTAT_MODULE, "failed to parse source var for stat");
            return 198;
        }
        p = end;
        specs[i].source_var_idx = src_idx;
    }

    #ifdef _OPENMP
    num_threads = omp_get_max_threads();
    #endif

    /*
     * === Load Phase ===
     * Use ctools_data_load for optimized parallel data loading with if/in filtering.
     *
     * Variable layout in loaded data:
     *   [0..nby-1]           = by-variables (for grouping)
     *   [nby]                = key variable (for window bounds)
     *   [nby+1..nby+nsource] = source variables (for stat computation)
     */
    double load_start = ctools_timer_seconds();

    size_t nvars_to_load = nby + 1 + nsource;
    load_indices = (int *)malloc(nvars_to_load * sizeof(int));
    if (!load_indices) {
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    /* Build load index array: [by_vars..., key, source_vars...] */
    size_t idx = 0;
    for (size_t b = 0; b < nby; b++) {
        load_indices[idx++] = by_indices[b];
    }
    load_indices[idx++] = key_idx;
    for (size_t v = 0; v < nsource; v++) {
        load_indices[idx++] = source_indices[v];
    }

    /* Load data using filtered loading (handles if/in at load time) */
    ctools_filtered_data_init(&filtered);

    stata_retcode rc = ctools_data_load(&filtered, load_indices, nvars_to_load, 0, 0, 0);
    if (rc != STATA_OK) {
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error(CRANGESTAT_MODULE, "failed to load data");
        return 920;
    }

    /* Get filtered observation count and obs_map */
    size_t nobs = filtered.data.nobs;
    obs_map = filtered.obs_map;

    if (nobs == 0) {
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error(CRANGESTAT_MODULE, "no observations");
        return 2000;
    }

    t_load = ctools_timer_seconds() - load_start;

    /*
     * === Sort Phase ===
     * Sort by (by_vars, key) using optimized algorithm selection.
     *
     * PERFORMANCE: Use SORT_ALG_COUNTING which:
     *   1. Tries counting sort first (O(n+k), optimal for small-range keys like dates)
     *   2. Falls back to LSD radix if range is too large (still O(n*d), faster than comparison sorts)
     *
     * This provides massive speedups for typical rangestat workloads where the key
     * variable (e.g., date, time period) has a small number of unique values.
     * Benchmark: 20M obs with 20 unique dates: IPS4o=3.1s, Counting=0.01s (280x faster)
     */
    double sort_start = ctools_timer_seconds();

    /* Build sort variable indices: [1, 2, ..., nby, nby+1] (1-based into loaded data) */
    size_t nsort = nby + 1;
    int *sort_vars = (int *)malloc(nsort * sizeof(int));
    if (!sort_vars) {
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    for (size_t i = 0; i < nsort; i++) {
        sort_vars[i] = (int)(i + 1);  /* 1-based indices */
    }

    /* Sort using counting sort (with LSD radix fallback for large ranges) */
    rc = ctools_sort_dispatch(&filtered.data, sort_vars, nsort, SORT_ALG_COUNTING);
    if (rc != STATA_OK) {
        free(sort_vars);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error(CRANGESTAT_MODULE, "sort failed");
        return 920;
    }

    /*
     * Compute inverse permutation before applying.
     * sort_order[sorted_idx] = original_idx
     * inv_perm[original_idx] = sorted_idx
     *
     * This lets us map results (in sorted order) back to original observation order.
     */
    inv_perm = (size_t *)malloc(nobs * sizeof(size_t));
    if (!inv_perm) {
        free(sort_vars);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    /* Compute inverse: inv_perm[original_idx] = sorted_idx */
    for (size_t i = 0; i < nobs; i++) {
        inv_perm[filtered.data.sort_order[i]] = i;
    }

    /* Save forward permutation for sequential store-phase reads.
     * sort_order[sorted_idx] = original_idx, so iterating sorted_idx gives
     * sequential reads from data[] and scattered writes via obs_map[]. */
    fwd_perm = (perm_idx_t *)malloc(nobs * sizeof(perm_idx_t));
    if (fwd_perm) {
        memcpy(fwd_perm, filtered.data.sort_order, nobs * sizeof(perm_idx_t));
    }
    /* Non-fatal: falls back to inv_perm-based store if alloc fails */

    /* Apply permutation to sort data in place */
    rc = ctools_apply_permutation(&filtered.data);
    if (rc != STATA_OK) {
        free(inv_perm);
        free(sort_vars);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error(CRANGESTAT_MODULE, "permutation failed");
        return 920;
    }

    free(sort_vars);

    t_sort = ctools_timer_seconds() - sort_start;

    /*
     * Set up pointers to sorted data for computation.
     * Data layout after loading:
     *   filtered.data.vars[0..nby-1] = by-variables
     *   filtered.data.vars[nby] = key variable
     *   filtered.data.vars[nby+1..nby+nsource] = source variables
     */
    key_data = filtered.data.vars[nby].data.dbl;

    /* Set up source_data pointers (don't allocate, just point into stata_data) */
    source_data = (double **)malloc(nsource * sizeof(double *));
    if (!source_data) {
        free(inv_perm);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }
    for (size_t v = 0; v < nsource; v++) {
        source_data[v] = filtered.data.vars[nby + 1 + v].data.dbl;
    }

    /* Set up by_data pointers */
    if (nby > 0) {
        by_data = (double **)malloc(nby * sizeof(double *));
        if (!by_data) {
            free(source_data);
            free(inv_perm);
            ctools_filtered_data_free(&filtered);
            free(load_indices);
            free(source_indices);
            free(result_indices);
            free(specs);
            if (by_indices) free(by_indices);
            ctools_error_alloc(CRANGESTAT_MODULE);
            return 920;
        }
        for (size_t b = 0; b < nby; b++) {
            by_data[b] = filtered.data.vars[b].data.dbl;
        }
    }

    /* Allocate result arrays (these are separate, not in stata_data) */
    result_data = (double **)malloc(nstats * sizeof(double *));
    if (!result_data) {
        if (by_data) free(by_data);
        free(source_data);
        free(inv_perm);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    for (size_t s = 0; s < nstats; s++) {
        result_data[s] = (double *)malloc(nobs * sizeof(double));
        if (!result_data[s]) {
            for (size_t j = 0; j < s; j++) free(result_data[j]);
            free(result_data);
            if (by_data) free(by_data);
            free(source_data);
            free(inv_perm);
            ctools_filtered_data_free(&filtered);
            free(load_indices);
            free(source_indices);
            free(result_indices);
            free(specs);
            if (by_indices) free(by_indices);
            ctools_error_alloc(CRANGESTAT_MODULE);
            return 920;
        }
        /* Initialize to missing */
        for (size_t i = 0; i < nobs; i++) {
            result_data[s][i] = SV_missval;
        }
    }

    /* === Group Detection Phase === */
    double groups_start = ctools_timer_seconds();

    size_t ngroups = 1;
    groups = (group_info *)malloc((nobs + 1) * sizeof(group_info));
    if (!groups) {
        free(inv_perm);
        for (size_t s = 0; s < nstats; s++) free(result_data[s]);
        free(result_data);
        if (by_data) free(by_data);
        free(source_data);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    if (nby > 0) {
        detect_groups(by_data, nobs, nby, groups, &ngroups);
    } else {
        groups[0].start = 0;
        groups[0].count = nobs;
        ngroups = 1;
    }

    /* Shrink groups array from (nobs+1) to actual ngroups.
     * For 10M rows with 1000 groups: 160MB -> 16KB. */
    if (ngroups < nobs + 1) {
        group_info *shrunk = (group_info *)realloc(groups, ngroups * sizeof(group_info));
        if (shrunk) groups = shrunk;  /* realloc failure is non-fatal */
    }

    t_groups = ctools_timer_seconds() - groups_start;

    /* === Computation Phase === */
    double compute_start = ctools_timer_seconds();

    /* Compute max group size for work array allocation */
    size_t max_group_size = 0;
    for (size_t g = 0; g < ngroups; g++) {
        if (groups[g].count > max_group_size)
            max_group_size = groups[g].count;
    }

    /* Allocate work arrays for each thread (sized to max_group_size, not nobs) */
    work_arrays = (double **)malloc(num_threads * sizeof(double *));
    if (!work_arrays) {
        free(groups);
        free(inv_perm);
        for (size_t s = 0; s < nstats; s++) free(result_data[s]);
        free(result_data);
        if (by_data) free(by_data);
        free(source_data);
        ctools_filtered_data_free(&filtered);
        free(load_indices);
        free(source_indices);
        free(result_indices);
        free(specs);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CRANGESTAT_MODULE);
        return 920;
    }

    for (int t = 0; t < num_threads; t++) {
        work_arrays[t] = (double *)malloc(max_group_size * sizeof(double));
        if (!work_arrays[t]) {
            for (int j = 0; j < t; j++) free(work_arrays[j]);
            free(work_arrays);
            free(groups);
            free(inv_perm);
            for (size_t s = 0; s < nstats; s++) free(result_data[s]);
            free(result_data);
            if (by_data) free(by_data);
            free(source_data);
            ctools_filtered_data_free(&filtered);
            free(load_indices);
            free(source_indices);
            free(result_indices);
            free(specs);
            if (by_indices) free(by_indices);
            ctools_error_alloc(CRANGESTAT_MODULE);
            return 920;
        }
    }

    /*
     * OPTIMIZATION: Cross-group parallelism for many small groups
     *
     * When there are many by-groups with small average size, parallelizing
     * across groups is more efficient than parallelizing within each group.
     * This avoids OpenMP region setup/teardown overhead per group.
     */
    size_t avg_group_size = nobs / ngroups;
    int use_cross_group_parallel = (ngroups >= CRANGESTAT_PARALLEL_GROUPS_MIN &&
                                    avg_group_size < CRANGESTAT_SMALL_GROUP_THRESHOLD);

    if (use_cross_group_parallel) {
        /*
         * Cross-group parallel processing: single parallel region over all groups.
         * Each thread processes complete groups independently.
         *
         * Pre-allocate per-thread O(1) buffers (prefix arrays reused across groups,
         * sparse tables built/freed per group since groups are small).
         */

        /* Determine which source variables need O(1) structures */
        int *a_needs_prefix = (int *)calloc(nsource, sizeof(int));
        int *a_needs_sparse = (int *)calloc(nsource, sizeof(int));
        int a_has_o1 = (a_needs_prefix && a_needs_sparse &&
                        max_group_size >= CRANGESTAT_PREFIX_SUM_THRESHOLD);

        if (a_has_o1) {
            for (size_t s = 0; s < nstats; s++) {
                int src_idx = specs[s].source_var_idx;
                if (src_idx < 0 || src_idx >= (int)nsource) continue;
                if (stat_uses_prefix_sums(specs[s].type)) a_needs_prefix[src_idx] = 1;
                if (stat_uses_sparse_table(specs[s].type)) a_needs_sparse[src_idx] = 1;
            }
        }

        /* Per-thread prefix arrays: thread t, source var v → a_prefix[t * nsource + v] */
        prefix_arrays *a_prefix = NULL;
        if (a_has_o1) {
            a_prefix = (prefix_arrays *)calloc((size_t)num_threads * nsource, sizeof(prefix_arrays));
            if (a_prefix) {
                for (int t = 0; t < num_threads; t++) {
                    for (size_t v = 0; v < nsource; v++) {
                        if (!a_needs_prefix[v]) continue;
                        size_t idx = (size_t)t * nsource + v;
                        a_prefix[idx].prefix_sum = (double *)malloc((max_group_size + 1) * sizeof(double));
                        a_prefix[idx].prefix_sum2 = (double *)malloc((max_group_size + 1) * sizeof(double));
                        a_prefix[idx].prefix_count = (size_t *)malloc((max_group_size + 1) * sizeof(size_t));
                        if (!a_prefix[idx].prefix_sum || !a_prefix[idx].prefix_sum2 ||
                            !a_prefix[idx].prefix_count) {
                            /* Disable prefix for this var on alloc failure */
                            a_needs_prefix[v] = 0;
                            if (a_prefix[idx].prefix_sum) { free(a_prefix[idx].prefix_sum); a_prefix[idx].prefix_sum = NULL; }
                            if (a_prefix[idx].prefix_sum2) { free(a_prefix[idx].prefix_sum2); a_prefix[idx].prefix_sum2 = NULL; }
                            if (a_prefix[idx].prefix_count) { free(a_prefix[idx].prefix_count); a_prefix[idx].prefix_count = NULL; }
                        }
                    }
                }
            } else {
                a_has_o1 = 0;
            }
        }

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1)
        #endif
        for (size_t g = 0; g < ngroups; g++) {
            #ifdef _OPENMP
            int tid = omp_get_thread_num();
            #else
            int tid = 0;
            #endif

            double *work = work_arrays[tid];
            size_t g_start = groups[g].start;
            size_t g_count = groups[g].count;
            size_t g_end = g_start + g_count;

            /* Build O(1) structures for this group */
            int group_has_o1 = a_has_o1 && (g_count >= CRANGESTAT_PREFIX_SUM_THRESHOLD);
            sparse_table group_sparse[16]; /* Stack-allocated for small nsource */
            sparse_table *group_sparse_ptr = (nsource <= 16) ? group_sparse : NULL;
            int sparse_built = 0;

            if (group_has_o1) {
                /* Build prefix arrays (reuse per-thread buffers) */
                for (size_t v = 0; v < nsource; v++) {
                    if (a_needs_prefix[v]) {
                        size_t idx = (size_t)tid * nsource + v;
                        a_prefix[idx].nobs = g_count;
                        build_prefix_arrays(source_data[v], g_start, g_count, &a_prefix[idx]);
                    }
                }

                /* Build sparse tables per group (groups are small, so alloc overhead is minimal) */
                if (group_sparse_ptr) {
                    memset(group_sparse_ptr, 0, nsource * sizeof(sparse_table));
                    for (size_t v = 0; v < nsource; v++) {
                        if (a_needs_sparse[v]) {
                            if (build_sparse_table(source_data[v], g_start, g_count,
                                                   &group_sparse_ptr[v], SV_missval) != 0) {
                                /* Build failed, disable for this group */
                                group_sparse_ptr = NULL;
                                break;
                            }
                        }
                    }
                    if (group_sparse_ptr) sparse_built = 1;
                }
            }

            /* Process each observation in this group with sliding window.
             * Since data is sorted by key, window bounds only advance forward. */
            size_t win_lo = g_start;
            size_t win_hi = g_start;
            for (size_t i = g_start; i < g_end; i++) {
                double key_val = key_data[i];

                /* Missing keys sort to end — stop processing */
                if (key_val >= SV_missval) break;

                /* Compute window bounds */
                double low_bound = config.low_is_missing ? -DBL_MAX : key_val + config.interval_low;
                double high_bound = config.high_is_missing ? DBL_MAX : key_val + config.interval_high;

                /* Slide window pointers forward (O(1) amortized) */
                while (win_lo < g_end && key_data[win_lo] < low_bound) win_lo++;
                while (win_hi < g_end && key_data[win_hi] <= high_bound) win_hi++;

                size_t win_start = win_lo;
                size_t win_end = win_hi;

                /* Compute each statistic for this observation */
                for (size_t s = 0; s < nstats; s++) {
                    int src_idx = specs[s].source_var_idx;
                    if (src_idx < 0 || src_idx >= (int)nsource) continue;

                    double val;
                    stat_type stype = specs[s].type;

                    if (stat_needs_work_array(stype)) {
                        val = compute_window_stat_complex(
                            stype,
                            source_data[src_idx],
                            win_start,
                            win_end,
                            i,
                            config.exclude_self,
                            work
                        );
                    } else if (group_has_o1 && a_needs_prefix[src_idx] &&
                               stat_uses_prefix_sums(stype)) {
                        size_t pidx = (size_t)tid * nsource + (size_t)src_idx;
                        val = compute_stat_with_prefix(
                            stype,
                            &a_prefix[pidx],
                            source_data[src_idx],
                            g_start,
                            win_start,
                            win_end,
                            i,
                            config.exclude_self
                        );
                    } else if (group_has_o1 && sparse_built && a_needs_sparse[src_idx] &&
                               stat_uses_sparse_table(stype)) {
                        val = compute_minmax_with_sparse(
                            stype,
                            &group_sparse_ptr[src_idx],
                            source_data[src_idx],
                            g_start,
                            win_start,
                            win_end,
                            i,
                            config.exclude_self
                        );
                    } else {
                        val = compute_window_stat_simple(
                            stype,
                            source_data[src_idx],
                            win_start,
                            win_end,
                            i,
                            config.exclude_self
                        );
                    }

                    result_data[s][i] = val;
                }
            }

            /* Free per-group sparse tables */
            if (sparse_built) {
                for (size_t v = 0; v < nsource; v++) {
                    if (a_needs_sparse[v]) {
                        free_sparse_table(&group_sparse_ptr[v]);
                    }
                }
            }
        }

        /* Cleanup per-thread prefix arrays */
        if (a_prefix) {
            for (int t = 0; t < num_threads; t++) {
                for (size_t v = 0; v < nsource; v++) {
                    size_t idx = (size_t)t * nsource + v;
                    if (a_prefix[idx].prefix_sum) free(a_prefix[idx].prefix_sum);
                    if (a_prefix[idx].prefix_sum2) free(a_prefix[idx].prefix_sum2);
                    if (a_prefix[idx].prefix_count) free(a_prefix[idx].prefix_count);
                }
            }
            free(a_prefix);
        }
        if (a_needs_prefix) free(a_needs_prefix);
        if (a_needs_sparse) free(a_needs_sparse);
    } else {
        /* Standard path: process groups sequentially with internal parallelism */

        /*
         * OPTIMIZATION: O(1) range queries using prefix sums and sparse tables
         *
         * Pre-allocate O(1) data structures once at max_group_size, then
         * rebuild per group to avoid per-group alloc/free overhead.
         */
        int use_o1_queries = (max_group_size >= CRANGESTAT_PREFIX_SUM_THRESHOLD);

        /* Determine which source variables need which data structures */
        int *needs_prefix = NULL;
        int *needs_sparse = NULL;
        prefix_arrays *prefix_arr = NULL;
        sparse_table *sparse_tbl = NULL;

        if (use_o1_queries) {
            needs_prefix = (int *)calloc(nsource, sizeof(int));
            needs_sparse = (int *)calloc(nsource, sizeof(int));
            prefix_arr = (prefix_arrays *)calloc(nsource, sizeof(prefix_arrays));
            sparse_tbl = (sparse_table *)calloc(nsource, sizeof(sparse_table));

            if (!needs_prefix || !needs_sparse || !prefix_arr || !sparse_tbl) {
                use_o1_queries = 0;
                if (needs_prefix) free(needs_prefix);
                if (needs_sparse) free(needs_sparse);
                if (prefix_arr) free(prefix_arr);
                if (sparse_tbl) free(sparse_tbl);
                needs_prefix = NULL;
                needs_sparse = NULL;
                prefix_arr = NULL;
                sparse_tbl = NULL;
            }
        }

        if (use_o1_queries) {
            for (size_t s = 0; s < nstats; s++) {
                int src_idx = specs[s].source_var_idx;
                if (src_idx < 0 || src_idx >= (int)nsource) continue;

                if (stat_uses_prefix_sums(specs[s].type)) {
                    needs_prefix[src_idx] = 1;
                }
                if (stat_uses_sparse_table(specs[s].type)) {
                    needs_sparse[src_idx] = 1;
                }
            }

            /* Pre-allocate buffers at max_group_size */
            for (size_t v = 0; v < nsource; v++) {
                if (needs_prefix[v]) {
                    prefix_arr[v].prefix_sum = (double *)malloc((max_group_size + 1) * sizeof(double));
                    prefix_arr[v].prefix_sum2 = (double *)malloc((max_group_size + 1) * sizeof(double));
                    prefix_arr[v].prefix_count = (size_t *)malloc((max_group_size + 1) * sizeof(size_t));
                    prefix_arr[v].nobs = max_group_size;

                    if (!prefix_arr[v].prefix_sum || !prefix_arr[v].prefix_sum2 ||
                        !prefix_arr[v].prefix_count) {
                        needs_prefix[v] = 0;
                        if (prefix_arr[v].prefix_sum) free(prefix_arr[v].prefix_sum);
                        if (prefix_arr[v].prefix_sum2) free(prefix_arr[v].prefix_sum2);
                        if (prefix_arr[v].prefix_count) free(prefix_arr[v].prefix_count);
                        prefix_arr[v].prefix_sum = NULL;
                        prefix_arr[v].prefix_sum2 = NULL;
                        prefix_arr[v].prefix_count = NULL;
                    }
                }

                if (needs_sparse[v]) {
                    if (alloc_sparse_table(&sparse_tbl[v], max_group_size) != 0) {
                        needs_sparse[v] = 0;
                    }
                }
            }
        }

        for (size_t g = 0; g < ngroups; g++) {
            size_t g_start = groups[g].start;
            size_t g_count = groups[g].count;
            size_t g_end = g_start + g_count;

            /*
             * OPTIMIZATION: When excludeself is false, observations with the same
             * key value get the same window and same result. We process in "runs"
             * of same-key observations, computing once per unique key.
             *
             * This can provide massive speedups when there are many duplicate keys.
             * For example, panel data with 100K IDs and 20 time periods:
             * - Without optimization: 2M window computations
             * - With optimization: 20 window computations
             *
             * For small groups (< 1000 obs), skip run detection as overhead exceeds benefit.
             */
            if (!config.exclude_self && g_count >= CRANGESTAT_RUN_OPT_THRESHOLD) {
                /* First pass: identify run boundaries (starts of each unique key) */
                size_t *run_starts = (size_t *)malloc((g_count + 1) * sizeof(size_t));
                if (!run_starts) {
                    /* Fall back to standard path on allocation failure */
                    goto standard_obs_loop;
                }

                size_t num_runs = 0;
                run_starts[num_runs++] = g_start;

                for (size_t i = g_start + 1; i < g_end; i++) {
                    double prev_key = key_data[i - 1];
                    double curr_key = key_data[i];

                    /* Check for run boundary (different key or missing) */
                    int prev_miss = (prev_key >= SV_missval);
                    int curr_miss = (curr_key >= SV_missval);

                    if (prev_miss != curr_miss || (!prev_miss && prev_key != curr_key)) {
                        run_starts[num_runs++] = i;
                    }
                }
                run_starts[num_runs] = g_end;  /* Sentinel */

                /* Build O(1) structures for this group (reuse pre-allocated buffers) */
                int group_has_o1 = use_o1_queries && (g_count >= CRANGESTAT_PREFIX_SUM_THRESHOLD);
                if (group_has_o1) {
                    for (size_t v = 0; v < nsource; v++) {
                        if (needs_prefix[v]) {
                            prefix_arr[v].nobs = g_count;
                            build_prefix_arrays(source_data[v], g_start, g_count, &prefix_arr[v]);
                        }
                        if (needs_sparse[v]) {
                            rebuild_sparse_table(&sparse_tbl[v], source_data[v], g_start,
                                                 g_count, SV_missval);
                        }
                    }
                }

                /* Process runs in parallel */
                #ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic, 1)
                #endif
                for (size_t r = 0; r < num_runs; r++) {
                    #ifdef _OPENMP
                    int tid = omp_get_thread_num();
                    #else
                    int tid = 0;
                    #endif

                    double *work = work_arrays[tid];
                    size_t run_start = run_starts[r];
                    size_t run_end = run_starts[r + 1];
                    double key_val = key_data[run_start];

                    /* Skip missing keys */
                    if (key_val >= SV_missval) {
                        continue;
                    }

                    /* Compute window bounds */
                    double low_bound = config.low_is_missing ? -DBL_MAX : key_val + config.interval_low;
                    double high_bound = config.high_is_missing ? DBL_MAX : key_val + config.interval_high;

                    /* Find window using binary search within group */
                    size_t win_start = lower_bound(key_data, g_start, g_end, low_bound);
                    size_t win_end = upper_bound(key_data, g_start, g_end, high_bound);

                    /* Compute each statistic ONCE for this run */
                    for (size_t s = 0; s < nstats; s++) {
                        int src_idx = specs[s].source_var_idx;
                        if (src_idx < 0 || src_idx >= (int)nsource) continue;

                        double val;
                        stat_type stype = specs[s].type;

                        if (stat_needs_work_array(stype)) {
                            val = compute_window_stat_complex(
                                stype,
                                source_data[src_idx],
                                win_start,
                                win_end,
                                0,  /* self_idx unused when exclude_self=0 */
                                0,  /* exclude_self=0 */
                                work
                            );
                        } else if (group_has_o1 && needs_prefix[src_idx] &&
                                   stat_uses_prefix_sums(stype)) {
                            val = compute_stat_with_prefix(
                                stype,
                                &prefix_arr[src_idx],
                                source_data[src_idx],
                                g_start,
                                win_start,
                                win_end,
                                0,
                                0
                            );
                        } else if (group_has_o1 && needs_sparse[src_idx] &&
                                   stat_uses_sparse_table(stype)) {
                            val = compute_minmax_with_sparse(
                                stype,
                                &sparse_tbl[src_idx],
                                source_data[src_idx],
                                g_start,
                                win_start,
                                win_end,
                                0,
                                0
                            );
                        } else {
                            val = compute_window_stat_simple(
                                stype,
                                source_data[src_idx],
                                win_start,
                                win_end,
                                0,
                                0
                            );
                        }

                        /* Apply result to ALL observations in this run */
                        for (size_t i = run_start; i < run_end; i++) {
                            result_data[s][i] = val;
                        }
                    }
                }

                free(run_starts);
            } else {
standard_obs_loop:;
                /* Rebuild O(1) structures for this group (reuse pre-allocated buffers) */
                int group_has_o1 = use_o1_queries && (g_count >= CRANGESTAT_PREFIX_SUM_THRESHOLD);
                if (group_has_o1) {
                    for (size_t v = 0; v < nsource; v++) {
                        if (needs_prefix[v]) {
                            prefix_arr[v].nobs = g_count;
                            build_prefix_arrays(source_data[v], g_start, g_count, &prefix_arr[v]);
                        }
                        if (needs_sparse[v]) {
                            rebuild_sparse_table(&sparse_tbl[v], source_data[v], g_start,
                                                 g_count, SV_missval);
                        }
                    }
                }

                /*
                 * OPTIMIZATION: Chunked parallel sliding window
                 *
                 * Each thread gets a contiguous chunk of non-missing observations.
                 * Each thread initializes its window with ONE binary search for its
                 * first observation, then slides the pointers forward (O(1) amortized).
                 * Total: O(P * log N + N) instead of O(N * log N).
                 */

                /* Find boundary before missing keys (they sort to end) */
                size_t g_end_nonmiss = g_end;
                while (g_end_nonmiss > g_start && key_data[g_end_nonmiss - 1] >= SV_missval)
                    g_end_nonmiss--;

                size_t nonmiss_count = g_end_nonmiss - g_start;

                #ifdef _OPENMP
                #pragma omp parallel if(nonmiss_count > CRANGESTAT_PARALLEL_THRESHOLD)
                #endif
                {
                    #ifdef _OPENMP
                    int tid = omp_get_thread_num();
                    int nthreads_actual = omp_get_num_threads();
                    #else
                    int tid = 0;
                    int nthreads_actual = 1;
                    #endif

                    double *work = work_arrays[tid];

                    /* Each thread gets a contiguous chunk */
                    size_t chunk = (nonmiss_count + (size_t)nthreads_actual - 1) / (size_t)nthreads_actual;
                    size_t my_start = g_start + (size_t)tid * chunk;
                    size_t my_end = my_start + chunk;
                    if (my_start > g_end_nonmiss) my_start = g_end_nonmiss;
                    if (my_end > g_end_nonmiss) my_end = g_end_nonmiss;

                    if (my_start < my_end) {
                        /* Initialize window with ONE binary search for first obs */
                        double first_key = key_data[my_start];
                        double first_low = config.low_is_missing ? -DBL_MAX : first_key + config.interval_low;
                        double first_high = config.high_is_missing ? DBL_MAX : first_key + config.interval_high;

                        size_t win_lo = lower_bound(key_data, g_start, g_end_nonmiss, first_low);
                        size_t win_hi = upper_bound(key_data, g_start, g_end_nonmiss, first_high);

                        for (size_t i = my_start; i < my_end; i++) {
                            double key_val = key_data[i];
                            double low_bound = config.low_is_missing ? -DBL_MAX : key_val + config.interval_low;
                            double high_bound = config.high_is_missing ? DBL_MAX : key_val + config.interval_high;

                            /* Slide window pointers forward */
                            while (win_lo < g_end_nonmiss && key_data[win_lo] < low_bound) win_lo++;
                            while (win_hi < g_end_nonmiss && key_data[win_hi] <= high_bound) win_hi++;

                            size_t win_start = win_lo;
                            size_t win_end = win_hi;

                            /* Compute each statistic for this observation */
                            for (size_t s = 0; s < nstats; s++) {
                                int src_idx = specs[s].source_var_idx;
                                if (src_idx < 0 || src_idx >= (int)nsource) continue;

                                double val;
                                stat_type stype = specs[s].type;

                                if (stat_needs_work_array(stype)) {
                                    val = compute_window_stat_complex(
                                        stype,
                                        source_data[src_idx],
                                        win_start,
                                        win_end,
                                        i,
                                        config.exclude_self,
                                        work
                                    );
                                } else if (group_has_o1 && needs_prefix[src_idx] &&
                                           stat_uses_prefix_sums(stype)) {
                                    val = compute_stat_with_prefix(
                                        stype,
                                        &prefix_arr[src_idx],
                                        source_data[src_idx],
                                        g_start,
                                        win_start,
                                        win_end,
                                        i,
                                        config.exclude_self
                                    );
                                } else if (group_has_o1 && needs_sparse[src_idx] &&
                                           stat_uses_sparse_table(stype)) {
                                    val = compute_minmax_with_sparse(
                                        stype,
                                        &sparse_tbl[src_idx],
                                        source_data[src_idx],
                                        g_start,
                                        win_start,
                                        win_end,
                                        i,
                                        config.exclude_self
                                    );
                                } else {
                                    val = compute_window_stat_simple(
                                        stype,
                                        source_data[src_idx],
                                        win_start,
                                        win_end,
                                        i,
                                        config.exclude_self
                                    );
                                }

                                result_data[s][i] = val;
                            }
                        }
                    }
                }
            }
        }

        /* Cleanup O(1) data structures (once, after all groups) */
        if (use_o1_queries) {
            for (size_t v = 0; v < nsource; v++) {
                if (needs_prefix[v]) {
                    free(prefix_arr[v].prefix_sum);
                    free(prefix_arr[v].prefix_sum2);
                    free(prefix_arr[v].prefix_count);
                }
                if (needs_sparse[v]) {
                    free_sparse_table(&sparse_tbl[v]);
                }
            }
            free(needs_prefix);
            free(needs_sparse);
            free(prefix_arr);
            free(sparse_tbl);
        }
    }

    t_compute = ctools_timer_seconds() - compute_start;

    /* === Store Phase === */
    double store_start = ctools_timer_seconds();

    /*
     * OPTIMIZATION: Parallel result storage
     * Store multiple result variables in parallel - each variable is independent.
     *
     * When fwd_perm is available, iterate in sorted order for sequential data[] reads.
     * fwd_perm[sorted_idx] = original_idx, so data[sorted_idx] is sequential.
     */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) if(nstats > 1)
    #endif
    for (size_t s = 0; s < nstats; s++) {
        int result_idx = specs[s].result_var_idx;
        double *data = result_data[s];

        if (fwd_perm) {
            /* Sequential reads from data[sorted_idx], scattered writes via obs_map */
            for (size_t sorted_idx = 0; sorted_idx < nobs; sorted_idx++) {
                perm_idx_t orig_idx = fwd_perm[sorted_idx];
                SF_vstore(result_idx, (ST_int)obs_map[orig_idx], data[sorted_idx]);
            }
        } else {
            for (size_t i = 0; i < nobs; i++) {
                SF_vstore(result_idx, (ST_int)obs_map[i], data[inv_perm[i]]);
            }
        }
    }

    t_store = ctools_timer_seconds() - store_start;

    /* === Cleanup === */
    for (int t = 0; t < num_threads; t++) {
        free(work_arrays[t]);
    }
    free(work_arrays);
    free(groups);
    free(inv_perm);
    if (fwd_perm) free(fwd_perm);

    /* Free result arrays (allocated separately from stata_data) */
    for (size_t s = 0; s < nstats; s++) free(result_data[s]);
    free(result_data);

    /* Free pointer arrays (not the underlying data, which is in stata_data) */
    if (by_data) free(by_data);
    free(source_data);

    /* Free the main data structure (includes key_data, source vars, by vars) */
    ctools_filtered_data_free(&filtered);
    free(load_indices);

    /* Store timing scalars */
    double total = ctools_timer_seconds() - t_start;
    SF_scal_save("_crangestat_time_load", t_load);
    SF_scal_save("_crangestat_time_sort", t_sort);
    SF_scal_save("_crangestat_time_groups", t_groups);
    SF_scal_save("_crangestat_time_compute", t_compute);
    SF_scal_save("_crangestat_time_store", t_store);
    SF_scal_save("_crangestat_time_total", total);
    SF_scal_save("_crangestat_nobs", (double)nobs);
    SF_scal_save("_crangestat_nstats", (double)nstats);
    SF_scal_save("_crangestat_ngroups", (double)ngroups);

    /* Free remaining allocations */
    free(source_indices);
    free(result_indices);
    free(specs);
    if (by_indices) free(by_indices);

    CTOOLS_SAVE_THREAD_INFO("_crangestat");

    /* Verbose output handled by .ado file using stored scalars */

    return 0;
}

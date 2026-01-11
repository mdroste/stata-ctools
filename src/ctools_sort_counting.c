/*
    ctools_sort_counting.c
    Optimized Parallel Counting Sort Module

    Counting sort is optimal for data with a small range of distinct values.
    This is common in Stata for categorical variables (e.g., foreign, rep78,
    state codes, year, month).

    Algorithm Overview:
    1. Find min/max values to determine range
    2. Parallel histogram: each thread counts values in its chunk
    3. Combine histograms into global counts
    4. Compute prefix sums (bucket offsets)
    5. Parallel scatter: each thread writes to pre-computed offsets

    Time Complexity: O(n + k) where k is the range of values
    Space Complexity: O(n + k)

    Optimizations in this version:
    - OpenMP for efficient parallelism (replaces pthreads)
    - Prefetching in hot loops
    - Loop unrolling for better ILP
    - Cache-line padded per-thread arrays to avoid false sharing
    - Pre-allocated scatter offset copies (no malloc in parallel regions)
    - Parallel min/max and histogram phases
*/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"

/* Maximum range for counting sort (larger ranges use radix) */
#define COUNTING_SORT_MAX_RANGE 1000000

/* Minimum observations to use parallel counting sort */
#define COUNTING_PARALLEL_THRESHOLD 50000

/* Maximum threads */
#define COUNTING_MAX_THREADS 16

/* Use PREFETCH_DISTANCE from ctools_config.h */
#define COUNTING_PREFETCH_DISTANCE PREFETCH_DISTANCE

/* Use centralized macros from ctools_config.h */
#define COUNTING_RESTRICT CTOOLS_RESTRICT
#define COUNTING_PREFETCH(addr) CTOOLS_PREFETCH(addr)
#define COUNTING_PREFETCH_W(addr) CTOOLS_PREFETCH_W(addr)

/* ============================================================================
   Utility Functions
   ============================================================================ */

static void *counting_aligned_alloc(size_t alignment, size_t size)
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

/* ============================================================================
   Sequential Counting Sort (for small datasets)
   ============================================================================ */

static stata_retcode counting_sort_numeric_seq(size_t * COUNTING_RESTRICT order,
                                               const double * COUNTING_RESTRICT data,
                                               size_t nobs)
{
    double min_val_d = data[order[0]];
    double max_val_d = data[order[0]];
    int64_t min_val;
    size_t range;
    size_t *counts = NULL;
    size_t *offsets = NULL;
    size_t *temp_order = NULL;
    size_t i;
    size_t missing_bucket;
    stata_retcode rc = STATA_OK;

    /* Find min/max */
    for (i = 0; i < nobs; i++) {
        double v = data[order[i]];
        if (SF_is_missing(v)) {
            continue;
        }
        if (v != floor(v)) {
            return STATA_ERR_UNSUPPORTED_TYPE;
        }
        if (v < min_val_d) min_val_d = v;
        if (v > max_val_d) max_val_d = v;
    }

    min_val = (int64_t)min_val_d;
    range = (size_t)((int64_t)max_val_d - min_val + 1);
    missing_bucket = range;

    if (range > COUNTING_SORT_MAX_RANGE) {
        return STATA_ERR_UNSUPPORTED_TYPE;
    }

    counts = (size_t *)calloc(range + 1, sizeof(size_t));
    offsets = (size_t *)malloc((range + 1) * sizeof(size_t));
    temp_order = (size_t *)malloc(nobs * sizeof(size_t));

    if (!counts || !offsets || !temp_order) {
        rc = STATA_ERR_MEMORY;
        goto cleanup;
    }

    /* Count with unrolling */
    i = 0;
    for (; i + 4 <= nobs; i += 4) {
        double v0 = data[order[i]];
        double v1 = data[order[i+1]];
        double v2 = data[order[i+2]];
        double v3 = data[order[i+3]];

        if (SF_is_missing(v0)) counts[missing_bucket]++;
        else counts[(size_t)((int64_t)v0 - min_val)]++;

        if (SF_is_missing(v1)) counts[missing_bucket]++;
        else counts[(size_t)((int64_t)v1 - min_val)]++;

        if (SF_is_missing(v2)) counts[missing_bucket]++;
        else counts[(size_t)((int64_t)v2 - min_val)]++;

        if (SF_is_missing(v3)) counts[missing_bucket]++;
        else counts[(size_t)((int64_t)v3 - min_val)]++;
    }
    for (; i < nobs; i++) {
        double v = data[order[i]];
        if (SF_is_missing(v)) {
            counts[missing_bucket]++;
        } else {
            counts[(size_t)((int64_t)v - min_val)]++;
        }
    }

    /* Prefix sum */
    offsets[0] = 0;
    for (i = 1; i <= missing_bucket; i++) {
        offsets[i] = offsets[i-1] + counts[i-1];
    }

    /* Scatter */
    for (i = 0; i < nobs; i++) {
        size_t idx = order[i];
        double v = data[idx];
        size_t bucket;

        if (SF_is_missing(v)) {
            bucket = missing_bucket;
        } else {
            bucket = (size_t)((int64_t)v - min_val);
        }

        temp_order[offsets[bucket]++] = idx;
    }

    memcpy(order, temp_order, nobs * sizeof(size_t));

cleanup:
    free(counts);
    free(offsets);
    free(temp_order);

    return rc;
}

/* ============================================================================
   Parallel Counting Sort Implementation
   ============================================================================ */

static stata_retcode counting_sort_numeric_parallel(size_t * COUNTING_RESTRICT order,
                                                    const double * COUNTING_RESTRICT data,
                                                    size_t nobs, int num_threads)
{
    stata_retcode rc = STATA_OK;

    if (num_threads > COUNTING_MAX_THREADS) num_threads = COUNTING_MAX_THREADS;

    /* Thread-local min/max results */
    double *thread_min = NULL;
    double *thread_max = NULL;
    int *thread_has_non_integer = NULL;

    /* Histogram data */
    size_t **local_counts = NULL;
    size_t *global_counts = NULL;
    size_t *global_offsets = NULL;
    size_t **thread_offsets = NULL;
    size_t *temp_order = NULL;

    double global_min, global_max;
    int64_t min_val;
    size_t range, missing_bucket;
    size_t chunk_size;
    int has_non_integer = 0;
    int t;

    chunk_size = (nobs + num_threads - 1) / num_threads;

    /* Allocate thread-local min/max arrays (cache-line aligned) */
    thread_min = (double *)counting_aligned_alloc(64, num_threads * sizeof(double));
    thread_max = (double *)counting_aligned_alloc(64, num_threads * sizeof(double));
    thread_has_non_integer = (int *)calloc(num_threads, sizeof(int));

    if (!thread_min || !thread_max || !thread_has_non_integer) {
        rc = STATA_ERR_MEMORY;
        goto cleanup;
    }

    /* Phase 1: Parallel min/max finding */
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        size_t start = (size_t)tid * chunk_size;
        size_t end = (size_t)(tid + 1) * chunk_size;
        if (end > nobs) end = nobs;

        double local_min = data[order[start < nobs ? start : 0]];
        double local_max = local_min;
        int local_non_int = 0;

        for (size_t i = start; i < end; i++) {
            double v = data[order[i]];

            if (SF_is_missing(v)) continue;

            if (v != floor(v)) {
                local_non_int = 1;
            }

            if (v < local_min) local_min = v;
            if (v > local_max) local_max = v;
        }

        thread_min[tid] = local_min;
        thread_max[tid] = local_max;
        thread_has_non_integer[tid] = local_non_int;
    }

    /* Combine thread results */
    global_min = thread_min[0];
    global_max = thread_max[0];

    for (t = 0; t < num_threads; t++) {
        if (thread_min[t] < global_min) global_min = thread_min[t];
        if (thread_max[t] > global_max) global_max = thread_max[t];
        if (thread_has_non_integer[t]) has_non_integer = 1;
    }

    if (has_non_integer) {
        rc = STATA_ERR_UNSUPPORTED_TYPE;
        goto cleanup;
    }

    min_val = (int64_t)global_min;
    range = (size_t)((int64_t)global_max - min_val + 1);
    missing_bucket = range;

    if (range > COUNTING_SORT_MAX_RANGE) {
        rc = STATA_ERR_UNSUPPORTED_TYPE;
        goto cleanup;
    }

    /* Allocate histogram arrays */
    /* Use cache-line padded range to avoid false sharing */
    size_t padded_range = ((range + 1 + 7) / 8) * 8;

    local_counts = (size_t **)malloc(num_threads * sizeof(size_t *));
    global_counts = (size_t *)calloc(range + 1, sizeof(size_t));
    global_offsets = (size_t *)malloc((range + 1) * sizeof(size_t));
    thread_offsets = (size_t **)malloc(num_threads * sizeof(size_t *));
    temp_order = (size_t *)counting_aligned_alloc(64, nobs * sizeof(size_t));

    if (!local_counts || !global_counts || !global_offsets ||
        !thread_offsets || !temp_order) {
        rc = STATA_ERR_MEMORY;
        goto cleanup;
    }

    /* Pre-allocate all per-thread count arrays */
    for (t = 0; t < num_threads; t++) {
        local_counts[t] = (size_t *)counting_aligned_alloc(64, padded_range * sizeof(size_t));
        thread_offsets[t] = (size_t *)malloc((range + 1) * sizeof(size_t));
        if (!local_counts[t] || !thread_offsets[t]) {
            rc = STATA_ERR_MEMORY;
            goto cleanup;
        }
        memset(local_counts[t], 0, padded_range * sizeof(size_t));
    }

    /* Phase 2: Parallel histogram */
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        size_t start = (size_t)tid * chunk_size;
        size_t end = (size_t)(tid + 1) * chunk_size;
        if (end > nobs) end = nobs;

        size_t *my_counts = local_counts[tid];

        /* Unrolled counting loop with prefetching */
        size_t i = start;
        for (; i + 4 <= end; i += 4) {
            if (i + COUNTING_PREFETCH_DISTANCE < end) {
                COUNTING_PREFETCH(&data[order[i + COUNTING_PREFETCH_DISTANCE]]);
            }

            double v0 = data[order[i]];
            double v1 = data[order[i+1]];
            double v2 = data[order[i+2]];
            double v3 = data[order[i+3]];

            if (SF_is_missing(v0)) my_counts[missing_bucket]++;
            else my_counts[(size_t)((int64_t)v0 - min_val)]++;

            if (SF_is_missing(v1)) my_counts[missing_bucket]++;
            else my_counts[(size_t)((int64_t)v1 - min_val)]++;

            if (SF_is_missing(v2)) my_counts[missing_bucket]++;
            else my_counts[(size_t)((int64_t)v2 - min_val)]++;

            if (SF_is_missing(v3)) my_counts[missing_bucket]++;
            else my_counts[(size_t)((int64_t)v3 - min_val)]++;
        }
        for (; i < end; i++) {
            double v = data[order[i]];
            if (SF_is_missing(v)) {
                my_counts[missing_bucket]++;
            } else {
                my_counts[(size_t)((int64_t)v - min_val)]++;
            }
        }
    }

    /* Phase 3: Combine histograms (can be parallelized for large ranges) */
    if (range > 10000) {
        /* Parallel reduction for large ranges */
        #pragma omp parallel for num_threads(num_threads)
        for (size_t b = 0; b <= missing_bucket; b++) {
            size_t sum = 0;
            for (int tt = 0; tt < num_threads; tt++) {
                sum += local_counts[tt][b];
            }
            global_counts[b] = sum;
        }
    } else {
        /* Sequential for small ranges */
        for (size_t b = 0; b <= missing_bucket; b++) {
            global_counts[b] = 0;
            for (t = 0; t < num_threads; t++) {
                global_counts[b] += local_counts[t][b];
            }
        }
    }

    /* Prefix sum */
    global_offsets[0] = 0;
    for (size_t b = 1; b <= missing_bucket; b++) {
        global_offsets[b] = global_offsets[b-1] + global_counts[b-1];
    }

    /* Compute per-thread offsets within each bucket */
    for (size_t b = 0; b <= missing_bucket; b++) {
        size_t offset = global_offsets[b];
        for (t = 0; t < num_threads; t++) {
            thread_offsets[t][b] = offset;
            offset += local_counts[t][b];
        }
    }

    /* Phase 4: Parallel scatter */
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        size_t start = (size_t)tid * chunk_size;
        size_t end = (size_t)(tid + 1) * chunk_size;
        if (end > nobs) end = nobs;

        /* Use pre-allocated thread_offsets - make local copy */
        size_t *my_offsets = thread_offsets[tid];

        /* Scatter with prefetching */
        for (size_t i = start; i < end; i++) {
            if (i + COUNTING_PREFETCH_DISTANCE < end) {
                COUNTING_PREFETCH(&data[order[i + COUNTING_PREFETCH_DISTANCE]]);
            }

            size_t idx = order[i];
            double v = data[idx];
            size_t bucket;

            if (SF_is_missing(v)) {
                bucket = missing_bucket;
            } else {
                bucket = (size_t)((int64_t)v - min_val);
            }

            temp_order[my_offsets[bucket]++] = idx;
        }
    }

    /* Copy result */
    memcpy(order, temp_order, nobs * sizeof(size_t));

cleanup:
    free(thread_min);
    free(thread_max);
    free(thread_has_non_integer);
    free(global_counts);
    free(global_offsets);
    free(temp_order);

    if (local_counts) {
        for (t = 0; t < num_threads; t++) {
            free(local_counts[t]);
        }
        free(local_counts);
    }

    if (thread_offsets) {
        for (t = 0; t < num_threads; t++) {
            free(thread_offsets[t]);
        }
        free(thread_offsets);
    }

    return rc;
}

/* ============================================================================
   Variable-Level Sort Function
   ============================================================================ */

static stata_retcode counting_sort_by_numeric_var(stata_data *data, int var_idx)
{
    const double *dbl_data = data->vars[var_idx].data.dbl;
    int num_threads;

    /* Determine thread count */
    num_threads = omp_get_max_threads();
    if (num_threads > COUNTING_MAX_THREADS) num_threads = COUNTING_MAX_THREADS;

    if (data->nobs < COUNTING_PARALLEL_THRESHOLD) {
        return counting_sort_numeric_seq(data->sort_order, dbl_data, data->nobs);
    }

    if (data->nobs < (size_t)MIN_OBS_PER_THREAD * (size_t)num_threads) {
        num_threads = (int)(data->nobs / MIN_OBS_PER_THREAD);
        if (num_threads < 2) num_threads = 2;
    }

    return counting_sort_numeric_parallel(data->sort_order, dbl_data,
                                          data->nobs, num_threads);
}

/* ============================================================================
   Permutation Application (Parallel)
   ============================================================================ */

static stata_retcode counting_apply_permutation(stata_data *data)
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
            double *new_data = (double *)counting_aligned_alloc(64, nobs * sizeof(double));
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

/*
    Check if counting sort is appropriate for this variable.
    Returns 1 if suitable, 0 otherwise.
*/
int ctools_counting_sort_suitable(stata_data *data, int var_idx)
{
    double *dbl_data;
    double min_val, max_val;
    size_t i;
    size_t range;

    if (data->vars[var_idx].type != STATA_TYPE_DOUBLE) {
        return 0;  /* Strings not supported */
    }

    dbl_data = data->vars[var_idx].data.dbl;

    /* Sample to check suitability (check every 1000th element) */
    min_val = max_val = dbl_data[0];
    for (i = 0; i < data->nobs; i += 1000) {
        double v = dbl_data[i];
        if (SF_is_missing(v)) continue;
        if (v != floor(v)) return 0;  /* Non-integer */
        if (v < min_val) min_val = v;
        if (v > max_val) max_val = v;
    }

    range = (size_t)((int64_t)max_val - (int64_t)min_val + 1);

    return (range <= COUNTING_SORT_MAX_RANGE);
}

/*
    Sort data using parallel counting sort.

    This algorithm is optimal for integer data with a small range of values.
    It achieves O(n + k) time complexity where k is the range.

    Returns STATA_ERR_UNSUPPORTED_TYPE if data is not suitable for counting sort
    (non-integer values or range too large). Caller should fall back to radix sort.
*/
stata_retcode ctools_sort_counting(stata_data *data, int *sort_vars, size_t nsort)
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

        if (data->vars[var_idx].type != STATA_TYPE_DOUBLE) {
            return STATA_ERR_UNSUPPORTED_TYPE;
        }

        rc = counting_sort_by_numeric_var(data, var_idx);
        if (rc != STATA_OK) {
            return rc;
        }
    }

    rc = counting_apply_permutation(data);
    if (rc != STATA_OK) {
        return rc;
    }

    return STATA_OK;
}

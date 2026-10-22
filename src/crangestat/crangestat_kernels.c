#include "crangestat_internal.h"

size_t crangestat_lower_bound(const double *key, size_t start, size_t end, double target)
{
    while (start < end)
    {
        size_t mid = start + (end - start) / 2;
        if (key[mid] < target)
        {
            start = mid + 1;
        }
        else
        {
            end = mid;
        }
    }
    return start;
}

size_t crangestat_upper_bound(const double *key, size_t start, size_t end, double target)
{
    while (start < end)
    {
        size_t mid = start + (end - start) / 2;
        if (key[mid] <= target)
        {
            start = mid + 1;
        }
        else
        {
            end = mid;
        }
    }
    return start;
}

double crangestat_compute_percentile(double *work, size_t n, double pctl)
{
    if (n == 0)
        return SV_missval;
    if (n == 1)
        return work[0];

    double pos = (pctl / 100.0) * (double)(n - 1);
    size_t k = (size_t)floor(pos);
    if (k >= n)
        k = n - 1;

    double val = ctools_quickselect_double(work, n, k);

    double frac = pos - (double)k;
    if (frac > 0.0 && k + 1 < n)
    {
        /* work[k+1..n-1] >= val after the select; its min is sorted[k+1] */
        double next_val = DBL_MAX;
        for (size_t i = k + 1; i < n; i++)
        {
            if (work[i] < next_val)
                next_val = work[i];
        }
        if (next_val < DBL_MAX)
        {
            val = val + frac * (next_val - val);
        }
    }
    return val;
}

void crangestat_build_prefix_arrays(const double *data, size_t start, size_t count,
                                    prefix_arrays *pa)
{
    const double miss = SV_missval;
    pa->origin = 0.0;
    for (size_t i = 0; i < count; i++)
    {
        if (data[start + i] < miss)
        {
            pa->origin = data[start + i];
            break;
        }
    }

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
     * The chunks are a fixed number of logical work items shared out with
     * "omp for", never one chunk per thread id: this function also runs inside
     * the cross-group parallel loop (where a nested region has one thread), and
     * a runtime may supply a smaller team than requested. Fixed chunking also
     * makes the sums independent of the thread count.
     *
     * For small groups, sequential is faster due to less overhead.
     */
    enum
    {
        PREFIX_CHUNKS = 64
    };
    size_t parallel_threshold = 10000;

    if (count >= parallel_threshold)
    {
        size_t chunk_size = (count + PREFIX_CHUNKS - 1) / PREFIX_CHUNKS;

        /* Per-chunk totals */
        double chunk_sum[PREFIX_CHUNKS];
        double chunk_sum2[PREFIX_CHUNKS];
        size_t chunk_count[PREFIX_CHUNKS];

/* Phase 1: Compute local prefix sums within each chunk (parallel) */
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (!omp_in_parallel())
#endif
        for (int c = 0; c < PREFIX_CHUNKS; c++)
        {
            size_t my_start = (size_t)c * chunk_size;
            size_t my_end = my_start + chunk_size;
            if (my_start > count)
                my_start = count;
            if (my_end > count)
                my_end = count;

            double local_sum = 0.0, local_sum2 = 0.0;
            size_t local_count = 0;

            for (size_t i = my_start; i < my_end; i++)
            {
                double val = data[start + i];
                size_t idx = i + 1;

                if (val < miss)
                {
                    val -= pa->origin;
                    local_sum += val;
                    local_sum2 += val * val;
                    local_count++;
                }

                pa->prefix_sum[idx] = local_sum;
                pa->prefix_sum2[idx] = local_sum2;
                pa->prefix_count[idx] = local_count;
            }

            chunk_sum[c] = local_sum;
            chunk_sum2[c] = local_sum2;
            chunk_count[c] = local_count;
        }

        /* Phase 2: Compute prefix of chunk totals (sequential) */
        double prefix_chunk_sum = 0.0, prefix_chunk_sum2 = 0.0;
        size_t prefix_chunk_count = 0;

        for (int c = 0; c < PREFIX_CHUNKS; c++)
        {
            double old_sum = chunk_sum[c];
            double old_sum2 = chunk_sum2[c];
            size_t old_count = chunk_count[c];

            chunk_sum[c] = prefix_chunk_sum;
            chunk_sum2[c] = prefix_chunk_sum2;
            chunk_count[c] = prefix_chunk_count;

            prefix_chunk_sum += old_sum;
            prefix_chunk_sum2 += old_sum2;
            prefix_chunk_count += old_count;
        }

/* Phase 3: Add chunk prefix to each element (parallel);
 * chunk 0 is already correct */
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (!omp_in_parallel())
#endif
        for (int c = 1; c < PREFIX_CHUNKS; c++)
        {
            size_t my_start = (size_t)c * chunk_size;
            size_t my_end = my_start + chunk_size;
            if (my_start > count)
                my_start = count;
            if (my_end > count)
                my_end = count;

            double add_sum = chunk_sum[c];
            double add_sum2 = chunk_sum2[c];
            size_t add_count = chunk_count[c];

            for (size_t i = my_start; i < my_end; i++)
            {
                size_t idx = i + 1;
                pa->prefix_sum[idx] += add_sum;
                pa->prefix_sum2[idx] += add_sum2;
                pa->prefix_count[idx] += add_count;
            }
        }
    }
    else
    {
        /* Sequential prefix sum for small groups */
        double sum = 0.0, sum2 = 0.0;
        size_t cnt = 0;

        for (size_t i = 0; i < count; i++)
        {
            double val = data[start + i];
            size_t idx = i + 1;

            if (val < miss)
            {
                val -= pa->origin;
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

void crangestat_query_prefix_sums(const prefix_arrays *pa, size_t g_start, size_t win_start,
                                  size_t win_end, double *out_sum, double *out_sum2,
                                  size_t *out_count)
{
    size_t a = win_start - g_start; /* Convert to group-relative indices */
    size_t b = win_end - g_start;

    *out_sum = pa->prefix_sum[b] - pa->prefix_sum[a];
    *out_sum2 = pa->prefix_sum2[b] - pa->prefix_sum2[a];
    *out_count = pa->prefix_count[b] - pa->prefix_count[a];
}

int crangestat_build_sparse_table(const double *data, size_t start, size_t count, sparse_table *st,
                                  const double miss)
{
    memset(st, 0, sizeof(*st));
    if (count == 0)
        return 0;

    /* Compute max_log = floor(log2(count)) */
    int max_log = 0;
    size_t temp = count;
    while (temp > 1)
    {
        temp >>= 1;
        max_log++;
    }
    if (max_log >= CRANGESTAT_SPARSE_TABLE_LOG)
    {
        max_log = CRANGESTAT_SPARSE_TABLE_LOG - 1;
    }
    st->max_log = max_log;
    st->nobs = count;

    /* Allocate log table */
    st->log_table = (int *)malloc((count + 1) * sizeof(int));
    if (!st->log_table)
        return -1;

    st->log_table[0] = 0;
    st->log_table[1] = 0;
    for (size_t i = 2; i <= count; i++)
    {
        st->log_table[i] = st->log_table[i / 2] + 1;
    }

    /* Allocate pointer arrays for sparse tables */
    int nlevels = max_log + 1;
    st->sparse_min = (double **)malloc(nlevels * sizeof(double *));
    st->sparse_max = (double **)malloc(nlevels * sizeof(double *));
    if (!st->sparse_min || !st->sparse_max)
    {
        free(st->log_table);
        if (st->sparse_min)
            free(st->sparse_min);
        if (st->sparse_max)
            free(st->sparse_max);
        memset(st, 0, sizeof(*st));
        return -1;
    }

    /* Single contiguous allocation for all sparse table data:
     * 2 * nlevels * count doubles (min and max for each level) */
    double *block = (double *)malloc(2 * (size_t)nlevels * count * sizeof(double));
    if (!block)
    {
        free(st->sparse_min);
        free(st->sparse_max);
        free(st->log_table);
        memset(st, 0, sizeof(*st));
        return -1;
    }

    for (int j = 0; j < nlevels; j++)
    {
        st->sparse_min[j] = block + (size_t)j * count;
        st->sparse_max[j] = block + (size_t)(nlevels + j) * count;
    }

    /* Initialize level 0: individual elements */
    for (size_t i = 0; i < count; i++)
    {
        double val = data[start + i];
        if (val >= miss)
        {
            /* Missing values: use sentinel values */
            st->sparse_min[0][i] = DBL_MAX;
            st->sparse_max[0][i] = -DBL_MAX;
        }
        else
        {
            st->sparse_min[0][i] = val;
            st->sparse_max[0][i] = val;
        }
    }

    /* Build higher levels */
    for (int j = 1; j <= max_log; j++)
    {
        size_t range_size = (size_t)1 << j;
        size_t half_range = (size_t)1 << (j - 1);
        size_t limit = (range_size <= count) ? count - range_size + 1 : 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (count > 10000)
#endif
        for (size_t i = 0; i < limit; i++)
        {
            double min1 = st->sparse_min[j - 1][i];
            double min2 = st->sparse_min[j - 1][i + half_range];
            double max1 = st->sparse_max[j - 1][i];
            double max2 = st->sparse_max[j - 1][i + half_range];

            st->sparse_min[j][i] = (min1 < min2) ? min1 : min2;
            st->sparse_max[j][i] = (max1 > max2) ? max1 : max2;
        }
    }

    return 0;
}

double crangestat_query_range_min(const sparse_table *st, size_t l, size_t r)
{
    if (l >= r)
        return DBL_MAX;
    int k = st->log_table[r - l];
    if (k > st->max_log)
        k = st->max_log;
    size_t range = (size_t)1 << k;
    const double *level = st->sparse_min[k];

    double min1 = level[l];
    for (size_t p = l + range; p + range < r; p += range)
    {
        if (level[p] < min1)
            min1 = level[p];
    }
    double min2 = level[r - range];

    return (min1 < min2) ? min1 : min2;
}

double crangestat_query_range_max(const sparse_table *st, size_t l, size_t r)
{
    if (l >= r)
        return -DBL_MAX;
    int k = st->log_table[r - l];
    if (k > st->max_log)
        k = st->max_log;
    size_t range = (size_t)1 << k;
    const double *level = st->sparse_max[k];

    double max1 = level[l];
    for (size_t p = l + range; p + range < r; p += range)
    {
        if (level[p] > max1)
            max1 = level[p];
    }
    double max2 = level[r - range];

    return (max1 > max2) ? max1 : max2;
}

void crangestat_free_sparse_table(sparse_table *st)
{
    if (st->sparse_min)
    {
        /* All level data is in a single contiguous block starting at sparse_min[0] */
        free(st->sparse_min[0]);
        free(st->sparse_min);
    }
    if (st->sparse_max)
    {
        /* sparse_max pointers point into the same block; just free the pointer array */
        free(st->sparse_max);
    }
    if (st->log_table)
        free(st->log_table);
    memset(st, 0, sizeof(*st));
}

int crangestat_alloc_sparse_table(sparse_table *st, size_t max_count)
{
    if (max_count == 0)
        return 0;

    /* Compute max_log for the largest possible group */
    int max_log = 0;
    size_t temp = max_count;
    while (temp > 1)
    {
        temp >>= 1;
        max_log++;
    }
    if (max_log >= CRANGESTAT_SPARSE_TABLE_LOG)
        max_log = CRANGESTAT_SPARSE_TABLE_LOG - 1;

    st->max_log = max_log;
    st->nobs = max_count;

    /* Allocate log table for [0, max_count] */
    st->log_table = (int *)malloc((max_count + 1) * sizeof(int));
    if (!st->log_table)
        return -1;

    /* Pre-fill log table */
    st->log_table[0] = 0;
    if (max_count >= 1)
        st->log_table[1] = 0;
    for (size_t i = 2; i <= max_count; i++)
    {
        st->log_table[i] = st->log_table[i / 2] + 1;
    }

    /* Allocate pointer arrays */
    int nlevels = max_log + 1;
    st->sparse_min = (double **)malloc(nlevels * sizeof(double *));
    st->sparse_max = (double **)malloc(nlevels * sizeof(double *));
    if (!st->sparse_min || !st->sparse_max)
    {
        free(st->log_table);
        if (st->sparse_min)
            free(st->sparse_min);
        if (st->sparse_max)
            free(st->sparse_max);
        st->log_table = NULL;
        st->sparse_min = NULL;
        st->sparse_max = NULL;
        return -1;
    }

    /* Single contiguous allocation: stride = max_count per level */
    double *block = (double *)malloc(2 * (size_t)nlevels * max_count * sizeof(double));
    if (!block)
    {
        free(st->sparse_min);
        free(st->sparse_max);
        free(st->log_table);
        st->sparse_min = NULL;
        st->sparse_max = NULL;
        st->log_table = NULL;
        return -1;
    }

    for (int j = 0; j < nlevels; j++)
    {
        st->sparse_min[j] = block + (size_t)j * max_count;
        st->sparse_max[j] = block + (size_t)(nlevels + j) * max_count;
    }

    return 0;
}

void crangestat_rebuild_sparse_table(sparse_table *st, const double *data, size_t start,
                                     size_t count, double miss)
{
    if (count == 0)
    {
        st->nobs = 0;
        return;
    }

    /* Recompute max_log for this group's count */
    int max_log = 0;
    size_t temp = count;
    while (temp > 1)
    {
        temp >>= 1;
        max_log++;
    }
    if (max_log >= CRANGESTAT_SPARSE_TABLE_LOG)
        max_log = CRANGESTAT_SPARSE_TABLE_LOG - 1;

    st->max_log = max_log;
    st->nobs = count;

    /* Rebuild log table for [0, count] */
    st->log_table[0] = 0;
    if (count >= 1)
        st->log_table[1] = 0;
    for (size_t i = 2; i <= count; i++)
    {
        st->log_table[i] = st->log_table[i / 2] + 1;
    }

    /* Initialize level 0 */
    for (size_t i = 0; i < count; i++)
    {
        double val = data[start + i];
        if (val >= miss)
        {
            st->sparse_min[0][i] = DBL_MAX;
            st->sparse_max[0][i] = -DBL_MAX;
        }
        else
        {
            st->sparse_min[0][i] = val;
            st->sparse_max[0][i] = val;
        }
    }

    /* Build higher levels */
    for (int j = 1; j <= max_log; j++)
    {
        size_t range_size = (size_t)1 << j;
        size_t half_range = (size_t)1 << (j - 1);
        size_t limit = (range_size <= count) ? count - range_size + 1 : 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (count > 10000)
#endif
        for (size_t i = 0; i < limit; i++)
        {
            double min1 = st->sparse_min[j - 1][i];
            double min2 = st->sparse_min[j - 1][i + half_range];
            double max1 = st->sparse_max[j - 1][i];
            double max2 = st->sparse_max[j - 1][i + half_range];

            st->sparse_min[j][i] = (min1 < min2) ? min1 : min2;
            st->sparse_max[j][i] = (max1 > max2) ? max1 : max2;
        }
    }
}

int crangestat_stat_needs_work_array(stat_type type)
{
    switch (type)
    {
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

int crangestat_stat_uses_prefix_sums(stat_type type)
{
    switch (type)
    {
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

int crangestat_stat_uses_sparse_table(stat_type type)
{
    switch (type)
    {
    case STAT_MIN:
    case STAT_MAX:
        return 1;
    default:
        return 0;
    }
}

double crangestat_stable_window_variance(const double *data, size_t start, size_t end, size_t self,
                                         int exclude_self)
{
    size_t n = 0;
    long double origin = 0, sum = 0;
    for (size_t i = start; i < end; i++)
    {
        if ((exclude_self && i == self) || !(data[i] < SV_missval))
            continue;
        if (n == 0)
            origin = data[i];
        sum += (long double)data[i] - origin;
        n++;
    }
    if (n < 2)
        return SV_missval;
    long double mean = sum / n, ss = 0, correction = 0;
    for (size_t i = start; i < end; i++)
    {
        if ((exclude_self && i == self) || !(data[i] < SV_missval))
            continue;
        long double d = ((long double)data[i] - origin) - mean;
        ss += d * d;
        correction += d;
    }
    return (double)((ss - correction * correction / n) / (n - 1));
}

double crangestat_compute_stat_with_prefix(stat_type type, const prefix_arrays *pa,
                                           const double *data, size_t g_start, size_t win_start,
                                           size_t win_end, size_t self_idx, int exclude_self)
{
    double sum, sum2;
    size_t n;

    /* A window without rows (including an inverted interval, low > high)
     * gives missing for every statistic, as in rangestat */
    if (win_end <= win_start)
        return SV_missval;

    /* Query prefix sums - O(1) */
    crangestat_query_prefix_sums(pa, g_start, win_start, win_end, &sum, &sum2, &n);

    /* Handle excludeself by subtracting self contribution */
    if (exclude_self && self_idx >= win_start && self_idx < win_end)
    {
        double self_val = data[self_idx];
        if (self_val < SV_missval)
        {
            self_val -= pa->origin;
            sum -= self_val;
            sum2 -= self_val * self_val;
            n -= 1;
        }
    }

    if (n == 0)
    {
        /* rangestat: a window with rows but no nonmissing value has count 0
         * and sum 0; only a window without rows leaves them missing. */
        if (win_end > win_start && (type == STAT_COUNT || type == STAT_SUM))
            return 0.0;
        return SV_missval;
    }

    if (type == STAT_SUM || type == STAT_MEAN) {
        double total = sum + (double)n * pa->origin;
        double scale = fabs(pa->prefix_sum[win_end-g_start]) +
                       fabs(pa->prefix_sum[win_start-g_start]) + fabs((double)n*pa->origin);
        if (!isfinite(total) || fabs(total) <= sqrt(DBL_EPSILON)*scale) {
            /* Local Neumaier sum avoids both prefix subtraction and the
             * cancellation from adding back a distant group origin. */
            double local=0.0, correction=0.0;
            for (size_t i=win_start;i<win_end;i++) {
                if ((exclude_self && i==self_idx) || data[i]>=SV_missval) continue;
                double x=data[i], next=local+x;
                correction += fabs(local)>=fabs(x) ? (local-next)+x : (x-next)+local;
                local=next;
            }
            total=local+correction;
        }
        return type==STAT_MEAN ? total/(double)n : total;
    }
    switch (type)
    {
    case STAT_COUNT:
        return (double)n;

    case STAT_SUM:
        return sum + (double)n * pa->origin;

    case STAT_MEAN:
        return pa->origin + sum / (double)n;

    case STAT_SD:
    case STAT_VARIANCE:
    {
        if (n < 2)
            return SV_missval;
        double mean = sum / (double)n;
        double ss = sum2 - sum * mean;
        double scale = fabs(pa->prefix_sum2[win_end - g_start]) +
                       fabs(pa->prefix_sum2[win_start - g_start]) + fabs(sum * mean);
        double var = ss / (double)(n - 1);
        /* Fall back before cancellation consumes half the mantissa. */
        if (!isfinite(var) || ss <= sqrt(DBL_EPSILON) * scale)
            var =
                crangestat_stable_window_variance(data, win_start, win_end, self_idx, exclude_self);
        return (type == STAT_SD) ? sqrt(var) : var;
    }

    default:
        return SV_missval;
    }
}

double crangestat_compute_minmax_with_sparse(stat_type type, const sparse_table *st,
                                             const double *data, size_t g_start, size_t win_start,
                                             size_t win_end, size_t self_idx, int exclude_self)
{
    if (win_start >= win_end)
        return SV_missval;

    /* Convert to group-relative indices */
    size_t l = win_start - g_start;
    size_t r = win_end - g_start;

    double result;
    if (type == STAT_MIN)
    {
        result = crangestat_query_range_min(st, l, r);
        if (result >= DBL_MAX - 1.0)
            return SV_missval; /* All missing */
    }
    else
    {
        result = crangestat_query_range_max(st, l, r);
        if (result <= -DBL_MAX + 1.0)
            return SV_missval; /* All missing */
    }

    /* Handle excludeself */
    if (exclude_self && self_idx >= win_start && self_idx < win_end)
    {
        double self_val = data[self_idx];
        if (self_val < SV_missval)
        {
            if ((type == STAT_MIN && self_val == result) ||
                (type == STAT_MAX && self_val == result))
            {
                /*
                    Self equals the min/max - need to find second-best.
                    Query left and right of self separately.
                */
                size_t self_rel = self_idx - g_start;
                double left_result, right_result;

                if (type == STAT_MIN)
                {
                    left_result =
                        (self_rel > l) ? crangestat_query_range_min(st, l, self_rel) : DBL_MAX;
                    right_result = (self_rel + 1 < r)
                                       ? crangestat_query_range_min(st, self_rel + 1, r)
                                       : DBL_MAX;
                    result = (left_result < right_result) ? left_result : right_result;
                    if (result >= DBL_MAX - 1.0)
                        return SV_missval;
                }
                else
                {
                    left_result =
                        (self_rel > l) ? crangestat_query_range_max(st, l, self_rel) : -DBL_MAX;
                    right_result = (self_rel + 1 < r)
                                       ? crangestat_query_range_max(st, self_rel + 1, r)
                                       : -DBL_MAX;
                    result = (left_result > right_result) ? left_result : right_result;
                    if (result <= -DBL_MAX + 1.0)
                        return SV_missval;
                }
            }
        }
    }

    return result;
}

void crangestat_simd_sum_minmax(const double *data, size_t start, size_t end, double *out_sum,
                                double *out_sum2, double *out_min, double *out_max, size_t *out_n)
{
    const double miss = SV_missval;
    double sum = 0.0, sum2 = 0.0;
    double minval = DBL_MAX, maxval = -DBL_MAX;
    size_t n = 0;
    size_t i = start;

#if CTOOLS_HAS_AVX2
    /* AVX2: process 4 doubles at a time */
    if (end - start >= 8)
    {
        __m256d vsum = _mm256_setzero_pd();
        __m256d vsum2 = _mm256_setzero_pd();
        __m256d vmin = _mm256_set1_pd(DBL_MAX);
        __m256d vmax = _mm256_set1_pd(-DBL_MAX);
        __m256d vmiss = _mm256_set1_pd(miss);

        for (; i + 4 <= end; i += 4)
        {
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
    if (end - start >= 4)
    {
        float64x2_t vsum = vdupq_n_f64(0.0);
        float64x2_t vsum2 = vdupq_n_f64(0.0);
        float64x2_t vmin = vdupq_n_f64(DBL_MAX);
        float64x2_t vmax = vdupq_n_f64(-DBL_MAX);
        float64x2_t vmiss = vdupq_n_f64(miss);

        for (; i + 2 <= end; i += 2)
        {
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
    for (; i < end; i++)
    {
        double val = data[i];
        if (val < miss)
        {
            sum += val;
            sum2 += val * val;
            if (val < minval)
                minval = val;
            if (val > maxval)
                maxval = val;
            n++;
        }
    }

    *out_sum = sum;
    *out_sum2 = sum2;
    *out_min = minval;
    *out_max = maxval;
    *out_n = n;
}

double crangestat_compute_window_stat_simple(stat_type type, const double *data, size_t win_start,
                                             size_t win_end, size_t self_idx, int exclude_self)
{
    /* Window without rows (or inverted interval): missing, as in rangestat */
    if (win_end <= win_start)
        return SV_missval;
    /* rangestat's first/last retain endpoint missing values. excludeself
     * replaces the self row with ordinary missing; it does not remove it. */
    if (type == STAT_FIRST || type == STAT_LAST)
    {
        size_t endpoint = type == STAT_FIRST ? win_start : win_end - 1;
        return exclude_self && endpoint == self_idx ? SV_missval : data[endpoint];
    }
    if (type == STAT_SD || type == STAT_VARIANCE)
    {
        double var =
            crangestat_stable_window_variance(data, win_start, win_end, self_idx, exclude_self);
        return type == STAT_SD && var < SV_missval ? sqrt(var) : var;
    }
    const double miss = SV_missval;
    size_t n = 0;
    double sum = 0.0, sum2 = 0.0;
    double minval = DBL_MAX, maxval = -DBL_MAX;
    double first_val = miss, last_val = miss;

    /* Check if self is in window and non-missing (for excludeself optimization) */
    int self_in_window = 0;
    int self_is_missing = 1;
    double self_val = 0.0;
    if (exclude_self && self_idx >= win_start && self_idx < win_end)
    {
        self_in_window = 1;
        self_val = data[self_idx];
        self_is_missing = (self_val >= miss);
    }

    /* Use SIMD for large windows */
    if ((win_end - win_start) >= 8)
    {
        crangestat_simd_sum_minmax(data, win_start, win_end, &sum, &sum2, &minval, &maxval, &n);

        /* SIMD-then-subtract optimization for excludeself */
        if (exclude_self && self_in_window && !self_is_missing)
        {
            /* Subtract self's contribution from aggregates */
            n -= 1;
            sum -= self_val;
            sum2 -= self_val * self_val;

            /* For min/max: if self equals min/max, need to rescan */
            /* Only rescan if we actually need min or max */
            if (type == STAT_MIN && self_val == minval)
            {
                minval = DBL_MAX;
                for (size_t i = win_start; i < win_end; i++)
                {
                    if (i == self_idx)
                        continue;
                    double val = data[i];
                    if (val < miss && val < minval)
                        minval = val;
                }
            }
            if (type == STAT_MAX && self_val == maxval)
            {
                maxval = -DBL_MAX;
                for (size_t i = win_start; i < win_end; i++)
                {
                    if (i == self_idx)
                        continue;
                    double val = data[i];
                    if (val < miss && val > maxval)
                        maxval = val;
                }
            }
        }

        /* Find first and last non-missing (still need sequential scan) */
        if (type == STAT_FIRST || type == STAT_FIRSTNM || type == STAT_LAST || type == STAT_LASTNM)
        {
            for (size_t i = win_start; i < win_end; i++)
            {
                if (exclude_self && i == self_idx)
                    continue;
                if (data[i] < miss)
                {
                    first_val = data[i];
                    break;
                }
            }
            for (size_t i = win_end; i > win_start; i--)
            {
                if (exclude_self && (i - 1) == self_idx)
                    continue;
                if (data[i - 1] < miss)
                {
                    last_val = data[i - 1];
                    break;
                }
            }
        }
    }
    else
    {
        /* Scalar path for small windows */
        for (size_t i = win_start; i < win_end; i++)
        {
            if (exclude_self && i == self_idx)
                continue;
            double val = data[i];
            if (val < miss)
            {
                if (n == 0)
                    first_val = val;
                last_val = val;
                sum += val;
                sum2 += val * val;
                if (val < minval)
                    minval = val;
                if (val > maxval)
                    maxval = val;
                n++;
            }
        }
    }

    if (n == 0)
    {
        /* rangestat: count and sum are 0 for a window whose rows are all
         * missing, and missing only when the window has no rows. */
        if (win_end > win_start && (type == STAT_COUNT || type == STAT_SUM))
            return 0.0;
        return SV_missval;
    }

    switch (type)
    {
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
    case STAT_VARIANCE:
    {
        if (n < 2)
            return SV_missval;
        double mean = sum / (double)n;
        double var = (sum2 - sum * mean) / (double)(n - 1);
        if (var < 0.0)
            var = 0.0;
        return (type == STAT_SD) ? sqrt(var) : var;
    }

    case STAT_FIRST:
    case STAT_FIRSTNM:
        return first_val;

    case STAT_LAST:
    case STAT_LASTNM:
        return last_val;

    case STAT_SKEWNESS:
    {
        /* rangestat uses standardized population moments, which are defined
         * for two distinct values (no small-sample bias correction). */
        if (n < 2)
            return SV_missval;
        double mean = sum / (double)n;
        /* Need second pass for higher moments */
        double m2 = 0.0, m3 = 0.0;
        for (size_t i = win_start; i < win_end; i++)
        {
            if (exclude_self && i == self_idx)
                continue;
            double val = data[i];
            if (val < miss)
            {
                double d = val - mean;
                m2 += d * d;
                m3 += d * d * d;
            }
        }
        m2 /= (double)n;
        m3 /= (double)n;
        if (m2 <= 0.0)
            return SV_missval;
        double sd = sqrt(m2);
        return m3 / (sd * sd * sd);
    }

    case STAT_KURTOSIS:
    {
        if (n < 2)
            return SV_missval;
        double mean = sum / (double)n;
        /* Need second pass for higher moments */
        double m2 = 0.0, m4 = 0.0;
        for (size_t i = win_start; i < win_end; i++)
        {
            if (exclude_self && i == self_idx)
                continue;
            double val = data[i];
            if (val < miss)
            {
                double d = val - mean;
                double d2 = d * d;
                m2 += d2;
                m4 += d2 * d2;
            }
        }
        m2 /= (double)n;
        m4 /= (double)n;
        if (m2 <= 0.0)
            return SV_missval;
        return m4 / (m2 * m2);
    }

    default:
        return SV_missval;
    }
}

double crangestat_compute_window_stat_complex(stat_type type, const double *data, size_t win_start,
                                              size_t win_end, size_t self_idx, int exclude_self,
                                              double *work)
{
    const double miss = SV_missval;

    /* Extract non-missing values into work array */
    size_t n = 0;
    for (size_t i = win_start; i < win_end; i++)
    {
        if (exclude_self && i == self_idx)
            continue;
        double val = data[i];
        if (val < miss)
        {
            work[n++] = val;
        }
    }

    if (n == 0)
        return SV_missval;

    switch (type)
    {
    case STAT_MEDIAN:
        return crangestat_compute_percentile(work, n, 50.0);

    case STAT_IQR:
    {
        double p25 = crangestat_compute_percentile(work, n, 25.0);
        double p75 = crangestat_compute_percentile(work, n, 75.0);
        if (p25 >= miss || p75 >= miss)
            return SV_missval;
        return p75 - p25;
    }

    case STAT_P1:
        return crangestat_compute_percentile(work, n, 1.0);
    case STAT_P5:
        return crangestat_compute_percentile(work, n, 5.0);
    case STAT_P10:
        return crangestat_compute_percentile(work, n, 10.0);
    case STAT_P25:
        return crangestat_compute_percentile(work, n, 25.0);
    case STAT_P75:
        return crangestat_compute_percentile(work, n, 75.0);
    case STAT_P90:
        return crangestat_compute_percentile(work, n, 90.0);
    case STAT_P95:
        return crangestat_compute_percentile(work, n, 95.0);
    case STAT_P99:
        return crangestat_compute_percentile(work, n, 99.0);

    default:
        return SV_missval;
    }
}

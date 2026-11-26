/*
    cegen_stats.c
    Statistics of one source within by-groups

    count, mean, total, sd, min, max, mdev, skew, kurt, std, pc, median,
    pctile, iqr, mad and seq.

    egen computes these two ways, and so does this file:

    - Without by, egen takes the moments from summarize. Here the sample is
      one group: sums run over fixed blocks of rows in parallel, each with
      compensated (Neumaier) summation, and the block sums are combined in
      block order, so the result does not depend on the number of threads.

    - With by, egen accumulates running sums over each group in observation
      order. Here each group is summed in plain double precision in the same
      order, so data sorted by the by-variables give the bits egen gives.
      Sorted data are runs of rows, evaluated in parallel. Unsorted data are
      not sorted for the moment functions: the groups are numbered, and one
      pass over the rows per moment adds each value to its group's sum, in
      observation order again. Only the order statistics gather each
      group's values.

    The percentile arithmetic differs the same way (_pctile without by, a
    formula on the sorted values with by) and both rules are kept.

    Expressions are evaluated in Stata's order, where it shows in the last
    bit: division binds tighter than multiplication, so a*b/c is a*(b/c),
    and an integer power is repeated multiplication from the left.

    License: MIT
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ctools_alloc.h"
#include "ctools_config.h"
#include "ctools_select.h"
#include "cegen_internal.h"

/* Rows per block of the single-group sums */
#define CEGEN_SUM_BLOCK 16384

/* ---------------------------------------------------------------------------
   Arithmetic as Stata does it
   --------------------------------------------------------------------------- */

/* a / b: missing when either is missing or b is zero */
static inline double st_div(double a, double b, double miss)
{
    if (!(a < miss) || !(b < miss) || b == 0.0)
        return miss;
    return a / b;
}

/* sqrt(a): missing for a missing or negative argument */
static inline double st_sqrt(double a, double miss)
{
    return (a < miss && a >= 0.0) ? sqrt(a) : miss;
}

/* d^3 and d^4 as Stata takes an integer power: multiplied from the left */
static inline double st_cube(double d)
{
    return d * d * d;
}

static inline double st_fourth(double d)
{
    return d * d * d * d;
}

/* scale*value/divisor, which Stata evaluates as scale*(value/divisor) */
static inline double st_scale_div(double scale, double value, double divisor, double miss)
{
    const double q = st_div(value, divisor, miss);
    return q < miss ? scale * q : miss;
}

/* ---------------------------------------------------------------------------
   Percentiles
   --------------------------------------------------------------------------- */

/* Mean of the i-th and (i+1)-th smallest of work (one-based, 1 <= i < n) */
static double order_pair_mean(double *work, size_t n, size_t i)
{
    double lower = ctools_quickselect_double(work, n, i - 1);
    /* quickselect leaves work[i..n-1] >= lower: the next one is their minimum */
    double upper = work[i];
    for (size_t j = i + 1; j < n; j++)
        if (work[j] < upper)
            upper = work[j];
    return (lower + upper) / 2.0;
}

double cegen_median(double *work, size_t n)
{
    if (n == 0)
        return SV_missval;
    if (n & 1)
        return ctools_quickselect_double(work, n, (n - 1) / 2);
    return order_pair_mean(work, n, n / 2);
}

double cegen_percentile(double *work, size_t n, double p, int by_formula)
{
    if (n == 0)
        return SV_missval;
    if (by_formula)
    {
        /* egen with by, on the sorted values x[1..n]:
         *   k = round(n*p/100)
         *   (x[k] + x[k+1])/2  if 100*k == n*p,  else  x[floor(n*p/100) + 1]
         * evaluated in double precision as Stata evaluates it: n*p/100 is
         * n*(p/100). */
        const double hb = (double)n * (p / 100.0);
        const double k = floor(hb + 0.5);
        if (100.0 * k == (double)n * p)
        {
            if (k < 1.0 || k >= (double)n)
                return SV_missval;
            return order_pair_mean(work, n, (size_t)k);
        }
        const double j = floor(hb);
        if (j >= (double)n)
            return SV_missval;
        return ctools_quickselect_double(work, n, (size_t)j);
    }
    /* _pctile: with h = (n*p)/100, x[floor(h) + 1] unless h is an integer,
     * when it is the mean of x[h] and x[h+1] */
    const double h = ((double)n * p) / 100.0;
    const size_t j = (size_t)floor(h);
    if (j >= n)
        return ctools_quickselect_double(work, n, n - 1);
    if (h - (double)j > 0.0 || j == 0)
        return ctools_quickselect_double(work, n, j);
    return order_pair_mean(work, n, j);
}

/* Copy the nonmissing values of v[0..count) to work; returns how many */
static size_t compact_nonmissing(const double *v, size_t count, double *work, double miss)
{
    size_t n = 0;
    for (size_t i = 0; i < count; i++)
    {
        work[n] = v[i];
        n += (v[i] < miss);
    }
    return n;
}

/* median, pctile, iqr or mad of the nonmissing values of v; work holds
 * count values and may be v itself when v can be reordered */
static double order_statistic(const cegen_job *job, const double *v, size_t count, double *work,
                              double miss, int by_formula)
{
    size_t n = compact_nonmissing(v, count, work, miss);
    switch (job->fn)
    {
    case CEGEN_FN_MEDIAN:
        return cegen_median(work, n);
    case CEGEN_FN_PCTILE:
        return cegen_percentile(work, n, job->p, by_formula);
    case CEGEN_FN_IQR:
    {
        double q1 = cegen_percentile(work, n, 25.0, by_formula);
        double q3 = cegen_percentile(work, n, 75.0, by_formula);
        return (q1 < miss && q3 < miss) ? q3 - q1 : miss;
    }
    default: /* CEGEN_FN_MAD */
    {
        double center = cegen_median(work, n);
        if (!(center < miss))
            return miss;
        for (size_t i = 0; i < n; i++)
            work[i] = fabs(work[i] - center);
        return cegen_median(work, n);
    }
    }
}

/* ---------------------------------------------------------------------------
   seq()
   --------------------------------------------------------------------------- */

typedef struct
{
    double from, to, block, range;
    int swapped;
} seq_rule;

static seq_rule seq_rule_of(const cegen_job *job)
{
    seq_rule r;
    r.from = job->seq_from;
    r.to = job->seq_to;
    r.block = job->seq_block;
    r.swapped = r.from > r.to;
    if (r.swapped)
    {
        double t = r.from;
        r.from = r.to;
        r.to = t;
    }
    r.range = r.to - r.from + 1.0;
    return r;
}

/* Value of the j-th row (zero based) of its group:
 * from + int(mod(j / block, to - from + 1)), counted down when from > to */
static inline double seq_value(const seq_rule *r, size_t j)
{
    const double q = (double)j / r->block;
    const double m = q - r->range * floor(q / r->range);
    const double v = r->from + trunc(m);
    return r->swapped ? r->to + r->from - v : v;
}

/* ---------------------------------------------------------------------------
   One group, the whole sample (egen without by)
   --------------------------------------------------------------------------- */

typedef enum
{
    SUM_VALUE = 0, /* x */
    SUM_ABSDEV,    /* |x - c| */
    SUM_DEV2,      /* (x - c)^2 */
    SUM_DEV3,      /* (x - c)^3 */
    SUM_DEV4       /* (x - c)^4 */
} sum_term;

typedef struct
{
    double sum, comp;
} sum_acc;

static inline void sum_add(sum_acc *a, double v)
{
    const double t = a->sum + v;
    if (fabs(a->sum) >= fabs(v))
        a->comp += (a->sum - t) + v;
    else
        a->comp += (v - t) + a->sum;
    a->sum = t;
}

static inline double sum_term_of(double x, sum_term term, double center)
{
    switch (term)
    {
    case SUM_VALUE:
        return x;
    case SUM_ABSDEV:
        return fabs(x - center);
    case SUM_DEV2:
    {
        const double d = x - center;
        return d * d;
    }
    case SUM_DEV3:
        return st_cube(x - center);
    default:
        return st_fourth(x - center);
    }
}

/*
    Sum of term(x[i]) over the nonmissing x and their number. Blocks are
    fixed, each is summed with compensation, and the block results are
    combined in order: the same bits for any number of threads. Returns 0
    when the block table cannot be allocated.
*/
static int blocked_sum(const double *x, size_t n, sum_term term, double center, double miss,
                       double *sum_out, size_t *count_out)
{
    const size_t nblocks = (n + CEGEN_SUM_BLOCK - 1) / CEGEN_SUM_BLOCK;
    *sum_out = 0.0;
    *count_out = 0;
    if (!nblocks)
        return 1;
    sum_acc *blocks = (sum_acc *)ctools_safe_malloc2(nblocks, sizeof(sum_acc));
    if (!blocks)
        return 0;
    size_t count = 0;
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n)) reduction(+ : count)
    for (size_t b = 0; b < nblocks; b++)
    {
        const size_t first = b * CEGEN_SUM_BLOCK;
        const size_t stop = first + CEGEN_SUM_BLOCK < n ? first + CEGEN_SUM_BLOCK : n;
        sum_acc acc = {0.0, 0.0};
        size_t present = 0;
        for (size_t i = first; i < stop; i++)
        {
            if (x[i] < miss)
            {
                sum_add(&acc, sum_term_of(x[i], term, center));
                present++;
            }
        }
        blocks[b] = acc;
        count += present;
    }
    sum_acc total = {0.0, 0.0};
    for (size_t b = 0; b < nblocks; b++)
    {
        sum_add(&total, blocks[b].sum);
        sum_add(&total, blocks[b].comp);
    }
    free(blocks);
    *sum_out = total.sum + total.comp;
    *count_out = count;
    return 1;
}

static void fill_all(double *out, size_t n, double value)
{
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n))
    for (size_t i = 0; i < n; i++)
        out[i] = value;
}

static stata_retcode stats_single(const cegen_job *job, const stata_data *data, double *out)
{
    const size_t n = data->nobs;
    const double miss = SV_missval;
    const int threads = cegen_threads_for(n);
    const stata_variable *src = job->nsrc ? &data->vars[0] : NULL;
    const double *x = (src && src->type == STATA_TYPE_DOUBLE) ? src->data.dbl : NULL;
    double sum = 0.0, value = miss;
    size_t count = 0;

    switch (job->fn)
    {
    case CEGEN_FN_COUNT:
    {
        if (x)
        {
#pragma omp parallel for schedule(static) num_threads(threads) reduction(+ : count)
            for (size_t i = 0; i < n; i++)
                count += (x[i] < miss);
        }
        else
        {
            char *const *s = src->data.str;
#pragma omp parallel for schedule(static) num_threads(threads) reduction(+ : count)
            for (size_t i = 0; i < n; i++)
                count += (s[i] && s[i][0] != '\0');
        }
        fill_all(out, n, (double)count);
        return STATA_OK;
    }
    case CEGEN_FN_SEQ:
    {
        const seq_rule rule = seq_rule_of(job);
#pragma omp parallel for schedule(static) num_threads(threads)
        for (size_t i = 0; i < n; i++)
            out[i] = seq_value(&rule, i);
        return STATA_OK;
    }
    case CEGEN_FN_MIN:
    case CEGEN_FN_MAX:
    {
        /* option missing compares the missing codes like other values */
        const int all = job->missing;
        double lo = INFINITY, hi = -INFINITY;
#pragma omp parallel for schedule(static) num_threads(threads) reduction(min : lo) reduction(max : hi) reduction(+ : count)
        for (size_t i = 0; i < n; i++)
        {
            const double v = x[i];
            if (all || v < miss)
            {
                if (v < lo)
                    lo = v;
                if (v > hi)
                    hi = v;
                count++;
            }
        }
        if (count)
            value = job->fn == CEGEN_FN_MIN ? lo : hi;
        fill_all(out, n, value);
        return STATA_OK;
    }
    case CEGEN_FN_MEDIAN:
    case CEGEN_FN_PCTILE:
    case CEGEN_FN_IQR:
    case CEGEN_FN_MAD:
    {
        double *work = (double *)ctools_safe_malloc2(n ? n : 1, sizeof(double));
        if (!work)
            return STATA_ERR_MEMORY;
        /* mad() uses the sorted-value formula with and without by */
        value = order_statistic(job, x, n, work, miss, 0);
        free(work);
        fill_all(out, n, value);
        return STATA_OK;
    }
    default:
        break;
    }

    /* The remaining functions start from the sum and the count */
    if (!blocked_sum(x, n, SUM_VALUE, 0.0, miss, &sum, &count))
        return STATA_ERR_MEMORY;
    const double mean = count ? sum / (double)count : miss;
    double m2 = 0.0, mk = 0.0;
    size_t unused = 0;

    switch (job->fn)
    {
    case CEGEN_FN_MEAN:
        value = mean;
        break;
    case CEGEN_FN_TOTAL:
        value = (job->missing && !count) ? miss : sum;
        break;
    case CEGEN_FN_SD:
    case CEGEN_FN_STD:
        if (count >= 2)
        {
            if (!blocked_sum(x, n, SUM_DEV2, mean, miss, &m2, &unused))
                return STATA_ERR_MEMORY;
            value = st_sqrt(m2 / (double)(count - 1), miss);
        }
        break;
    case CEGEN_FN_MDEV:
        if (count)
        {
            if (!blocked_sum(x, n, SUM_ABSDEV, mean, miss, &mk, &unused))
                return STATA_ERR_MEMORY;
            value = mk / (double)count;
        }
        break;
    case CEGEN_FN_SKEW:
    case CEGEN_FN_KURT:
        if (count)
        {
            const sum_term high = job->fn == CEGEN_FN_SKEW ? SUM_DEV3 : SUM_DEV4;
            if (!blocked_sum(x, n, SUM_DEV2, mean, miss, &m2, &unused) ||
                !blocked_sum(x, n, high, mean, miss, &mk, &unused))
                return STATA_ERR_MEMORY;
            m2 /= (double)count;
            mk /= (double)count;
            if (job->fn == CEGEN_FN_SKEW)
            {
                const double root = st_sqrt(m2, miss);
                value = root < miss ? st_div(mk, st_cube(root), miss) : miss;
            }
            else
                value = st_div(mk, m2 * m2, miss);
        }
        break;
    default:
        break;
    }

    if (job->fn == CEGEN_FN_STD)
    {
        /* sd()*(x - mean)/sd + mean() */
        const double sd = value, scale = job->std_sd, shift = job->std_mean;
#pragma omp parallel for schedule(static) num_threads(threads)
        for (size_t i = 0; i < n; i++)
        {
            const double z = x[i] < miss ? st_scale_div(scale, x[i] - mean, sd, miss) : miss;
            out[i] = z < miss ? z + shift : miss;
        }
        return STATA_OK;
    }
    if (job->fn == CEGEN_FN_PC)
    {
        const double scale = job->prop ? 1.0 : 100.0;
#pragma omp parallel for schedule(static) num_threads(threads)
        for (size_t i = 0; i < n; i++)
            out[i] = x[i] < miss ? st_scale_div(scale, x[i], sum, miss) : miss;
        return STATA_OK;
    }
    fill_all(out, n, value);
    return STATA_OK;
}

/* ---------------------------------------------------------------------------
   By-groups
   --------------------------------------------------------------------------- */

typedef struct
{
    double *values; /* the group's source values */
    double *work;   /* results by position, or a reorderable copy */
    size_t cap;
} group_scratch;

static int scratch_reserve(group_scratch *s, size_t count)
{
    if (count <= s->cap)
        return 1;
    size_t cap = s->cap ? s->cap : 256;
    while (cap < count)
        cap *= 2;
    double *values = (double *)realloc(s->values, cap * sizeof(double));
    if (values)
        s->values = values;
    double *work = (double *)realloc(s->work, cap * sizeof(double));
    if (work)
        s->work = work;
    if (!values || !work)
        return 0;
    s->cap = cap;
    return 1;
}

/* Sum and number of the nonmissing values, in order */
static inline double group_sum(const double *v, size_t count, double miss, size_t *present)
{
    double sum = 0.0;
    size_t n = 0;
    for (size_t i = 0; i < count; i++)
    {
        if (v[i] < miss)
        {
            sum += v[i];
            n++;
        }
    }
    *present = n;
    return sum;
}

/* Sum of (v - mean)^2 over the nonmissing values, in order */
static inline double group_ss(const double *v, size_t count, double mean, double miss)
{
    double ss = 0.0;
    for (size_t i = 0; i < count; i++)
    {
        if (v[i] < miss)
        {
            const double d = v[i] - mean;
            ss += d * d;
        }
    }
    return ss;
}

/* One value for the whole group */
static double group_scalar(const cegen_job *job, const double *v, size_t count, double miss)
{
    size_t n = 0;
    switch (job->fn)
    {
    case CEGEN_FN_COUNT:
        (void)group_sum(v, count, miss, &n);
        return (double)n;
    case CEGEN_FN_MEAN:
    {
        const double sum = group_sum(v, count, miss, &n);
        return n ? sum / (double)n : miss;
    }
    case CEGEN_FN_TOTAL:
    {
        const double sum = group_sum(v, count, miss, &n);
        return (job->missing && !n) ? miss : sum;
    }
    case CEGEN_FN_SD:
    {
        const double sum = group_sum(v, count, miss, &n);
        if (n < 2)
            return miss;
        return st_sqrt(group_ss(v, count, sum / (double)n, miss) / (double)(n - 1), miss);
    }
    case CEGEN_FN_MIN:
    case CEGEN_FN_MAX:
    {
        /* Without option missing: the extreme of the nonmissing values, or
         * system missing. With it: the extreme of all values, the missing
         * codes ordered above every number. */
        const int all = job->missing, is_min = job->fn == CEGEN_FN_MIN;
        double best = miss;
        int have = 0;
        for (size_t i = 0; i < count; i++)
        {
            const double x = v[i];
            if (!all && !(x < miss))
                continue;
            if (!have || (is_min ? x < best : x > best))
                best = x;
            have = 1;
        }
        return best;
    }
    case CEGEN_FN_MDEV:
    {
        const double sum = group_sum(v, count, miss, &n);
        if (!n)
            return miss;
        const double mean = sum / (double)n;
        double dev = 0.0;
        for (size_t i = 0; i < count; i++)
            if (v[i] < miss)
                dev += fabs(v[i] - mean);
        return dev / (double)n;
    }
    case CEGEN_FN_SKEW:
    case CEGEN_FN_KURT:
    {
        const double sum = group_sum(v, count, miss, &n);
        if (!n)
            return miss;
        const double mean = sum / (double)n;
        const int skew = job->fn == CEGEN_FN_SKEW;
        double m2 = group_ss(v, count, mean, miss) / (double)n, mk = 0.0;
        for (size_t i = 0; i < count; i++)
            if (v[i] < miss)
                mk += skew ? st_cube(v[i] - mean) : st_fourth(v[i] - mean);
        mk /= (double)n;
        if (skew)
        {
            const double root = st_sqrt(m2, miss);
            return root < miss ? st_div(mk, st_cube(root), miss) : miss;
        }
        return st_div(mk, m2 * m2, miss);
    }
    default:
        return miss;
    }
}

/* One value per row of the group, written to res[0..count) */
static void group_rows(const cegen_job *job, const seq_rule *rule, const double *v, size_t count,
                       double *res, double miss)
{
    size_t n = 0;
    if (job->fn == CEGEN_FN_SEQ)
    {
        for (size_t j = 0; j < count; j++)
            res[j] = seq_value(rule, j);
        return;
    }
    const double sum = group_sum(v, count, miss, &n);
    if (job->fn == CEGEN_FN_PC)
    {
        const double scale = job->prop ? 1.0 : 100.0;
        for (size_t j = 0; j < count; j++)
            res[j] = v[j] < miss ? st_scale_div(scale, v[j], sum, miss) : miss;
        return;
    }
    /* std */
    const double mean = n ? sum / (double)n : miss;
    const double sd =
        n >= 2 ? st_sqrt(group_ss(v, count, mean, miss) / (double)(n - 1), miss) : miss;
    const double scale = job->std_sd, shift = job->std_mean;
    for (size_t j = 0; j < count; j++)
    {
        const double z = v[j] < miss ? st_scale_div(scale, v[j] - mean, sd, miss) : miss;
        res[j] = z < miss ? z + shift : miss;
    }
}

static int fn_is_rowwise(cegen_fn fn)
{
    return fn == CEGEN_FN_STD || fn == CEGEN_FN_PC || fn == CEGEN_FN_SEQ;
}

static int fn_is_order_statistic(cegen_fn fn)
{
    return fn == CEGEN_FN_MEDIAN || fn == CEGEN_FN_PCTILE || fn == CEGEN_FN_IQR ||
           fn == CEGEN_FN_MAD;
}

/* Evaluate one group: rows at positions first .. first+count-1 of order */
static int group_eval(const cegen_job *job, const seq_rule *rule, const stata_variable *src,
                      const perm_idx_t *order, size_t first, size_t count, double *out,
                      group_scratch *s, double miss)
{
    const cegen_fn fn = job->fn;
    const double *x = (src && src->type == STATA_TYPE_DOUBLE) ? src->data.dbl : NULL;

    /* count() of a string source: the nonempty values */
    if (fn == CEGEN_FN_COUNT && !x)
    {
        char *const *str = src->data.str;
        size_t n = 0;
        for (size_t j = 0; j < count; j++)
        {
            const char *t = str[order ? order[first + j] : first + j];
            n += (t && t[0] != '\0');
        }
        for (size_t j = 0; j < count; j++)
            out[order ? order[first + j] : first + j] = (double)n;
        return 1;
    }

    /* The group's values as one contiguous run: in place when the rows are
     * already grouped, gathered through the order otherwise */
    const double *v = NULL;
    if (x)
    {
        if (order)
        {
            if (!scratch_reserve(s, count))
                return 0;
            for (size_t j = 0; j < count; j++)
                s->values[j] = x[order[first + j]];
            v = s->values;
        }
        else
            v = x + first;
    }

    if (fn_is_rowwise(fn))
    {
        if (!order)
        {
            group_rows(job, rule, v, count, out + first, miss);
            return 1;
        }
        if (!scratch_reserve(s, count))
            return 0;
        group_rows(job, rule, v, count, s->work, miss);
        for (size_t j = 0; j < count; j++)
            out[order[first + j]] = s->work[j];
        return 1;
    }

    double value;
    if (fn_is_order_statistic(fn))
    {
        if (!scratch_reserve(s, count))
            return 0;
        value = order_statistic(job, v, count, s->work, miss, 1);
    }
    else
        value = group_scalar(job, v, count, miss);

    if (order)
        for (size_t j = 0; j < count; j++)
            out[order[first + j]] = value;
    else
        for (size_t j = 0; j < count; j++)
            out[first + j] = value;
    return 1;
}

static stata_retcode stats_by(const cegen_job *job, const stata_data *data,
                              const cegen_groups *groups, double *out)
{
    const size_t n = data->nobs;
    const double miss = SV_missval;
    const stata_variable *src = job->nsrc ? &data->vars[0] : NULL;
    const seq_rule rule = seq_rule_of(job);
    int threads = cegen_threads_for(n);
    if ((size_t)threads > groups->ngroups)
        threads = groups->ngroups ? (int)groups->ngroups : 1;

    /* Chunks of consecutive groups with about equal rows, several per
     * thread so that uneven groups still balance */
    const size_t max_chunks = (size_t)threads * 8;
    size_t *cuts = (size_t *)malloc((max_chunks + 1) * sizeof(size_t));
    if (!cuts)
        return STATA_ERR_MEMORY;
    const size_t nchunks = cegen_groups_chunks(groups, max_chunks, cuts);
    int failed = 0;

#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        group_scratch scratch = {NULL, NULL, 0};
#pragma omp for schedule(dynamic, 1)
        for (size_t c = 0; c < nchunks; c++)
        {
            int stop;
#pragma omp atomic read
            stop = failed;
            if (stop)
                continue;
            for (size_t g = cuts[c]; g < cuts[c + 1]; g++)
            {
                const size_t first = groups->starts[g];
                const size_t count = groups->starts[g + 1] - first;
                if (!group_eval(job, &rule, src, groups->order, first, count, out, &scratch,
                                miss))
                {
#pragma omp atomic write
                    failed = 1;
                    break;
                }
            }
        }
        free(scratch.values);
        free(scratch.work);
    }
    free(cuts);
    return failed ? STATA_ERR_MEMORY : STATA_OK;
}

/* ---------------------------------------------------------------------------
   By-groups known by number (numeric keys, unsorted data)
   --------------------------------------------------------------------------- */

static int fn_is_moment(cegen_fn fn)
{
    return !fn_is_order_statistic(fn);
}

/*
    The functions that need no order within a group: each pass visits the
    rows once, in observation order, and adds to the row's group. The result
    of a group is then spread to its rows in parallel.
*/
static stata_retcode stats_by_ids(const cegen_job *job, const stata_data *data,
                                  const cegen_groups *groups, double *out)
{
    const size_t n = data->nobs, ngroups = groups->ngroups;
    const double miss = SV_missval;
    const ST_int *ids = groups->ids;
    const cegen_fn fn = job->fn;
    const stata_variable *src = job->nsrc ? &data->vars[0] : NULL;
    const double *x = (src && src->type == STATA_TYPE_DOUBLE) ? src->data.dbl : NULL;
    const int threads = cegen_threads_for(n);

    double *a = (double *)ctools_safe_calloc2(ngroups ? ngroups : 1, sizeof(double));
    double *b = (double *)ctools_safe_calloc2(ngroups ? ngroups : 1, sizeof(double));
    double *c = (double *)ctools_safe_calloc2(ngroups ? ngroups : 1, sizeof(double));
    uint32_t *count = (uint32_t *)ctools_safe_calloc2(ngroups ? ngroups : 1, sizeof(uint32_t));
    if (!a || !b || !c || !count)
    {
        free(a);
        free(b);
        free(c);
        free(count);
        return STATA_ERR_MEMORY;
    }

    if (fn == CEGEN_FN_SEQ)
    {
        const seq_rule rule = seq_rule_of(job);
        for (size_t i = 0; i < n; i++)
            out[i] = seq_value(&rule, count[ids[i] - 1]++);
        goto done;
    }
    if (fn == CEGEN_FN_COUNT && !x)
    {
        char *const *str = src->data.str;
        for (size_t i = 0; i < n; i++)
            count[ids[i] - 1] += (str[i] && str[i][0] != '\0');
    }
    else if (fn == CEGEN_FN_MIN || fn == CEGEN_FN_MAX)
    {
        /* a: the extreme so far; count: whether the group has one */
        const int all = job->missing, is_min = fn == CEGEN_FN_MIN;
        for (size_t i = 0; i < n; i++)
        {
            const double v = x[i];
            const size_t g = (size_t)ids[i] - 1;
            if (!all && !(v < miss))
                continue;
            if (!count[g] || (is_min ? v < a[g] : v > a[g]))
                a[g] = v;
            count[g] = 1;
        }
        for (size_t g = 0; g < ngroups; g++)
            if (!count[g])
                a[g] = miss;
    }
    else
    {
        /* a: the sum; count: the nonmissing values */
        for (size_t i = 0; i < n; i++)
        {
            if (x[i] < miss)
            {
                const size_t g = (size_t)ids[i] - 1;
                a[g] += x[i];
                count[g]++;
            }
        }
    }

    /* a becomes the group's result (for std: its mean, with the sd in b) */
    switch (fn)
    {
    case CEGEN_FN_COUNT:
        for (size_t g = 0; g < ngroups; g++)
            a[g] = (double)count[g];
        break;
    case CEGEN_FN_TOTAL:
        if (job->missing)
            for (size_t g = 0; g < ngroups; g++)
                if (!count[g])
                    a[g] = miss;
        break;
    case CEGEN_FN_MEAN:
    case CEGEN_FN_SD:
    case CEGEN_FN_STD:
    case CEGEN_FN_MDEV:
    case CEGEN_FN_SKEW:
    case CEGEN_FN_KURT:
    {
        for (size_t g = 0; g < ngroups; g++)
            a[g] = count[g] ? a[g] / (double)count[g] : miss;
        if (fn == CEGEN_FN_MEAN)
            break;
        /* b: sum of squared (mdev: absolute) deviations; c: of the 3rd or 4th power */
        const int absolute = fn == CEGEN_FN_MDEV, cube = fn == CEGEN_FN_SKEW,
                  fourth = fn == CEGEN_FN_KURT;
        for (size_t i = 0; i < n; i++)
        {
            if (!(x[i] < miss))
                continue;
            const size_t g = (size_t)ids[i] - 1;
            const double d = x[i] - a[g];
            b[g] += absolute ? fabs(d) : d * d;
            if (cube)
                c[g] += st_cube(d);
            else if (fourth)
                c[g] += st_fourth(d);
        }
        for (size_t g = 0; g < ngroups; g++)
        {
            const double m = (double)count[g];
            if (fn == CEGEN_FN_MDEV)
                a[g] = count[g] ? b[g] / m : miss;
            else if (fn == CEGEN_FN_SKEW || fn == CEGEN_FN_KURT)
            {
                if (!count[g])
                    a[g] = miss;
                else if (fn == CEGEN_FN_SKEW)
                {
                    const double root = st_sqrt(b[g] / m, miss);
                    a[g] = root < miss ? st_div(c[g] / m, st_cube(root), miss) : miss;
                }
                else
                {
                    const double m2 = b[g] / m;
                    a[g] = st_div(c[g] / m, m2 * m2, miss);
                }
            }
            else
            {
                /* sd and std: b becomes the standard deviation */
                b[g] = count[g] >= 2 ? st_sqrt(b[g] / (m - 1.0), miss) : miss;
                if (fn == CEGEN_FN_SD)
                    a[g] = b[g];
            }
        }
        break;
    }
    default:
        break;
    }

    if (fn == CEGEN_FN_STD)
    {
        const double scale = job->std_sd, shift = job->std_mean;
#pragma omp parallel for schedule(static) num_threads(threads)
        for (size_t i = 0; i < n; i++)
        {
            const size_t g = (size_t)ids[i] - 1;
            const double z = x[i] < miss ? st_scale_div(scale, x[i] - a[g], b[g], miss) : miss;
            out[i] = z < miss ? z + shift : miss;
        }
    }
    else if (fn == CEGEN_FN_PC)
    {
        const double scale = job->prop ? 1.0 : 100.0;
#pragma omp parallel for schedule(static) num_threads(threads)
        for (size_t i = 0; i < n; i++)
            out[i] = x[i] < miss ? st_scale_div(scale, x[i], a[(size_t)ids[i] - 1], miss) : miss;
    }
    else
    {
#pragma omp parallel for schedule(static) num_threads(threads)
        for (size_t i = 0; i < n; i++)
            out[i] = a[(size_t)ids[i] - 1];
    }

done:
    free(a);
    free(b);
    free(c);
    free(count);
    return STATA_OK;
}

stata_retcode cegen_stats_run(const cegen_job *job, stata_data *data, cegen_result *result,
                              cegen_report *report)
{
    const size_t n = data->nobs;
    if (job->fn != CEGEN_FN_SEQ && job->nsrc != 1)
        return STATA_ERR_INVALID_INPUT;
    if (job->nsrc && job->fn != CEGEN_FN_COUNT && data->vars[0].type != STATA_TYPE_DOUBLE)
        return STATA_ERR_UNSUPPORTED_TYPE;

    result->num = (double *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, n ? n : 1, sizeof(double));
    if (!result->num)
        return STATA_ERR_MEMORY;

    if (!job->nby)
    {
        report->ngroups = n ? 1 : 0;
        return stats_single(job, data, result->num);
    }

    /* Only the order statistics need each group's rows together */
    cegen_groups groups;
    stata_retcode rc = cegen_groups_build(data, job->nsrc, job->nby,
                                          fn_is_moment(job->fn) ? 0 : CEGEN_GROUPS_RUNS, &groups);
    if (rc == STATA_OK)
    {
        report->ngroups = groups.ngroups;
        if (groups.starts)
            rc = stats_by(job, data, &groups, result->num);
        else
            rc = stats_by_ids(job, data, &groups, result->num);
    }
    cegen_groups_free(&groups);
    return rc;
}

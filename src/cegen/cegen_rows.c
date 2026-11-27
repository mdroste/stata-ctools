/*
    cegen_rows.c
    Functions across the variables of one observation

    rowtotal, rowmean, rowsd, rowmin, rowmax, rowfirst, rowlast, rowmiss,
    rownonmiss, rowmedian, rowpctile, anycount, anymatch, anyvalue and diff.

    Each observation is independent, so the rows run in parallel over
    contiguous ranges. The arithmetic follows egen's own order of
    operations, including the ones that show in the result: rowtotal()
    accumulates in the new variable, so a float total is rounded to float
    after every addition, and rowmax(), rowfirst() and rowlast() carry a
    missing code through when every variable is missing.

    License: MIT
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ctools_alloc.h"
#include "ctools_config.h"
#include "cegen_internal.h"

/* ---------------------------------------------------------------------------
   values(): the integers to match, sorted, with how often each was listed
   --------------------------------------------------------------------------- */

typedef struct
{
    double *value;
    double *times; /* a value listed twice counts twice in anycount() */
    size_t n;
} value_set;

static int value_cmp(const void *a, const void *b)
{
    const double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static int value_set_init(value_set *set, const double *list, size_t nlist)
{
    set->n = 0;
    set->value = (double *)malloc((nlist ? nlist : 1) * sizeof(double));
    set->times = (double *)malloc((nlist ? nlist : 1) * sizeof(double));
    if (!set->value || !set->times)
        return 0;
    if (nlist)
        memcpy(set->value, list, nlist * sizeof(double));
    qsort(set->value, nlist, sizeof(double), value_cmp);
    for (size_t i = 0; i < nlist; i++)
    {
        if (set->n && set->value[set->n - 1] == set->value[i])
            set->times[set->n - 1] += 1.0;
        else
        {
            set->value[set->n] = set->value[i];
            set->times[set->n++] = 1.0;
        }
    }
    return 1;
}

static void value_set_free(value_set *set)
{
    free(set->value);
    free(set->times);
}

/* How often x is listed (0 when it is not) */
static inline double value_set_times(const value_set *set, double x)
{
    if (set->n <= 8)
    {
        for (size_t i = 0; i < set->n; i++)
            if (set->value[i] == x)
                return set->times[i];
        return 0.0;
    }
    size_t lo = 0, hi = set->n;
    while (lo < hi)
    {
        size_t mid = lo + (hi - lo) / 2;
        if (set->value[mid] < x)
            lo = mid + 1;
        else
            hi = mid;
    }
    return (lo < set->n && set->value[lo] == x) ? set->times[lo] : 0.0;
}

/* ---------------------------------------------------------------------------
   One row
   --------------------------------------------------------------------------- */

/* rowpctile()'s rule on the n sorted values: with t = n*p/100, the mean of
 * x[t] and x[t+1] when t is an integer (within 1e-8), else x[ceil(t)] */
static double row_percentile(double *sorted, size_t n, double p, double miss)
{
    if (n == 0)
        return miss;
    const double threshold = (double)n * p / 100.0;
    double index = ceil(threshold);
    if (fabs(index - threshold) < 1e-8)
        index += 1.0;
    if (index < 1.0 || index > (double)n)
        return miss;
    const size_t i = (size_t)index; /* one based */
    if (fabs(threshold - (index - 1.0)) < 1e-8)
    {
        if (i < 2)
            return miss;
        return (sorted[i - 2] + sorted[i - 1]) / 2.0;
    }
    return sorted[i - 1];
}

static inline void insertion_sort(double *a, size_t n)
{
    for (size_t i = 1; i < n; i++)
    {
        const double key = a[i];
        size_t j = i;
        while (j > 0 && a[j - 1] > key)
        {
            a[j] = a[j - 1];
            j--;
        }
        a[j] = key;
    }
}

typedef struct
{
    const cegen_job *job;
    const double **num; /* numeric columns, NULL where the column is a string */
    char ***str;        /* string columns, NULL where the column is numeric */
    size_t ncols;
    const value_set *values;
    double miss;
    int stepwise; /* rowtotal: round the running total to float */
} row_context;

static inline int cell_missing(const row_context *c, size_t col, size_t row)
{
    if (c->num[col])
        return !(c->num[col][row] < c->miss);
    const char *s = c->str[col][row];
    return !s || s[0] == '\0';
}

static double row_value(const row_context *c, size_t row, double *buffer)
{
    const size_t k = c->ncols;
    const double miss = c->miss;
    const double **x = c->num;

    switch (c->job->fn)
    {
    case CEGEN_FN_ROWTOTAL:
    {
        double g;
        if (!c->job->missing)
        {
            /* missing values add 0 */
            g = x[0][row] < miss ? x[0][row] : 0.0;
            if (c->stepwise)
                g = cegen_type_round(g, CEGEN_T_FLOAT);
            for (size_t j = 1; j < k; j++)
            {
                if (!(g < miss))
                    break;
                g += x[j][row] < miss ? x[j][row] : 0.0;
                if (c->stepwise)
                    g = cegen_type_round(g, CEGEN_T_FLOAT);
            }
            return g < miss ? g : miss;
        }
        /* option missing: missing when every variable is missing */
        g = x[0][row] < miss ? x[0][row] : miss;
        if (c->stepwise)
            g = cegen_type_round(g, CEGEN_T_FLOAT);
        for (size_t j = 1; j < k; j++)
        {
            const double v = x[j][row];
            if (!(v < miss))
                continue;
            g = g < miss ? g + v : v;
            if (c->stepwise)
                g = cegen_type_round(g, CEGEN_T_FLOAT);
        }
        return g < miss ? g : miss;
    }
    case CEGEN_FN_ROWMEAN:
    case CEGEN_FN_ROWSD:
    {
        double sum = 0.0;
        size_t n = 0;
        for (size_t j = 0; j < k; j++)
        {
            if (x[j][row] < miss)
            {
                sum += x[j][row];
                n++;
            }
        }
        if (!n)
            return miss;
        const double mean = sum / (double)n;
        if (c->job->fn == CEGEN_FN_ROWMEAN)
            return mean;
        if (n < 2)
            return miss;
        double ss = 0.0;
        for (size_t j = 0; j < k; j++)
        {
            if (x[j][row] < miss)
            {
                const double d = x[j][row] - mean;
                ss += d * d;
            }
        }
        return sqrt(ss / (double)(n - 1));
    }
    case CEGEN_FN_ROWMIN:
    {
        double g = miss;
        for (size_t j = 0; j < k; j++)
            if (x[j][row] < g)
                g = x[j][row];
        return g;
    }
    case CEGEN_FN_ROWMAX:
    {
        /* a missing result takes the code of the last variable */
        double g = x[0][row];
        for (size_t j = 1; j < k; j++)
        {
            const double v = x[j][row];
            if (!(g < miss) || (v > g && v < miss))
                g = v;
        }
        return g;
    }
    case CEGEN_FN_ROWFIRST:
    {
        /* the first nonmissing value; else the last variable's missing code */
        for (size_t j = 0; j < k; j++)
            if (x[j][row] < miss)
                return x[j][row];
        return x[k - 1][row];
    }
    case CEGEN_FN_ROWLAST:
    {
        for (size_t j = k; j-- > 0;)
            if (x[j][row] < miss)
                return x[j][row];
        return x[0][row];
    }
    case CEGEN_FN_ROWMISS:
    case CEGEN_FN_ROWNONMISS:
    {
        size_t n = 0;
        for (size_t j = 0; j < k; j++)
            n += (size_t)cell_missing(c, j, row);
        return (double)(c->job->fn == CEGEN_FN_ROWMISS ? n : k - n);
    }
    case CEGEN_FN_ROWMEDIAN:
    case CEGEN_FN_ROWPCTILE:
    {
        size_t n = 0;
        for (size_t j = 0; j < k; j++)
        {
            buffer[n] = x[j][row];
            n += (x[j][row] < miss);
        }
        if (n > 32)
            qsort(buffer, n, sizeof(double), value_cmp);
        else
            insertion_sort(buffer, n);
        return row_percentile(buffer, n, c->job->p, miss);
    }
    case CEGEN_FN_ANYCOUNT:
    {
        double n = 0.0;
        for (size_t j = 0; j < k; j++)
            n += value_set_times(c->values, x[j][row]);
        return n;
    }
    case CEGEN_FN_ANYMATCH:
    {
        for (size_t j = 0; j < k; j++)
            if (value_set_times(c->values, x[j][row]) > 0.0)
                return 1.0;
        return 0.0;
    }
    case CEGEN_FN_ANYVALUE:
        return value_set_times(c->values, x[0][row]) > 0.0 ? x[0][row] : miss;
    case CEGEN_FN_DIFF:
    {
        if (x[0])
        {
            for (size_t j = 1; j < k; j++)
                if (x[j][row] != x[0][row])
                    return 1.0;
            return 0.0;
        }
        const char *first = c->str[0][row] ? c->str[0][row] : "";
        for (size_t j = 1; j < k; j++)
            if (strcmp(first, c->str[j][row] ? c->str[j][row] : "") != 0)
                return 1.0;
        return 0.0;
    }
    default:
        return miss;
    }
}

/* ---------------------------------------------------------------------------
   Integer storage types

   egen builds rowtotal(), rowfirst(), rowlast(), rowmax() and rowmin() in
   the new variable itself: generate stores the first term, then one replace
   per variable updates it. In an integer type that makes every step show:
   generate truncates, and a replace that meets a value the type cannot hold
   widens the whole variable (to float or double for a fraction, to the next
   integer type for a large integer), after which later steps round to the
   new type. The same steps are taken here, one pass over the rows per
   variable, so the values and the final type are the ones egen arrives at.
   --------------------------------------------------------------------------- */

static int chain_function(cegen_fn fn)
{
    return fn == CEGEN_FN_ROWTOTAL || fn == CEGEN_FN_ROWFIRST || fn == CEGEN_FN_ROWLAST ||
           fn == CEGEN_FN_ROWMAX || fn == CEGEN_FN_ROWMIN;
}

/* The type of an integer variable after replace has stored values with
 * this range (fractional: some value is not an integer) */
cegen_type cegen_chain_promote(cegen_type type, int any, int fractional, double lo, double hi)
{
    static const double low[] = {-127.0, -32767.0, -2147483647.0};
    static const double high[] = {100.0, 32740.0, 2147483620.0};
    if (!any || type > CEGEN_T_LONG)
        return type;
    if (fractional)
        return type == CEGEN_T_LONG ? CEGEN_T_DOUBLE : CEGEN_T_FLOAT;
    for (int t = (int)type; t <= (int)CEGEN_T_LONG; t++)
        if (lo >= low[t] && hi <= high[t])
            return (cegen_type)t;
    return CEGEN_T_DOUBLE;
}

/* The value replace assigns at one step, from the running value g and the
 * step's variable v */
static inline double chain_step(const cegen_job *job, double g, double v, double miss)
{
    switch (job->fn)
    {
    case CEGEN_FN_ROWTOTAL:
        if (!job->missing)
            return g < miss ? g + (v < miss ? v : 0.0) : g;
        if (!(v < miss))
            return g;
        return g < miss ? g + v : v;
    case CEGEN_FN_ROWMAX:
        return (!(g < miss) || (v > g && v < miss)) ? v : g;
    case CEGEN_FN_ROWMIN:
        return v < g ? v : g;
    default: /* rowfirst, rowlast */
        return v < miss ? v : g;
    }
}

static void rows_chain(const cegen_job *job, const double **x, size_t k, size_t n, double *out,
                       cegen_type *final_type)
{
    const double miss = SV_missval;
    const int threads = cegen_threads_for(n);
    cegen_type type = job->type;

    /* generate: the first term, truncated to the type (rowmin starts from
     * missing and takes every variable in a step; rowfirst runs backwards) */
    const int backwards = job->fn == CEGEN_FN_ROWFIRST;
    const int from_missing = job->fn == CEGEN_FN_ROWMIN;
    const double *first = x[backwards ? k - 1 : 0];
    const int zero_missing = job->fn == CEGEN_FN_ROWTOTAL && !job->missing;
    const int plain_missing = job->fn == CEGEN_FN_ROWTOTAL && job->missing;
#pragma omp parallel for schedule(static) num_threads(threads)
    for (size_t i = 0; i < n; i++)
    {
        double v = from_missing ? miss : first[i];
        if (zero_missing && !(v < miss))
            v = 0.0;
        if (plain_missing && !(v < miss))
            v = miss;
        out[i] = cegen_type_round(v, type);
    }

    for (size_t step = from_missing ? 0 : 1; step < k; step++)
    {
        const double *v = x[backwards ? k - 1 - step : step];
        int fractional = 0, any = 0;
        double lo = INFINITY, hi = -INFINITY;
#pragma omp parallel for schedule(static) num_threads(threads) reduction(| : fractional, any) reduction(min : lo) reduction(max : hi)
        for (size_t i = 0; i < n; i++)
        {
            const double value = chain_step(job, out[i], v[i], miss);
            out[i] = value;
            if (value < miss)
            {
                any = 1;
                fractional |= (value != trunc(value));
                if (value < lo)
                    lo = value;
                if (value > hi)
                    hi = value;
            }
        }
        type = cegen_chain_promote(type, any, fractional, lo, hi);
        if (type == CEGEN_T_FLOAT)
        {
#pragma omp parallel for schedule(static) num_threads(threads)
            for (size_t i = 0; i < n; i++)
                out[i] = cegen_type_round(out[i], CEGEN_T_FLOAT);
        }
    }
    *final_type = type;
}

/* ---------------------------------------------------------------------------
   Driver
   --------------------------------------------------------------------------- */

stata_retcode cegen_rows_run(const cegen_job *job, stata_data *data, cegen_result *result,
                             cegen_report *report)
{
    const size_t n = data->nobs, k = job->nsrc;
    if (!k || k > data->nvars)
        return STATA_ERR_INVALID_INPUT;

    /* Strings are read only where egen reads them: for missingness, and
     * in diff() when every variable is a string */
    const int strings_ok = job->fn == CEGEN_FN_ROWMISS || job->fn == CEGEN_FN_ROWNONMISS ||
                           job->fn == CEGEN_FN_DIFF;
    size_t nstring = 0;
    for (size_t j = 0; j < k; j++)
        nstring += data->vars[j].type == STATA_TYPE_STRING;
    if (nstring && (!strings_ok || (job->fn == CEGEN_FN_DIFF && nstring != k)))
        return STATA_ERR_UNSUPPORTED_TYPE;
    if (job->fn == CEGEN_FN_ANYVALUE && k != 1)
        return STATA_ERR_INVALID_INPUT;

    stata_retcode rc = STATA_OK;
    const double **num = (const double **)calloc(k, sizeof(double *));
    char ***str = (char ***)calloc(k, sizeof(char **));
    value_set values = {NULL, NULL, 0};
    result->num = (double *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, n ? n : 1, sizeof(double));
    if (!num || !str || !result->num || !value_set_init(&values, job->list, job->nlist))
    {
        rc = STATA_ERR_MEMORY;
        goto done;
    }
    for (size_t j = 0; j < k; j++)
    {
        if (data->vars[j].type == STATA_TYPE_STRING)
            str[j] = data->vars[j].data.str;
        else
            num[j] = data->vars[j].data.dbl;
    }

    if (job->type <= CEGEN_T_LONG && chain_function(job->fn))
    {
        rows_chain(job, num, k, n, result->num, &report->final_type);
        report->has_final_type = 1;
        report->ngroups = n;
        goto done;
    }

    row_context context;
    context.job = job;
    context.num = num;
    context.str = str;
    context.ncols = k;
    context.values = &values;
    context.miss = SV_missval;
    context.stepwise = job->fn == CEGEN_FN_ROWTOTAL && job->type == CEGEN_T_FLOAT;

    const int needs_buffer = job->fn == CEGEN_FN_ROWMEDIAN || job->fn == CEGEN_FN_ROWPCTILE;
    const int threads = cegen_threads_for(n);
    double *out = result->num;
    int failed = 0;
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        double *buffer = needs_buffer ? (double *)malloc(k * sizeof(double)) : NULL;
        if (needs_buffer && !buffer)
        {
#pragma omp atomic write
            failed = 1;
        }
#pragma omp for schedule(static)
        for (size_t i = 0; i < n; i++)
        {
            if (needs_buffer && !buffer)
                continue;
            out[i] = row_value(&context, i, buffer);
        }
        free(buffer);
    }
    if (failed)
        rc = STATA_ERR_MEMORY;
    report->ngroups = n;

done:
    free(num);
    free(str);
    value_set_free(&values);
    return rc;
}

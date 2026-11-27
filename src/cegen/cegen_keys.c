/*
    cegen_keys.c
    group(), tag(), cut(), fill() and ma()

    group() and tag() number or mark the distinct values of a varlist.
    Sorted data are runs of equal keys, found and filled in parallel. For
    unsorted numeric keys the groups are numbered without a sort; tag() is
    then done, and group() sorts only one row per distinct value to put the
    numbers in key order. String keys sort the sample.
    cut() bins a variable at given cut points or at quantiles it selects
    itself; fill() extends a pattern and ma() is a centered moving average.

    License: MIT
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ctools_alloc.h"
#include "ctools_config.h"
#include "ctools_select.h"
#include "cegen_internal.h"

/* egen group(, label) defines at most this many label values */
#define CEGEN_MAX_LABELS 65536

/* ---------------------------------------------------------------------------
   group() and tag()
   --------------------------------------------------------------------------- */

/* 1 when a key of the row is missing (a string key: empty) */
static int key_missing(const stata_data *data, size_t nkeys, size_t row, double miss)
{
    for (size_t k = 0; k < nkeys; k++)
    {
        const stata_variable *v = &data->vars[k];
        if (v->type == STATA_TYPE_DOUBLE)
        {
            if (!(v->data.dbl[row] < miss))
                return 1;
        }
        else if (!v->data.str[row] || v->data.str[row][0] == '\0')
            return 1;
    }
    return 0;
}

/* Remember the observation of a group's row for group(, label) */
static stata_retcode keep_rep(cegen_report *report, size_t ngroups, perm_idx_t obs)
{
    if (report->reps_overflow)
        return STATA_OK;
    if (ngroups > CEGEN_MAX_LABELS)
    {
        report->reps_overflow = 1;
        return STATA_OK;
    }
    if (!report->reps)
    {
        report->reps = (perm_idx_t *)malloc((size_t)CEGEN_MAX_LABELS * sizeof(perm_idx_t));
        if (!report->reps)
            return STATA_ERR_MEMORY;
    }
    report->reps[report->nreps++] = obs;
    return STATA_OK;
}

/*
    The numbered groups in ascending key order: sorts one row of each group
    (first[g], its earliest row) by the keys. Returns the permutation of the
    groups in *order_out (aligned; the caller frees it).
*/
static stata_retcode ids_key_order(const stata_data *data, size_t nkeys, const perm_idx_t *first,
                                   size_t ngroups, perm_idx_t **order_out)
{
    stata_data small;
    stata_data_init(&small);
    stata_retcode rc = STATA_OK;
    int *sort_vars = (int *)malloc(nkeys * sizeof(int));
    small.vars = (stata_variable *)ctools_safe_cacheline_alloc2(nkeys, sizeof(stata_variable));
    small.sort_order = (perm_idx_t *)ctools_safe_cacheline_alloc2(ngroups, sizeof(perm_idx_t));
    if (!sort_vars || !small.vars || !small.sort_order)
        rc = STATA_ERR_MEMORY;
    if (small.vars)
        memset(small.vars, 0, nkeys * sizeof(stata_variable));
    small.nobs = ngroups;
    small.nvars = small.vars ? nkeys : 0;
    for (size_t k = 0; k < nkeys && rc == STATA_OK; k++)
    {
        const double *column = data->vars[k].data.dbl;
        double *values = (double *)ctools_safe_cacheline_alloc2(ngroups, sizeof(double));
        if (!values)
        {
            rc = STATA_ERR_MEMORY;
            break;
        }
        small.vars[k].type = STATA_TYPE_DOUBLE;
        small.vars[k].nobs = ngroups;
        small.vars[k].data.dbl = values;
        sort_vars[k] = (int)k + 1;
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(ngroups))
        for (size_t g = 0; g < ngroups; g++)
            values[g] = column[first[g]];
    }
    if (rc == STATA_OK)
    {
        for (size_t g = 0; g < ngroups; g++)
            small.sort_order[g] = (perm_idx_t)g;
        if (ngroups > 1)
            rc = ctools_sort_dispatch(&small, sort_vars, nkeys, SORT_ALG_AUTO);
    }
    if (rc == STATA_OK)
    {
        *order_out = small.sort_order;
        small.sort_order = NULL;
    }
    free(sort_vars);
    stata_data_free(&small);
    return rc;
}

/*
    group() and tag() over numbered groups. Returns STATA_ERR_UNSUPPORTED_TYPE
    when the groups are too many for sorting their keys to pay (the caller
    sorts the sample instead).
*/
static stata_retcode keys_from_ids(const cegen_job *job, const stata_data *data,
                                   const perm_idx_t *obs_map, const cegen_groups *groups,
                                   double *out, cegen_report *report)
{
    const size_t n = data->nobs, nkeys = job->nsrc, total = groups->ngroups;
    const double miss = SV_missval;
    const ST_int *ids = groups->ids;
    const int tag = job->fn == CEGEN_FN_TAG;
    const int threads = cegen_threads_for(n);
    if (!tag && total > n / 4)
        return STATA_ERR_UNSUPPORTED_TYPE;

    /* The earliest row of each group: descending, so the lowest row wins */
    perm_idx_t *first = (perm_idx_t *)ctools_safe_malloc2(total ? total : 1, sizeof(perm_idx_t));
    double *number = (double *)ctools_safe_malloc2(total ? total : 1, sizeof(double));
    perm_idx_t *key_order = NULL;
    stata_retcode rc = STATA_OK;
    if (!first || !number)
    {
        free(first);
        free(number);
        return STATA_ERR_MEMORY;
    }
    for (size_t i = n; i-- > 0;)
        first[ids[i] - 1] = (perm_idx_t)i;

    size_t ngroups = 0;
    if (tag)
    {
#pragma omp parallel for schedule(static) num_threads(threads)
        for (size_t i = 0; i < n; i++)
            out[i] = 0.0;
        for (size_t g = 0; g < total; g++)
        {
            if (!job->missing && key_missing(data, nkeys, first[g], miss))
                continue;
            out[first[g]] = 1.0;
            ngroups++;
        }
    }
    else
    {
        rc = ids_key_order(data, nkeys, first, total, &key_order);
        for (size_t pos = 0; pos < total && rc == STATA_OK; pos++)
        {
            const size_t g = key_order[pos];
            if (!job->missing && key_missing(data, nkeys, first[g], miss))
            {
                number[g] = miss;
                continue;
            }
            number[g] = (double)(++ngroups);
            if (job->label)
                rc = keep_rep(report, ngroups, obs_map[first[g]]);
        }
        if (rc == STATA_OK)
        {
#pragma omp parallel for schedule(static) num_threads(threads)
            for (size_t i = 0; i < n; i++)
                out[i] = number[ids[i] - 1];
        }
    }
    report->ngroups = ngroups;
    ctools_aligned_free(key_order);
    free(first);
    free(number);
    return rc;
}

static stata_retcode keys_group_tag(const cegen_job *job, stata_data *data,
                                    const perm_idx_t *obs_map, double *out, cegen_report *report)
{
    const size_t n = data->nobs, nkeys = job->nsrc;
    const double miss = SV_missval;
    const int tag = job->fn == CEGEN_FN_TAG;

    if (!nkeys)
    {
        /* group() of no variables: one group */
        for (size_t i = 0; i < n; i++)
            out[i] = 1.0;
        report->ngroups = n ? 1 : 0;
        return STATA_OK;
    }

    /* Runs when the data are sorted (or the keys are strings), else numbers */
    cegen_groups runs;
    stata_retcode rc = cegen_groups_build(data, 0, nkeys, 0, &runs);
    if (rc == STATA_OK && !runs.starts)
    {
        rc = keys_from_ids(job, data, obs_map, &runs, out, report);
        cegen_groups_free(&runs);
        if (rc != STATA_ERR_UNSUPPORTED_TYPE)
            return rc;
        rc = cegen_groups_build(data, 0, nkeys, CEGEN_GROUPS_RUNS | CEGEN_GROUPS_KEY_ORDER,
                                &runs);
    }
    if (rc != STATA_OK)
    {
        cegen_groups_free(&runs);
        return rc;
    }
    const perm_idx_t *order = runs.order;

    /* The number of each run: 0 for a run of missing keys, which stay out
     * unless option missing was given. Runs are in key order. */
    double *id = (double *)ctools_safe_malloc2(runs.ngroups ? runs.ngroups : 1, sizeof(double));
    int threads = cegen_threads_for(n);
    if ((size_t)threads > runs.ngroups)
        threads = runs.ngroups ? (int)runs.ngroups : 1;
    const size_t max_chunks = (size_t)threads * 8;
    size_t *cuts = (size_t *)malloc((max_chunks + 1) * sizeof(size_t));
    if (!id || !cuts)
    {
        free(id);
        free(cuts);
        cegen_groups_free(&runs);
        return STATA_ERR_MEMORY;
    }
    size_t ngroups = 0;
    for (size_t g = 0; g < runs.ngroups; g++)
    {
        const size_t row = order ? order[runs.starts[g]] : runs.starts[g];
        if (!job->missing && key_missing(data, nkeys, row, miss))
        {
            id[g] = 0.0;
            continue;
        }
        id[g] = (double)(++ngroups);
        if (job->label && !tag)
        {
            rc = keep_rep(report, ngroups, obs_map[row]);
            if (rc != STATA_OK)
                break;
        }
    }
    report->ngroups = ngroups;

    if (rc == STATA_OK)
    {
        const size_t nchunks = cegen_groups_chunks(&runs, max_chunks, cuts);
        if (tag)
        {
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n))
            for (size_t i = 0; i < n; i++)
                out[i] = 0.0;
        }
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if (threads > 1)
        for (size_t c = 0; c < nchunks; c++)
        {
            for (size_t g = cuts[c]; g < cuts[c + 1]; g++)
            {
                const size_t first = runs.starts[g], stop = runs.starts[g + 1];
                if (tag)
                {
                    if (id[g] == 0.0)
                        continue;
                    /* the earliest observation of the run */
                    size_t lowest = order ? order[first] : first;
                    if (order)
                        for (size_t pos = first + 1; pos < stop; pos++)
                            if (order[pos] < lowest)
                                lowest = order[pos];
                    out[lowest] = 1.0;
                    continue;
                }
                const double value = id[g] == 0.0 ? miss : id[g];
                if (order)
                    for (size_t pos = first; pos < stop; pos++)
                        out[order[pos]] = value;
                else
                    for (size_t pos = first; pos < stop; pos++)
                        out[pos] = value;
            }
        }
    }
    free(id);
    free(cuts);
    cegen_groups_free(&runs);
    return rc;
}

/* ---------------------------------------------------------------------------
   cut()
   --------------------------------------------------------------------------- */

/*
    Quantile k of nq of the m sorted values, as pctile, nq() takes it: with
    h = m*k/nq, the mean of x[h] and x[h+1] when h is an integer, else
    x[floor(h) + 1] (one based). Returns the zero-based position of the
    value and sets *pair when the next value is averaged in.
*/
static size_t cut_quantile_position(size_t m, size_t k, size_t nq, int *pair)
{
    const unsigned long long t = (unsigned long long)m * (unsigned long long)k;
    const size_t whole = (size_t)(t / nq);
    *pair = (t % nq == 0);
    if (*pair)
        return whole ? whole - 1 : 0;
    return whole < m ? whole : m - 1;
}

static stata_retcode keys_cut(const cegen_job *job, const stata_data *data, double *out,
                              cegen_report *report)
{
    const size_t n = data->nobs;
    const double miss = SV_missval;
    const double *x = data->vars[0].data.dbl;
    const int xfloat = job->xfloat;
    double *cuts = NULL;
    size_t ncuts = 0;

    if (job->nq > 0)
    {
        /* group(#): the minimum, the # - 1 quantiles, and the maximum + 1 */
        const size_t nq = (size_t)job->nq;
        double *work = (double *)ctools_safe_malloc2(n ? n : 1, sizeof(double));
        cuts = (double *)ctools_safe_malloc2(nq + 1, sizeof(double));
        if (!work || !cuts)
        {
            free(work);
            free(cuts);
            return STATA_ERR_MEMORY;
        }
        size_t m = 0;
        double lo = INFINITY, hi = -INFINITY;
        for (size_t i = 0; i < n; i++)
        {
            const double v = xfloat ? cegen_type_round(x[i], CEGEN_T_FLOAT) : x[i];
            if (!(v < miss))
                continue;
            work[m++] = v;
            if (v < lo)
                lo = v;
            if (v > hi)
                hi = v;
        }
        if ((double)nq > job->nobs_total + 1.0)
        {
            /* pctile refuses more quantiles than _N + 1 (the ado checks
             * this first; the limit is not the size of the sample) */
            report->bad_groups = 1;
            m = 0;
        }
        else if (!m)
            report->all_missing = 1; /* pctile: no observations, r(2000) */
        if (m)
        {
            cuts[ncuts++] = lo;
            /* Ascending quantiles: each selection leaves everything before
             * its position no larger, so the next one only works to its right */
            size_t done = 0;
            for (size_t k = 1; k < nq; k++)
            {
                int pair = 0;
                const size_t pos = cut_quantile_position(m, k, nq, &pair);
                if (pos < done)
                {
                    /* the same value as the quantile before */
                    cuts[ncuts] = cuts[ncuts - 1];
                    ncuts++;
                    continue;
                }
                double q = ctools_quickselect_double(work + done, m - done, pos - done);
                if (pair && pos + 1 < m)
                {
                    double next = work[pos + 1];
                    for (size_t j = pos + 2; j < m; j++)
                        if (work[j] < next)
                            next = work[j];
                    q = (q + next) / 2.0;
                }
                cuts[ncuts++] = q;
                done = pos;
            }
            cuts[ncuts++] = hi + 1.0;
        }
        free(work);
    }
    else
    {
        ncuts = job->nlist;
        cuts = (double *)ctools_safe_malloc2(ncuts ? ncuts : 1, sizeof(double));
        if (!cuts)
            return STATA_ERR_MEMORY;
        if (ncuts)
            memcpy(cuts, job->list, ncuts * sizeof(double));
    }

    /* at(): values at or above the last cut point are left out */
    const double top = (job->nq > 0 || !ncuts) ? miss : cuts[ncuts - 1];
    const int icodes = job->icodes;
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n))
    for (size_t i = 0; i < n; i++)
    {
        const double v = xfloat ? cegen_type_round(x[i], CEGEN_T_FLOAT) : x[i];
        if (!(v < miss) || v >= top || !ncuts || v < cuts[0])
        {
            out[i] = miss;
            continue;
        }
        /* the last cut point at or below v */
        size_t lo = 0, hi = ncuts;
        while (hi - lo > 1)
        {
            const size_t mid = lo + (hi - lo) / 2;
            if (cuts[mid] <= v)
                lo = mid;
            else
                hi = mid;
        }
        out[i] = icodes ? (double)lo : cuts[lo];
    }
    report->cuts = cuts;
    report->ncuts = ncuts;
    report->ngroups = ncuts;
    return STATA_OK;
}

/* ---------------------------------------------------------------------------
   fill()
   --------------------------------------------------------------------------- */

static stata_retcode keys_fill(const cegen_job *job, size_t n, double *out)
{
    const double miss = SV_missval;
    const size_t k = job->nlist;
    if (k < 2 || k > n)
        return STATA_ERR_INVALID_INPUT;
    /* An integer type that cannot hold a value is promoted by the ado, so
     * only float rounds the running values */
    const int round_out = job->type == CEGEN_T_FLOAT;
    const int round_step = job->xfloat;
    double *step = (double *)ctools_safe_malloc2(n, sizeof(double));
    if (!step)
        return STATA_ERR_MEMORY;

    for (size_t i = 0; i < n; i++)
        out[i] = i < k ? (round_out ? cegen_type_round(job->list[i], CEGEN_T_FLOAT) : job->list[i])
                       : miss;
    /* The differences of the pattern, repeated from its second half on */
    for (size_t i = 0; i < n; i++)
    {
        double d = miss;
        if (i + 1 < n && out[i] < miss && out[i + 1] < miss)
        {
            d = out[i + 1] - out[i];
            if (round_step)
                d = cegen_type_round(d, CEGEN_T_FLOAT);
        }
        step[i] = d;
    }
    const size_t half = (size_t)floor((double)k / 2.0 + 0.5);
    for (size_t i = half; i < n; i++)
        step[i] = step[i - half];
    for (size_t i = half; i < n; i++)
    {
        double v = miss;
        if (out[i - 1] < miss && step[i - 1] < miss)
        {
            v = out[i - 1] + step[i - 1];
            if (round_out)
                v = cegen_type_round(v, CEGEN_T_FLOAT);
        }
        out[i] = v;
    }
    free(step);
    return STATA_OK;
}

/* ---------------------------------------------------------------------------
   ma()
   --------------------------------------------------------------------------- */

/* After one of ma()'s replace steps: widen the type as replace does for the
 * values now held, and round them when the type has become float */
static void ma_settle(double *out, size_t n, int threads, cegen_type *type)
{
    const double miss = SV_missval;
    int fractional = 0, any = 0;
    double lo = INFINITY, hi = -INFINITY;
#pragma omp parallel for schedule(static) num_threads(threads) reduction(| : fractional, any) reduction(min : lo) reduction(max : hi)
    for (size_t i = 0; i < n; i++)
    {
        const double value = out[i];
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
    *type = cegen_chain_promote(*type, any, fractional, lo, hi);
    if (*type == CEGEN_T_FLOAT)
    {
#pragma omp parallel for schedule(static) num_threads(threads)
        for (size_t i = 0; i < n; i++)
            out[i] = cegen_type_round(out[i], CEGEN_T_FLOAT);
    }
}

/*
    ma() in an integer type. egen creates the variable empty in the type
    asked for and fills it by replace: the values, then one replace per
    distance adding both neighbours, then the division by their count. Each
    step can widen the type (cegen_chain_promote), and later steps round to
    the wider type, so the steps are taken one pass at a time here and the
    type reached is the type of the result.
*/
static void ma_chain(const cegen_job *job, const double *x, size_t n, double *out,
                     cegen_type *final_type)
{
    const double miss = SV_missval;
    const int threads = cegen_threads_for(n);
    const size_t half = (size_t)((job->window - 1) / 2);
    cegen_type type = job->type;

#pragma omp parallel for schedule(static) num_threads(threads)
    for (size_t i = 0; i < n; i++)
        out[i] = x[i];
    ma_settle(out, n, threads, &type);

    for (size_t d = 1; d <= half; d++)
    {
#pragma omp parallel for schedule(static) num_threads(threads)
        for (size_t i = 0; i < n; i++)
        {
            /* neighbours outside the range are absent, not missing */
            const double lag = i >= d ? x[i - d] : 0.0, lead = i + d < n ? x[i + d] : 0.0;
            if (out[i] < miss)
                out[i] = out[i] + lag + lead;
        }
        ma_settle(out, n, threads, &type);
    }

#pragma omp parallel for schedule(static) num_threads(threads)
    for (size_t i = 0; i < n; i++)
    {
        const size_t count = 1 + (i < half ? i : half) + (n - 1 - i < half ? n - 1 - i : half);
        if (out[i] < miss)
            out[i] /= (double)count;
    }
    ma_settle(out, n, threads, &type);

    if (!job->nomiss)
    {
#pragma omp parallel for schedule(static) num_threads(threads)
        for (size_t i = 0; i < n; i++)
        {
            const size_t count = 1 + (i < half ? i : half) + (n - 1 - i < half ? n - 1 - i : half);
            if (count != (size_t)job->window)
                out[i] = miss;
        }
    }
    *final_type = type;
}

static stata_retcode keys_ma(const cegen_job *job, const stata_data *data, double *out,
                             cegen_report *report)
{
    const size_t n = data->nobs;
    const double miss = SV_missval;
    const double *x = data->vars[0].data.dbl;
    const size_t half = (size_t)((job->window - 1) / 2);
    const int round_out = job->type == CEGEN_T_FLOAT;

    int has_missing = 0;
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n)) reduction(| : has_missing)
    for (size_t i = 0; i < n; i++)
        has_missing |= !(x[i] < miss);
    if (has_missing)
    {
        /* egen stops: "missing values encountered" (its ma() runs under
         * version 6, where every missing code equals missing) */
        report->all_missing = 1;
        return STATA_OK;
    }

    if (job->type <= CEGEN_T_LONG)
    {
        ma_chain(job, x, n, out, &report->final_type);
        report->has_final_type = 1;
        return STATA_OK;
    }

#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n))
    for (size_t i = 0; i < n; i++)
    {
        double v = x[i];
        size_t count = 1;
        if (round_out)
            v = cegen_type_round(v, CEGEN_T_FLOAT);
        for (size_t d = 1; d <= half; d++)
        {
            /* neighbours outside the range are absent, not missing */
            const int has_lag = i >= d, has_lead = i + d < n;
            const double lag = has_lag ? x[i - d] : 0.0, lead = has_lead ? x[i + d] : 0.0;
            if (!(v < miss))
                break; /* the running value overflowed its type */
            v = v + lag + lead;
            if (round_out)
                v = cegen_type_round(v, CEGEN_T_FLOAT);
            count += (size_t)has_lag + (size_t)has_lead;
        }
        if (v < miss)
        {
            v /= (double)count;
            if (round_out)
                v = cegen_type_round(v, CEGEN_T_FLOAT);
            if (!job->nomiss && count != (size_t)job->window)
                v = miss;
        }
        out[i] = v;
    }
    return STATA_OK;
}

/* ---------------------------------------------------------------------------
   Driver
   --------------------------------------------------------------------------- */

stata_retcode cegen_keys_run(const cegen_job *job, stata_data *data, const perm_idx_t *obs_map,
                             cegen_result *result, cegen_report *report)
{
    const size_t n = data->nobs;
    result->num = (double *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, n ? n : 1, sizeof(double));
    if (!result->num)
        return STATA_ERR_MEMORY;

    switch (job->fn)
    {
    case CEGEN_FN_GROUP:
    case CEGEN_FN_TAG:
        return keys_group_tag(job, data, obs_map, result->num, report);
    case CEGEN_FN_CUT:
        if (job->nsrc != 1 || data->vars[0].type != STATA_TYPE_DOUBLE)
            return STATA_ERR_UNSUPPORTED_TYPE;
        return keys_cut(job, data, result->num, report);
    case CEGEN_FN_FILL:
        return keys_fill(job, n, result->num);
    case CEGEN_FN_MA:
        if (job->nsrc != 1 || data->vars[0].type != STATA_TYPE_DOUBLE || job->window < 3)
            return STATA_ERR_UNSUPPORTED_TYPE;
        return keys_ma(job, data, result->num, report);
    default:
        return STATA_ERR_INVALID_INPUT;
    }
}

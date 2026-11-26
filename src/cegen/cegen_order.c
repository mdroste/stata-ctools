/*
    cegen_order.c
    Functions of the sorted values within by-groups: rank() and mode()

    Both need each group's values in order. Two ways to get them:

    - By-groups of moderate size: the groups are found first (runs of sorted
      data, or numbered without a sort), and each group's values are sorted
      on their own, the groups in parallel. Sorting many short runs is far
      cheaper than one sort of the whole sample by (by-variables, value).

    - No by-variables, a few very large groups, or a string source: one
      stable sort of the sample by the by-variables and the source together
      (the parallel sorts csort uses; no sort when the rows are already in
      that order).

    Either way a group is then a run of positions whose values ascend,
    missing values last, and a tie or a distinct value is a run of equal
    neighbours: ranks and frequencies come from one pass.

    License: MIT
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ctools_alloc.h"
#include "ctools_config.h"
#include "ctools_order.h"
#include "cegen_internal.h"

/* Groups larger than this are left to the parallel sort of the sample */
#define CEGEN_ORDER_GROUP_MAX 1000000

/* ---------------------------------------------------------------------------
   rank()
   --------------------------------------------------------------------------- */

/* Ranks of one group: positions first .. first+count-1 of order, ascending */
static void rank_group(const cegen_job *job, const double *x, const perm_idx_t *order,
                       size_t first, size_t count, double *out, double miss)
{
    const size_t stop = first + count;
    size_t pos = first;
    while (pos < stop)
    {
        const size_t row = order ? order[pos] : pos;
        const double v = x[row];
        if (!(v < miss))
            break; /* missing values sort last and get no rank */
        if (job->unique)
        {
            out[row] = (double)(pos - first + 1);
            pos++;
            continue;
        }
        /* the run of values tied with v */
        size_t end = pos + 1;
        while (end < stop && x[order ? order[end] : end] == v)
            end++;
        /* lowest rank of the tie; the default adds (ties - 1)/2 for the mean */
        double rank = (double)(pos - first + 1);
        if (!job->track && !job->field)
            rank += (double)(end - pos - 1) / 2.0;
        for (size_t j = pos; j < end; j++)
            out[order ? order[j] : j] = rank;
        pos = end;
    }
    for (; pos < stop; pos++)
        out[order ? order[pos] : pos] = miss;
}

/* ---------------------------------------------------------------------------
   mode()
   --------------------------------------------------------------------------- */

typedef struct
{
    const stata_variable *src;
    const perm_idx_t *order;
    const size_t *runs; /* starts of the runs of equal values */
    double miss;
} mode_view;

static inline size_t mode_run_row(const mode_view *v, size_t run)
{
    const size_t pos = v->runs[run];
    return v->order ? v->order[pos] : pos;
}

static inline int mode_run_missing(const mode_view *v, size_t run)
{
    const size_t row = mode_run_row(v, run);
    if (v->src->type == STATA_TYPE_DOUBLE)
        return !(v->src->data.dbl[row] < v->miss);
    const char *s = v->src->data.str[row];
    return !s || s[0] == '\0';
}

typedef struct
{
    int group_missing, multiple, any_value;
    size_t max_modes, bad_groups, strlen_max;
} mode_tally;

static void mode_tally_merge(mode_tally *into, const mode_tally *from)
{
    into->group_missing |= from->group_missing;
    into->multiple |= from->multiple;
    into->any_value |= from->any_value;
    if (from->max_modes > into->max_modes)
        into->max_modes = from->max_modes;
    if (from->strlen_max > into->strlen_max)
        into->strlen_max = from->strlen_max;
    into->bad_groups += from->bad_groups;
}

/*
    The mode of one group, whose distinct values are runs run_first ..
    run_stop-1 in ascending order. Returns the run holding the mode, or
    SIZE_MAX when the group gets a missing result.
*/
static size_t mode_group(const cegen_job *job, const mode_view *v, size_t run_first,
                         size_t run_stop, mode_tally *tally)
{
    /* Highest frequency among the values that count: missing values only
     * with option missing */
    size_t best = 0, nmodes = 0, first_mode = SIZE_MAX, last_mode = SIZE_MAX;
    for (size_t r = run_first; r < run_stop; r++)
    {
        if (!job->missing && mode_run_missing(v, r))
            continue;
        const size_t freq = v->runs[r + 1] - v->runs[r];
        if (freq > best)
        {
            best = freq;
            nmodes = 0;
            first_mode = r;
        }
        if (freq == best)
        {
            nmodes++;
            last_mode = r;
        }
    }
    const long long want = job->nummode;
    if (!nmodes)
    {
        /* Only missing values, which egen gives frequency 0: every distinct
         * missing code ties, and the result is a missing value either way */
        tally->group_missing = 1;
        nmodes = run_stop - run_first;
        if (nmodes > tally->max_modes)
            tally->max_modes = nmodes;
        if (!job->minmode && !job->maxmode && !want)
        {
            tally->multiple = 1;
            return SIZE_MAX;
        }
        if (job->maxmode)
            return run_stop - 1;
        if (want > 1)
        {
            if ((size_t)want > nmodes)
            {
                tally->bad_groups++;
                return SIZE_MAX;
            }
            return run_first + (size_t)want - 1;
        }
        return run_first;
    }
    tally->any_value = 1;
    if (nmodes > tally->max_modes)
        tally->max_modes = nmodes;
    if (job->maxmode)
        return last_mode;
    if (job->minmode || want == 1)
        return first_mode;
    if (want > 1)
    {
        if ((size_t)want > nmodes)
        {
            tally->bad_groups++;
            return SIZE_MAX;
        }
        /* the want-th mode in ascending order */
        size_t seen = 0;
        for (size_t r = first_mode; r <= last_mode; r++)
        {
            if (!job->missing && mode_run_missing(v, r))
                continue;
            if (v->runs[r + 1] - v->runs[r] == best && ++seen == (size_t)want)
                return r;
        }
        return SIZE_MAX;
    }
    if (nmodes > 1)
    {
        tally->multiple = 1;
        return SIZE_MAX;
    }
    return first_mode;
}

/* ---------------------------------------------------------------------------
   Groups sorted one by one (numeric source)
   --------------------------------------------------------------------------- */

typedef struct
{
    double value;
    perm_idx_t row;
} order_item;

/* Ascending by value, the missing codes last; ties in observation order */
static int order_item_cmp(const void *a, const void *b)
{
    const order_item *x = (const order_item *)a, *y = (const order_item *)b;
    if (x->value < y->value)
        return -1;
    if (x->value > y->value)
        return 1;
    return (x->row > y->row) - (x->row < y->row);
}

static void order_items_sort(order_item *items, size_t count)
{
    if (count > 24)
    {
        qsort(items, count, sizeof(order_item), order_item_cmp);
        return;
    }
    for (size_t i = 1; i < count; i++)
    {
        const order_item key = items[i];
        size_t j = i;
        while (j > 0 && order_item_cmp(&items[j - 1], &key) > 0)
        {
            items[j] = items[j - 1];
            j--;
        }
        items[j] = key;
    }
}

typedef struct
{
    order_item *items;
    perm_idx_t *rows; /* the group's rows in value order */
    size_t *runs;     /* mode(): starts of the runs of equal values */
    size_t cap;
} order_scratch;

static int order_scratch_reserve(order_scratch *s, size_t count, int with_runs)
{
    if (count <= s->cap)
        return 1;
    size_t cap = s->cap ? s->cap : 256;
    while (cap < count)
        cap *= 2;
    order_item *items = (order_item *)realloc(s->items, cap * sizeof(order_item));
    if (items)
        s->items = items;
    perm_idx_t *rows = (perm_idx_t *)realloc(s->rows, cap * sizeof(perm_idx_t));
    if (rows)
        s->rows = rows;
    size_t *runs = s->runs;
    if (with_runs)
    {
        runs = (size_t *)realloc(s->runs, (cap + 1) * sizeof(size_t));
        if (runs)
            s->runs = runs;
    }
    if (!items || !rows || (with_runs && !runs))
        return 0;
    s->cap = cap;
    return 1;
}

static stata_retcode order_by_groups(const cegen_job *job, const stata_data *data,
                                     const cegen_groups *groups, double *out, mode_tally *total)
{
    const size_t n = data->nobs;
    const double miss = SV_missval;
    const stata_variable *src = &data->vars[0];
    const double *x = src->data.dbl;
    const perm_idx_t *order = groups->order;
    const int is_mode = job->fn == CEGEN_FN_MODE;

    int threads = cegen_threads_for(n);
    if ((size_t)threads > groups->ngroups)
        threads = groups->ngroups ? (int)groups->ngroups : 1;
    const size_t max_chunks = (size_t)threads * 8;
    size_t *cuts = (size_t *)malloc((max_chunks + 1) * sizeof(size_t));
    if (!cuts)
        return STATA_ERR_MEMORY;
    const size_t nchunks = cegen_groups_chunks(groups, max_chunks, cuts);
    int failed = 0;

#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        order_scratch s = {NULL, NULL, NULL, 0};
        mode_tally tally = {0, 0, 0, 0, 0, 0};
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
                if (!order_scratch_reserve(&s, count, is_mode))
                {
#pragma omp atomic write
                    failed = 1;
                    break;
                }
                for (size_t j = 0; j < count; j++)
                {
                    const perm_idx_t row = order ? order[first + j] : (perm_idx_t)(first + j);
                    s.items[j].value = x[row];
                    s.items[j].row = row;
                }
                order_items_sort(s.items, count);
                for (size_t j = 0; j < count; j++)
                    s.rows[j] = s.items[j].row;
                if (!is_mode)
                {
                    rank_group(job, x, s.rows, 0, count, out, miss);
                    continue;
                }
                size_t nruns = 0;
                for (size_t j = 0; j < count; j++)
                    if (j == 0 || s.items[j].value != s.items[j - 1].value)
                        s.runs[nruns++] = j;
                s.runs[nruns] = count;
                const mode_view view = {src, s.rows, s.runs, miss};
                const size_t chosen = mode_group(job, &view, 0, nruns, &tally);
                const double value = chosen != SIZE_MAX ? s.items[s.runs[chosen]].value : miss;
                for (size_t j = 0; j < count; j++)
                    out[s.rows[j]] = value;
            }
        }
#pragma omp critical(cegen_order_tally)
        mode_tally_merge(total, &tally);
        free(s.items);
        free(s.rows);
        free(s.runs);
    }
    free(cuts);
    return failed ? STATA_ERR_MEMORY : STATA_OK;
}

/* ---------------------------------------------------------------------------
   The sample sorted at once
   --------------------------------------------------------------------------- */

static stata_retcode order_by_sort(const cegen_job *job, stata_data *data, cegen_result *result,
                                   cegen_report *report, mode_tally *total)
{
    const size_t n = data->nobs, nby = job->nby, nkeys = nby + 1;
    const double miss = SV_missval;
    const stata_variable *src = &data->vars[0];
    const int is_string = src->type == STATA_TYPE_STRING;

    stata_retcode rc = STATA_OK;
    int *keys = (int *)malloc(nkeys * sizeof(int));
    int *sort_vars = (int *)malloc(nkeys * sizeof(int));
    size_t *outer = NULL, *runs = NULL, *first_run = NULL, *cuts = NULL;
    if (!keys || !sort_vars)
    {
        rc = STATA_ERR_MEMORY;
        goto done;
    }
    /* by-variables, then the source */
    for (size_t k = 0; k < nby; k++)
        keys[k] = (int)(job->nsrc + k);
    keys[nby] = 0;
    for (size_t k = 0; k < nkeys; k++)
        sort_vars[k] = keys[k] + 1;

    const perm_idx_t *order = NULL;
    if (!cegen_rows_sorted(data, keys, nkeys))
    {
        if (!data->sort_order)
        {
            rc = STATA_ERR_INVALID_INPUT;
            goto done;
        }
        rc = ctools_sort_dispatch(data, sort_vars, nkeys, SORT_ALG_AUTO);
        if (rc != STATA_OK)
            goto done;
        order = data->sort_order;
    }

    /* The by-groups */
    cegen_groups groups;
    memset(&groups, 0, sizeof(groups));
    outer = (size_t *)ctools_safe_malloc2(nby ? n + 1 : 2, sizeof(size_t));
    if (!outer)
    {
        rc = STATA_ERR_MEMORY;
        goto done;
    }
    if (nby)
        groups.ngroups = ctools_group_starts_parallel(data, order, keys, nby, outer);
    else
    {
        groups.ngroups = 1;
        outer[0] = 0;
        outer[1] = n;
    }
    groups.starts = outer;
    groups.order = order;
    report->ngroups = groups.ngroups;

    int threads = cegen_threads_for(n);
    if ((size_t)threads > groups.ngroups)
        threads = (int)groups.ngroups;
    const size_t max_chunks = (size_t)threads * 8;
    cuts = (size_t *)malloc((max_chunks + 1) * sizeof(size_t));
    if (!cuts)
    {
        rc = STATA_ERR_MEMORY;
        goto done;
    }
    const size_t nchunks = cegen_groups_chunks(&groups, max_chunks, cuts);

    if (job->fn == CEGEN_FN_RANK)
    {
        const double *x = src->data.dbl;
        double *out = result->num;
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if (threads > 1)
        for (size_t c = 0; c < nchunks; c++)
            for (size_t g = cuts[c]; g < cuts[c + 1]; g++)
                rank_group(job, x, order, outer[g], outer[g + 1] - outer[g], out, miss);
        goto done;
    }

    /* mode(): the distinct values of each group */
    runs = (size_t *)ctools_safe_malloc2(n + 1, sizeof(size_t));
    first_run = (size_t *)ctools_safe_malloc2(groups.ngroups + 1, sizeof(size_t));
    if (!runs || !first_run)
    {
        rc = STATA_ERR_MEMORY;
        goto done;
    }
    const size_t nruns = ctools_group_starts_parallel(data, order, keys, nkeys, runs);
    /* Every group start is a run start: walk the two lists together */
    for (size_t g = 0, r = 0; g < groups.ngroups; g++)
    {
        while (runs[r] < outer[g])
            r++;
        first_run[g] = r;
    }
    first_run[groups.ngroups] = nruns;

    const mode_view view = {src, order, runs, miss};
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        mode_tally tally = {0, 0, 0, 0, 0, 0};
#pragma omp for schedule(dynamic, 1)
        for (size_t c = 0; c < nchunks; c++)
        {
            for (size_t g = cuts[c]; g < cuts[c + 1]; g++)
            {
                const size_t chosen =
                    mode_group(job, &view, first_run[g], first_run[g + 1], &tally);
                if (is_string)
                {
                    char *text = NULL;
                    if (chosen != SIZE_MAX)
                    {
                        text = src->data.str[mode_run_row(&view, chosen)];
                        if (text && text[0] == '\0')
                            text = NULL;
                    }
                    if (text)
                    {
                        const size_t len = strlen(text);
                        if (len > tally.strlen_max)
                            tally.strlen_max = len;
                    }
                    for (size_t pos = outer[g]; pos < outer[g + 1]; pos++)
                        result->str[order ? order[pos] : pos] = text;
                }
                else
                {
                    const double value =
                        chosen != SIZE_MAX ? src->data.dbl[mode_run_row(&view, chosen)] : miss;
                    for (size_t pos = outer[g]; pos < outer[g + 1]; pos++)
                        result->num[order ? order[pos] : pos] = value;
                }
            }
        }
#pragma omp critical(cegen_order_tally)
        mode_tally_merge(total, &tally);
    }

done:
    free(keys);
    free(sort_vars);
    free(outer);
    free(runs);
    free(first_run);
    free(cuts);
    return rc;
}

/* ---------------------------------------------------------------------------
   Driver
   --------------------------------------------------------------------------- */

stata_retcode cegen_order_run(const cegen_job *job, stata_data *data, cegen_result *result,
                              cegen_report *report)
{
    const size_t n = data->nobs;
    const double miss = SV_missval;
    if (job->nsrc != 1)
        return STATA_ERR_INVALID_INPUT;
    stata_variable *src = &data->vars[0];
    const int is_string = src->type == STATA_TYPE_STRING;
    if (job->fn == CEGEN_FN_RANK && is_string)
        return STATA_ERR_UNSUPPORTED_TYPE;

    if (is_string)
    {
        result->str = (char **)ctools_safe_calloc2(n ? n : 1, sizeof(char *));
        if (!result->str)
            return STATA_ERR_MEMORY;
    }
    else
    {
        result->num =
            (double *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, n ? n : 1, sizeof(double));
        if (!result->num)
            return STATA_ERR_MEMORY;
    }
    if (n == 0)
        return STATA_OK;

    /* field ranks the negated values, highest first */
    if (job->fn == CEGEN_FN_RANK && job->field)
    {
        double *x = src->data.dbl;
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n))
        for (size_t i = 0; i < n; i++)
            x[i] = x[i] < miss ? -x[i] : miss;
    }

    mode_tally total = {0, 0, 0, 0, 0, 0};
    stata_retcode rc = STATA_OK;
    int done = 0;
    if (job->nby && !is_string)
    {
        cegen_groups groups;
        rc = cegen_groups_build(data, job->nsrc, job->nby, CEGEN_GROUPS_RUNS, &groups);
        if (rc == STATA_OK && groups.max_count <= CEGEN_ORDER_GROUP_MAX)
        {
            report->ngroups = groups.ngroups;
            rc = order_by_groups(job, data, &groups, result->num, &total);
            done = 1;
        }
        cegen_groups_free(&groups);
    }
    if (rc == STATA_OK && !done)
        rc = order_by_sort(job, data, result, report, &total);
    if (rc != STATA_OK || job->fn != CEGEN_FN_MODE)
        return rc;

    report->group_missing = total.group_missing;
    report->multiple = total.multiple;
    report->max_modes = total.max_modes;
    report->bad_groups = total.bad_groups;
    report->strlen_max = total.strlen_max;
    if (!job->missing && !total.any_value)
    {
        /* Every value is missing: egen warns and generates missing values */
        report->all_missing = 1;
        if (is_string)
            memset(result->str, 0, n * sizeof(char *));
        else
            for (size_t i = 0; i < n; i++)
                result->num[i] = miss;
    }
    return rc;
}

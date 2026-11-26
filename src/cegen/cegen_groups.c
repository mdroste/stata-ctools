/*
    cegen_groups.c
    By-groups of the loaded sample

    The groups are Stata's by-groups: runs of equal key values, each missing
    code its own group, strings compared by their bytes. Three ways to find
    them, cheapest first:

    - data already in key order (the usual panel layout): one parallel pass
      over neighbouring rows, no permutation at all;
    - numeric keys in any order: the shared parallel remap of the HDFE
      commands numbers the groups without sorting (dense integer codes by a
      direct table, anything else by hash partitions), and one stable
      counting pass turns the numbers into runs when a function needs them;
    - string keys, or groups that must come in key order: the stable sorts
      csort uses.

    In every case the rows of a group come in observation order, so sums
    accumulate in the same order on every run and for any thread count.

    License: MIT
*/

#include <stdlib.h>
#include <string.h>

#include "ctools_alloc.h"
#include "ctools_config.h"
#include "ctools_order.h"
#include "hdfe/hdfe_remap.h"
#include "cegen_internal.h"

int cegen_threads_for(size_t n)
{
    int threads = 1;
#ifdef _OPENMP
    if (n >= MIN_OBS_PER_THREAD * 2)
        threads = ctools_get_max_threads();
#else
    (void)n;
#endif
    return threads < 1 ? 1 : threads;
}

int cegen_rows_sorted(const stata_data *data, const int *keys, size_t nkeys)
{
    const size_t n = data->nobs;
    if (n < 2 || !nkeys)
        return 1;
    int unsorted = 0;
    /* One numeric key is the common case and vectorizes as a plain scan */
    if (nkeys == 1 && data->vars[keys[0]].type == STATA_TYPE_DOUBLE)
    {
        const double *x = data->vars[keys[0]].data.dbl;
#pragma omp parallel for schedule(static) reduction(| : unsorted) num_threads(cegen_threads_for(n))
        for (size_t i = 1; i < n; i++)
            unsorted |= (x[i - 1] > x[i]);
        return !unsorted;
    }
    /* Blocks let a worker stop comparing once any block has found disorder */
    const size_t block = 8192, nblocks = (n - 1 + block - 1) / block;
#pragma omp parallel for schedule(dynamic, 4) reduction(| : unsorted) num_threads(cegen_threads_for(n))
    for (size_t b = 0; b < nblocks; b++)
    {
        if (unsorted)
            continue;
        const size_t first = 1 + b * block, stop = first + block < n ? first + block : n;
        for (size_t i = first; i < stop; i++)
        {
            if (ctools_compare_rows(data, i - 1, keys, data, i, keys, nkeys) > 0)
            {
                unsorted = 1;
                break;
            }
        }
    }
    return !unsorted;
}

void cegen_groups_free(cegen_groups *groups)
{
    if (!groups)
        return;
    free(groups->starts);
    free(groups->owned_order);
    free(groups->ids);
    memset(groups, 0, sizeof(*groups));
}

/*
    Number the groups of numeric keys without sorting: ids[row] in
    1 .. *ngroups. One key is one remap; further keys are folded in one at a
    time, the pair (group so far, level of the next key) being remapped as
    one exact integer. Returns STATA_ERR_UNSUPPORTED_TYPE when that integer
    would not be exact (the caller sorts instead).
*/
static stata_retcode groups_number(const stata_data *data, size_t first_key, size_t nkeys,
                                   ST_int **ids_out, size_t *ngroups)
{
    const size_t n = data->nobs;
    const ST_int rows = (ST_int)n;
    ST_int *ids = (ST_int *)ctools_safe_malloc2(n, sizeof(ST_int));
    ST_int *next = NULL;
    double *pair = NULL;
    ST_int levels = 0;
    stata_retcode rc = STATA_OK;

    if (!ids)
        return STATA_ERR_MEMORY;
    if (hdfe_remap_levels(data->vars[first_key].data.dbl, rows, ids, &levels) != 0)
        rc = STATA_ERR_MEMORY;
    for (size_t k = 1; k < nkeys && rc == STATA_OK; k++)
    {
        ST_int levels_k = 0;
        if (!next)
        {
            next = (ST_int *)ctools_safe_malloc2(n, sizeof(ST_int));
            pair = (double *)ctools_safe_malloc2(n, sizeof(double));
            if (!next || !pair)
            {
                rc = STATA_ERR_MEMORY;
                break;
            }
        }
        if (hdfe_remap_levels(data->vars[first_key + k].data.dbl, rows, next, &levels_k) != 0)
        {
            rc = STATA_ERR_MEMORY;
            break;
        }
        if ((double)levels * (double)levels_k > 0x1p53)
        {
            rc = STATA_ERR_UNSUPPORTED_TYPE;
            break;
        }
        const double width = (double)levels_k;
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n))
        for (size_t i = 0; i < n; i++)
            pair[i] = (double)(ids[i] - 1) * width + (double)(next[i] - 1);
        if (hdfe_remap_levels(pair, rows, ids, &levels) != 0)
            rc = STATA_ERR_MEMORY;
    }
    free(next);
    free(pair);
    if (rc != STATA_OK)
    {
        free(ids);
        return rc;
    }
    *ids_out = ids;
    *ngroups = (size_t)levels;
    return STATA_OK;
}

/* Runs from group numbers: a stable counting sort of the rows */
static stata_retcode groups_runs_from_ids(const ST_int *ids, size_t n, size_t ngroups,
                                          cegen_groups *groups)
{
    size_t *starts = (size_t *)ctools_safe_calloc2(ngroups + 1, sizeof(size_t));
    size_t *cursor = (size_t *)ctools_safe_malloc2(ngroups ? ngroups : 1, sizeof(size_t));
    perm_idx_t *order = (perm_idx_t *)ctools_safe_malloc2(n ? n : 1, sizeof(perm_idx_t));
    if (!starts || !cursor || !order)
    {
        free(starts);
        free(cursor);
        free(order);
        return STATA_ERR_MEMORY;
    }
    for (size_t i = 0; i < n; i++)
        starts[ids[i]]++;
    size_t max_count = 0;
    for (size_t g = 0; g < ngroups; g++)
    {
        if (starts[g + 1] > max_count)
            max_count = starts[g + 1];
        starts[g + 1] += starts[g];
        cursor[g] = starts[g];
    }
    /* rows in ascending order within each group */
    for (size_t i = 0; i < n; i++)
        order[cursor[ids[i] - 1]++] = (perm_idx_t)i;
    free(cursor);
    groups->starts = starts;
    groups->owned_order = order;
    groups->order = order;
    groups->max_count = max_count;
    return STATA_OK;
}

stata_retcode cegen_groups_build(stata_data *data, size_t first_key, size_t nkeys, int flags,
                                 cegen_groups *groups)
{
    const size_t n = data->nobs;
    memset(groups, 0, sizeof(*groups));

    if (!nkeys || n == 0)
    {
        groups->starts = (size_t *)malloc(2 * sizeof(size_t));
        if (!groups->starts)
            return STATA_ERR_MEMORY;
        groups->ngroups = n ? 1 : 0;
        groups->starts[0] = 0;
        groups->starts[1] = n;
        groups->max_count = n;
        return STATA_OK;
    }

    int *keys = (int *)malloc(nkeys * sizeof(int));
    int *sort_vars = (int *)malloc(nkeys * sizeof(int));
    size_t *starts = (size_t *)ctools_safe_malloc2(n + 1, sizeof(size_t));
    if (!keys || !sort_vars || !starts)
    {
        free(keys);
        free(sort_vars);
        free(starts);
        return STATA_ERR_MEMORY;
    }
    for (size_t k = 0; k < nkeys; k++)
    {
        keys[k] = (int)(first_key + k);          /* zero-based column */
        sort_vars[k] = (int)(first_key + k) + 1; /* one-based sort key */
    }

    stata_retcode rc = STATA_OK;
    const perm_idx_t *order = NULL;
    const int sorted = cegen_rows_sorted(data, keys, nkeys);
    if (!sorted && !(flags & CEGEN_GROUPS_KEY_ORDER))
    {
        /* Numeric keys: number the groups instead of sorting the rows */
        int numeric = 1;
        for (size_t k = 0; k < nkeys; k++)
            numeric &= data->vars[first_key + k].type == STATA_TYPE_DOUBLE;
        if (numeric)
        {
            rc = groups_number(data, first_key, nkeys, &groups->ids, &groups->ngroups);
            if (rc == STATA_OK && (flags & CEGEN_GROUPS_RUNS))
                rc = groups_runs_from_ids(groups->ids, n, groups->ngroups, groups);
            if (rc != STATA_ERR_UNSUPPORTED_TYPE)
            {
                free(keys);
                free(sort_vars);
                free(starts);
                if (rc != STATA_OK)
                    cegen_groups_free(groups);
                return rc;
            }
            rc = STATA_OK; /* too many key combinations to number: sort */
        }
    }
    if (!sorted)
    {
        if (!data->sort_order)
            rc = STATA_ERR_INVALID_INPUT;
        else
            rc = ctools_sort_dispatch(data, sort_vars, nkeys, SORT_ALG_AUTO);
        order = data->sort_order;
    }
    if (rc == STATA_OK)
    {
        size_t ngroups = ctools_group_starts_parallel(data, order, keys, nkeys, starts);
        size_t max_count = 0;
        for (size_t g = 0; g < ngroups; g++)
        {
            size_t count = starts[g + 1] - starts[g];
            if (count > max_count)
                max_count = count;
        }
        size_t *trimmed = (size_t *)realloc(starts, (ngroups + 1) * sizeof(size_t));
        groups->starts = trimmed ? trimmed : starts;
        starts = NULL;
        groups->ngroups = ngroups;
        groups->order = order;
        groups->max_count = max_count;
    }
    free(keys);
    free(sort_vars);
    free(starts);
    return rc;
}

size_t cegen_groups_chunks(const cegen_groups *groups, size_t max_chunks, size_t *cuts)
{
    const size_t ngroups = groups->ngroups;
    if (!ngroups || !max_chunks)
    {
        cuts[0] = 0;
        return 0;
    }
    const size_t n = groups->starts[ngroups];
    size_t want = max_chunks < ngroups ? max_chunks : ngroups;
    size_t nchunks = 0;
    cuts[0] = 0;
    for (size_t c = 1; c < want; c++)
    {
        /* the first group that starts at or after row c * n / want */
        const size_t target = (size_t)(((unsigned long long)n * c) / want);
        size_t lo = cuts[nchunks] + 1, hi = ngroups;
        while (lo < hi)
        {
            size_t mid = lo + (hi - lo) / 2;
            if (groups->starts[mid] < target)
                lo = mid + 1;
            else
                hi = mid;
        }
        if (lo >= ngroups)
            break;
        cuts[++nchunks] = lo;
    }
    cuts[++nchunks] = ngroups;
    return nchunks;
}

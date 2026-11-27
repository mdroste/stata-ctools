/*
 * ccollapse_quantile.c
 *
 * Percentiles (p1 ... p99, median) and the interquartile range.
 *
 * The definition is collapse's, which is that of summarize, detail and
 * _pctile: with the nonmissing values in ascending order, W their total
 * weight (their number when unweighted) and P = W p / 100, percentile p is
 * the first value whose cumulative weight exceeds P, or the mean of that
 * value and the next one when a cumulative weight equals P. iqr is p75 - p25.
 *
 * Groups are independent: threads take groups, gather a group's nonmissing
 * values into a buffer of their own and select from it. Unweighted values
 * need no sort, only the order statistics the percentiles ask for
 * (quickselect, narrowed from one requested rank to the next). Weighted
 * values are sorted with their weights, equal values by weight as collapse
 * orders them (the order sets the running total when weights are negative
 * or too small to change it).
 *
 * A group of a million rows or more would keep one thread busy while the
 * others wait, so it is selected from with all of them: a sorted sample of
 * the group brackets each requested rank between two values, one parallel
 * pass counts the values below each bracket and collects the few inside it,
 * and the order statistics are selected from those. The result is the exact
 * order statistic; a bracket that misses its rank (a sample that happens to
 * misrepresent the group) sends the group to the ordinary path.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ccollapse_internal.h"
#include "ctools_cancel.h"
#include "ctools_config.h"
#include "ctools_select.h"

typedef struct
{
    double x, w;
} cc_pair;

/* The percentile targets of one source */
typedef struct
{
    size_t source;
    int pct[99];      /* distinct percentiles, ascending */
    int npct;
    double value[99]; /* scratch: the value of pct[k] in the current group */
} cc_quantile_source;

static int cc_pct_slot(const cc_quantile_source *qs, int pct)
{
    for (int k = 0; k < qs->npct; k++)
        if (qs->pct[k] == pct)
            return k;
    return -1;
}

static void cc_pct_add(cc_quantile_source *qs, int pct)
{
    int k = 0;
    while (k < qs->npct && qs->pct[k] < pct)
        k++;
    if (k < qs->npct && qs->pct[k] == pct)
        return;
    memmove(&qs->pct[k + 1], &qs->pct[k], (size_t)(qs->npct - k) * sizeof(qs->pct[0]));
    qs->pct[k] = pct;
    qs->npct++;
}

/* Put the order statistics of a[lo, hi) at ranks[0 .. nranks) (ascending,
 * distinct, within [lo, hi)) in place, with smaller values before and larger
 * values after each of them */
static void cc_multiselect(double *a, size_t lo, size_t hi, const size_t *ranks, size_t nranks)
{
    while (nranks > 0 && hi - lo > 1)
    {
        size_t mid = nranks / 2, k = ranks[mid];
        ctools_quickselect_double(a + lo, hi - lo, k - lo);
        /* recurse into the side with fewer ranks, continue with the other */
        size_t nleft = mid, nright = nranks - mid - 1;
        if (nleft <= nright)
        {
            cc_multiselect(a, lo, k, ranks, nleft);
            lo = k + 1;
            ranks += mid + 1;
            nranks = nright;
        }
        else
        {
            cc_multiselect(a, k + 1, hi, ranks + mid + 1, nright);
            hi = k;
            nranks = nleft;
        }
    }
}

/* Unweighted percentiles of the n values in a (reordered) */
static void cc_percentiles(double *a, size_t n, cc_quantile_source *qs)
{
    size_t ranks[2 * 99];
    size_t nranks = 0;
    for (int k = 0; k < qs->npct; k++)
    {
        double P = (double)n * (double)qs->pct[k] / 100.0;
        size_t i = (size_t)floor(P);
        if (P == (double)i && i > 0 && (nranks == 0 || ranks[nranks - 1] != i - 1))
            ranks[nranks++] = i - 1;
        if (i < n && (nranks == 0 || ranks[nranks - 1] != i))
            ranks[nranks++] = i;
    }
    cc_multiselect(a, 0, n, ranks, nranks);
    for (int k = 0; k < qs->npct; k++)
    {
        double P = (double)n * (double)qs->pct[k] / 100.0;
        size_t i = (size_t)floor(P);
        if (i >= n)
            qs->value[k] = a[n - 1];
        else if (P == (double)i && i > 0)
            qs->value[k] = (a[i - 1] + a[i]) / 2.0;
        else
            qs->value[k] = a[i];
    }
}

/* collapse's order of the pairs: by value, and equal values by weight */
static inline int cc_pair_before(const cc_pair *a, const cc_pair *b)
{
    return a->x < b->x || (a->x == b->x && a->w < b->w);
}

/* Stable merge sort of pairs by value, then weight; tmp holds n pairs */
static void cc_sort_pairs(cc_pair *a, cc_pair *tmp, size_t n)
{
    enum { CC_RUN = 16 };
    for (size_t lo = 0; lo < n; lo += CC_RUN)
    {
        size_t hi = lo + CC_RUN < n ? lo + CC_RUN : n;
        for (size_t i = lo + 1; i < hi; i++)
        {
            cc_pair key = a[i];
            size_t j = i;
            while (j > lo && cc_pair_before(&key, &a[j - 1]))
            {
                a[j] = a[j - 1];
                j--;
            }
            a[j] = key;
        }
    }
    cc_pair *src = a, *dst = tmp;
    for (size_t width = CC_RUN; width < n; width *= 2)
    {
        for (size_t lo = 0; lo < n; lo += 2 * width)
        {
            size_t mid = lo + width < n ? lo + width : n;
            size_t hi = lo + 2 * width < n ? lo + 2 * width : n;
            size_t l = lo, r = mid, out = lo;
            while (l < mid && r < hi)
                dst[out++] = cc_pair_before(&src[r], &src[l]) ? src[r++] : src[l++];
            while (l < mid)
                dst[out++] = src[l++];
            while (r < hi)
                dst[out++] = src[r++];
        }
        cc_pair *swap = src;
        src = dst;
        dst = swap;
    }
    if (src != a)
        memcpy(a, src, n * sizeof(*a));
}

/*
    Weighted percentiles of the n pairs in a (sorted here). The running
    total reaches P = W p / 100 at the first pair whose cumulative weight is
    P or more. When it equals P, the percentile is the mean of that value and
    the next. Weights too small to change the total can leave it at P over
    several pairs; collapse then averages the last of them with the next
    when every weight is positive, and the first of them when a weight is
    negative (iweights).

    A total of zero or less (iweights only) is collapse's rescaled case:
    W = 0 gives missing (w * n / W is missing in Stata), and with W negative
    each weight becomes w * n / W, as in collapse's ado implementation
    (collapse_16's remakew), which turns the running total around. Equal
    values are then taken in the order of their rescaled weights (by
    descending weight), percentile p is the first value whose rescaled
    running total reaches G = (rescaled total) p / 100 to within 1e-6, and
    the mean of it and the next value when that total is within 1e-6 of G.
    This reproduced collapse in every such cell probed (32,000 cells, exact
    and inexact weights). W is the total in ascending order.
*/
static void cc_weighted_percentiles_rescaled(cc_pair *a, size_t n, double W,
                                             cc_quantile_source *qs)
{
    const double missval = SV_missval;
    if (W == 0.0)
    {
        for (int k = 0; k < qs->npct; k++)
            qs->value[k] = missval;
        return;
    }
    if (W < 0.0)
    {
        /* each run of equal values, ascending by weight, turned around */
        for (size_t lo = 0; lo < n;)
        {
            size_t hi = lo + 1;
            while (hi < n && a[hi].x == a[lo].x)
                hi++;
            for (size_t i = lo, j = hi - 1; i < j; i++, j--)
            {
                cc_pair t = a[i];
                a[i] = a[j];
                a[j] = t;
            }
            lo = hi;
        }
    }
    const double m = (double)n;
    double total = 0.0;
    for (size_t i = 0; i < n; i++)
        total += a[i].w * m / W;
    for (int k = 0; k < qs->npct; k++)
    {
        const double G = total * (double)qs->pct[k] / 100.0;
        double cumulative = 0.0;
        size_t i = 0;
        for (; i < n; i++)
        {
            cumulative += a[i].w * m / W;
            if (cumulative >= G - 1e-6)
                break;
        }
        if (i == n)
            i = n - 1;
        qs->value[k] = fabs(cumulative - G) <= 1e-6 && i + 1 < n ? (a[i].x + a[i + 1].x) / 2.0
                                                                 : a[i].x;
    }
}

static void cc_weighted_percentiles(cc_pair *a, cc_pair *tmp, size_t n, cc_quantile_source *qs)
{
    cc_sort_pairs(a, tmp, n);
    double W = 0.0;
    int negative = 0;
    for (size_t i = 0; i < n; i++)
    {
        W += a[i].w;
        negative |= a[i].w < 0;
    }
    if (negative && !(W > 0.0))
    {
        cc_weighted_percentiles_rescaled(a, n, W, qs);
        return;
    }
    size_t i = 0;
    double cumulative = a[0].w;
    for (int k = 0; k < qs->npct; k++)
    {
        double P = W * (double)qs->pct[k] / 100.0;
        while (i + 1 < n && cumulative < P)
            cumulative += a[++i].w;
        if (cumulative == P && !negative)
            while (i + 1 < n && cumulative + a[i + 1].w == P)
                cumulative += a[++i].w;
        qs->value[k] = cumulative == P && i + 1 < n ? (a[i].x + a[i + 1].x) / 2.0 : a[i].x;
    }
}

/* The targets of qs in group g, from qs->value; n nonmissing values */
static void cc_assign(const cc_job *job, const cc_quantile_source *qs, size_t g, size_t n,
                      double missval)
{
    for (size_t t = 0; t < job->ntargets; t++)
    {
        const cc_target *target = &job->targets[t];
        if ((size_t)target->source != qs->source)
            continue;
        if (target->stat == CC_STAT_PCTILE)
            job->values[t][g] = n ? qs->value[cc_pct_slot(qs, target->pct)] : missval;
        else if (target->stat == CC_STAT_IQR)
        {
            /* a missing quartile (rescaled iweights) makes the range missing */
            double hi = qs->value[cc_pct_slot(qs, 75)], lo = qs->value[cc_pct_slot(qs, 25)];
            job->values[t][g] =
                n && hi < missval && lo < missval ? cc_storable(hi - lo, missval) : missval;
        }
    }
}

/* One source in one group, on the calling thread. buffer holds the group's
 * values, or twice its pairs when weighted: the second half, half pairs
 * away, is the sorter's scratch. */
static void cc_group_quantiles(const cc_job *job, cc_quantile_source *qs, size_t g,
                               void *buffer, size_t half, double missval)
{
    const cc_groups *groups = job->groups;
    const perm_idx_t *perm = groups->perm;
    const double *w = job->weight;
    const double *x = job->data->vars[job->nby + qs->source].data.dbl;
    size_t lo = groups->starts[g], hi = groups->starts[g + 1], n = 0;
    if (!w)
    {
        double *values = buffer;
        if (perm)
        {
            for (size_t p = lo; p < hi; p++)
            {
                double v = x[perm[p]];
                if (v < missval)
                    values[n++] = v;
            }
        }
        else
        {
            for (size_t p = lo; p < hi; p++)
            {
                double v = x[p];
                if (v < missval)
                    values[n++] = v;
            }
        }
        if (n)
            cc_percentiles(values, n, qs);
    }
    else
    {
        cc_pair *pairs = buffer;
        for (size_t p = lo; p < hi; p++)
        {
            size_t row = perm ? perm[p] : p;
            double v = x[row];
            if (v < missval)
            {
                pairs[n].x = v;
                pairs[n].w = w[row];
                n++;
            }
        }
        if (n)
            cc_weighted_percentiles(pairs, pairs + half, n, qs);
    }
    cc_assign(job, qs, g, n, missval);
}

/* ===========================================================================
   Large groups: selection with every thread
   =========================================================================== */

#define CC_LARGE_GROUP ((size_t)1 << 20)
#define CC_SAMPLE_SIZE ((size_t)1 << 16)

static int cc_compare_doubles(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* The region of v among nint ascending, disjoint intervals [a[k], b[k]]:
 * 2k + 1 inside interval k, 2k in the gap below it, 2 nint above them all */
static inline size_t cc_region(double v, const double *a, const double *b, size_t nint)
{
    size_t lo = 0, hi = nint;
    while (lo < hi)
    {
        size_t mid = lo + (hi - lo) / 2;
        if (v > b[mid])
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo < nint && v >= a[lo] ? 2 * lo + 1 : 2 * lo;
}

/*
    Unweighted percentiles of the large group at positions [lo, hi), with
    team threads. Returns 1 with qs->value set and *n_out the number of
    nonmissing values, 0 when the group must take the ordinary path, -1 when
    memory ran out.
*/
static int cc_percentiles_sampled(const double *x, const perm_idx_t *perm, size_t lo, size_t hi,
                                  double missval, cc_quantile_source *qs, int team, size_t *n_out)
{
    enum { CC_MAX_RANKS = 2 * 99 };
    size_t ranks[CC_MAX_RANKS], nranks = 0;
    double rank_value[CC_MAX_RANKS];
    double a[CC_MAX_RANKS], b[CC_MAX_RANKS];
    size_t of_interval[CC_MAX_RANKS];
    size_t below_interval[CC_MAX_RANKS], inside_start[CC_MAX_RANKS + 1];
    const size_t len = hi - lo;
    const size_t nregions = 2 * CC_MAX_RANKS + 1;
    double *sample = NULL, *inside = NULL;
    size_t *counts = NULL;
    int result = 0;

    size_t n = 0;
    #pragma omp parallel for schedule(static) num_threads(team) reduction(+: n)
    for (size_t p = lo; p < hi; p++)
        n += x[perm ? perm[p] : p] < missval;
    *n_out = n;
    if (n < CC_LARGE_GROUP / 2)
        return 0;

    for (int k = 0; k < qs->npct; k++)
    {
        double P = (double)n * (double)qs->pct[k] / 100.0;
        size_t i = (size_t)floor(P);
        if (P == (double)i && i > 0 && (nranks == 0 || ranks[nranks - 1] != i - 1))
            ranks[nranks++] = i - 1;
        if (i < n && (nranks == 0 || ranks[nranks - 1] != i))
            ranks[nranks++] = i;
    }

    sample = malloc(CC_SAMPLE_SIZE * sizeof(*sample));
    counts = calloc((size_t)team * nregions, sizeof(*counts));
    if (!sample || !counts)
    {
        result = -1;
        goto done;
    }
    /* Evenly spaced rows: representative of ordered and of shuffled data */
    size_t stride = len / CC_SAMPLE_SIZE, ns = 0;
    for (size_t j = 0; j < CC_SAMPLE_SIZE; j++)
    {
        size_t p = lo + j * stride + stride / 2;
        double v = x[perm ? perm[p] : p];
        if (v < missval)
            sample[ns++] = v;
    }
    if (ns < CC_SAMPLE_SIZE / 2)
        goto done;
    qsort(sample, ns, sizeof(*sample), cc_compare_doubles);

    /* A bracket of sample values around each rank, four standard errors of
     * its place in the sample wide; overlapping brackets become one interval */
    const size_t margin = 2 * (size_t)sqrt((double)ns) + 16;
    size_t nint = 0;
    for (size_t r = 0; r < nranks; r++)
    {
        size_t at = (size_t)((double)ranks[r] / (double)n * (double)ns);
        double below = at > margin ? sample[at - margin] : -HUGE_VAL;
        double above = at + margin < ns ? sample[at + margin] : HUGE_VAL;
        if (nint && below <= b[nint - 1])
        {
            if (above > b[nint - 1])
                b[nint - 1] = above;
        }
        else
        {
            a[nint] = below;
            b[nint] = above;
            nint++;
        }
        of_interval[r] = nint - 1;
    }

    /* Pass 1: every share of the rows counts its values by region */
    #pragma omp parallel for schedule(static, 1) num_threads(team)
    for (int t = 0; t < team; t++)
    {
        size_t *mine = counts + (size_t)t * nregions;
        size_t from = lo + len * (size_t)t / (size_t)team;
        size_t to = lo + len * (size_t)(t + 1) / (size_t)team;
        for (size_t p = from; p < to; p++)
        {
            double v = x[perm ? perm[p] : p];
            if (v < missval)
                mine[cc_region(v, a, b, nint)]++;
        }
    }
    /* The values below each interval, the values inside it, and where each
     * share writes its part of them */
    size_t running = 0, kept = 0;
    for (size_t region = 0; region < 2 * nint + 1; region++)
    {
        size_t total = 0;
        for (int t = 0; t < team; t++)
        {
            size_t count = counts[(size_t)t * nregions + region];
            if (region & 1)
                counts[(size_t)t * nregions + region] = kept + total;
            total += count;
        }
        if (region & 1)
        {
            below_interval[region / 2] = running;
            inside_start[region / 2] = kept;
            kept += total;
        }
        running += total;
    }
    inside_start[nint] = kept;
    for (size_t r = 0; r < nranks; r++)
    {
        size_t k = of_interval[r];
        if (ranks[r] < below_interval[k] ||
            ranks[r] - below_interval[k] >= inside_start[k + 1] - inside_start[k])
            goto done; /* the bracket missed its rank */
    }
    inside = malloc((kept ? kept : 1) * sizeof(*inside));
    if (!inside)
    {
        result = -1;
        goto done;
    }
    /* Pass 2: collect the values inside the intervals */
    #pragma omp parallel for schedule(static, 1) num_threads(team)
    for (int t = 0; t < team; t++)
    {
        size_t *mine = counts + (size_t)t * nregions;
        size_t from = lo + len * (size_t)t / (size_t)team;
        size_t to = lo + len * (size_t)(t + 1) / (size_t)team;
        for (size_t p = from; p < to; p++)
        {
            double v = x[perm ? perm[p] : p];
            if (v < missval)
            {
                size_t region = cc_region(v, a, b, nint);
                if (region & 1)
                    inside[mine[region]++] = v;
            }
        }
    }
    for (size_t r = 0; r < nranks;)
    {
        size_t k = of_interval[r], first = r, local[CC_MAX_RANKS];
        while (r < nranks && of_interval[r] == k)
        {
            local[r - first] = inside_start[k] + (ranks[r] - below_interval[k]);
            r++;
        }
        cc_multiselect(inside, inside_start[k], inside_start[k + 1], local, r - first);
        for (size_t q = first; q < r; q++)
            rank_value[q] = inside[local[q - first]];
    }
    for (int k = 0; k < qs->npct; k++)
    {
        double P = (double)n * (double)qs->pct[k] / 100.0;
        size_t i = (size_t)floor(P);
        size_t want = i < n ? i : n - 1, at = 0;
        while (at + 1 < nranks && ranks[at] < want)
            at++;
        if (i < n && P == (double)i && i > 0)
            qs->value[k] = (rank_value[at - 1] + rank_value[at]) / 2.0;
        else
            qs->value[k] = rank_value[at];
    }
    result = 1;

done:
    free(sample);
    free(counts);
    free(inside);
    return result;
}

/* ===========================================================================
   Entry point
   =========================================================================== */

stata_retcode cc_compute_quantiles(const cc_job *job)
{
    const cc_groups *groups = job->groups;
    const size_t ngroups = groups->ngroups;
    const size_t *starts = groups->starts;
    const perm_idx_t *perm = groups->perm;
    const double *w = job->weight;
    const double missval = SV_missval;

    size_t nq = 0;
    for (size_t t = 0; t < job->ntargets; t++)
        nq += cc_stat_is_quantile(job->targets[t].stat);
    if (!nq)
        return STATA_OK;
    if (!starts || (!perm && !groups->contiguous))
        return STATA_ERR_INVALID_INPUT;

    /* The sources that have percentile targets, with their percentiles */
    cc_quantile_source *proto = calloc(job->nsources, sizeof(*proto));
    if (!proto)
        return STATA_ERR_MEMORY;
    size_t nqs = 0;
    for (size_t s = 0; s < job->nsources; s++)
    {
        cc_quantile_source *qs = &proto[nqs];
        qs->source = s;
        qs->npct = 0;
        for (size_t t = 0; t < job->ntargets; t++)
        {
            const cc_target *target = &job->targets[t];
            if ((size_t)target->source != s)
                continue;
            if (target->stat == CC_STAT_PCTILE)
                cc_pct_add(qs, target->pct);
            else if (target->stat == CC_STAT_IQR)
            {
                cc_pct_add(qs, 25);
                cc_pct_add(qs, 75);
            }
        }
        if (qs->npct)
            nqs++;
    }

    const size_t per_value = w ? 2 * sizeof(cc_pair) : sizeof(double);
    int failure = STATA_OK;

    /* Large groups first, each with every thread (unweighted values only) */
    const int sampled = !w && job->threads > 1;
    size_t largest = 0;
    for (size_t g = 0; g < ngroups && failure == STATA_OK; g++)
    {
        size_t rows = starts[g + 1] - starts[g];
        if (!sampled || rows < CC_LARGE_GROUP)
        {
            if (rows > largest)
                largest = rows;
            continue;
        }
        if (ctools_cancel_requested())
        {
            failure = (int)STATA_ERR_CANCELLED;
            break;
        }
        void *buffer = NULL;
        for (size_t q = 0; q < nqs && failure == STATA_OK; q++)
        {
            cc_quantile_source *qs = &proto[q];
            const double *x = job->data->vars[job->nby + qs->source].data.dbl;
            size_t n = 0;
            int done = cc_percentiles_sampled(x, perm, starts[g], starts[g + 1], missval, qs,
                                              job->threads, &n);
            if (done > 0)
                cc_assign(job, qs, g, n, missval);
            else if (done < 0)
                failure = (int)STATA_ERR_MEMORY;
            else
            {
                if (!buffer)
                    buffer = malloc(rows * per_value);
                if (!buffer)
                    failure = (int)STATA_ERR_MEMORY;
                else
                    cc_group_quantiles(job, qs, g, buffer, rows, missval);
            }
        }
        free(buffer);
    }
    if (failure != STATA_OK)
    {
        free(proto);
        return (stata_retcode)failure;
    }

    if (largest == 0)
        largest = 1;
    if (largest > SIZE_MAX / per_value)
    {
        free(proto);
        return STATA_ERR_MEMORY;
    }
    int team = ctools_workspace_threads(job->threads, ngroups, largest * per_value);
    size_t chunk = ngroups / ((size_t)team * 16);
    if (chunk < 1)
        chunk = 1;
    if (chunk > 256)
        chunk = 256;

    #pragma omp parallel num_threads(team) if (team > 1)
    {
        /* each thread: its own value buffer and its own copy of the lists */
        void *buffer = malloc(largest * per_value);
        cc_quantile_source *mine = malloc((nqs ? nqs : 1) * sizeof(*mine));
        if (!buffer || !mine)
        {
            #pragma omp atomic write
            failure = (int)STATA_ERR_MEMORY;
        }
        else
            memcpy(mine, proto, nqs * sizeof(*mine));

        #pragma omp for schedule(dynamic, chunk)
        for (size_t g = 0; g < ngroups; g++)
        {
            int seen;
            #pragma omp atomic read
            seen = failure;
            if (seen != STATA_OK)
                continue;
            if (sampled && starts[g + 1] - starts[g] >= CC_LARGE_GROUP)
                continue; /* done above */
            if ((g & 1023) == 0 && ctools_cancel_requested())
            {
                #pragma omp atomic write
                failure = (int)STATA_ERR_CANCELLED;
                continue;
            }
            for (size_t q = 0; q < nqs; q++)
                cc_group_quantiles(job, &mine[q], g, buffer, largest, missval);
        }
        free(buffer);
        free(mine);
    }
    free(proto);
    return (stata_retcode)failure;
}

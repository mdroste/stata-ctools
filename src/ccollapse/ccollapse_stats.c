/*
 * ccollapse_stats.c
 *
 * Every collapse statistic that needs no ordering of the values: mean, sum,
 * rawsum, count, percent, sd, semean, sebinomial, sepoisson, min, max,
 * first, last, firstnm and lastnm.
 *
 * Each statistic is a pass (or two) over a source column that adds row i to
 * the accumulators of group gid[i]. A group's rows are visited in their
 * order in the dataset and by one thread, so the sums are the sequential
 * sums collapse forms on data sorted by the by() variables, and they do not
 * depend on the number of threads. Two ways to use the threads keep that:
 *
 *   rows in group order   the rows are cut at group boundaries and every
 *                         thread takes a run of whole groups;
 *   otherwise             every thread takes a source variable.
 *
 * The formulas are collapse's (W: sum of weights, n: nonmissing values):
 *
 *   mean        sum(w x) / W
 *   sum         sum(w x); with aweights sum(w x) n / W
 *   sd          sqrt(sum(w (x - mean)^2) / (W - 1));
 *               with aweights sqrt(sum(w (x - mean)^2) (n / W) / (n - 1))
 *   semean      sd / sqrt(N),  N = n, or W with fweights; collapse forms it
 *               that way when sd of the same variable is also asked for,
 *               and otherwise as sqrt(sum(w (x - mean)^2) / ((W - 1) N))
 *   sebinomial  sqrt(mean (1 - mean) / N), missing unless all values are 0/1
 *   sepoisson   sqrt(int(mean N + .5)) / N, missing if a value is negative
 *   percent     100 W / (W summed over all groups)
 *   count       n, or W with fweights, iweights and pweights
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

/* Least rows for which the rows of ordered groups are shared out */
#define CC_STATS_PARALLEL_MIN ((size_t)262144)

/* What the targets of one source ask for */
typedef struct
{
    int any;      /* the source has a target handled here */
    int moments;  /* n, W, sum(w x) */
    int rawsum;   /* unweighted sum next to a weighted one */
    int spread;   /* sum(w (x - mean)^2) */
    int binary;   /* all values 0 or 1? */
    int negative; /* any value negative? */
    int min, max;
    int last;     /* last row of each group */
    int firstnm, lastnm;
    int percent;
    int sd;       /* sd itself is a target */
} cc_needs;

/* Per-group accumulators of one source; unused ones stay NULL */
typedef struct
{
    uint32_t *n;         /* nonmissing values */
    double *sw;          /* sum of their weights (weighted only) */
    double *swx;         /* sum(w x), or sum(x) */
    double *rs;          /* sum(x) when swx is weighted */
    double *mean;
    double *b;           /* sum(w (x - mean)^2) */
    int with_sd;         /* sd is a target too: semean is formed from it */
    double *mn;          /* least value, missing values included */
    double *mx;          /* greatest nonmissing value */
    double *mxall;       /* greatest value, missing values included */
    unsigned char *flag; /* bit 0: a value other than 0/1; bit 1: a negative value */
    perm_idx_t *last;
    perm_idx_t *fnm, *lnm;
    double total;        /* W over all groups, for percent */
} cc_acc;

static void cc_acc_free(cc_acc *a)
{
    free(a->n);
    free(a->sw);
    free(a->swx);
    free(a->rs);
    free(a->mean);
    free(a->b);
    free(a->mn);
    free(a->mx);
    free(a->mxall);
    free(a->flag);
    free(a->last);
    free(a->fnm);
    free(a->lnm);
    memset(a, 0, sizeof(*a));
}

static cc_needs cc_source_needs(const cc_job *job, size_t source)
{
    cc_needs need;
    memset(&need, 0, sizeof(need));
    int weighted = job->weight != NULL;
    for (size_t t = 0; t < job->ntargets; t++)
    {
        if ((size_t)job->targets[t].source != source)
            continue;
        switch (job->targets[t].stat)
        {
        case CC_STAT_MEAN:
        case CC_STAT_SUM:
        case CC_STAT_COUNT:
            need.moments = 1;
            break;
        case CC_STAT_RAWSUM:
            need.moments = 1;
            need.rawsum = weighted;
            break;
        case CC_STAT_PERCENT:
            need.moments = need.percent = 1;
            break;
        case CC_STAT_SD:
            need.sd = 1;
            need.moments = need.spread = 1;
            break;
        case CC_STAT_SEMEAN:
            need.moments = need.spread = 1;
            break;
        case CC_STAT_SEBINOMIAL:
            need.moments = need.binary = 1;
            break;
        case CC_STAT_SEPOISSON:
            need.moments = need.negative = 1;
            break;
        case CC_STAT_MIN:
            need.min = 1;
            break;
        case CC_STAT_MAX:
            need.max = 1;
            break;
        case CC_STAT_FIRST:
            break;
        case CC_STAT_LAST:
            need.last = 1;
            break;
        case CC_STAT_FIRSTNM:
            need.firstnm = 1;
            break;
        case CC_STAT_LASTNM:
            need.lastnm = need.last = 1;
            break;
        default:
            continue; /* percentiles: ccollapse_quantile.c */
        }
        need.any = 1;
    }
    return need;
}

static stata_retcode cc_acc_alloc(cc_acc *a, const cc_needs *need, size_t ngroups, int weighted)
{
    size_t count = ngroups ? ngroups : 1;
    int ok = 1;
    memset(a, 0, sizeof(*a));
    if (need->moments)
    {
        ok &= (a->n = calloc(count, sizeof(*a->n))) != NULL;
        ok &= (a->swx = calloc(count, sizeof(*a->swx))) != NULL;
        if (weighted)
            ok &= (a->sw = calloc(count, sizeof(*a->sw))) != NULL;
        if (need->rawsum)
            ok &= (a->rs = calloc(count, sizeof(*a->rs))) != NULL;
    }
    a->with_sd = need->sd;
    if (need->spread)
    {
        ok &= (a->mean = malloc(count * sizeof(*a->mean))) != NULL;
        ok &= (a->b = calloc(count, sizeof(*a->b))) != NULL;
    }
    if (need->binary || need->negative)
        ok &= (a->flag = calloc(count, 1)) != NULL;
    if (need->min)
        ok &= (a->mn = malloc(count * sizeof(*a->mn))) != NULL;
    if (need->max)
    {
        ok &= (a->mx = malloc(count * sizeof(*a->mx))) != NULL;
        ok &= (a->mxall = malloc(count * sizeof(*a->mxall))) != NULL;
    }
    if (need->last)
        ok &= (a->last = malloc(count * sizeof(*a->last))) != NULL;
    if (need->firstnm)
        ok &= (a->fnm = malloc(count * sizeof(*a->fnm))) != NULL;
    if (need->lastnm)
        ok &= (a->lnm = malloc(count * sizeof(*a->lnm))) != NULL;
    if (!ok)
    {
        cc_acc_free(a);
        return STATA_ERR_MEMORY;
    }
    for (size_t g = 0; g < ngroups; g++)
    {
        if (a->mn)
            a->mn[g] = HUGE_VAL;
        if (a->mx)
        {
            a->mx[g] = -HUGE_VAL;
            a->mxall[g] = -HUGE_VAL;
        }
        if (a->fnm)
            a->fnm[g] = CC_NONE;
        if (a->lnm)
            a->lnm[g] = CC_NONE;
    }
    return STATA_OK;
}

/* ===========================================================================
   Passes over rows [lo, hi) of a numeric source
   =========================================================================== */

static void cc_pass_moments(const double *x, const double *w, const uint32_t *gid, size_t lo,
                            size_t hi, double missval, cc_acc *a)
{
    uint32_t *n = a->n;
    double *swx = a->swx;
    if (!w)
    {
        for (size_t i = lo; i < hi; i++)
        {
            double v = x[i];
            if (v < missval)
            {
                uint32_t g = gid[i];
                n[g]++;
                swx[g] += v;
            }
        }
        return;
    }
    double *sw = a->sw, *rs = a->rs;
    if (!rs)
    {
        for (size_t i = lo; i < hi; i++)
        {
            double v = x[i];
            if (v < missval)
            {
                uint32_t g = gid[i];
                n[g]++;
                sw[g] += w[i];
                swx[g] += w[i] * v;
            }
        }
        return;
    }
    for (size_t i = lo; i < hi; i++)
    {
        double v = x[i];
        if (v < missval)
        {
            uint32_t g = gid[i];
            n[g]++;
            sw[g] += w[i];
            swx[g] += w[i] * v;
            rs[g] += v;
        }
    }
}

static void cc_pass_spread(const double *x, const double *w, const uint32_t *gid, size_t lo,
                           size_t hi, double missval, cc_acc *a)
{
    const double *mean = a->mean;
    double *b = a->b;
    if (!w)
    {
        for (size_t i = lo; i < hi; i++)
        {
            double v = x[i];
            if (v < missval)
            {
                uint32_t g = gid[i];
                double d = v - mean[g];
                b[g] += d * d;
            }
        }
        return;
    }
    for (size_t i = lo; i < hi; i++)
    {
        double v = x[i];
        if (v < missval)
        {
            uint32_t g = gid[i];
            double d = v - mean[g];
            b[g] += (w[i] * d) * d;
        }
    }
}

static void cc_pass_flags(const double *x, const uint32_t *gid, size_t lo, size_t hi,
                          double missval, cc_acc *a)
{
    unsigned char *flag = a->flag;
    for (size_t i = lo; i < hi; i++)
    {
        double v = x[i];
        if (v < missval)
            flag[gid[i]] |= (unsigned char)((v != 0.0 && v != 1.0) | ((v < 0.0) << 1));
    }
}

static void cc_pass_extremes(const double *x, const uint32_t *gid, size_t lo, size_t hi,
                             double missval, cc_acc *a)
{
    double *mn = a->mn, *mx = a->mx, *mxall = a->mxall;
    if (mn)
    {
        /* missing values sort above every number, so the least value of a
         * group is a number whenever the group has one */
        for (size_t i = lo; i < hi; i++)
        {
            uint32_t g = gid[i];
            if (x[i] < mn[g])
                mn[g] = x[i];
        }
    }
    if (mx)
    {
        for (size_t i = lo; i < hi; i++)
        {
            uint32_t g = gid[i];
            double v = x[i];
            if (v < missval)
            {
                if (v > mx[g])
                    mx[g] = v;
            }
            else if (v > mxall[g])
                mxall[g] = v;
        }
    }
}

static void cc_pass_rows(const stata_variable *var, const uint32_t *gid, size_t lo, size_t hi,
                         double missval, cc_acc *a)
{
    perm_idx_t *last = a->last, *fnm = a->fnm, *lnm = a->lnm;
    if (last)
    {
        for (size_t i = lo; i < hi; i++)
            last[gid[i]] = (perm_idx_t)i;
    }
    if (!fnm && !lnm)
        return;
    if (var->type == STATA_TYPE_DOUBLE)
    {
        const double *x = var->data.dbl;
        for (size_t i = lo; i < hi; i++)
        {
            if (x[i] < missval)
            {
                uint32_t g = gid[i];
                if (fnm && fnm[g] == CC_NONE)
                    fnm[g] = (perm_idx_t)i;
                if (lnm)
                    lnm[g] = (perm_idx_t)i;
            }
        }
    }
    else
    {
        char **s = var->data.str;
        for (size_t i = lo; i < hi; i++)
        {
            if (s[i] && s[i][0])
            {
                uint32_t g = gid[i];
                if (fnm && fnm[g] == CC_NONE)
                    fnm[g] = (perm_idx_t)i;
                if (lnm)
                    lnm[g] = (perm_idx_t)i;
            }
        }
    }
}

/* ===========================================================================
   One source
   =========================================================================== */

/* Rows [0, nobs) cut at group boundaries into nparts runs of whole groups;
 * nparts == 1 visits the rows in one piece */
typedef struct
{
    const size_t *cut; /* [nparts + 1] */
    int nparts;
} cc_parts;

static inline double cc_group_weight(const cc_acc *a, size_t g)
{
    return a->sw ? a->sw[g] : (double)a->n[g];
}

/* Accumulate everything the targets of source need */
static void cc_accumulate(const cc_job *job, size_t source, const cc_needs *need,
                          const cc_parts *parts, cc_acc *a)
{
    const stata_variable *var = &job->data->vars[job->nby + source];
    const uint32_t *gid = job->groups->gid;
    const double *w = job->weight;
    const double missval = SV_missval;
    const size_t ngroups = job->groups->ngroups;
    const int numeric = var->type == STATA_TYPE_DOUBLE;
    const double *x = numeric ? var->data.dbl : NULL;

    #pragma omp parallel for schedule(static, 1) num_threads(parts->nparts) if (parts->nparts > 1)
    for (int p = 0; p < parts->nparts; p++)
    {
        size_t lo = parts->cut[p], hi = parts->cut[p + 1];
        if (numeric && need->moments)
            cc_pass_moments(x, w, gid, lo, hi, missval, a);
        if (numeric && (need->binary || need->negative))
            cc_pass_flags(x, gid, lo, hi, missval, a);
        if (numeric && (need->min || need->max))
            cc_pass_extremes(x, gid, lo, hi, missval, a);
        if (need->last || need->firstnm || need->lastnm)
            cc_pass_rows(var, gid, lo, hi, missval, a);
    }

    if (numeric && need->spread)
    {
        for (size_t g = 0; g < ngroups; g++)
            a->mean[g] = a->n[g] ? a->swx[g] / cc_group_weight(a, g) : 0.0;
        #pragma omp parallel for schedule(static, 1) num_threads(parts->nparts) if (parts->nparts > 1)
        for (int p = 0; p < parts->nparts; p++)
            cc_pass_spread(x, w, gid, parts->cut[p], parts->cut[p + 1], missval, a);
    }

    if (numeric && need->percent)
    {
        double total = 0.0;
        if (w && job->groups->contiguous)
        {
            /* the weights of the rows in group order, summed as one sequence */
            for (size_t i = 0; i < job->groups->nobs; i++)
                if (x[i] < missval)
                    total += w[i];
        }
        else
        {
            for (size_t g = 0; g < ngroups; g++)
                total += cc_group_weight(a, g);
        }
        a->total = total;
    }
}

/* The value of one numeric target in group g */
static inline double cc_finish(const cc_job *job, const cc_target *target, const cc_acc *a,
                               const double *x, size_t g, double missval)
{
    const cc_wtype wtype = job->wtype;
    const perm_idx_t *first = job->groups->first;
    double n = a->n ? (double)a->n[g] : 0.0;
    double W = a->n ? cc_group_weight(a, g) : 0.0;
    /* the number of observations a standard error or a count refers to */
    double N = (wtype == CC_W_NONE || wtype == CC_W_AWEIGHT) ? n : W;

    switch (target->stat)
    {
    case CC_STAT_MEAN:
        return n > 0 ? cc_storable(a->swx[g] / W, missval) : missval;
    case CC_STAT_SUM:
        if (n == 0)
            return 0.0;
        return cc_storable(wtype == CC_W_AWEIGHT ? a->swx[g] * n / W : a->swx[g], missval);
    case CC_STAT_RAWSUM:
        if (n == 0)
            return 0.0;
        return cc_storable(a->rs ? a->rs[g] : a->swx[g], missval);
    case CC_STAT_COUNT:
        return cc_storable(N, missval);
    case CC_STAT_PERCENT:
        return cc_storable(W / a->total * 100.0, missval);
    case CC_STAT_SD:
    case CC_STAT_SEMEAN:
    {
        if (n == 0)
            return missval;
        /* aweights are scaled to sum to n */
        double ss = wtype == CC_W_AWEIGHT ? a->b[g] * (n / W) : a->b[g];
        double df = wtype == CC_W_AWEIGHT ? n - 1.0 : W - 1.0;
        if (target->stat == CC_STAT_SD)
            return cc_storable(sqrt(ss / df), missval);
        if (a->with_sd)
            return cc_storable(sqrt(ss / df) / sqrt(N), missval);
        return cc_storable(sqrt(ss / (df * N)), missval);
    }
    case CC_STAT_SEBINOMIAL:
    {
        if (n == 0 || (a->flag[g] & 1))
            return missval;
        double mean = a->swx[g] / W;
        return cc_storable(sqrt(mean * (1.0 - mean) / N), missval);
    }
    case CC_STAT_SEPOISSON:
    {
        if (n == 0 || (a->flag[g] & 2))
            return missval;
        double mean = a->swx[g] / W;
        return cc_storable(sqrt(floor(mean * N + 0.5)) / N, missval);
    }
    case CC_STAT_MIN:
        return a->mn[g];
    case CC_STAT_MAX:
        /* without a number in the group, the greatest of its missing values */
        return a->mx[g] > -HUGE_VAL ? a->mx[g] : a->mxall[g];
    case CC_STAT_FIRST:
        return x[first[g]];
    case CC_STAT_LAST:
        return x[a->last[g]];
    case CC_STAT_FIRSTNM:
        return x[a->fnm[g] != CC_NONE ? a->fnm[g] : first[g]];
    case CC_STAT_LASTNM:
        return x[a->lnm[g] != CC_NONE ? a->lnm[g] : a->last[g]];
    default:
        return missval;
    }
}

/* The row whose string group g takes for a string target */
static inline perm_idx_t cc_finish_row(const cc_job *job, const cc_target *target,
                                       const cc_acc *a, size_t g)
{
    const perm_idx_t *first = job->groups->first;
    switch (target->stat)
    {
    case CC_STAT_LAST:
        return a->last[g];
    case CC_STAT_FIRSTNM:
        return a->fnm[g] != CC_NONE ? a->fnm[g] : first[g];
    case CC_STAT_LASTNM:
        return a->lnm[g] != CC_NONE ? a->lnm[g] : a->last[g];
    default:
        return first[g];
    }
}

static stata_retcode cc_source(const cc_job *job, size_t source, const cc_parts *parts)
{
    cc_needs need = cc_source_needs(job, source);
    if (!need.any)
        return STATA_OK;
    const stata_variable *var = &job->data->vars[job->nby + source];
    const size_t ngroups = job->groups->ngroups;
    const double missval = SV_missval;
    cc_acc a;
    stata_retcode rc = cc_acc_alloc(&a, &need, ngroups, job->weight != NULL);
    if (rc != STATA_OK)
        return rc;

    cc_accumulate(job, source, &need, parts, &a);

    for (size_t t = 0; t < job->ntargets; t++)
    {
        const cc_target *target = &job->targets[t];
        if ((size_t)target->source != source || cc_stat_is_quantile(target->stat))
            continue;
        if (var->type == STATA_TYPE_DOUBLE)
        {
            double *out = job->values[t];
            const double *x = var->data.dbl;
            #pragma omp parallel for schedule(static) num_threads(parts->nparts) \
                if (parts->nparts > 1 && ngroups >= 65536)
            for (size_t g = 0; g < ngroups; g++)
                out[g] = cc_finish(job, target, &a, x, g, missval);
        }
        else
        {
            perm_idx_t *out = job->string_rows[t];
            for (size_t g = 0; g < ngroups; g++)
                out[g] = cc_finish_row(job, target, &a, g);
        }
    }
    cc_acc_free(&a);
    return STATA_OK;
}

stata_retcode cc_compute_stats(const cc_job *job)
{
    const cc_groups *groups = job->groups;
    const size_t nobs = groups->nobs, ngroups = groups->ngroups;
    int threads = job->threads > 0 ? job->threads : 1;
    stata_retcode rc = STATA_OK;

    size_t whole[2] = {0, nobs};
    cc_parts one = {whole, 1};

    /* Rows in group order: every thread takes a run of whole groups */
    if (groups->contiguous && groups->starts && ngroups > 1 && threads > 1 &&
        nobs >= CC_STATS_PARALLEL_MIN)
    {
        size_t *cut = malloc(((size_t)threads + 1) * sizeof(*cut));
        if (!cut)
            return STATA_ERR_MEMORY;
        int nparts = 0;
        cut[0] = 0;
        for (int t = 1; t < threads; t++)
        {
            size_t row = nobs * (size_t)t / (size_t)threads;
            size_t at = groups->starts[groups->gid[row]];
            if (at > cut[nparts])
                cut[++nparts] = at;
        }
        cut[++nparts] = nobs;
        cc_parts parts = {cut, nparts};
        for (size_t s = 0; s < job->nsources && rc == STATA_OK; s++)
        {
            if (ctools_cancel_requested())
                rc = STATA_ERR_CANCELLED;
            else
                rc = cc_source(job, s, &parts);
        }
        free(cut);
        return rc;
    }

    /* Otherwise every thread takes a source */
    int team = (size_t)threads < job->nsources ? threads : (int)job->nsources;
    if (team <= 1 || nobs < CC_STATS_PARALLEL_MIN / 4)
    {
        for (size_t s = 0; s < job->nsources && rc == STATA_OK; s++)
        {
            if (ctools_cancel_requested())
                rc = STATA_ERR_CANCELLED;
            else
                rc = cc_source(job, s, &one);
        }
        return rc;
    }
    int failure = STATA_OK;
    #pragma omp parallel for schedule(dynamic, 1) num_threads(team)
    for (size_t s = 0; s < job->nsources; s++)
    {
        int seen;
        #pragma omp atomic read
        seen = failure;
        if (seen != STATA_OK)
            continue;
        stata_retcode src_rc =
            ctools_cancel_requested() ? STATA_ERR_CANCELLED : cc_source(job, s, &one);
        if (src_rc != STATA_OK)
        {
            #pragma omp atomic write
            failure = (int)src_rc;
        }
    }
    return (stata_retcode)failure;
}

#include <math.h>
#include <string.h>
#include "ctools_ols.h"
#include "cpsmatch_stats.h"

typedef struct
{
    size_t n;
    double anchor, mean;
    dd_real shifted_sum, squares;
} cpsmatch_moment;

static int cpsmatch_stat_nonmissing(double x)
{
    return isfinite(x) && x < SV_missval;
}

static double cpsmatch_stat_value(double x)
{
    return cpsmatch_stat_nonmissing(x) ? x : SV_missval;
}

static double cpsmatch_moment_variance(const cpsmatch_moment *moment)
{
    if (moment->n < 2)
        return SV_missval;
    /* Match summarize's squared deviations from its double-precision mean.
     * Re-centering that sum changes results for large-offset outcomes. */
    double variance = (moment->squares.hi + moment->squares.lo) / (moment->n - 1);
    return cpsmatch_stat_value(variance < 0.0 ? 0.0 : variance);
}

int ctools_get_max_threads(void);

/* (a.hi + a.lo) + (b.hi + b.lo) */
static inline dd_real cpsmatch_dd_add(dd_real a, dd_real b)
{
    dd_real s = two_sum(a.hi, b.hi), t = two_sum(a.lo, b.lo);
    s.lo += t.hi;
    double h = s.hi + s.lo;
    s.lo = s.lo - (h - s.hi);
    s.hi = h;
    s.lo += t.lo;
    h = s.hi + s.lo;
    s.lo = s.lo - (h - s.hi);
    s.hi = h;
    return s;
}

/* One contiguous chunk's part of the sums */
typedef struct
{
    size_t n[4], n_weighted;
    dd_real shifted[4], squares[4], weighted_sum, sum_weight, sum_weight_squared;
    double n0_on, n0_off, n1_on, n1_off;
} cpsmatch_stat_part;

enum { CPSMATCH_STAT_CHUNK = 262144, CPSMATCH_STAT_PARTS = 64 };

void cpsmatch_compute_statistics(size_t n, const double *treatment, const double *outcome,
                                 const double *weight, const double *support,
                                 cpsmatch_statistics *result)
{
    /* Full control/treated arms, supported treated, and used controls. The
     * control variance is unweighted; the matched control mean is weighted.
     * Each sum is shifted by its anchor, the first value it takes in row
     * order. Rows split into chunks whose double-double sums add up in
     * chunk order: one chunk (fewer rows than two chunks' worth) is the
     * sequential sum, and more give it unless an exact sum lies within
     * some 1e-25 of a rounding boundary of its double result. */
    cpsmatch_moment moments[4] = {0};
    dd_real weighted_sum = {0}, sum_weight = {0}, sum_weight_squared = {0};
    double weighted_anchor = 0.0;
    size_t n_weighted = 0;
    memset(result, 0, sizeof(*result));

    int have[5] = {0}, found = 0;
    for (size_t i = 0; i < n && outcome && found < 5; i++)
    {
        double t = treatment[i];
        if (t != 0.0 && t != 1.0)
            continue;
        double s = support ? support[i] : 1.0;
        double w = weight ? weight[i] : 0.0;
        double y = outcome[i];
        if (!cpsmatch_stat_nonmissing(y))
            continue;
        int slot[3] = {(int)t, (t == 1.0 && s == 1.0) ? 2 : -1,
                       (t == 0.0 && w > 0.0 && cpsmatch_stat_nonmissing(w)) ? 3 : -1};
        for (int k = 0; k < 3; k++)
            if (slot[k] >= 0 && !have[slot[k]])
            {
                have[slot[k]] = 1;
                moments[slot[k]].anchor = y;
                found++;
            }
        if (slot[2] == 3 && s == 1.0 && !have[4])
        {
            have[4] = 1;
            weighted_anchor = y;
            found++;
        }
    }

    int chunks = 1;
    if (n >= 2 * (size_t)CPSMATCH_STAT_CHUNK)
    {
        chunks = ctools_get_max_threads();
        if (chunks < 1)
            chunks = 1;
        if (chunks > CPSMATCH_STAT_PARTS)
            chunks = CPSMATCH_STAT_PARTS;
    }
    cpsmatch_stat_part parts[CPSMATCH_STAT_PARTS];
    memset(parts, 0, sizeof(parts));
#pragma omp parallel for schedule(static, 1) num_threads(chunks) if (chunks > 1)
    for (int c = 0; c < chunks; c++)
    {
        cpsmatch_stat_part *part = &parts[c];
        const size_t lo = n * (size_t)c / (size_t)chunks, hi = n * ((size_t)c + 1) / (size_t)chunks;
        for (size_t i = lo; i < hi; i++)
        {
            double t = treatment[i];
            if (t != 0.0 && t != 1.0)
                continue;
            double s = support ? support[i] : 1.0;
            double w = weight ? weight[i] : 0.0;
            if (t == 0.0)
            {
                part->n0_on += s == 1.0;
                part->n0_off += s == 0.0;
            }
            else
            {
                part->n1_on += s == 1.0;
                part->n1_off += s == 0.0;
            }
            if (!outcome)
                continue;
            if (t == 0.0 && cpsmatch_stat_nonmissing(w))
            {
                /* gen double weight_squared = _weight^2 rounds each product
                 * to double and turns overflow into missing before
                 * summarize. */
                double square = w * w;
                if (cpsmatch_stat_nonmissing(square))
                    part->sum_weight_squared = dd_add_d(part->sum_weight_squared, square);
            }
            double y = outcome[i];
            if (!cpsmatch_stat_nonmissing(y))
                continue;
            int k = (int)t;
            part->shifted[k] = dd_add_d(part->shifted[k], y - moments[k].anchor);
            part->n[k]++;
            if (t == 1.0 && s == 1.0)
            {
                part->shifted[2] = dd_add_d(part->shifted[2], y - moments[2].anchor);
                part->n[2]++;
            }
            if (t == 0.0 && w > 0.0 && cpsmatch_stat_nonmissing(w))
            {
                part->shifted[3] = dd_add_d(part->shifted[3], y - moments[3].anchor);
                part->n[3]++;
                if (s == 1.0)
                {
                    part->weighted_sum = dd_add_d(part->weighted_sum, w * (y - weighted_anchor));
                    part->sum_weight = dd_add_d(part->sum_weight, w);
                    part->n_weighted++;
                }
            }
        }
    }
    for (int c = 0; c < chunks; c++)
    {
        const cpsmatch_stat_part *part = &parts[c];
        result->n0_on += part->n0_on;
        result->n0_off += part->n0_off;
        result->n1_on += part->n1_on;
        result->n1_off += part->n1_off;
        for (int k = 0; k < 4; k++)
        {
            moments[k].shifted_sum =
                c ? cpsmatch_dd_add(moments[k].shifted_sum, part->shifted[k]) : part->shifted[k];
            moments[k].n += part->n[k];
        }
        weighted_sum = c ? cpsmatch_dd_add(weighted_sum, part->weighted_sum) : part->weighted_sum;
        sum_weight = c ? cpsmatch_dd_add(sum_weight, part->sum_weight) : part->sum_weight;
        sum_weight_squared = c ? cpsmatch_dd_add(sum_weight_squared, part->sum_weight_squared)
                               : part->sum_weight_squared;
        n_weighted += part->n_weighted;
    }
    for (size_t j = 0; j < 4; j++)
    {
        cpsmatch_moment *moment = &moments[j];
        moment->mean =
            moment->n
                ? cpsmatch_stat_value(moment->anchor +
                                      (moment->shifted_sum.hi + moment->shifted_sum.lo) / moment->n)
                : SV_missval;
    }
    if (outcome)
    {
#pragma omp parallel for schedule(static, 1) num_threads(chunks) if (chunks > 1)
        for (int c = 0; c < chunks; c++)
        {
            cpsmatch_stat_part *part = &parts[c];
            const size_t lo = n * (size_t)c / (size_t)chunks, hi = n * ((size_t)c + 1) / (size_t)chunks;
            for (size_t i = lo; i < hi; i++)
            {
                double t = treatment[i], y = outcome[i];
                if ((t != 0.0 && t != 1.0) || !cpsmatch_stat_nonmissing(y))
                    continue;
                double s = support ? support[i] : 1.0;
                double w = weight ? weight[i] : 0.0;
                int slot[3] = {(int)t, (t == 1.0 && s == 1.0) ? 2 : -1,
                               (t == 0.0 && w > 0.0 && cpsmatch_stat_nonmissing(w)) ? 3 : -1};
                for (int k = 0; k < 3; k++)
                {
                    if (slot[k] < 0)
                        continue;
                    double d = y - moments[slot[k]].mean;
                    dd_real square = two_prod(d, d);
                    part->squares[slot[k]] = dd_add_d(part->squares[slot[k]], square.hi);
                    part->squares[slot[k]] = dd_add_d(part->squares[slot[k]], square.lo);
                }
            }
        }
        for (int c = 0; c < chunks; c++)
            for (int k = 0; k < 4; k++)
                moments[k].squares =
                    c ? cpsmatch_dd_add(moments[k].squares, parts[c].squares[k]) : parts[c].squares[k];
    }

    result->y1 = moments[2].mean;
    result->y0 = n_weighted
                     ? cpsmatch_stat_value(weighted_anchor + (weighted_sum.hi + weighted_sum.lo) /
                                                                 (sum_weight.hi + sum_weight.lo))
                     : SV_missval;
    result->n1 = (double)moments[2].n;
    result->n0 = (double)n_weighted;
    result->var1 = cpsmatch_moment_variance(&moments[2]);
    result->var0 = cpsmatch_moment_variance(&moments[3]);
    result->weight_squared = cpsmatch_stat_value(sum_weight_squared.hi + sum_weight_squared.lo);
    result->y1_u = moments[1].mean;
    result->y0_u = moments[0].mean;
    result->n1_u = (double)moments[1].n;
    result->n0_u = (double)moments[0].n;
    result->var1_u = cpsmatch_moment_variance(&moments[1]);
    result->var0_u = cpsmatch_moment_variance(&moments[0]);
}

ST_retcode cpsmatch_save_statistics(const cpsmatch_statistics *result)
{
    const char *names[] = {"_cpsmatch_y1",     "_cpsmatch_y0",    "_cpsmatch_n1",
                           "_cpsmatch_n0",     "_cpsmatch_var1",  "_cpsmatch_var0",
                           "_cpsmatch_w2",     "_cpsmatch_y1_u",  "_cpsmatch_y0_u",
                           "_cpsmatch_n1_u",   "_cpsmatch_n0_u",  "_cpsmatch_var1_u",
                           "_cpsmatch_var0_u", "_cpsmatch_n0_on", "_cpsmatch_n0_off",
                           "_cpsmatch_n1_on",  "_cpsmatch_n1_off"};
    const double values[] = {result->y1,
                             result->y0,
                             result->n1,
                             result->n0,
                             result->var1,
                             result->var0,
                             result->weight_squared,
                             result->y1_u,
                             result->y0_u,
                             result->n1_u,
                             result->n0_u,
                             result->var1_u,
                             result->var0_u,
                             result->n0_on,
                             result->n0_off,
                             result->n1_on,
                             result->n1_off};
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++)
    {
        ST_retcode rc = SF_scal_save((char *)names[i], values[i]);
        if (rc)
            return rc;
    }
    return 0;
}

/* ========================================================================
 * Exact sums and the summarize/generate passes of _cpsmatch_effects
 * ======================================================================== */

#ifdef _OPENMP
#include <omp.h>
#endif

/* Digit k holds bits [32k, 32k + 32) of the sum times 2^1074, so every
 * double is an integer there. An add puts less than 2^32 into a digit, and
 * carries are pushed up every 2^30 adds, before an int64 digit can
 * overflow; normalized digits below the top one lie in [0, 2^32). */
static void cpsmatch_xsum_norm(cpsmatch_xsum *a)
{
    for (int k = 0; k < CPSMATCH_XDIGITS - 1; k++)
    {
        int64_t carry = a->d[k] >> 32; /* arithmetic shift: floor */
        a->d[k] = (int64_t)((uint64_t)a->d[k] & 0xffffffffu);
        a->d[k + 1] += carry;
    }
    a->pending = 0;
}

void cpsmatch_xsum_add(cpsmatch_xsum *a, double x)
{
    uint64_t bits;
    memcpy(&bits, &x, sizeof(bits));
    unsigned e = (unsigned)(bits >> 52) & 0x7ffu;
    uint64_t m = bits & 0xfffffffffffffULL;
    if (e)
        m |= (uint64_t)1 << 52;
    else if (!m)
        return;
    else
        e = 1;
    const unsigned b = e - 1, k = b >> 5, s = b & 31u;
    const unsigned __int128 v = (unsigned __int128)m << s;
    int64_t lo = (int64_t)((uint64_t)v & 0xffffffffu);
    int64_t mid = (int64_t)((uint64_t)(v >> 32) & 0xffffffffu);
    int64_t hi = (int64_t)(uint64_t)(v >> 64);
    if (bits >> 63)
    {
        lo = -lo;
        mid = -mid;
        hi = -hi;
    }
    a->d[k] += lo;
    a->d[k + 1] += mid;
    a->d[k + 2] += hi;
    if (++a->pending == ((uint32_t)1 << 30))
        cpsmatch_xsum_norm(a);
}

void cpsmatch_xsum_merge(cpsmatch_xsum *into, const cpsmatch_xsum *from)
{
    cpsmatch_xsum f = *from;
    cpsmatch_xsum_norm(&f);
    cpsmatch_xsum_norm(into);
    for (int k = 0; k < CPSMATCH_XDIGITS; k++)
        into->d[k] += f.d[k];
    into->pending = 1;
}

static int cpsmatch_clz128(unsigned __int128 w)
{
    uint64_t hi = (uint64_t)(w >> 64), lo = (uint64_t)w;
    return hi ? __builtin_clzll(hi) : 64 + __builtin_clzll(lo);
}

/* Round to nearest, ties to even */
double cpsmatch_xsum_value(const cpsmatch_xsum *in)
{
    cpsmatch_xsum a = *in;
    cpsmatch_xsum_norm(&a);
    const int neg = a.d[CPSMATCH_XDIGITS - 1] < 0;
    if (neg)
    {
        for (int k = 0; k < CPSMATCH_XDIGITS; k++)
            a.d[k] = -a.d[k];
        cpsmatch_xsum_norm(&a);
    }
    int t = CPSMATCH_XDIGITS - 1;
    while (t >= 0 && a.d[t] == 0)
        t--;
    if (t < 0)
        return 0.0;
    unsigned __int128 w = 0;
    for (int k = t; k >= t - 2; k--)
        w = (w << 32) | (k >= 0 ? (uint64_t)a.d[k] : 0);
    int sticky = 0;
    for (int k = t - 3; k >= 0 && !sticky; k--)
        sticky = a.d[k] != 0;
    const int base = 32 * (t - 2) - 1074; /* exponent of w's lowest bit */
    const int bits = 128 - cpsmatch_clz128(w);
    int low = base + bits - 53; /* exponent of the result's lowest bit */
    if (low < -1074)
        low = -1074;
    const int drop = low - base;
    unsigned __int128 q;
    if (drop > 0)
    {
        const unsigned __int128 rem = w & (((unsigned __int128)1 << drop) - 1);
        const unsigned __int128 half = (unsigned __int128)1 << (drop - 1);
        q = w >> drop;
        if (rem > half || (rem == half && (sticky || (q & 1))))
            q++;
    }
    else
    {
        q = w << (-drop);
    }
    const double r = ldexp((double)(uint64_t)q, low);
    return neg ? -r : r;
}

/* Stata's arithmetic on possibly missing values */
static inline int sm_mis(double x)
{
    return !(x < SV_missval);
}
static inline double sm_fin(double r)
{
    return (isfinite(r) && r < SV_missval) ? r : SV_missval;
}
static inline double sm_add(double a, double b)
{
    return (sm_mis(a) || sm_mis(b)) ? SV_missval : sm_fin(a + b);
}
static inline double sm_sub(double a, double b)
{
    return (sm_mis(a) || sm_mis(b)) ? SV_missval : sm_fin(a - b);
}
static inline double sm_mul(double a, double b)
{
    return (sm_mis(a) || sm_mis(b)) ? SV_missval : sm_fin(a * b);
}
static inline double sm_div(double a, double b)
{
    return (sm_mis(a) || sm_mis(b) || b == 0.0) ? SV_missval : sm_fin(a / b);
}
/* x^2 is x*x in Stata */
static inline double sm_sq(double a)
{
    return sm_mul(a, a);
}
/* max() ignores missing arguments */
static inline double sm_max(double a, double b)
{
    if (sm_mis(a))
        return b;
    if (sm_mis(b))
        return a;
    return a > b ? a : b;
}
/* generate's default storage type: float, missing beyond its range */
static inline double sm_float(double x)
{
    if (sm_mis(x))
        return SV_missval;
    const float f = (float)x;
    if (!isfinite(f) || fabsf(f) > 0x1.fffffep+126f)
        return SV_missval;
    return (double)f;
}

static double cpsmatch_xsum_total(const cpsmatch_xsum *work, int threads, int j)
{
    cpsmatch_xsum total;
    memset(&total, 0, sizeof(total));
    for (int t = 0; t < threads; t++)
        cpsmatch_xsum_merge(&total, &work[(size_t)t * CPSMATCH_EFFECT_SUMS + j]);
    return sm_fin(cpsmatch_xsum_value(&total));
}

static double cpsmatch_xmean(const cpsmatch_xsum *work, int threads, int j, size_t n)
{
    return n ? sm_fin(cpsmatch_xsum_total(work, threads, j) / (double)n) : SV_missval;
}

static double cpsmatch_xvar(const cpsmatch_xsum *work, int threads, int j, size_t n)
{
    return n >= 2 ? sm_fin(cpsmatch_xsum_total(work, threads, j) / (double)(n - 1)) : SV_missval;
}

/* Static row ranges: [lo, hi) of thread tid */
#define CPSMATCH_RANGE(n, tid, team, lo, hi)                                                       \
    const size_t lo = (n) * (tid) / (team), hi = (n) * ((tid) + 1) / (team)

void cpsmatch_compute_effects(size_t n, const double *arm, const double *support,
                              const double *weight, const double *y, const double *matched,
                              const double *self, int ai, int ate, int samplevar, int altvariance,
                              int neighbors, int threads, cpsmatch_xsum *work,
                              cpsmatch_effects *out)
{
    const double miss = SV_missval;
    if (threads < 1)
        threads = 1;
    memset(work, 0, (size_t)threads * CPSMATCH_EFFECT_SUMS * sizeof(*work));
    /* summarize v if _treated == 1 (A) and == 0 (B), each also on support
     * (C, D); summarize _v on support, treated (E) and controls (F) */
    size_t nA = 0, nB = 0, nC = 0, nD = 0, nE = 0, nF = 0;
#pragma omp parallel num_threads(threads) reduction(+ : nA, nB, nC, nD, nE, nF)
    {
        size_t tid = 0, team = 1;
#ifdef _OPENMP
        tid = (size_t)omp_get_thread_num();
        team = (size_t)omp_get_num_threads();
#endif
        cpsmatch_xsum *w = work + tid * CPSMATCH_EFFECT_SUMS;
        CPSMATCH_RANGE(n, tid, team, lo, hi);
        for (size_t i = lo; i < hi; i++)
        {
            const double t = arm[i], v = y[i];
            const int on = support[i] == 1.0;
            if (!sm_mis(v))
            {
                if (t == 1.0)
                {
                    cpsmatch_xsum_add(&w[0], v);
                    nA++;
                    if (on)
                    {
                        cpsmatch_xsum_add(&w[2], v);
                        nC++;
                    }
                }
                else if (t == 0.0)
                {
                    cpsmatch_xsum_add(&w[1], v);
                    nB++;
                    if (on)
                    {
                        cpsmatch_xsum_add(&w[3], v);
                        nD++;
                    }
                }
            }
            const double mv = matched ? matched[i] : miss;
            if (on && !sm_mis(mv))
            {
                if (t == 1.0)
                {
                    cpsmatch_xsum_add(&w[4], mv);
                    nE++;
                }
                else if (t == 0.0)
                {
                    cpsmatch_xsum_add(&w[5], mv);
                    nF++;
                }
            }
        }
    }
    out->u1u = cpsmatch_xmean(work, threads, 0, nA);
    out->n1u = (double)nA;
    out->u0u = cpsmatch_xmean(work, threads, 1, nB);
    out->n0u = (double)nB;
    out->m1t = cpsmatch_xmean(work, threads, 2, nC);
    out->n1 = (double)nC;
    out->m0t = cpsmatch_xmean(work, threads, 4, nE);
    out->m1u = cpsmatch_xmean(work, threads, 5, nF);
    out->m0u = cpsmatch_xmean(work, threads, 3, nD);
    out->n0 = (double)nD;

    /* The scalars _cpsmatch_effects forms from these */
    const double att = sm_sub(out->m1t, out->m0t);
    const double atu = ate ? sm_sub(out->m1u, out->m0u) : miss;
    const double n1 = out->n1, n0 = ate ? out->n0 : 0.0;
    const double effect =
        ate ? sm_add(sm_div(sm_mul(att, n1), n0 + n1), sm_div(sm_mul(atu, n0), n0 + n1)) : miss;
    const double frac = (double)ai / ((double)ai + 1.0), nb = (double)neighbors;
    const double u1u = out->u1u, u0u = out->u0u, m1t = out->m1t, m0u = out->m0u;

    memset(work, 0, (size_t)threads * CPSMATCH_EFFECT_SUMS * sizeof(*work));
#pragma omp parallel num_threads(threads)
    {
        size_t tid = 0, team = 1;
#ifdef _OPENMP
        tid = (size_t)omp_get_thread_num();
        team = (size_t)omp_get_num_threads();
#endif
        cpsmatch_xsum *w = work + tid * CPSMATCH_EFFECT_SUMS;
        CPSMATCH_RANGE(n, tid, team, lo, hi);
        for (size_t i = lo; i < hi; i++)
        {
            const double t = arm[i], v = y[i];
            const int on = support[i] == 1.0;
            if (!sm_mis(v))
            {
                double d;
                if (t == 1.0)
                {
                    d = v - u1u;
                    cpsmatch_xsum_add(&w[0], d * d);
                    if (on)
                    {
                        d = v - m1t;
                        cpsmatch_xsum_add(&w[2], d * d);
                    }
                }
                else if (t == 0.0)
                {
                    d = v - u0u;
                    cpsmatch_xsum_add(&w[1], d * d);
                    if (on)
                    {
                        d = v - m0u;
                        cpsmatch_xsum_add(&w[3], d * d);
                    }
                }
            }
            if (!ai || !on)
                continue;
            /* generate w = max(_weight, 0), shat, vt, vu and ve, all float */
            const double mv = matched ? matched[i] : miss, sv = self ? self[i] : miss;
            const double wt = sm_float(sm_max(weight[i], 0.0));
            const double shat =
                altvariance ? sm_float(sv) : sm_float(sm_mul(frac, sm_sq(sm_sub(v, sv))));
            const double one_t = sm_sub(1.0, t);
            double vt, vu = miss, ve = miss;
            if (samplevar)
            {
                vt = sm_float(sm_mul(shat, sm_sq(sm_sub(t, sm_mul(one_t, wt)))));
                if (ate)
                {
                    vu = sm_float(sm_mul(shat, sm_sq(sm_sub(one_t, sm_mul(t, wt)))));
                    ve = sm_float(sm_mul(shat, sm_sq(sm_add(1.0, wt))));
                }
            }
            else
            {
                const double wpart = sm_sub(sm_sq(wt), sm_div(wt, nb));
                vt = sm_float(sm_add(sm_max(0.0, sm_mul(t, sm_sq(sm_sub(sm_sub(v, mv), att)))),
                                     sm_mul(sm_mul(one_t, wpart), shat)));
                if (ate)
                {
                    vu = sm_float(
                        sm_add(sm_max(0.0, sm_mul(one_t, sm_sq(sm_sub(sm_sub(mv, v), atu)))),
                               sm_mul(sm_mul(t, wpart), shat)));
                    const double inner = sm_sub(
                        sm_add(sm_mul(t, sm_sub(v, mv)), sm_mul(one_t, sm_sub(mv, v))), effect);
                    ve = sm_float(sm_add(
                        sm_max(0.0, sm_sq(inner)),
                        sm_mul(sm_sub(sm_add(sm_sq(wt), sm_mul(2.0, wt)), sm_div(wt, nb)), shat)));
                }
            }
            if (!sm_mis(vt))
                cpsmatch_xsum_add(&w[4], vt);
            if (!sm_mis(vu))
                cpsmatch_xsum_add(&w[5], vu);
            if (!sm_mis(ve))
                cpsmatch_xsum_add(&w[6], ve);
        }
    }
    out->var1u = cpsmatch_xvar(work, threads, 0, nA);
    out->var0u = cpsmatch_xvar(work, threads, 1, nB);
    out->var1 = cpsmatch_xvar(work, threads, 2, nC);
    out->var0 = cpsmatch_xvar(work, threads, 3, nD);
    out->sumvt = ai ? cpsmatch_xsum_total(work, threads, 4) : miss;
    out->sumvu = ai && ate ? cpsmatch_xsum_total(work, threads, 5) : miss;
    out->sumve = ai && ate ? cpsmatch_xsum_total(work, threads, 6) : miss;
}

void cpsmatch_compute_effect_weights(size_t n, const double *arm, const double *support,
                                     const double *weight, int threads, cpsmatch_xsum *work,
                                     cpsmatch_effect_weights *out)
{
    if (threads < 1)
        threads = 1;
    memset(work, 0, (size_t)threads * CPSMATCH_EFFECT_SUMS * sizeof(*work));
#pragma omp parallel num_threads(threads)
    {
        size_t tid = 0, team = 1;
#ifdef _OPENMP
        tid = (size_t)omp_get_thread_num();
        team = (size_t)omp_get_num_threads();
#endif
        cpsmatch_xsum *w = work + tid * CPSMATCH_EFFECT_SUMS;
        CPSMATCH_RANGE(n, tid, team, lo, hi);
        for (size_t i = lo; i < hi; i++)
        {
            const double t = arm[i];
            if (support[i] != 1.0 || (t != 0.0 && t != 1.0))
                continue;
            /* generate double w2 = _weight^2 if _support == 1, later
             * replace w2 = (1 + _weight)^2 if _support == 1 */
            const double w2 = sm_sq(weight[i]), w1 = sm_sq(sm_add(1.0, weight[i]));
            const int k = t == 0.0 ? 0 : 1;
            if (!sm_mis(w2))
                cpsmatch_xsum_add(&w[k], w2);
            if (!sm_mis(w1))
                cpsmatch_xsum_add(&w[2 + k], w1);
        }
    }
    out->w2c = cpsmatch_xsum_total(work, threads, 0);
    out->w2t = cpsmatch_xsum_total(work, threads, 1);
    out->w1c = cpsmatch_xsum_total(work, threads, 2);
    out->w1t = cpsmatch_xsum_total(work, threads, 3);
}

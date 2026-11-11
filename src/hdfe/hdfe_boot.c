/*
 * hdfe_boot.c
 *
 * Pairs-bootstrap engine shared by the HDFE estimators (see hdfe_boot.h and
 * docs/development/BOOTSTRAP_PLAN.md). No SPI calls.
 * Part of the ctools Stata plugin suite
 */

#include "hdfe_boot.h"
#include "../ctools_cancel.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "../ctools_rng.h"

/* ========================================================================
 * Resampling frame and draws
 * ======================================================================== */

static ST_int stratum_lo(const hdfe_boot_design *d, ST_int s)
{
    return d->stratum_start ? d->stratum_start[s] : 0;
}

static ST_int stratum_hi(const hdfe_boot_design *d, ST_int s)
{
    return d->stratum_start ? d->stratum_start[s + 1] : d->U;
}

ST_retcode hdfe_boot_design_check(const hdfe_boot_design *d)
{
    if (!d || d->U < 1 || d->S < 1 || (!d->stratum_start && d->S != 1))
        return 198;
    if (d->stratum_start && (d->stratum_start[0] != 0 || d->stratum_start[d->S] != d->U))
        return 198;
    for (ST_int s = 0; s < d->S; s++)
    {
        ST_int lo = stratum_lo(d, s), hi = stratum_hi(d, s);
        if (hi < lo)
            return 198;
        if (hi == lo)
            return 2001;
        if (d->unit_freq)
        {
            double total = 0.0;
            for (ST_int u = lo; u < hi; u++)
            {
                double f = d->unit_freq[u];
                /* at most INT_MAX / 2, so that a unit's drawn count, within
                 * a few sqrt(f) of f, fits an ST_int */
                if (!isfinite(f) || f < 0.0 || f != floor(f) || f > (double)(INT_MAX / 2))
                    return 198;
                total += f;
            }
            if (total < 1.0)
                return 2001;
            /* exact in double, and in draws(fast)'s 64-bit counts */
            if (total > 9007199254740992.0)
                return 198;
        }
    }
    return 0;
}

void hdfe_boot_draw_fast(const hdfe_boot_design *d, ST_int r, ST_int *unit_mult)
{
    ctools_rng rng;
    /* One generator per replicate, keyed on (seed, r): the draws do not
     * depend on which replicates run together or in what order. */
    ctools_rng_seed(&rng, d->seed ^ (0x9E3779B97F4A7C15ULL * ((uint64_t)r + 1u)));
    memset(unit_mult, 0, (size_t)d->U * sizeof(*unit_mult));
    for (ST_int s = 0; s < d->S; s++)
    {
        ST_int lo = stratum_lo(d, s), n = stratum_hi(d, s) - lo;
        if (!d->unit_freq)
        {
            for (ST_int t = 0; t < n; t++)
                unit_mult[lo + (ST_int)ctools_rng_uniform(&rng, (size_t)n)]++;
            continue;
        }
        /* Multinomial(F, freq / F) by conditional binomials */
        double mass = 0.0;
        for (ST_int u = lo; u < lo + n; u++)
            mass += d->unit_freq[u];
        int64_t remaining = (int64_t)mass;
        for (ST_int u = lo; u < lo + n && remaining > 0; u++)
        {
            double f = d->unit_freq[u];
            int64_t k;
            if (u == lo + n - 1 || mass <= f)
                k = remaining;
            else
            {
                double p = f / mass;
                k = ctools_rng_binomial(&rng, remaining, p < 1.0 ? p : 1.0);
            }
            /* a unit's count beyond INT_MAX would need a draw some 10^4
             * standard deviations above its frequency (at most INT_MAX / 2) */
            unit_mult[u] = k > INT_MAX ? INT_MAX : (ST_int)k;
            remaining -= k;
            mass -= f;
        }
    }
}

typedef struct
{
    ST_int draw; /* drawn unit, 1-based within the stratum */
    double tie;
    ST_int idx; /* unit, 0-based within the stratum */
} boot_draw_key;

static int boot_draw_key_cmp(const void *a, const void *b)
{
    const boot_draw_key *x = a, *y = b;
    if (x->draw != y->draw)
        return x->draw < y->draw ? -1 : 1;
    if (x->tie != y->tie)
        return x->tie < y->tie ? -1 : 1;
    return (x->idx > y->idx) - (x->idx < y->idx);
}

ST_retcode hdfe_boot_map_bsample(const hdfe_boot_design *d, const ST_double *u_draw,
                                 const ST_double *u_tie, int cluster_mode, ST_int *unit_mult)
{
    for (ST_int u = 0; u < d->U; u++)
    {
        if (!(u_draw[u] >= 0.0 && u_draw[u] < 1.0) || !(u_tie[u] >= 0.0 && u_tie[u] < 1.0))
            return 198;
    }
    ST_int widest = 0;
    for (ST_int s = 0; s < d->S; s++)
    {
        ST_int n = stratum_hi(d, s) - stratum_lo(d, s);
        if (n > widest)
            widest = n;
    }
    boot_draw_key *keys = malloc((size_t)widest * sizeof(*keys));
    ST_int *runs = malloc((size_t)widest * sizeof(*runs));
    if (!keys || !runs)
    {
        free(keys);
        free(runs);
        return 920;
    }
    for (ST_int s = 0; s < d->S; s++)
    {
        ST_int lo = stratum_lo(d, s), n = stratum_hi(d, s) - lo;
        for (ST_int k = 0; k < n; k++)
        {
            /* bsample's int(uniform() * n + 1), rounded as Stata rounds it:
             * the product, then the sum, each in double */
            volatile double t = u_draw[lo + k] * (double)n;
            t = t + 1.0;
            keys[k].draw = (ST_int)t;
            keys[k].tie = u_tie[lo + k];
            keys[k].idx = k;
        }
        /* The cluster mapping sorts on the unit too; the observation
         * mapping's sort has no third key, but a tie in both uniforms is a
         * probability-zero event, so the index only makes the order total. */
        qsort(keys, (size_t)n, sizeof(*keys), boot_draw_key_cmp);
        /* runs[q]: the count of a run of equal draws, at its last position */
        ST_int run = 0;
        for (ST_int q = 0; q < n; q++)
        {
            run = (q > 0 && keys[q].draw == keys[q - 1].draw) ? run + 1 : 1;
            runs[q] = (q == n - 1 || keys[q + 1].draw != keys[q].draw) ? run : 0;
        }
        for (ST_int q = 0; q < n; q++)
        {
            if (cluster_mode)
                /* bsample maps the counts back through the permutation: the
                 * unit numbered q takes the count at position p[q] */
                unit_mult[lo + q] = runs[keys[q].idx];
            else
                /* the last row of each run carries the run's count */
                unit_mult[lo + keys[q].idx] = runs[q];
        }
    }
    free(keys);
    free(runs);
    return 0;
}

double hdfe_boot_draw_total(const hdfe_boot_design *d)
{
    if (!d->unit_freq)
        return (double)d->U;
    double total = 0.0;
    for (ST_int u = 0; u < d->U; u++)
        total += d->unit_freq[u];
    return total;
}

ST_int hdfe_boot_draw_units(const hdfe_boot_design *d)
{
    const double total = hdfe_boot_draw_total(d);
    return total > (double)INT_MAX ? -1 : (ST_int)total;
}

ST_retcode hdfe_boot_map_stata(const hdfe_boot_design *d, const ST_double *u_draw,
                               const ST_double *u_tie, int cluster_mode, ST_int *unit_mult)
{
    if (!d->unit_freq || cluster_mode)
        return hdfe_boot_map_bsample(d, u_draw, u_tie, cluster_mode, unit_mult);
    const ST_int Ue = hdfe_boot_draw_units(d);
    if (Ue < 1)
        return 198;
    ST_int *owner = malloc((size_t)Ue * sizeof(*owner));
    ST_int *start = malloc((size_t)(d->S + 1) * sizeof(*start));
    ST_int *mult = malloc((size_t)Ue * sizeof(*mult));
    if (!owner || !start || !mult)
    {
        free(owner);
        free(start);
        free(mult);
        return 920;
    }
    /* The expanded frame: within a stratum, each unit's first copy, then
     * the further copies, unit by unit */
    ST_int pos = 0;
    for (ST_int s = 0; s < d->S; s++)
    {
        ST_int lo = stratum_lo(d, s), hi = stratum_hi(d, s);
        start[s] = pos;
        for (ST_int u = lo; u < hi; u++)
        {
            if (d->unit_freq[u] >= 1.0)
                owner[pos++] = u;
        }
        for (ST_int u = lo; u < hi; u++)
        {
            for (ST_int c = 1; c < (ST_int)d->unit_freq[u]; c++)
                owner[pos++] = u;
        }
    }
    start[d->S] = pos;
    hdfe_boot_design e = {Ue, d->S, start, NULL, 0};
    ST_retcode rc = pos == Ue ? hdfe_boot_design_check(&e) : 198;
    if (!rc)
        rc = hdfe_boot_map_bsample(&e, u_draw, u_tie, 0, mult);
    if (!rc)
    {
        memset(unit_mult, 0, (size_t)d->U * sizeof(*unit_mult));
        for (ST_int p = 0; p < Ue; p++)
            unit_mult[owner[p]] += mult[p];
    }
    free(owner);
    free(start);
    free(mult);
    return rc;
}

/* ========================================================================
 * Accumulation
 * ======================================================================== */

ST_retcode hdfe_boot_acc_init(hdfe_boot_acc *acc, ST_int K_b, ST_int L, ST_int Nd,
                              ST_int R_cap)
{
    memset(acc, 0, sizeof(*acc));
    if (K_b < 0 || L < 0 || Nd < 0 || R_cap < 1)
        return 198;
    acc->K_b = K_b;
    acc->L = L;
    acc->Nd = Nd;
    acc->R_cap = R_cap;
    acc->b = malloc((size_t)R_cap * (size_t)(K_b > 0 ? K_b : 1) * sizeof(*acc->b));
    acc->status = malloc((size_t)R_cap * sizeof(*acc->status));
    if (L > 0)
    {
        acc->fe_mean = calloc((size_t)L, sizeof(*acc->fe_mean));
        acc->fe_m2 = calloc((size_t)L, sizeof(*acc->fe_m2));
        acc->fe_n = calloc((size_t)L, sizeof(*acc->fe_n));
    }
    if (Nd > 0)
    {
        acc->d_mean = calloc((size_t)Nd, sizeof(*acc->d_mean));
        acc->d_m2 = calloc((size_t)Nd, sizeof(*acc->d_m2));
        acc->d_n = calloc((size_t)Nd, sizeof(*acc->d_n));
    }
    if (!acc->b || !acc->status || (L > 0 && (!acc->fe_mean || !acc->fe_m2 || !acc->fe_n)) ||
        (Nd > 0 && (!acc->d_mean || !acc->d_m2 || !acc->d_n)))
    {
        hdfe_boot_acc_free(acc);
        return 920;
    }
    return 0;
}

void hdfe_boot_acc_free(hdfe_boot_acc *acc)
{
    if (!acc)
        return;
    free(acc->b);
    free(acc->status);
    free(acc->fe_mean);
    free(acc->fe_m2);
    free(acc->fe_n);
    free(acc->d_mean);
    free(acc->d_m2);
    free(acc->d_n);
    free(acc->fe_draws);
    memset(acc, 0, sizeof(*acc));
}

ST_retcode hdfe_boot_acc_keep_fe(hdfe_boot_acc *acc, size_t max_bytes)
{
    if (!acc || acc->L < 1 || acc->R != 0)
        return 198;
    const size_t n = (size_t)acc->L * (size_t)acc->R_cap;
    if (n > max_bytes / sizeof(float))
        return 909;
    free(acc->fe_draws);
    if (!(acc->fe_draws = malloc(n * sizeof(float))))
        return 920;
    for (size_t i = 0; i < n; i++)
        acc->fe_draws[i] = NAN;
    return 0;
}

static int fe_cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

ST_retcode hdfe_boot_fe_ci(const hdfe_boot_acc *acc, ST_double level, ST_double *lb,
                           ST_double *ub)
{
    if (!acc || !acc->fe_draws || !(level > 0.0 && level < 100.0))
        return 198;
    const ST_int R = acc->R;
    double *x = malloc((size_t)(R > 0 ? R : 1) * sizeof(double));
    if (!x)
        return 920;
    for (ST_int i = 0; i < acc->L; i++)
    {
        const float *f = acc->fe_draws + (size_t)i * (size_t)acc->R_cap;
        ST_int n = 0;
        for (ST_int r = 0; r < R; r++)
        {
            if (!isnan(f[r]))
                x[n++] = (double)f[r];
        }
        if (n < 2)
        {
            lb[i] = ub[i] = NAN;
            continue;
        }
        qsort(x, (size_t)n, sizeof(double), fe_cmp_double);
        lb[i] = hdfe_boot_pctile(x, n, (100.0 - level) / 2.0);
        ub[i] = hdfe_boot_pctile(x, n, (100.0 + level) / 2.0);
    }
    free(x);
    return 0;
}

/*
 * Welford over replicates reps[0 .. n) for each of len elements (the fixed
 * effects, or d when use_d), skipping the replicates whose status[j] has
 * no fixed effects and NaN values. Each element folds the replicates in
 * order, so the sums are those of folding one replicate at a time,
 * whatever the number of threads over the elements.
 */
static void welford_batch(const hdfe_boot_rep *reps, const ST_int *status, ST_int n, int use_d,
                          ST_int len, ST_double *mean, ST_double *m2, ST_int *count,
                          ST_int threads)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)                                    \
    if (threads > 1 && (size_t)len * (size_t)n > 100000)
#endif
    for (ST_int i = 0; i < len; i++)
    {
        ST_double mu = mean[i], s = m2[i];
        ST_int c = count[i];
        for (ST_int j = 0; j < n; j++)
        {
            /* The fixed effects come from the replicate's own fit, which is
             * valid when only the omitted-coefficient pattern differs (plan
             * 5.5). */
            if (status[j] != HDFE_BOOT_OK && status[j] != HDFE_BOOT_COLLINEAR)
                continue;
            const ST_double *x = use_d ? reps[j].d : reps[j].fe;
            if (!x || isnan(x[i]))
                continue;
            ST_double delta = x[i] - mu;
            c++;
            mu += delta / c;
            s += delta * (x[i] - mu);
        }
        mean[i] = mu;
        m2[i] = s;
        count[i] = c;
    }
    (void)threads;
}

/* Fold replicates reps[0 .. n) into acc, in replicate order */
static void acc_add_batch(hdfe_boot_acc *acc, const hdfe_boot_rep *reps, ST_int n,
                          ST_int threads)
{
    const ST_int R0 = acc->R;
    for (ST_int j = 0; j < n; j++)
    {
        const hdfe_boot_rep *rep = &reps[j];
        ST_double *row = acc->b + (size_t)acc->R * (size_t)acc->K_b;
        ST_int status = (rep->status >= 0 && rep->status < HDFE_BOOT_NSTATUS) ? rep->status
                                                                             : HDFE_BOOT_OTHER;
        for (ST_int k = 0; k < acc->K_b; k++)
            row[k] = status == HDFE_BOOT_OK ? rep->b[k] : NAN;
        acc->status[acc->R] = status;
        acc->status_count[status]++;
        if (acc->fe_draws && rep->fe &&
            (status == HDFE_BOOT_OK || status == HDFE_BOOT_COLLINEAR))
        {
            for (ST_int i = 0; i < acc->L; i++)
                acc->fe_draws[(size_t)i * (size_t)acc->R_cap + (size_t)acc->R] =
                    (float)rep->fe[i];
        }
        acc->iterations += rep->iterations;
        acc->fe_split_reps += rep->fe_split != 0;
        acc->R++;
    }
    if (acc->L > 0)
        welford_batch(reps, acc->status + R0, n, 0, acc->L, acc->fe_mean, acc->fe_m2, acc->fe_n,
                      threads);
    if (acc->Nd > 0)
        welford_batch(reps, acc->status + R0, n, 1, acc->Nd, acc->d_mean, acc->d_m2, acc->d_n,
                      threads);
}

/* ========================================================================
 * Replicate scheduling
 * ======================================================================== */

/* Memory of one replicate in flight: the estimator's workspace plus the
 * multiplicities and result buffers */
static size_t slot_bytes(const hdfe_boot_design *d, const hdfe_boot_model_ops *ops,
                         const void *model, size_t outputs)
{
    return ops->workspace_bytes(model) + (size_t)d->U * sizeof(ST_int) +
           outputs * sizeof(ST_double);
}

ST_int hdfe_boot_parallel(const hdfe_boot_design *d, const hdfe_boot_model_ops *ops,
                          const void *model, ST_int K_b, ST_int L, ST_int Nd, ST_int threads,
                          size_t budget)
{
    if (!d || !ops || !ops->workspace_bytes || K_b < 0 || L < 0 || Nd < 0)
        return 0;
    const size_t slot = slot_bytes(d, ops, model, (size_t)K_b + (size_t)L + (size_t)Nd);
    ST_int T = threads < 1 ? 1 : threads;
    if (slot > 0 && budget / slot < (size_t)T)
        T = (ST_int)(budget / slot);
    return T;
}

ST_retcode hdfe_boot_run(const hdfe_boot_design *d, const ST_int *supplied,
                         const hdfe_boot_model_ops *ops, const void *model, ST_int first,
                         ST_int count, ST_int threads, size_t budget, hdfe_boot_acc *acc)
{
    if (!d || !ops || !ops->fit || !ops->workspace_bytes || !ops->workspace_new ||
        !ops->workspace_free || !acc || count < 0 || first != acc->R ||
        count > acc->R_cap - acc->R)
        return 198;
    if (count == 0)
        return 0;

    const ST_int all = threads < 1 ? 1 : threads;
    const size_t outputs = (size_t)acc->K_b + (size_t)acc->L + (size_t)acc->Nd;
    const size_t out_bytes = (outputs > 0 ? outputs : 1) * sizeof(ST_double);
    const size_t slot = slot_bytes(d, ops, model, outputs);
    /* T replicates in flight, each with its own workspace */
    ST_int T = hdfe_boot_parallel(d, ops, model, acc->K_b, acc->L, acc->Nd, all, budget);
    if (T < 1)
        return 909;
    if (T > count)
        T = count;

    /* The threads left over go to each replicate's own column parallelism,
     * as far as the budget holds their buffers (plan 6.5) */
    ST_int inner = all / T;
    size_t extra = 0;
    if (inner > 1 && ops->thread_bytes)
    {
        extra = ops->thread_bytes(model);
        const size_t per = budget / (size_t)T;
        if (per <= slot)
            inner = 1;
        else if (extra > 0 && (per - slot) / extra < (size_t)(inner - 1))
            inner = 1 + (ST_int)((per - slot) / extra);
    }
    if (inner < 1)
        inner = 1;

    /* Replicates are taken dynamically from a window of W >= T, all of the
     * chunk when their results fit: a slow replicate then holds up no other
     * thread. Large fixed-effect outputs limit the window to a multiple of
     * T within the budget (and 256 MB). */
    ST_int W = T;
    {
        const size_t used = (size_t)T * (slot + (size_t)(inner - 1) * extra);
        size_t room = budget > used ? budget - used : 0;
        if (room > ((size_t)256 << 20))
            room = (size_t)256 << 20;
        const size_t more = room / out_bytes;
        if (more >= (size_t)(count - T))
            W = count;
        else
            W = T + (ST_int)(more / (size_t)T) * T;
    }

    ST_retcode rc = 0;
    void **ws = calloc((size_t)T, sizeof(*ws));
    ST_int *mult = malloc((size_t)T * (size_t)d->U * sizeof(*mult));
    ST_double *buffers = malloc((size_t)W * out_bytes);
    hdfe_boot_rep *reps = calloc((size_t)W, sizeof(*reps));
    ST_retcode *fit_rc = calloc((size_t)W, sizeof(*fit_rc));
    if (!ws || !mult || !buffers || !reps || !fit_rc)
    {
        rc = 920;
        goto cleanup;
    }
    for (ST_int t = 0; t < T; t++)
    {
        if (!(ws[t] = ops->workspace_new(model)))
        {
            rc = 920;
            goto cleanup;
        }
    }
    for (ST_int j = 0; j < W; j++)
    {
        ST_double *buf = buffers + (size_t)j * outputs;
        reps[j].b = buf;
        reps[j].fe = acc->L > 0 ? buf + acc->K_b : NULL;
        reps[j].d = acc->Nd > 0 ? buf + acc->K_b + acc->L : NULL;
    }

    for (ST_int start = 0; start < count && !rc; start += W)
    {
        if (ctools_cancel_requested()) { rc = 1; break; }
        ctools_cancel_team team = {0};
        const ST_int n = count - start < W ? count - start : W;
        const ST_int nthr = n < T ? n : T;
        (void)nthr;
        for (ST_int j = 0; j < n; j++)
        {
            reps[j].status = HDFE_BOOT_OK;
            reps[j].iterations = 0;
            reps[j].fe_split = 0;
            reps[j].threads = inner;
            fit_rc[j] = 0;
        }
#ifdef _OPENMP
        /* A replicate's own parallel regions nest in this one */
        const int levels = omp_get_max_active_levels();
        if (nthr > 1 && inner > 1 && levels < 2)
            omp_set_max_active_levels(2);
#pragma omp parallel num_threads(nthr) if (nthr > 1)
#endif
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
            /* Regions the fit opens without a thread count get its share */
            omp_set_num_threads(inner);
#endif
            ST_int *m = mult + (size_t)tid * (size_t)d->U;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1) nowait
#endif
            for (ST_int j = 0; j < n; j++)
            {
                if (ctools_cancel_requested()) { fit_rc[j] = 1; continue; }
                if (supplied)
                    memcpy(m, supplied + (size_t)(start + j) * (size_t)d->U,
                           (size_t)d->U * sizeof(*m));
                else
                    hdfe_boot_draw_fast(d, first + start + j, m);
                fit_rc[j] = ops->fit(model, ws[tid], m, &reps[j]);
            }
            ctools_cancel_wait(&team);
        }
#ifdef _OPENMP
        omp_set_max_active_levels(levels);
#endif
        if (ctools_cancel_requested()) { rc = 1; break; }
        /* Fold in replicate order, so the sums do not depend on threads */
        for (ST_int j = 0; j < n && !rc; j++)
            rc = fit_rc[j];
        if (!rc)
            acc_add_batch(acc, reps, n, all);
    }

cleanup:
    if (ws)
    {
        for (ST_int t = 0; t < T; t++)
        {
            if (ws[t])
                ops->workspace_free(ws[t]);
        }
    }
    free(ws);
    free(mult);
    free(buffers);
    free(reps);
    free(fit_rc);
    return rc;
}

/* ========================================================================
 * Summary statistics
 * ======================================================================== */

double hdfe_boot_normal(double x)
{
    return 0.5 * erfc(-x / sqrt(2.0));
}

/* Acklam's rational approximation, refined by one Halley step */
double hdfe_boot_invnormal(double p)
{
    static const double a[6] = {-3.969683028665376e+01, 2.209460984245205e+02,
                                -2.759285104469687e+02, 1.383577518672690e+02,
                                -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[5] = {-5.447609879822406e+01, 1.615858368580409e+02,
                                -1.556989798598866e+02, 6.680131188771972e+01,
                                -1.328068155288572e+01};
    static const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
                                -2.400758277161838e+00, -2.549732539343734e+00,
                                4.374664141464968e+00,  2.938163982698783e+00};
    static const double e[4] = {7.784695709041462e-03, 3.224671290700398e-01,
                                2.445134137142996e+00, 3.754408661907416e+00};
    if (!(p > 0.0 && p < 1.0))
        return p == 0.0 ? -INFINITY : (p == 1.0 ? INFINITY : NAN);
    double x, q, r;
    if (p < 0.02425)
    {
        q = sqrt(-2.0 * log(p));
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((e[0] * q + e[1]) * q + e[2]) * q + e[3]) * q + 1.0);
    }
    else if (p > 1.0 - 0.02425)
    {
        q = sqrt(-2.0 * log1p(-p));
        x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((e[0] * q + e[1]) * q + e[2]) * q + e[3]) * q + 1.0);
    }
    else
    {
        q = p - 0.5;
        r = q * q;
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    }
    double err = hdfe_boot_normal(x) - p;
    double u = err * sqrt(2.0 * 3.14159265358979323846) * exp(0.5 * x * x);
    return x - u / (1.0 + 0.5 * x * u);
}

double hdfe_boot_pctile(const double *sorted, ST_int n, double p)
{
    if (n < 1 || !(p >= 0.0 && p <= 100.0))
        return NAN;
    /* P = n p / 100, i = floor(P): x(i + 1) when P is not an integer,
     * otherwise the mean of x(i) and x(i + 1), with x(0) = x(1) and
     * x(n + 1) = x(n) (1-based order statistics) */
    double P = (double)n * p / 100.0;
    double i = floor(P);
    ST_int lo = (ST_int)i;
    if (P > i)
        return sorted[lo < n ? lo : n - 1];
    ST_int a = lo < 1 ? 0 : lo - 1, b = lo < n ? lo : n - 1;
    return 0.5 * (sorted[a] + sorted[b]);
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

void hdfe_boot_summary_free(hdfe_boot_summary *s)
{
    if (!s)
        return;
    free(s->b_bs);
    free(s->se);
    free(s->bias);
    free(s->z0);
    free(s->V);
    free(s->ci_normal);
    free(s->ci_percentile);
    free(s->ci_bc);
    memset(s, 0, sizeof(*s));
}

ST_retcode hdfe_boot_summarize(const hdfe_boot_acc *acc, const ST_double *b_hat, int mse,
                               ST_double level, hdfe_boot_summary *out)
{
    memset(out, 0, sizeof(*out));
    if (!acc || !b_hat || !(level > 0.0 && level < 100.0))
        return 198;
    const ST_int K = acc->K_b, R = acc->R;
    out->K_b = K;

    /* Complete replicates (casewise) */
    ST_int n = 0;
    ST_int *rows = malloc((size_t)(R > 0 ? R : 1) * sizeof(*rows));
    if (!rows)
        return 920;
    for (ST_int r = 0; r < R; r++)
    {
        const ST_double *row = acc->b + (size_t)r * (size_t)K;
        int complete = acc->status[r] == HDFE_BOOT_OK;
        for (ST_int k = 0; k < K && complete; k++)
            complete = isfinite(row[k]);
        if (complete)
            rows[n++] = r;
    }
    out->N_reps = n;
    out->N_misreps = R - n;
    if (n < 2)
    {
        free(rows);
        return 2000;
    }

    size_t KK = (size_t)(K > 0 ? K : 1);
    out->b_bs = calloc(KK, sizeof(double));
    out->se = calloc(KK, sizeof(double));
    out->bias = calloc(KK, sizeof(double));
    out->z0 = calloc(KK, sizeof(double));
    out->V = calloc(KK * KK, sizeof(double));
    out->ci_normal = calloc(2 * KK, sizeof(double));
    out->ci_percentile = calloc(2 * KK, sizeof(double));
    out->ci_bc = calloc(2 * KK, sizeof(double));
    double *x = malloc((size_t)n * sizeof(double));
    if (!out->b_bs || !out->se || !out->bias || !out->z0 || !out->V || !out->ci_normal ||
        !out->ci_percentile || !out->ci_bc || !x)
    {
        free(rows);
        free(x);
        hdfe_boot_summary_free(out);
        return 920;
    }

#define REP(i, k) acc->b[(size_t)rows[i] * (size_t)K + (size_t)(k)]
    for (ST_int k = 0; k < K; k++)
    {
        double s = 0.0;
        for (ST_int i = 0; i < n; i++)
            s += REP(i, k);
        out->b_bs[k] = s / n;
    }
    /* V about the replicate mean (n - 1), or with mse about b_hat (n) */
    for (ST_int j = 0; j < K; j++)
    {
        for (ST_int k = 0; k <= j; k++)
        {
            double cj = mse ? b_hat[j] : out->b_bs[j];
            double ck = mse ? b_hat[k] : out->b_bs[k];
            double s = 0.0;
            for (ST_int i = 0; i < n; i++)
                s += (REP(i, j) - cj) * (REP(i, k) - ck);
            s /= mse ? (double)n : (double)(n - 1);
            out->V[(size_t)j * KK + (size_t)k] = out->V[(size_t)k * KK + (size_t)j] = s;
        }
    }

    const double zalpha = hdfe_boot_invnormal((100.0 + level) / 200.0);
    for (ST_int k = 0; k < K; k++)
    {
        /* se and the normal interval use the standard deviation about the
         * replicate mean, also under mse, as _bs_sum does */
        double s = 0.0;
        for (ST_int i = 0; i < n; i++)
        {
            double dev = REP(i, k) - out->b_bs[k];
            s += dev * dev;
            x[i] = REP(i, k);
        }
        double sd = sqrt(s / (n - 1));
        out->se[k] = sd;
        out->bias[k] = out->b_bs[k] - b_hat[k];
        out->ci_normal[k] = b_hat[k] - zalpha * sd;
        out->ci_normal[KK + (size_t)k] = b_hat[k] + zalpha * sd;

        qsort(x, (size_t)n, sizeof(double), cmp_double);
        out->ci_percentile[k] = hdfe_boot_pctile(x, n, (100.0 - level) / 2.0);
        out->ci_percentile[KK + (size_t)k] = hdfe_boot_pctile(x, n, (100.0 + level) / 2.0);

        /* Bias-corrected: z0 from #{b^r <= b_hat}, missing at 0 or n */
        ST_int below = 0;
        for (ST_int i = 0; i < n; i++)
            below += x[i] <= b_hat[k];
        out->z0[k] = out->ci_bc[k] = out->ci_bc[KK + (size_t)k] = NAN;
        if (below > 0 && below < n)
        {
            double z0 = hdfe_boot_invnormal((double)below / n);
            out->z0[k] = z0;
            out->ci_bc[k] = hdfe_boot_pctile(x, n, 100.0 * hdfe_boot_normal(2.0 * z0 - zalpha));
            out->ci_bc[KK + (size_t)k] =
                hdfe_boot_pctile(x, n, 100.0 * hdfe_boot_normal(2.0 * z0 + zalpha));
        }
    }
#undef REP
    free(rows);
    free(x);
    return 0;
}

void hdfe_boot_fe_se(const hdfe_boot_acc *acc, ST_double *fe_se, ST_double *d_se)
{
    if (fe_se)
    {
        for (ST_int l = 0; l < acc->L; l++)
            fe_se[l] = acc->fe_n[l] >= 2 ? sqrt(acc->fe_m2[l] / (acc->fe_n[l] - 1)) : NAN;
    }
    if (d_se)
    {
        for (ST_int i = 0; i < acc->Nd; i++)
            d_se[i] = acc->d_n[i] >= 2 ? sqrt(acc->d_m2[i] / (acc->d_n[i] - 1)) : NAN;
    }
}

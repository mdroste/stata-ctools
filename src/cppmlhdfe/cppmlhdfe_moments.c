#include "cppmlhdfe_internal.h"

/* Transform the intercept of a centered sandwich, V <- T V T'. These
 * differences can cancel to zero even when the slopes are well conditioned.
 * Keep product and sum roundoff until the final conversion to double. */
static dd_real ppml_vce_add_product(dd_real sum, double a, double b)
{
    dd_real product = two_prod(a, b);
    return dd_add_d(dd_add_d(sum, product.hi), product.lo);
}

void ppml_uncenter_vce(double *V, const double *means, int K)
{
    const int stride = K + 1;
    dd_real corner = {V[K * stride + K], 0};
    for (int i = 0; i < K; i++)
    {
        corner = ppml_vce_add_product(corner, -2 * means[i], V[i * stride + K]);
        for (int j = 0; j < K; j++)
        {
            dd_real product = two_prod(means[i], means[j]);
            corner = ppml_vce_add_product(corner, product.hi, V[i * stride + j]);
            corner = ppml_vce_add_product(corner, product.lo, V[i * stride + j]);
        }
    }
    for (int i = 0; i < K; i++)
    {
        dd_real side = {V[i * stride + K], 0};
        for (int j = 0; j < K; j++)
            side = ppml_vce_add_product(side, -V[i * stride + j], means[j]);
        V[i * stride + K] = V[K * stride + i] = side.hi + side.lo;
    }
    V[K * stride + K] = corner.hi + corner.lo;
}

/* The sum of x over fixed blocks of rows, each summed in order and the blocks
 * then added in order: the same value at every thread count. */
#define PPML_SUM_BLOCK 32768
static ST_double ppml_block_sum(const ST_double *x, ST_int N)
{
    const ST_int nblocks = N / PPML_SUM_BLOCK + (N % PPML_SUM_BLOCK != 0);
    ST_double total = 0.0;
    ST_double *part =
        nblocks > 8 ? (ST_double *)malloc((size_t)nblocks * sizeof(ST_double)) : NULL;
    ST_int b;
    /* without the scratch the blocks run on this thread and add up directly */
    #pragma omp parallel for schedule(static) if (part != NULL)
    for (b = 0; b < nblocks; b++)
    {
        const ST_int lo = b * PPML_SUM_BLOCK;
        const ST_int hi = N - lo > PPML_SUM_BLOCK ? lo + PPML_SUM_BLOCK : N;
        ST_double sum = 0.0;
        for (ST_int row = lo; row < hi; row++)
            sum += x[row];
        if (part)
            part[b] = sum;
        else
            total += sum;
    }
    if (part)
    {
        for (b = 0; b < nblocks; b++)
            total += part[b];
        free(part);
    }
    return total;
}

void ppml_update_fe_weights(HDFE_State *S, const ST_double *irls_w, ST_int N)
{
    S->sum_weights = ppml_block_sum(irls_w, N);
    ppml_update_fe_counts(S, irls_w, N);
}

void ppml_update_fe_counts(HDFE_State *S, const ST_double *irls_w, ST_int N)
{
    ST_int g, i, lev;

    /* each factor's counts are an independent pass over the rows */
#pragma omp parallel for schedule(static) private(i, lev) if (S->G > 1 && N > 100000)
    for (g = 0; g < S->G; g++)
    {
        FE_Factor *f = &S->factors[g];
        if (f->slope_center)
        {
            memset(f->weighted_counts, 0, f->num_levels * sizeof(ST_double));
            memset(f->slope_center, 0, f->num_levels * sizeof(ST_double));
            for (i = 0; i < N; i++)
            {
                lev = f->levels[i] - 1;
                f->weighted_counts[lev] += irls_w[i];
                f->slope_center[lev] += irls_w[i] * f->slope[i];
            }
            for (i = 0; i < f->num_levels; i++)
                if (f->weighted_counts[i] > 0)
                    f->slope_center[i] /= f->weighted_counts[i];
        }
        memset(f->weighted_counts, 0, f->num_levels * sizeof(ST_double));
        for (i = 0; i < N; i++)
        {
            lev = f->levels[i] - 1;
            ST_double loading =
                f->slope ? f->slope[i] - (f->slope_center ? f->slope_center[lev] : 0.0) : 1.0;
            f->weighted_counts[lev] += irls_w[i] * loading * loading;
        }
        if (f->inv_weighted_counts)
        {
            for (i = 0; i < f->num_levels; i++)
            {
                f->inv_weighted_counts[i] =
                    (f->weighted_counts[i] > 0.0) ? 1.0 / f->weighted_counts[i] : 0.0;
            }
        }
    }
}

/* y*log(y/mu) - (y - mu) for mu = exp(eta_total). When y is close to a
   large mean the two terms nearly cancel and y*(log(y) - eta) keeps only
   about eps*|log y|/|y/mu - 1| of its digits; mu*((1+r)*log1p(r) - r) with
   r = y/mu - 1 is exact to rounding. y = 0 gives mu. The part for y > 0
   can remember log(y), which one branch reads: *log_y is NaN until that
   branch has evaluated it and holds it afterwards (NULL: nothing kept). */
static inline ST_double ppml_unit_deviance_pos(ST_double y, ST_double mu, ST_double eta_total,
                                               ST_double *log_y)
{
    ST_double r = (y - mu) / mu;
    if (fabs(r) < 1e-2) {
        /* (1+r)log1p(r) - r = sum_{k>=2} (-1)^k r^k / (k(k-1)) */
        ST_double term = r * r, sum = 0.0;
        for (int k = 2; k <= 12; k++) {
            sum += ((k % 2) ? -term : term) / (ST_double)(k * (k - 1));
            term *= r;
        }
        return mu * sum;
    }
    if (fabs(r) < 1.0) return mu * ((1.0 + r) * log1p(r) - r);
    if (!log_y) return y * (log(y) - eta_total) - (y - mu);
    /* NaN: not evaluated yet */
    if (*log_y != *log_y) *log_y = log(y);
    return y * (*log_y - eta_total) - (y - mu);
}

static ST_double ppml_unit_deviance_mu(ST_double y, ST_double mu, ST_double eta_total)
{
    if (y <= 0.0) return mu;
    return ppml_unit_deviance_pos(y, mu, eta_total, NULL);
}

static ST_double ppml_unit_deviance(ST_double y, ST_double eta_total)
{
    return ppml_unit_deviance_mu(y, exp(eta_total), eta_total);
}

/* y*log(y) - y - lngamma(y+1), the saturated Poisson log likelihood of one
   outcome; Stirling's series avoids the cancellation of its large terms */
static ST_double ppml_saturated_loglik(ST_double y)
{
    if (y <= 0.0) return 0.0;
    if (y < 10.0) return y * log(y) - y - lgamma(y + 1.0);
    ST_double inv = 1.0 / y, inv2 = inv * inv;
    return -0.5 * (1.8378770664093454836 + log(y))   /* log(2*pi) */
           - inv * (1.0 / 12.0 - inv2 * (1.0 / 360.0 - inv2 * (1.0 / 1260.0 - inv2 / 1680.0)));
}

/* The deviance and log likelihood are double-double sums over fixed blocks
 * of rows, each summed in row order, and then the block sums in block order.
 * The blocks run in parallel, and since they do not depend on the thread
 * count neither does the result. */
#define PPML_SUM_BLOCK 32768

static dd_real ppml_dd_add_dd(dd_real a, dd_real b)
{
    return dd_add_d(dd_add_d(a, b.hi), b.lo);
}

static dd_real ppml_deviance_block(const ST_double *y, const ST_double *eta, const ST_double *offset,
                                   const ST_double *w_user, ST_int lo, ST_int hi)
{
    dd_real dev = {0.0, 0.0};
    for (ST_int i = lo; i < hi; i++) {
        ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        ST_double eta_total = eta[i] + (offset != NULL ? offset[i] : 0.0);
        dev = dd_add_d(dev, w * ppml_unit_deviance(y[i], eta_total));
    }
    return dev;
}

static dd_real ppml_loglik_block(const ST_double *y, const ST_double *eta, const ST_double *offset,
                                 ST_double log_scale, const ST_double *w_user, ST_int lo, ST_int hi)
{
    dd_real ll = {0.0, 0.0};
    for (ST_int i = lo; i < hi; i++) {
        ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        ST_double eta_total = eta[i] + (offset != NULL ? offset[i] : 0.0) + log_scale;
        ll = dd_add_d(ll, w * (ppml_saturated_loglik(y[i]) - ppml_unit_deviance(y[i], eta_total)));
    }
    return ll;
}

ST_double ppml_compute_deviance(const ST_double *y, const ST_double *eta, const ST_double *offset,
                                const ST_double *w_user, ST_int N)
{
    ST_int nblocks = (N + PPML_SUM_BLOCK - 1) / PPML_SUM_BLOCK;
    dd_real *part = nblocks > 1 ? (dd_real *)malloc((size_t)nblocks * sizeof(dd_real)) : NULL;
    dd_real dev = {0.0, 0.0};
    if (!part) {
        for (ST_int b = 0; b < nblocks; b++) {
            ST_int lo = b * PPML_SUM_BLOCK, hi = lo + PPML_SUM_BLOCK < N ? lo + PPML_SUM_BLOCK : N;
            dev = ppml_dd_add_dd(dev, ppml_deviance_block(y, eta, offset, w_user, lo, hi));
        }
        return 2.0 * (dev.hi + dev.lo);
    }
    #pragma omp parallel for schedule(static)
    for (ST_int b = 0; b < nblocks; b++) {
        ST_int lo = b * PPML_SUM_BLOCK, hi = lo + PPML_SUM_BLOCK < N ? lo + PPML_SUM_BLOCK : N;
        part[b] = ppml_deviance_block(y, eta, offset, w_user, lo, hi);
    }
    for (ST_int b = 0; b < nblocks; b++) dev = ppml_dd_add_dd(dev, part[b]);
    free(part);
    return 2.0 * (dev.hi + dev.lo);
}

/* ------------------------------------------------------------------------
 * The passes of one IRLS update. Each replaces several passes over the rows
 * that computed the same numbers one array at a time: the arithmetic of every
 * row, and the order of every sum, are those of the separate passes.
 * ---------------------------------------------------------------------- */

/* Rows per block of the linear-predictor pass: its K + 3 arrays stay in the
 * cache while the regressors are added one after another */
#define PPML_STEP_BLOCK 4096

void ppml_irls_step_eta(const ST_double *X_part, ST_int N, ST_int K, const ST_double *beta,
                        const ST_double *z, const ST_double *z_part, const ST_double *offset,
                        ST_double *eta_new, int *eta_not_finite, int *mu_overflow)
{
    const ST_int nblocks = N / PPML_STEP_BLOCK + (N % PPML_STEP_BLOCK != 0);
    int bad = 0, over = 0;
    ST_int b;
    #pragma omp parallel for schedule(static) reduction(| : bad, over) if (N > 100000)
    for (b = 0; b < nblocks; b++) {
        const ST_int lo = b * PPML_STEP_BLOCK;
        const ST_int hi = N - lo > PPML_STEP_BLOCK ? lo + PPML_STEP_BLOCK : N;
        /* X~ b, each row accumulated over the regressors in order */
        for (ST_int i = lo; i < hi; i++) eta_new[i] = 0.0;
        for (ST_int k = 0; k < K; k++) {
            const ST_double *x = X_part + (size_t)k * N;
            const ST_double bk = beta[k];
            for (ST_int i = lo; i < hi; i++) eta_new[i] += x[i] * bk;
        }
        for (ST_int i = lo; i < hi; i++) {
            eta_new[i] += z[i] - z_part[i];
            if (!isfinite(eta_new[i])) bad |= 1;
            else {
                /* exp() overflows only above log(DBL_MAX) = 709.78; below 709
                 * it is finite, so only those rare rows need the call */
                const ST_double eta_off = eta_new[i] + (offset ? offset[i] : 0.0);
                if (eta_off > 709.0 && !isfinite(exp(eta_off))) over |= 1;
            }
        }
    }
    *eta_not_finite |= bad;
    *mu_overflow |= over;
}

/* The row loop of ppml_irls_step_means is a function of its own that starts
 * on a 4 KB page boundary. Measured on a Ryzen 9 9950X (2026-11-26): where
 * this loop sits within its page decides whether a hardware thread running
 * it slows the other hardware thread of its core. At one position in eight
 * of the loop as it was (inlined in the parallel region, with a call of fmax
 * and a branch on the outcome after the sum), rows that call log1p cost the
 * sibling thread 30% more; one thread per core, or rows without log1p,
 * showed nothing. The loop sat wherever the link put it, after all the code
 * before it, so the step's speed changed with any build of unrelated code.
 * Aligned to a page it is at the same place in every build. In this form
 * the loop is also less exposed: of 64 places tried, the worst two cost the
 * sibling thread 8%. The once-per-fit passes in this file were timed the
 * same way and are left where the link puts them: they have at most a few
 * places that cost 5 to 9%, and aligning them showed no gain in a plugin. */
#if defined(__GNUC__) && !defined(__clang__)
#define PPML_ONE_PAGE __attribute__((noinline, noclone, aligned(4096)))
#elif defined(__clang__)
#define PPML_ONE_PAGE __attribute__((noinline, aligned(4096)))
#else
#define PPML_ONE_PAGE
#endif

#define PPML_ROWS_NOT_FINITE 1
#define PPML_ROWS_OVERFLOW 2
#define PPML_ROWS_FLOORED 4

/* Rows lo..hi-1 of ppml_irls_step_means: returns their deviance sum (before
 * the factor 2), the largest change of eta in *change and PPML_ROWS_* bits
 * in *flags. */
PPML_ONE_PAGE
static dd_real ppml_step_means_rows(const ST_double *y, const ST_double *eta_new,
                                    const ST_double *offset, const ST_double *w_user,
                                    ST_int lo, ST_int hi, ST_double *eta, ST_double *mu,
                                    ST_double *change, int *flags)
{
    dd_real dev = {0.0, 0.0};
    ST_double max_change = 0.0;
    int bad = 0, over = 0, floored = 0;
    for (ST_int i = lo; i < hi; i++) {
        const ST_double z_hat = eta_new[i];
        const ST_double off = (offset != NULL) ? offset[i] : 0.0;
        if (!isfinite(z_hat)) { bad = 1; continue; }
        const ST_double yi = y[i];
        /* Zero outcomes with negligible means may have no finite FE
         * coefficient after a tolerance-limited separation check. They
         * must not prevent convergence of the identified parameters.
         * Compares instead of a call of fmax and of branches on the
         * outcome: the largest change is the same number (a NaN difference
         * is skipped as fmax skips it, and no change is below +0). */
        const ST_double moved = ((yi > 0) | (mu[i] > 1e-12)) ? fabs(z_hat - eta[i]) : 0.0;
        if (moved > max_change) max_change = moved;
        eta[i] = z_hat;
        const ST_double eta_total = z_hat + off;
        ST_double m = exp(eta_total);
        /* the deviance takes the mean before the floor */
        const ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        dev = dd_add_d(dev, w * ppml_unit_deviance_mu(yi, m, eta_total));
        /* As ppmlhdfe: floor mu for the weights only, never cap it (a
           capped mean no longer agrees with eta and can stall the
           deviance at a wrong point); an overflowing mean stops the fit */
        over |= !isfinite(m);
        if (m < PPML_MU_FLOOR) m = PPML_MU_FLOOR;
        mu[i] = m;
        floored |= (yi > 0) & (m <= PPML_MU_FLOOR);
    }
    *change = max_change;
    *flags = (bad ? PPML_ROWS_NOT_FINITE : 0) | (over ? PPML_ROWS_OVERFLOW : 0) |
             (floored ? PPML_ROWS_FLOORED : 0);
    return dev;
}

ST_retcode ppml_irls_step_means(const ST_double *y, const ST_double *eta_new,
                                const ST_double *offset, const ST_double *w_user, ST_int N,
                                ST_double *eta, ST_double *mu, ppml_step_means *out)
{
    const ST_int nblocks = (N + PPML_SUM_BLOCK - 1) / PPML_SUM_BLOCK;
    dd_real one = {0.0, 0.0};
    dd_real *part = nblocks > 1 ? (dd_real *)malloc((size_t)nblocks * sizeof(dd_real)) : &one;
    if (!part) return 920;
    ST_double max_change = 0.0;
    int flags = 0;
    ST_int b;
    #pragma omp parallel for schedule(static) reduction(max : max_change)                          \
        reduction(| : flags) if (N > 100000)
    for (b = 0; b < nblocks; b++) {
        const ST_int lo = b * PPML_SUM_BLOCK;
        const ST_int hi = lo + PPML_SUM_BLOCK < N ? lo + PPML_SUM_BLOCK : N;
        ST_double block_change = 0.0;
        int block_flags = 0;
        part[b] = ppml_step_means_rows(y, eta_new, offset, w_user, lo, hi, eta, mu,
                                       &block_change, &block_flags);
        if (block_change > max_change) max_change = block_change;
        flags |= block_flags;
    }
    dd_real total = {0.0, 0.0};
    for (b = 0; b < nblocks; b++) total = ppml_dd_add_dd(total, part[b]);
    if (part != &one) free(part);
    out->max_eta_change = max_change;
    out->eta_not_finite = (flags & PPML_ROWS_NOT_FINITE) != 0;
    out->mu_overflow = (flags & PPML_ROWS_OVERFLOW) != 0;
    out->floored_positive = (flags & PPML_ROWS_FLOORED) != 0;
    out->deviance = 2.0 * (total.hi + total.lo);
    return 0;
}

ST_retcode ppml_irls_step_weights(const ST_double *mu, const ST_double *w_user,
                                  ST_double user_weight_scale, const ST_int *excluded, ST_int N,
                                  ST_double *irls_w, ST_double *sum_out)
{
    const ST_int nblocks = N / PPML_SUM_BLOCK + (N % PPML_SUM_BLOCK != 0);
    ST_double one = 0.0;
    ST_double *part = nblocks > 1 ? (ST_double *)malloc((size_t)nblocks * sizeof(ST_double)) : &one;
    if (!part) return 920;
    ST_int b;
    #pragma omp parallel for schedule(static) if (N > 100000)
    for (b = 0; b < nblocks; b++) {
        const ST_int lo = b * PPML_SUM_BLOCK;
        const ST_int hi = N - lo > PPML_SUM_BLOCK ? lo + PPML_SUM_BLOCK : N;
        ST_double sum = 0.0;
        for (ST_int i = lo; i < hi; i++) {
            const ST_double w_u = (w_user != NULL) ? w_user[i] : 1.0;
            ST_double w = mu[i] * (w_u / user_weight_scale);
            if (w < DBL_MIN) w = DBL_MIN;
            /* A row separation(mu) has flagged leaves the fit: its mean is
             * held at a marker of 2^-46 and the fit is redone without it
             * once this pass ends. Its weight here is exactly zero rather
             * than that marker's: a positive weight sixteen orders below the
             * others is what a projection can neither use nor certify (the
             * iteration on the rows once reported a breakdown for it), while
             * a zero one simply takes the row out. */
            if (excluded && excluded[i]) w = 0.0;
            irls_w[i] = w;
            sum += w;
        }
        part[b] = sum;
    }
    ST_double total = 0.0;
    for (b = 0; b < nblocks; b++) total += part[b];
    if (part != &one) free(part);
    *sum_out = total;
    return 0;
}

/* ll = sum of y*log(mu) - mu - lngamma(y+1) = saturated ll - deviance/2 */
ST_double ppml_compute_loglik(const ST_double *y, const ST_double *eta, const ST_double *offset,
                              ST_double log_scale, const ST_double *w_user, ST_int N)
{
    ST_int nblocks = (N + PPML_SUM_BLOCK - 1) / PPML_SUM_BLOCK;
    dd_real *part = nblocks > 1 ? (dd_real *)malloc((size_t)nblocks * sizeof(dd_real)) : NULL;
    dd_real ll = {0.0, 0.0};
    if (!part) {
        for (ST_int b = 0; b < nblocks; b++) {
            ST_int lo = b * PPML_SUM_BLOCK, hi = lo + PPML_SUM_BLOCK < N ? lo + PPML_SUM_BLOCK : N;
            ll = ppml_dd_add_dd(ll, ppml_loglik_block(y, eta, offset, log_scale, w_user, lo, hi));
        }
        return ll.hi + ll.lo;
    }
    #pragma omp parallel for schedule(static)
    for (ST_int b = 0; b < nblocks; b++) {
        ST_int lo = b * PPML_SUM_BLOCK, hi = lo + PPML_SUM_BLOCK < N ? lo + PPML_SUM_BLOCK : N;
        part[b] = ppml_loglik_block(y, eta, offset, log_scale, w_user, lo, hi);
    }
    for (ST_int b = 0; b < nblocks; b++) ll = ppml_dd_add_dd(ll, part[b]);
    free(part);
    return ll.hi + ll.lo;
}

/* The terms of ppml_compute_loglik for eta and, in the same pass, for the
 * constant eta0 (the null model, no offset or scale), each summed as
 * ppml_compute_loglik sums them: the two results are those of two separate
 * calls, while the saturated term and exp(eta0) are evaluated once */
static void ppml_loglik_pair_block(const ST_double *y, const ST_double *eta, const ST_double *offset,
                                   ST_double log_scale, const ST_double *w_user, ST_double eta0,
                                   ST_double mu0, ST_int lo, ST_int hi, dd_real *ll, dd_real *ll0)
{
    dd_real a = {0.0, 0.0}, b = {0.0, 0.0};
    for (ST_int i = lo; i < hi; i++) {
        ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        ST_double eta_total = eta[i] + (offset != NULL ? offset[i] : 0.0) + log_scale;
        ST_double saturated = ppml_saturated_loglik(y[i]);
        a = dd_add_d(a, w * (saturated - ppml_unit_deviance(y[i], eta_total)));
        b = dd_add_d(b, w * (saturated - ppml_unit_deviance_mu(y[i], mu0, eta0)));
    }
    *ll = a;
    *ll0 = b;
}

void ppml_compute_loglik_pair(const ST_double *y, const ST_double *eta, const ST_double *offset,
                              ST_double log_scale, const ST_double *w_user, ST_double eta0, ST_int N,
                              ST_double *ll_out, ST_double *ll0_out)
{
    /* the null model's eta_total as ppml_compute_loglik forms it */
    eta0 = eta0 + 0.0 + 0.0;
    const ST_double mu0 = exp(eta0);
    ST_int nblocks = (N + PPML_SUM_BLOCK - 1) / PPML_SUM_BLOCK;
    dd_real *part = nblocks > 1 ? (dd_real *)malloc((size_t)nblocks * 2 * sizeof(dd_real)) : NULL;
    dd_real ll = {0.0, 0.0}, ll0 = {0.0, 0.0};
    if (!part) {
        for (ST_int b = 0; b < nblocks; b++) {
            ST_int lo = b * PPML_SUM_BLOCK, hi = lo + PPML_SUM_BLOCK < N ? lo + PPML_SUM_BLOCK : N;
            dd_real a, c;
            ppml_loglik_pair_block(y, eta, offset, log_scale, w_user, eta0, mu0, lo, hi, &a, &c);
            ll = ppml_dd_add_dd(ll, a);
            ll0 = ppml_dd_add_dd(ll0, c);
        }
    } else {
        #pragma omp parallel for schedule(static)
        for (ST_int b = 0; b < nblocks; b++) {
            ST_int lo = b * PPML_SUM_BLOCK, hi = lo + PPML_SUM_BLOCK < N ? lo + PPML_SUM_BLOCK : N;
            ppml_loglik_pair_block(y, eta, offset, log_scale, w_user, eta0, mu0, lo, hi,
                                   &part[2 * (size_t)b], &part[2 * (size_t)b + 1]);
        }
        for (ST_int b = 0; b < nblocks; b++) {
            ll = ppml_dd_add_dd(ll, part[2 * (size_t)b]);
            ll0 = ppml_dd_add_dd(ll0, part[2 * (size_t)b + 1]);
        }
        free(part);
    }
    *ll_out = ll.hi + ll.lo;
    *ll0_out = ll0.hi + ll0.lo;
}

/* ------------------------------------------------------------------------
 * Passes of the steps after the IRLS loop
 * ---------------------------------------------------------------------- */

/* One block of ppml_means_deviance: ppml_deviance_block, which also stores
 * every mean it evaluates (floored for mu, as it is for the weights) */
static dd_real ppml_means_deviance_block(const ST_double *y, const ST_double *eta,
                                         const ST_double *offset, const ST_double *w_user,
                                         ST_int lo, ST_int hi, ST_double *mu)
{
    dd_real dev = {0.0, 0.0};
    for (ST_int i = lo; i < hi; i++) {
        ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        ST_double eta_total = eta[i] + (offset != NULL ? offset[i] : 0.0);
        const ST_double m = exp(eta_total);
        dev = dd_add_d(dev, w * ppml_unit_deviance_mu(y[i], m, eta_total));
        mu[i] = m < PPML_MU_FLOOR ? PPML_MU_FLOOR : m;
    }
    return dev;
}

ST_double ppml_means_deviance(const ST_double *y, const ST_double *eta, const ST_double *offset,
                              const ST_double *w_user, ST_int N, ST_double *mu)
{
    ST_int nblocks = (N + PPML_SUM_BLOCK - 1) / PPML_SUM_BLOCK;
    dd_real *part = nblocks > 1 ? (dd_real *)malloc((size_t)nblocks * sizeof(dd_real)) : NULL;
    dd_real dev = {0.0, 0.0};
    if (!part) {
        for (ST_int b = 0; b < nblocks; b++) {
            ST_int lo = b * PPML_SUM_BLOCK, hi = lo + PPML_SUM_BLOCK < N ? lo + PPML_SUM_BLOCK : N;
            dev = ppml_dd_add_dd(dev, ppml_means_deviance_block(y, eta, offset, w_user, lo, hi, mu));
        }
        return 2.0 * (dev.hi + dev.lo);
    }
    #pragma omp parallel for schedule(static)
    for (ST_int b = 0; b < nblocks; b++) {
        ST_int lo = b * PPML_SUM_BLOCK, hi = lo + PPML_SUM_BLOCK < N ? lo + PPML_SUM_BLOCK : N;
        part[b] = ppml_means_deviance_block(y, eta, offset, w_user, lo, hi, mu);
    }
    for (ST_int b = 0; b < nblocks; b++) dev = ppml_dd_add_dd(dev, part[b]);
    free(part);
    return 2.0 * (dev.hi + dev.lo);
}

/* Rows per block of ppml_fixed_effect_part: its K + 3 arrays stay in the
 * cache while the regressors are added one after another */
#define PPML_PART_BLOCK 4096

void ppml_fixed_effect_part(const ST_double *X, ST_int N, ST_int K, const ST_double *beta,
                            const ST_double *eta, ST_double *xb, ST_double *fe)
{
    const ST_int nblocks = N / PPML_PART_BLOCK + (N % PPML_PART_BLOCK != 0);
    ST_int b;
    #pragma omp parallel for schedule(static) if (N > 100000)
    for (b = 0; b < nblocks; b++) {
        const ST_int lo = b * PPML_PART_BLOCK;
        const ST_int hi = N - lo > PPML_PART_BLOCK ? lo + PPML_PART_BLOCK : N;
        for (ST_int i = lo; i < hi; i++) xb[i] = 0.0;
        for (ST_int k = 0; k < K; k++) {
            const ST_double *x = X + (size_t)k * N;
            const ST_double bk = beta[k];
            for (ST_int i = lo; i < hi; i++) xb[i] += x[i] * bk;
        }
        for (ST_int i = lo; i < hi; i++) fe[i] = eta[i] - xb[i];
    }
}

dd_real ppml_block_sum_dd(const ST_double *x, ST_int N)
{
    const ST_int nblocks = (N + PPML_SUM_BLOCK - 1) / PPML_SUM_BLOCK;
    dd_real *part = nblocks > 8 ? (dd_real *)malloc((size_t)nblocks * sizeof(dd_real)) : NULL;
    dd_real total = {0.0, 0.0};
    ST_int b;
    /* without the scratch the blocks run on this thread and add up directly */
    #pragma omp parallel for schedule(static) if (part != NULL)
    for (b = 0; b < nblocks; b++) {
        const ST_int lo = b * PPML_SUM_BLOCK;
        const ST_int hi = lo + PPML_SUM_BLOCK < N ? lo + PPML_SUM_BLOCK : N;
        dd_real acc = {0.0, 0.0};
        for (ST_int row = lo; row < hi; row++) acc = dd_add_d(acc, x[row]);
        if (part) part[b] = acc;
        else total = ppml_dd_add_dd(total, acc);
    }
    if (part) {
        for (b = 0; b < nblocks; b++) total = ppml_dd_add_dd(total, part[b]);
        free(part);
    }
    return total;
}

/* What a row's two terms of ppml_compute_loglik_pair_scaled need from its
 * outcome alone: the saturated term (a logarithm, and a log-gamma below 10),
 * the null model's whole term, and log(y) for the model's deviance far from
 * the mean. Outcomes repeat (counts take few values), so a block of rows
 * remembers them by value in a small table and evaluates them once per
 * distinct value instead of once per row. What the table returns is the
 * number the evaluation returned for that very value, so every row adds
 * what it adds without the table. */
#define PPML_Y_MEMO 512 /* slots, a power of two */
typedef struct
{
    unsigned long long key; /* the bits of y; 0 (y = +0 never asks) when empty */
    ST_double saturated;    /* ppml_saturated_loglik(y) */
    ST_double null_term;    /* saturated less the null model's unit deviance */
    ST_double log_y;        /* log(y), NaN until the deviance has needed it */
} ppml_y_memo;

/* ppml_loglik_pair_block for the outcomes y[i] * y_scale */
static void ppml_loglik_pair_scaled_block(const ST_double *y, ST_double y_scale,
                                          const ST_double *eta, const ST_double *offset,
                                          ST_double log_scale, const ST_double *w_user,
                                          ST_double eta0, ST_double mu0, ST_int lo, ST_int hi,
                                          dd_real *ll, dd_real *ll0)
{
    dd_real a = {0.0, 0.0}, b = {0.0, 0.0};
    ppml_y_memo memo[PPML_Y_MEMO];
    for (int slot = 0; slot < PPML_Y_MEMO; slot++) memo[slot].key = 0;
    for (ST_int i = lo; i < hi; i++) {
        const ST_double yi = y[i] * y_scale;
        ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        ST_double eta_total = eta[i] + (offset != NULL ? offset[i] : 0.0) + log_scale;
        const ST_double mu = exp(eta_total);
        if (yi <= 0.0) {
            /* no saturated term; the unit deviances are the means */
            a = dd_add_d(a, w * (0.0 - mu));
            b = dd_add_d(b, w * (0.0 - mu0));
            continue;
        }
        unsigned long long key;
        memcpy(&key, &yi, sizeof(key));
        /* the high bits of the product mix every bit of the value */
        ppml_y_memo *m = &memo[(key * 0x9E3779B97F4A7C15ULL) >> 55];
        if (m->key != key) {
            m->key = key;
            m->saturated = ppml_saturated_loglik(yi);
            m->null_term = m->saturated - ppml_unit_deviance_mu(yi, mu0, eta0);
            m->log_y = NAN;
        }
        a = dd_add_d(a, w * (m->saturated - ppml_unit_deviance_pos(yi, mu, eta_total, &m->log_y)));
        b = dd_add_d(b, w * m->null_term);
    }
    *ll = a;
    *ll0 = b;
}

void ppml_compute_loglik_pair_scaled(const ST_double *y, ST_double y_scale, const ST_double *eta,
                                     const ST_double *offset, ST_double log_scale,
                                     const ST_double *w_user, ST_double eta0, ST_int N,
                                     ST_double *ll_out, ST_double *ll0_out)
{
    /* the null model's eta_total as ppml_compute_loglik forms it */
    eta0 = eta0 + 0.0 + 0.0;
    const ST_double mu0 = exp(eta0);
    ST_int nblocks = (N + PPML_SUM_BLOCK - 1) / PPML_SUM_BLOCK;
    dd_real *part = nblocks > 1 ? (dd_real *)malloc((size_t)nblocks * 2 * sizeof(dd_real)) : NULL;
    dd_real ll = {0.0, 0.0}, ll0 = {0.0, 0.0};
    if (!part) {
        for (ST_int b = 0; b < nblocks; b++) {
            ST_int lo = b * PPML_SUM_BLOCK, hi = lo + PPML_SUM_BLOCK < N ? lo + PPML_SUM_BLOCK : N;
            dd_real a, c;
            ppml_loglik_pair_scaled_block(y, y_scale, eta, offset, log_scale, w_user, eta0, mu0, lo, hi,
                                          &a, &c);
            ll = ppml_dd_add_dd(ll, a);
            ll0 = ppml_dd_add_dd(ll0, c);
        }
    } else {
        #pragma omp parallel for schedule(static)
        for (ST_int b = 0; b < nblocks; b++) {
            ST_int lo = b * PPML_SUM_BLOCK, hi = lo + PPML_SUM_BLOCK < N ? lo + PPML_SUM_BLOCK : N;
            ppml_loglik_pair_scaled_block(y, y_scale, eta, offset, log_scale, w_user, eta0, mu0, lo, hi,
                                          &part[2 * (size_t)b], &part[2 * (size_t)b + 1]);
        }
        for (ST_int b = 0; b < nblocks; b++) {
            ll = ppml_dd_add_dd(ll, part[2 * (size_t)b]);
            ll0 = ppml_dd_add_dd(ll0, part[2 * (size_t)b + 1]);
        }
        free(part);
    }
    *ll_out = ll.hi + ll.lo;
    *ll0_out = ll0.hi + ll0.lo;
}

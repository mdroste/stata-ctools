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
    ST_int g, i, lev;

    S->sum_weights = ppml_block_sum(irls_w, N);

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
   r = y/mu - 1 is exact to rounding. y = 0 gives mu. */
static ST_double ppml_unit_deviance_mu(ST_double y, ST_double mu, ST_double eta_total)
{
    if (y <= 0.0) return mu;
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
    return y * (log(y) - eta_total) - (y - mu);
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

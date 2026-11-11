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

void ppml_update_fe_weights(HDFE_State *S, const ST_double *irls_w, ST_int N)
{
    ST_int g, i, lev;
    ST_double sum_w = 0.0;

    for (i = 0; i < N; i++)
        sum_w += irls_w[i];
    S->sum_weights = sum_w;

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
static ST_double ppml_unit_deviance(ST_double y, ST_double eta_total)
{
    ST_double mu = exp(eta_total);
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

ST_double ppml_compute_deviance(const ST_double *y, const ST_double *eta, const ST_double *offset,
                                const ST_double *w_user, ST_int N)
{
    dd_real dev = {0.0, 0.0};
    for (ST_int i = 0; i < N; i++) {
        ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        ST_double eta_total = eta[i] + (offset != NULL ? offset[i] : 0.0);
        dev = dd_add_d(dev, w * ppml_unit_deviance(y[i], eta_total));
    }
    return 2.0 * (dev.hi + dev.lo);
}

/* ll = sum of y*log(mu) - mu - lngamma(y+1) = saturated ll - deviance/2 */
ST_double ppml_compute_loglik(const ST_double *y, const ST_double *eta, const ST_double *offset,
                              ST_double log_scale, const ST_double *w_user, ST_int N)
{
    dd_real ll = {0.0, 0.0};
    for (ST_int i = 0; i < N; i++) {
        ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        ST_double eta_total = eta[i] + (offset != NULL ? offset[i] : 0.0) + log_scale;
        ll = dd_add_d(ll, w * (ppml_saturated_loglik(y[i]) - ppml_unit_deviance(y[i], eta_total)));
    }
    return ll.hi + ll.lo;
}

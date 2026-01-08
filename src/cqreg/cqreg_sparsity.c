/*
 * cqreg_sparsity.c
 *
 * Sparsity estimation for quantile regression variance computation.
 * Part of the ctools suite.
 */

#include "cqreg_sparsity.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Mathematical constants */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

#ifndef M_SQRT_2PI
#define M_SQRT_2PI 2.50662827463100050242
#endif

/* Default alpha for bandwidth selection */
#define CQREG_BW_ALPHA 0.05

/* ============================================================================
 * Statistical Functions
 * ============================================================================ */

ST_double cqreg_dnorm(ST_double x)
{
    return exp(-0.5 * x * x) / M_SQRT_2PI;
}

/*
 * Standard normal CDF using Abramowitz & Stegun approximation.
 * Maximum error: 7.5e-8
 */
ST_double cqreg_pnorm(ST_double x)
{
    const ST_double a1 =  0.254829592;
    const ST_double a2 = -0.284496736;
    const ST_double a3 =  1.421413741;
    const ST_double a4 = -1.453152027;
    const ST_double a5 =  1.061405429;
    const ST_double p  =  0.3275911;

    int sign = 1;
    if (x < 0) {
        sign = -1;
        x = -x;
    }

    ST_double t = 1.0 / (1.0 + p * x);
    ST_double t2 = t * t;
    ST_double t3 = t2 * t;
    ST_double t4 = t3 * t;
    ST_double t5 = t4 * t;

    ST_double y = 1.0 - (a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t5) * exp(-x*x/2.0);

    return 0.5 * (1.0 + sign * y);
}

/*
 * Standard normal quantile function using rational approximation.
 * Accurate to about 1e-9.
 */
ST_double cqreg_qnorm(ST_double p)
{
    /* Coefficients for rational approximation */
    static const ST_double a[] = {
        -3.969683028665376e+01,
         2.209460984245205e+02,
        -2.759285104469687e+02,
         1.383577518672690e+02,
        -3.066479806614716e+01,
         2.506628277459239e+00
    };
    static const ST_double b[] = {
        -5.447609879822406e+01,
         1.615858368580409e+02,
        -1.556989798598866e+02,
         6.680131188771972e+01,
        -1.328068155288572e+01
    };
    static const ST_double c[] = {
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
         4.374664141464968e+00,
         2.938163982698783e+00
    };
    static const ST_double d[] = {
         7.784695709041462e-03,
         3.224671290700398e-01,
         2.445134137142996e+00,
         3.754408661907416e+00
    };

    static const ST_double p_low  = 0.02425;
    static const ST_double p_high = 1.0 - p_low;

    ST_double q, r;

    if (p <= 0.0) return -INFINITY;
    if (p >= 1.0) return INFINITY;

    if (p < p_low) {
        /* Rational approximation for lower region */
        q = sqrt(-2.0 * log(p));
        return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
               ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    }
    else if (p <= p_high) {
        /* Rational approximation for central region */
        q = p - 0.5;
        r = q * q;
        return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q /
               (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0);
    }
    else {
        /* Rational approximation for upper region */
        q = sqrt(-2.0 * log(1.0 - p));
        return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
                ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    }
}

/* ============================================================================
 * Kernel Functions
 * ============================================================================ */

ST_double cqreg_kernel_epanechnikov(ST_double u)
{
    if (fabs(u) >= 1.0) return 0.0;
    return 0.75 * (1.0 - u * u);
}

ST_double cqreg_kernel_gaussian(ST_double u)
{
    return cqreg_dnorm(u);
}

ST_double cqreg_kernel_triangular(ST_double u)
{
    ST_double abs_u = fabs(u);
    if (abs_u >= 1.0) return 0.0;
    return 1.0 - abs_u;
}

/* ============================================================================
 * Bandwidth Selection
 * ============================================================================ */

ST_double cqreg_bandwidth_hsheather(ST_int N, ST_double q, ST_double alpha)
{
    /*
     * Hall-Sheather (1988) optimal bandwidth for quantile regression.
     *
     * Formula (from R's quantreg package):
     * h = n^{-1/3} * z_{1-alpha/2}^{2/3} * ((1.5 * phi(z_q)^2) / (2*z_q^2 + 1))^{1/3}
     *
     * where:
     *   z_q = Phi^{-1}(q)     - quantile of standard normal at q
     *   phi(z_q) = N(z_q; 0,1) - standard normal density at z_q
     *   z_{1-alpha/2}         - critical value for confidence level (default alpha=0.05)
     *
     * Reference: Hall, P. and Sheather, S.J. (1988), JRSS(B), 50, 381-391.
     */
    ST_double z_alpha = cqreg_qnorm(1.0 - alpha / 2.0);
    ST_double z_q = cqreg_qnorm(q);
    ST_double phi_zq = cqreg_dnorm(z_q);

    ST_double numer = 1.5 * phi_zq * phi_zq;
    ST_double denom = 2.0 * z_q * z_q + 1.0;  /* NOTE: uses z_q, not z_alpha! */

    ST_double h = pow(z_alpha, 2.0/3.0) * pow(numer / denom, 1.0/3.0) * pow((ST_double)N, -1.0/3.0);

    return h;
}

ST_double cqreg_bandwidth_bofinger(ST_int N, ST_double q)
{
    /*
     * Bofinger (1975) bandwidth:
     * h = (9/2 * phi(z_q)^4 * (2*z_q^2 + 1)^2 / n)^{1/5}
     */
    ST_double z_q = cqreg_qnorm(q);
    ST_double phi_zq = cqreg_dnorm(z_q);

    ST_double phi4 = pow(phi_zq, 4);
    ST_double term = 2.0 * z_q * z_q + 1.0;

    ST_double h = pow(4.5 * phi4 * term * term / (ST_double)N, 0.2);

    return h;
}

ST_double cqreg_bandwidth_chamberlain(ST_int N, ST_double q, ST_double alpha)
{
    /*
     * Chamberlain (1994) bandwidth:
     * h = z_{1-alpha/2} * sqrt(q*(1-q) / n) / phi(z_q)
     */
    ST_double z_alpha = cqreg_qnorm(1.0 - alpha / 2.0);
    ST_double z_q = cqreg_qnorm(q);
    ST_double phi_zq = cqreg_dnorm(z_q);

    if (phi_zq < 1e-10) {
        phi_zq = 1e-10;  /* Prevent division by zero for extreme quantiles */
    }

    ST_double h = z_alpha * sqrt(q * (1.0 - q) / (ST_double)N) / phi_zq;

    return h;
}

ST_double cqreg_compute_bandwidth(ST_int N, ST_double q, cqreg_bw_method method)
{
    switch (method) {
        case CQREG_BW_HSHEATHER:
            return cqreg_bandwidth_hsheather(N, q, CQREG_BW_ALPHA);
        case CQREG_BW_BOFINGER:
            return cqreg_bandwidth_bofinger(N, q);
        case CQREG_BW_CHAMBERLAIN:
            return cqreg_bandwidth_chamberlain(N, q, CQREG_BW_ALPHA);
        default:
            return cqreg_bandwidth_hsheather(N, q, CQREG_BW_ALPHA);
    }
}

/* ============================================================================
 * Sorting
 * ============================================================================ */

static int compare_double(const void *a, const void *b)
{
    ST_double da = *(const ST_double *)a;
    ST_double db = *(const ST_double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

void cqreg_sort_ascending(ST_double *data, ST_int N)
{
    qsort(data, N, sizeof(ST_double), compare_double);
}

/* ============================================================================
 * Sample Quantile
 * ============================================================================ */

ST_double cqreg_sample_quantile(const ST_double *sorted_data, ST_int N, ST_double q)
{
    /*
     * Stata-compatible quantile calculation.
     *
     * Stata's _pctile uses ceiling(N * q) as the index (1-based).
     * For 0-based indexing: index = ceiling(N * q) - 1
     *
     * This returns the smallest value x such that at least q*100% of
     * the observations are <= x.
     */

    if (N <= 0) return 0.0;
    if (N == 1) return sorted_data[0];

    /* Compute the 1-based index using ceiling */
    ST_double raw_index = (ST_double)N * q;
    ST_int idx = (ST_int)ceil(raw_index);

    /* Convert to 0-based and clamp */
    idx = idx - 1;
    if (idx < 0) idx = 0;
    if (idx >= N) idx = N - 1;

    return sorted_data[idx];
}

/* ============================================================================
 * Density Estimation
 * ============================================================================ */

ST_double cqreg_kernel_density(ST_double x,
                               const ST_double *data,
                               ST_int N,
                               ST_double bandwidth)
{
    ST_double sum = 0.0;
    ST_double h_inv = 1.0 / bandwidth;

    #pragma omp parallel for reduction(+:sum) if(N > 10000)
    for (ST_int i = 0; i < N; i++) {
        ST_double u = (x - data[i]) * h_inv;
        sum += cqreg_kernel_epanechnikov(u);
    }

    return sum * h_inv / (ST_double)N;
}

ST_double cqreg_density_at_quantile(const ST_double *data,
                                    ST_int N,
                                    ST_double q,
                                    ST_double bandwidth)
{
    /* Find the sample quantile */
    ST_double x_q = cqreg_sample_quantile(data, N, q);

    /* Estimate density at the quantile */
    return cqreg_kernel_density(x_q, data, N, bandwidth);
}

/* ============================================================================
 * Fitted Density Estimation (Stata's default method)
 * ============================================================================ */

/*
 * Estimate per-observation densities using the "fitted" method.
 *
 * This is Stata's default density estimation for qreg:
 * 1. Compute fitted values yhat = X * beta
 * 2. For each observation, use local kernel density estimation
 * 3. The density is estimated at residual = 0 (the quantile hyperplane)
 *
 * Parameters:
 *   obs_density   - Output: per-observation density estimates (N)
 *   residuals     - Regression residuals (N)
 *   N             - Number of observations
 *   q             - Quantile
 *   bw_method     - Bandwidth selection method
 *
 * Returns:
 *   Average sparsity (for compatibility with IID method)
 */
ST_double cqreg_estimate_fitted_density(ST_double *obs_density,
                                        const ST_double *residuals,
                                        ST_int N,
                                        ST_double q,
                                        cqreg_bw_method bw_method)
{
    if (obs_density == NULL || residuals == NULL || N <= 0) {
        return 1.0;
    }

    /*
     * The fitted method estimates local density using kernel smoothing.
     * For quantile regression, the density at the quantile is estimated
     * by counting observations in a bandwidth around residual = 0.
     *
     * Bandwidth is computed in probability units, then scaled by the
     * local dispersion of residuals.
     */

    /* Compute bandwidth in probability units */
    ST_double h_prob = cqreg_compute_bandwidth(N, q, bw_method);

    /* Sort residuals to find quantiles */
    ST_double *sorted = (ST_double *)malloc(N * sizeof(ST_double));
    if (sorted == NULL) {
        /* Fallback: uniform density */
        for (ST_int i = 0; i < N; i++) {
            obs_density[i] = 1.0;
        }
        return 1.0;
    }
    memcpy(sorted, residuals, N * sizeof(ST_double));
    cqreg_sort_ascending(sorted, N);

    /* Find the residual scale from the bandwidth region around 0 */
    ST_double q_lo = q - h_prob;
    ST_double q_hi = q + h_prob;

    /* Clamp to valid range */
    if (q_lo < 0.01) q_lo = 0.01;
    if (q_hi > 0.99) q_hi = 0.99;

    /* Get corresponding residual quantiles */
    ST_double r_lo = cqreg_sample_quantile(sorted, N, q_lo);
    ST_double r_hi = cqreg_sample_quantile(sorted, N, q_hi);

    /* Bandwidth in residual scale */
    ST_double h_resid = (r_hi - r_lo) / 2.0;
    if (h_resid < 1e-10) {
        h_resid = 1.0;  /* Fallback for degenerate cases */
    }

    /*
     * Use kernel density estimation for each observation.
     * The density at residual=0 is estimated using Epanechnikov kernel.
     */
    ST_double total_density = 0.0;

    #pragma omp parallel for reduction(+:total_density) if(N > 10000)
    for (ST_int i = 0; i < N; i++) {
        ST_double u = residuals[i] / h_resid;

        /* Epanechnikov kernel at u */
        ST_double K = 0.0;
        if (fabs(u) < 1.0) {
            K = 0.75 * (1.0 - u * u);
        }

        /* Density estimate: K(u) / h */
        obs_density[i] = K / h_resid;

        /* For numerical stability, use a minimum density */
        if (obs_density[i] < 1e-10) {
            obs_density[i] = 1e-10;
        }

        total_density += obs_density[i];
    }

    free(sorted);

    /* Return average sparsity (1/f) */
    ST_double avg_density = total_density / N;
    if (avg_density < 1e-10) avg_density = 1e-10;

    return 1.0 / avg_density;
}

/*
 * Estimate sparsity using kernel density estimation on residuals.
 * This is an alternative to the difference quotient method.
 */
ST_double cqreg_estimate_kernel_density_sparsity(cqreg_sparsity_state *sp,
                                                  const ST_double *residuals)
{
    if (sp == NULL || residuals == NULL) {
        return 0.0;
    }

    ST_int N = sp->N;
    ST_double q = sp->quantile;

    /* Copy residuals and sort */
    memcpy(sp->sorted_resid, residuals, N * sizeof(ST_double));
    cqreg_sort_ascending(sp->sorted_resid, N);

    /* Compute bandwidth */
    ST_double h_prob = cqreg_compute_bandwidth(N, q, sp->bw_method);
    sp->bandwidth = h_prob;

    /* Get the residual at the quantile (should be near 0 for QR) */
    ST_double r_q = cqreg_sample_quantile(sp->sorted_resid, N, q);

    /* Convert bandwidth to residual scale using IQR */
    ST_double r_25 = cqreg_sample_quantile(sp->sorted_resid, N, 0.25);
    ST_double r_75 = cqreg_sample_quantile(sp->sorted_resid, N, 0.75);
    ST_double iqr = r_75 - r_25;
    if (iqr < 1e-10) iqr = 1.0;

    /* Silverman's rule of thumb for bandwidth in residual units */
    ST_double h_resid = 0.9 * iqr / 1.34 * pow((ST_double)N, -0.2);

    /* Estimate density at r_q using kernel density estimation */
    ST_double density = cqreg_kernel_density(r_q, sp->sorted_resid, N, h_resid);

    /* Sparsity = 1 / density */
    if (density < 1e-10) density = 1e-10;
    sp->sparsity = 1.0 / density;

    return sp->sparsity;
}

/* ============================================================================
 * Main Sparsity Estimation
 * ============================================================================ */

/* Threshold for identifying LP basis observations (residuals ≈ 0) */
#define CQREG_BASIS_THRESHOLD 1e-8

ST_double cqreg_estimate_sparsity(cqreg_sparsity_state *sp,
                                  const ST_double *residuals)
{
    if (sp == NULL || residuals == NULL) {
        return 0.0;
    }

    ST_int N = sp->N;
    ST_double q = sp->quantile;

    /*
     * Stata's qreg vce(iid, residual) method:
     * 1. Exclude LP basis observations (those with |residual| ≈ 0)
     * 2. Compute bandwidth using N_adj = N - n_basis
     * 3. Compute difference quotient sparsity from non-basis residuals
     *
     * The basis observations are exactly on the quantile hyperplane and
     * create discontinuities in the empirical quantile function.
     */

    /* Count and copy non-basis residuals */
    ST_int n_basis = 0;
    ST_int N_adj = 0;

    for (ST_int i = 0; i < N; i++) {
        if (fabs(residuals[i]) < CQREG_BASIS_THRESHOLD) {
            n_basis++;
        } else {
            sp->sorted_resid[N_adj++] = residuals[i];
        }
    }

    /* Handle edge case: all observations are basis (shouldn't happen) */
    if (N_adj < 2) {
        /* Fall back to using all residuals */
        memcpy(sp->sorted_resid, residuals, N * sizeof(ST_double));
        N_adj = N;
    }

    /* Sort the non-basis residuals */
    cqreg_sort_ascending(sp->sorted_resid, N_adj);

    /* Compute bandwidth using adjusted N (in probability units, 0-1 scale) */
    ST_double h_prob = cqreg_compute_bandwidth(N_adj, q, sp->bw_method);
    sp->bandwidth = h_prob;

    /*
     * Use the difference quotient (Siddiqui) method:
     * sparsity = (F^{-1}(q+h) - F^{-1}(q-h)) / (2h)
     *
     * where F^{-1} is the empirical quantile function of the residuals.
     */

    /* Compute indices for bandwidth region */
    ST_double q_lo = q - h_prob;
    ST_double q_hi = q + h_prob;

    /* Check for valid range - if out of bounds, cannot compute VCE */
    if (q_lo <= 0.0 || q_hi >= 1.0) {
        /*
         * This matches Stata's behavior: VCE computation fails if
         * tau ± h goes outside (0, 1). Return a fallback value.
         */
        ST_double r_25 = cqreg_sample_quantile(sp->sorted_resid, N_adj, 0.25);
        ST_double r_75 = cqreg_sample_quantile(sp->sorted_resid, N_adj, 0.75);
        sp->sparsity = (r_75 - r_25) / 0.5;  /* IQR-based fallback */
        if (sp->sparsity < 1e-10) sp->sparsity = 1.0;
        return sp->sparsity;
    }

    /* Get corresponding order statistics from non-basis residuals */
    ST_double x_lo = cqreg_sample_quantile(sp->sorted_resid, N_adj, q_lo);
    ST_double x_hi = cqreg_sample_quantile(sp->sorted_resid, N_adj, q_hi);

    /* Difference quotient sparsity */
    ST_double dq = 2.0 * h_prob;
    sp->sparsity = (x_hi - x_lo) / dq;

    /* Ensure positive sparsity */
    if (sp->sparsity < 1e-10) {
        sp->sparsity = 1e-10;
    }

    return sp->sparsity;
}

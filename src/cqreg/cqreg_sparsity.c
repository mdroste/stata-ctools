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
     * Hall-Sheather (1988) optimal bandwidth:
     * h = z_{1-alpha/2}^{2/3} * (1.5 * phi(z_q)^2 / (2*z_{1-alpha/2}^2 + 1))^{1/3} * n^{-1/3}
     *
     * Note: Stata's implementation multiplies by an additional factor of about 2
     * to account for the Epanechnikov kernel efficiency.
     */
    ST_double z_alpha = cqreg_qnorm(1.0 - alpha / 2.0);
    ST_double z_q = cqreg_qnorm(q);
    ST_double phi_zq = cqreg_dnorm(z_q);

    ST_double numer = 1.5 * phi_zq * phi_zq;
    ST_double denom = 2.0 * z_alpha * z_alpha + 1.0;

    ST_double h = pow(z_alpha, 2.0/3.0) * pow(numer / denom, 1.0/3.0) * pow((ST_double)N, -1.0/3.0);

    /* Stata's qreg uses approximately 2x this bandwidth for IID VCE */
    h *= 2.0;

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
    /* Type 7 quantile (R default): linear interpolation */
    ST_double index = (N - 1) * q;
    ST_int lo = (ST_int)floor(index);
    ST_int hi = (ST_int)ceil(index);

    if (lo < 0) lo = 0;
    if (hi >= N) hi = N - 1;

    if (lo == hi) {
        return sorted_data[lo];
    }

    ST_double frac = index - (ST_double)lo;
    return sorted_data[lo] * (1.0 - frac) + sorted_data[hi] * frac;
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
 * Main Sparsity Estimation
 * ============================================================================ */

ST_double cqreg_estimate_sparsity(cqreg_sparsity_state *sp,
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

    /* Compute bandwidth (as a proportion/quantile) */
    ST_double h = cqreg_compute_bandwidth(N, q, sp->bw_method);
    sp->bandwidth = h;

    /*
     * Estimate sparsity using difference quotient method (like Stata's qreg).
     *
     * The sparsity s = 1/f(0) is estimated as:
     * s = (x_{ceil(n*(q+h))} - x_{floor(n*(q-h))}) / (2*h)
     *
     * where x are ordered residuals and h is the bandwidth.
     * This estimates the reciprocal density at the conditional quantile.
     */

    /* Compute indices */
    ST_double q_lo = q - h;
    ST_double q_hi = q + h;

    /* Clamp to valid range */
    if (q_lo < 0.0) q_lo = 0.0;
    if (q_hi > 1.0) q_hi = 1.0;

    /* Get corresponding order statistics */
    ST_double x_lo = cqreg_sample_quantile(sp->sorted_resid, N, q_lo);
    ST_double x_hi = cqreg_sample_quantile(sp->sorted_resid, N, q_hi);

    /* Sparsity = difference in residuals / difference in quantiles
     * This is the slope of the quantile function, which is 1/density
     */
    ST_double dq = q_hi - q_lo;
    if (dq < 1e-10) dq = 2.0 * h;  /* Fallback */

    sp->sparsity = (x_hi - x_lo) / dq;

    /* Ensure positive sparsity */
    if (sp->sparsity < 1e-10) {
        sp->sparsity = 1e-10;
    }

    return sp->sparsity;
}

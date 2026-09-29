/*
 * cqreg_sparsity.c
 *
 * Sparsity estimation for quantile regression variance computation.
 * Part of the ctools suite.
 */

#include "cqreg_sparsity.h"
#include "cqreg_linalg.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include "../ctools_runtime.h"

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

static ST_double cqreg_dnorm(ST_double x)
{
    return exp(-0.5 * x * x) / M_SQRT_2PI;
}

/*
 * Standard normal quantile function using rational approximation.
 * Accurate to about 1e-9.
 */
static ST_double cqreg_qnorm(ST_double p)
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

    if (p <= 0.0) return -1e308;  /* Approximate -INFINITY (avoids -ffast-math issues) */
    if (p >= 1.0) return 1e308;   /* Approximate INFINITY */

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

/* The kernels of qreg's vce(, kernel()) method, as in _qreg_kernel
 * (qreg.ado), evaluated at u = r/kbwidth; Stata's epanechnikov kernel has
 * support |u| < sqrt(5). */
static ST_double cqreg_kernel_weight(cqreg_kernel_type kernel, ST_double u)
{
    ST_double au = fabs(u);

    switch (kernel) {
        case CQREG_KERNEL_EPAN2:
            return (au < 1.0) ? 0.75 * (1.0 - u * u) : 0.0;
        case CQREG_KERNEL_BIWEIGHT:
            return (au < 1.0) ? (15.0 / 16.0) * ((1.0 - u * u) * (1.0 - u * u)) : 0.0;
        case CQREG_KERNEL_COSINE:
            return (au < 0.5) ? 1.0 + cos(2.0 * M_PI * u) : 0.0;
        case CQREG_KERNEL_GAUSSIAN:
            return cqreg_dnorm(u);
        case CQREG_KERNEL_PARZEN:
            if (au <= 0.5) return 4.0 / 3.0 - 8.0 * (u * u) + 8.0 * (au * au * au);
            return (au <= 1.0) ? (8.0 / 3.0) * ((1.0 - au) * (1.0 - au) * (1.0 - au)) : 0.0;
        case CQREG_KERNEL_RECTANGLE:
            return (au < 1.0) ? 0.5 : 0.0;
        case CQREG_KERNEL_TRIANGLE:
            return (au < 1.0) ? 1.0 - au : 0.0;
        case CQREG_KERNEL_EPANECHNIKOV:
        default:
            return (au < sqrt(5.0)) ? 0.75 * (1.0 - 0.2 * (u * u)) / sqrt(5.0) : 0.0;
    }
}

/* ============================================================================
 * Bandwidth Selection
 * ============================================================================ */

static ST_double cqreg_bandwidth_hsheather(ST_int N, ST_double q, ST_double alpha)
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

static ST_double cqreg_bandwidth_bofinger(ST_int N, ST_double q)
{
    /*
     * Bofinger (1975) bandwidth:
     * h = (9/2 * phi(z_q)^4 * (2*z_q^2 + 1)^2 / n)^{1/5}
     */
    ST_double z_q = cqreg_qnorm(q);
    ST_double phi_zq = cqreg_dnorm(z_q);

    ST_double phi4 = pow(phi_zq, 4);
    ST_double term = 2.0 * z_q * z_q + 1.0;

    ST_double h = pow(4.5 * phi4 / (term * term) / (ST_double)N, 0.2);

    return h;
}

static ST_double cqreg_bandwidth_chamberlain(ST_int N, ST_double q, ST_double alpha)
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

    ST_double h = z_alpha * sqrt(q * (1.0 - q) / (ST_double)N);

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
 * Stata Percentiles
 * ============================================================================ */

static int compare_st_int(const void *a, const void *b)
{
    ST_int ia = *(const ST_int *)a;
    ST_int ib = *(const ST_int *)b;
    return (ia > ib) - (ia < ib);
}

/*
 * Percentile by Stata's default definition (_pctile; summarize, detail):
 * with P = n*pct/100, x_(floor(P)+1) if P is not an integer, and the average
 * (x_(P) + x_(P+1))/2 if it is. Selection with O(n) expected work; reorders
 * values (the multiset is unchanged, so repeated calls are allowed).
 */
static ST_double cqreg_stata_pctile(ST_double *values, ST_int n, ST_double pct)
{
    ST_double P = (ST_double)n * pct / 100.0;
    ST_double fl = floor(P);

    if (!(fl >= 0.0)) fl = 0.0;
    if (fl >= (ST_double)n) return cqreg_select(values, n, n - 1);

    ST_int i = (ST_int)fl;
    if (P == fl && i >= 1) {
        /* After selecting x_(i), x_(i+1) is the smallest value to its right */
        ST_double lo = cqreg_select(values, n, i - 1);
        ST_double hi = values[i];
        for (ST_int j = i + 1; j < n; j++) {
            if (values[j] < hi) hi = values[j];
        }
        return (lo + hi) / 2.0;
    }
    return cqreg_select(values, n, i);
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

    /*
     * Stata's qreg vce(iid, residual) method (GetVCE in qreg.ado):
     * 1. Drop exactly the K observations of the linear programming basis
     *    (r(basis)); other zero residuals, common with tied or discrete
     *    data, are kept
     * 2. Compute bandwidth using N_adj = N - K
     * 3. sparsity = (F^{-1}(q+h) - F^{-1}(q-h)) / (2h), where F^{-1} is
     *    Stata's percentile (_pctile) of the remaining residuals
     *
     * OPTIMIZATION: Use quickselect O(N) instead of full sort O(N log N)
     * since we only need order statistics at q-h and q+h.
     */
    sp->status = 0;
    sp->sparsity = 0.0;

    ST_int *drop = NULL;
    ST_int n_drop = (sp->basis != NULL && sp->n_basis > 0) ? sp->n_basis : 0;
    if (n_drop > 0) {
        drop = (ST_int *)malloc((size_t)n_drop * sizeof(ST_int));
        if (drop == NULL) {
            ctools_error("cqreg", "memory allocation failed");
            sp->status = 920;
            return sp->sparsity;
        }
        memcpy(drop, sp->basis, (size_t)n_drop * sizeof(ST_int));
        qsort(drop, (size_t)n_drop, sizeof(ST_int), compare_st_int);
    }

    /* Copy the residuals of the non-basis observations */
    ST_int N_adj = 0, next = 0;
    for (ST_int i = 0; i < N; i++) {
        while (next < n_drop && drop[next] < i) next++;
        if (next < n_drop && drop[next] == i) continue;
        sp->sorted_resid[N_adj++] = residuals[i];
    }
    free(drop);

    /* qreg: an estimate from zero non-basis residuals is an error */
    if (N_adj == 0) {
        ctools_error("cqreg", "no observations remain after dropping residuals "
                     "associated with the linear programming basis");
        sp->status = 2000;
        return sp->sparsity;
    }

    /* Bandwidth for N_adj observations (in probability units, 0-1 scale).
     * h is proportional to n^(-1/3) (Hall-Sheather), n^(-1/5) (Bofinger)
     * or n^(-1/2) (Chamberlain), so the bandwidth that Stata computed for
     * all N observations is rescaled, keeping Stata's invnormal() values. */
    ST_double h_prob;
    if (sp->full_bandwidth > 0.0) {
        ST_double power = (sp->bw_method == CQREG_BW_BOFINGER) ? 0.2 :
                          (sp->bw_method == CQREG_BW_CHAMBERLAIN) ? 0.5 : 1.0 / 3.0;
        h_prob = sp->full_bandwidth * pow((ST_double)N / (ST_double)N_adj, power);
    } else {
        h_prob = cqreg_compute_bandwidth(N_adj, q, sp->bw_method);
    }
    sp->bandwidth = h_prob;

    /* Compute bounds of the bandwidth region */
    ST_double q_lo = q - h_prob;
    ST_double q_hi = q + h_prob;

    /* qreg stops when tau -/+ h leaves [0, 1]: restricting the interval
     * would understate the standard errors. */
    if (q_lo < 0.0 || q_hi > 1.0) {
        ctools_error("cqreg", "VCE computation failed; try a different bandwidth or bsqreg");
        sp->status = 498;
        return sp->sparsity;
    }

    ST_double x_lo = cqreg_stata_pctile(sp->sorted_resid, N_adj, 100.0 * q_lo);
    ST_double x_hi = cqreg_stata_pctile(sp->sorted_resid, N_adj, 100.0 * q_hi);

    /* Difference quotient sparsity */
    ST_double dq = 2.0 * h_prob;
    sp->sparsity = (x_hi - x_lo) / dq;

    /* qreg: a sparsity below machine epsilon cannot give standard errors */
    if (!(sp->sparsity >= DBL_EPSILON)) {
        ctools_error("cqreg", "sparsity estimate of %g is too small; computation "
                     "of coefficient standard errors cannot be completed", sp->sparsity);
        sp->status = 498;
    }

    return sp->sparsity;
}

/* ============================================================================
 * Kernel Density Method
 * ============================================================================ */

ST_double cqreg_estimate_kernel_sparsity(cqreg_sparsity_state *sp,
                                         const ST_double *residuals,
                                         ST_double *kval)
{
    if (sp == NULL || residuals == NULL || kval == NULL) {
        return 0.0;
    }

    ST_int N = sp->N;
    ST_int i;

    /*
     * qreg's vce(, kernel()) method (VCE_kernel, KernelBWidth and
     * _qreg_kernel in qreg.ado), on the residuals of all N observations:
     * 1. sd and interquartile range as from summarize, detail
     * 2. kband = min(sd, IQR/1.34) * (invnormal(q+h) - invnormal(q-h)),
     *    where kfactor holds the difference of the normal quantiles
     * 3. kval_i = K(r_i/kband); the iid sparsity is kband / mean(kval)
     * Sums are compensated (Kahan) to track Stata's accurate sums.
     */
    sp->status = 0;
    sp->sparsity = 0.0;
    sp->kbwidth = 0.0;

    if (!(sp->kfactor > 0.0) || !isfinite(sp->kfactor) || N < 2) {
        ctools_error("cqreg", "VCE computation failed; try a different bandwidth or bsqreg");
        sp->status = 498;
        return sp->sparsity;
    }

    /* Mean and standard deviation (denominator N-1) */
    ST_double sum = 0.0, comp = 0.0;
    for (i = 0; i < N; i++) {
        ST_double yk = residuals[i] - comp;
        ST_double t = sum + yk;
        comp = (t - sum) - yk;
        sum = t;
    }
    ST_double mean = sum / (ST_double)N;
    sum = 0.0; comp = 0.0;
    for (i = 0; i < N; i++) {
        ST_double dev = residuals[i] - mean;
        ST_double yk = dev * dev - comp;
        ST_double t = sum + yk;
        comp = (t - sum) - yk;
        sum = t;
    }
    ST_double sd = sqrt(sum / (ST_double)(N - 1));

    /* Interquartile range by Stata's percentile definition */
    memcpy(sp->sorted_resid, residuals, (size_t)N * sizeof(ST_double));
    ST_double p25 = cqreg_stata_pctile(sp->sorted_resid, N, 25.0);
    ST_double p75 = cqreg_stata_pctile(sp->sorted_resid, N, 75.0);
    ST_double sr = (p75 - p25) / 1.34;

    ST_double kband = ((sd < sr) ? sd : sr) * sp->kfactor;
    sp->kbwidth = kband;

    /* A zero bandwidth (e.g. more than half the residuals tied at zero)
     * leaves the kernel weights undefined */
    if (!(kband > 0.0) || !isfinite(kband)) {
        ctools_error("cqreg", "kernel bandwidth of %g is not positive; computation "
                     "of coefficient standard errors cannot be completed", kband);
        sp->status = 498;
        return sp->sparsity;
    }

    sum = 0.0; comp = 0.0;
    for (i = 0; i < N; i++) {
        kval[i] = cqreg_kernel_weight(sp->kernel, residuals[i] / kband);
        ST_double yk = kval[i] - comp;
        ST_double t = sum + yk;
        comp = (t - sum) - yk;
        sum = t;
    }

    /* Infinite when no residual has positive kernel weight; the iid VCE
     * then fails, while the robust Hessian is simply singular (as in qreg) */
    sp->sparsity = kband / (sum / (ST_double)N);
    return sp->sparsity;
}

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
 * Sorting and Selection
 * ============================================================================ */

static int compare_double(const void *a, const void *b)
{
    ST_double da = *(const ST_double *)a;
    ST_double db = *(const ST_double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

/* ============================================================================
 * Quickselect Algorithm - O(N) expected time for finding k-th element
 * Much faster than full sort O(N log N) when we only need a few order statistics
 * ============================================================================ */

/* Partition helper for quickselect */
static ST_int partition(ST_double *arr, ST_int lo, ST_int hi)
{
    /* Use median-of-three pivot selection for better performance */
    ST_int mid = lo + (hi - lo) / 2;

    /* Sort lo, mid, hi to find median */
    if (arr[mid] < arr[lo]) {
        ST_double tmp = arr[lo]; arr[lo] = arr[mid]; arr[mid] = tmp;
    }
    if (arr[hi] < arr[lo]) {
        ST_double tmp = arr[lo]; arr[lo] = arr[hi]; arr[hi] = tmp;
    }
    if (arr[mid] < arr[hi]) {
        ST_double tmp = arr[mid]; arr[mid] = arr[hi]; arr[hi] = tmp;
    }

    /* Pivot is now at hi */
    ST_double pivot = arr[hi];
    ST_int i = lo - 1;

    for (ST_int j = lo; j < hi; j++) {
        if (arr[j] <= pivot) {
            i++;
            ST_double tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }
    }

    ST_double tmp = arr[i + 1];
    arr[i + 1] = arr[hi];
    arr[hi] = tmp;

    return i + 1;
}

/*
 * Quickselect: Find the k-th smallest element in O(N) expected time.
 * Modifies array in place (partial sorting).
 * After return, arr[k] contains the k-th smallest element.
 */
static void quickselect(ST_double *arr, ST_int lo, ST_int hi, ST_int k)
{
    while (lo < hi) {
        ST_int pivot_idx = partition(arr, lo, hi);

        if (pivot_idx == k) {
            return;
        } else if (k < pivot_idx) {
            hi = pivot_idx - 1;
        } else {
            lo = pivot_idx + 1;
        }
    }
}

/*
 * Find two order statistics efficiently.
 * Uses quickselect twice, leveraging the partial ordering from the first call.
 * Returns the values at the specified indices (0-based).
 * Modifies array in place.
 */
static void find_two_order_statistics(ST_double *arr, ST_int N,
                                       ST_int idx_lo, ST_int idx_hi,
                                       ST_double *val_lo, ST_double *val_hi)
{
    if (idx_lo > idx_hi) {
        /* Swap if needed */
        ST_int tmp = idx_lo;
        idx_lo = idx_hi;
        idx_hi = tmp;
    }

    /* Clamp indices */
    if (idx_lo < 0) idx_lo = 0;
    if (idx_hi >= N) idx_hi = N - 1;

    /* Find the lower order statistic first */
    quickselect(arr, 0, N - 1, idx_lo);
    *val_lo = arr[idx_lo];

    /* Now find the higher order statistic.
     * After the first quickselect, all elements in arr[0..idx_lo] are <= arr[idx_lo],
     * and all elements in arr[idx_lo+1..N-1] are >= arr[idx_lo].
     * So we only need to search in arr[idx_lo+1..N-1] for idx_hi.
     */
    if (idx_hi > idx_lo) {
        quickselect(arr, idx_lo + 1, N - 1, idx_hi);
    }
    *val_hi = arr[idx_hi];
}

/* ============================================================================
 * Parallel Merge Sort for double arrays
 * Uses OpenMP tasks for parallel divide-and-conquer
 * ============================================================================ */

#ifdef _OPENMP
#include <omp.h>

/* Threshold below which we use sequential sort */
#define PARALLEL_SORT_THRESHOLD 100000

/* Merge two sorted subarrays: data[lo..mid] and data[mid+1..hi] */
static void merge_subarrays(ST_double *data, ST_double *temp,
                            ST_int lo, ST_int mid, ST_int hi)
{
    ST_int i = lo, j = mid + 1, k = lo;

    /* Copy to temp */
    for (ST_int t = lo; t <= hi; t++) {
        temp[t] = data[t];
    }

    /* Merge back */
    while (i <= mid && j <= hi) {
        if (temp[i] <= temp[j]) {
            data[k++] = temp[i++];
        } else {
            data[k++] = temp[j++];
        }
    }

    /* Copy remaining elements */
    while (i <= mid) {
        data[k++] = temp[i++];
    }
    while (j <= hi) {
        data[k++] = temp[j++];
    }
}

/* Recursive parallel merge sort */
static void parallel_merge_sort_recursive(ST_double *data, ST_double *temp,
                                          ST_int lo, ST_int hi, int depth)
{
    if (lo >= hi) return;

    ST_int n = hi - lo + 1;

    /* Use qsort for small arrays or deep recursion */
    if (n < PARALLEL_SORT_THRESHOLD || depth > 4) {
        qsort(&data[lo], n, sizeof(ST_double), compare_double);
        return;
    }

    ST_int mid = lo + (hi - lo) / 2;

    /* Parallel recursive calls */
    #pragma omp task shared(data, temp) if(depth < 3)
    parallel_merge_sort_recursive(data, temp, lo, mid, depth + 1);

    #pragma omp task shared(data, temp) if(depth < 3)
    parallel_merge_sort_recursive(data, temp, mid + 1, hi, depth + 1);

    #pragma omp taskwait

    /* Merge the sorted halves */
    merge_subarrays(data, temp, lo, mid, hi);
}

static void parallel_merge_sort(ST_double *data, ST_int N)
{
    if (N <= 1) return;

    /* Allocate temporary buffer */
    ST_double *temp = (ST_double *)malloc(N * sizeof(ST_double));
    if (temp == NULL) {
        /* Fallback to qsort if allocation fails */
        qsort(data, N, sizeof(ST_double), compare_double);
        return;
    }

    #pragma omp parallel
    {
        #pragma omp single
        {
            parallel_merge_sort_recursive(data, temp, 0, N - 1, 0);
        }
    }

    free(temp);
}
#endif /* _OPENMP */

void cqreg_sort_ascending(ST_double *data, ST_int N)
{
#ifdef _OPENMP
    if (N > PARALLEL_SORT_THRESHOLD) {
        parallel_merge_sort(data, N);
    } else {
        qsort(data, N, sizeof(ST_double), compare_double);
    }
#else
    qsort(data, N, sizeof(ST_double), compare_double);
#endif
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
                                        const ST_double *y,
                                        const ST_double *residuals,
                                        ST_int N,
                                        ST_double q,
                                        cqreg_bw_method bw_method)
{
    if (obs_density == NULL || y == NULL || residuals == NULL || N <= 0) {
        return 1.0;
    }

    /*
     * Fitted method density estimation for robust VCE (Powell sandwich).
     *
     * For each observation i, estimate the conditional density of Y at the
     * fitted value ŷ_i:
     *   f_i = (1/(N*h)) * sum_j K((y_j - ŷ_i) / h)
     *
     * where:
     *   ŷ_i = y_i - r_i (fitted value at observation i)
     *   h = bandwidth in Y scale
     *
     * This is Stata's "fitted" method for vce(robust).
     */

    /* Compute bandwidth in probability units */
    ST_double h_prob = cqreg_compute_bandwidth(N, q, bw_method);

    /*
     * Compute bandwidth in Y scale using the fitted values.
     * The bandwidth should be proportional to the spread of fitted values.
     */
    ST_double *fitted = (ST_double *)malloc(N * sizeof(ST_double));
    if (fitted == NULL) {
        /* Fallback: uniform density */
        for (ST_int i = 0; i < N; i++) {
            obs_density[i] = 1.0 / N;
        }
        return (ST_double)N;
    }

    /* Compute fitted values: ŷ_i = y_i - r_i */
    for (ST_int i = 0; i < N; i++) {
        fitted[i] = y[i] - residuals[i];
    }
    cqreg_sort_ascending(fitted, N);

    /* Get bandwidth in fitted value scale from probability bandwidth */
    ST_double q_lo = q - h_prob;
    ST_double q_hi = q + h_prob;
    if (q_lo < 0.01) q_lo = 0.01;
    if (q_hi > 0.99) q_hi = 0.99;

    ST_double fit_lo = cqreg_sample_quantile(fitted, N, q_lo);
    ST_double fit_hi = cqreg_sample_quantile(fitted, N, q_hi);

    /* Bandwidth h in fitted value scale */
    ST_double h = (fit_hi - fit_lo) / 2.0;
    if (h < 1e-10) {
        h = 1.0;
    }

    free(fitted);

    /*
     * Compute per-observation density:
     * f_i = (1/(N*h)) * sum_j K((y_j - ŷ_i) / h)
     *
     * where ŷ_i = y_i - r_i (fitted value)
     */
    ST_double total_density = 0.0;
    ST_double h_inv = 1.0 / h;
    ST_double N_h_inv = 1.0 / ((ST_double)N * h);

    for (ST_int i = 0; i < N; i++) {
        ST_double yhat_i = y[i] - residuals[i];  /* Fitted value at i */
        ST_double sum_K = 0.0;

        /* Sum kernel contributions from all observations */
        for (ST_int j = 0; j < N; j++) {
            ST_double u = (y[j] - yhat_i) * h_inv;

            /* Epanechnikov kernel */
            if (fabs(u) < 1.0) {
                sum_K += 0.75 * (1.0 - u * u);
            }
        }

        obs_density[i] = sum_K * N_h_inv;

        /* Ensure minimum density for numerical stability */
        if (obs_density[i] < 1e-10) {
            obs_density[i] = 1e-10;
        }

        total_density += obs_density[i];
    }

    /* Return average sparsity */
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
     *
     * OPTIMIZATION: Use quickselect O(N) instead of full sort O(N log N)
     * since we only need two order statistics at q-h and q+h.
     */

    /* Copy non-basis residuals (basis obs have |residual| < threshold) */
    ST_int N_adj = 0;

    for (ST_int i = 0; i < N; i++) {
        if (fabs(residuals[i]) >= CQREG_BASIS_THRESHOLD) {
            sp->sorted_resid[N_adj++] = residuals[i];
        }
    }

    /* Handle edge case: all observations are basis (shouldn't happen) */
    if (N_adj < 2) {
        /* Fall back to using all residuals */
        memcpy(sp->sorted_resid, residuals, N * sizeof(ST_double));
        N_adj = N;
    }

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

    /* Check for valid range - if out of bounds, use IQR fallback */
    if (q_lo <= 0.0 || q_hi >= 1.0) {
        /*
         * This matches Stata's behavior: VCE computation fails if
         * tau ± h goes outside (0, 1). Return a fallback value.
         * Use quickselect for IQR computation too.
         */
        ST_int idx_25 = (ST_int)(0.25 * N_adj);
        ST_int idx_75 = (ST_int)(0.75 * N_adj);
        if (idx_25 >= N_adj) idx_25 = N_adj - 1;
        if (idx_75 >= N_adj) idx_75 = N_adj - 1;

        ST_double r_25, r_75;
        find_two_order_statistics(sp->sorted_resid, N_adj, idx_25, idx_75, &r_25, &r_75);

        sp->sparsity = (r_75 - r_25) / 0.5;  /* IQR-based fallback */
        if (sp->sparsity < 1e-10) sp->sparsity = 1.0;
        return sp->sparsity;
    }

    /*
     * OPTIMIZATION: Use quickselect to find order statistics in O(N) time
     * instead of full sort O(N log N). For N=5M, this saves ~0.5-1 second.
     *
     * Compute 0-based indices for the order statistics.
     * Stata uses ceiling(N * q) - 1 for the quantile index.
     */
    ST_int idx_lo = (ST_int)ceil(N_adj * q_lo) - 1;
    ST_int idx_hi = (ST_int)ceil(N_adj * q_hi) - 1;

    /* Clamp indices to valid range */
    if (idx_lo < 0) idx_lo = 0;
    if (idx_hi >= N_adj) idx_hi = N_adj - 1;

    ST_double x_lo, x_hi;
    find_two_order_statistics(sp->sorted_resid, N_adj, idx_lo, idx_hi, &x_lo, &x_hi);

    /* Difference quotient sparsity */
    ST_double dq = 2.0 * h_prob;
    sp->sparsity = (x_hi - x_lo) / dq;

    /* Ensure positive sparsity */
    if (sp->sparsity < 1e-10) {
        sp->sparsity = 1e-10;
    }

    return sp->sparsity;
}

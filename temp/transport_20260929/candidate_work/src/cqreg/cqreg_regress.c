/*
 * cqreg_regress.c
 *
 * Full quantile regression command orchestrator.
 * Part of the ctools suite.
 */

#include "cqreg_regress.h"
#include "cqreg_types.h"
#include "cqreg_ipm.h"
#include "cqreg_fn.h"
#include "cqreg_sparsity.h"
#include "cqreg_vce.h"
#include "cqreg_linalg.h"
#include "../ctools_runtime.h"
#include "../ctools_config.h"
#include "../ctools_types.h"  /* For ctools_data_load */
#include "../ctools_spi.h"
#include "../ctools_ols.h"    /* For detect_collinearity, fast_dot */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdarg.h>

/* Debug logging to file */
#define MAIN_DEBUG 0

#if MAIN_DEBUG
static FILE *main_debug_file = NULL;

static void main_debug_open(void) {
    if (main_debug_file == NULL) {
        main_debug_file = fopen("/tmp/cqreg_main_debug.log", "a");
        if (main_debug_file) {
            fprintf(main_debug_file, "\n=== cqreg_main Session Started ===\n");
            fflush(main_debug_file);
        }
    }
}

static void main_debug_close(void) {
    if (main_debug_file) {
        fprintf(main_debug_file, "=== cqreg_main Session Ended ===\n\n");
        fflush(main_debug_file);
        fclose(main_debug_file);
        main_debug_file = NULL;
    }
}

static void main_debug_log(const char *fmt, ...) {
    if (main_debug_file) {
        va_list args;
        va_start(args, fmt);
        vfprintf(main_debug_file, fmt, args);
        va_end(args);
        fflush(main_debug_file);
    }
}
#else
#define main_debug_open()
#define main_debug_close()
#define main_debug_log(...)
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/* ============================================================================
 * Helper: Read scalar from Stata
 * ============================================================================ */

static ST_double read_scalar(const char *name, ST_double default_val)
{
    ST_double val;
    if (SF_scal_use((char *)name, &val) != 0) {
        return default_val;
    }
    return val;
}

static ST_int read_scalar_int(const char *name, ST_int default_val)
{
    ST_double val;
    if (SF_scal_use((char *)name, &val) != 0) {
        return default_val;
    }
    return (ST_int)val;
}

/* ============================================================================
 * Helper: Store scalar to Stata
 * ============================================================================ */

static void store_scalar(const char *name, ST_double val)
{
    SF_scal_save((char *)name, val);
}

/* ============================================================================
 * Helper: Load data from Stata using ctools_data_load (parallel loading)
 *
 * This function uses the common ctools_data_load which:
 * - Loads all variables in parallel (one thread per variable)
 * - Uses 8x loop unrolling for throughput
 * - Returns data in stata_data structure
 *
 * After loading:
 * - state->y points to data->vars[0].data.dbl (depvar, view, not copied)
 * - state->X is allocated as contiguous N*K array with X columns + constant
 * ============================================================================ */

/*
 * Load data using ctools_data_load for efficient parallel loading.
 *
 * Uses ctools_data_load() which:
 * - Handles if/in filtering at load time
 * - Loads only filtered observations (memory efficient)
 * - Returns obs_map for write-back operations
 *
 * The filtered data is already compact, so we just copy directly.
 */
static ST_int load_data_filtered(cqreg_state *state,
                                  ctools_filtered_data *filtered,
                                  ST_int K_x)
{
    ST_int K = K_x + 1;  /* +1 for constant */
    ST_int i, k;
    ST_int N = (ST_int)filtered->data.nobs;

    main_debug_log("load_data_filtered: N=%d, K_x=%d\n", N, K_x);

    /* Allocate final y and X arrays */
    size_t y_size;
    if (ctools_safe_mul_size((size_t)N, sizeof(ST_double), &y_size) != 0) {
        ctools_error("cqreg", "Overflow computing y array size");
        return -1;
    }
    state->y = (ST_double *)ctools_cacheline_alloc(y_size);
    if (state->y == NULL) {
        ctools_error("cqreg", "Failed to allocate y array");
        return -1;
    }
    state->y_owned = 1;

    size_t x_size;
    if (ctools_safe_alloc_size((size_t)N, (size_t)K, sizeof(ST_double), &x_size) != 0) {
        ctools_error("cqreg", "Overflow computing X matrix size");
        return -1;
    }
    state->X = (ST_double *)ctools_cacheline_alloc(x_size);
    if (state->X == NULL) {
        ctools_error("cqreg", "Failed to allocate X matrix");
        return -1;
    }

    /* Direct copy - data is already filtered, no indirection needed */
    memcpy(state->y, filtered->data.vars[0].data.dbl, (size_t)N * sizeof(ST_double));

    /* Copy X columns in parallel - direct memcpy */
    #pragma omp parallel for if(K_x >= 2)
    for (k = 0; k < K_x; k++) {
        memcpy(&state->X[(size_t)k * N], filtered->data.vars[1 + k].data.dbl,
               (size_t)N * sizeof(ST_double));
    }

    /* Fill constant column */
    ST_double *const_col = &state->X[(size_t)K_x * N];
    for (i = 0; i < N; i++) {
        const_col[i] = 1.0;
    }

    main_debug_log("load_data_filtered: X matrix populated, constant added\n");

    return 0;
}

/* qsort comparator for ST_int */
static int cmp_st_int(const void *a, const void *b)
{
    ST_int va = *(const ST_int *)a;
    ST_int vb = *(const ST_int *)b;
    return (va > vb) - (va < vb);
}

/* ============================================================================
 * Helper: Load cluster variable from pre-loaded data
 * ============================================================================ */

static ST_int load_clusters_from_data(cqreg_state *state,
                                      const ST_double *cluster_data,
                                      ST_int N)
{
    ST_int idx;

    /* Use ctools_cacheline_alloc to match ctools_aligned_free in cqreg_state_free */
    state->cluster_ids = (ST_int *)ctools_cacheline_alloc(N * sizeof(ST_int));
    if (state->cluster_ids == NULL) {
        return -1;
    }

    /* Convert pre-loaded doubles to integer cluster IDs */
    for (idx = 0; idx < N; idx++) {
        state->cluster_ids[idx] = (ST_int)cluster_data[idx];
    }

    /* Count unique clusters via sort: O(N log N) instead of O(N * G) */
    ST_int *sorted = (ST_int *)malloc((size_t)N * sizeof(ST_int));
    if (sorted == NULL) {
        return -1;
    }
    memcpy(sorted, state->cluster_ids, (size_t)N * sizeof(ST_int));
    qsort(sorted, (size_t)N, sizeof(ST_int), cmp_st_int);

    ST_int num_unique = 1;
    for (idx = 1; idx < N; idx++) {
        if (sorted[idx] != sorted[idx - 1]) {
            num_unique++;
        }
    }
    free(sorted);

    state->num_clusters = num_unique;

    return 0;
}

/* ============================================================================
 * Helper: Siddiqui (Fitted) Sparsity Estimation
 * ============================================================================ */

/* qreg's r(498) messages when the density/sparsity step cannot be completed */
static void cqreg_vce_bandwidth_error(void)
{
    ctools_error("cqreg", "VCE computation failed; try a different bandwidth or bsqreg");
}

static void cqreg_vce_iterations_error(void)
{
    ctools_error("cqreg", "VCE computation failed; try increasing the maximum "
                 "number of iterations or try bsqreg");
}

/*
 * Estimate sparsity using the Siddiqui difference quotient method.
 * This is Stata's default "fitted" method for vce(iid).
 *
 * The method:
 * 1. Fit quantile regression at τ-h and τ+h
 * 2. Evaluate predicted quantiles at X̄ (mean of X)
 * 3. Sparsity = (Q̂(τ+h|X̄) - Q̂(τ-h|X̄)) / (2h)
 *
 * OPTIMIZATION: The auxiliary solves use relaxed tolerance (1e-6 vs 1e-12)
 * and early stopping since we only need ~3 decimal places of accuracy for
 * the sparsity estimate. This typically reduces iterations by 30-50%.
 *
 * Reference: Hendricks and Koenker (1992), Koenker (2005, sec. 3.4)
 */
static ST_retcode estimate_siddiqui_sparsity(const ST_double *y,
                                             const ST_double *X,
                                             ST_int N, ST_int K,
                                             ST_double quantile,
                                             cqreg_bw_method bw_method,
                                             cqreg_ipm_state *main_ipm,
                                             ST_double precomputed_bw,
                                             ST_double *out_sparsity)
{
    ST_double h = (precomputed_bw > 0) ? precomputed_bw : cqreg_compute_bandwidth(N, quantile, bw_method);
    ST_double q_lo = quantile - h;
    ST_double q_hi = quantile + h;

    /* qreg stops when tau -/+ h leaves (0, 1): restricting the interval
     * would understate the standard errors. */
    if (q_lo <= 0.0 || q_hi >= 1.0) {
        cqreg_vce_bandwidth_error();
        return 498;
    }

    /* Allocate temporary arrays for beta coefficients */
    ST_double *beta_lo = (ST_double *)malloc(K * sizeof(ST_double));
    ST_double *beta_hi = (ST_double *)malloc(K * sizeof(ST_double));

    if (beta_lo == NULL || beta_hi == NULL) {
        free(beta_lo);
        free(beta_hi);
        ctools_error_alloc("cqreg");
        return 920;
    }

    /*
     * OPTIMIZATION: Run both auxiliary solves in parallel.
     * - Reuse the main IPM state for one solve (zero allocation cost).
     * - Create a lightweight FN state for the other solve, without the
     *   legacy IPM workspace or N-sized per-thread scratch buffers.
     * - skip_crossover=1 skips the expensive crossover + residual
     *   computation since we only need beta for sparsity estimation.
     */
    ST_int skip_cross = (N > 1000) ? 1 : 0;

    cqreg_ipm_config aux_config;
    cqreg_ipm_config_init(&aux_config);
    aux_config.maxiter = main_ipm->config.maxiter;  /* qreg passes iterate() to the auxiliary fits */
    aux_config.tol_primal = 1e-6;
    aux_config.tol_dual = 1e-6;
    aux_config.tol_gap = 1e-6;
    aux_config.verbose = 0;
    aux_config.skip_crossover = skip_cross;

    cqreg_ipm_state *ipm_hi = cqreg_ipm_create_lite(N, K, &aux_config);
    if (ipm_hi == NULL) {
        free(beta_lo);
        free(beta_hi);
        ctools_error_alloc("cqreg");
        return 920;
    }

    /* Temporarily reconfigure the main state for the lo solve */
    cqreg_ipm_config saved_config = main_ipm->config;
    main_ipm->config = aux_config;

    ST_int rc_lo = 0, rc_hi = 0;
    double t_aux_start = ctools_timer_seconds();

#ifdef _OPENMP
    #pragma omp parallel sections
    {
        #pragma omp section
        { rc_lo = cqreg_fn_solve(main_ipm, y, X, q_lo, beta_lo); }
        #pragma omp section
        { rc_hi = cqreg_fn_solve(ipm_hi, y, X, q_hi, beta_hi); }
    }
#else
    rc_lo = cqreg_fn_solve(main_ipm, y, X, q_lo, beta_lo);
    rc_hi = cqreg_fn_solve(ipm_hi, y, X, q_hi, beta_hi);
#endif

    double t_aux_elapsed = ctools_timer_seconds() - t_aux_start;
    SF_scal_save("_cqreg_time_sparsity_aux", t_aux_elapsed);

    /* Restore original config and free lite state */
    main_ipm->config = saved_config;
    cqreg_ipm_free(ipm_hi);

    main_debug_log("Siddiqui: q_lo=%.4f, q_hi=%.4f, rc_lo=%d, rc_hi=%d\n",
                   q_lo, q_hi, rc_lo, rc_hi);
    main_debug_log("Siddiqui: beta_lo[0]=%.4f, beta_lo[1]=%.4f\n",
                   beta_lo[0], K > 1 ? beta_lo[1] : 0.0);
    main_debug_log("Siddiqui: beta_hi[0]=%.4f, beta_hi[1]=%.4f\n",
                   beta_hi[0], K > 1 ? beta_hi[1] : 0.0);

    /* Check convergence */
    if (rc_lo <= 0 || rc_hi <= 0) {
        main_debug_log("Siddiqui: FAILED - auxiliary solves did not converge\n");
        free(beta_lo);
        free(beta_hi);
        cqreg_vce_iterations_error();
        return 498;
    }

    /* Compute X̄ (mean of each column of X)
     * OPTIMIZATION: Parallel column summation with 8x unrolling */
    ST_double *x_bar = (ST_double *)calloc(K, sizeof(ST_double));
    if (x_bar == NULL) {
        free(beta_lo);
        free(beta_hi);
        ctools_error_alloc("cqreg");
        return 920;
    }

    ST_double inv_N = 1.0 / (ST_double)N;
    ST_int N8 = N - (N & 7);

    #pragma omp parallel for if(K >= 4 && N > 10000)
    for (ST_int k = 0; k < K; k++) {
        const ST_double *Xk = &X[(size_t)k * N];
        ST_double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
        ST_double s4 = 0.0, s5 = 0.0, s6 = 0.0, s7 = 0.0;
        ST_int i;
        for (i = 0; i < N8; i += 8) {
            s0 += Xk[i];     s1 += Xk[i + 1];
            s2 += Xk[i + 2]; s3 += Xk[i + 3];
            s4 += Xk[i + 4]; s5 += Xk[i + 5];
            s6 += Xk[i + 6]; s7 += Xk[i + 7];
        }
        ST_double sum = ((s0 + s4) + (s1 + s5)) + ((s2 + s6) + (s3 + s7));
        for (; i < N; i++) {
            sum += Xk[i];
        }
        x_bar[k] = sum * inv_N;
    }

    /* Compute Q̂(τ-h|X̄) = X̄'β_lo and Q̂(τ+h|X̄) = X̄'β_hi */
    ST_double Qhat_lo = 0.0;
    ST_double Qhat_hi = 0.0;
    for (ST_int k = 0; k < K; k++) {
        Qhat_lo += x_bar[k] * beta_lo[k];
        Qhat_hi += x_bar[k] * beta_hi[k];
    }

    /* Siddiqui sparsity = (Q̂(τ+h|X̄) - Q̂(τ-h|X̄)) / (2h) */
    ST_double sparsity = (Qhat_hi - Qhat_lo) / (2.0 * h);

    main_debug_log("Siddiqui: Qhat_lo=%.4f, Qhat_hi=%.4f, h=%.6f\n",
                   Qhat_lo, Qhat_hi, h);
    main_debug_log("Siddiqui: sparsity = (%.4f - %.4f) / (2*%.6f) = %.4f\n",
                   Qhat_hi, Qhat_lo, h, sparsity);

    /* Reported as computed: qreg's fitted iid VCE uses sparsity^2. */

    /* Cleanup */
    free(beta_lo);
    free(beta_hi);
    free(x_bar);

    main_debug_log("Siddiqui: returning sparsity=%.4f\n", sparsity);
    *out_sparsity = sparsity;
    return 0;
}

/*
 * Estimate per-observation densities for robust VCE using Stata's fitted method.
 *
 * This implements Stata's qreg vce(robust) approach:
 * 1. Fit QR at τ-h to get β_lo
 * 2. Fit QR at τ+h to get β_hi
 * 3. For each observation i:
 *    f_i = (2*h) / (x_i' * β_hi - x_i' * β_lo)
 *
 * This is the "fitted" method density estimator - it uses predicted values
 * from auxiliary quantile regressions rather than kernel density estimation.
 *
 * Parameters:
 *   obs_density  - Output: per-observation density estimates (N)
 *   y            - Response variable (N)
 *   X            - Design matrix (N x K, column-major)
 *   N            - Number of observations
 *   K            - Number of regressors
 *   quantile     - Target quantile (0 < τ < 1)
 *   bw_method    - Bandwidth selection method
 *   ipm_config   - IPM solver configuration
 *   main_beta    - Beta from main solve (for warm start)
 *   out_bandwidth - Output: bandwidth used
 *
 * Returns:
 *   Average sparsity (for compatibility with IID method)
 */
static ST_retcode estimate_fitted_per_obs_density(ST_double *obs_density,
                                                  const ST_double *y,
                                                  const ST_double *X,
                                                  ST_int N, ST_int K,
                                                  ST_double quantile,
                                                  cqreg_bw_method bw_method,
                                                  cqreg_ipm_state *main_ipm,
                                                  ST_double *out_bandwidth,
                                                  ST_double precomputed_bw,
                                                  ST_double *out_sparsity)
{
    ST_double h = (precomputed_bw > 0) ? precomputed_bw : cqreg_compute_bandwidth(N, quantile, bw_method);
    ST_double q_lo = quantile - h;
    ST_double q_hi = quantile + h;

    *out_bandwidth = h;

    /* qreg stops when tau -/+ h leaves (0, 1) */
    if (q_lo <= 0.0 || q_hi >= 1.0) {
        cqreg_vce_bandwidth_error();
        return 498;
    }

    /* Allocate temporary arrays for beta coefficients */
    ST_double *beta_lo = (ST_double *)malloc(K * sizeof(ST_double));
    ST_double *beta_hi = (ST_double *)malloc(K * sizeof(ST_double));

    if (beta_lo == NULL || beta_hi == NULL) {
        free(beta_lo);
        free(beta_hi);
        ctools_error_alloc("cqreg");
        return 920;
    }

    /*
     * OPTIMIZATION: Parallel auxiliary solves — same approach as
     * estimate_siddiqui_sparsity. Reuse main IPM state for one solve,
     * lightweight state for the other.
     */
    ST_int skip_cross = (N > 1000) ? 1 : 0;

    cqreg_ipm_config aux_config;
    cqreg_ipm_config_init(&aux_config);
    aux_config.maxiter = main_ipm->config.maxiter;  /* qreg passes iterate() to the auxiliary fits */
    aux_config.tol_primal = 1e-6;
    aux_config.tol_dual = 1e-6;
    aux_config.tol_gap = 1e-6;
    aux_config.verbose = 0;
    aux_config.skip_crossover = skip_cross;

    cqreg_ipm_state *ipm_hi = cqreg_ipm_create_lite(N, K, &aux_config);
    if (ipm_hi == NULL) {
        free(beta_lo);
        free(beta_hi);
        ctools_error_alloc("cqreg");
        return 920;
    }

    cqreg_ipm_config saved_config = main_ipm->config;
    main_ipm->config = aux_config;

    ST_int rc_lo = 0, rc_hi = 0;
    double t_aux_start2 = ctools_timer_seconds();

#ifdef _OPENMP
    #pragma omp parallel sections
    {
        #pragma omp section
        { rc_lo = cqreg_fn_solve(main_ipm, y, X, q_lo, beta_lo); }
        #pragma omp section
        { rc_hi = cqreg_fn_solve(ipm_hi, y, X, q_hi, beta_hi); }
    }
#else
    rc_lo = cqreg_fn_solve(main_ipm, y, X, q_lo, beta_lo);
    rc_hi = cqreg_fn_solve(ipm_hi, y, X, q_hi, beta_hi);
#endif

    SF_scal_save("_cqreg_time_sparsity_aux", ctools_timer_seconds() - t_aux_start2);

    main_ipm->config = saved_config;
    cqreg_ipm_free(ipm_hi);

    /* Check convergence */
    if (rc_lo <= 0 || rc_hi <= 0) {
        free(beta_lo);
        free(beta_hi);
        cqreg_vce_iterations_error();
        return 498;
    }

    /*
     * Compute per-observation densities using Stata's formula:
     * f_i = (2*h) / (xb_hi - xb_lo)
     *
     * where xb_lo = x_i' * beta_lo, xb_hi = x_i' * beta_hi
     *
     * If (xb_hi - xb_lo) is too small, set f_i = 0 (Stata uses sqrt(c(epsdouble)))
     *
     * OPTIMIZATION: Parallelize over observations with reduction for total_density.
     * Each observation's density is independent.
     */
    ST_double eps_sqrt = sqrt(2.2e-16);  /* sqrt of machine epsilon */
    ST_double total_density = 0.0;
    ST_double two_h = 2.0 * h;

    #pragma omp parallel for reduction(+:total_density) if(N > 1000)
    for (ST_int i = 0; i < N; i++) {
        /* Compute xb_lo = x_i' * beta_lo and xb_hi = x_i' * beta_hi
         * OPTIMIZATION: Compute both dot products in a single pass with 4x unrolling */
        ST_double xb_lo = 0.0;
        ST_double xb_hi = 0.0;

        ST_int K4 = K - (K & 3);
        ST_int k;
        for (k = 0; k < K4; k += 4) {
            ST_double x0 = X[(size_t)k * N + i];
            ST_double x1 = X[(size_t)(k + 1) * N + i];
            ST_double x2 = X[(size_t)(k + 2) * N + i];
            ST_double x3 = X[(size_t)(k + 3) * N + i];
            xb_lo += x0 * beta_lo[k] + x1 * beta_lo[k + 1] + x2 * beta_lo[k + 2] + x3 * beta_lo[k + 3];
            xb_hi += x0 * beta_hi[k] + x1 * beta_hi[k + 1] + x2 * beta_hi[k + 2] + x3 * beta_hi[k + 3];
        }
        for (; k < K; k++) {
            ST_double x_ik = X[(size_t)k * N + i];
            xb_lo += x_ik * beta_lo[k];
            xb_hi += x_ik * beta_hi[k];
        }

        ST_double diff = xb_hi - xb_lo;

        if (diff > eps_sqrt) {
            obs_density[i] = two_h / diff;
        } else {
            obs_density[i] = 0.0;
        }

        total_density += obs_density[i];
    }

    /* Cleanup */
    free(beta_lo);
    free(beta_hi);

    /* Return average sparsity (reported only; the VCE uses obs_density) */
    ST_double avg_density = total_density / N;
    if (avg_density < 1e-10) avg_density = 1e-10;

    main_debug_log("Fitted per-obs density: avg_density=%.6f, sparsity=%.6f\n",
                   avg_density, 1.0 / avg_density);

    *out_sparsity = 1.0 / avg_density;
    return 0;
}

/*
 * Set to exactly zero the residuals that are zero to working precision:
 * |r_i| <= 8 (K+1) eps (|y_i| + sum_j |x_ij b_j|), the rounding error of
 * y_i - x_i'beta. With tied or discrete data many observations besides the
 * basis lie on the quantile hyperplane; that rounding would otherwise stand
 * in for zero-width percentile gaps or kernel bandwidths, where qreg's
 * exactly computed residuals give zero (and its r(498)). The bound must stay
 * at rounding level: a larger one also zeroes genuine residuals whenever
 * the level of y or of cancelling x_ij b_j terms is large (e.g. an uncentered
 * year trend).
 */
static void cqreg_snap_zero_residuals(ST_double *r, const ST_double *y,
                                      const ST_double *X, const ST_double *beta,
                                      ST_int N, ST_int K)
{
    const ST_double tol = 8.0 * (ST_double)(K + 1) * DBL_EPSILON;
    #pragma omp parallel for if(N > 100000)
    for (ST_int i = 0; i < N; i++) {
        ST_double m = fabs(y[i]);
        for (ST_int j = 0; j < K; j++) {
            m += fabs(X[(size_t)j * N + i] * beta[j]);
        }
        if (fabs(r[i]) <= tol * m) r[i] = 0.0;
    }
}

/* ============================================================================
 * Helper: Store results to Stata
 * ============================================================================ */

static ST_retcode store_results(cqreg_state *state)
{
    ST_int K = state->K;
    ST_int k, j;
    char scalar_name[64];

    /* Store scalar results */
    store_scalar("__cqreg_N", (ST_double)state->N);
    store_scalar("__cqreg_K_keep", (ST_double)(K - 1));  /* Exclude constant for reporting */
    store_scalar("__cqreg_sum_adev", state->sum_adev);
    store_scalar("__cqreg_sum_rdev", state->sum_rdev);
    store_scalar("__cqreg_q_v", state->q_v);
    store_scalar("__cqreg_sparsity", state->sparsity);
    store_scalar("__cqreg_bandwidth", state->bandwidth);
    store_scalar("__cqreg_iterations", (ST_double)state->iterations);
    store_scalar("__cqreg_converged", (ST_double)state->converged);
    if (state->density_method == CQREG_DENSITY_KERNEL) {
        store_scalar("__cqreg_kbwidth", state->kbwidth);
    }

    /* Store coefficients (exclude constant for now, stored separately) */
    for (k = 0; k < K - 1; k++) {
        snprintf(scalar_name, sizeof(scalar_name), "__cqreg_beta_%d", k + 1);
        store_scalar(scalar_name, state->beta[k]);
    }

    /* Store constant */
    store_scalar("__cqreg_cons", state->beta[K - 1]);

    /* Store cluster info if applicable */
    if (state->vce_type == CQREG_VCE_CLUSTER) {
        store_scalar("__cqreg_N_clust", (ST_double)state->num_clusters);
    }

    /* Store VCE matrix */
    for (k = 0; k < K; k++) {
        for (j = 0; j < K; j++) {
            /* Stata matrices are 1-indexed */
            ST_retcode rc = ctools_mat_store("__cqreg_V", k + 1, j + 1,
                                            state->V[(size_t)k * K + j]);
            if (rc) return rc;
        }
    }
    return 0;
}

/* ============================================================================
 * Main Entry Point
 * ============================================================================ */

ST_retcode cqreg_full_regression(const char *args)
{
    (void)args;  /* Unused - arguments parsed via Stata scalars */
    main_debug_open();
    main_debug_log("cqreg_full_regression: ENTRY\n");

    cqreg_state *state = NULL;
    cqreg_ipm_config ipm_config;
    ST_retcode rc = 0;


    /* Timer */
    double time_start = ctools_timer_seconds();


    /* ========================================================================
     * Step 1: Read parameters from Stata scalars
     * ======================================================================== */

    ST_double quantile = read_scalar("__cqreg_quantile", 0.5);
    ST_int K_total = read_scalar_int("__cqreg_K", 2);  /* depvar + indepvars */
    ST_int G = read_scalar_int("__cqreg_G", 0);
    ST_int vce_type = read_scalar_int("__cqreg_vce_type", 0);
    ST_int bw_method = read_scalar_int("__cqreg_bw_method", 0);
    ST_int density_method = read_scalar_int("__cqreg_density_method", 1);  /* Default: fitted */
    ST_double precomputed_bw = read_scalar("__cqreg_bandwidth", -1.0);  /* Pre-computed in Stata for precision */
    ST_int verbose = read_scalar_int("__cqreg_verbose", 0);
    ST_double tolerance = read_scalar("__cqreg_tolerance", 1e-12);
    ST_int maxiter = read_scalar_int("__cqreg_maxiter", 200);
    ST_int nopreprocess = read_scalar_int("__cqreg_nopreprocess", 0);
    /* Kernel method: kernel and invnormal(q+h) - invnormal(q-h) from Stata */
    ST_int kernel = read_scalar_int("__cqreg_kernel", CQREG_KERNEL_EPANECHNIKOV);
    ST_double kfactor = read_scalar("__cqreg_kfactor", -1.0);

    if (G > 0 || maxiter < 1) {
        ctools_error("cqreg", "absorb() is unsupported and maxiter must be positive");
        main_debug_close();
        return 198;
    }

    /* Validate quantile */
    if (quantile <= 0.0 || quantile >= 1.0) {
        ctools_error("cqreg", "Quantile must be between 0 and 1");
        return 198;
    }

    /* qreg defines the residual method for the iid VCE only (r(184)) */
    if (density_method == CQREG_DENSITY_RESIDUAL && vce_type != CQREG_VCE_IID) {
        ctools_error("cqreg", "denmethod(residual) may only be combined with vce(iid)");
        return 184;
    }
    if (density_method == CQREG_DENSITY_KERNEL &&
        (kernel < CQREG_KERNEL_EPANECHNIKOV || kernel > CQREG_KERNEL_TRIANGLE)) {
        ctools_error("cqreg", "invalid kernel");
        return 198;
    }

    /* Offsets into the K x K matrices are ST_int (see CQREG_MAX_K) */
    if (K_total < 1 || K_total > CQREG_MAX_K) {
        ctools_error("cqreg", "at most %d regressors are supported", CQREG_MAX_K - 1);
        return 103;
    }

    /* Build variable index array for filtered loading.
     * Load ALL variables in one shot: depvar + indepvars + cluster var.
     * This avoids a redundant SPI pass for cluster loading later. */
    ST_int depvar_idx = 1;  /* First variable */
    ST_int K_x = K_total - 1;  /* Number of indepvars */
    ST_int K = K_x + 1;  /* Including constant */

    ST_int has_cluster = (vce_type == CQREG_VCE_CLUSTER) ? 1 : 0;
    ST_int nvars_to_load = 1 + K_x + has_cluster;
    int *var_indices = (int *)malloc((size_t)nvars_to_load * sizeof(int));
    if (var_indices == NULL) {
        ctools_error("cqreg", "Memory allocation failed for var_indices");
        return 920;
    }
    var_indices[0] = depvar_idx;
    for (ST_int k = 0; k < K_x; k++) {
        var_indices[1 + k] = k + 2;  /* Variables 2 to K_total */
    }
    if (has_cluster) {
        var_indices[1 + K_x] = K_total + 1;  /* Cluster variable */
    }

    /* Load data with if/in filtering using ctools_data_load.
     * This handles SF_ifobs internally and returns only filtered observations.
     * All vars (depvar, indepvars, cluster) are loaded in one parallel pass. */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    stata_retcode load_rc = ctools_data_load(&filtered, var_indices,
                                                       (size_t)nvars_to_load, 0, 0, CTOOLS_LOAD_CHECK_IF | CTOOLS_LOAD_NO_SORT_ORDER);
    free(var_indices);
    var_indices = NULL;

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        ctools_error("cqreg", "Failed to load data");
        return 920;
    }

    ST_int N = (ST_int)filtered.data.nobs;

    if (N <= 0) {
        ctools_filtered_data_free(&filtered);
        ctools_error("cqreg", "No observations");
        return 2000;
    }

    main_debug_log("cqreg_full_regression: N_filtered=%d\n", N);

    /* Number of variables passed to plugin (for validation if needed) */
    (void)SF_nvars();

    /* Cluster data is already loaded in filtered.data.vars[].
     * Layout: [depvar, x1..xK, cluster]
     */
    main_debug_log("Step 1 done: N=%d, K=%d, q=%.3f\n", N, K, quantile);

    /* ========================================================================
     * Step 2: Create state and load data
     * ======================================================================== */

    main_debug_log("Step 2: Creating state...\n");
    state = cqreg_state_create(N, K);
    if (state == NULL) {
        ctools_filtered_data_free(&filtered);
        ctools_error("cqreg", "Failed to create state");
        return 920;
    }


    state->quantile = quantile;
    state->vce_type = (cqreg_vce_type)vce_type;
    state->bw_method = (cqreg_bw_method)bw_method;
    state->density_method = (cqreg_density_method)density_method;

    main_debug_log("State created. Loading data...\n");

    /* Load data from filtered structure (already loaded with if/in filtering) */
    if (load_data_filtered(state, &filtered, K_x) != 0) {
        main_debug_log("ERROR: load_data_filtered failed\n");
        cqreg_state_free(state);
        ctools_filtered_data_free(&filtered);
        main_debug_close();
        return 198;
    }

    main_debug_log("Data loaded successfully\n");

    /* Record data load time (matches creghdfe - pure data transfer) */
    state->time_load = ctools_timer_seconds() - time_start;

    /* Check for degenerate cases: constant y or constant x columns */
    /* Constant y: all values of y are the same */
    {
        ST_double y_min = state->y[0];
        ST_double y_max = state->y[0];
        for (ST_int i = 1; i < N; i++) {
            if (state->y[i] < y_min) y_min = state->y[i];
            if (state->y[i] > y_max) y_max = state->y[i];
        }
        if (y_max - y_min < 1e-12) {
            ctools_error("cqreg", "Dependent variable has no variation (constant)");
            cqreg_state_free(state);
            ctools_filtered_data_free(&filtered);
            main_debug_close();
            return 198;
        }
    }

    /* Compute sample quantile and raw sum of deviations for pseudo R^2 */
    state->q_v = cqreg_compute_quantile(state->y, N, quantile);
    state->sum_rdev = cqreg_sum_raw_deviations(state->y, N, state->q_v, quantile);
    main_debug_log("Computed q_v=%.6f, sum_rdev=%.4f\n", state->q_v, state->sum_rdev);

    /* Load cluster variable from pre-loaded data (no SPI re-read) */
    if (vce_type == CQREG_VCE_CLUSTER) {
        ST_double *cluster_data = filtered.data.vars[1 + K_x].data.dbl;
        if (load_clusters_from_data(state, cluster_data, N) != 0) {
            cqreg_state_free(state);
            ctools_filtered_data_free(&filtered);
            return 198;
        }
    }

    /* ========================================================================
     * Step 3b: Collinearity detection
     * ======================================================================== */

    ST_int *is_collinear = NULL;

    if (K_x > 0) {
        /* Compute X'X for collinearity check */
        ST_double *xtx = (ST_double *)malloc((size_t)K_x * K_x * sizeof(ST_double));
        is_collinear = (ST_int *)calloc(K_x, sizeof(ST_int));

        if (!xtx || !is_collinear) {
            free(xtx);
            free(is_collinear);
            cqreg_state_free(state);
            ctools_filtered_data_free(&filtered);
            ctools_error("cqreg", "Memory allocation failed for collinearity check");
            main_debug_close();
            return 920;
        }

        for (ST_int i = 0; i < K_x; i++) {
            const ST_double *Xi = &state->X[(size_t)i * N];
            for (ST_int j = i; j < K_x; j++) {
                const ST_double *Xj = &state->X[(size_t)j * N];
                ST_double mean_i = 0.0, mean_j = 0.0;
                for (ST_int row = 0; row < N; row++) {
                    mean_i += Xi[row] / N;
                    mean_j += Xj[row] / N;
                }
                ST_double val = 0.0;
                for (ST_int row = 0; row < N; row++)
                    val += (Xi[row] - mean_i) * (Xj[row] - mean_j);
                xtx[i * K_x + j] = val;
                xtx[j * K_x + i] = val;
            }
        }

        /* Detect numerical collinearity via Cholesky */
        ST_int num_collinear = detect_collinearity(xtx, K_x, is_collinear, verbose);
        free(xtx);

        if (num_collinear < 0) {
            free(is_collinear);
            cqreg_state_free(state);
            ctools_filtered_data_free(&filtered);
            ctools_error("cqreg", "Collinearity detection failed");
            main_debug_close();
            return 920;
        }

        /* Store collinearity flags to Stata */
        char scalar_name[64];
        for (ST_int k = 0; k < K_x; k++) {
            snprintf(scalar_name, sizeof(scalar_name), "__cqreg_collinear_%d", k + 1);
            SF_scal_save(scalar_name, (ST_double)is_collinear[k]);
        }

        main_debug_log("Collinearity: %d of %d variables collinear\n", num_collinear, K_x);

        /* Compact X if any collinear variables found */
        if (num_collinear > 0) {
            ST_int K_keep_x = K_x - num_collinear;
            ST_int K_new = K_keep_x + 1;  /* +1 for constant */

            /* Compact X in-place: shift non-collinear columns to front */
            ST_int dst = 0;
            for (ST_int k = 0; k < K_x; k++) {
                if (!is_collinear[k]) {
                    if (dst != k) {
                        memcpy(&state->X[(size_t)dst * N], &state->X[(size_t)k * N],
                               (size_t)N * sizeof(ST_double));
                    }
                    dst++;
                }
            }
            /* Move constant column to new position */
            memcpy(&state->X[(size_t)K_keep_x * N], &state->X[(size_t)K_x * N],
                   (size_t)N * sizeof(ST_double));

            /* Update state dimensions */
            K_x = K_keep_x;
            K = K_new;
            state->K = K_new;

            /* Reallocate beta and V for new K */
            ctools_aligned_free(state->beta);
            state->beta = (ST_double *)ctools_cacheline_alloc(K_new * sizeof(ST_double));
            if (state->beta == NULL) {
                free(is_collinear);
                cqreg_state_free(state);
                ctools_filtered_data_free(&filtered);
                ctools_error("cqreg", "Memory allocation failed for beta");
                main_debug_close();
                return 920;
            }
            memset(state->beta, 0, K_new * sizeof(ST_double));

            ctools_aligned_free(state->V);
            state->V = (ST_double *)ctools_cacheline_alloc((size_t)K_new * K_new * sizeof(ST_double));
            if (state->V == NULL) {
                free(is_collinear);
                cqreg_state_free(state);
                ctools_filtered_data_free(&filtered);
                ctools_error("cqreg", "Memory allocation failed for V");
                main_debug_close();
                return 920;
            }
            memset(state->V, 0, (size_t)K_new * K_new * sizeof(ST_double));

            main_debug_log("Compacted X: K_x=%d, K=%d\n", K_x, K);
        }
    }

    /* ========================================================================
     * Step 4: Run IPM solver
     * ======================================================================== */


    double ipm_start = ctools_timer_seconds();

    /* Configure IPM */
    cqreg_ipm_config_init(&ipm_config);
    ipm_config.maxiter = maxiter;
    ipm_config.tol_primal = tolerance;
    ipm_config.tol_dual = tolerance;
    ipm_config.tol_gap = tolerance;
    ipm_config.verbose = verbose;
    ipm_config.use_mehrotra = 1;


    main_debug_log("Step 4: Creating IPM state...\n");

    /* The direct FN solver does not use the legacy IPM directions or its
     * N-sized per-thread buffers. Keep those only for experimental preprocessing. */
    ST_int use_preprocess = (nopreprocess == -1);
    state->ipm = use_preprocess ? cqreg_ipm_create(N, K, &ipm_config)
                                : cqreg_ipm_create_lite(N, K, &ipm_config);
    if (state->ipm == NULL) {
        main_debug_log("ERROR: cqreg_ipm_create failed\n");
        ctools_error("cqreg", "Failed to create IPM solver");
        rc = 920;
        goto cleanup;
    }

    main_debug_log("IPM state created. Choosing solver...\n");

    /* Solve quantile regression.
     * The preprocessing algorithm (Chernozhukov et al. 2020) is available but
     * currently DISABLED by default due to performance issues with the current
     * implementation. The direct Frisch-Newton solver is fast enough for most
     * use cases. To enable preprocessing (experimental), set nopreprocess=-1.
     */
    ST_int ipm_result;

    if (use_preprocess) {
        main_debug_log("Using preprocessing solver (N=%d > 5000)\n", N);
        ipm_result = cqreg_preprocess_solve(state->ipm, state->y, state->X, quantile, state->beta);
    } else {
        main_debug_log("Using direct Frisch-Newton solver\n");
        ipm_result = cqreg_fn_solve(state->ipm, state->y, state->X, quantile, state->beta);
    }

    main_debug_log("Solver returned %d\n", ipm_result);

    state->iterations = state->ipm->iterations;
    state->converged = ipm_result > 0 && state->ipm->converged;

    main_debug_log("Getting objective value...\n");
    /* Get objective value */
    state->sum_adev = cqreg_ipm_get_objective(state->ipm, quantile);

    main_debug_log("Getting residuals...\n");
    /* Get residuals */
    cqreg_ipm_get_residuals(state->ipm, state->residuals);

    main_debug_log("IPM phase complete. sum_adev=%.4f\n", state->sum_adev);
    state->time_ipm = ctools_timer_seconds() - ipm_start;

    /* Save IPM sub-timers now, before aux solves can overwrite them */
    double time_ipm_init = state->ipm->time_init;
    double time_ipm_iterate = state->ipm->time_iterate;
    double time_ipm_crossover = state->ipm->time_crossover;




    if (!state->converged) {
        ctools_error("cqreg", "IPM solver stopped without convergence after %d iterations", state->iterations);
        rc = 430;
        goto cleanup;
    }

    /* ========================================================================
     * Step 5: Estimate sparsity and compute VCE
     * ======================================================================== */

    main_debug_log("Step 5: Sparsity and VCE estimation...\n");
    double vce_start = ctools_timer_seconds();
    double time_sparsity = 0.0;
    double time_vce_matrix = 0.0;

    main_debug_log("Creating sparsity state...\n");
    /* Create sparsity state */
    state->sparsity_state = cqreg_sparsity_create(N, quantile, state->bw_method);
    if (state->sparsity_state == NULL) {
        main_debug_log("ERROR: cqreg_sparsity_create failed\n");
        ctools_error("cqreg", "Failed to create sparsity estimator");
        rc = 920;
        goto cleanup;
    }

    main_debug_log("Estimating sparsity/density (method=%d, vce=%d)...\n",
                   state->density_method, state->vce_type);

    /*
     * Estimate sparsity based on density method and VCE type:
     *
     * For IID VCE:
     * - RESIDUAL: Difference quotient of residual percentiles without the LP
     *   basis observations (matches qreg vce(iid, residual))
     * - FITTED: Siddiqui difference quotient at mean(X) (matches qreg default)
     * - KERNEL: kband / mean(K(r/kband)) (matches qreg vce(iid, kernel()))
     *
     * For Robust VCE:
     * - FITTED: Per-observation densities with Powell sandwich
     * - KERNEL: Kernel weights with Powell sandwich
     * (RESIDUAL was rejected above, as in qreg.)
     */
    if (state->density_method == CQREG_DENSITY_FITTED && state->vce_type == CQREG_VCE_ROBUST) {
        /* Fitted method for robust: compute per-observation densities using
         * Stata's approach - difference quotients from QR at τ±h.
         * f_i = (2h) / (Q̂(τ+h|x_i) - Q̂(τ-h|x_i))
         */
        rc = estimate_fitted_per_obs_density(state->obs_density,
                                             state->y, state->X,
                                             N, K, quantile,
                                             state->bw_method,
                                             state->ipm,
                                             &state->bandwidth,
                                             precomputed_bw,
                                             &state->sparsity);
    } else if (state->density_method == CQREG_DENSITY_FITTED) {
        /* Fitted method for IID: use Siddiqui difference quotient method.
         * This fits QR at τ±h and computes sparsity at mean(X).
         * Reuses main IPM state workspace (avoids ~2.8 GB allocation for N=5M).
         * Matches Stata's default vce(iid) "fitted" method exactly.
         */
        rc = estimate_siddiqui_sparsity(state->y, state->X,
                                        N, K, quantile,
                                        state->bw_method,
                                        state->ipm,
                                        precomputed_bw,
                                        &state->sparsity);
        state->bandwidth = (precomputed_bw > 0) ? precomputed_bw : cqreg_compute_bandwidth(N, quantile, state->bw_method);
    } else if (state->density_method == CQREG_DENSITY_KERNEL) {
        /* Kernel method (qreg's VCE_kernel): the bandwidth h is that of all
         * N observations; qreg stops when tau -/+ h leaves (0, 1). */
        ST_double h = (precomputed_bw > 0) ? precomputed_bw : cqreg_compute_bandwidth(N, quantile, state->bw_method);
        state->bandwidth = h;
        if (!(quantile - h > 0.0) || !(quantile + h < 1.0) ||
            SF_is_missing(kfactor) || !(kfactor > 0.0)) {
            cqreg_vce_bandwidth_error();
            rc = 498;
        } else {
            cqreg_snap_zero_residuals(state->residuals, state->y, state->X,
                                      state->beta, N, K);
            state->sparsity_state->kernel = (cqreg_kernel_type)kernel;
            state->sparsity_state->kfactor = kfactor;
            state->sparsity = cqreg_estimate_kernel_sparsity(state->sparsity_state,
                                                             state->residuals,
                                                             state->obs_density);
            state->kbwidth = state->sparsity_state->kbwidth;
            rc = state->sparsity_state->status;
            if (rc == 0 && state->vce_type != CQREG_VCE_IID) {
                /* qreg reports no sparsity for the robust kernel VCE */
                state->sparsity = SV_missval;
            } else if (rc == 0 && !isfinite(state->sparsity)) {
                ctools_error("cqreg", "no residual lies within the kernel bandwidth; "
                             "computation of coefficient standard errors cannot be completed");
                rc = 498;
            }
        }
    } else if (!state->ipm->basis_ok) {
        /* The residual method drops the LP basis observations, which only a
         * crossover to a basic solution identifies */
        ctools_error("cqreg", "denmethod(residual) requires the linear programming "
                     "basis, which the solver did not identify; use denmethod(fitted)");
        rc = 498;
    } else {
        /* Residual method: difference quotient sparsity without the K basis
         * observations (the IPM state has not been reused since the main solve) */
        cqreg_snap_zero_residuals(state->residuals, state->y, state->X,
                                  state->beta, N, K);
        state->sparsity_state->basis = state->ipm->basis;
        state->sparsity_state->n_basis = K;
        state->sparsity_state->full_bandwidth = precomputed_bw;
        state->sparsity = cqreg_estimate_sparsity(state->sparsity_state, state->residuals);
        state->bandwidth = state->sparsity_state->bandwidth;
        rc = state->sparsity_state->status;
    }

    /* Standard errors from a failed density estimate would be arbitrary:
     * stop like qreg (the message has been displayed). */
    if (rc != 0) goto cleanup;

    main_debug_log("Sparsity=%.4f, bandwidth=%.6f\n", state->sparsity, state->bandwidth);

    time_sparsity = ctools_timer_seconds() - vce_start;

    main_debug_log("Computing VCE...\n");
    double vce_matrix_start = ctools_timer_seconds();
    /* Compute VCE */
    if (cqreg_compute_vce(state, state->X) != 0) {
        main_debug_log("ERROR: VCE computation failed\n");
        ctools_error("cqreg", "VCE computation failed");
        rc = 498;
        goto cleanup;
    }
    main_debug_log("VCE computed\n");
    time_vce_matrix = ctools_timer_seconds() - vce_matrix_start;

    state->time_vce = ctools_timer_seconds() - vce_start;


    /* ========================================================================
     * Step 6: Store results to Stata
     * ======================================================================== */

    main_debug_log("Step 6: Storing results to Stata...\n");
    rc = store_results(state);
    if (rc) goto cleanup;
    main_debug_log("Results stored\n");

    state->time_total = ctools_timer_seconds() - time_start;

    /* Store timing scalars for verbose display in Stata */
    SF_scal_save("_cqreg_time_load", state->time_load);
    SF_scal_save("_cqreg_time_ipm", state->time_ipm);
    SF_scal_save("_cqreg_time_ipm_init", time_ipm_init);
    SF_scal_save("_cqreg_time_ipm_iterate", time_ipm_iterate);
    SF_scal_save("_cqreg_time_ipm_crossover", time_ipm_crossover);
    SF_scal_save("_cqreg_time_vce", state->time_vce);
    SF_scal_save("_cqreg_time_sparsity", time_sparsity);
    SF_scal_save("_cqreg_time_vce_matrix", time_vce_matrix);
    SF_scal_save("_cqreg_time_total", state->time_total);
    SF_scal_save("_cqreg_ipm_iterations", (ST_double)state->iterations);
    CTOOLS_SAVE_THREAD_INFO("_cqreg");

    /* ========================================================================
     * Cleanup
     * ======================================================================== */

cleanup:
    main_debug_log("Cleanup: freeing state...\n");
    cqreg_state_free(state);
    main_debug_log("State freed\n");
    ctools_filtered_data_free(&filtered);
    main_debug_log("indepvar_idx freed\n");
    free(is_collinear);

    main_debug_log("cqreg_full_regression: EXIT (rc=%d)\n", rc);
    main_debug_close();
    return rc;
}

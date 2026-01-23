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
#include "cqreg_hdfe.h"
#include "cqreg_linalg.h"
#include "../ctools_error.h"
#include "../ctools_timer.h"
#include "../ctools_config.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
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

static ST_int load_data(cqreg_state *state,
                        ST_int depvar_idx,
                        const ST_int *indepvar_idx,
                        ST_int K_x,
                        ST_int N,
                        ST_int in1, ST_int in2)
{
    ST_int K = K_x + 1;  /* +1 for constant */
    ST_int k, i, idx, obs;
    ST_double val;

    main_debug_log("load_data: N=%d, K_x=%d, in1=%d, in2=%d\n", N, K_x, in1, in2);

    /* Allocate y array directly (owned by state) - with overflow check */
    size_t y_size;
    if (ctools_safe_mul_size((size_t)N, sizeof(ST_double), &y_size) != 0) {
        ctools_error("cqreg", "Overflow computing y array size");
        return -1;
    }
    state->y = (ST_double *)cqreg_aligned_alloc(y_size, CQREG_CACHE_LINE);
    if (state->y == NULL) {
        ctools_error("cqreg", "Failed to allocate y array");
        return -1;
    }
    state->y_owned = 1;

    /* Allocate X matrix: column-major, N rows x K columns - with overflow check */
    size_t x_size;
    if (ctools_safe_alloc_size((size_t)N, (size_t)K, sizeof(ST_double), &x_size) != 0) {
        ctools_error("cqreg", "Overflow computing X matrix size");
        return -1;
    }
    state->X = (ST_double *)cqreg_aligned_alloc(x_size, CQREG_CACHE_LINE);
    if (state->X == NULL) {
        ctools_error("cqreg", "Failed to allocate X matrix");
        return -1;
    }

    /*
     * Load data from Stata, filtering by if condition.
     * Only include observations where SF_ifobs() returns true.
     */
    idx = 0;
    for (obs = in1; obs <= in2; obs++) {
        if (!SF_ifobs(obs)) continue;

        /* Load dependent variable */
        if (SF_vdata(depvar_idx, obs, &val) != 0) {
            ctools_error("cqreg", "Failed to read depvar at obs %d", obs);
            return -1;
        }
        state->y[idx] = val;

        /* Load independent variables */
        for (k = 0; k < K_x; k++) {
            if (SF_vdata(indepvar_idx[k], obs, &val) != 0) {
                ctools_error("cqreg", "Failed to read indepvar %d at obs %d", k+1, obs);
                return -1;
            }
            state->X[k * N + idx] = val;
        }

        idx++;
    }

    main_debug_log("load_data: loaded %d observations (expected %d)\n", idx, N);

    /* Add constant as last column */
    ST_double *const_col = &state->X[K_x * N];
    ST_int N8 = N - (N % 8);
    for (i = 0; i < N8; i += 8) {
        const_col[i]     = 1.0;
        const_col[i + 1] = 1.0;
        const_col[i + 2] = 1.0;
        const_col[i + 3] = 1.0;
        const_col[i + 4] = 1.0;
        const_col[i + 5] = 1.0;
        const_col[i + 6] = 1.0;
        const_col[i + 7] = 1.0;
    }
    for (; i < N; i++) {
        const_col[i] = 1.0;
    }

    main_debug_log("load_data: X matrix populated, constant added\n");

    return 0;
}

/* ============================================================================
 * Helper: Load cluster variable
 * ============================================================================ */

static ST_int load_clusters(cqreg_state *state,
                            ST_int cluster_var_idx,
                            ST_int N,
                            ST_int in1, ST_int in2)
{
    ST_int idx, obs;

    /* Use cqreg_aligned_alloc to match cqreg_aligned_free in cqreg_state_free */
    state->cluster_ids = (ST_int *)cqreg_aligned_alloc(N * sizeof(ST_int), CQREG_CACHE_LINE);
    if (state->cluster_ids == NULL) {
        return -1;
    }

    /* Load cluster IDs, filtering by if condition */
    idx = 0;
    for (obs = in1; obs <= in2; obs++) {
        if (!SF_ifobs(obs)) continue;

        ST_double val;
        if (SF_vdata(cluster_var_idx, obs, &val) != 0) {
            ctools_error("cqreg", "Failed to read cluster var at obs %d", obs);
            return -1;
        }
        state->cluster_ids[idx] = (ST_int)val;
        idx++;
    }

    /* Count unique clusters */
    ST_int max_clusters = N;
    ST_int *seen = (ST_int *)calloc(max_clusters, sizeof(ST_int));
    if (seen == NULL) {
        return -1;
    }

    ST_int num_unique = 0;
    for (idx = 0; idx < N; idx++) {
        ST_int cid = state->cluster_ids[idx];
        ST_int found = 0;
        for (ST_int j = 0; j < num_unique; j++) {
            if (seen[j] == cid) {
                found = 1;
                break;
            }
        }
        if (!found) {
            seen[num_unique++] = cid;
        }
    }
    free(seen);

    state->num_clusters = num_unique;

    return 0;
}

/* ============================================================================
 * Helper: Siddiqui (Fitted) Sparsity Estimation
 * ============================================================================ */

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
static ST_double estimate_siddiqui_sparsity(const ST_double *y,
                                            const ST_double *X,
                                            ST_int N, ST_int K,
                                            ST_double quantile,
                                            cqreg_bw_method bw_method,
                                            const cqreg_ipm_config *ipm_config,
                                            const ST_double *main_beta)
{
    ST_double h = cqreg_compute_bandwidth(N, quantile, bw_method);
    ST_double q_lo = quantile - h;
    ST_double q_hi = quantile + h;

    /* Check bounds */
    if (q_lo <= 0.0 || q_hi >= 1.0) {
        /* Fallback to IQR-based estimate */
        return 1.0;
    }

    /* Allocate temporary arrays for beta coefficients */
    ST_double *beta_lo = (ST_double *)malloc(K * sizeof(ST_double));
    ST_double *beta_hi = (ST_double *)malloc(K * sizeof(ST_double));

    if (beta_lo == NULL || beta_hi == NULL) {
        free(beta_lo);
        free(beta_hi);
        return 1.0;
    }

    /*
     * Configure auxiliary solves for sparsity estimation.
     * These need to match Stata's precision to get accurate VCE.
     * Using same tolerance as main solve ensures sparsity estimate
     * is accurate enough for standard error computation.
     */
    cqreg_ipm_config aux_config;
    cqreg_ipm_config_init(&aux_config);
    aux_config.maxiter = 200;       /* Same as main solve */
    aux_config.tol_primal = 1e-8;   /* Tight tolerance for accurate sparsity */
    aux_config.tol_dual = 1e-8;     /* Tight tolerance for accurate sparsity */
    aux_config.tol_gap = 1e-8;      /* Tight tolerance for accurate sparsity */
    aux_config.verbose = 0;         /* Suppress output for aux solves */
    aux_config.use_mehrotra = ipm_config->use_mehrotra;

    /* Create temporary IPM states for the auxiliary regressions */
    cqreg_ipm_state *ipm_lo = cqreg_ipm_create(N, K, &aux_config);
    cqreg_ipm_state *ipm_hi = cqreg_ipm_create(N, K, &aux_config);

    if (ipm_lo == NULL || ipm_hi == NULL) {
        cqreg_ipm_free(ipm_lo);
        cqreg_ipm_free(ipm_hi);
        free(beta_lo);
        free(beta_hi);
        return 1.0;
    }

    /*
     * OPTIMIZATION: Initialize auxiliary betas from main_beta.
     * Even though τ±h differs from τ, the coefficients are typically close.
     * This provides a better starting point than OLS, reducing iterations.
     */
    if (main_beta != NULL) {
        memcpy(beta_lo, main_beta, K * sizeof(ST_double));
        memcpy(beta_hi, main_beta, K * sizeof(ST_double));
    }

    /* Solve quantile regression at q_lo and q_hi.
     * Run the two solves in parallel if OpenMP is available.
     */
    ST_int rc_lo = 0, rc_hi = 0;

#ifdef _OPENMP
    #pragma omp parallel sections
    {
        #pragma omp section
        { rc_lo = cqreg_fn_solve(ipm_lo, y, X, q_lo, beta_lo); }
        #pragma omp section
        { rc_hi = cqreg_fn_solve(ipm_hi, y, X, q_hi, beta_hi); }
    }
#else
    rc_lo = cqreg_fn_solve(ipm_lo, y, X, q_lo, beta_lo);
    rc_hi = cqreg_fn_solve(ipm_hi, y, X, q_hi, beta_hi);
#endif

    main_debug_log("Siddiqui: q_lo=%.4f, q_hi=%.4f, rc_lo=%d, rc_hi=%d\n",
                   q_lo, q_hi, rc_lo, rc_hi);
    main_debug_log("Siddiqui: beta_lo[0]=%.4f, beta_lo[1]=%.4f\n",
                   beta_lo[0], K > 1 ? beta_lo[1] : 0.0);
    main_debug_log("Siddiqui: beta_hi[0]=%.4f, beta_hi[1]=%.4f\n",
                   beta_hi[0], K > 1 ? beta_hi[1] : 0.0);

    /* Check convergence */
    if (rc_lo <= 0 || rc_hi <= 0) {
        main_debug_log("Siddiqui: FAILED - auxiliary solves did not converge\n");
        cqreg_ipm_free(ipm_lo);
        cqreg_ipm_free(ipm_hi);
        free(beta_lo);
        free(beta_hi);
        return 1.0;
    }

    /* Compute X̄ (mean of each column of X) */
    ST_double *x_bar = (ST_double *)calloc(K, sizeof(ST_double));
    if (x_bar == NULL) {
        cqreg_ipm_free(ipm_lo);
        cqreg_ipm_free(ipm_hi);
        free(beta_lo);
        free(beta_hi);
        return 1.0;
    }

    for (ST_int k = 0; k < K; k++) {
        const ST_double *Xk = &X[k * N];
        ST_double sum = 0.0;
        for (ST_int i = 0; i < N; i++) {
            sum += Xk[i];
        }
        x_bar[k] = sum / (ST_double)N;
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

    /* Ensure positive sparsity */
    if (sparsity < 1e-10) {
        main_debug_log("Siddiqui: sparsity too small, clamping to 1e-10\n");
        sparsity = 1e-10;
    }

    /* Cleanup */
    cqreg_ipm_free(ipm_lo);
    cqreg_ipm_free(ipm_hi);
    free(beta_lo);
    free(beta_hi);
    free(x_bar);

    main_debug_log("Siddiqui: returning sparsity=%.4f\n", sparsity);
    return sparsity;
}

/* ============================================================================
 * Helper: Store results to Stata
 * ============================================================================ */

static void store_results(cqreg_state *state)
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

    /* Store coefficients (exclude constant for now, stored separately) */
    for (k = 0; k < K - 1; k++) {
        snprintf(scalar_name, sizeof(scalar_name), "__cqreg_beta_%d", k + 1);
        store_scalar(scalar_name, state->beta[k]);
    }

    /* Store constant */
    store_scalar("__cqreg_cons", state->beta[K - 1]);

    /* Store HDFE info if applicable */
    if (state->has_hdfe) {
        store_scalar("__cqreg_df_a", (ST_double)state->df_a);
    }

    /* Store cluster info if applicable */
    if (state->vce_type == CQREG_VCE_CLUSTER) {
        store_scalar("__cqreg_N_clust", (ST_double)state->num_clusters);
    }

    /* Store VCE matrix */
    for (k = 0; k < K; k++) {
        for (j = 0; j < K; j++) {
            /* Stata matrices are 1-indexed */
            SF_mat_store("__cqreg_V", k + 1, j + 1, state->V[k * K + j]);
        }
    }
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
    ST_int verbose = read_scalar_int("__cqreg_verbose", 0);
    ST_double tolerance = read_scalar("__cqreg_tolerance", 1e-12);
    ST_int maxiter = read_scalar_int("__cqreg_maxiter", 200);
    ST_int nopreprocess = read_scalar_int("__cqreg_nopreprocess", 0);

    /* Validate quantile */
    if (quantile <= 0.0 || quantile >= 1.0) {
        ctools_error("cqreg", "Quantile must be between 0 and 1");
        return 198;
    }

    /* Get observation range and count valid observations (respecting if condition) */
    ST_int in1 = SF_in1();
    ST_int in2 = SF_in2();

    /* Count observations that pass the if condition */
    ST_int N = 0;
    for (ST_int obs = in1; obs <= in2; obs++) {
        if (SF_ifobs(obs)) N++;
    }

    if (N <= 0) {
        ctools_error("cqreg", "No observations");
        return 2000;
    }

    main_debug_log("cqreg_full_regression: N_range=%d, N_valid=%d\n", N_range, N);

    /* Number of variables passed to plugin (for validation if needed) */
    (void)SF_nvars();

    /* Parse variable indices:
     *   vars 1 to K_total: depvar, indepvars
     *   vars K_total+1 to K_total+G: FE variables
     *   var K_total+G+1: cluster variable (if vce_type == cluster)
     */
    ST_int depvar_idx = 1;  /* First variable */
    ST_int K_x = K_total - 1;  /* Number of indepvars */
    ST_int K = K_x + 1;  /* Including constant */

    ST_int *indepvar_idx = NULL;
    ST_int *fe_var_idx = NULL;
    ST_int cluster_var_idx = 0;

    /* Allocate index arrays */
    if (K_x > 0) {
        indepvar_idx = (ST_int *)malloc(K_x * sizeof(ST_int));
        if (indepvar_idx == NULL) {
            ctools_error("cqreg", "Memory allocation failed");
            return 920;
        }
        for (ST_int k = 0; k < K_x; k++) {
            indepvar_idx[k] = k + 2;  /* Variables 2 to K_total */
        }
    }

    if (G > 0) {
        fe_var_idx = (ST_int *)malloc(G * sizeof(ST_int));
        if (fe_var_idx == NULL) {
            free(indepvar_idx);
            ctools_error("cqreg", "Memory allocation failed");
            return 920;
        }
        for (ST_int g = 0; g < G; g++) {
            fe_var_idx[g] = K_total + g + 1;
        }
    }

    if (vce_type == CQREG_VCE_CLUSTER) {
        cluster_var_idx = K_total + G + 1;
    }



    main_debug_log("Step 1 done: N=%d, K=%d, G=%d, q=%.3f\n", N, K, G, quantile);

    /* ========================================================================
     * Step 2: Create state and load data
     * ======================================================================== */

    main_debug_log("Step 2: Creating state...\n");
    state = cqreg_state_create(N, K);
    if (state == NULL) {
        free(indepvar_idx);
        free(fe_var_idx);
        ctools_error("cqreg", "Failed to create state");
        return 920;
    }


    state->quantile = quantile;
    state->vce_type = (cqreg_vce_type)vce_type;
    state->bw_method = (cqreg_bw_method)bw_method;
    state->density_method = (cqreg_density_method)density_method;

    main_debug_log("State created. Loading data...\n");

    /* Load data (filtering by if condition) */
    if (load_data(state, depvar_idx, indepvar_idx, K_x, N, in1, in2) != 0) {
        main_debug_log("ERROR: load_data failed\n");
        cqreg_state_free(state);
        free(indepvar_idx);
        free(fe_var_idx);
        main_debug_close();
        return 198;
    }

    main_debug_log("Data loaded successfully\n");

    /* Compute sample quantile and raw sum of deviations for pseudo R^2 */
    state->q_v = cqreg_compute_quantile(state->y, N, quantile);
    state->sum_rdev = cqreg_sum_raw_deviations(state->y, N, state->q_v, quantile);
    main_debug_log("Computed q_v=%.6f, sum_rdev=%.4f\n", state->q_v, state->sum_rdev);

    state->time_load = ctools_timer_seconds() - time_start;

    /* Load cluster variable if needed (filtering by if condition) */
    if (vce_type == CQREG_VCE_CLUSTER) {
        if (load_clusters(state, cluster_var_idx, N, in1, in2) != 0) {
            cqreg_state_free(state);
            free(indepvar_idx);
            free(fe_var_idx);
            return 198;
        }
    }

    /* ========================================================================
     * Step 3: HDFE projection (if applicable)
     * ======================================================================== */

    if (G > 0) {
        double hdfe_start = ctools_timer_seconds();

        /* Initialize HDFE with filtered observation count */
        if (cqreg_hdfe_init(state, fe_var_idx, G, N, in1, in2, 10000, 1e-8) != 0) {
            cqreg_state_free(state);
            free(indepvar_idx);
            free(fe_var_idx);
            ctools_error("cqreg", "HDFE initialization failed");
            return 198;
        }

        /* Partial out fixed effects */
        if (cqreg_hdfe_partial_out(state, state->y, state->X, N, K) != 0) {
            cqreg_hdfe_cleanup(state);
            cqreg_state_free(state);
            free(indepvar_idx);
            free(fe_var_idx);
            ctools_error("cqreg", "HDFE partial out failed");
            return 198;
        }

        state->df_a = cqreg_hdfe_get_df_absorbed(state);
        state->time_hdfe = ctools_timer_seconds() - hdfe_start;

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

    /* Create IPM state */
    state->ipm = cqreg_ipm_create(N, K, &ipm_config);
    if (state->ipm == NULL) {
        main_debug_log("ERROR: cqreg_ipm_create failed\n");
        cqreg_hdfe_cleanup(state);
        cqreg_state_free(state);
        free(indepvar_idx);
        free(fe_var_idx);
        ctools_error("cqreg", "Failed to create IPM solver");
        main_debug_close();
        return 920;
    }

    main_debug_log("IPM state created. Choosing solver...\n");

    /* Solve quantile regression.
     * The preprocessing algorithm (Chernozhukov et al. 2020) is available but
     * currently DISABLED by default due to performance issues with the current
     * implementation. The direct Frisch-Newton solver is fast enough for most
     * use cases. To enable preprocessing (experimental), set nopreprocess=-1.
     */
    ST_int ipm_result;
    ST_int use_preprocess = (nopreprocess == -1);  /* Only if explicitly enabled */

    if (use_preprocess) {
        main_debug_log("Using preprocessing solver (N=%d > 5000)\n", N);
        ipm_result = cqreg_preprocess_solve(state->ipm, state->y, state->X, quantile, state->beta);
    } else {
        main_debug_log("Using direct Frisch-Newton solver\n");
        ipm_result = cqreg_fn_solve(state->ipm, state->y, state->X, quantile, state->beta);
    }

    main_debug_log("Solver returned %d\n", ipm_result);

    state->iterations = abs(ipm_result);
    state->converged = (ipm_result > 0) ? 1 : 0;

    main_debug_log("Getting objective value...\n");
    /* Get objective value */
    state->sum_adev = cqreg_ipm_get_objective(state->ipm, quantile);

    main_debug_log("Getting residuals...\n");
    /* Get residuals */
    cqreg_ipm_get_residuals(state->ipm, state->residuals);

    main_debug_log("IPM phase complete. sum_adev=%.4f\n", state->sum_adev);
    state->time_ipm = ctools_timer_seconds() - ipm_start;




    if (!state->converged) {
        ctools_error("cqreg", "IPM solver did not converge in %d iterations", maxiter);
        /* Continue anyway to return partial results */
    }

    /* ========================================================================
     * Step 5: Estimate sparsity and compute VCE
     * ======================================================================== */

    main_debug_log("Step 5: Sparsity and VCE estimation...\n");
    double vce_start = ctools_timer_seconds();

    main_debug_log("Creating sparsity state...\n");
    /* Create sparsity state */
    state->sparsity_state = cqreg_sparsity_create(N, quantile, state->bw_method);
    if (state->sparsity_state == NULL) {
        main_debug_log("ERROR: cqreg_sparsity_create failed\n");
        cqreg_hdfe_cleanup(state);
        cqreg_state_free(state);
        free(indepvar_idx);
        free(fe_var_idx);
        ctools_error("cqreg", "Failed to create sparsity estimator");
        main_debug_close();
        return 920;
    }

    main_debug_log("Estimating sparsity/density (method=%d, vce=%d)...\n",
                   state->density_method, state->vce_type);

    /*
     * Estimate sparsity based on density method and VCE type:
     *
     * For IID VCE:
     * - RESIDUAL: Difference quotient on residuals (matches qreg vce(iid, residual))
     * - FITTED: Kernel density on residuals for scalar sparsity (matches qreg default)
     *
     * For Robust VCE:
     * - RESIDUAL: Difference quotient sparsity with score-based sandwich
     * - FITTED: Per-observation densities with Powell sandwich
     */
    if (state->density_method == CQREG_DENSITY_FITTED && state->vce_type == CQREG_VCE_ROBUST) {
        /* Fitted method for robust: compute per-observation densities */
        state->sparsity = cqreg_estimate_fitted_density(state->obs_density,
                                                        state->y,
                                                        state->residuals,
                                                        N, quantile,
                                                        state->bw_method);
        state->bandwidth = cqreg_compute_bandwidth(N, quantile, state->bw_method);
    } else if (state->density_method == CQREG_DENSITY_FITTED) {
        /* Fitted method for IID: use Siddiqui difference quotient method.
         * This fits QR at τ±h and computes sparsity at mean(X).
         * Uses warm-start from main solve's beta for speed.
         * Matches Stata's default vce(iid) "fitted" method exactly.
         */
        state->sparsity = estimate_siddiqui_sparsity(state->y, state->X,
                                                     N, K, quantile,
                                                     state->bw_method,
                                                     &ipm_config,
                                                     state->beta);
        state->bandwidth = cqreg_compute_bandwidth(N, quantile, state->bw_method);
    } else {
        /* Residual method: use difference quotient sparsity */
        state->sparsity = cqreg_estimate_sparsity(state->sparsity_state, state->residuals);
        state->bandwidth = state->sparsity_state->bandwidth;
    }

    main_debug_log("Sparsity=%.4f, bandwidth=%.6f\n", state->sparsity, state->bandwidth);


    main_debug_log("Computing VCE...\n");
    /* Compute VCE */
    if (cqreg_compute_vce(state, state->X) != 0) {
        main_debug_log("WARNING: VCE computation failed\n");
        ctools_error("cqreg", "VCE computation failed");
        /* Continue with zero VCE */
    }
    main_debug_log("VCE computed\n");

    state->time_vce = ctools_timer_seconds() - vce_start;


    /* ========================================================================
     * Step 6: Store results to Stata
     * ======================================================================== */

    main_debug_log("Step 6: Storing results to Stata...\n");
    store_results(state);
    main_debug_log("Results stored\n");

    state->time_total = ctools_timer_seconds() - time_start;

    /* Store timing scalars for verbose display in Stata */
    SF_scal_save("_cqreg_time_load", state->time_load);
    SF_scal_save("_cqreg_time_hdfe", state->time_hdfe);
    SF_scal_save("_cqreg_time_ipm", state->time_ipm);
    SF_scal_save("_cqreg_time_vce", state->time_vce);
    SF_scal_save("_cqreg_time_total", state->time_total);
    CTOOLS_SAVE_THREAD_INFO("_cqreg");

    /* ========================================================================
     * Cleanup
     * ======================================================================== */

    main_debug_log("Cleanup: freeing state...\n");
    cqreg_hdfe_cleanup(state);
    main_debug_log("HDFE cleaned up\n");
    cqreg_state_free(state);
    main_debug_log("State freed\n");
    free(indepvar_idx);
    main_debug_log("indepvar_idx freed\n");
    free(fe_var_idx);
    main_debug_log("fe_var_idx freed\n");

    main_debug_log("cqreg_full_regression: EXIT (rc=%d)\n", rc);
    main_debug_close();
    return rc;
}

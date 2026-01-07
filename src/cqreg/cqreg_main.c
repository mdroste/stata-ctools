/*
 * cqreg_main.c
 *
 * Full quantile regression command orchestrator.
 * Part of the ctools suite.
 */

#include "cqreg_main.h"
#include "cqreg_types.h"
#include "cqreg_ipm.h"
#include "cqreg_sparsity.h"
#include "cqreg_vce.h"
#include "cqreg_hdfe.h"
#include "cqreg_linalg.h"
#include "../ctools_error.h"
#include "../ctools_timer.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

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
 * Helper: Load data from Stata
 * ============================================================================ */

static ST_int load_data(cqreg_state *state,
                        ST_int depvar_idx,
                        const ST_int *indepvar_idx,
                        ST_int K_x,
                        ST_int in1, ST_int in2)
{
    ST_int N = in2 - in1 + 1;
    ST_int i, k, obs;

    /* Load dependent variable */
    for (i = 0; i < N; i++) {
        obs = in1 + i;
        if (SF_vdata(depvar_idx, obs, &state->y[i]) != 0) {
            ctools_error("cqreg", "Failed to read y at obs %d", obs);
            return -1;
        }
    }

    /* Load independent variables (column-major) */
    for (k = 0; k < K_x; k++) {
        ST_double *col = &state->X[k * N];
        for (i = 0; i < N; i++) {
            obs = in1 + i;
            if (SF_vdata(indepvar_idx[k], obs, &col[i]) != 0) {
                ctools_error("cqreg", "Failed to read X[%d] at obs %d", k, obs);
                return -1;
            }
        }
    }

    /* Add constant as last column */
    ST_double *const_col = &state->X[K_x * N];
    for (i = 0; i < N; i++) {
        const_col[i] = 1.0;
    }

    return 0;
}

/* ============================================================================
 * Helper: Load cluster variable
 * ============================================================================ */

static ST_int load_clusters(cqreg_state *state,
                            ST_int cluster_var_idx,
                            ST_int in1, ST_int in2)
{
    ST_int N = in2 - in1 + 1;
    ST_int i, obs;

    state->cluster_ids = (ST_int *)malloc(N * sizeof(ST_int));
    if (state->cluster_ids == NULL) {
        return -1;
    }

    for (i = 0; i < N; i++) {
        obs = in1 + i;
        ST_double val;
        if (SF_vdata(cluster_var_idx, obs, &val) != 0) {
            ctools_error("cqreg", "Failed to read cluster var at obs %d", obs);
            return -1;
        }
        state->cluster_ids[i] = (ST_int)val;
    }

    /* Count unique clusters */
    ST_int max_clusters = N;
    ST_int *seen = (ST_int *)calloc(max_clusters, sizeof(ST_int));
    if (seen == NULL) {
        return -1;
    }

    ST_int num_unique = 0;
    for (i = 0; i < N; i++) {
        ST_int cid = state->cluster_ids[i];
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
    ST_int verbose = read_scalar_int("__cqreg_verbose", 0);
    ST_double tolerance = read_scalar("__cqreg_tolerance", 1e-8);
    ST_int maxiter = read_scalar_int("__cqreg_maxiter", 200);

    /* Validate quantile */
    if (quantile <= 0.0 || quantile >= 1.0) {
        ctools_error("cqreg", "Quantile must be between 0 and 1");
        return 198;
    }

    /* Get observation range */
    ST_int in1 = SF_in1();
    ST_int in2 = SF_in2();
    ST_int N = in2 - in1 + 1;

    if (N <= 0) {
        ctools_error("cqreg", "No observations");
        return 2000;
    }

    /* Number of variables passed to plugin */
    ST_int nvars = SF_nvars();

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


    if (verbose) {
        ctools_msg("cqreg", "N=%d, K=%d, G=%d, quantile=%.3f", N, K, G, quantile);
    }

    /* ========================================================================
     * Step 2: Create state and load data
     * ======================================================================== */

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

    /* Load data */
    if (load_data(state, depvar_idx, indepvar_idx, K_x, in1, in2) != 0) {
        cqreg_state_free(state);
        free(indepvar_idx);
        free(fe_var_idx);
        return 198;
    }

    state->time_load = ctools_timer_seconds() - time_start;
    double time_after_load = ctools_timer_seconds();



    /* Load cluster variable if needed */
    if (vce_type == CQREG_VCE_CLUSTER) {
        if (load_clusters(state, cluster_var_idx, in1, in2) != 0) {
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

        /* Initialize HDFE */
        if (cqreg_hdfe_init(state, fe_var_idx, G, in1, in2, 10000, 1e-8) != 0) {
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

        if (verbose) {
            ctools_msg("cqreg", "HDFE: df_absorbed=%d, time=%.3fs",
                       state->df_a, state->time_hdfe);
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


    /* Create IPM state */
    state->ipm = cqreg_ipm_create(N, K, &ipm_config);
    if (state->ipm == NULL) {
        cqreg_hdfe_cleanup(state);
        cqreg_state_free(state);
        free(indepvar_idx);
        free(fe_var_idx);
        ctools_error("cqreg", "Failed to create IPM solver");
        return 920;
    }


    /* Solve quantile regression using smoothed IPM
     * Note: Preprocessing disabled - smoothed IPM gives approximate solutions
     * which don't work well with LP-based preprocessing methods.
     * For large N, consider implementing Frisch-Newton exact LP solver.
     */
    ST_int ipm_result = cqreg_ipm_solve(state->ipm, state->y, state->X,
                                         quantile, state->beta);

    state->iterations = abs(ipm_result);
    state->converged = (ipm_result > 0) ? 1 : 0;

    /* Get objective value */
    state->sum_adev = cqreg_ipm_get_objective(state->ipm, quantile);

    /* Get residuals */
    cqreg_ipm_get_residuals(state->ipm, state->residuals);

    state->time_ipm = ctools_timer_seconds() - ipm_start;



    if (verbose) {
        ctools_msg("cqreg", "IPM: iterations=%d, converged=%d, objective=%.4f, time=%.3fs",
                   state->iterations, state->converged, state->sum_adev, state->time_ipm);
    }

    if (!state->converged) {
        ctools_error("cqreg", "IPM solver did not converge in %d iterations", maxiter);
        /* Continue anyway to return partial results */
    }

    /* ========================================================================
     * Step 5: Estimate sparsity and compute VCE
     * ======================================================================== */

    double vce_start = ctools_timer_seconds();

    /* Create sparsity state */
    state->sparsity_state = cqreg_sparsity_create(N, quantile, state->bw_method);
    if (state->sparsity_state == NULL) {
        cqreg_hdfe_cleanup(state);
        cqreg_state_free(state);
        free(indepvar_idx);
        free(fe_var_idx);
        ctools_error("cqreg", "Failed to create sparsity estimator");
        return 920;
    }

    /* Estimate sparsity */
    state->sparsity = cqreg_estimate_sparsity(state->sparsity_state, state->residuals);
    state->bandwidth = state->sparsity_state->bandwidth;

    if (verbose) {
        ctools_msg("cqreg", "Sparsity: %.4f, bandwidth: %.6f", state->sparsity, state->bandwidth);
    }

    /* Compute VCE */
    if (cqreg_compute_vce(state, state->X) != 0) {
        ctools_error("cqreg", "VCE computation failed");
        /* Continue with zero VCE */
    }

    state->time_vce = ctools_timer_seconds() - vce_start;


    /* ========================================================================
     * Step 6: Store results to Stata
     * ======================================================================== */

    store_results(state);

    state->time_total = ctools_timer_seconds() - time_start;

    if (verbose) {
        ctools_msg("cqreg", "Total time: %.3fs (load: %.3f, hdfe: %.3f, ipm: %.3f, vce: %.3f)",
                   state->time_total, state->time_load, state->time_hdfe,
                   state->time_ipm, state->time_vce);
    }

    /* ========================================================================
     * Cleanup
     * ======================================================================== */

    cqreg_hdfe_cleanup(state);
    cqreg_state_free(state);
    free(indepvar_idx);
    free(fe_var_idx);

    return rc;
}

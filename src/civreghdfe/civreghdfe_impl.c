/*
    civreghdfe_impl.c
    Instrumental Variables Regression with High-Dimensional Fixed Effects

    Implements 2SLS/IV regression with HDFE absorption. Reuses the creghdfe
    infrastructure for HDFE demeaning, singleton detection, and VCE computation.

    Algorithm:
    1. Read data: y, X_endog, X_exog, Z (instruments), FE vars, weights, cluster
    2. Detect and remove singletons
    3. Partial out FEs from all variables using CG solver
    4. First stage: Regress each endogenous var on [exogenous + instruments]
    5. Second stage: Regress y on [exogenous + predicted endogenous]
    6. Compute VCE (corrected for 2SLS)
    7. Store results to Stata
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "civreghdfe_impl.h"
#include "../ctools_config.h"
#include "../ctools_runtime.h"
#include "../ctools_ols.h"
#include "../ctools_types.h"  /* For ctools_data_load */
#include "../ctools_spi.h"  /* Error-checking SPI wrappers */
#include "civreghdfe_matrix.h"
#include "civreghdfe_estimate.h"
#include "civreghdfe_vce.h"
#include "civreghdfe_tests.h"
#include "../ctools_hdfe_utils.h"
#include "../ctools_unroll.h"

/* Shared OLS functions */
#define cholesky ctools_cholesky
#define invert_from_cholesky ctools_invert_from_cholesky

/* Debug flag */
#define CIVREGHDFE_DEBUG 0


/*
    Full IV regression with HDFE.

    Variable layout from Stata:
    [y, X_endog_1, ..., X_endog_Ke, X_exog_1, ..., X_exog_Kx, Z_1, ..., Z_Kz, FE_1, ..., FE_G, cluster, weight]

    Scalars from Stata:
    - __civreghdfe_K_endog: number of endogenous regressors
    - __civreghdfe_K_exog: number of exogenous regressors
    - __civreghdfe_K_iv: number of instruments (including exogenous)
    - __civreghdfe_G: number of FE factors
    - __civreghdfe_has_cluster, __civreghdfe_has_weights, __civreghdfe_weight_type
    - __civreghdfe_vce_type: 0=unadjusted, 1=robust, 2=cluster
    - __civreghdfe_maxiter, __civreghdfe_tolerance, __civreghdfe_verbose
*/
static ST_retcode do_iv_regression(void)
{
    /* Start total timing */
    double t_total_start = ctools_timer_seconds();
    double t_load = 0, t_extract = 0, t_singleton = 0, t_hdfe_setup = 0, t_fwl = 0, t_partial_out = 0;
    double t_postproc = 0, t_estimate = 0, t_store = 0;

    ST_int in1 = SF_in1();
    ST_int in2 = SF_in2();
    ST_int N_total = in2 - in1 + 1;

    /* Read configuration scalars */
    ST_double dval;
    ST_int K_endog, K_exog, K_iv, G;
    ST_int has_cluster, has_weights, weight_type;
    ST_int vce_type, maxiter, verbose;
    ST_double tolerance;

    SF_scal_use("__civreghdfe_K_endog", &dval); K_endog = (ST_int)dval;
    SF_scal_use("__civreghdfe_K_exog", &dval); K_exog = (ST_int)dval;
    SF_scal_use("__civreghdfe_K_iv", &dval); K_iv = (ST_int)dval;
    SF_scal_use("__civreghdfe_G", &dval); G = (ST_int)dval;
    SF_scal_use("__civreghdfe_has_cluster", &dval); has_cluster = (ST_int)dval;
    ST_int has_cluster2 = 0;
    SF_scal_use("__civreghdfe_has_cluster2", &dval); has_cluster2 = (ST_int)dval;
    SF_scal_use("__civreghdfe_has_weights", &dval); has_weights = (ST_int)dval;
    SF_scal_use("__civreghdfe_weight_type", &dval); weight_type = (ST_int)dval;
    SF_scal_use("__civreghdfe_vce_type", &dval); vce_type = (ST_int)dval;
    SF_scal_use("__civreghdfe_maxiter", &dval); maxiter = (ST_int)dval;
    SF_scal_use("__civreghdfe_tolerance", &dval); tolerance = dval;
    SF_scal_use("__civreghdfe_verbose", &dval); verbose = (ST_int)dval;
    ST_int nested_fe_index;
    SF_scal_use("__civreghdfe_nested_fe_index", &dval); nested_fe_index = (ST_int)dval;

    /* New scalars for estimation method */
    ST_int est_method = 0;
    ST_double kclass_user = 0, fuller_alpha = 0;
    SF_scal_use("__civreghdfe_est_method", &dval); est_method = (ST_int)dval;
    SF_scal_use("__civreghdfe_kclass", &dval); kclass_user = dval;
    SF_scal_use("__civreghdfe_fuller", &dval); fuller_alpha = dval;

    /* HAC parameters */
    ST_int kernel_type = 0, bw = 0, kiefer = 0, dkraay_T = 0, hac_panel = 0;
    SF_scal_use("__civreghdfe_kernel", &dval); kernel_type = (ST_int)dval;
    SF_scal_use("__civreghdfe_bw", &dval); bw = (ST_int)dval;
    SF_scal_use("__civreghdfe_kiefer", &dval); kiefer = (ST_int)dval;
    SF_scal_use("__civreghdfe_dkraay_T", &dval); dkraay_T = (ST_int)dval;
    SF_scal_use("__civreghdfe_hac_panel", &dval); hac_panel = (ST_int)dval;

    /* DOF adjustment parameters */
    ST_int dofminus = 0, sdofminus_opt = 0, nopartialsmall = 0, center = 0;
    SF_scal_use("__civreghdfe_dofminus", &dval); dofminus = (ST_int)dval;
    SF_scal_use("__civreghdfe_sdofminus", &dval); sdofminus_opt = (ST_int)dval;
    SF_scal_use("__civreghdfe_nopartialsmall", &dval); nopartialsmall = (ST_int)dval;
    SF_scal_use("__civreghdfe_center", &dval); center = (ST_int)dval;

    /* Calculate variable positions */
    /* Layout: [y, X_endog (Ke), X_exog (Kx), Z (Kz), FE (G), cluster?, cluster2?, weight?] */
    ST_int var_y_idx = 0;  /* Index in loaded_data.vars[] */
    ST_int var_endog_start_idx = 1;
    ST_int var_exog_start_idx = var_endog_start_idx + K_endog;
    ST_int var_iv_start_idx = var_exog_start_idx + K_exog;
    ST_int var_fe_start_idx = var_iv_start_idx + K_iv;
    ST_int var_cluster_idx = has_cluster ? (var_fe_start_idx + G) : -1;
    ST_int var_cluster2_idx = has_cluster2 ? (var_fe_start_idx + G + (has_cluster ? 1 : 0)) : -1;
    ST_int var_weight_idx = has_weights ? (var_fe_start_idx + G + (has_cluster ? 1 : 0) + (has_cluster2 ? 1 : 0)) : -1;

    ST_int K_total = K_exog + K_endog;  /* Total regressors */

    /* Start data load timing */
    double t_load_start = ctools_timer_seconds();

    /* ================================================================
     * PARALLEL DATA LOADING using ctools_data_load
     * This handles if/in filtering at load time, loading only filtered observations.
     * ================================================================ */

    /* Calculate total variables to load */
    ST_int total_vars = 1 + K_endog + K_exog + K_iv + G +
                        (has_cluster ? 1 : 0) + (has_cluster2 ? 1 : 0) + (has_weights ? 1 : 0);

    /* Build variable indices array (1-based for Stata) */
    int *var_indices = (int *)malloc(total_vars * sizeof(int));
    if (!var_indices) {
        SF_error("civreghdfe: Memory allocation failed\n");
        return 920;
    }
    for (ST_int i = 0; i < total_vars; i++) {
        var_indices[i] = i + 1;  /* 1-based Stata variable indices */
    }

    /* Load all variables in parallel with if/in filtering */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);

    stata_retcode load_rc = ctools_data_load(&filtered, var_indices, total_vars, 0, 0, 0);
    free(var_indices);

    /* Capture pure SPI data load time */
    t_load = ctools_timer_seconds() - t_load_start;

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        SF_error("civreghdfe: Parallel data load failed\n");
        return 920;
    }

    /* Start extraction timer */
    double t_extract_start = ctools_timer_seconds();

    /* Get filtered observation count - data is already filtered */
    N_total = (ST_int)filtered.data.nobs;

    if (N_total == 0) {
        SF_error("civreghdfe: No observations selected\n");
        ctools_filtered_data_free(&filtered);
        return 2001;
    }

    /* Allocate output arrays (with overflow-safe multiplication) */
    ST_double *y = (ST_double *)ctools_safe_malloc2((size_t)N_total, sizeof(ST_double));
    ST_double *X_endog = (K_endog > 0) ? (ST_double *)ctools_safe_malloc3((size_t)N_total, (size_t)K_endog, sizeof(ST_double)) : NULL;
    ST_double *X_exog = (K_exog > 0) ? (ST_double *)ctools_safe_malloc3((size_t)N_total, (size_t)K_exog, sizeof(ST_double)) : NULL;
    ST_double *Z = (ST_double *)ctools_safe_malloc3((size_t)N_total, (size_t)K_iv, sizeof(ST_double));
    ST_double *weights = has_weights ? (ST_double *)ctools_safe_malloc2((size_t)N_total, sizeof(ST_double)) : NULL;
    ST_int *cluster_ids = has_cluster ? (ST_int *)ctools_safe_malloc2((size_t)N_total, sizeof(ST_int)) : NULL;
    ST_int *cluster2_ids = has_cluster2 ? (ST_int *)ctools_safe_malloc2((size_t)N_total, sizeof(ST_int)) : NULL;

    /* Allocate FE level arrays (with overflow-safe multiplication) */
    ST_int **fe_levels = (ST_int **)ctools_safe_malloc2((size_t)G, sizeof(ST_int *));
    int fe_alloc_failed = 0;
    if (fe_levels) {
        for (ST_int g = 0; g < G; g++) {
            fe_levels[g] = (ST_int *)ctools_safe_malloc2((size_t)N_total, sizeof(ST_int));
            if (!fe_levels[g]) fe_alloc_failed = 1;
        }
    } else {
        fe_alloc_failed = 1;
    }

    if (!y || !Z || (K_endog > 0 && !X_endog) || (K_exog > 0 && !X_exog) ||
        (has_weights && !weights) || (has_cluster && !cluster_ids) ||
        (has_cluster2 && !cluster2_ids) || fe_alloc_failed) {
        SF_error("civreghdfe: Memory allocation failed\n");
        free(y); free(X_endog); free(X_exog); free(Z);
        free(weights); free(cluster_ids); free(cluster2_ids);
        if (fe_levels) {
            for (ST_int g = 0; g < G; g++) free(fe_levels[g]);
            free(fe_levels);
        }
        ctools_filtered_data_free(&filtered);
        return 920;
    }

    /* Extract data from filtered.data - data is already filtered, use direct copy */
    ST_int i;
    memcpy(y, filtered.data.vars[var_y_idx].data.dbl, N_total * sizeof(ST_double));

    for (ST_int k = 0; k < K_endog; k++) {
        memcpy(X_endog + k * N_total, filtered.data.vars[var_endog_start_idx + k].data.dbl,
               N_total * sizeof(ST_double));
    }
    for (ST_int k = 0; k < K_exog; k++) {
        memcpy(X_exog + k * N_total, filtered.data.vars[var_exog_start_idx + k].data.dbl,
               N_total * sizeof(ST_double));
    }
    for (ST_int k = 0; k < K_iv; k++) {
        memcpy(Z + k * N_total, filtered.data.vars[var_iv_start_idx + k].data.dbl,
               N_total * sizeof(ST_double));
    }
    for (ST_int g = 0; g < G; g++) {
        double *src = filtered.data.vars[var_fe_start_idx + g].data.dbl;
        for (i = 0; i < N_total; i++) {
            /* Use -1 as sentinel for missing FE values to avoid undefined cast behavior */
            fe_levels[g][i] = SF_is_missing(src[i]) ? -1 : (ST_int)src[i];
        }
    }
    if (has_cluster) {
        double *src = filtered.data.vars[var_cluster_idx].data.dbl;
        for (i = 0; i < N_total; i++) {
            /* Use -1 as sentinel for missing cluster values */
            cluster_ids[i] = SF_is_missing(src[i]) ? -1 : (ST_int)src[i];
        }
    }
    if (has_cluster2) {
        double *src = filtered.data.vars[var_cluster2_idx].data.dbl;
        for (i = 0; i < N_total; i++) {
            /* Use -1 as sentinel for missing cluster2 values */
            cluster2_ids[i] = SF_is_missing(src[i]) ? -1 : (ST_int)src[i];
        }
    }
    if (has_weights) {
        memcpy(weights, filtered.data.vars[var_weight_idx].data.dbl, N_total * sizeof(ST_double));
    }

    /* Free filtered data - we've extracted what we need */
    ctools_filtered_data_free(&filtered);

    /* End extraction, start missing value check timing */
    t_extract = ctools_timer_seconds() - t_extract_start;
    double t_missing_start = ctools_timer_seconds();

    /* Drop observations with missing values */
    ST_int *valid_mask = (ST_int *)calloc(N_total, sizeof(ST_int));
    ST_int N_valid = 0;

    for (i = 0; i < N_total; i++) {
        int is_valid = 1;

        /* Check y */
        if (SF_is_missing(y[i])) is_valid = 0;

        /* Check X_endog */
        for (ST_int k = 0; k < K_endog && is_valid; k++) {
            if (SF_is_missing(X_endog[k * N_total + i])) is_valid = 0;
        }

        /* Check X_exog */
        for (ST_int k = 0; k < K_exog && is_valid; k++) {
            if (SF_is_missing(X_exog[k * N_total + i])) is_valid = 0;
        }

        /* Check Z */
        for (ST_int k = 0; k < K_iv && is_valid; k++) {
            if (SF_is_missing(Z[k * N_total + i])) is_valid = 0;
        }

        /* Check FE variables for missing values (marked as -1 during extraction) */
        for (ST_int g = 0; g < G && is_valid; g++) {
            if (fe_levels[g][i] < 0) is_valid = 0;
        }

        /* Check weights */
        if (has_weights && is_valid) {
            if (SF_is_missing(weights[i]) || weights[i] <= 0) is_valid = 0;
        }

        /* Check cluster variables for missing values (marked as -1 during extraction) */
        if (has_cluster && is_valid) {
            if (cluster_ids[i] < 0) is_valid = 0;
        }
        if (has_cluster2 && is_valid) {
            if (cluster2_ids[i] < 0) is_valid = 0;
        }

        valid_mask[i] = is_valid;
        if (is_valid) N_valid++;
    }

    if (N_valid < K_total + K_iv + 1) {
        SF_error("civreghdfe: Insufficient observations\n");
        free(y); free(X_endog); free(X_exog); free(Z);
        free(weights); free(cluster_ids); free(cluster2_ids); free(valid_mask);
        for (ST_int g = 0; g < G; g++) free(fe_levels[g]);
        free(fe_levels);
        return 2001;
    }

    /* Fused compaction: compact FE levels first (needed for singleton detection),
     * then run singleton detection, then compact all data in one pass.
     * This avoids two separate O(N×K) data compaction passes. */

    /* Step 1: Compact only FE levels for singleton detection */
    ST_int **fe_levels_c = (ST_int **)malloc((size_t)G * sizeof(ST_int *));
    int fe_c_alloc_failed = 0;
    if (fe_levels_c) {
        for (ST_int g = 0; g < G; g++) {
            fe_levels_c[g] = (ST_int *)malloc((size_t)N_valid * sizeof(ST_int));
            if (!fe_levels_c[g]) fe_c_alloc_failed = 1;
        }
    } else {
        fe_c_alloc_failed = 1;
    }
    if (fe_c_alloc_failed) {
        SF_error("civreghdfe: Memory allocation failed for FE compaction\n");
        if (fe_levels_c) {
            for (ST_int g = 0; g < G; g++) if (fe_levels_c[g]) free(fe_levels_c[g]);
            free(fe_levels_c);
        }
        free(y); free(X_endog); free(X_exog); free(Z);
        free(weights); free(cluster_ids); free(cluster2_ids); free(valid_mask);
        for (ST_int g = 0; g < G; g++) free(fe_levels[g]);
        free(fe_levels);
        return 920;
    }
    {
        ST_int idx = 0;
        for (ST_int ii = 0; ii < N_total; ii++) {
            if (!valid_mask[ii]) continue;
            for (ST_int g = 0; g < G; g++)
                fe_levels_c[g][idx] = fe_levels[g][ii];
            idx++;
        }
    }

    /* End missing value check/compact timing, start singleton timing */
    double t_missing = ctools_timer_seconds() - t_missing_start;
    double t_singleton_start = ctools_timer_seconds();

    /* Step 2: Singleton detection on compacted FE levels */
    ST_int *singleton_mask = (ST_int *)malloc((size_t)N_valid * sizeof(ST_int));
    if (!singleton_mask) {
        SF_error("civreghdfe: Memory allocation failed for singleton mask\n");
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
        free(y); free(X_endog); free(X_exog); free(Z);
        free(weights); free(cluster_ids); free(cluster2_ids); free(valid_mask);
        for (ST_int g = 0; g < G; g++) free(fe_levels[g]);
        free(fe_levels);
        return 920;
    }
    ST_int num_singletons_total = ctools_remove_singletons(
        fe_levels_c, G, N_valid, singleton_mask, 100, verbose
    );

    /* Step 3: Build combined index mapping original obs → final position.
     * Combine valid_mask and singleton_mask into one compaction pass. */
    ST_int N = N_valid - num_singletons_total;

    /* Allocate final arrays */
    ST_double *y_c = (ST_double *)malloc((size_t)N * sizeof(ST_double));
    ST_double *X_endog_c = (K_endog > 0) ? (ST_double *)malloc((size_t)N * K_endog * sizeof(ST_double)) : NULL;
    ST_double *X_exog_c = (K_exog > 0) ? (ST_double *)malloc((size_t)N * K_exog * sizeof(ST_double)) : NULL;
    ST_double *Z_c = (ST_double *)malloc((size_t)N * K_iv * sizeof(ST_double));
    ST_double *weights_c = has_weights ? (ST_double *)malloc((size_t)N * sizeof(ST_double)) : NULL;
    ST_int *cluster_ids_c = has_cluster ? (ST_int *)malloc((size_t)N * sizeof(ST_int)) : NULL;
    ST_int *cluster2_ids_c = has_cluster2 ? (ST_int *)malloc((size_t)N * sizeof(ST_int)) : NULL;

    if (!y_c || !Z_c || (K_endog > 0 && !X_endog_c) || (K_exog > 0 && !X_exog_c) ||
        (has_weights && !weights_c) || (has_cluster && !cluster_ids_c) ||
        (has_cluster2 && !cluster2_ids_c)) {
        SF_error("civreghdfe: Memory allocation failed for compacted arrays\n");
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        free(singleton_mask);
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
        free(y); free(X_endog); free(X_exog); free(Z);
        free(weights); free(cluster_ids); free(cluster2_ids); free(valid_mask);
        for (ST_int g = 0; g < G; g++) free(fe_levels[g]);
        free(fe_levels);
        return 920;
    }

    /* Single fused compaction pass: skip invalid AND singleton observations */
    {
        ST_int valid_idx = 0;  /* Index into singleton_mask (N_valid-sized) */
        ST_int out_idx = 0;    /* Index into final arrays (N-sized) */
        for (ST_int ii = 0; ii < N_total; ii++) {
            if (!valid_mask[ii]) continue;
            /* valid_idx tracks position in the N_valid-sized arrays */
            if (!singleton_mask[valid_idx]) {
                valid_idx++;
                continue;
            }
            y_c[out_idx] = y[ii];
            for (ST_int k = 0; k < K_endog; k++)
                X_endog_c[k * N + out_idx] = X_endog[k * N_total + ii];
            for (ST_int k = 0; k < K_exog; k++)
                X_exog_c[k * N + out_idx] = X_exog[k * N_total + ii];
            for (ST_int k = 0; k < K_iv; k++)
                Z_c[k * N + out_idx] = Z[k * N_total + ii];
            if (has_weights) weights_c[out_idx] = weights[ii];
            if (has_cluster) cluster_ids_c[out_idx] = cluster_ids[ii];
            if (has_cluster2) cluster2_ids_c[out_idx] = cluster2_ids[ii];
            out_idx++;
            valid_idx++;
        }
    }

    /* Compact FE levels from the already-compacted fe_levels_c using singleton_mask */
    if (num_singletons_total > 0) {
        for (ST_int g = 0; g < G; g++) {
            ST_int *new_levels = (ST_int *)malloc((size_t)N * sizeof(ST_int));
            if (new_levels) {
                ctools_compact_array_int(fe_levels_c[g], new_levels, singleton_mask, N_valid, N);
                free(fe_levels_c[g]);
                fe_levels_c[g] = new_levels;
            }
        }
    }

    free(singleton_mask);

    /* Free original arrays */
    free(y); free(X_endog); free(X_exog); free(Z);
    free(weights); free(cluster_ids); free(cluster2_ids); free(valid_mask);
    for (ST_int g = 0; g < G; g++) free(fe_levels[g]);
    free(fe_levels);

    /* Check we still have enough observations */
    if (N < K_total + K_iv + 1) {
        SF_error("civreghdfe: Insufficient observations after singleton removal\n");
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
        return 2001;
    }

    /* End singleton timing, start HDFE setup timing */
    t_singleton = ctools_timer_seconds() - t_singleton_start;
    double t_hdfe_setup_start = ctools_timer_seconds();

    /* Initialize HDFE state (reuse creghdfe infrastructure) */
    /* This requires setting up FE_Factor structures */
    HDFE_State *state = (HDFE_State *)calloc(1, sizeof(HDFE_State));
    state->G = G;
    state->N = N;
    state->K = 1 + K_endog + K_exog + K_iv;  /* Total columns to demean */
    state->in1 = 1;
    state->in2 = N;
    state->has_weights = has_weights;
    state->weight_type = weight_type;
    state->weights = weights_c;
    state->maxiter = maxiter;
    state->tolerance = tolerance;
    state->verbose = verbose;

    /* Set up factors */
    state->factors = (FE_Factor *)calloc(G, sizeof(FE_Factor));
    for (ST_int g = 0; g < G; g++) {
        /* Find max level value for remapping */
        ST_int max_level = 0;
        for (ST_int i = 0; i < N; i++) {
            if (fe_levels_c[g][i] > max_level) max_level = fe_levels_c[g][i];
        }

        /* Remap FE levels to contiguous 1-based indices
           The CG solver expects levels to be 1, 2, 3, ..., num_levels
           so it can use levels[i] - 1 as array indices */
        ST_int *remap = (ST_int *)calloc(max_level + 1, sizeof(ST_int));
        if (!remap) {
            SF_error("civreghdfe: Memory allocation failed for level remap\n");
            /* Clean up the g factors already initialized; free remaining fe_levels */
            state->G = g;  /* Only cleanup factors 0..g-1 */
            ctools_hdfe_state_cleanup(state);
            free(state);
            free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
            free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
            for (ST_int fg = g; fg < G; fg++) free(fe_levels_c[fg]);
            free(fe_levels_c);
            return 920;
        }

        /* First pass: assign contiguous level IDs starting from 1 */
        ST_int next_level = 1;
        for (ST_int i = 0; i < N; i++) {
            ST_int old_level = fe_levels_c[g][i];
            if (remap[old_level] == 0) {
                remap[old_level] = next_level++;
            }
        }

        ST_int num_levels = next_level - 1;  /* Total unique levels */

        /* Allocate factor arrays with correct size */
        state->factors[g].has_intercept = 1;
        state->factors[g].levels = fe_levels_c[g];  /* Will be remapped in place */
        state->factors[g].num_levels = num_levels;
        state->factors[g].max_level = num_levels;  /* After remapping, max_level == num_levels */
        state->factors[g].counts = (ST_double *)calloc(num_levels, sizeof(ST_double));
        state->factors[g].weighted_counts = has_weights ? (ST_double *)calloc(num_levels, sizeof(ST_double)) : NULL;
        state->factors[g].means = (ST_double *)calloc(num_levels, sizeof(ST_double));

        /* Second pass: remap levels in place and count */
        for (ST_int i = 0; i < N; i++) {
            ST_int old_level = fe_levels_c[g][i];
            ST_int new_level = remap[old_level];  /* 1-based */
            fe_levels_c[g][i] = new_level;        /* Remap in place */

            state->factors[g].counts[new_level - 1] += 1.0;  /* 0-based index */
            if (has_weights) {
                state->factors[g].weighted_counts[new_level - 1] += weights_c[i];
            }
        }

        free(remap);

        /* Build sorted permutation for cache-friendly scatter-gather */
        ctools_build_sorted_permutation(&state->factors[g], N);

    }

    /* Allocate inv_counts, inv_weighted_counts, and CG solver buffers */
    ST_int num_threads = omp_get_max_threads();
    if (num_threads > 8) num_threads = 8;
    state->num_threads = num_threads;

    /* Max columns partialled out: y + endogenous + exogenous + instruments */
    ST_int max_partial_cols = 1 + K_endog + K_exog + K_iv;
    if (ctools_hdfe_alloc_buffers(state, 1, max_partial_cols) != 0) {
        SF_error("civreghdfe: Memory allocation failed for HDFE buffers\n");
        ctools_hdfe_state_cleanup(state);
        free(state);
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        free(fe_levels_c);
        return 920;
    }

    state->factors_initialized = 1;

    /* Set global state for creghdfe solver functions */
    g_state = state;

    /* End HDFE setup timing, start FWL timing */
    t_hdfe_setup = ctools_timer_seconds() - t_hdfe_setup_start;
    double t_fwl_start = ctools_timer_seconds();

    /* FWL partialling: read partial indices for iterative absorption with FE */
    ST_double dval_partial;
    ST_int n_partial = 0;
    ST_int *partial_indices = NULL;
    ST_int *is_partial = NULL;

    if (SF_scal_use("__civreghdfe_n_partial", &dval_partial) == 0) {
        n_partial = (ST_int)dval_partial;
    }

    if (n_partial > 0 && K_exog > 0) {
        /* Read partial variable indices (1-based indices into X_exog) */
        partial_indices = (ST_int *)malloc(n_partial * sizeof(ST_int));
        is_partial = (ST_int *)calloc(K_exog, sizeof(ST_int));  /* Mask for partial vars */
        ST_double dval_idx;

        for (ST_int pi = 0; pi < n_partial; pi++) {
            char scal_name[64];
            snprintf(scal_name, sizeof(scal_name), "__civreghdfe_partial_%d", (int)(pi + 1));
            if (SF_scal_use(scal_name, &dval_idx) == 0) {
                partial_indices[pi] = (ST_int)dval_idx;
                if (partial_indices[pi] >= 1 && partial_indices[pi] <= K_exog) {
                    is_partial[partial_indices[pi] - 1] = 1;
                }
            }
        }
    }

    /* End FWL timing, start partial out timing */
    t_fwl = ctools_timer_seconds() - t_fwl_start;
    double t_partial_start = ctools_timer_seconds();

    /* Partial out FEs from all variables using CG solver */
    ST_int total_cols = 1 + K_endog + K_exog + K_iv;

    /* Combine all data into one array for parallel processing */
    ST_double *all_data = (ST_double *)malloc((size_t)N * total_cols * sizeof(ST_double));
    if (!all_data) {
        SF_error("civreghdfe: Memory allocation failed for demeaning buffer\n");
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        free(fe_levels_c);  /* levels freed by state cleanup */
        ctools_hdfe_state_cleanup(state);
        free(state);
        g_state = NULL;
        return 920;
    }

    /* Copy y */
    memcpy(all_data, y_c, N * sizeof(ST_double));

    /* Copy X_endog */
    if (K_endog > 0) {
        memcpy(all_data + N, X_endog_c, N * K_endog * sizeof(ST_double));
    }

    /* Copy X_exog */
    if (K_exog > 0) {
        memcpy(all_data + N * (1 + K_endog), X_exog_c, N * K_exog * sizeof(ST_double));
    }

    /* Copy Z */
    memcpy(all_data + N * (1 + K_endog + K_exog), Z_c, N * K_iv * sizeof(ST_double));

    /* Iterative absorption: alternate between FE demeaning and partial variable
       residualization until convergence. This matches reghdfe's approach.
       Key insight: use the CURRENT partial variable columns from all_data for
       each projection, not a fixed precomputed P. */
    if (n_partial > 0 && K_exog > 0 && partial_indices && is_partial) {
        /* Allocate workspace */
        ST_double *P_cur = (ST_double *)malloc((size_t)N * n_partial * sizeof(ST_double));
        ST_double *PtP = (ST_double *)calloc(n_partial * n_partial, sizeof(ST_double));
        ST_double *PtP_inv = (ST_double *)calloc(n_partial * n_partial, sizeof(ST_double));
        ST_double *Ptx = (ST_double *)calloc(n_partial, sizeof(ST_double));
        ST_double *coef = (ST_double *)calloc(n_partial, sizeof(ST_double));
        ST_double *old_data = (ST_double *)malloc((size_t)N * total_cols * sizeof(ST_double));

        /* Get column offsets for partial variables in all_data */
        ST_int *partial_col_offsets = (ST_int *)malloc(n_partial * sizeof(ST_int));
        for (ST_int pi = 0; pi < n_partial; pi++) {
            ST_int exog_idx = partial_indices[pi] - 1;  /* 0-based index in X_exog */
            partial_col_offsets[pi] = 1 + K_endog + exog_idx;  /* Column in all_data */
        }

        const int max_iter = 100;
        const double tol = 1e-10;

        for (int iter = 0; iter < max_iter; iter++) {
            /* Save current state for convergence check */
            memcpy(old_data, all_data, (size_t)N * total_cols * sizeof(ST_double));

            /* Step 1: Demean by FE */
            partial_out_columns(state, all_data, N, total_cols, num_threads);

            /* Step 2: Extract current P from all_data and compute (P'P)^{-1} */
            for (ST_int pi = 0; pi < n_partial; pi++) {
                ST_int col = partial_col_offsets[pi];
                memcpy(P_cur + pi * N, all_data + col * N, N * sizeof(ST_double));
            }

            /* Compute P'P */
            memset(PtP, 0, n_partial * n_partial * sizeof(ST_double));
            ctools_matmul_atb(P_cur, P_cur, N, n_partial, n_partial, PtP);

            /* Invert P'P */
            memcpy(PtP_inv, PtP, n_partial * n_partial * sizeof(ST_double));
            int ptp_ok = (cholesky(PtP_inv, n_partial) == 0 &&
                          invert_from_cholesky(PtP_inv, n_partial, PtP_inv) == 0);

            if (!ptp_ok) {
                break;
            }

            /* Step 3: Residualize non-partial columns by current P */
            for (ST_int col = 0; col < total_cols; col++) {
                /* Skip the partial variable columns themselves */
                int is_partial_col = 0;
                for (ST_int pi = 0; pi < n_partial; pi++) {
                    if (col == partial_col_offsets[pi]) {
                        is_partial_col = 1;
                        break;
                    }
                }
                if (is_partial_col) continue;

                ST_double *x = all_data + col * N;

                /* Compute P'x */
                for (ST_int p = 0; p < n_partial; p++) {
                    ST_double sum = 0.0;
                    for (ST_int i = 0; i < N; i++) {
                        sum += P_cur[p * N + i] * x[i];
                    }
                    Ptx[p] = sum;
                }

                /* Compute coef = (P'P)^{-1} P'x */
                for (ST_int p = 0; p < n_partial; p++) {
                    ST_double sum = 0.0;
                    for (ST_int q = 0; q < n_partial; q++) {
                        sum += PtP_inv[q * n_partial + p] * Ptx[q];
                    }
                    coef[p] = sum;
                }

                /* Residualize: x = x - P * coef */
                for (ST_int i = 0; i < N; i++) {
                    ST_double fitted = 0.0;
                    for (ST_int p = 0; p < n_partial; p++) {
                        fitted += P_cur[p * N + i] * coef[p];
                    }
                    x[i] -= fitted;
                }
            }

            /* Check convergence: max change in any element */
            ST_double max_change = 0.0;
            for (ST_int j = 0; j < N * total_cols; j++) {
                ST_double diff = fabs(all_data[j] - old_data[j]);
                if (diff > max_change) max_change = diff;
            }

            if (max_change < tol) {
                break;
            }
        }

        free(P_cur);
        free(PtP);
        free(PtP_inv);
        free(Ptx);
        free(coef);
        free(old_data);
        free(partial_col_offsets);
    } else {
        /* No partial variables - just demean by FE once */
        partial_out_columns(state, all_data, N, total_cols, num_threads);
    }

    /* End partial out timing, start post-processing timing */
    t_partial_out = ctools_timer_seconds() - t_partial_start;
    double t_postproc_start = ctools_timer_seconds();

    /* Extract demeaned data */
    ST_double *y_dem = all_data;
    ST_double *X_endog_dem = (K_endog > 0) ? all_data + N : NULL;
    ST_double *X_exog_dem = (K_exog > 0) ? all_data + N * (1 + K_endog) : NULL;
    ST_double *Z_dem = all_data + N * (1 + K_endog + K_exog);

    /* Remove partial variables from X_exog AND from Z (instruments) after convergence.
       Z contains [exog_vars, excluded_instruments], so partial exog vars must be removed
       from both X_exog and the first K_exog columns of Z. */
    if (n_partial > 0 && K_exog > 0 && partial_indices && is_partial) {
        ST_int K_exog_new = K_exog - n_partial;
        ST_int K_excl = K_iv - K_exog;  /* Number of excluded instruments */

        /* Remove partial columns from X_exog */
        if (K_exog_new > 0) {
            ST_int new_idx = 0;
            for (ST_int k = 0; k < K_exog; k++) {
                if (!is_partial[k]) {
                    if (new_idx != k) {
                        memcpy(X_exog_dem + new_idx * N, X_exog_dem + k * N, N * sizeof(ST_double));
                    }
                    new_idx++;
                }
            }
        }

        /* Remove partial columns from Z (which starts with K_exog columns of exog vars) */
        /* After removal, Z should have K_exog_new + K_excl columns */
        ST_int new_idx = 0;
        for (ST_int k = 0; k < K_exog; k++) {
            if (!is_partial[k]) {
                if (new_idx != k) {
                    memcpy(Z_dem + new_idx * N, Z_dem + k * N, N * sizeof(ST_double));
                }
                new_idx++;
            }
        }
        /* Shift excluded instruments down */
        for (ST_int k = 0; k < K_excl; k++) {
            memcpy(Z_dem + new_idx * N, Z_dem + (K_exog + k) * N, N * sizeof(ST_double));
            new_idx++;
        }

        K_exog = K_exog_new;
        K_iv = K_exog_new + K_excl;
        X_exog_dem = (K_exog > 0) ? all_data + N * (1 + K_endog) : NULL;

        free(partial_indices);
        free(is_partial);
        partial_indices = NULL;
        is_partial = NULL;
    } else if (partial_indices) {
        free(partial_indices);
        free(is_partial);
        partial_indices = NULL;
        is_partial = NULL;
    }

    /* Compute df_a (absorbed degrees of freedom) and mobility groups */
    ST_int df_a = 0;
    ST_int df_a_nested = 0;  /* Levels from FE nested within cluster */
    ST_int mobility_groups = 0;

    ctools_compute_hdfe_dof(state->factors, G, N, &df_a, &mobility_groups);

    for (ST_int g = 0; g < G; g++) {
        /* Check if this FE is nested in cluster (1-indexed) */
        if (nested_fe_index > 0 && g == (nested_fe_index - 1)) {
            df_a_nested = state->factors[g].num_levels;
        }
        /* Save per-FE num_levels to Stata scalars for absorbed DOF table */
        char scalar_name[64];
        snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_num_levels_%d", (int)(g + 1));
        ctools_scal_save(scalar_name, (ST_double)state->factors[g].num_levels);
    }

    /* Save mobility groups for absorbed DOF table display */
    ctools_scal_save("__civreghdfe_mobility_groups", (ST_double)mobility_groups);

    /* For VCE calculation when FE is nested in cluster:
       - nested FE levels don't contribute to df_a for VCE
       - nested_adj = 1 to account for the constant */
    ST_int df_a_for_vce = df_a - df_a_nested;
    ST_int nested_adj = (df_a_nested > 0) ? 1 : 0;

    /* Remap cluster IDs to contiguous 1-based indices using shared utility */
    ST_int num_clusters = 0;
    if (has_cluster) {
        if (ctools_remap_cluster_ids(cluster_ids_c, N, &num_clusters) != 0) {
            SF_error("civreghdfe: Failed to remap cluster IDs (memory allocation error)\n");
            free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
            free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
            free(fe_levels_c);  /* levels freed by state cleanup */
            ctools_hdfe_state_cleanup(state);
            free(state);
            g_state = NULL;
            return 920;
        }
    }

    /* Remap cluster2 IDs for two-way clustering */
    ST_int num_clusters2 = 0;
    if (has_cluster2) {
        if (ctools_remap_cluster_ids(cluster2_ids_c, N, &num_clusters2) != 0) {
            SF_error("civreghdfe: Failed to remap cluster2 IDs (memory allocation error)\n");
            free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
            free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
            free(fe_levels_c);  /* levels freed by state cleanup */
            ctools_hdfe_state_cleanup(state);
            free(state);
            g_state = NULL;
            return 920;
        }
    }

    /* Detect FEs nested within cluster variable */
    ST_int *fe_nested = (ST_int *)calloc(G, sizeof(ST_int));
    if (has_cluster && fe_nested) {
        for (ST_int g = 0; g < G; g++) {
            ST_int is_nested = ctools_fe_nested_in_cluster(
                state->factors[g].levels, state->factors[g].num_levels,
                cluster_ids_c, N);
            if (is_nested < 0) is_nested = 0;
            fe_nested[g] = is_nested;

            char scalar_name[64];
            snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_fe_nested_%d", (int)(g + 1));
            ctools_scal_save(scalar_name, (ST_double)is_nested);
        }

    }
    if (fe_nested) free(fe_nested);

    /* Update K_total after possible FWL reduction of K_exog */
    K_total = K_exog + K_endog;

    /* ================================================================
     * COLLINEARITY DETECTION for X (regressors) and Z (instruments)
     * Detect after HDFE partialling, before estimation.
     * Same approach as creghdfe: FE-absorbed variance check + Cholesky.
     * ================================================================ */
    ST_int K_exog_orig = K_exog;
    ST_int K_endog_orig = K_endog;
    ST_int K_total_orig = K_total;

    /* is_collinear_x: flags for [exog, endog] in internal order */
    ST_int *is_collinear_x = (ST_int *)calloc(K_total, sizeof(ST_int));
    ST_int num_collinear_x = 0;

    if (!is_collinear_x) {
        SF_error("civreghdfe: Memory allocation failed for collinearity arrays\n");
        free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        ctools_hdfe_state_cleanup(state);
        free(state);
        g_state = NULL;
        return 920;
    }

    /* Stage 1: FE-absorbed variance check — mark columns with near-zero
       variance after partialling as collinear (same as creghdfe) */
    {
        for (ST_int k = 0; k < K_exog; k++) {
            ST_double xx_partial = fast_dot(X_exog_dem + k * N, X_exog_dem + k * N, N);
            if (xx_partial < 1e-30) {
                is_collinear_x[k] = 1;
                num_collinear_x++;
            }
        }
        for (ST_int k = 0; k < K_endog; k++) {
            ST_double xx_partial = fast_dot(X_endog_dem + k * N, X_endog_dem + k * N, N);
            if (xx_partial < 1e-30) {
                is_collinear_x[K_exog + k] = 1;
                num_collinear_x++;
            }
        }
    }

    /* Stage 2: Numerical collinearity via Cholesky on X'X
       Build concatenated X = [X_exog, X_endog] in column-major order */
    if (K_total > 1) {
        ST_double *XtX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        if (!XtX) {
            free(is_collinear_x);
            SF_error("civreghdfe: Memory allocation failed for XtX\n");
            free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
            free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
            ctools_hdfe_state_cleanup(state);
            free(state);
            g_state = NULL;
            return 920;
        }

        /* Compute X'X where X = [X_exog_dem, X_endog_dem] using block matmul.
         * XtX (K_total x K_total) has blocks:
         *   [Xe'Xe  Xe'Xn]
         *   [Xn'Xe  Xn'Xn]
         * Use ctools_matmul_atb for each block (OpenMP + SIMD) instead of
         * K_total² individual fast_dot calls. */
        if (K_exog > 0 && K_endog > 0) {
            /* Both exog and endog present — compute 3 blocks (4th by symmetry) */
            ST_double *blk_ee = (ST_double *)malloc(K_exog * K_exog * sizeof(ST_double));
            ST_double *blk_en = (ST_double *)malloc(K_exog * K_endog * sizeof(ST_double));
            ST_double *blk_nn = (ST_double *)malloc(K_endog * K_endog * sizeof(ST_double));
            if (blk_ee && blk_en && blk_nn) {
                ctools_matmul_atb(X_exog_dem, X_exog_dem, N, K_exog, K_exog, blk_ee);
                ctools_matmul_atb(X_exog_dem, X_endog_dem, N, K_exog, K_endog, blk_en);
                ctools_matmul_atb(X_endog_dem, X_endog_dem, N, K_endog, K_endog, blk_nn);

                /* Place blocks into K_total x K_total XtX (column-major) */
                for (ST_int jj = 0; jj < K_exog; jj++)
                    for (ST_int ii = 0; ii < K_exog; ii++)
                        XtX[jj * K_total + ii] = blk_ee[jj * K_exog + ii];
                for (ST_int jj = 0; jj < K_endog; jj++)
                    for (ST_int ii = 0; ii < K_exog; ii++) {
                        XtX[(K_exog + jj) * K_total + ii] = blk_en[jj * K_exog + ii];
                        XtX[ii * K_total + (K_exog + jj)] = blk_en[jj * K_exog + ii];
                    }
                for (ST_int jj = 0; jj < K_endog; jj++)
                    for (ST_int ii = 0; ii < K_endog; ii++)
                        XtX[(K_exog + jj) * K_total + (K_exog + ii)] = blk_nn[jj * K_endog + ii];
            } else {
                /* Fallback: element-wise */
                for (ST_int ii = 0; ii < K_total; ii++) {
                    const ST_double *xi = (ii < K_exog) ? X_exog_dem + ii * N
                                                        : X_endog_dem + (ii - K_exog) * N;
                    for (ST_int jj = ii; jj < K_total; jj++) {
                        const ST_double *xj = (jj < K_exog) ? X_exog_dem + jj * N
                                                            : X_endog_dem + (jj - K_exog) * N;
                        ST_double val = fast_dot(xi, xj, N);
                        XtX[jj * K_total + ii] = val;
                        XtX[ii * K_total + jj] = val;
                    }
                }
            }
            free(blk_ee); free(blk_en); free(blk_nn);
        } else if (K_exog > 0) {
            ctools_matmul_atb(X_exog_dem, X_exog_dem, N, K_exog, K_exog, XtX);
        } else {
            ctools_matmul_atb(X_endog_dem, X_endog_dem, N, K_endog, K_endog, XtX);
        }

        ST_int num_chol_collinear = detect_collinearity(XtX, K_total, is_collinear_x, verbose);
        free(XtX);

        if (num_chol_collinear < 0) {
            free(is_collinear_x);
            SF_error("civreghdfe: Collinearity detection failed\n");
            free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
            free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
            ctools_hdfe_state_cleanup(state);
            free(state);
            g_state = NULL;
            return 920;
        }

        /* Recount total collinear */
        num_collinear_x = 0;
        for (ST_int k = 0; k < K_total; k++) {
            if (is_collinear_x[k]) num_collinear_x++;
        }
    }

    /* Store collinearity flags to Stata scalars (in [exog, endog] order) */
    {
        char scalar_name[64];
        ctools_scal_save("__civreghdfe_num_collinear", (ST_double)num_collinear_x);
        for (ST_int k = 0; k < K_total_orig; k++) {
            snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_collinear_%d", (int)(k + 1));
            ctools_scal_save(scalar_name, (ST_double)is_collinear_x[k]);
        }
    }

    /* Compact X arrays in-place if collinear columns found */
    if (num_collinear_x > 0) {
        /* Compact X_exog_dem */
        ST_int new_K_exog = 0;
        for (ST_int k = 0; k < K_exog; k++) {
            if (!is_collinear_x[k]) {
                if (new_K_exog != k) {
                    memcpy(X_exog_dem + new_K_exog * N, X_exog_dem + k * N, N * sizeof(ST_double));
                }
                new_K_exog++;
            }
        }

        /* Compact X_endog_dem */
        ST_int new_K_endog = 0;
        for (ST_int k = 0; k < K_endog; k++) {
            if (!is_collinear_x[K_exog + k]) {
                if (new_K_endog != k) {
                    memcpy(X_endog_dem + new_K_endog * N, X_endog_dem + k * N, N * sizeof(ST_double));
                }
                new_K_endog++;
            }
        }

        /* Compact exogenous portion of Z_dem (first K_exog columns match X_exog) */
        ST_int K_excl = K_iv - K_exog;
        ST_int new_z_idx = 0;
        for (ST_int k = 0; k < K_exog; k++) {
            if (!is_collinear_x[k]) {
                if (new_z_idx != k) {
                    memcpy(Z_dem + new_z_idx * N, Z_dem + k * N, N * sizeof(ST_double));
                }
                new_z_idx++;
            }
        }
        /* Shift excluded instruments down after compacted exog portion */
        for (ST_int k = 0; k < K_excl; k++) {
            if (new_z_idx != K_exog + k) {
                memcpy(Z_dem + new_z_idx * N, Z_dem + (K_exog + k) * N, N * sizeof(ST_double));
            }
            new_z_idx++;
        }

        K_exog = new_K_exog;
        K_endog = new_K_endog;
        K_total = K_exog + K_endog;
        K_iv = new_z_idx;  /* new_K_exog + K_excl */

        if (verbose) {
            char msg[128];
            snprintf(msg, sizeof(msg), "civreghdfe: Dropped %d collinear regressor(s)\n",
                     (int)num_collinear_x);
            SF_display(msg);
        }
    }

    /* Z collinearity detection */
    ST_int *is_collinear_z = (ST_int *)calloc(K_iv, sizeof(ST_int));
    ST_int num_collinear_z = 0;

    if (is_collinear_z && K_iv > 1) {
        ST_double *ZtZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
        if (ZtZ) {
            ctools_matmul_atb(Z_dem, Z_dem, N, K_iv, K_iv, ZtZ);
            ST_int num_z_collinear = detect_collinearity(ZtZ, K_iv, is_collinear_z, verbose);
            free(ZtZ);

            if (num_z_collinear > 0) {
                /* Compact Z_dem in-place */
                ST_int new_K_iv = 0;
                for (ST_int k = 0; k < K_iv; k++) {
                    if (!is_collinear_z[k]) {
                        if (new_K_iv != k) {
                            memcpy(Z_dem + new_K_iv * N, Z_dem + k * N, N * sizeof(ST_double));
                        }
                        new_K_iv++;
                    }
                }
                num_collinear_z = K_iv - new_K_iv;
                K_iv = new_K_iv;

                if (verbose) {
                    char msg[128];
                    snprintf(msg, sizeof(msg), "civreghdfe: Dropped %d collinear instrument(s)\n",
                             (int)num_collinear_z);
                    SF_display(msg);
                }
            }
        }
    }
    free(is_collinear_z);

    /* Post-compaction identification check */
    if (K_iv < K_total) {
        free(is_collinear_x);
        SF_error("civreghdfe: Model is underidentified after removing collinear variables\n");
        free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        ctools_hdfe_state_cleanup(state);
        free(state);
        g_state = NULL;
        return 481;
    }

    /* Store compacted dimensions */
    ctools_scal_save("__civreghdfe_K_keep", (ST_double)K_total);

    /* Allocate output arrays */
    ST_double *beta = (ST_double *)calloc(K_total, sizeof(ST_double));
    ST_double *V = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *first_stage_F = (ST_double *)calloc(K_endog, sizeof(ST_double));

    /* Check allocations - critical for preventing NULL pointer dereference */
    if (!beta || !V || !first_stage_F) {
        SF_error("civreghdfe: Memory allocation failed for output arrays\n");
        free(beta); free(V); free(first_stage_F);
        free(is_collinear_x);
        free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        ctools_hdfe_state_cleanup(state);
        free(state);
        g_state = NULL;
        return 920;  /* Memory allocation error */
    }

    /* End post-processing timing, start estimation timing */
    t_postproc = ctools_timer_seconds() - t_postproc_start;
    double t_estimate_start = ctools_timer_seconds();

    /* Compute k-class IV estimation (2SLS, LIML, Fuller, etc.) */
    /* For cluster VCE, pass df_a_for_vce (excluding nested FE levels) and nested_adj */
    /* For panel-aware HAC, pass first FE variable as panel IDs */
    ST_int *hac_panel_ids = NULL;
    ST_int num_hac_panels = 0;
    if (hac_panel && G > 0) {
        hac_panel_ids = fe_levels_c[0];  /* First absorb variable is panel ID */
        num_hac_panels = state->factors[0].num_levels;
    }

    ST_double lambda = 1.0;
    ST_retcode rc = ivest_compute_2sls(
        y_dem, X_exog_dem, X_endog_dem, Z_dem,
        weights_c, weight_type,
        N, K_exog, K_endog, K_iv,
        beta, V, first_stage_F,
        vce_type, cluster_ids_c, num_clusters,
        cluster2_ids_c, num_clusters2,
        df_a_for_vce, nested_adj, verbose,
        est_method, kclass_user, fuller_alpha, &lambda,
        kernel_type, bw, kiefer,
        hac_panel_ids, num_hac_panels,
        sdofminus_opt, center,
        kiefer ? Z_c : NULL  /* Pass original Z for Kiefer VCE */
    );

    if (rc != STATA_OK) {
        /* Cleanup and return */
        free(is_collinear_x);
        free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        free(beta); free(V); free(first_stage_F);
        ctools_hdfe_state_cleanup(state);
        free(state);
        g_state = NULL;
        return rc;
    }

    /* End estimation timing (includes VCE), start stats computation */
    t_estimate = ctools_timer_seconds() - t_estimate_start;
    double t_stats_start = ctools_timer_seconds();

    /* Compute RSS, TSS, R-squared, Root MSE, and F-statistic */
    /* These are computed on demeaned (partialled-out) data */
    ST_double rss = 0.0;
    ST_double tss = 0.0;
    ST_double y_mean = 0.0;
    ST_double sum_w = 0.0;

    /* Compute weighted mean of y_dem */
    for (ST_int i = 0; i < N; i++) {
        ST_double w = (weights_c && weight_type != 0) ? weights_c[i] : 1.0;
        y_mean += w * y_dem[i];
        sum_w += w;
    }
    y_mean /= sum_w;

    /* Compute TSS (centered) = sum(w * (y - ybar)^2) */
    for (ST_int i = 0; i < N; i++) {
        ST_double w = (weights_c && weight_type != 0) ? weights_c[i] : 1.0;
        ST_double dev = y_dem[i] - y_mean;
        tss += w * dev * dev;
    }

    /* Compute RSS = sum(w * resid^2) */
    /* Compute fitted values: yhat = X_exog * beta[0:K_exog-1] + X_endog * beta[K_exog:K_total-1] */
    for (ST_int i = 0; i < N; i++) {
        ST_double fitted = 0.0;
        for (ST_int j = 0; j < K_exog; j++) {
            fitted += X_exog_dem[j * N + i] * beta[j];
        }
        for (ST_int j = 0; j < K_endog; j++) {
            fitted += X_endog_dem[j * N + i] * beta[K_exog + j];
        }
        ST_double resid = y_dem[i] - fitted;
        ST_double w = (weights_c && weight_type != 0) ? weights_c[i] : 1.0;
        rss += w * resid * resid;
    }

    ST_double r2 = (tss > 0) ? 1.0 - rss / tss : 0.0;

    /* Compute df_r:
       - For Driscoll-Kraay (vce_type == 4): df_r = T - 1 (time periods minus 1)
       - For clustered VCE (vce_type == 2, 3): df_r = num_clusters - 1 (Stata convention)
       - Otherwise: df_r = N - K - df_a - dofminus
       - nopartialsmall: exclude partialled variables from K for DOF calc */
    ST_int K_for_dof = nopartialsmall ? (K_total - n_partial) : K_total;
    ST_int df_r_val;
    if (vce_type == 4 && dkraay_T > 0) {
        /* Driscoll-Kraay: df_r = T - 1 */
        df_r_val = dkraay_T - 1;
    } else if (has_cluster) {
        df_r_val = num_clusters - 1;
    } else {
        df_r_val = N - K_for_dof - df_a - dofminus;
    }
    if (df_r_val <= 0) df_r_val = 1;
    /* For rmse, use sdofminus = max(1, df_a_for_vce) when nested, df_a otherwise
       (matches ivreghdfe line 670: if HDFE.df_a=0, force absorb_ct to 1) */
    ST_int sdofminus = (has_cluster && df_a_nested > 0) ?
                       (df_a_for_vce > 0 ? df_a_for_vce : 1) : df_a;
    ST_double rmse = sqrt(rss / (N - K_total - sdofminus > 0 ? N - K_total - sdofminus : 1));

    /* Compute model F-statistic: (R2 / K) / ((1 - R2) / df_r) */
    ST_double f_stat = 0.0;
    if (K_total > 0 && r2 < 1.0) {
        f_stat = (r2 / K_total) / ((1.0 - r2) / df_r_val);
    }

    /* End stats timing, start store timing */
    double t_stats = ctools_timer_seconds() - t_stats_start;
    double t_store_start = ctools_timer_seconds();

    /* Store results to Stata using checked wrappers */
    /* Scalars - use error-checking wrappers to prevent silent failures */
    ctools_scal_save("__civreghdfe_N", (ST_double)N);
    ctools_scal_save("__civreghdfe_df_r", (ST_double)df_r_val);
    ctools_scal_save("__civreghdfe_df_a", (ST_double)df_a);
    ctools_scal_save("__civreghdfe_df_a_for_vce", (ST_double)df_a_for_vce);
    ctools_scal_save("__civreghdfe_K", (ST_double)K_total);
    ctools_scal_save("__civreghdfe_rss", rss);
    ctools_scal_save("__civreghdfe_tss", tss);
    ctools_scal_save("__civreghdfe_r2", r2);
    ctools_scal_save("__civreghdfe_rmse", rmse);
    ctools_scal_save("__civreghdfe_F", f_stat);

    if (vce_type == 4 && dkraay_T > 0) {
        /* Driscoll-Kraay: N_clust = number of time periods */
        ctools_scal_save("__civreghdfe_N_clust", (ST_double)dkraay_T);
    } else if (has_cluster) {
        ctools_scal_save("__civreghdfe_N_clust", (ST_double)num_clusters);
    }

    if (has_cluster2) {
        ctools_scal_save("__civreghdfe_N_clust2", (ST_double)num_clusters2);
    }

    /* Store lambda for LIML/Fuller */
    if (est_method == 1 || est_method == 2) {
        ctools_scal_save("__civreghdfe_lambda", lambda);
    }

    /* Matrices: e(b) and e(V) - expand from compact to full dimensions,
       then reorder from [exog, endog] to [endog, exog] for Stata.
       Compact: beta[0..K_exog-1] = exog, beta[K_exog..K_total-1] = endog
       Full: K_exog_orig exog + K_endog_orig endog (zeros for collinear)
       Stata: [endog_orig, exog_orig] to match ivreghdfe convention */

    /* Build mapping from full index to compact index (-1 if collinear) */
    /* is_collinear_x[0..K_exog_orig-1] = exog flags,
       is_collinear_x[K_exog_orig..K_total_orig-1] = endog flags */
    ST_int *full_to_compact = (ST_int *)malloc(K_total_orig * sizeof(ST_int));
    if (full_to_compact) {
        ST_int compact_idx = 0;
        for (ST_int k = 0; k < K_total_orig; k++) {
            if (is_collinear_x[k]) {
                full_to_compact[k] = -1;
            } else {
                full_to_compact[k] = compact_idx++;
            }
        }
    }

    /* Store beta in [endog_orig, exog_orig] order with zeros for collinear */
    for (ST_int e = 0; e < K_endog_orig; e++) {
        ST_int full_idx = K_exog_orig + e;  /* Position in full [exog, endog] */
        ST_double val = 0.0;
        if (full_to_compact && full_to_compact[full_idx] >= 0) {
            val = beta[full_to_compact[full_idx]];
        }
        ctools_mat_store("__civreghdfe_b", 1, e + 1, val);
    }
    for (ST_int x = 0; x < K_exog_orig; x++) {
        ST_double val = 0.0;
        if (full_to_compact && full_to_compact[x] >= 0) {
            val = beta[full_to_compact[x]];
        }
        ctools_mat_store("__civreghdfe_b", 1, K_endog_orig + x + 1, val);
    }

    /* Store V in [endog_orig, exog_orig] order */
    /* new_order: [endog_0..endog_{Ke_orig-1}, exog_0..exog_{Kx_orig-1}]
       For new index i in [0, K_endog_orig): full index = K_exog_orig + i
       For new index i in [K_endog_orig, K_total_orig): full index = i - K_endog_orig */
    for (ST_int new_i = 0; new_i < K_total_orig; new_i++) {
        ST_int full_i = (new_i < K_endog_orig) ? (K_exog_orig + new_i) : (new_i - K_endog_orig);
        ST_int compact_i = full_to_compact ? full_to_compact[full_i] : full_i;
        for (ST_int new_j = 0; new_j < K_total_orig; new_j++) {
            ST_int full_j = (new_j < K_endog_orig) ? (K_exog_orig + new_j) : (new_j - K_endog_orig);
            ST_int compact_j = full_to_compact ? full_to_compact[full_j] : full_j;
            ST_double val = 0.0;
            if (compact_i >= 0 && compact_j >= 0) {
                val = V[compact_j * K_total + compact_i];
            }
            ctools_mat_store("__civreghdfe_V", new_i + 1, new_j + 1, val);
        }
    }
    free(full_to_compact);

    /* First stage F-stats */
    for (ST_int e = 0; e < K_endog; e++) {
        char name[64];
        snprintf(name, sizeof(name), "__civreghdfe_F1_%d", (int)(e + 1));
        ctools_scal_save(name, first_stage_F[e]);
    }

    /* End store timing, compute total */
    t_store = ctools_timer_seconds() - t_store_start;
    double t_total = ctools_timer_seconds() - t_total_start;

    /* Save timing scalars - less critical but still use checked wrappers */
    ctools_scal_save("_civreghdfe_time_load", t_load);
    ctools_scal_save("_civreghdfe_time_extract", t_extract);
    ctools_scal_save("_civreghdfe_time_missing", t_missing);
    ctools_scal_save("_civreghdfe_time_singleton", t_singleton);
    ctools_scal_save("_civreghdfe_time_remap", t_hdfe_setup);
    ctools_scal_save("_civreghdfe_time_fwl", t_fwl);
    ctools_scal_save("_civreghdfe_time_partial", t_partial_out);
    ctools_scal_save("_civreghdfe_time_dof", t_postproc);
    ctools_scal_save("_civreghdfe_time_estimate", t_estimate);
    ctools_scal_save("_civreghdfe_time_stats", t_stats);
    ctools_scal_save("_civreghdfe_time_store", t_store);
    ctools_scal_save("_civreghdfe_time_total", t_total);
    CTOOLS_SAVE_THREAD_INFO("_civreghdfe");

    /* Cleanup */
    free(is_collinear_x);
    free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
    free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
    free(beta); free(V); free(first_stage_F);

    /* Cleanup state */
    ctools_hdfe_state_cleanup(state);
    free(fe_levels_c);
    free(state);
    g_state = NULL;

    return STATA_OK;
}

/*
    Main entry point for civreghdfe plugin.
*/
ST_retcode civreghdfe_main(const char *args)
{
    if (args == NULL || strlen(args) == 0) {
        SF_error("civreghdfe: No subcommand specified\n");
        return 198;
    }

    if (strcmp(args, "iv_regression") == 0) {
        return do_iv_regression();
    }

    SF_error("civreghdfe: Unknown subcommand\n");
    return 198;
}

/*
 * Global flag to skip saving results to Stata (for debugging).
 * When set to 1, all SF_scal_save and SF_mat_store calls are skipped.
 */
int g_civreghdfe_noreturn = 0;

/*
 * Cleanup function for civreghdfe persistent state.
 * Frees the global HDFE state if allocated.
 * Safe to call multiple times (idempotent).
 */
void civreghdfe_cleanup_state(void)
{
    /* civreghdfe cleans up g_state at the end of do_iv_regression(),
       so this is just a safety measure for interrupted execution */
    if (g_state != NULL) {
        /* Note: Full cleanup would require knowing G and num_threads,
           which we don't have here. The main function handles full cleanup.
           This just nulls the pointer to prevent double-free issues. */
        g_state = NULL;
    }
}

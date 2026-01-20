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
#include "../ctools_timer.h"
#include "../ctools_ols.h"
#include "../ctools_types.h"  /* For ctools_data_load_selective */
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
    double t_load = 0, t_singleton = 0, t_hdfe_setup = 0, t_fwl = 0, t_partial_out = 0;
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
    ST_int kernel_type = 0, bw = 0, dkraay = 0;
    SF_scal_use("__civreghdfe_kernel", &dval); kernel_type = (ST_int)dval;
    SF_scal_use("__civreghdfe_bw", &dval); bw = (ST_int)dval;
    SF_scal_use("__civreghdfe_dkraay", &dval); dkraay = (ST_int)dval;
    (void)dkraay;  /* Reserved for Driscoll-Kraay panel HAC */

    if (verbose) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "civreghdfe: N=%d, K_endog=%d, K_exog=%d, K_iv=%d, G=%d\n",
                 (int)N_total, (int)K_endog, (int)K_exog, (int)K_iv, (int)G);
        SF_display(buf);
    }

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
     * PARALLEL DATA LOADING using ctools_data_load_selective
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

    /* Load all variables in parallel using thread pool */
    stata_data loaded_data;
    stata_data_init(&loaded_data);

    stata_retcode load_rc = ctools_data_load_selective(&loaded_data, var_indices, total_vars, 0, 0);
    free(var_indices);

    if (load_rc != STATA_OK) {
        stata_data_free(&loaded_data);
        SF_error("civreghdfe: Parallel data load failed\n");
        return 920;
    }

    /* Build selection index for if/in filtering */
    ST_int N_range = in2 - in1 + 1;
    ST_int *sel_idx = (ST_int *)malloc(N_range * sizeof(ST_int));
    if (!sel_idx) {
        stata_data_free(&loaded_data);
        SF_error("civreghdfe: Memory allocation failed\n");
        return 920;
    }

    ST_int N_ifobs = 0;
    for (ST_int obs = in1; obs <= in2; obs++) {
        if (SF_ifobs(obs)) {
            sel_idx[N_ifobs] = obs - in1;  /* 0-based index into loaded_data */
            N_ifobs++;
        }
    }

    if (N_ifobs == 0) {
        SF_error("civreghdfe: No observations selected\n");
        free(sel_idx);
        stata_data_free(&loaded_data);
        return 2001;
    }

    /* Check if selection is identity (no filtering needed) */
    ST_int sel_is_identity = (N_ifobs == N_range) ? 1 : 0;

    N_total = N_ifobs;

    /* Allocate output arrays */
    ST_double *y = (ST_double *)malloc((size_t)N_total * sizeof(ST_double));
    ST_double *X_endog = (K_endog > 0) ? (ST_double *)malloc((size_t)N_total * K_endog * sizeof(ST_double)) : NULL;
    ST_double *X_exog = (K_exog > 0) ? (ST_double *)malloc((size_t)N_total * K_exog * sizeof(ST_double)) : NULL;
    ST_double *Z = (ST_double *)malloc((size_t)N_total * K_iv * sizeof(ST_double));
    ST_double *weights = has_weights ? (ST_double *)malloc((size_t)N_total * sizeof(ST_double)) : NULL;
    ST_int *cluster_ids = has_cluster ? (ST_int *)malloc((size_t)N_total * sizeof(ST_int)) : NULL;
    ST_int *cluster2_ids = has_cluster2 ? (ST_int *)malloc((size_t)N_total * sizeof(ST_int)) : NULL;

    /* Allocate FE level arrays */
    ST_int **fe_levels = (ST_int **)malloc((size_t)G * sizeof(ST_int *));
    int fe_alloc_failed = 0;
    if (fe_levels) {
        for (ST_int g = 0; g < G; g++) {
            fe_levels[g] = (ST_int *)malloc((size_t)N_total * sizeof(ST_int));
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
        free(sel_idx);
        stata_data_free(&loaded_data);
        return 920;
    }

    /* Extract data from loaded_data into output arrays with if/in filtering */
    ST_int i;
    if (sel_is_identity) {
        /* Fast path: direct copy */
        memcpy(y, loaded_data.vars[var_y_idx].data.dbl, N_total * sizeof(ST_double));

        for (ST_int k = 0; k < K_endog; k++) {
            memcpy(X_endog + k * N_total, loaded_data.vars[var_endog_start_idx + k].data.dbl,
                   N_total * sizeof(ST_double));
        }
        for (ST_int k = 0; k < K_exog; k++) {
            memcpy(X_exog + k * N_total, loaded_data.vars[var_exog_start_idx + k].data.dbl,
                   N_total * sizeof(ST_double));
        }
        for (ST_int k = 0; k < K_iv; k++) {
            memcpy(Z + k * N_total, loaded_data.vars[var_iv_start_idx + k].data.dbl,
                   N_total * sizeof(ST_double));
        }
        for (ST_int g = 0; g < G; g++) {
            double *src = loaded_data.vars[var_fe_start_idx + g].data.dbl;
            for (i = 0; i < N_total; i++) {
                fe_levels[g][i] = (ST_int)src[i];
            }
        }
        if (has_cluster) {
            double *src = loaded_data.vars[var_cluster_idx].data.dbl;
            for (i = 0; i < N_total; i++) {
                cluster_ids[i] = (ST_int)src[i];
            }
        }
        if (has_cluster2) {
            double *src = loaded_data.vars[var_cluster2_idx].data.dbl;
            for (i = 0; i < N_total; i++) {
                cluster2_ids[i] = (ST_int)src[i];
            }
        }
        if (has_weights) {
            memcpy(weights, loaded_data.vars[var_weight_idx].data.dbl, N_total * sizeof(ST_double));
        }
    } else {
        /* Slow path: selective copy using sel_idx */
        double *src_y = loaded_data.vars[var_y_idx].data.dbl;
        for (i = 0; i < N_total; i++) {
            y[i] = src_y[sel_idx[i]];
        }

        for (ST_int k = 0; k < K_endog; k++) {
            double *src = loaded_data.vars[var_endog_start_idx + k].data.dbl;
            for (i = 0; i < N_total; i++) {
                X_endog[k * N_total + i] = src[sel_idx[i]];
            }
        }
        for (ST_int k = 0; k < K_exog; k++) {
            double *src = loaded_data.vars[var_exog_start_idx + k].data.dbl;
            for (i = 0; i < N_total; i++) {
                X_exog[k * N_total + i] = src[sel_idx[i]];
            }
        }
        for (ST_int k = 0; k < K_iv; k++) {
            double *src = loaded_data.vars[var_iv_start_idx + k].data.dbl;
            for (i = 0; i < N_total; i++) {
                Z[k * N_total + i] = src[sel_idx[i]];
            }
        }
        for (ST_int g = 0; g < G; g++) {
            double *src = loaded_data.vars[var_fe_start_idx + g].data.dbl;
            for (i = 0; i < N_total; i++) {
                fe_levels[g][i] = (ST_int)src[sel_idx[i]];
            }
        }
        if (has_cluster) {
            double *src = loaded_data.vars[var_cluster_idx].data.dbl;
            for (i = 0; i < N_total; i++) {
                cluster_ids[i] = (ST_int)src[sel_idx[i]];
            }
        }
        if (has_cluster2) {
            double *src = loaded_data.vars[var_cluster2_idx].data.dbl;
            for (i = 0; i < N_total; i++) {
                cluster2_ids[i] = (ST_int)src[sel_idx[i]];
            }
        }
        if (has_weights) {
            double *src = loaded_data.vars[var_weight_idx].data.dbl;
            for (i = 0; i < N_total; i++) {
                weights[i] = src[sel_idx[i]];
            }
        }
    }

    /* Free loaded_data and sel_idx - we've extracted what we need */
    stata_data_free(&loaded_data);
    free(sel_idx);

    /* End data load, start missing value check timing */
    t_load = ctools_timer_seconds() - t_load_start;
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

        /* Check FE - already loaded and converted to int, check for very large negative values
           (which could indicate missing - Stata missing values become very large when cast to int) */
        /* Note: FE values were already validated as non-missing during the data load */

        /* Check weights */
        if (has_weights && is_valid) {
            if (SF_is_missing(weights[i]) || weights[i] <= 0) is_valid = 0;
        }

        /* Check cluster - already loaded as int, was already validated during data load */

        valid_mask[i] = is_valid;
        if (is_valid) N_valid++;
    }

    if (verbose) {
        char buf[256];
        snprintf(buf, sizeof(buf), "civreghdfe: %d valid observations (dropped %d)\n",
                 (int)N_valid, (int)(N_total - N_valid));
        SF_display(buf);
    }

    if (N_valid < K_total + K_iv + 1) {
        SF_error("civreghdfe: Insufficient observations\n");
        free(y); free(X_endog); free(X_exog); free(Z);
        free(weights); free(cluster_ids); free(cluster2_ids); free(valid_mask);
        for (ST_int g = 0; g < G; g++) free(fe_levels[g]);
        free(fe_levels);
        return 2001;
    }

    /* Compact data to remove invalid observations - cast to size_t to prevent 32-bit overflow */
    ST_double *y_c = (ST_double *)malloc((size_t)N_valid * sizeof(ST_double));
    ST_double *X_endog_c = (K_endog > 0) ? (ST_double *)malloc((size_t)N_valid * K_endog * sizeof(ST_double)) : NULL;
    ST_double *X_exog_c = (K_exog > 0) ? (ST_double *)malloc((size_t)N_valid * K_exog * sizeof(ST_double)) : NULL;
    ST_double *Z_c = (ST_double *)malloc((size_t)N_valid * K_iv * sizeof(ST_double));
    ST_double *weights_c = has_weights ? (ST_double *)malloc((size_t)N_valid * sizeof(ST_double)) : NULL;
    ST_int *cluster_ids_c = has_cluster ? (ST_int *)malloc((size_t)N_valid * sizeof(ST_int)) : NULL;
    ST_int *cluster2_ids_c = has_cluster2 ? (ST_int *)malloc((size_t)N_valid * sizeof(ST_int)) : NULL;
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

    /* Check all allocations */
    if (!y_c || !Z_c || (K_endog > 0 && !X_endog_c) || (K_exog > 0 && !X_exog_c) ||
        (has_weights && !weights_c) || (has_cluster && !cluster_ids_c) ||
        (has_cluster2 && !cluster2_ids_c) || fe_c_alloc_failed) {
        SF_error("civreghdfe: Memory allocation failed for compacted arrays\n");
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        if (fe_levels_c) {
            for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
            free(fe_levels_c);
        }
        free(y); free(X_endog); free(X_exog); free(Z);
        free(weights); free(cluster_ids); free(cluster2_ids); free(valid_mask);
        for (ST_int g = 0; g < G; g++) free(fe_levels[g]);
        free(fe_levels);
        return 920;
    }

    ST_int idx = 0;
    for (ST_int i = 0; i < N_total; i++) {
        if (!valid_mask[i]) continue;

        y_c[idx] = y[i];

        for (ST_int k = 0; k < K_endog; k++) {
            X_endog_c[k * N_valid + idx] = X_endog[k * N_total + i];
        }
        for (ST_int k = 0; k < K_exog; k++) {
            X_exog_c[k * N_valid + idx] = X_exog[k * N_total + i];
        }
        for (ST_int k = 0; k < K_iv; k++) {
            Z_c[k * N_valid + idx] = Z[k * N_total + i];
        }
        for (ST_int g = 0; g < G; g++) {
            fe_levels_c[g][idx] = fe_levels[g][i];
        }
        if (has_weights) weights_c[idx] = weights[i];
        if (has_cluster) cluster_ids_c[idx] = cluster_ids[i];
        if (has_cluster2) cluster2_ids_c[idx] = cluster2_ids[i];

        idx++;
    }

    /* Free original arrays */
    free(y); free(X_endog); free(X_exog); free(Z);
    free(weights); free(cluster_ids); free(cluster2_ids); free(valid_mask);
    for (ST_int g = 0; g < G; g++) free(fe_levels[g]);
    free(fe_levels);

    ST_int N = N_valid;

    /* End missing value check/compact timing, start singleton timing */
    double t_missing = ctools_timer_seconds() - t_missing_start;
    double t_singleton_start = ctools_timer_seconds();

    /* Singleton detection using shared utility from ctools_hdfe_utils */
    ST_int *singleton_mask = (ST_int *)malloc((size_t)N * sizeof(ST_int));
    if (!singleton_mask) {
        SF_error("civreghdfe: Memory allocation failed for singleton mask\n");
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
        return 920;
    }
    ST_int num_singletons_total = ctools_remove_singletons(
        fe_levels_c, G, N, singleton_mask, 100, verbose
    );

    /* Compact data if singletons were found */
    if (num_singletons_total > 0) {
        ST_int N_after = N - num_singletons_total;

        /* Allocate compacted arrays - cast to size_t to prevent 32-bit overflow */
        ST_double *y_new = (ST_double *)malloc((size_t)N_after * sizeof(ST_double));
        ST_double *X_endog_new = (K_endog > 0) ? (ST_double *)malloc((size_t)N_after * K_endog * sizeof(ST_double)) : NULL;
        ST_double *X_exog_new = (K_exog > 0) ? (ST_double *)malloc((size_t)N_after * K_exog * sizeof(ST_double)) : NULL;
        ST_double *Z_new = (ST_double *)malloc((size_t)N_after * K_iv * sizeof(ST_double));
        ST_double *weights_new = has_weights ? (ST_double *)malloc((size_t)N_after * sizeof(ST_double)) : NULL;
        ST_int *cluster_ids_new = has_cluster ? (ST_int *)malloc((size_t)N_after * sizeof(ST_int)) : NULL;
        ST_int *cluster2_ids_new = has_cluster2 ? (ST_int *)malloc((size_t)N_after * sizeof(ST_int)) : NULL;
        ST_int **fe_levels_new = (ST_int **)malloc((size_t)G * sizeof(ST_int *));
        for (ST_int g = 0; g < G; g++) {
            fe_levels_new[g] = (ST_int *)malloc((size_t)N_after * sizeof(ST_int));
        }

        /* Compact using shared utilities */
        ctools_compact_array_double(y_c, y_new, singleton_mask, N, N_after);
        if (K_endog > 0) ctools_compact_matrix_double(X_endog_c, X_endog_new, singleton_mask, N, N_after, K_endog);
        if (K_exog > 0) ctools_compact_matrix_double(X_exog_c, X_exog_new, singleton_mask, N, N_after, K_exog);
        ctools_compact_matrix_double(Z_c, Z_new, singleton_mask, N, N_after, K_iv);
        if (has_weights) ctools_compact_array_double(weights_c, weights_new, singleton_mask, N, N_after);
        if (has_cluster) ctools_compact_array_int(cluster_ids_c, cluster_ids_new, singleton_mask, N, N_after);
        if (has_cluster2) ctools_compact_array_int(cluster2_ids_c, cluster2_ids_new, singleton_mask, N, N_after);
        for (ST_int g = 0; g < G; g++) {
            ctools_compact_array_int(fe_levels_c[g], fe_levels_new[g], singleton_mask, N, N_after);
        }

        /* Free old arrays and swap in new ones */
        free(y_c); y_c = y_new;
        if (X_endog_c) free(X_endog_c); X_endog_c = X_endog_new;
        if (X_exog_c) free(X_exog_c); X_exog_c = X_exog_new;
        free(Z_c); Z_c = Z_new;
        if (weights_c) free(weights_c); weights_c = weights_new;
        if (cluster_ids_c) free(cluster_ids_c); cluster_ids_c = cluster_ids_new;
        if (cluster2_ids_c) free(cluster2_ids_c); cluster2_ids_c = cluster2_ids_new;
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
        fe_levels_c = fe_levels_new;

        N = N_after;
    } else if (verbose) {
        SF_display("civreghdfe: No singletons found\n");
    }
    free(singleton_mask);

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
            /* Clean up previously allocated factors */
            for (ST_int fg = 0; fg < g; fg++) {
                free(state->factors[fg].counts);
                if (state->factors[fg].weighted_counts) free(state->factors[fg].weighted_counts);
                free(state->factors[fg].means);
            }
            free(state->factors);
            free(state);
            free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
            free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
            for (ST_int fg = 0; fg < G; fg++) free(fe_levels_c[fg]);
            free(fe_levels_c);
            /* Note: singleton_mask already freed after singleton removal */
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

        /* Compute inverse counts for fast division in projection */
        state->factors[g].inv_counts = (ST_double *)malloc(num_levels * sizeof(ST_double));
        if (state->factors[g].inv_counts) {
            for (ST_int lev = 0; lev < num_levels; lev++) {
                state->factors[g].inv_counts[lev] =
                    (state->factors[g].counts[lev] > 0) ? 1.0 / state->factors[g].counts[lev] : 0.0;
            }
        }
        if (has_weights && state->factors[g].weighted_counts) {
            state->factors[g].inv_weighted_counts = (ST_double *)malloc(num_levels * sizeof(ST_double));
            if (state->factors[g].inv_weighted_counts) {
                for (ST_int lev = 0; lev < num_levels; lev++) {
                    state->factors[g].inv_weighted_counts[lev] =
                        (state->factors[g].weighted_counts[lev] > 0) ? 1.0 / state->factors[g].weighted_counts[lev] : 0.0;
                }
            }
        } else {
            state->factors[g].inv_weighted_counts = NULL;
        }
    }

    /* Allocate CG solver buffers */
    ST_int num_threads = omp_get_max_threads();
    if (num_threads > 8) num_threads = 8;

    state->thread_cg_r = (ST_double **)malloc(num_threads * sizeof(ST_double *));
    state->thread_cg_u = (ST_double **)malloc(num_threads * sizeof(ST_double *));
    state->thread_cg_v = (ST_double **)malloc(num_threads * sizeof(ST_double *));
    state->thread_proj = (ST_double **)malloc(num_threads * sizeof(ST_double *));

    /* thread_fe_means needs num_threads * G arrays, one per (thread, factor) combination */
    state->thread_fe_means = (ST_double **)calloc(num_threads * G, sizeof(ST_double *));

    for (ST_int t = 0; t < num_threads; t++) {
        state->thread_cg_r[t] = (ST_double *)malloc(N * sizeof(ST_double));
        state->thread_cg_u[t] = (ST_double *)malloc(N * sizeof(ST_double));
        state->thread_cg_v[t] = (ST_double *)malloc(N * sizeof(ST_double));
        state->thread_proj[t] = (ST_double *)malloc(N * sizeof(ST_double));

        /* Allocate mean arrays for each factor - use num_levels (contiguous after remapping) */
        for (ST_int g = 0; g < G; g++) {
            state->thread_fe_means[t * G + g] = (ST_double *)malloc(
                state->factors[g].num_levels * sizeof(ST_double));
        }
    }

    state->factors_initialized = 1;
    state->num_threads = num_threads;  /* Required for partial_out_columns */

    /* Set global state for creghdfe solver functions */
    g_state = state;

    /* End HDFE setup timing, start FWL timing */
    t_hdfe_setup = ctools_timer_seconds() - t_hdfe_setup_start;
    double t_fwl_start = ctools_timer_seconds();

    /* FWL partialling: partial out specified exogenous variables before FE absorption */
    ST_double dval_partial;
    ST_int n_partial = 0;
    if (SF_scal_use("__civreghdfe_n_partial", &dval_partial) == 0) {
        n_partial = (ST_int)dval_partial;
    }

    if (n_partial > 0 && K_exog > 0) {
        if (verbose) {
            char buf[256];
            snprintf(buf, sizeof(buf), "civreghdfe: Partialling out %d exogenous variable(s) via FWL...\n", (int)n_partial);
            SF_display(buf);
        }

        /* Read partial variable indices (1-based indices into X_exog) */
        ST_int *partial_indices = (ST_int *)malloc(n_partial * sizeof(ST_int));
        ST_int *is_partial = (ST_int *)calloc(K_exog, sizeof(ST_int));  /* Mask for partial vars */
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

        /* Build matrix P of partial variables (N x n_partial) */
        ST_double *P = (ST_double *)malloc((size_t)N * n_partial * sizeof(ST_double));
        for (ST_int pi = 0; pi < n_partial; pi++) {
            ST_int idx = partial_indices[pi] - 1;  /* Convert to 0-based */
            memcpy(P + pi * N, X_exog_c + idx * N, N * sizeof(ST_double));
        }

        /* Compute P'P and invert */
        ST_double *PtP = (ST_double *)calloc(n_partial * n_partial, sizeof(ST_double));
        ST_double *PtP_inv = (ST_double *)calloc(n_partial * n_partial, sizeof(ST_double));

        civreghdfe_matmul_atb(P, P, N, n_partial, n_partial, PtP);
        memcpy(PtP_inv, PtP, n_partial * n_partial * sizeof(ST_double));

        if (cholesky(PtP_inv, n_partial) == 0 && invert_from_cholesky(PtP_inv, n_partial, PtP_inv) == 0) {
            /* Residualize y: y = y - P(P'P)^{-1}P'y */
            ST_double *Pty = (ST_double *)calloc(n_partial, sizeof(ST_double));
            for (ST_int p = 0; p < n_partial; p++) {
                ST_double sum = 0.0;
                for (ST_int i = 0; i < N; i++) {
                    sum += P[p * N + i] * y_c[i];
                }
                Pty[p] = sum;
            }

            ST_double *coef = (ST_double *)calloc(n_partial, sizeof(ST_double));
            for (ST_int p = 0; p < n_partial; p++) {
                ST_double sum = 0.0;
                for (ST_int q = 0; q < n_partial; q++) {
                    sum += PtP_inv[q * n_partial + p] * Pty[q];
                }
                coef[p] = sum;
            }

            for (ST_int i = 0; i < N; i++) {
                ST_double fitted = 0.0;
                for (ST_int p = 0; p < n_partial; p++) {
                    fitted += P[p * N + i] * coef[p];
                }
                y_c[i] -= fitted;
            }

            /* Residualize X_endog */
            for (ST_int k = 0; k < K_endog; k++) {
                ST_double *Ptx = Pty;  /* Reuse */
                for (ST_int p = 0; p < n_partial; p++) {
                    ST_double sum = 0.0;
                    for (ST_int i = 0; i < N; i++) {
                        sum += P[p * N + i] * X_endog_c[k * N + i];
                    }
                    Ptx[p] = sum;
                }

                for (ST_int p = 0; p < n_partial; p++) {
                    ST_double sum = 0.0;
                    for (ST_int q = 0; q < n_partial; q++) {
                        sum += PtP_inv[q * n_partial + p] * Ptx[q];
                    }
                    coef[p] = sum;
                }

                for (ST_int i = 0; i < N; i++) {
                    ST_double fitted = 0.0;
                    for (ST_int p = 0; p < n_partial; p++) {
                        fitted += P[p * N + i] * coef[p];
                    }
                    X_endog_c[k * N + i] -= fitted;
                }
            }

            /* Residualize remaining X_exog (those not being partialled) */
            for (ST_int k = 0; k < K_exog; k++) {
                if (is_partial[k]) continue;  /* Skip partial vars themselves */

                ST_double *Ptx = Pty;
                for (ST_int p = 0; p < n_partial; p++) {
                    ST_double sum = 0.0;
                    for (ST_int i = 0; i < N; i++) {
                        sum += P[p * N + i] * X_exog_c[k * N + i];
                    }
                    Ptx[p] = sum;
                }

                for (ST_int p = 0; p < n_partial; p++) {
                    ST_double sum = 0.0;
                    for (ST_int q = 0; q < n_partial; q++) {
                        sum += PtP_inv[q * n_partial + p] * Ptx[q];
                    }
                    coef[p] = sum;
                }

                for (ST_int i = 0; i < N; i++) {
                    ST_double fitted = 0.0;
                    for (ST_int p = 0; p < n_partial; p++) {
                        fitted += P[p * N + i] * coef[p];
                    }
                    X_exog_c[k * N + i] -= fitted;
                }
            }

            /* Residualize Z */
            for (ST_int k = 0; k < K_iv; k++) {
                ST_double *Ptx = Pty;
                for (ST_int p = 0; p < n_partial; p++) {
                    ST_double sum = 0.0;
                    for (ST_int i = 0; i < N; i++) {
                        sum += P[p * N + i] * Z_c[k * N + i];
                    }
                    Ptx[p] = sum;
                }

                for (ST_int p = 0; p < n_partial; p++) {
                    ST_double sum = 0.0;
                    for (ST_int q = 0; q < n_partial; q++) {
                        sum += PtP_inv[q * n_partial + p] * Ptx[q];
                    }
                    coef[p] = sum;
                }

                for (ST_int i = 0; i < N; i++) {
                    ST_double fitted = 0.0;
                    for (ST_int p = 0; p < n_partial; p++) {
                        fitted += P[p * N + i] * coef[p];
                    }
                    Z_c[k * N + i] -= fitted;
                }
            }

            free(Pty);
            free(coef);

            /* Build reduced X_exog with partialled variables removed */
            ST_int K_exog_new = K_exog - n_partial;
            if (K_exog_new > 0) {
                ST_double *X_exog_reduced = (ST_double *)malloc((size_t)N * K_exog_new * sizeof(ST_double));
                ST_int new_idx = 0;
                for (ST_int k = 0; k < K_exog; k++) {
                    if (!is_partial[k]) {
                        memcpy(X_exog_reduced + new_idx * N, X_exog_c + k * N, N * sizeof(ST_double));
                        new_idx++;
                    }
                }
                free(X_exog_c);
                X_exog_c = X_exog_reduced;
            } else {
                free(X_exog_c);
                X_exog_c = NULL;
            }
            K_exog = K_exog_new;

            if (verbose) {
                char buf[256];
                snprintf(buf, sizeof(buf), "civreghdfe: FWL partialling complete, K_exog reduced to %d\n", (int)K_exog);
                SF_display(buf);
            }
        } else {
            SF_error("civreghdfe: Warning - FWL partialling failed (P'P singular)\n");
        }

        free(partial_indices);
        free(is_partial);
        free(P);
        free(PtP);
        free(PtP_inv);
    }

    /* End FWL timing, start partial out timing */
    t_fwl = ctools_timer_seconds() - t_fwl_start;
    double t_partial_start = ctools_timer_seconds();

    /* Partial out FEs from all variables using CG solver */
    if (verbose) {
        SF_display("civreghdfe: Partialling out fixed effects...\n");
    }

    ST_int total_cols = 1 + K_endog + K_exog + K_iv;

    /* Combine all data into one array for parallel processing */
    ST_double *all_data = (ST_double *)malloc((size_t)N * total_cols * sizeof(ST_double));
    if (!all_data) {
        SF_error("civreghdfe: Memory allocation failed for demeaning buffer\n");
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
        /* Cleanup state */
        for (ST_int t = 0; t < num_threads; t++) {
            free(state->thread_cg_r[t]); free(state->thread_cg_u[t]);
            free(state->thread_cg_v[t]); free(state->thread_proj[t]);
            for (ST_int fg = 0; fg < G; fg++) {
                if (state->thread_fe_means[t * G + fg])
                    free(state->thread_fe_means[t * G + fg]);
            }
        }
        free(state->thread_cg_r); free(state->thread_cg_u);
        free(state->thread_cg_v); free(state->thread_proj);
        free(state->thread_fe_means);
        for (ST_int fg = 0; fg < G; fg++) {
            free(state->factors[fg].counts);
            if (state->factors[fg].inv_counts) free(state->factors[fg].inv_counts);
            if (state->factors[fg].weighted_counts) free(state->factors[fg].weighted_counts);
            if (state->factors[fg].inv_weighted_counts) free(state->factors[fg].inv_weighted_counts);
            free(state->factors[fg].means);
        }
        free(state->factors);
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

    /* Demean all columns in parallel using shared helper */
    partial_out_columns(state, all_data, N, total_cols, num_threads);

    /* End partial out timing, start post-processing timing */
    t_partial_out = ctools_timer_seconds() - t_partial_start;
    double t_postproc_start = ctools_timer_seconds();

    if (verbose) {
        SF_display("civreghdfe: Fixed effects partialled out\n");
    }

    /* Extract demeaned data */
    ST_double *y_dem = all_data;
    ST_double *X_endog_dem = (K_endog > 0) ? all_data + N : NULL;
    ST_double *X_exog_dem = (K_exog > 0) ? all_data + N * (1 + K_endog) : NULL;
    ST_double *Z_dem = all_data + N * (1 + K_endog + K_exog);

    /* Compute df_a (absorbed degrees of freedom) */
    /* Single FE: df_a = num_levels (constant absorbed by FE)
       Multi-FE: df_a = sum(num_levels) - (G - 1) for redundant levels */
    ST_int df_a = 0;
    ST_int df_a_nested = 0;  /* Levels from FE nested within cluster */
    for (ST_int g = 0; g < G; g++) {
        df_a += state->factors[g].num_levels;
        /* Check if this FE is nested in cluster (1-indexed) */
        if (nested_fe_index > 0 && g == (nested_fe_index - 1)) {
            df_a_nested = state->factors[g].num_levels;
        }
        /* Save per-FE num_levels to Stata scalars for absorbed DOF table */
        char scalar_name[64];
        snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_num_levels_%d", (int)(g + 1));
        SF_scal_save(scalar_name, (ST_double)state->factors[g].num_levels);
    }
    if (G > 1) df_a -= (G - 1);  /* Subtract redundant levels for multi-way FE */

    /* For VCE calculation when FE is nested in cluster:
       - nested FE levels don't contribute to df_a for VCE
       - nested_adj = 1 to account for the constant */
    ST_int df_a_for_vce = df_a - df_a_nested;
    ST_int nested_adj = (df_a_nested > 0) ? 1 : 0;

    if (verbose && has_cluster) {
        char buf[512];
        snprintf(buf, sizeof(buf),
            "civreghdfe: nested_fe_index=%d, df_a=%d, df_a_nested=%d, df_a_for_vce=%d, nested_adj=%d\n",
            (int)nested_fe_index, (int)df_a, (int)df_a_nested, (int)df_a_for_vce, (int)nested_adj);
        SF_display(buf);
    }

    /* Remap cluster IDs to contiguous 1-based indices using shared utility */
    ST_int num_clusters = 0;
    if (has_cluster) {
        if (ctools_remap_cluster_ids(cluster_ids_c, N, &num_clusters) != 0) {
            SF_error("civreghdfe: Failed to remap cluster IDs (memory allocation error)\n");
            free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
            free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
            for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
            free(fe_levels_c);
            /* Cleanup state (singleton_mask already freed after singleton removal) */
            for (ST_int t = 0; t < num_threads; t++) {
                free(state->thread_cg_r[t]); free(state->thread_cg_u[t]);
                free(state->thread_cg_v[t]); free(state->thread_proj[t]);
                for (ST_int fg = 0; fg < G; fg++) {
                    if (state->thread_fe_means[t * G + fg])
                        free(state->thread_fe_means[t * G + fg]);
                }
            }
            free(state->thread_cg_r); free(state->thread_cg_u);
            free(state->thread_cg_v); free(state->thread_proj);
            free(state->thread_fe_means);
            for (ST_int fg = 0; fg < G; fg++) {
                free(state->factors[fg].counts);
                if (state->factors[fg].inv_counts) free(state->factors[fg].inv_counts);
                if (state->factors[fg].weighted_counts) free(state->factors[fg].weighted_counts);
                if (state->factors[fg].inv_weighted_counts) free(state->factors[fg].inv_weighted_counts);
                free(state->factors[fg].means);
            }
            free(state->factors);
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
            for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
            free(fe_levels_c);
            /* Cleanup state (singleton_mask already freed after singleton removal) */
            for (ST_int t = 0; t < num_threads; t++) {
                free(state->thread_cg_r[t]); free(state->thread_cg_u[t]);
                free(state->thread_cg_v[t]); free(state->thread_proj[t]);
                for (ST_int fg = 0; fg < G; fg++) {
                    if (state->thread_fe_means[t * G + fg])
                        free(state->thread_fe_means[t * G + fg]);
                }
            }
            free(state->thread_cg_r); free(state->thread_cg_u);
            free(state->thread_cg_v); free(state->thread_proj);
            free(state->thread_fe_means);
            for (ST_int fg = 0; fg < G; fg++) {
                free(state->factors[fg].counts);
                if (state->factors[fg].inv_counts) free(state->factors[fg].inv_counts);
                if (state->factors[fg].weighted_counts) free(state->factors[fg].weighted_counts);
                if (state->factors[fg].inv_weighted_counts) free(state->factors[fg].inv_weighted_counts);
                free(state->factors[fg].means);
            }
            free(state->factors);
            free(state);
            g_state = NULL;
            return 920;
        }
    }

    /* Detect FEs nested within cluster variable
       An FE is nested if every FE level maps to exactly one cluster.
       This is done by checking if for each FE level, all observations have the same cluster_id.
       We use a simple approach: for each FE, track first seen cluster_id per level. */
    ST_int *fe_nested = (ST_int *)calloc(G, sizeof(ST_int));  /* 1 if nested, 0 otherwise */
    if (has_cluster && fe_nested) {
        for (ST_int g = 0; g < G; g++) {
            ST_int num_levels = state->factors[g].num_levels;
            /* Allocate array to store first-seen cluster for each level */
            ST_int *level_cluster = (ST_int *)malloc(num_levels * sizeof(ST_int));
            if (!level_cluster) continue;

            /* Initialize to -1 (not yet seen) */
            for (ST_int lev = 0; lev < num_levels; lev++) {
                level_cluster[lev] = -1;
            }

            /* Check each observation */
            int is_nested = 1;
            for (ST_int i = 0; i < N && is_nested; i++) {
                ST_int lev = fe_levels_c[g][i] - 1;  /* 0-based level index */
                ST_int clust = cluster_ids_c[i];
                if (level_cluster[lev] == -1) {
                    /* First time seeing this level */
                    level_cluster[lev] = clust;
                } else if (level_cluster[lev] != clust) {
                    /* This level maps to multiple clusters - not nested */
                    is_nested = 0;
                }
            }

            fe_nested[g] = is_nested;
            free(level_cluster);

            /* Store nested status to Stata scalar */
            char scalar_name[64];
            snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_fe_nested_%d", (int)(g + 1));
            SF_scal_save(scalar_name, (ST_double)is_nested);
        }

        if (verbose) {
            char buf[256];
            snprintf(buf, sizeof(buf), "civreghdfe: FE nested detection complete\n");
            SF_display(buf);
            for (ST_int g = 0; g < G; g++) {
                snprintf(buf, sizeof(buf), "  FE %d: %s\n", (int)(g + 1),
                         fe_nested[g] ? "nested" : "not nested");
                SF_display(buf);
            }
        }
    }
    if (fe_nested) free(fe_nested);

    /* Update K_total after possible FWL reduction of K_exog */
    K_total = K_exog + K_endog;

    /* Allocate output arrays */
    ST_double *beta = (ST_double *)calloc(K_total, sizeof(ST_double));
    ST_double *V = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *first_stage_F = (ST_double *)calloc(K_endog, sizeof(ST_double));

    /* End post-processing timing, start estimation timing */
    t_postproc = ctools_timer_seconds() - t_postproc_start;
    double t_estimate_start = ctools_timer_seconds();

    /* Compute k-class IV estimation (2SLS, LIML, Fuller, etc.) */
    /* For cluster VCE, pass df_a_for_vce (excluding nested FE levels) and nested_adj */
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
        kernel_type, bw
    );

    if (rc != STATA_OK) {
        /* Cleanup and return */
        free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        free(beta); free(V); free(first_stage_F);
        /* Cleanup state */
        for (ST_int t = 0; t < num_threads; t++) {
            free(state->thread_cg_r[t]); free(state->thread_cg_u[t]);
            free(state->thread_cg_v[t]); free(state->thread_proj[t]);
            for (ST_int g = 0; g < G; g++) {
                if (state->thread_fe_means[t * G + g])
                    free(state->thread_fe_means[t * G + g]);
            }
        }
        free(state->thread_cg_r); free(state->thread_cg_u);
        free(state->thread_cg_v); free(state->thread_proj);
        free(state->thread_fe_means);
        for (ST_int g = 0; g < G; g++) {
            free(state->factors[g].levels);  /* Same as fe_levels_c[g] */
            free(state->factors[g].counts);
            if (state->factors[g].inv_counts) free(state->factors[g].inv_counts);
            if (state->factors[g].weighted_counts) free(state->factors[g].weighted_counts);
            if (state->factors[g].inv_weighted_counts) free(state->factors[g].inv_weighted_counts);
            free(state->factors[g].means);
        }
        free(state->factors);
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
    ST_int df_r_val = N - K_total - df_a;
    if (df_r_val <= 0) df_r_val = 1;
    ST_double rmse = sqrt(rss / df_r_val);

    /* Compute model F-statistic: (R2 / K) / ((1 - R2) / df_r) */
    ST_double f_stat = 0.0;
    if (K_total > 0 && r2 < 1.0) {
        f_stat = (r2 / K_total) / ((1.0 - r2) / df_r_val);
    }

    /* End stats timing, start store timing */
    double t_stats = ctools_timer_seconds() - t_stats_start;
    double t_store_start = ctools_timer_seconds();

    /* Store results to Stata */
    /* Scalars */
    SF_scal_save("__civreghdfe_N", (ST_double)N);
    SF_scal_save("__civreghdfe_df_r", (ST_double)df_r_val);
    SF_scal_save("__civreghdfe_df_a", (ST_double)df_a);
    SF_scal_save("__civreghdfe_K", (ST_double)K_total);
    SF_scal_save("__civreghdfe_rss", rss);
    SF_scal_save("__civreghdfe_tss", tss);
    SF_scal_save("__civreghdfe_r2", r2);
    SF_scal_save("__civreghdfe_rmse", rmse);
    SF_scal_save("__civreghdfe_F", f_stat);

    if (has_cluster) {
        SF_scal_save("__civreghdfe_N_clust", (ST_double)num_clusters);
    }

    if (has_cluster2) {
        SF_scal_save("__civreghdfe_N_clust2", (ST_double)num_clusters2);
    }

    /* Store lambda for LIML/Fuller */
    if (est_method == 1 || est_method == 2) {
        SF_scal_save("__civreghdfe_lambda", lambda);
    }

    /* Matrices: e(b) and e(V) - reorder from [exog, endog] to [endog, exog] */
    /* Internal order: beta[0..K_exog-1] = exog, beta[K_exog..K_total-1] = endog
       Stata order: [endog, exog] to match ivreghdfe convention */

    /* Store beta in [endog, exog] order */
    for (ST_int e = 0; e < K_endog; e++) {
        SF_mat_store("__civreghdfe_b", 1, e + 1, beta[K_exog + e]);
    }
    for (ST_int x = 0; x < K_exog; x++) {
        SF_mat_store("__civreghdfe_b", 1, K_endog + x + 1, beta[x]);
    }

    /* Store V in [endog, exog] order - need to reorder rows and columns */
    /* new_order: [endog_0, ..., endog_{Ke-1}, exog_0, ..., exog_{Kx-1}]
       For new index i in [0, K_endog): old index = K_exog + i
       For new index i in [K_endog, K_total): old index = i - K_endog */
    for (ST_int new_i = 0; new_i < K_total; new_i++) {
        ST_int old_i = (new_i < K_endog) ? (K_exog + new_i) : (new_i - K_endog);
        for (ST_int new_j = 0; new_j < K_total; new_j++) {
            ST_int old_j = (new_j < K_endog) ? (K_exog + new_j) : (new_j - K_endog);
            SF_mat_store("__civreghdfe_V", new_i + 1, new_j + 1, V[old_j * K_total + old_i]);
        }
    }

    /* First stage F-stats */
    for (ST_int e = 0; e < K_endog; e++) {
        char name[64];
        snprintf(name, sizeof(name), "__civreghdfe_F1_%d", (int)(e + 1));
        SF_scal_save(name, first_stage_F[e]);
    }

    /* End store timing, compute total */
    t_store = ctools_timer_seconds() - t_store_start;
    double t_total = ctools_timer_seconds() - t_total_start;

    /* Save timing scalars */
    SF_scal_save("_civreghdfe_time_load", t_load + t_missing);  /* Combine load and missing check */
    SF_scal_save("_civreghdfe_time_singleton", t_singleton);
    SF_scal_save("_civreghdfe_time_setup", t_hdfe_setup);
    SF_scal_save("_civreghdfe_time_fwl", t_fwl);
    SF_scal_save("_civreghdfe_time_partial", t_partial_out);
    SF_scal_save("_civreghdfe_time_postproc", t_postproc);
    SF_scal_save("_civreghdfe_time_estimate", t_estimate);
    SF_scal_save("_civreghdfe_time_stats", t_stats);
    SF_scal_save("_civreghdfe_time_store", t_store);
    SF_scal_save("_civreghdfe_time_total", t_total);
    CTOOLS_SAVE_THREAD_INFO("_civreghdfe");

    /* Cleanup */
    free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
    free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
    free(beta); free(V); free(first_stage_F);

    /* Cleanup state */
    for (ST_int t = 0; t < num_threads; t++) {
        free(state->thread_cg_r[t]); free(state->thread_cg_u[t]);
        free(state->thread_cg_v[t]); free(state->thread_proj[t]);
        for (ST_int g = 0; g < G; g++) {
            if (state->thread_fe_means[t * G + g])
                free(state->thread_fe_means[t * G + g]);
        }
    }
    free(state->thread_cg_r); free(state->thread_cg_u);
    free(state->thread_cg_v); free(state->thread_proj);
    free(state->thread_fe_means);
    for (ST_int g = 0; g < G; g++) {
        free(state->factors[g].levels);  /* Same as fe_levels_c[g] */
        free(state->factors[g].counts);
        if (state->factors[g].inv_counts) free(state->factors[g].inv_counts);
        if (state->factors[g].weighted_counts) free(state->factors[g].weighted_counts);
        if (state->factors[g].inv_weighted_counts) free(state->factors[g].inv_weighted_counts);
        free(state->factors[g].means);
    }
    free(state->factors);
    free(fe_levels_c);  /* Free the pointer array */
    free(state);
    g_state = NULL;

    if (verbose) {
        SF_display("civreghdfe: Done\n");
    }

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

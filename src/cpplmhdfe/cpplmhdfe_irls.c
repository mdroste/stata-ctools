/*
 * cpplmhdfe_irls.c
 *
 * IRLS (Iteratively Reweighted Least Squares) loop for PPML estimation
 * with high-dimensional fixed effects.
 *
 * Algorithm:
 *   1. Load data, set up FE, remove singletons
 *   2. Detect separation (FE-level: all y=0 in a group)
 *   3. Initialize beta from OLS of log(max(y,1)) on partialled X
 *   4. IRLS loop:
 *      a. Compute IRLS weights: w_irls = mu * w_user
 *      b. Compute working depvar: z = eta + (y - mu)/mu
 *      c. Partial out [z, X] via CG solver with IRLS weights
 *      d. Solve WLS: beta = (X~'X~)^{-1} X~'z~
 *      e. Update eta, mu; check deviance convergence
 *   5. VCE: sandwich with score residuals (y - mu)
 *   6. Store results
 *
 * Part of the ctools Stata plugin suite
 */

#include "cpplmhdfe_irls.h"
#include "cpplmhdfe_separation.h"
#include "../creghdfe/creghdfe_utils.h"
#include "../creghdfe/creghdfe_solver.h"
#include "../creghdfe/creghdfe_vce.h"
#include "../ctools_matrix.h"
#include "../ctools_hdfe_utils.h"
#include "../ctools_config.h"
#include "../ctools_types.h"
#include "../ctools_spi.h"
#include "../ctools_ols.h"
#include "../ctools_matrix.h"

/* Use backward compat names */
#define cholesky ctools_cholesky
#define invert_from_cholesky ctools_invert_from_cholesky

/* Forward declarations for helpers defined at bottom of file */
static void ppml_update_fe_weights(HDFE_State *S, const ST_double *irls_w, ST_int N);
static ST_double ppml_compute_deviance(const ST_double *y, const ST_double *mu,
                                        const ST_double *w_user, ST_int N);
static ST_double ppml_compute_loglik(const ST_double *y, const ST_double *mu,
                                      const ST_double *w_user, ST_int N);

/* Global state pointer for CG solver (same pattern as creghdfe) */
static HDFE_State *g_ppml_state = NULL;

static void cleanup_ppml_state(void)
{
    if (g_ppml_state) {
        ctools_hdfe_state_cleanup(g_ppml_state);
        free(g_ppml_state);
        g_ppml_state = NULL;
    }
}

/* ========================================================================
 * Main PPML regression
 * ======================================================================== */

ST_retcode do_ppml_regression(int argc, char *argv[])
{
    ST_int K, G, N_orig, N, in1, in2;
    ST_int k, g, i, j, idx;
    ST_double val;
    ST_int verbose;
    ST_int num_singletons;
    ST_int *mask = NULL;
    PPMLFactorData *factors = NULL;
    char scalar_name[64];
    double t_start, t_load, t_remap, t_singleton, t_dof, t_irls, t_vce;
    ST_int max_iter_singleton = 100;
    ST_int mobility_groups = 1;
    ST_int df_a = 0;
    ST_int num_threads;

    /* Data arrays */
    ST_double *data = NULL;     /* N x K matrix (column-major) */
    ST_int *cluster_ids = NULL;
    ST_int num_clusters = 0;
    perm_idx_t *obs_map = NULL;

    (void)argc;
    (void)argv;

    t_start = get_time_sec();

    in1 = SF_in1();
    in2 = SF_in2();

    /* Read parameters */
    SF_scal_use("__cpplmhdfe_K", &val); K = (ST_int)val;
    SF_scal_use("__cpplmhdfe_G", &val); G = (ST_int)val;
    SF_scal_use("__cpplmhdfe_verbose", &val); verbose = (ST_int)val;

    ST_int maxiter;
    ST_double tolerance;
    SF_scal_use("__cpplmhdfe_maxiter", &val); maxiter = (ST_int)val;
    SF_scal_use("__cpplmhdfe_tolerance", &val); tolerance = val;

    ST_int vcetype;
    SF_scal_use("__cpplmhdfe_vce_type", &val); vcetype = (ST_int)val;

    ST_int has_weights = 0, weight_type = 0;
    if (SF_scal_use("__cpplmhdfe_has_weights", &val) == 0)
        has_weights = (ST_int)val;
    if (SF_scal_use("__cpplmhdfe_weight_type", &val) == 0)
        weight_type = (ST_int)val;

    ST_int has_offset = 0;
    if (SF_scal_use("__cpplmhdfe_has_offset", &val) == 0)
        has_offset = (ST_int)val;

    /* IRLS parameters */
    ST_int irls_maxiter = 1000;
    ST_double irls_tol = 1e-8;
    ST_double sep_tol = 1e-8;
    if (SF_scal_use("__cpplmhdfe_irls_maxiter", &val) == 0)
        irls_maxiter = (ST_int)val;
    if (SF_scal_use("__cpplmhdfe_irls_tol", &val) == 0)
        irls_tol = val;
    if (SF_scal_use("__cpplmhdfe_sep_tol", &val) == 0)
        sep_tol = val;

    /* DOF adjustment */
    ST_int dof_adjust_type = 0;
    if (SF_scal_use("__cpplmhdfe_dof_adjust_type", &val) == 0)
        dof_adjust_type = (ST_int)val;

    /* Validation */
    if (G < 1 || G > 10) {
        SF_error("cpplmhdfe: invalid number of FE groups (must be 1-10)\n");
        return 198;
    }
    if (K < 2) {
        SF_error("cpplmhdfe: need at least depvar and one indepvar (K < 2)\n");
        return 198;
    }

    /* Determine threads */
#ifdef _OPENMP
    num_threads = ctools_get_max_threads();
    if (num_threads > K) num_threads = K;
    if (num_threads < 1) num_threads = 1;
#else
    num_threads = 1;
#endif

    /* ================================================================
     * STEP 1: Data loading
     * varlist: depvar indepvars... fe1..feG [cluster] [weight] [offset]
     * ================================================================ */
    ST_int total_vars = K + G + (vcetype == 2 ? 1 : 0) + (has_weights ? 1 : 0) + (has_offset ? 1 : 0);

    int *var_indices = (int *)malloc(total_vars * sizeof(int));
    if (!var_indices) {
        SF_error("cpplmhdfe: memory allocation failed\n");
        return 920;
    }
    for (i = 0; i < total_vars; i++) {
        var_indices[i] = i + 1;
    }

    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);

    stata_retcode load_rc = ctools_data_load(&filtered, var_indices, total_vars, 0, 0, 0);
    free(var_indices);

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        SF_error("cpplmhdfe: parallel data load failed\n");
        return 920;
    }

    N_orig = (ST_int)filtered.data.nobs;
    obs_map = filtered.obs_map;

    if (N_orig <= 0) {
        ctools_filtered_data_free(&filtered);
        SF_error("cpplmhdfe: no observations\n");
        return 2000;
    }

    t_load = get_time_sec();

    /* Allocate working arrays */
    factors = (PPMLFactorData *)calloc(G, sizeof(PPMLFactorData));
    mask = (ST_int *)malloc(N_orig * sizeof(ST_int));
    data = (ST_double *)malloc((size_t)N_orig * K * sizeof(ST_double));

    ST_double *weights = NULL;
    if (has_weights) {
        weights = (ST_double *)malloc(N_orig * sizeof(ST_double));
    }

    ST_double *offset_arr = NULL;
    if (has_offset) {
        offset_arr = (ST_double *)malloc(N_orig * sizeof(ST_double));
    }

    if (!factors || !mask || !data || (has_weights && !weights) || (has_offset && !offset_arr)) {
        if (factors) free(factors);
        if (mask) free(mask);
        if (data) free(data);
        if (weights) free(weights);
        if (offset_arr) free(offset_arr);
        ctools_filtered_data_free(&filtered);
        SF_error("cpplmhdfe: memory allocation failed\n");
        return 920;
    }

    for (i = 0; i < N_orig; i++) mask[i] = 1;

    /* Per-factor level arrays */
    for (g = 0; g < G; g++) {
        factors[g].levels = (ST_int *)malloc(N_orig * sizeof(ST_int));
        factors[g].num_obs = N_orig;
        factors[g].counts = NULL;
        factors[g].num_levels = 0;
        if (!factors[g].levels) {
            for (i = 0; i < g; i++) free(factors[i].levels);
            free(factors); free(mask); free(data);
            if (weights) free(weights);
            if (offset_arr) free(offset_arr);
            ctools_filtered_data_free(&filtered);
            return 920;
        }
    }

    /* Copy numeric data to column-major */
    for (k = 0; k < K; k++) {
        double *src = filtered.data.vars[k].data.dbl;
        double *dst = &data[k * N_orig];
        memcpy(dst, src, N_orig * sizeof(double));
    }

    /* Copy weights */
    if (has_weights) {
        ST_int weight_pos = K + G + (vcetype == 2 ? 1 : 0);
        memcpy(weights, filtered.data.vars[weight_pos].data.dbl, N_orig * sizeof(ST_double));
    }

    /* Copy offset */
    if (has_offset) {
        ST_int offset_pos = K + G + (vcetype == 2 ? 1 : 0) + (has_weights ? 1 : 0);
        memcpy(offset_arr, filtered.data.vars[offset_pos].data.dbl, N_orig * sizeof(ST_double));
    }

    /* Validate y >= 0 */
    {
        ST_double *y = data;  /* First column */
        for (i = 0; i < N_orig; i++) {
            if (y[i] < 0.0) {
                SF_error("cpplmhdfe: dependent variable must be non-negative\n");
                for (g = 0; g < G; g++) free(factors[g].levels);
                free(factors); free(mask); free(data);
                if (weights) free(weights);
                if (offset_arr) free(offset_arr);
                ctools_filtered_data_free(&filtered);
                return 198;
            }
        }
    }

    /* FE remap + counting */
    ST_double *weighted_counts_orig[10] = {NULL};
    for (g = 0; g < G; g++) {
        ST_int *counts_g = NULL;
        ST_double *wcounts_g = NULL;
        if (remap_and_count(filtered.data.vars[K + g].data.dbl, N_orig,
                            factors[g].levels, &factors[g].num_levels,
                            &counts_g,
                            has_weights ? weights : NULL,
                            has_weights ? &wcounts_g : NULL) == 0) {
            factors[g].counts = counts_g;
            if (has_weights && wcounts_g)
                weighted_counts_orig[g] = wcounts_g;
        }
    }

    /* Save cluster raw values before freeing filtered data */
    double *cluster_raw_values = NULL;
    ST_int cluster_matches_fe = -1;
    if (vcetype == 2) {
        cluster_raw_values = (double *)malloc(N_orig * sizeof(double));
        if (cluster_raw_values) {
            memcpy(cluster_raw_values, filtered.data.vars[K + G].data.dbl, N_orig * sizeof(double));
            /* Check for FE match */
            double *cluster_data = filtered.data.vars[K + G].data.dbl;
            for (g = 0; g < G; g++) {
                double *fe_data = filtered.data.vars[K + g].data.dbl;
                if ((ST_int)cluster_data[0] != (ST_int)fe_data[0] ||
                    (ST_int)cluster_data[N_orig/2] != (ST_int)fe_data[N_orig/2] ||
                    (ST_int)cluster_data[N_orig-1] != (ST_int)fe_data[N_orig-1])
                    continue;
                ST_int all_match = 1;
                for (idx = 0; idx < N_orig && all_match; idx++) {
                    if ((ST_int)cluster_data[idx] != (ST_int)fe_data[idx])
                        all_match = 0;
                }
                if (all_match) { cluster_matches_fe = g; break; }
            }
        }
    }

    /* Free filtered data, keep obs_map */
    stata_data_free(&filtered.data);

    /* Verify remapping */
    for (g = 0; g < G; g++) {
        if (factors[g].num_levels == 0) {
            for (i = 0; i < G; i++) {
                free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
            }
            free(factors); free(mask); free(data);
            if (weights) free(weights);
            if (offset_arr) free(offset_arr);
            if (cluster_raw_values) free(cluster_raw_values);
            SF_error("cpplmhdfe: FE remapping failed\n");
            return 920;
        }
    }

    t_remap = get_time_sec();

    /* ================================================================
     * STEP 2: Singleton removal
     * ================================================================ */
    num_singletons = 0;
    N = N_orig;

    {
        ST_int **fe_levels = (ST_int **)malloc(G * sizeof(ST_int *));
        if (fe_levels) {
            for (g = 0; g < G; g++)
                fe_levels[g] = factors[g].levels;
            num_singletons = ctools_remove_singletons(fe_levels, G, N_orig, mask, max_iter_singleton, (verbose >= 1));
            free(fe_levels);
            N = 0;
            for (i = 0; i < N_orig; i++)
                if (mask[i]) N++;
        }
    }

    if (N == 0) {
        ctools_scal_save("__cpplmhdfe_N", 0.0);
        ctools_scal_save("__cpplmhdfe_num_singletons", (ST_double)num_singletons);
        SF_error("cpplmhdfe: all observations are singletons\n");
        for (g = 0; g < G; g++) {
            free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask); free(data);
        if (weights) free(weights);
        if (offset_arr) free(offset_arr);
        if (cluster_raw_values) free(cluster_raw_values);
        if (obs_map) free(obs_map);
        return 2001;
    }

    t_singleton = get_time_sec();

    /* Recount levels after singleton removal */
    ST_int orig_num_levels[10];
    for (g = 0; g < G; g++) {
        orig_num_levels[g] = factors[g].num_levels;
        memset(factors[g].counts, 0, orig_num_levels[g] * sizeof(ST_int));
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) {
                ST_int level = factors[g].levels[i] - 1;
                if (level >= 0 && level < orig_num_levels[g])
                    factors[g].counts[level]++;
            }
        }
        ST_int nlev = 0;
        for (i = 0; i < orig_num_levels[g]; i++)
            if (factors[g].counts[i] > 0) nlev++;
        factors[g].num_levels = nlev;
    }

    /* ================================================================
     * STEP 3: DOF computation
     * ================================================================ */
    df_a = 0;
    for (g = 0; g < G; g++)
        df_a += factors[g].num_levels;

    if (dof_adjust_type == 1) {
        if (G >= 2) {
            mobility_groups = 1;
            df_a -= 1;
        }
    } else if (G >= 2) {
        ST_int *fe1_c = (ST_int *)malloc(N * sizeof(ST_int));
        ST_int *fe2_c = (ST_int *)malloc(N * sizeof(ST_int));
        if (fe1_c && fe2_c) {
            ST_int *remap1 = (ST_int *)calloc(orig_num_levels[0] + 1, sizeof(ST_int));
            ST_int *remap2 = (ST_int *)calloc(orig_num_levels[1] + 1, sizeof(ST_int));
            if (remap1 && remap2) {
                ST_int next1 = 1, next2 = 1;
                idx = 0;
                for (i = 0; i < N_orig; i++) {
                    if (mask[i]) {
                        ST_int lev1 = factors[0].levels[i];
                        ST_int lev2 = factors[1].levels[i];
                        if (remap1[lev1] == 0) remap1[lev1] = next1++;
                        if (remap2[lev2] == 0) remap2[lev2] = next2++;
                        fe1_c[idx] = remap1[lev1];
                        fe2_c[idx] = remap2[lev2];
                        idx++;
                    }
                }
                mobility_groups = count_connected_components(
                    fe1_c, fe2_c, N, factors[0].num_levels, factors[1].num_levels);
                if (mobility_groups < 0) mobility_groups = 1;
                free(remap1); free(remap2);
            }
        }
        if (fe1_c) free(fe1_c);
        if (fe2_c) free(fe2_c);
        df_a -= mobility_groups;
        if (G > 2 && (dof_adjust_type == 0 || dof_adjust_type == 3)) {
            ST_int extra = G - 2;
            df_a -= extra;
            mobility_groups += extra;
        }
    } else {
        mobility_groups = 0;
    }

    t_dof = get_time_sec();

    /* ================================================================
     * STEP 4: Build HDFE state (compacted)
     * ================================================================ */
    cleanup_ppml_state();
    g_ppml_state = (HDFE_State *)calloc(1, sizeof(HDFE_State));
    if (!g_ppml_state) {
        for (g = 0; g < G; g++) {
            free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask); free(data);
        if (weights) free(weights);
        if (offset_arr) free(offset_arr);
        if (cluster_raw_values) free(cluster_raw_values);
        if (obs_map) free(obs_map);
        return 920;
    }

    g_ppml_state->G = G;
    g_ppml_state->N = N;
    g_ppml_state->K = K;
    g_ppml_state->in1 = in1;
    g_ppml_state->in2 = in2;
    g_ppml_state->maxiter = maxiter;
    g_ppml_state->tolerance = tolerance;
    g_ppml_state->verbose = verbose;
    g_ppml_state->num_threads = num_threads;
    g_ppml_state->factors_initialized = 1;
    g_ppml_state->df_a = df_a;
    g_ppml_state->mobility_groups = mobility_groups;
    g_ppml_state->has_weights = 1;  /* IRLS always uses weights */
    g_ppml_state->weight_type = 1;  /* Treat as aweight for CG solver */
    g_ppml_state->weights = NULL;
    g_ppml_state->sum_weights = 0.0;

    /* Allocate and remap factors */
    g_ppml_state->factors = (FE_Factor *)calloc(G, sizeof(FE_Factor));
    if (!g_ppml_state->factors) {
        cleanup_ppml_state();
        for (g = 0; g < G; g++) {
            free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask); free(data);
        if (weights) free(weights);
        if (offset_arr) free(offset_arr);
        if (cluster_raw_values) free(cluster_raw_values);
        if (obs_map) free(obs_map);
        return 920;
    }

    for (g = 0; g < G; g++) {
        g_ppml_state->factors[g].num_levels = factors[g].num_levels;
        g_ppml_state->factors[g].max_level = factors[g].num_levels - 1;
        g_ppml_state->factors[g].has_intercept = 1;
        g_ppml_state->factors[g].levels = (ST_int *)malloc(N * sizeof(ST_int));
        g_ppml_state->factors[g].counts = (ST_double *)calloc(factors[g].num_levels, sizeof(ST_double));
        g_ppml_state->factors[g].weighted_counts = (ST_double *)calloc(factors[g].num_levels, sizeof(ST_double));
        g_ppml_state->factors[g].means = NULL;

        if (!g_ppml_state->factors[g].levels || !g_ppml_state->factors[g].counts ||
            !g_ppml_state->factors[g].weighted_counts) {
            cleanup_ppml_state();
            for (i = 0; i < G; i++) {
                free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
            }
            free(factors); free(mask); free(data);
            if (weights) free(weights);
            if (offset_arr) free(offset_arr);
            for (i = 0; i < G; i++)
                if (weighted_counts_orig[i]) free(weighted_counts_orig[i]);
            if (cluster_raw_values) free(cluster_raw_values);
            if (obs_map) free(obs_map);
            return 920;
        }

        /* Remap to contiguous levels */
        ST_int *remap = (ST_int *)calloc(orig_num_levels[g] + 1, sizeof(ST_int));
        if (!remap) {
            cleanup_ppml_state();
            for (i = 0; i < G; i++) {
                free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
            }
            free(factors); free(mask); free(data);
            if (weights) free(weights);
            if (offset_arr) free(offset_arr);
            for (i = 0; i < G; i++)
                if (weighted_counts_orig[i]) free(weighted_counts_orig[i]);
            if (cluster_raw_values) free(cluster_raw_values);
            if (obs_map) free(obs_map);
            return 920;
        }
        ST_int next_level = 1;
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) {
                ST_int old_level = factors[g].levels[i];
                if (remap[old_level] == 0)
                    remap[old_level] = next_level++;
                g_ppml_state->factors[g].levels[idx] = remap[old_level];
                g_ppml_state->factors[g].counts[remap[old_level] - 1] += 1.0;
                idx++;
            }
        }
        free(remap);

        ctools_build_sorted_permutation(&g_ppml_state->factors[g], N);
    }

    /* Free original weighted counts */
    for (g = 0; g < G; g++) {
        if (weighted_counts_orig[g]) free(weighted_counts_orig[g]);
    }

    /* ================================================================
     * STEP 5: Compact data (remove singleton rows)
     * ================================================================ */
    ST_double *data_compact = (ST_double *)malloc((size_t)N * K * sizeof(ST_double));
    ST_double *w_user_compact = NULL;
    ST_double *offset_compact = NULL;

    if (has_weights) w_user_compact = (ST_double *)malloc(N * sizeof(ST_double));
    if (has_offset) offset_compact = (ST_double *)malloc(N * sizeof(ST_double));

    if (!data_compact || (has_weights && !w_user_compact) || (has_offset && !offset_compact)) {
        if (data_compact) free(data_compact);
        if (w_user_compact) free(w_user_compact);
        if (offset_compact) free(offset_compact);
        cleanup_ppml_state();
        for (g = 0; g < G; g++) {
            free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask); free(data);
        if (weights) free(weights);
        if (offset_arr) free(offset_arr);
        if (cluster_raw_values) free(cluster_raw_values);
        if (obs_map) free(obs_map);
        return 920;
    }

    /* Compact data columns */
    for (k = 0; k < K; k++) {
        const ST_double *src_col = data + k * N_orig;
        ST_double *dst_col = data_compact + k * N;
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) {
                dst_col[idx++] = src_col[i];
            }
        }
    }

    /* Compact weights */
    if (has_weights) {
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) w_user_compact[idx++] = weights[i];
        }
    }

    /* Compact offset */
    if (has_offset) {
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) offset_compact[idx++] = offset_arr[i];
        }
    }

    free(data);
    data = data_compact;
    if (weights) free(weights);
    weights = NULL;
    if (offset_arr) free(offset_arr);
    offset_arr = NULL;

    /* ================================================================
     * STEP 6: Separation detection (FE-level)
     * ================================================================ */
    ST_int num_separated = 0;
    {
        ST_int *sep_mask = (ST_int *)calloc(N, sizeof(ST_int));
        if (sep_mask) {
            num_separated = ppml_detect_separation_fe(
                data, g_ppml_state->factors, G, N, sep_mask);

            if (num_separated > 0 && verbose) {
                char msg[128];
                snprintf(msg, sizeof(msg), "(cpplmhdfe: %d separated observations detected and dropped)\n", num_separated);
                SF_display(msg);
            }

            /* Remove separated observations if any */
            if (num_separated > 0) {
                /* Build new mask: invert sep_mask */
                ST_int *keep = (ST_int *)malloc(N * sizeof(ST_int));
                if (keep) {
                    for (i = 0; i < N; i++) keep[i] = !sep_mask[i];
                    ST_int N_new = N - num_separated;

                    /* Compact data */
                    ST_double *data_new = (ST_double *)malloc((size_t)N_new * K * sizeof(ST_double));
                    if (data_new) {
                        for (k = 0; k < K; k++) {
                            idx = 0;
                            for (i = 0; i < N; i++) {
                                if (keep[i])
                                    data_new[k * N_new + idx++] = data[k * N + i];
                            }
                        }
                        free(data);
                        data = data_new;
                    }

                    /* Compact w_user */
                    if (w_user_compact) {
                        ST_double *w_new = (ST_double *)malloc(N_new * sizeof(ST_double));
                        if (w_new) {
                            idx = 0;
                            for (i = 0; i < N; i++)
                                if (keep[i]) w_new[idx++] = w_user_compact[i];
                            free(w_user_compact);
                            w_user_compact = w_new;
                        }
                    }

                    /* Compact offset */
                    if (offset_compact) {
                        ST_double *o_new = (ST_double *)malloc(N_new * sizeof(ST_double));
                        if (o_new) {
                            idx = 0;
                            for (i = 0; i < N; i++)
                                if (keep[i]) o_new[idx++] = offset_compact[i];
                            free(offset_compact);
                            offset_compact = o_new;
                        }
                    }

                    /* Compact FE levels */
                    for (g = 0; g < G; g++) {
                        ST_int *lev_new = (ST_int *)malloc(N_new * sizeof(ST_int));
                        if (lev_new) {
                            idx = 0;
                            for (i = 0; i < N; i++)
                                if (keep[i]) lev_new[idx++] = g_ppml_state->factors[g].levels[i];
                            free(g_ppml_state->factors[g].levels);
                            g_ppml_state->factors[g].levels = lev_new;
                            /* Rebuild sorted permutation */
                            ctools_build_sorted_permutation(&g_ppml_state->factors[g], N_new);
                        }
                    }

                    /* Compact obs_map */
                    if (obs_map) {
                        /* obs_map maps compacted index → original Stata obs.
                         * But after singleton removal, obs_map is still N_orig size.
                         * We need to build a new mapping for the post-separation N_new obs.
                         * The current obs_map was already compacted by singleton mask.
                         * We need to re-compact using the keep array. */
                        /* First, build the obs_map for post-singleton N observations */
                        perm_idx_t *obs_map_n = (perm_idx_t *)malloc(N * sizeof(perm_idx_t));
                        if (obs_map_n) {
                            idx = 0;
                            for (i = 0; i < N_orig; i++) {
                                if (mask[i]) {
                                    obs_map_n[idx++] = obs_map[i];
                                }
                            }
                            /* Now compact obs_map_n by keep */
                            perm_idx_t *obs_map_new = (perm_idx_t *)malloc(N_new * sizeof(perm_idx_t));
                            if (obs_map_new) {
                                idx = 0;
                                for (i = 0; i < N; i++) {
                                    if (keep[i]) obs_map_new[idx++] = obs_map_n[i];
                                }
                                free(obs_map);
                                obs_map = obs_map_new;
                                /* Update mask to reflect the new mapping */
                                free(mask);
                                mask = (ST_int *)malloc(N_new * sizeof(ST_int));
                                if (mask) {
                                    for (i = 0; i < N_new; i++) mask[i] = 1;
                                }
                            }
                            free(obs_map_n);
                        }
                    }

                    N = N_new;
                    g_ppml_state->N = N;
                    free(keep);
                }
            }
            free(sep_mask);
        }
    }

    /* ================================================================
     * STEP 7: Allocate IRLS arrays and HDFE buffers
     * ================================================================ */
    ST_double *mu = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *eta = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *irls_w = (ST_double *)malloc(N * sizeof(ST_double));
    /* z_orig: saves pre-partialled working depvar for FE reconstruction.
     * aug_copy: [z, X1, ..., Xk] for CG partialling (modified in place). */
    ST_int K_x = K - 1;
    ST_int K_aug = K_x + 1;  /* working depvar + X vars */
    ST_double *z_orig = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *aug_copy = (ST_double *)malloc((size_t)N * K_aug * sizeof(ST_double));

    /* Pre-allocate IRLS working buffers (reused across iterations) */
    ST_double *xtx = (ST_double *)malloc(K_x * K_x * sizeof(ST_double));
    ST_double *xty = (ST_double *)malloc(K_x * sizeof(ST_double));
    ST_double *irls_beta = (ST_double *)malloc(K_x * sizeof(ST_double));
    ST_double *inv_xx = (ST_double *)malloc(K_x * K_x * sizeof(ST_double));
    ST_int *sep_m = (ST_int *)malloc(N * sizeof(ST_int));

    if (!mu || !eta || !irls_w || !z_orig || !aug_copy ||
        !xtx || !xty || !irls_beta || !inv_xx || !sep_m) {
        if (mu) free(mu);
        if (eta) free(eta);
        if (irls_w) free(irls_w);
        if (z_orig) free(z_orig);
        if (aug_copy) free(aug_copy);
        if (xtx) free(xtx);
        if (xty) free(xty);
        if (irls_beta) free(irls_beta);
        if (inv_xx) free(inv_xx);
        if (sep_m) free(sep_m);
        cleanup_ppml_state();
        for (g = 0; g < G; g++) {
            free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); if (mask) free(mask); free(data);
        if (w_user_compact) free(w_user_compact);
        if (offset_compact) free(offset_compact);
        if (cluster_raw_values) free(cluster_raw_values);
        if (obs_map) free(obs_map);
        return 920;
    }

    /* Set IRLS weights for HDFE state */
    g_ppml_state->weights = irls_w;

    /* Allocate HDFE buffers (alloc_proj=0, max_columns=K_aug) */
    if (ctools_hdfe_alloc_buffers(g_ppml_state, 0, K_aug) != 0) {
        free(mu); free(eta);
        free(z_orig); free(aug_copy);
        free(xtx); free(xty); free(irls_beta); free(inv_xx); free(sep_m);
        /* irls_w is owned by g_ppml_state->weights, freed by cleanup */
        cleanup_ppml_state();
        for (g = 0; g < G; g++) {
            free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); if (mask) free(mask); free(data);
        if (w_user_compact) free(w_user_compact);
        if (offset_compact) free(offset_compact);
        if (cluster_raw_values) free(cluster_raw_values);
        if (obs_map) free(obs_map);
        return 920;
    }

    /* ================================================================
     * STEP 8: Initialize IRLS
     * Initialize eta from OLS of log(max(y,1)) on X with unit weights
     * ================================================================ */
    ST_double *y = data;  /* First column of data */
    ST_double *X = data + N;  /* Columns 1..K-1 */

    /* Standardize y and X to match ppmlhdfe's convergence path.
     * ppmlhdfe divides by sample stdev (with floor 1e-3) before IRLS.
     * The VCE is invariant to standardization algebraically, but matching
     * the convergence trajectory gives the same terminal state, eliminating
     * systematic VCE precision differences. */
    ST_double stdev_y = 1.0;
    ST_double *stdev_x = (ST_double *)malloc(K_x * sizeof(ST_double));
    {
        ST_double sum_y = 0.0, sum_y2 = 0.0;
        for (i = 0; i < N; i++) {
            sum_y += y[i];
            sum_y2 += y[i] * y[i];
        }
        ST_double var_y = (sum_y2 - sum_y * sum_y / N) / (N - 1);
        stdev_y = (var_y > 0.0) ? sqrt(var_y) : 1.0;
        if (stdev_y < 1e-3) stdev_y = 1e-3;
        for (i = 0; i < N; i++)
            y[i] /= stdev_y;
    }
    if (stdev_x) {
        for (k = 0; k < K_x; k++) {
            ST_double *xk = &X[k * N];
            ST_double sx = 0.0, sx2 = 0.0;
            for (i = 0; i < N; i++) {
                sx += xk[i];
                sx2 += xk[i] * xk[i];
            }
            ST_double var = (sx2 - sx * sx / N) / (N - 1);
            ST_double sd = (var > 0.0) ? sqrt(var) : 1.0;
            if (sd < 1e-3) sd = 1e-3;
            stdev_x[k] = sd;
            for (i = 0; i < N; i++)
                xk[i] /= sd;
        }
    }

    /* Initial mu: match ppmlhdfe's "simple" initialization:
     * mu = (y + mean(y, w)) / 2
     * This produces the same IRLS trajectory as ppmlhdfe. */
    {
        ST_double y_wmean = 0.0, w_sum_init = 0.0;
        for (i = 0; i < N; i++) {
            ST_double w_u = (w_user_compact != NULL) ? w_user_compact[i] : 1.0;
            y_wmean += w_u * y[i];
            w_sum_init += w_u;
        }
        y_wmean /= w_sum_init;
        for (i = 0; i < N; i++) {
            mu[i] = 0.5 * (y[i] + y_wmean);
            if (mu[i] < 1e-18) mu[i] = 1e-18;
            eta[i] = log(mu[i]);
            if (has_offset && offset_compact) {
                eta[i] -= offset_compact[i];
            }
        }
    }

    /* Initialize IRLS weights: w_irls = mu * w_user */
    for (i = 0; i < N; i++) {
        ST_double w_u = (w_user_compact != NULL) ? w_user_compact[i] : 1.0;
        irls_w[i] = mu[i] * w_u;
        if (irls_w[i] < 1e-10) irls_w[i] = 1e-10;  /* Floor to prevent zero weights */
    }

    /* Set weighted_counts for initial weights */
    ppml_update_fe_weights(g_ppml_state, irls_w, N);

    /* ================================================================
     * STEP 8b: Pre-IRLS collinearity detection
     * Partial out FE from X with initial weights, then detect collinearity
     * in X'X to identify variables collinear with the absorbed FEs.
     * ================================================================ */
    ST_int *irls_keep_idx = NULL;
    ST_int K_irls = K_x;  /* number of non-collinear X vars for IRLS */
    {
        /* Set up [z_dummy, X] for partialling - z doesn't matter, just X */
        for (i = 0; i < N; i++)
            aug_copy[i] = 0.0;  /* dummy z */
        for (k = 0; k < K_x; k++)
            memcpy(&aug_copy[(k + 1) * N], &X[k * N], N * sizeof(ST_double));

        partial_out_columns(g_ppml_state, aug_copy, N, K_aug, num_threads);

        /* Build X'X from partialled X columns */
        ST_double *xtx_pre = (ST_double *)calloc(K_x * K_x, sizeof(ST_double));
        if (xtx_pre) {
            for (i = 0; i < K_x; i++) {
                for (j = 0; j <= i; j++) {
                    ST_double d = fast_dot(&aug_copy[(i+1)*N], &aug_copy[(j+1)*N], N);
                    xtx_pre[i * K_x + j] = d;
                    xtx_pre[j * K_x + i] = d;
                }
            }
            ST_int *pre_collinear = (ST_int *)calloc(K_x, sizeof(ST_int));
            if (pre_collinear) {
                detect_collinearity(xtx_pre, K_x, pre_collinear, verbose);
                ST_int n_coll = 0;
                for (k = 0; k < K_x; k++)
                    if (pre_collinear[k]) n_coll++;
                K_irls = K_x - n_coll;
                if (n_coll > 0 && K_irls > 0) {
                    irls_keep_idx = (ST_int *)malloc(K_irls * sizeof(ST_int));
                    if (irls_keep_idx) {
                        idx = 0;
                        for (k = 0; k < K_x; k++)
                            if (!pre_collinear[k])
                                irls_keep_idx[idx++] = k;
                    }
                    if (verbose) {
                        char msg[128];
                        snprintf(msg, sizeof(msg),
                            "(cpplmhdfe: %d variable(s) collinear with FE, excluded from IRLS)\n", n_coll);
                        SF_display(msg);
                    }
                }
                free(pre_collinear);
            }
            free(xtx_pre);
        }
    }

    /* Allocate compact IRLS arrays if collinearity detected */
    ST_double *xtx_k = NULL, *xty_k = NULL, *beta_k = NULL, *inv_xx_k = NULL;
    if (irls_keep_idx && K_irls < K_x) {
        xtx_k = (ST_double *)malloc(K_irls * K_irls * sizeof(ST_double));
        xty_k = (ST_double *)malloc(K_irls * sizeof(ST_double));
        beta_k = (ST_double *)malloc(K_irls * sizeof(ST_double));
        inv_xx_k = (ST_double *)malloc(K_irls * K_irls * sizeof(ST_double));
    }

    /* ================================================================
     * STEP 9: IRLS LOOP
     * ================================================================ */
    ST_double deviance = 1e30, deviance_old;
    ST_int irls_iter;
    ST_int irls_converged = 0;

    /* Adaptive CG tolerance matching ppmlhdfe:
     * Start at 1e-4 (fast early iterations), tighten toward 1e-9
     * (precise final iterations) based on relative deviance change. */
    ST_double target_inner_tol = fmax(1e-12, fmin(1e-9, 0.1 * tolerance));
    ST_double start_inner_tol = 1e-4;
    g_ppml_state->tolerance = start_inner_tol;

    for (irls_iter = 0; irls_iter < irls_maxiter; irls_iter++) {
        deviance_old = deviance;

        /* (a) Compute working depvar: z = eta + (y - mu)/mu
         * Save in z_orig for FE reconstruction; copy into aug_copy for CG. */
        for (i = 0; i < N; i++) {
            ST_double z = eta[i] + (y[i] - mu[i]) / mu[i];
            z_orig[i] = z;
            aug_copy[i] = z;
        }

        /* (b) Copy X directly into aug_copy for CG (CG modifies in place) */
        for (k = 0; k < K_x; k++) {
            memcpy(&aug_copy[(k + 1) * N], &X[k * N], N * sizeof(ST_double));
        }

        /* (d) Partial out via CG solver using current IRLS weights */
        ST_int cg_iters = partial_out_columns(g_ppml_state, aug_copy, N, K_aug, num_threads);
        if (cg_iters < 0) cg_iters = -cg_iters;

        /* (e) Solve WLS: beta = (X~'W X~)^{-1} X~'W z~ using partialled data */
        memset(irls_beta, 0, K_x * sizeof(ST_double));

        if (irls_keep_idx && K_irls < K_x && xtx_k && xty_k && beta_k && inv_xx_k) {
            /* Compact solve: only non-collinear columns.
             * Use dd_real (double-double) precision for X'WX and X'Wz to match
             * ppmlhdfe/reghdfe's quadcross(). Small precision differences in
             * each IRLS iteration accumulate in the mu trajectory, which then
             * gets amplified through cond(X'WX) in the final VCE. */
            ST_int ki, kj;
            memset(xtx_k, 0, K_irls * K_irls * sizeof(ST_double));
            memset(xty_k, 0, K_irls * sizeof(ST_double));

            for (ki = 0; ki < K_irls; ki++) {
                ST_int ci = irls_keep_idx[ki];
                const ST_double *xi = &aug_copy[(ci + 1) * N];
                const ST_double *z_par = &aug_copy[0];
                {
                    dd_real acc = {0.0, 0.0};
                    for (i = 0; i < N; i++) {
                        ST_double wxi = irls_w[i] * xi[i];
                        dd_real prod = two_prod(wxi, z_par[i]);
                        acc = dd_add_d(acc, prod.hi);
                        acc = dd_add_d(acc, prod.lo);
                    }
                    xty_k[ki] = acc.hi + acc.lo;
                }
                for (kj = 0; kj <= ki; kj++) {
                    ST_int cj = irls_keep_idx[kj];
                    const ST_double *xj = &aug_copy[(cj + 1) * N];
                    dd_real acc = {0.0, 0.0};
                    for (i = 0; i < N; i++) {
                        ST_double wxi = irls_w[i] * xi[i];
                        dd_real prod = two_prod(wxi, xj[i]);
                        acc = dd_add_d(acc, prod.hi);
                        acc = dd_add_d(acc, prod.lo);
                    }
                    xtx_k[ki * K_irls + kj] = acc.hi + acc.lo;
                    xtx_k[kj * K_irls + ki] = acc.hi + acc.lo;
                }
            }

            memcpy(inv_xx_k, xtx_k, K_irls * K_irls * sizeof(ST_double));
            if (cholesky(inv_xx_k, K_irls) != 0) {
                break;
            }
            invert_from_cholesky(inv_xx_k, K_irls, inv_xx_k);

            for (ki = 0; ki < K_irls; ki++) {
                beta_k[ki] = 0.0;
                for (kj = 0; kj < K_irls; kj++)
                    beta_k[ki] += inv_xx_k[ki * K_irls + kj] * xty_k[kj];
            }
            /* Expand to full beta */
            for (ki = 0; ki < K_irls; ki++)
                irls_beta[irls_keep_idx[ki]] = beta_k[ki];
        } else {
            /* Full solve: all columns (no collinearity detected).
             * Use dd_real for X'WX and X'Wz to match quadcross() precision. */
            memset(xtx, 0, K_x * K_x * sizeof(ST_double));
            memset(xty, 0, K_x * sizeof(ST_double));
            for (j = 0; j < K_x; j++) {
                const ST_double *xj_ptr = &aug_copy[(j + 1) * N];
                const ST_double *z_par = &aug_copy[0];
                {
                    dd_real acc = {0.0, 0.0};
                    for (i = 0; i < N; i++) {
                        ST_double wxi = irls_w[i] * xj_ptr[i];
                        dd_real prod = two_prod(wxi, z_par[i]);
                        acc = dd_add_d(acc, prod.hi);
                        acc = dd_add_d(acc, prod.lo);
                    }
                    xty[j] = acc.hi + acc.lo;
                }
                for (k = 0; k <= j; k++) {
                    const ST_double *xk_ptr = &aug_copy[(k + 1) * N];
                    dd_real acc = {0.0, 0.0};
                    for (i = 0; i < N; i++) {
                        ST_double wxi = irls_w[i] * xj_ptr[i];
                        dd_real prod = two_prod(wxi, xk_ptr[i]);
                        acc = dd_add_d(acc, prod.hi);
                        acc = dd_add_d(acc, prod.lo);
                    }
                    xtx[j * K_x + k] = acc.hi + acc.lo;
                    xtx[k * K_x + j] = acc.hi + acc.lo;
                }
            }

            memcpy(inv_xx, xtx, K_x * K_x * sizeof(ST_double));
            if (cholesky(inv_xx, K_x) != 0) {
                for (k = 0; k < K_x; k++)
                    inv_xx[k * K_x + k] = xtx[k * K_x + k] + 1e-10;
                if (cholesky(inv_xx, K_x) != 0) {
                    break;
                }
            }
            invert_from_cholesky(inv_xx, K_x, inv_xx);

            for (i = 0; i < K_x; i++) {
                irls_beta[i] = 0.0;
                for (j = 0; j < K_x; j++)
                    irls_beta[i] += inv_xx[i * K_x + j] * xty[j];
            }
        }

        /* (f) Update eta and mu.
         * z_pred[i] = X_tilde * beta + FE_component
         * where FE_component = z_orig - z_partialled (original z minus CG result)
         */
        for (i = 0; i < N; i++) {
            ST_double xb = 0.0;
            for (k = 0; k < K_x; k++) {
                xb += aug_copy[(k + 1) * N + i] * irls_beta[k];
            }
            /* FE component: original z minus partialled z */
            ST_double fe_part = z_orig[i] - aug_copy[i];
            ST_double z_hat = xb + fe_part;
            ST_double off = (offset_compact != NULL) ? offset_compact[i] : 0.0;
            /* eta = X*beta + FE (without offset); z didn't include offset */
            eta[i] = z_hat;
            mu[i] = exp(eta[i] + off);
            /* Clamp mu to prevent overflow/underflow */
            if (mu[i] > 1e18) mu[i] = 1e18;
            if (mu[i] < 1e-18) mu[i] = 1e-18;
        }

        /* (g) Compute deviance */
        deviance = ppml_compute_deviance(y, mu, w_user_compact, N);

        /* (h) Check convergence (matching ppmlhdfe's criterion) */
        if (irls_iter > 0) {
            ST_double denom_eps = deviance < deviance_old ? deviance : deviance_old;
            if (denom_eps < 0.1) denom_eps = 0.1;
            ST_double rel_change = fabs(deviance - deviance_old) / denom_eps;
            if (verbose) {
                char msg[256];
                snprintf(msg, sizeof(msg),
                    "  iter %d: deviance=%.10e eps=%.4e\n",
                    irls_iter + 1, deviance, rel_change);
                SF_display(msg);
            }
            if (rel_change < irls_tol) {
                irls_converged = 1;
                irls_iter++;
                break;
            }

            /* Adaptive CG tolerance: tighten as IRLS converges.
             * Matches ppmlhdfe's scheme: new_tol = max(target, eps * 0.01) */
            if (rel_change < 0.1) {
                ST_double new_tol = fmax(target_inner_tol, rel_change * 0.01);
                if (new_tol < g_ppml_state->tolerance) {
                    g_ppml_state->tolerance = new_tol;
                }
            }
        } else if (verbose) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                "  iter %d: deviance=%.10e eps=.\n",
                irls_iter + 1, deviance);
            SF_display(msg);
        }

        /* (i) Update IRLS weights for next iteration */
        for (i = 0; i < N; i++) {
            ST_double w_u = (w_user_compact != NULL) ? w_user_compact[i] : 1.0;
            irls_w[i] = mu[i] * w_u;
            if (irls_w[i] < 1e-10) irls_w[i] = 1e-10;
        }
        ppml_update_fe_weights(g_ppml_state, irls_w, N);

        /* (j) Post-iteration separation check */
        {
            memset(sep_m, 0, N * sizeof(ST_int));
            ST_int new_sep = ppml_detect_separation_mu(y, mu, sep_tol, N, sep_m);
            if (new_sep > 0) {
                num_separated += new_sep;
                if (verbose) {
                    char msg[128];
                    snprintf(msg, sizeof(msg),
                        "(cpplmhdfe: %d additional separated obs detected at iter %d)\n",
                        new_sep, irls_iter + 1);
                    SF_display(msg);
                }
                /* For simplicity in v1: set mu to a small value rather than
                 * removing obs mid-loop (would require full re-compaction) */
                for (i = 0; i < N; i++) {
                    if (sep_m[i]) {
                        mu[i] = 1e-18;
                        irls_w[i] = 1e-18;
                    }
                }
                ppml_update_fe_weights(g_ppml_state, irls_w, N);
            }
        }
    }

    if (verbose) {
        char msg[256];
        if (irls_converged) {
            snprintf(msg, sizeof(msg), "(IRLS converged in %d iterations, deviance = %.8g)\n",
                     irls_iter, deviance);
        } else {
            snprintf(msg, sizeof(msg), "(IRLS did NOT converge in %d iterations, deviance = %.8g)\n",
                     irls_maxiter, deviance);
        }
        SF_display(msg);
    }

    t_irls = get_time_sec();

    /* ================================================================
     * STEP 10: Final OLS on converged data for beta, VCE
     * ================================================================ */

    /* Ensure final partialled data uses tightest CG tolerance.
     * If the adaptive scheme hasn't tightened to target yet (e.g. fast
     * convergence in few iterations), do one extra pass at target_inner_tol. */
    if (target_inner_tol < g_ppml_state->tolerance) {
        g_ppml_state->tolerance = target_inner_tol;
        memcpy(aug_copy, z_orig, N * sizeof(ST_double));
        for (k = 0; k < K_x; k++)
            memcpy(&aug_copy[(k + 1) * N], &X[k * N], N * sizeof(ST_double));
        partial_out_columns(g_ppml_state, aug_copy, N, K_aug, num_threads);
    }

    /* Collinearity detection */
    ST_int *is_collinear = (ST_int *)calloc(K_x, sizeof(ST_int));
    ST_int num_collinear = 0;
    ST_int K_keep;

    if (is_collinear) {
        ST_double *xtx_final = (ST_double *)malloc(K_x * K_x * sizeof(ST_double));
        if (xtx_final) {
            for (i = 0; i < K_x; i++) {
                for (j = 0; j < K_x; j++) {
                    xtx_final[i * K_x + j] = fast_dot(&aug_copy[(i+1)*N], &aug_copy[(j+1)*N], N);
                }
            }
            detect_collinearity(xtx_final, K_x, is_collinear, verbose);
            free(xtx_final);
        }
        for (k = 0; k < K_x; k++)
            if (is_collinear[k]) num_collinear++;
    }
    K_keep = K_x - num_collinear;

    /* Store collinearity flags */
    ctools_scal_save("__cpplmhdfe_num_collinear", (ST_double)num_collinear);
    for (k = 0; k < K_x; k++) {
        snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_collinear_%d", k + 1);
        ctools_scal_save(scalar_name, (ST_double)(is_collinear ? is_collinear[k] : 0));
    }

    /* Build non-collinear index */
    ST_int *keep_idx = (ST_int *)malloc(K_keep * sizeof(ST_int));
    if (keep_idx && is_collinear) {
        idx = 0;
        for (k = 0; k < K_x; k++) {
            if (!is_collinear[k])
                keep_idx[idx++] = k;
        }
    }

    /* Final OLS: beta and inv(X'X)
     * ppmlhdfe calls reghdfe_solve_ols with w=1 (unweighted) for the final
     * beta and VCE computation. The IRLS loop uses weighted OLS, but the
     * final solve uses unweighted OLS on the partialled data. For fweights,
     * ppmlhdfe passes true_w to reghdfe which uses weighted OLS. */
    ST_double *beta_final = (ST_double *)calloc(K_keep, sizeof(ST_double));
    ST_double *inv_xx_final = NULL;
    ST_double *V_keep = NULL;
    ST_double *xtx_keep = NULL;
    ST_double *xty_keep = NULL;
    ST_int final_use_weights = (has_weights && weight_type == 2);

    if (K_keep > 0 && keep_idx && beta_final) {
        xtx_keep = (ST_double *)calloc(K_keep * K_keep, sizeof(ST_double));
        xty_keep = (ST_double *)calloc(K_keep, sizeof(ST_double));
        inv_xx_final = (ST_double *)malloc(K_keep * K_keep * sizeof(ST_double));
        V_keep = (ST_double *)calloc(K_keep * K_keep, sizeof(ST_double));

        if (xtx_keep && xty_keep && inv_xx_final && V_keep) {

            /* Build X'X (or X'WX for fweight) using dd_real precision. */
            for (i = 0; i < K_keep; i++) {
                for (j = 0; j <= i; j++) {
                    dd_real acc = {0.0, 0.0};
                    const ST_double *xi = &aug_copy[(keep_idx[i]+1)*N];
                    const ST_double *xj = &aug_copy[(keep_idx[j]+1)*N];
                    for (idx = 0; idx < N; idx++) {
                        ST_double val = final_use_weights ? irls_w[idx] * xi[idx] : xi[idx];
                        dd_real prod = two_prod(val, xj[idx]);
                        acc = dd_add_d(acc, prod.hi);
                        acc = dd_add_d(acc, prod.lo);
                    }
                    xtx_keep[i * K_keep + j] = acc.hi + acc.lo;
                    xtx_keep[j * K_keep + i] = acc.hi + acc.lo;
                }
                {
                    dd_real acc = {0.0, 0.0};
                    const ST_double *xi = &aug_copy[(keep_idx[i]+1)*N];
                    const ST_double *z_part = aug_copy;
                    for (idx = 0; idx < N; idx++) {
                        ST_double val = final_use_weights ? irls_w[idx] * xi[idx] : xi[idx];
                        dd_real prod = two_prod(val, z_part[idx]);
                        acc = dd_add_d(acc, prod.hi);
                        acc = dd_add_d(acc, prod.lo);
                    }
                    xty_keep[i] = acc.hi + acc.lo;
                }
            }

            /* Solve with iterative refinement for maximum beta precision. */
            memcpy(inv_xx_final, xtx_keep, K_keep * K_keep * sizeof(ST_double));
            if (cholesky(inv_xx_final, K_keep) == 0) {
                invert_from_cholesky(inv_xx_final, K_keep, inv_xx_final);
                for (i = 0; i < K_keep; i++) {
                    beta_final[i] = 0.0;
                    for (j = 0; j < K_keep; j++)
                        beta_final[i] += inv_xx_final[i * K_keep + j] * xty_keep[j];
                }
                /* Iterative refinement */
                for (ST_int refine = 0; refine < 2; refine++) {
                    ST_double r_keep[64];
                    if (K_keep > 64) break;
                    for (i = 0; i < K_keep; i++) {
                        const ST_double *xi = &aug_copy[(keep_idx[i]+1)*N];
                        const ST_double *z_part = aug_copy;
                        dd_real acc = {0.0, 0.0};
                        for (idx = 0; idx < N; idx++) {
                            ST_double val = final_use_weights ? irls_w[idx] * xi[idx] : xi[idx];
                            dd_real prod = two_prod(val, z_part[idx]);
                            acc = dd_add_d(acc, prod.hi);
                            acc = dd_add_d(acc, prod.lo);
                        }
                        for (j = 0; j < K_keep; j++) {
                            dd_real prod = two_prod(xtx_keep[i * K_keep + j], beta_final[j]);
                            acc = dd_add_d(acc, -prod.hi);
                            acc = dd_add_d(acc, -prod.lo);
                        }
                        r_keep[i] = acc.hi + acc.lo;
                    }
                    for (i = 0; i < K_keep; i++) {
                        ST_double delta = 0.0;
                        for (j = 0; j < K_keep; j++)
                            delta += inv_xx_final[i * K_keep + j] * r_keep[j];
                        beta_final[i] += delta;
                    }
                }
            }
        }
    }

    /* ================================================================
     * STEP 11: VCE computation
     * Matching creghdfe/reghdfe: add weighted means back to X, include
     * constant column, extend inv_xx via block partition formula.
     * This accounts for the intercept in the sandwich VCE.
     * ================================================================ */

    /* Compute N_eff (used for DOF, corner, etc.) */
    ST_double N_eff = (ST_double)N;
    if (has_weights && weight_type == 2 && w_user_compact) {
        N_eff = 0.0;
        for (i = 0; i < N; i++) N_eff += w_user_compact[i];
    }

    /* Compute IRLS-weighted means of original X columns (before partialling).
     * These means are added back to the partialled X for VCE computation.
     * Uses IRLS weights = mean(x, HDFE.weight) matching ppmlhdfe. */
    ST_int K_with_cons = K_keep + 1;
    ST_double *means_x = (ST_double *)calloc(K_keep, sizeof(ST_double));
    if (means_x && K_keep > 0) {
        ST_double sum_w_eff = 0.0;
        for (i = 0; i < N; i++) sum_w_eff += irls_w[i];
        for (k = 0; k < K_keep; k++) {
            ST_double wmean = 0.0;
            const ST_double *xk = &X[keep_idx[k] * N];
            for (i = 0; i < N; i++) {
                wmean += irls_w[i] * xk[i];
            }
            means_x[k] = wmean / sum_w_eff;
        }
    }

    /* Extend inv_xx_final to (K_keep+1)x(K_keep+1) using block partition formula.
     * The (K+1)th element is the constant/intercept.
     * inv([A b; b' c]) with A=X'WX, b=X'W1, c=1'W1.
     * Using: side = -inv(A) * mean_x, corner = 1/c - mean_x' * side */
    ST_double *inv_xx_ext = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));
    if (inv_xx_ext && inv_xx_final && means_x && K_keep > 0) {
        ST_double *side = (ST_double *)calloc(K_keep, sizeof(ST_double));
        if (side) {
            for (j = 0; j < K_keep; j++) {
                for (i = 0; i < K_keep; i++) {
                    side[j] -= means_x[i] * inv_xx_final[i * K_keep + j];
                }
            }

            /* N_corner = 1'W1 where W = diag(effective weights used in bread)
             * Non-fweight: sum(w_norm) = N. Fweight: sum(irls_w) = sum(mu*fw). */
            ST_double N_corner = (has_weights && weight_type == 2) ?
                N_eff : (ST_double)N;
            /* Actually for fweight, the bread uses raw irls_w = mu*fw,
             * so N_corner should be sum(irls_w), not sum(fw).
             * But N_eff = sum(fw). Let me compute correctly. */
            if (has_weights && weight_type == 2) {
                N_corner = 0.0;
                for (i = 0; i < N; i++) N_corner += irls_w[i];
            }
            ST_double corner = 1.0 / N_corner;
            for (i = 0; i < K_keep; i++) {
                corner -= means_x[i] * side[i];
            }

            /* Build extended matrix */
            for (i = 0; i < K_keep; i++) {
                for (j = 0; j < K_keep; j++) {
                    inv_xx_ext[i * K_with_cons + j] = inv_xx_final[i * K_keep + j];
                }
                inv_xx_ext[i * K_with_cons + K_keep] = side[i];
                inv_xx_ext[K_keep * K_with_cons + i] = side[i];
            }
            inv_xx_ext[K_keep * K_with_cons + K_keep] = corner;

            free(side);
        }
    }

    /* Build X_eff with means added back + constant column (K_with_cons columns) */
    ST_double *X_eff = NULL;
    if (K_keep > 0) {
        X_eff = (ST_double *)malloc((size_t)N * K_with_cons * sizeof(ST_double));
        if (X_eff) {
            for (k = 0; k < K_keep; k++) {
                const ST_double *xp = &aug_copy[(keep_idx[k] + 1) * N];
                ST_double mk = means_x ? means_x[k] : 0.0;
                for (i = 0; i < N; i++) {
                    X_eff[k * N + i] = xp[i] + mk;
                }
            }
            /* Constant column */
            for (i = 0; i < N; i++) {
                X_eff[K_keep * N + i] = 1.0;
            }
        }
    }

    /* Compute VCE residual.
     * For fweight: resid = (y - mu), pass fw separately to ctools_vce.
     * Otherwise: OLS residual from partialled data: z̃ - X̃β.
     * This matches reghdfe's internal VCE residual computation. */
    ST_double *vce_resid = (ST_double *)malloc(N * sizeof(ST_double));
    if (vce_resid) {
        if (has_weights && weight_type == 2 && w_user_compact) {
            /* fweight: resid = (y - mu), fw passed separately */
            for (i = 0; i < N; i++)
                vce_resid[i] = y[i] - mu[i];
        } else if (beta_final && K_keep > 0) {
            /* OLS residual from partialled data: z̃ - X̃β */
            for (i = 0; i < N; i++) {
                ST_double xb = 0.0;
                for (k = 0; k < K_keep; k++) {
                    xb += aug_copy[(keep_idx[k] + 1) * N + i] * beta_final[k];
                }
                vce_resid[i] = aug_copy[i] - xb;
            }
        } else {
            for (i = 0; i < N; i++)
                vce_resid[i] = (y[i] - mu[i]) / mu[i];
        }
    }

    /* Cluster setup */
    ST_int df_a_nested = 0;

    if (vcetype == 2) {
        cluster_ids = (ST_int *)malloc(N * sizeof(ST_int));
        if (cluster_ids) {
            if (cluster_matches_fe >= 0) {
                num_clusters = g_ppml_state->factors[cluster_matches_fe].num_levels;
                for (idx = 0; idx < N; idx++)
                    cluster_ids[idx] = g_ppml_state->factors[cluster_matches_fe].levels[idx] - 1;

                for (g = 0; g < G; g++) {
                    ST_int is_nested;
                    if (g == cluster_matches_fe) {
                        is_nested = 1;
                    } else {
                        is_nested = ctools_fe_nested_in_cluster(
                            g_ppml_state->factors[g].levels, g_ppml_state->factors[g].num_levels,
                            cluster_ids, N);
                        if (is_nested < 0) is_nested = 0;
                    }
                    if (is_nested) df_a_nested += g_ppml_state->factors[g].num_levels;
                    snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_fe_nested_%d", g + 1);
                    ctools_scal_save(scalar_name, (ST_double)is_nested);
                }
            } else if (cluster_raw_values) {
                double *clust_compact = (double *)malloc(N * sizeof(double));
                ST_int *clust_levels = (ST_int *)malloc(N * sizeof(ST_int));
                if (clust_compact && clust_levels) {
                    idx = 0;
                    for (i = 0; i < N_orig; i++) {
                        if (mask && mask[i]) {
                            if (idx < N) clust_compact[idx] = cluster_raw_values[i];
                            idx++;
                        }
                    }
                    ST_int ncl = 0;
                    if (remap_values_sorted(clust_compact, N, clust_levels, &ncl) == 0) {
                        num_clusters = ncl;
                        for (i = 0; i < N; i++)
                            cluster_ids[i] = clust_levels[i] - 1;
                    }
                    free(clust_levels);
                }
                if (clust_compact) free(clust_compact);

                for (g = 0; g < G; g++) {
                    ST_int is_nested = ctools_fe_nested_in_cluster(
                        g_ppml_state->factors[g].levels, g_ppml_state->factors[g].num_levels,
                        cluster_ids, N);
                    if (is_nested < 0) is_nested = 0;
                    if (is_nested) df_a_nested += g_ppml_state->factors[g].num_levels;
                    snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_fe_nested_%d", g + 1);
                    ctools_scal_save(scalar_name, (ST_double)is_nested);
                }
            }
        }
    } else {
        for (g = 0; g < G; g++) {
            snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_fe_nested_%d", g + 1);
            ctools_scal_save(scalar_name, 0.0);
        }
    }

    /* Compute VCE: PPML sandwich estimator matching creghdfe/reghdfe.
     * Uses extended (K_keep+1)x(K_keep+1) system with constant column,
     * then extracts K_keep×K_keep submatrix for reporting.
     *
     * For fweight: D = extended H^{-1} (raw), resid = (y-mu), weights = fw
     * For others: D extended & normalized, resid = OLS_resid, weights = irls_w as aweight */
    ST_double *V_ext = NULL;  /* (K_with_cons)x(K_with_cons) full VCE */
    if (K_keep > 0 && vce_resid && X_eff && inv_xx_ext) {
        V_ext = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));

        if (vcetype == 0) {
            /* Unadjusted VCE: V = H^{-1} (Poisson dispersion = 1) */
            for (i = 0; i < K_keep; i++)
                for (j = 0; j < K_keep; j++)
                    V_keep[i * K_keep + j] = inv_xx_final[i * K_keep + j];

        } else if (V_ext && ((vcetype == 1) || (vcetype == 2 && cluster_ids))) {

            if (has_weights && weight_type == 2 && w_user_compact) {
                /* Fweight path: D = extended H^{-1} (raw), resid = (y-mu), fw separate */
                ctools_vce_data d;
                d.X_eff = X_eff;
                d.D = inv_xx_ext;
                d.resid = vce_resid;
                d.weights = w_user_compact;
                d.weight_type = 2;
                d.N = N;
                d.K = K_with_cons;
                d.normalize_weights = 0;

                if (vcetype == 1) {
                    ST_double dof_adj = N_eff / (N_eff - 1);
                    ctools_vce_robust(&d, dof_adj, V_ext);
                } else {
                    ST_double dof_adj = (ST_double)num_clusters / (num_clusters - 1);
                    ctools_vce_cluster(&d, cluster_ids, num_clusters, dof_adj, V_ext);
                }
            } else {
                /* Non-fweight path: unweighted, matching ppmlhdfe.
                 * ppmlhdfe calls reghdfe_solve_ols with w=1 for the final VCE,
                 * so meat uses resid^2 without IRLS weight scaling. */
                ctools_vce_data d;
                d.X_eff = X_eff;
                d.D = inv_xx_ext;
                d.resid = vce_resid;
                d.weights = NULL;
                d.weight_type = 0;  /* no weights */
                d.N = N;
                d.K = K_with_cons;
                d.normalize_weights = 0;

                if (vcetype == 1) {
                    ST_double dof_adj = (ST_double)N / (N - 1);
                    ctools_vce_robust(&d, dof_adj, V_ext);
                } else {
                    ST_double dof_adj = (ST_double)num_clusters / (num_clusters - 1);
                    ctools_vce_cluster(&d, cluster_ids, num_clusters, dof_adj, V_ext);
                }
            }

            /* Extract K_keep×K_keep submatrix (dropping constant row/col) */
            for (i = 0; i < K_keep; i++)
                for (j = 0; j < K_keep; j++)
                    V_keep[i * K_keep + j] = V_ext[i * K_with_cons + j];
        }
    }
    if (V_ext) free(V_ext);

    if (cluster_raw_values) free(cluster_raw_values);

    t_vce = get_time_sec();

    /* ================================================================
     * STEP 12: Un-standardize and store results
     * ================================================================ */

    /* Un-standardize y and mu for log-likelihood computation.
     * y_actual = y_std * stdev_y, mu_actual = mu_std * stdev_y */
    for (i = 0; i < N; i++) {
        y[i] *= stdev_y;
        mu[i] *= stdev_y;
    }

    /* Un-standardize beta: beta_actual = beta_std / stdev_x */
    if (stdev_x && keep_idx && beta_final) {
        for (k = 0; k < K_keep; k++)
            beta_final[k] /= stdev_x[keep_idx[k]];
    }

    /* Un-standardize V: V_actual[i,j] = V_std[i,j] / (stdev_x[i] * stdev_x[j]) */
    if (stdev_x && keep_idx && V_keep) {
        for (i = 0; i < K_keep; i++)
            for (j = 0; j < K_keep; j++)
                V_keep[i * K_keep + j] /= stdev_x[keep_idx[i]] * stdev_x[keep_idx[j]];
    }

    ST_double ll = ppml_compute_loglik(y, mu, w_user_compact, N);

    /* Compute log-likelihood of null model (intercept only): ll_0 */
    ST_double y_mean = 0.0, w_sum = 0.0;
    for (i = 0; i < N; i++) {
        ST_double w_u = (w_user_compact != NULL) ? w_user_compact[i] : 1.0;
        y_mean += w_u * y[i];
        w_sum += w_u;
    }
    y_mean /= w_sum;
    if (y_mean < 1e-18) y_mean = 1e-18;
    ST_double ll_0 = 0.0;
    for (i = 0; i < N; i++) {
        ST_double w_u = (w_user_compact != NULL) ? w_user_compact[i] : 1.0;
        if (y[i] > 0.0) {
            ll_0 += w_u * (y[i] * log(y_mean) - y_mean - lgamma(y[i] + 1.0));
        } else {
            ll_0 += w_u * (-y_mean);
        }
    }

    /* N reporting: for fweights, report sum of weights; otherwise report obs count */
    {
        ST_double N_report = (ST_double)N;
        if (weight_type == 2 && has_weights && w_user_compact) {
            N_report = 0.0;
            for (i = 0; i < N; i++) N_report += w_user_compact[i];
        }
        ctools_scal_save("__cpplmhdfe_N", N_report);
    }

    ctools_scal_save("__cpplmhdfe_num_singletons", (ST_double)num_singletons);
    ctools_scal_save("__cpplmhdfe_num_separated", (ST_double)num_separated);
    ctools_scal_save("__cpplmhdfe_K_keep", (ST_double)K_keep);
    ctools_scal_save("__cpplmhdfe_df_a", (ST_double)df_a);
    ctools_scal_save("__cpplmhdfe_mobility_groups", (ST_double)mobility_groups);
    ctools_scal_save("__cpplmhdfe_deviance", deviance * stdev_y);
    ctools_scal_save("__cpplmhdfe_ll", ll);
    ctools_scal_save("__cpplmhdfe_ll_0", ll_0);
    ctools_scal_save("__cpplmhdfe_irls_iterations", (ST_double)irls_iter);
    ctools_scal_save("__cpplmhdfe_irls_converged", (ST_double)irls_converged);

    for (g = 0; g < G; g++) {
        snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_num_levels_%d", g + 1);
        ctools_scal_save(scalar_name, (ST_double)factors[g].num_levels);
    }

    /* Store betas */
    for (k = 0; k < K_keep; k++) {
        snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_beta_%d", k + 1);
        ctools_scal_save(scalar_name, beta_final ? beta_final[k] : 0.0);
    }

    /* Store VCE matrix */
    if (V_keep) {
        for (i = 0; i < K_keep; i++) {
            for (j = 0; j < K_keep; j++) {
                ctools_mat_store("__cpplmhdfe_V", i + 1, j + 1, V_keep[i * K_keep + j]);
            }
        }
    }


    /* Store cluster count */
    if (vcetype == 2 && num_clusters > 0) {
        ctools_scal_save("__cpplmhdfe_N_clust", (ST_double)num_clusters);
    }

    /* Timing */
    ctools_scal_save("_cpplmhdfe_time_load", t_load - t_start);
    ctools_scal_save("_cpplmhdfe_time_remap", t_remap - t_load);
    ctools_scal_save("_cpplmhdfe_time_singleton", t_singleton - t_remap);
    ctools_scal_save("_cpplmhdfe_time_dof", t_dof - t_singleton);
    ctools_scal_save("_cpplmhdfe_time_irls", t_irls - t_dof);
    ctools_scal_save("_cpplmhdfe_time_vce", t_vce - t_irls);
    ctools_scal_save("_cpplmhdfe_time_total", t_vce - t_start);
    CTOOLS_SAVE_THREAD_INFO("_cpplmhdfe");

    /* ================================================================
     * Cleanup
     * ================================================================ */
    free(mu); free(eta);
    free(z_orig); free(aug_copy);
    free(xtx); free(xty); free(irls_beta); free(inv_xx); free(sep_m);
    if (irls_keep_idx) free(irls_keep_idx);
    if (xtx_k) free(xtx_k);
    if (xty_k) free(xty_k);
    if (beta_k) free(beta_k);
    if (inv_xx_k) free(inv_xx_k);
    if (vce_resid) free(vce_resid);
    if (beta_final) free(beta_final);
    if (inv_xx_final) free(inv_xx_final);
    if (V_keep) free(V_keep);
    if (xtx_keep) free(xtx_keep);
    if (xty_keep) free(xty_keep);
    if (is_collinear) free(is_collinear);
    if (keep_idx) free(keep_idx);
    if (X_eff) free(X_eff);
    if (means_x) free(means_x);
    if (inv_xx_ext) free(inv_xx_ext);
    if (cluster_ids) free(cluster_ids);
    if (stdev_x) free(stdev_x);
    if (w_user_compact) free(w_user_compact);
    if (offset_compact) free(offset_compact);

    for (g = 0; g < G; g++) {
        free(factors[g].levels);
        if (factors[g].counts) free(factors[g].counts);
    }
    free(factors);
    if (mask) free(mask);
    free(data);
    if (obs_map) free(obs_map);

    /* Detach irls_w from state before cleanup (it's owned by us) */
    g_ppml_state->weights = NULL;
    cleanup_ppml_state();

    return 0;
}

/* ========================================================================
 * Helper: Update FE weighted counts from IRLS weights
 * ======================================================================== */
static void ppml_update_fe_weights(HDFE_State *S, const ST_double *irls_w, ST_int N)
{
    ST_int g, i, lev;
    ST_double sum_w = 0.0;

    for (i = 0; i < N; i++)
        sum_w += irls_w[i];
    S->sum_weights = sum_w;

    for (g = 0; g < S->G; g++) {
        FE_Factor *f = &S->factors[g];
        memset(f->weighted_counts, 0, f->num_levels * sizeof(ST_double));
        for (i = 0; i < N; i++) {
            lev = f->levels[i] - 1;
            f->weighted_counts[lev] += irls_w[i];
        }
        if (f->inv_weighted_counts) {
            for (i = 0; i < f->num_levels; i++) {
                f->inv_weighted_counts[i] = (f->weighted_counts[i] > 0.0) ?
                    1.0 / f->weighted_counts[i] : 0.0;
            }
        }
    }
}

/* ========================================================================
 * Helper: Compute Poisson deviance
 * D = 2 * sum[y*log(y/mu) - (y - mu)]  (y*log(y/mu) = 0 when y = 0)
 * ======================================================================== */
static ST_double ppml_compute_deviance(const ST_double *y, const ST_double *mu,
                                        const ST_double *w_user, ST_int N)
{
    ST_int i;
    ST_double dev = 0.0;

    for (i = 0; i < N; i++) {
        ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        ST_double d;
        if (y[i] > 0.0) {
            d = y[i] * log(y[i] / mu[i]) - (y[i] - mu[i]);
        } else {
            d = mu[i];  /* 0*log(0/mu) - (0 - mu) = mu */
        }
        dev += w * d;
    }

    return 2.0 * dev;
}

/* ========================================================================
 * Helper: Compute Poisson log-likelihood
 * ll = sum[y*log(mu) - mu - log(y!)]
 * ======================================================================== */
static ST_double ppml_compute_loglik(const ST_double *y, const ST_double *mu,
                                      const ST_double *w_user, ST_int N)
{
    ST_int i;
    ST_double ll = 0.0;

    for (i = 0; i < N; i++) {
        ST_double w = (w_user != NULL) ? w_user[i] : 1.0;
        ST_double contrib;
        if (y[i] > 0.0) {
            contrib = y[i] * log(mu[i]) - mu[i] - lgamma(y[i] + 1.0);
        } else {
            contrib = -mu[i];
        }
        ll += w * contrib;
    }

    return ll;
}

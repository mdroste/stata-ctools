/*
 * creghdfe_regress.c
 *
 * Full regression command for creghdfe
 * Orchestrates HDFE init, partial out, and OLS
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_regress.h"
#include "creghdfe_utils.h"
#include "creghdfe_hdfe.h"
#include "creghdfe_solver.h"
#include "creghdfe_ols.h"
#include "creghdfe_vce.h"
#include "../ctools_hdfe_utils.h"
#include "../ctools_config.h"
#include "../ctools_types.h"  /* For ctools_data_load */
#include "../ctools_spi.h"  /* Error-checking SPI wrappers */

/*
 * FULLY COMBINED: HDFE init + Partial out + OLS in one shot
 * This reads ALL variables exactly once, eliminating data transfer overhead
 *
 * Expected varlist: depvar indepvars... fe1_var fe2_var ... feG_var [cluster_var]
 * Expected scalars:
 *   __creghdfe_K: number of data variables (depvar + indepvars)
 *   __creghdfe_G: number of FE groups
 *   __creghdfe_drop_singletons: whether to drop singletons
 *   __creghdfe_maxiter, __creghdfe_tolerance: CG solver parameters
 *   __creghdfe_verbose, __creghdfe_standardize
 *   __creghdfe_vce_type (0=unadjusted, 1=robust, 2=cluster)
 *   __creghdfe_compute_dof: whether to compute DOF/mobility groups
 *
 * Returns via scalars:
 *   __creghdfe_N: final number of observations
 *   __creghdfe_num_singletons: total singletons dropped
 *   __creghdfe_num_levels_1, __creghdfe_num_levels_2, ...: levels per FE
 *   __creghdfe_df_a: degrees of freedom absorbed
 *   __creghdfe_mobility_groups: mobility groups (if G>=2)
 *   __creghdfe_ols_N, __creghdfe_K_keep, etc.: OLS results
 */
ST_retcode do_full_regression(int argc, char *argv[])
{
    ST_int K, G, N_orig, N, in1, in2;
    ST_int k, g, i, j, idx;
    ST_double val;
    ST_int drop_singletons, verbose, compute_dof;
    ST_int num_singletons;
    ST_int *mask = NULL;  /* 1 = keep, 0 = drop */
    FactorData *factors = NULL;
    char scalar_name[64];
    double t_start, t_load, t_copy, t_remap, t_singleton, t_dof, t_partial, t_ols, t_vce;
    ST_int max_iter_singleton = 100;
    ST_int mobility_groups = 1;
    ST_int df_a = 0;
    ST_int compute_resid = 0;
    ST_int resid_var_idx = 0;

    /* Data arrays */
    ST_double *data = NULL;  /* N x K matrix (column-major) */
    ST_double *xtx = NULL, *xtx_keep = NULL, *xty_keep = NULL;
    ST_double *beta_keep = NULL, *inv_xx_keep = NULL, *V_keep = NULL;
    ST_double *data_keep = NULL;
    ST_int *is_collinear = NULL, *keep_idx = NULL;
    ST_int num_collinear, K_keep, K_x, vcetype, df_r, df_a_nested;
    ST_double tss_within, rss;
    ST_double *means = NULL, *stdevs = NULL, *tss = NULL;
    ST_int maxiter, standardize;
    ST_int num_threads;
    ST_int *cluster_ids = NULL;
    ST_int num_clusters = 0;
    perm_idx_t *obs_map = NULL;  /* Maps filtered index to 1-based Stata obs */

    (void)argc;  /* Unused */
    (void)argv;  /* Unused */

    t_start = get_time_sec();

    in1 = SF_in1();
    in2 = SF_in2();

    /* Read parameters */
    SF_scal_use("__creghdfe_K", &val); K = (ST_int)val;
    SF_scal_use("__creghdfe_G", &val); G = (ST_int)val;
    SF_scal_use("__creghdfe_drop_singletons", &val); drop_singletons = (ST_int)val;
    SF_scal_use("__creghdfe_verbose", &val); verbose = (ST_int)val;
    SF_scal_use("__creghdfe_maxiter", &val); maxiter = (ST_int)val;
    SF_scal_use("__creghdfe_tolerance", &val);
    ST_double tolerance = val;
    SF_scal_use("__creghdfe_standardize", &val); standardize = (ST_int)val;
    SF_scal_use("__creghdfe_vce_type", &val); vcetype = (ST_int)val;
    SF_scal_use("__creghdfe_df_a_nested", &val); df_a_nested = (ST_int)val;

    /* Read weight parameters */
    ST_int has_weights = 0, weight_type = 0;
    if (SF_scal_use("__creghdfe_has_weights", &val) == 0) {
        has_weights = (ST_int)val;
    }
    if (SF_scal_use("__creghdfe_weight_type", &val) == 0) {
        weight_type = (ST_int)val;
    }

    compute_dof = 0;
    if (SF_scal_use("__creghdfe_compute_dof", &val) == 0) {
        compute_dof = (ST_int)val;
    }

    /* Read DOF adjustment type: 0=all, 1=none, 2=firstpair, 3=pairwise */
    ST_int dof_adjust_type = 0;
    if (SF_scal_use("__creghdfe_dof_adjust_type", &val) == 0) {
        dof_adjust_type = (ST_int)val;
    }

    /* Read savefe flag */
    ST_int savefe = 0;
    ST_int savefe_var_idx = 0;
    if (SF_scal_use("__creghdfe_savefe", &val) == 0) {
        savefe = (ST_int)val;
    }
    if (SF_scal_use("__creghdfe_savefe_idx", &val) == 0) {
        savefe_var_idx = (ST_int)val;
    }

    /* Read groupvar flag */
    ST_int compute_groupvar = 0;
    ST_int groupvar_var_idx = 0;
    if (SF_scal_use("__creghdfe_compute_groupvar", &val) == 0) {
        compute_groupvar = (ST_int)val;
    }
    if (SF_scal_use("__creghdfe_groupvar_idx", &val) == 0) {
        groupvar_var_idx = (ST_int)val;
    }

    /* Check if we should compute and store residuals */
    if (SF_scal_use("__creghdfe_compute_resid", &val) == 0) {
        compute_resid = (ST_int)val;
    }
    if (SF_scal_use("__creghdfe_resid_var_idx", &val) == 0) {
        resid_var_idx = (ST_int)val;
    }

    /* Check if quad-precision accumulation is requested */
    ST_int use_quad = 0;
    if (SF_scal_use("__creghdfe_use_quad", &val) == 0) {
        use_quad = (ST_int)val;
    }

    if (G < 1 || G > 10) {
        SF_error("creghdfe: invalid number of FE groups (must be 1-10)\n");
        return 198;
    }
    if (K < 2) {
        SF_error("creghdfe: need at least depvar and one indepvar (K < 2)\n");
        return 198;
    }
    if (K > 1000) {
        SF_error("creghdfe: too many variables (K > 1000)\n");
        return 198;
    }

    /* Determine threads */
#ifdef _OPENMP
    num_threads = omp_get_max_threads();
    if (num_threads > K) num_threads = K;
    if (num_threads < 1) num_threads = 1;
#else
    num_threads = 1;
#endif

    /* ================================================================
     * STEP 1: PARALLEL data loading using ctools_data_load
     * This handles if/in filtering at load time, loading only filtered observations.
     * ================================================================ */

    /* Calculate total variables to load:
     * K numeric (depvar + indepvars) + G FE vars + (1 cluster if vcetype==2) + (1 weight if has_weights) */
    ST_int total_vars = K + G + (vcetype == 2 ? 1 : 0) + (has_weights ? 1 : 0);

    /* Build variable indices array (1-based for Stata) */
    int *var_indices = (int *)malloc(total_vars * sizeof(int));
    if (!var_indices) {
        SF_error("creghdfe: memory allocation failed\n");
        return 1;
    }

    /* Variables 1..K are numeric data, K+1..K+G are FE vars */
    for (i = 0; i < total_vars; i++) {
        var_indices[i] = i + 1;
    }

    /* Load all variables in parallel with if/in filtering */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);

    stata_retcode load_rc = ctools_data_load(&filtered, var_indices, total_vars, 0, 0, 0);
    free(var_indices);

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        SF_error("creghdfe: parallel data load failed\n");
        return 920;
    }

    /* Get filtered observation count and obs_map */
    N_orig = (ST_int)filtered.data.nobs;
    obs_map = filtered.obs_map;

    if (N_orig <= 0) {
        ctools_filtered_data_free(&filtered);
        SF_error("creghdfe: no observations\n");
        return 198;
    }

    /* Verify all variables were loaded successfully (defensive check) */
    for (i = 0; i < total_vars; i++) {
        if (filtered.data.vars[i].nobs > 0) {
            if (filtered.data.vars[i].type == STATA_TYPE_DOUBLE && filtered.data.vars[i].data.dbl == NULL) {
                char errmsg[128];
                snprintf(errmsg, sizeof(errmsg), "creghdfe: variable %d failed to load (NULL data pointer, nobs=%zu)\n",
                         i + 1, filtered.data.vars[i].nobs);
                SF_error(errmsg);
                ctools_filtered_data_free(&filtered);
                return 920;
            }
            if (filtered.data.vars[i].type == STATA_TYPE_STRING && filtered.data.vars[i].data.str == NULL) {
                char errmsg[128];
                snprintf(errmsg, sizeof(errmsg), "creghdfe: string variable %d failed to load (NULL data pointer, nobs=%zu)\n",
                         i + 1, filtered.data.vars[i].nobs);
                SF_error(errmsg);
                ctools_filtered_data_free(&filtered);
                return 920;
            }
        }
    }

    t_load = get_time_sec();

    /* Allocate factor structures (no hash tables needed with sort-based remapping) */
    factors = (FactorData *)calloc(G, sizeof(FactorData));
    mask = (ST_int *)ctools_safe_malloc2((size_t)N_orig, sizeof(ST_int));

    /* Allocate data matrix (column-major) */
    data = (ST_double *)ctools_safe_malloc3((size_t)N_orig, (size_t)K, sizeof(ST_double));
    means = (ST_double *)malloc(K * sizeof(ST_double));
    stdevs = (ST_double *)malloc(K * sizeof(ST_double));
    tss = (ST_double *)malloc(K * sizeof(ST_double));

    /* Allocate weight array if using weights */
    ST_double *weights = NULL;
    if (has_weights) {
        weights = (ST_double *)ctools_safe_malloc2((size_t)N_orig, sizeof(ST_double));
    }

    if (!factors || !mask || !data || !means || !stdevs || !tss ||
        (has_weights && !weights)) {
        if (factors) free(factors);
        if (mask) free(mask);
        if (data) free(data);
        if (means) free(means);
        if (stdevs) free(stdevs);
        if (tss) free(tss);
        if (weights) free(weights);
        ctools_filtered_data_free(&filtered);
        SF_error("creghdfe: memory allocation failed\n");
        return 1;
    }

    /* Initialize */
    for (i = 0; i < N_orig; i++) mask[i] = 1;
    for (k = 0; k < K; k++) means[k] = 0.0;

    /* Allocate per-factor level arrays */
    for (g = 0; g < G; g++) {
        factors[g].levels = (ST_int *)malloc(N_orig * sizeof(ST_int));
        factors[g].num_obs = N_orig;
        factors[g].counts = NULL;
        factors[g].num_levels = 0;

        if (!factors[g].levels) {
            for (i = 0; i < g; i++) {
                if (factors[i].levels) free(factors[i].levels);
            }
            free(factors); free(mask);
            free(data); free(means); free(stdevs); free(tss);
            if (weights) free(weights);
            ctools_filtered_data_free(&filtered);
            return 1;
        }
    }

    /* Copy numeric data from filtered.data to column-major data array.
     * Data is already filtered, so use direct copy. */
    for (k = 0; k < K; k++) {
        double *src = filtered.data.vars[k].data.dbl;
        double *dst = &data[k * N_orig];
        double sum = 0.0;
        for (idx = 0; idx < N_orig; idx++) {
            dst[idx] = src[idx];
            sum += src[idx];
        }
        means[k] = sum;
    }
    if (has_weights) {
        ST_int weight_var_pos = K + G + (vcetype == 2 ? 1 : 0);
        memcpy(weights, filtered.data.vars[weight_var_pos].data.dbl, N_orig * sizeof(ST_double));
    }

    /* Compute means */
    for (k = 0; k < K; k++) {
        means[k] /= N_orig;
    }

    t_copy = get_time_sec();

    /* Use sort-based remapping for FE variables (much faster than hash tables).
     * Data is already filtered, use direct access. */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (g = 0; g < G; g++) {
        ST_int num_levels_g = 0;
        if (remap_values_sorted(filtered.data.vars[K + g].data.dbl, N_orig,
                                factors[g].levels, &num_levels_g) == 0) {
            factors[g].num_levels = num_levels_g;
        }
    }

    /* Save cluster variable raw values BEFORE freeing filtered data (if clustering).
     * Also check if cluster variable matches any FE variable (common optimization). */
    double *cluster_raw_values = NULL;
    ST_int cluster_matches_fe = -1;  /* -1 = no match, 0..G-1 = which FE it matches */

    if (vcetype == 2) {
        /* Cluster variable is at position K+G in filtered.data */
        cluster_raw_values = (double *)malloc(N_orig * sizeof(double));
        if (cluster_raw_values) {
            /* Extract cluster values - direct copy since data is filtered */
            memcpy(cluster_raw_values, filtered.data.vars[K + G].data.dbl, N_orig * sizeof(double));

            /* Check if cluster values match any FE variable (common case: vce(cluster i) with absorb(i t)). */
            double *cluster_data = filtered.data.vars[K + G].data.dbl;
            for (g = 0; g < G; g++) {
                double *fe_data = filtered.data.vars[K + g].data.dbl;

                /* Quick rejection: check first, middle, and last values */
                if ((ST_int)cluster_data[0] != (ST_int)fe_data[0] ||
                    (ST_int)cluster_data[N_orig/2] != (ST_int)fe_data[N_orig/2] ||
                    (ST_int)cluster_data[N_orig-1] != (ST_int)fe_data[N_orig-1]) {
                    continue;
                }

                /* Full comparison - data is already filtered, compare directly */
                ST_int all_match = 1;
                for (idx = 0; idx < N_orig && all_match; idx++) {
                    if ((ST_int)cluster_data[idx] != (ST_int)fe_data[idx]) {
                        all_match = 0;
                    }
                }

                if (all_match) {
                    cluster_matches_fe = g;
                    break;
                }
            }
        }
    }

    /* Free filtered.data but keep obs_map for write-back operations */
    stata_data_free(&filtered.data);
    /* obs_map kept for write-back of residuals/groupvar/savefe */

    /* Verify FE remapping succeeded */
    for (g = 0; g < G; g++) {
        if (factors[g].num_levels == 0) {
            for (i = 0; i < G; i++) {
                if (factors[i].levels) free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
            }
            free(factors); free(mask);
            free(data); free(means); free(stdevs); free(tss);
            if (weights) free(weights);
            SF_error("creghdfe: FE remapping failed\n");
            return 1;
        }
    }

    /* Allocate and compute counts per FE level */
    for (g = 0; g < G; g++) {
        factors[g].counts = (ST_int *)calloc(factors[g].num_levels, sizeof(ST_int));
        if (!factors[g].counts) {
            for (i = 0; i < G; i++) {
                if (factors[i].levels) free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
            }
            free(factors); free(mask);
            free(data); free(means); free(stdevs); free(tss);
            if (weights) free(weights);
            return 1;
        }
        for (i = 0; i < N_orig; i++) {
            factors[g].counts[factors[g].levels[i] - 1]++;
        }
    }

    /* Compute weighted counts per FE level if using weights */
    ST_double *weighted_counts_orig[10] = {NULL};  /* Max 10 FE groups */
    if (has_weights) {
        for (g = 0; g < G; g++) {
            weighted_counts_orig[g] = (ST_double *)calloc(factors[g].num_levels, sizeof(ST_double));
            if (!weighted_counts_orig[g]) {
                for (i = 0; i < g; i++) {
                    if (weighted_counts_orig[i]) free(weighted_counts_orig[i]);
                }
                for (i = 0; i < G; i++) {
                    if (factors[i].levels) free(factors[i].levels);
                    if (factors[i].counts) free(factors[i].counts);
                }
                free(factors); free(mask);
                free(data); free(means); free(stdevs); free(tss);
                free(weights);
                return 1;
            }
            for (i = 0; i < N_orig; i++) {
                weighted_counts_orig[g][factors[g].levels[i] - 1] += weights[i];
            }
        }
    }

    t_remap = get_time_sec();

    /* ================================================================
     * STEP 2: Iteratively drop singletons using shared utility
     * ================================================================ */
    num_singletons = 0;
    N = N_orig;

    if (drop_singletons) {
        /* Build array of level pointers for shared utility (cast to ST_int* for API compatibility) */
        ST_int **fe_levels = (ST_int **)malloc(G * sizeof(ST_int *));
        if (fe_levels) {
            for (g = 0; g < G; g++) {
                fe_levels[g] = (ST_int *)factors[g].levels;
            }

            num_singletons = ctools_remove_singletons(fe_levels, G, N_orig, mask, max_iter_singleton, (verbose >= 1));
            free(fe_levels);

            /* Count remaining observations */
            N = 0;
            for (i = 0; i < N_orig; i++) {
                if (mask[i]) N++;
            }
        }
    }

    /* Check if all observations were dropped as singletons */
    if (N == 0) {
        ctools_scal_save("__creghdfe_N", 0.0);
        ctools_scal_save("__creghdfe_num_singletons", (ST_double)num_singletons);
        ctools_scal_save("__creghdfe_K_keep", 0.0);
        SF_error("creghdfe: all observations are singletons\n");
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        if (weights) free(weights);
        return 2001;  /* Custom error code for singleton issue */
    }

    t_singleton = get_time_sec();

    /* Store original num_levels BEFORE recount (needed for remap array sizing) */
    ST_int orig_num_levels[10];  /* Max 10 FE groups */
    for (g = 0; g < G; g++) {
        orig_num_levels[g] = factors[g].num_levels;
    }

    /* Recount levels after singleton removal.
     * Re-accumulate counts from masked observations, then count non-zeros. */
    for (g = 0; g < G; g++) {
        memset(factors[g].counts, 0, orig_num_levels[g] * sizeof(ST_int));
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) {
                ST_int level = factors[g].levels[i] - 1;
                if (level >= 0 && level < orig_num_levels[g]) {
                    factors[g].counts[level]++;
                }
            }
        }
        ST_int num_levels_after = 0;
        for (i = 0; i < orig_num_levels[g]; i++) {
            if (factors[g].counts[i] > 0) num_levels_after++;
        }
        factors[g].num_levels = num_levels_after;
    }

    /* ================================================================
     * STEP 3: Compute DOF and mobility groups
     * Respects dof_adjust_type: 0=all, 1=none, 2=firstpair, 3=pairwise
     * ================================================================ */
    df_a = 0;
    for (g = 0; g < G; g++) {
        df_a += factors[g].num_levels;
    }

    /* Array to store group assignments if groupvar requested */
    ST_int *group_assignments = NULL;
    if (compute_groupvar) {
        group_assignments = (ST_int *)calloc(N, sizeof(ST_int));
    }

    if (dof_adjust_type == 1) {
        /* dofadjustments(none): skip mobility calculation entirely */
        mobility_groups = 1;
        /* All observations in group 1 if groupvar requested */
        if (group_assignments) {
            for (idx = 0; idx < N; idx++) {
                group_assignments[idx] = 1;
            }
        }
    } else if (compute_dof && G >= 2) {
        /* Build compact level arrays for remaining observations only */
        ST_int *fe1_compact = (ST_int *)malloc(N * sizeof(ST_int));
        ST_int *fe2_compact = (ST_int *)malloc(N * sizeof(ST_int));

        if (fe1_compact && fe2_compact) {
            /* Use orig_num_levels for remap array sizing since factors[g].levels
             * still contains original IDs from 1 to orig_num_levels[g] */
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
                        fe1_compact[idx] = remap1[lev1];
                        fe2_compact[idx] = remap2[lev2];
                        idx++;
                    }
                }

                /* Count connected components with optional group tracking */
                mobility_groups = count_connected_components(
                    fe1_compact, fe2_compact, N,
                    factors[0].num_levels, factors[1].num_levels
                );

                /* If groupvar requested, compute group assignments using union-find */
                if (group_assignments) {
                    /* Use union-find to assign each observation to a group */
                    ST_int total_nodes = factors[0].num_levels + factors[1].num_levels;
                    ST_int *parent = (ST_int *)malloc(total_nodes * sizeof(ST_int));
                    if (parent) {
                        /* Initialize: each node is its own parent */
                        for (i = 0; i < total_nodes; i++) parent[i] = i;

                        /* Find with path compression */
                        #define FIND(x) ({ ST_int _x = (x); while (parent[_x] != _x) { parent[_x] = parent[parent[_x]]; _x = parent[_x]; } _x; })

                        /* Union the FE levels */
                        for (idx = 0; idx < N; idx++) {
                            ST_int node1 = fe1_compact[idx] - 1;
                            ST_int node2 = factors[0].num_levels + fe2_compact[idx] - 1;
                            ST_int root1 = FIND(node1);
                            ST_int root2 = FIND(node2);
                            if (root1 != root2) parent[root1] = root2;
                        }

                        /* Assign group IDs (1-indexed) */
                        ST_int *root_to_group = (ST_int *)calloc(total_nodes, sizeof(ST_int));
                        ST_int next_group = 1;
                        if (root_to_group) {
                            for (idx = 0; idx < N; idx++) {
                                ST_int node1 = fe1_compact[idx] - 1;
                                ST_int root = FIND(node1);
                                if (root_to_group[root] == 0) {
                                    root_to_group[root] = next_group++;
                                }
                                group_assignments[idx] = root_to_group[root];
                            }
                            free(root_to_group);
                        }

                        #undef FIND
                        free(parent);
                    }
                }

                if (mobility_groups < 0) mobility_groups = 1;

                free(remap1);
                free(remap2);
            }
            free(fe1_compact);
            free(fe2_compact);
        }

        df_a -= mobility_groups;

        /* For G >= 3, add additional mobility groups for each FE beyond the second.
         * This approximates reghdfe's more complex multi-way FE DOF calculation.
         * Only apply if dof_adjust_type is 0 (all) or 3 (pairwise) */
        if (G > 2 && (dof_adjust_type == 0 || dof_adjust_type == 3)) {
            ST_int extra_mobility = G - 2;
            df_a -= extra_mobility;
            mobility_groups += extra_mobility;
        }
    } else if (G == 1) {
        mobility_groups = 0;
        /* All observations in group 1 if groupvar requested */
        if (group_assignments) {
            for (idx = 0; idx < N; idx++) {
                group_assignments[idx] = 1;
            }
        }
    }

    t_dof = get_time_sec();

    /* ================================================================
     * STEP 4: Set up global state for CG solver (compacted data)
     * ================================================================ */
    cleanup_state();
    g_state = (HDFE_State *)calloc(1, sizeof(HDFE_State));

    if (!g_state) {
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    g_state->G = G;
    g_state->N = N;
    g_state->K = K;
    g_state->in1 = in1;
    g_state->in2 = in2;
    g_state->maxiter = maxiter;
    g_state->tolerance = tolerance;
    g_state->verbose = verbose;
    g_state->num_threads = num_threads;
    g_state->factors_initialized = 1;
    g_state->df_a = df_a;
    g_state->mobility_groups = mobility_groups;
    g_state->has_weights = has_weights;
    g_state->weight_type = weight_type;
    g_state->weights = NULL;  /* Will be set after compaction */
    g_state->sum_weights = 0.0;

    /* Allocate and copy compacted factors to global state */
    g_state->factors = (FE_Factor *)calloc(G, sizeof(FE_Factor));
    if (!g_state->factors) {
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    for (g = 0; g < G; g++) {
        g_state->factors[g].num_levels = factors[g].num_levels;
        g_state->factors[g].max_level = factors[g].num_levels - 1;  /* Remapped levels are 0 to n-1 */
        g_state->factors[g].has_intercept = 1;
        g_state->factors[g].levels = (ST_int *)malloc(N * sizeof(ST_int));
        g_state->factors[g].counts = (ST_double *)calloc(factors[g].num_levels, sizeof(ST_double));
        g_state->factors[g].means = NULL;  /* Not used in creghdfe - CG solver uses thread_fe_means */
        g_state->factors[g].weighted_counts = NULL;

        /* Allocate weighted_counts if using weights */
        if (has_weights) {
            g_state->factors[g].weighted_counts = (ST_double *)calloc(factors[g].num_levels, sizeof(ST_double));
        }

        if (!g_state->factors[g].levels || !g_state->factors[g].counts ||
            (has_weights && !g_state->factors[g].weighted_counts)) {
            cleanup_state();
            for (i = 0; i < G; i++) {
                if (factors[i].levels) free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
            }
            free(factors); free(mask);
            free(data); free(means); free(stdevs); free(tss);
            if (weights) free(weights);
            for (i = 0; i < G; i++) {
                if (weighted_counts_orig[i]) free(weighted_counts_orig[i]);
            }
            return 1;
        }

        /* Remap to contiguous levels and copy.
         * Use orig_num_levels since factors[g].levels still has original IDs. */
        ST_int *remap = (ST_int *)calloc(orig_num_levels[g] + 1, sizeof(ST_int));
        if (!remap) {
            cleanup_state();
            for (i = 0; i < G; i++) {
                if (factors[i].levels) free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
            }
            free(factors); free(mask);
            free(data); free(means); free(stdevs); free(tss);
            if (weights) free(weights);
            for (i = 0; i < G; i++) {
                if (weighted_counts_orig[i]) free(weighted_counts_orig[i]);
            }
            return 1;
        }
        ST_int next_level = 1;
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) {
                ST_int old_level = factors[g].levels[i];
                if (remap[old_level] == 0) {
                    remap[old_level] = next_level++;
                }
                g_state->factors[g].levels[idx] = remap[old_level];
                g_state->factors[g].counts[remap[old_level] - 1] += 1.0;
                /* Accumulate weighted counts if using weights */
                if (has_weights) {
                    g_state->factors[g].weighted_counts[remap[old_level] - 1] += weights[i];
                }
                idx++;
            }
        }
        free(remap);

        /* Build sorted permutation for cache-friendly scatter-gather */
        ctools_build_sorted_permutation(&g_state->factors[g], N);
    }

    /* Allocate inv_counts, inv_weighted_counts, and thread buffers */
    g_state->num_threads = num_threads;
    if (ctools_hdfe_alloc_buffers(g_state, 0, K) != 0) {
        cleanup_state();
        for (i = 0; i < G; i++) {
            if (factors[i].levels) free(factors[i].levels);
            if (factors[i].counts) free(factors[i].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        if (weights) free(weights);
        for (i = 0; i < G; i++) {
            if (weighted_counts_orig[i]) free(weighted_counts_orig[i]);
        }
        return 1;
    }

    /* ================================================================
     * STEP 5: Compact data (remove singleton rows)
     * ================================================================ */
    ST_double *data_compact = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K, sizeof(ST_double));
    ST_double *means_compact = (ST_double *)ctools_safe_malloc2((size_t)K, sizeof(ST_double));

    if (!data_compact || !means_compact) {
        if (data_compact) free(data_compact);
        if (means_compact) free(means_compact);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    /* Compact data and recompute means (weighted if using weights).
     * Build compact index first, then parallelize column copy + means. */
    for (k = 0; k < K; k++) means_compact[k] = 0.0;

    /* Build compact index map: compact_idx[j] = destination index for j-th valid obs */
    ST_int *compact_map = (ST_int *)malloc(N * sizeof(ST_int));
    if (!compact_map) {
        free(data_compact); free(means_compact);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }
    idx = 0;
    for (i = 0; i < N_orig; i++) {
        if (mask[i]) {
            compact_map[idx] = i;
            idx++;
        }
    }

    /* Compute sum_w_compact */
    ST_double sum_w_compact = 0.0;
    if (has_weights) {
        for (idx = 0; idx < N; idx++) {
            sum_w_compact += weights[compact_map[idx]];
        }
    }

    /* Parallel copy + mean accumulation per column */
    #pragma omp parallel for schedule(static) if(K > 1)
    for (k = 0; k < K; k++) {
        ST_double col_mean = 0.0;
        const ST_double *src_col = data + k * N_orig;
        ST_double *dst_col = data_compact + k * N;
        ST_int ii;
        for (ii = 0; ii < N; ii++) {
            ST_int orig_i = compact_map[ii];
            ST_double val = src_col[orig_i];
            dst_col[ii] = val;
            ST_double w = (has_weights) ? weights[orig_i] : 1.0;
            col_mean += w * val;
        }
        means_compact[k] = col_mean;
    }
    free(compact_map);

    /* Divide by sum of weights (or N if unweighted) */
    ST_double mean_divisor = (has_weights) ? sum_w_compact : (ST_double)N;
    for (k = 0; k < K; k++) {
        means_compact[k] /= mean_divisor;
    }

    /* Free original data, use compacted */
    free(data);
    data = data_compact;
    free(means);
    means = means_compact;

    /* Compact weights and assign to g_state */
    if (has_weights) {
        ST_double *weights_compact = (ST_double *)malloc(N * sizeof(ST_double));
        if (!weights_compact) {
            cleanup_state();
            for (g = 0; g < G; g++) {
                if (factors[g].levels) free(factors[g].levels);
                if (factors[g].counts) free(factors[g].counts);
            }
            free(factors); free(mask);
            free(data); free(stdevs); free(tss);
            free(weights);
            for (g = 0; g < G; g++) {
                if (weighted_counts_orig[g]) free(weighted_counts_orig[g]);
            }
            return 1;
        }

        /* Copy weights for non-singleton observations */
        ST_double sum_w = 0.0;
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) {
                weights_compact[idx] = weights[i];
                sum_w += weights[i];
                idx++;
            }
        }

        /* For aweight/pweight: normalize weights so sum(w) = N
         * This is what reghdfe does (reghdfe.mata line 3598)
         * For fweight: keep raw weights (no normalization) */
        if (weight_type == 1 || weight_type == 3) {
            ST_double scale = (ST_double)N / sum_w;
            for (idx = 0; idx < N; idx++) {
                weights_compact[idx] *= scale;
            }
            /* Also normalize the weighted_counts in factors and recompute inv_weighted_counts */
            for (g = 0; g < G; g++) {
                if (g_state->factors[g].weighted_counts) {
                    for (i = 0; i < g_state->factors[g].num_levels; i++) {
                        g_state->factors[g].weighted_counts[i] *= scale;
                    }
                    /* Recompute inv_weighted_counts after scaling */
                    if (g_state->factors[g].inv_weighted_counts) {
                        for (i = 0; i < g_state->factors[g].num_levels; i++) {
                            g_state->factors[g].inv_weighted_counts[i] =
                                (g_state->factors[g].weighted_counts[i] > 0) ?
                                1.0 / g_state->factors[g].weighted_counts[i] : 0.0;
                        }
                    }
                }
            }
            /* After normalization, sum_w should be N for aw/pw */
            /* But we keep the original sum_w for fweight calculations */
        }

        free(weights);
        g_state->weights = weights_compact;
        g_state->sum_weights = sum_w;

        /* Free original weighted counts */
        for (g = 0; g < G; g++) {
            if (weighted_counts_orig[g]) free(weighted_counts_orig[g]);
        }
    }

    /* ================================================================
     * STEP 6: Compute TSS and stdevs on compacted data (weighted if using weights)
     * Note: for aweight/pweight, weights have been normalized so sum(w) = N
     * ================================================================ */
    for (k = 0; k < K; k++) {
        ST_double *col = &data[k * N];
        ST_double mean_k = means[k];

        /* Compute TSS = sum((x - mean)^2) [or weighted variant].
         * When quad option is specified, use double-double arithmetic
         * to match Mata's quadcross() precision. Otherwise use Kahan
         * compensated summation for near-quad precision without
         * blocking FMA/SIMD vectorization. */
        ST_double ss = 0.0;
        if (has_weights && g_state->weights != NULL) {
            if (use_quad) {
                ss = dd_sum_sq_dev_weighted(col, mean_k, g_state->weights, N);
            } else {
                /* Kahan compensated summation */
                ST_double kc = 0.0;
                for (idx = 0; idx < N; idx++) {
                    ST_double d = col[idx] - mean_k;
                    ST_double sq = d * d;
                    ST_double wsq = g_state->weights[idx] * sq;
                    ST_double t = ss + wsq;
                    kc += (ss - t) + wsq;
                    ss = t;
                }
                ss += kc;
            }
            ST_double df_stdev = (weight_type == 2) ? (g_state->sum_weights - 1.0) : (ST_double)(N - 1);
            tss[k] = ss;
            stdevs[k] = sqrt(ss / df_stdev);
        } else {
            if (use_quad) {
                ss = dd_sum_sq_dev(col, mean_k, N);
            } else {
                /* Kahan compensated summation */
                ST_double kc = 0.0;
                for (idx = 0; idx < N; idx++) {
                    ST_double d = col[idx] - mean_k;
                    ST_double sq = d * d;
                    ST_double t = ss + sq;
                    kc += (ss - t) + sq;
                    ss = t;
                }
                ss += kc;
            }
            tss[k] = ss;
            stdevs[k] = sqrt(ss / (N - 1));
        }
        if (stdevs[k] < 1e-30) stdevs[k] = 1.0;
    }

    /* ================================================================
     * STEP 7: Standardize and partial out
     * ================================================================ */
    if (standardize) {
        for (k = 0; k < K; k++) {
            ST_double *col = &data[k * N];
            ST_double inv_stdev = 1.0 / stdevs[k];
            #pragma omp simd
            for (idx = 0; idx < N; idx++) {
                col[idx] *= inv_stdev;
            }
        }
    }

    /* Partial out via CG solver using shared helper */
    ST_int max_iters = partial_out_columns(g_state, data, N, K, num_threads);
    if (max_iters < 0) max_iters = -max_iters;  /* Handle failure indicator */

    /* Save iteration count to Stata scalar */
    ctools_scal_save("__creghdfe_iterations", (ST_double)max_iters);

    t_partial = get_time_sec();

    /* ================================================================
     * STEP 8: OLS with collinearity detection
     * ================================================================ */
    K_x = K - 1;  /* Excluding y */
    /* For fweight, use sum(weights) as effective N in df calculation (reghdfe.mata line 3594) */
    ST_int N_eff_df = (weight_type == 2) ? (ST_int)g_state->sum_weights : N;
    df_r = N_eff_df - df_a;

    /* Allocate for collinearity detection - cast to size_t to prevent 32-bit overflow */
    xtx = (ST_double *)malloc((size_t)K_x * K_x * sizeof(ST_double));
    is_collinear = (ST_int *)malloc(K_x * sizeof(ST_int));

    if (!xtx || !is_collinear) {
        if (xtx) free(xtx);
        if (is_collinear) free(is_collinear);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    memset(is_collinear, 0, K_x * sizeof(ST_int));

    /* Check for FE-absorbed collinearity - use very small tolerance
     * Variables that are constant within FE groups but vary across groups
     * (like time-invariant covariates) should NOT be marked collinear.
     * Only mark as collinear if the variable has essentially zero variance
     * after partialling (relative tolerance 1e-14 or absolute < 1e-30). */
    ST_double collinear_tol = 1e-14;
    num_collinear = 0;
    for (k = 0; k < K_x; k++) {
        ST_double xx_partial = fast_dot(&data[(k+1) * N], &data[(k+1) * N], N);
        ST_double xx_orig = tss[k + 1];
        /* Only mark collinear if essentially zero (both relative and absolute) */
        if (xx_orig > 0 && (xx_partial / xx_orig) <= collinear_tol && xx_partial < 1e-30) {
            is_collinear[k] = 1;
            num_collinear++;
            snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_collinear_varnum_%d", k + 1);
            ctools_scal_save(scalar_name, 1.0);
        }
    }

    /* Compute X'X */
    memset(xtx, 0, K_x * K_x * sizeof(ST_double));
    for (i = 0; i < K_x; i++) {
        for (j = 0; j < K_x; j++) {
            xtx[i * K_x + j] = fast_dot(&data[(i+1) * N], &data[(j+1) * N], N);
        }
    }

    /* Detect numerical collinearity via Cholesky */
    ST_int num_numerical_collinear = detect_collinearity(xtx, K_x, is_collinear, verbose);
    if (num_numerical_collinear < 0) {
        free(xtx); free(is_collinear);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    /* Recount collinear */
    num_collinear = 0;
    for (k = 0; k < K_x; k++) {
        if (is_collinear[k]) num_collinear++;
    }

    K_keep = K_x - num_collinear;

    /* Store collinearity flags */
    ctools_scal_save("__creghdfe_num_collinear", (ST_double)num_collinear);
    for (k = 0; k < K_x; k++) {
        snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_collinear_%d", k + 1);
        ctools_scal_save(scalar_name, (ST_double)is_collinear[k]);
    }

    if (K_keep == 0) {
        /* All X variables are collinear with FE - report as omitted (like reghdfe) */
        ctools_scal_save("__creghdfe_K_keep", 0.0);
        ctools_scal_save("__creghdfe_ols_N", (ST_double)N);
        ctools_scal_save("__creghdfe_N", (ST_double)N);
        ctools_scal_save("__creghdfe_has_cons", 1.0);
        ctools_scal_save("__creghdfe_cons", means[0]);  /* Constant = mean(y) */
        ctools_scal_save("__creghdfe_rss", tss[0]);  /* RSS = TSS when no X vars */
        ctools_scal_save("__creghdfe_tss", tss[0]);  /* Total TSS */
        ctools_scal_save("__creghdfe_tss_within", tss[0]);
        ctools_scal_save("__creghdfe_df_a", (ST_double)df_a);
        ctools_scal_save("__creghdfe_df_a_nested_computed", 0.0);
        ctools_scal_save("__creghdfe_mobility_groups", (ST_double)mobility_groups);
        ctools_scal_save("__creghdfe_num_singletons", (ST_double)num_singletons);
        /* Save number of FE levels */
        for (g = 0; g < G; g++) {
            snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_num_levels_%d", g + 1);
            ctools_scal_save(scalar_name, (ST_double)factors[g].num_levels);
        }
        free(xtx); free(is_collinear);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 0;
    }

    /* Build index of non-collinear variables */
    keep_idx = (ST_int *)malloc(K_keep * sizeof(ST_int));
    if (!keep_idx) {
        free(xtx); free(is_collinear);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    idx = 0;
    for (k = 0; k < K_x; k++) {
        if (!is_collinear[k]) {
            keep_idx[idx++] = k + 1;  /* +1 to skip y */
        }
    }

    /* Allocate for OLS - cast to size_t to prevent 32-bit overflow */
    ST_int K_with_cons = K_keep + 1;
    data_keep = (ST_double *)malloc((size_t)N * (K_keep + 1) * sizeof(ST_double));
    xtx_keep = (ST_double *)malloc(K_keep * K_keep * sizeof(ST_double));
    xty_keep = (ST_double *)malloc(K_keep * sizeof(ST_double));
    beta_keep = (ST_double *)malloc(K_with_cons * sizeof(ST_double));
    inv_xx_keep = (ST_double *)malloc(K_with_cons * K_with_cons * sizeof(ST_double));
    V_keep = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));
    ST_double *means_x = (ST_double *)malloc(K_keep * sizeof(ST_double));

    if (!data_keep || !xtx_keep || !xty_keep || !beta_keep || !inv_xx_keep || !V_keep || !means_x) {
        if (data_keep) free(data_keep);
        if (xtx_keep) free(xtx_keep);
        if (xty_keep) free(xty_keep);
        if (beta_keep) free(beta_keep);
        if (inv_xx_keep) free(inv_xx_keep);
        if (V_keep) free(V_keep);
        if (means_x) free(means_x);
        free(keep_idx); free(xtx); free(is_collinear);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    /* Copy y and non-collinear X */
    for (idx = 0; idx < N; idx++) {
        data_keep[0 * N + idx] = data[0 * N + idx];
    }
    for (k = 0; k < K_keep; k++) {
        for (idx = 0; idx < N; idx++) {
            data_keep[(k+1) * N + idx] = data[keep_idx[k] * N + idx];
        }
        if (standardize) {
            means_x[k] = means[keep_idx[k]] / stdevs[keep_idx[k]];
        } else {
            means_x[k] = means[keep_idx[k]];
        }
    }

    /* Compute X'X and X'y on non-collinear data (weighted if using weights) */
    if (has_weights) {
        compute_xtx_xty_weighted(data_keep, g_state->weights, weight_type, N, K_keep + 1, xtx_keep, xty_keep);
        /* Weighted TSS_within = sum(w_i * y_i^2)
         * Use quad-precision or Kahan summation based on option. */
        if (use_quad) {
            tss_within = dd_sum_sq_weighted(data_keep, g_state->weights, N);
        } else {
            tss_within = 0.0;
            for (idx = 0; idx < N; idx++) {
                volatile ST_double sq = data_keep[idx] * data_keep[idx];
                tss_within += g_state->weights[idx] * sq;
            }
        }
    } else {
        compute_xtx_xty(data_keep, N, K_keep + 1, xtx_keep, xty_keep);
        /* Compute tss_within: use quad-precision or Kahan summation. */
        if (use_quad) {
            tss_within = dd_sum_sq(data_keep, N);
        } else {
            tss_within = 0.0;
            for (idx = 0; idx < N; idx++) {
                volatile ST_double sq = data_keep[idx] * data_keep[idx];
                tss_within += sq;
            }
        }
    }

    t_ols = get_time_sec();

    /* Cholesky and invert */
    ST_double *inv_xx_x = (ST_double *)malloc(K_keep * K_keep * sizeof(ST_double));
    if (!inv_xx_x) {
        free(data_keep); free(xtx_keep); free(xty_keep); free(beta_keep);
        free(inv_xx_keep); free(V_keep); free(means_x); free(keep_idx);
        free(xtx); free(is_collinear);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    memcpy(inv_xx_x, xtx_keep, K_keep * K_keep * sizeof(ST_double));
    if (cholesky(inv_xx_x, K_keep) != 0) {
        SF_error("creghdfe: X'X not positive definite\n");
        free(inv_xx_x); free(data_keep); free(xtx_keep); free(xty_keep); free(beta_keep);
        free(inv_xx_keep); free(V_keep); free(means_x); free(keep_idx);
        free(xtx); free(is_collinear);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 198;
    }

    if (invert_from_cholesky(inv_xx_x, K_keep, inv_xx_x) != 0) {
        free(inv_xx_x); free(data_keep); free(xtx_keep); free(xty_keep); free(beta_keep);
        free(inv_xx_keep); free(V_keep); free(means_x); free(keep_idx);
        free(xtx); free(is_collinear);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    /* Compute beta */
    for (i = 0; i < K_keep; i++) {
        beta_keep[i] = 0.0;
        for (j = 0; j < K_keep; j++) {
            beta_keep[i] += inv_xx_x[i * K_keep + j] * xty_keep[j];
        }
    }

    /* Compute RSS */
    rss = tss_within;
    for (k = 0; k < K_keep; k++) {
        rss -= beta_keep[k] * xty_keep[k];
    }

    /* Extend inv_xx using block partition formula */
    ST_double *side = (ST_double *)malloc(K_keep * sizeof(ST_double));
    if (!side) {
        free(inv_xx_x); free(data_keep); free(xtx_keep); free(xty_keep); free(beta_keep);
        free(inv_xx_keep); free(V_keep); free(means_x); free(keep_idx);
        free(xtx); free(is_collinear);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    for (j = 0; j < K_keep; j++) {
        side[j] = 0.0;
        for (i = 0; i < K_keep; i++) {
            side[j] -= means_x[i] * inv_xx_x[i * K_keep + j];
        }
    }

    /* For fweight, use sum(weights) as effective N in corner calculation
     * This is the 1'W1 term in the block partition formula */
    ST_double N_corner = (weight_type == 2 && has_weights) ? g_state->sum_weights : (ST_double)N;
    ST_double corner = 1.0 / N_corner;
    for (i = 0; i < K_keep; i++) {
        corner -= means_x[i] * side[i];
    }

    /* Build extended inv_xx matrix */
    for (i = 0; i < K_keep; i++) {
        for (j = 0; j < K_keep; j++) {
            inv_xx_keep[i * K_with_cons + j] = inv_xx_x[i * K_keep + j];
        }
        inv_xx_keep[i * K_with_cons + K_keep] = side[i];
        inv_xx_keep[K_keep * K_with_cons + i] = side[i];
    }
    inv_xx_keep[K_keep * K_with_cons + K_keep] = corner;

    /* Compute constant */
    ST_double y_mean = means[0];
    ST_double xb_mean = 0.0;
    for (k = 0; k < K_keep; k++) {
        xb_mean += means_x[k] * beta_keep[k];
    }
    beta_keep[K_keep] = y_mean - xb_mean;

    free(side);
    free(inv_xx_x);

    /* Build data with constant for VCE - cast to size_t to prevent 32-bit overflow */
    ST_double *data_with_cons = (ST_double *)malloc((size_t)N * (K_with_cons + 1) * sizeof(ST_double));
    if (!data_with_cons) {
        free(data_keep); free(xtx_keep); free(xty_keep);
        free(beta_keep); free(inv_xx_keep); free(V_keep); free(means_x);
        free(keep_idx); free(xtx); free(is_collinear);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
        }
        free(factors); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    /* Copy y (partialled) */
    for (idx = 0; idx < N; idx++) {
        data_with_cons[0 * N + idx] = data_keep[0 * N + idx];
    }
    /* Copy X and add back means */
    for (k = 0; k < K_keep; k++) {
        for (idx = 0; idx < N; idx++) {
            data_with_cons[(k + 1) * N + idx] = data_keep[(k + 1) * N + idx] + means_x[k];
        }
    }
    /* Add constant column */
    for (idx = 0; idx < N; idx++) {
        data_with_cons[K_with_cons * N + idx] = 1.0;
    }

    /* ================================================================
     * STEP 9: Compute VCE
     * ================================================================ */

    /* Read cluster variable if clustering */
    ST_int df_a_nested_computed = 0;  /* Track FEs nested within cluster */

    if (vcetype == 2) {
        cluster_ids = (ST_int *)malloc(N * sizeof(ST_int));
        if (!cluster_ids) {
            free(data_with_cons);
            free(data_keep); free(xtx_keep); free(xty_keep);
            free(beta_keep); free(inv_xx_keep); free(V_keep); free(means_x);
            free(keep_idx); free(xtx); free(is_collinear);
            cleanup_state();
            for (g = 0; g < G; g++) {
                if (factors[g].levels) free(factors[g].levels);
                if (factors[g].counts) free(factors[g].counts);
            }
            free(factors); free(mask);
            free(data); free(means); free(stdevs); free(tss);
            return 1;
        }

        /* FAST PATH: If cluster variable matches an FE variable (e.g., vce(cluster i) with absorb(i t)),
         * we can reuse the existing FE levels as cluster IDs - no hashing or SF_vdata needed! */
        if (cluster_matches_fe >= 0) {
            /* Use existing FE levels (converted to 0-based) as cluster_ids */
            ST_int *fe_levels_for_cluster = g_state->factors[cluster_matches_fe].levels;
            num_clusters = g_state->factors[cluster_matches_fe].num_levels;

            for (idx = 0; idx < N; idx++) {
                cluster_ids[idx] = fe_levels_for_cluster[idx] - 1;  /* Convert 1-based to 0-based */
            }

            /* Check FE nesting within cluster */
            for (g = 0; g < G; g++) {
                ST_int is_nested;
                if (g == cluster_matches_fe) {
                    is_nested = 1;  /* FE is the cluster variable - trivially nested */
                } else {
                    is_nested = ctools_fe_nested_in_cluster(
                        g_state->factors[g].levels, g_state->factors[g].num_levels,
                        cluster_ids, N);
                    if (is_nested < 0) is_nested = 0;  /* Treat alloc failure as not nested */
                }
                if (is_nested) df_a_nested_computed += g_state->factors[g].num_levels;
                snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_nested_%d", g + 1);
                ctools_scal_save(scalar_name, (ST_double)is_nested);
            }
        } else {
            /* Cluster variable is different from all FE variables.
             * Use saved cluster_raw_values and remap with sort (still much faster than SF_vdata + hash). */
            if (!cluster_raw_values) {
                free(cluster_ids);
                free(data_with_cons);
                free(data_keep); free(xtx_keep); free(xty_keep);
                free(beta_keep); free(inv_xx_keep); free(V_keep); free(means_x);
                free(keep_idx); free(xtx); free(is_collinear);
                cleanup_state();
                for (g = 0; g < G; g++) {
                    if (factors[g].levels) free(factors[g].levels);
                    if (factors[g].counts) free(factors[g].counts);
                }
                free(factors); free(mask);
                free(data); free(means); free(stdevs); free(tss);
                return 1;
            }

            /* Extract cluster values for non-singleton observations and remap using sort */
            double *cluster_values_compact = (double *)malloc(N * sizeof(double));
            ST_int *cluster_levels = (ST_int *)malloc(N * sizeof(ST_int));

            if (!cluster_values_compact || !cluster_levels) {
                if (cluster_values_compact) free(cluster_values_compact);
                if (cluster_levels) free(cluster_levels);
                free(cluster_raw_values);
                free(cluster_ids);
                free(data_with_cons);
                free(data_keep); free(xtx_keep); free(xty_keep);
                free(beta_keep); free(inv_xx_keep); free(V_keep); free(means_x);
                free(keep_idx); free(xtx); free(is_collinear);
                cleanup_state();
                for (g = 0; g < G; g++) {
                    if (factors[g].levels) free(factors[g].levels);
                    if (factors[g].counts) free(factors[g].counts);
                }
                free(factors); free(mask);
                free(data); free(means); free(stdevs); free(tss);
                return 1;
            }

            /* Compact cluster values (remove singletons) */
            idx = 0;
            for (i = 0; i < N_orig; i++) {
                if (mask[i]) {
                    cluster_values_compact[idx] = cluster_raw_values[i];
                    idx++;
                }
            }

            /* Remap using sort (much faster than hash table) */
            ST_int num_cluster_levels = 0;
            if (remap_values_sorted(cluster_values_compact, N, cluster_levels, &num_cluster_levels) != 0) {
                free(cluster_values_compact);
                free(cluster_levels);
                free(cluster_raw_values);
                free(cluster_ids);
                free(data_with_cons);
                free(data_keep); free(xtx_keep); free(xty_keep);
                free(beta_keep); free(inv_xx_keep); free(V_keep); free(means_x);
                free(keep_idx); free(xtx); free(is_collinear);
                cleanup_state();
                for (g = 0; g < G; g++) {
                    if (factors[g].levels) free(factors[g].levels);
                    if (factors[g].counts) free(factors[g].counts);
                }
                free(factors); free(mask);
                free(data); free(means); free(stdevs); free(tss);
                return 1;
            }

            num_clusters = num_cluster_levels;
            for (idx = 0; idx < N; idx++) {
                cluster_ids[idx] = cluster_levels[idx] - 1;  /* Convert 1-based to 0-based */
            }

            free(cluster_values_compact);
            free(cluster_levels);

            /* Check if any FE is nested within cluster */
            for (g = 0; g < G; g++) {
                ST_int is_nested = ctools_fe_nested_in_cluster(
                    g_state->factors[g].levels, g_state->factors[g].num_levels,
                    cluster_ids, N);
                if (is_nested < 0) is_nested = 0;  /* Treat alloc failure as not nested */
                if (is_nested) df_a_nested_computed += g_state->factors[g].num_levels;
                snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_nested_%d", g + 1);
                ctools_scal_save(scalar_name, (ST_double)is_nested);
            }
        }

        /* Free saved cluster raw values */
        if (cluster_raw_values) free(cluster_raw_values);
    } else {
        /* No clustering - save 0 for all FE nested status */
        for (g = 0; g < G; g++) {
            snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_nested_%d", g + 1);
            ctools_scal_save(scalar_name, 0.0);
        }
    }

    /* Use computed df_a_nested if we found nested FEs, otherwise use passed value */
    if (df_a_nested_computed > 0) {
        df_a_nested = df_a_nested_computed;
    }

    df_r -= K_keep;

    /* Compute residuals (needed for robust/cluster VCE and optionally stored) */
    ST_double *resid = (ST_double *)malloc(N * sizeof(ST_double));
    if (resid) {
        for (idx = 0; idx < N; idx++) {
            ST_double y_hat = 0.0;
            for (k = 0; k < K_keep; k++) {
                /* Use volatile to prevent FMA, matching Mata's matrix multiply */
                volatile ST_double prod = data_keep[(k + 1) * N + idx] * beta_keep[k];
                y_hat += prod;
            }
            resid[idx] = data_keep[idx] - y_hat;
        }

        /* Recompute RSS from residuals: use quad-precision or Kahan. */
        if (use_quad) {
            if (has_weights) {
                rss = dd_sum_sq_weighted(resid, g_state->weights, N);
            } else {
                rss = dd_sum_sq(resid, N);
            }
        } else {
            rss = 0.0;
            if (has_weights) {
                for (idx = 0; idx < N; idx++) {
                    volatile ST_double sq = resid[idx] * resid[idx];
                    rss += g_state->weights[idx] * sq;
                }
            } else {
                for (idx = 0; idx < N; idx++) {
                    volatile ST_double sq = resid[idx] * resid[idx];
                    rss += sq;
                }
            }
        }
    }

    if (vcetype == 0) {
        compute_vce_unadjusted(inv_xx_keep, rss, df_r, K_with_cons, V_keep);
    } else if (vcetype == 1 || vcetype == 2) {
        if (resid) {
            /* Pass weights for weighted VCE (NULL if no weights) */
            ST_double *vce_weights = has_weights ? g_state->weights : NULL;
            /* N_eff = sum(weights) for fweight, else N */
            ST_int N_eff = (weight_type == 2) ? (ST_int)g_state->sum_weights : N;
            if (vcetype == 1) {
                compute_vce_robust(data_with_cons, resid, inv_xx_keep, vce_weights, weight_type, N, N_eff, K_with_cons, df_a, V_keep);
            } else {
                ST_int df_m_cluster = K_keep;
                compute_vce_cluster(data_with_cons, resid, inv_xx_keep, vce_weights, weight_type, cluster_ids, N, N_eff, K_with_cons, num_clusters, V_keep, df_m_cluster, df_a, df_a_nested);
            }
        }
    }

    free(data_with_cons);

    /* Destandardize if needed */
    if (standardize) {
        ST_double stdev_y = stdevs[0];

        /* Destandardize residuals first */
        if (resid) {
            for (idx = 0; idx < N; idx++) {
                resid[idx] *= stdev_y;
            }
        }

        /* Match reghdfe: multiply standardized RSS/TSS_within by stdev_y^2.
         * reghdfe computes RSS from standardized residuals via quadcross(),
         * then does sol.rss = sol.rss * stdev_y ^ 2 (reghdfe.mata:3717). */
        rss = rss * stdev_y * stdev_y;
        tss_within = tss_within * stdev_y * stdev_y;

        for (k = 0; k < K_keep; k++) {
            ST_int orig_idx = keep_idx[k];
            ST_double stdev_x = stdevs[orig_idx];
            beta_keep[k] = beta_keep[k] * stdev_y / stdev_x;
        }

        for (i = 0; i < K_with_cons; i++) {
            ST_double stdev_x_i = (i < K_keep) ? (stdevs[keep_idx[i]] / stdev_y) : (1.0 / stdev_y);
            for (j = 0; j < K_with_cons; j++) {
                ST_double stdev_x_j = (j < K_keep) ? (stdevs[keep_idx[j]] / stdev_y) : (1.0 / stdev_y);
                V_keep[i * K_with_cons + j] = V_keep[i * K_with_cons + j] / (stdev_x_i * stdev_x_j);
            }
        }

        /* Re-compute constant */
        xb_mean = 0.0;
        for (k = 0; k < K_keep; k++) {
            xb_mean += means[keep_idx[k]] * beta_keep[k];
        }
        beta_keep[K_keep] = means[0] - xb_mean;
    }

    /* Store residuals back to Stata if requested.
     * Uses obs_map to write to correct Stata observations. */
    if (compute_resid && resid_var_idx > 0 && resid) {
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) {
                SF_vstore(resid_var_idx, (ST_int)obs_map[i], resid[idx]);
                idx++;
            }
        }
    }

    /* Store groupvar (mobility group assignments) back to Stata if requested */
    if (compute_groupvar && groupvar_var_idx > 0 && group_assignments) {
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) {
                SF_vstore(groupvar_var_idx, (ST_int)obs_map[i], (ST_double)group_assignments[idx]);
                idx++;
            }
        }
    }

    /* Compute and store savefe (fixed effect estimates) if requested */
    if (savefe && savefe_var_idx > 0 && resid) {
        /* For each FE g, compute FE coefficients as the mean of residuals within each level */
        for (g = 0; g < G; g++) {
            ST_int num_levels_g = g_state->factors[g].num_levels;
            ST_double *fe_coefs = (ST_double *)calloc(num_levels_g, sizeof(ST_double));
            ST_double *fe_counts = (ST_double *)calloc(num_levels_g, sizeof(ST_double));

            if (fe_coefs && fe_counts) {
                /* Sum residuals per FE level */
                for (idx = 0; idx < N; idx++) {
                    ST_int level = g_state->factors[g].levels[idx] - 1;
                    if (level >= 0 && level < num_levels_g) {
                        fe_coefs[level] += resid[idx];
                        fe_counts[level] += 1.0;
                    }
                }

                /* Compute mean and center (subtract grand mean) */
                ST_double grand_mean = 0.0;
                ST_double grand_count = 0.0;
                for (ST_int lev = 0; lev < num_levels_g; lev++) {
                    if (fe_counts[lev] > 0) {
                        fe_coefs[lev] /= fe_counts[lev];
                        grand_mean += fe_coefs[lev] * fe_counts[lev];
                        grand_count += fe_counts[lev];
                    }
                }
                if (grand_count > 0) {
                    grand_mean /= grand_count;
                }
                for (ST_int lev = 0; lev < num_levels_g; lev++) {
                    fe_coefs[lev] -= grand_mean;
                }

                /* Store FE estimates back to Stata variable using obs_map */
                ST_int fe_var_idx = savefe_var_idx + g;
                idx = 0;
                for (i = 0; i < N_orig; i++) {
                    if (mask[i]) {
                        ST_int level = g_state->factors[g].levels[idx] - 1;
                        if (level >= 0 && level < num_levels_g) {
                            SF_vstore(fe_var_idx, (ST_int)obs_map[i], fe_coefs[level]);
                        }
                        idx++;
                    }
                }

                free(fe_coefs);
                free(fe_counts);
            }
        }
    }

    /* Free residuals and group assignments */
    if (resid) free(resid);
    if (group_assignments) free(group_assignments);

    t_vce = get_time_sec();

    /* ================================================================
     * STEP 10: Store all results
     * ================================================================ */

    /* HDFE init results */
    /* For fweight, report N as sum of weights (like reghdfe) */
    if (weight_type == 2 && has_weights) {
        ctools_scal_save("__creghdfe_N", g_state->sum_weights);
    } else {
        ctools_scal_save("__creghdfe_N", (ST_double)N);
    }
    ctools_scal_save("__creghdfe_num_singletons", (ST_double)num_singletons);
    for (g = 0; g < G; g++) {
        snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_num_levels_%d", g + 1);
        ctools_scal_save(scalar_name, (ST_double)factors[g].num_levels);
    }
    ctools_scal_save("__creghdfe_df_a", (ST_double)df_a);
    ctools_scal_save("__creghdfe_df_a_nested_computed", (ST_double)df_a_nested);
    ctools_scal_save("__creghdfe_mobility_groups", (ST_double)mobility_groups);

    /* OLS results */
    /* For fweight, report N as sum of weights (like reghdfe) */
    if (weight_type == 2 && has_weights) {
        ctools_scal_save("__creghdfe_ols_N", g_state->sum_weights);
    } else {
        ctools_scal_save("__creghdfe_ols_N", (ST_double)N);
    }
    ctools_scal_save("__creghdfe_K_keep", (ST_double)K_keep);
    ctools_scal_save("__creghdfe_has_cons", 1.0);
    ctools_scal_save("__creghdfe_rss", rss);
    ctools_scal_save("__creghdfe_tss_within", tss_within);
    ctools_scal_save("__creghdfe_tss", tss[0]);

    /* Store betas */
    for (k = 0; k < K_keep; k++) {
        snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_beta_%d", k + 1);
        ctools_scal_save(scalar_name, beta_keep[k]);
    }
    ctools_scal_save("__creghdfe_cons", beta_keep[K_keep]);

    /* Store VCE matrix directly */
    for (i = 0; i < K_with_cons; i++) {
        for (j = 0; j < K_with_cons; j++) {
            ctools_mat_store("__creghdfe_V", i + 1, j + 1, V_keep[i * K_with_cons + j]);
        }
    }

    /* Store number of clusters */
    if (vcetype == 2 && num_clusters > 0) {
        ctools_scal_save("__creghdfe_N_clust", (ST_double)num_clusters);
    }

    /* Save timing results to Stata scalars */
    ctools_scal_save("_creghdfe_time_load", t_load - t_start);
    ctools_scal_save("_creghdfe_time_copy", t_copy - t_load);
    ctools_scal_save("_creghdfe_time_remap", t_remap - t_copy);
    ctools_scal_save("_creghdfe_time_singleton", t_singleton - t_remap);
    ctools_scal_save("_creghdfe_time_dof", t_dof - t_singleton);
    ctools_scal_save("_creghdfe_time_partial", t_partial - t_dof);
    ctools_scal_save("_creghdfe_time_ols", t_ols - t_partial);
    ctools_scal_save("_creghdfe_time_vce", t_vce - t_ols);
    ctools_scal_save("_creghdfe_time_total", t_vce - t_start);
    CTOOLS_SAVE_THREAD_INFO("_creghdfe");

    /* ================================================================
     * Cleanup
     * ================================================================ */
    free(means_x);
    free(data_keep); free(xtx_keep); free(xty_keep); free(beta_keep);
    free(inv_xx_keep); free(V_keep); free(keep_idx); free(xtx); free(is_collinear);
    if (cluster_ids) free(cluster_ids);
    if (obs_map) free(obs_map);

    /* Clean up local factor data */
    for (g = 0; g < G; g++) {
        if (factors[g].levels) free(factors[g].levels);
        if (factors[g].counts) free(factors[g].counts);
    }
    free(factors); free(mask);
    free(data); free(means); free(stdevs); free(tss);

    /* Clean up global state to prevent memory leaks between calls */
    cleanup_state();

    return 0;
}

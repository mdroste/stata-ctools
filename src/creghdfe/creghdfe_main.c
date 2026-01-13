/*
 * creghdfe_main.c
 *
 * Main entry point and full_regression command
 * This is the largest module - orchestrates HDFE init, partial out, and OLS
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_main.h"
#include "creghdfe_utils.h"
#include "creghdfe_hdfe.h"
#include "creghdfe_solver.h"
#include "creghdfe_ols.h"
#include "creghdfe_vce.h"

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
    ST_int k, g, i, j, obs, idx;
    ST_double val;
    ST_int drop_singletons, verbose, compute_dof;
    ST_int num_singletons, num_singletons_iter, iter;
    ST_int *mask = NULL;  /* 1 = keep, 0 = drop */
    FactorData *factors = NULL;
    IntHashTable **hash_tables = NULL;
    char msg[256];
    char scalar_name[64];
    double t_start, t_read, t_singleton, t_dof, t_partial, t_ols, t_vce;
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
    ST_int *iter_results = NULL;
    ST_double *means = NULL, *stdevs = NULL, *tss = NULL;
    ST_int maxiter, standardize;
    ST_int num_threads;
    ST_int *cluster_ids = NULL;
    ST_int num_clusters = 0;
    ST_int t;

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

    /* Check if we should compute and store residuals */
    if (SF_scal_use("__creghdfe_compute_resid", &val) == 0) {
        compute_resid = (ST_int)val;
    }
    if (SF_scal_use("__creghdfe_resid_var_idx", &val) == 0) {
        resid_var_idx = (ST_int)val;
    }

    if (G < 1 || G > 10) {
        SF_error("creghdfe: invalid number of FE groups\n");
        return 198;
    }
    if (K < 2) {
        SF_error("creghdfe: need at least depvar and one indepvar\n");
        return 198;
    }

    /* Count valid observations */
    N_orig = 0;
    for (obs = in1; obs <= in2; obs++) {
        if (SF_ifobs(obs)) N_orig++;
    }
    if (N_orig <= 0) {
        SF_error("creghdfe: no observations\n");
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

    if (verbose >= 1) {
        sprintf(msg, "{txt}   C plugin: full_regression N=%d, K=%d, G=%d, threads=%d\n",
                N_orig, K, G, num_threads);
        SF_display(msg);
    }

    /* ================================================================
     * STEP 1: Read ALL data and FE variables in a SINGLE PASS
     * This is the key optimization - everything read once
     * ================================================================ */

    /* Allocate factor structures */
    factors = (FactorData *)calloc(G, sizeof(FactorData));
    hash_tables = (IntHashTable **)calloc(G, sizeof(IntHashTable *));
    mask = (ST_int *)malloc(N_orig * sizeof(ST_int));

    /* Allocate data matrix */
    data = (ST_double *)malloc(N_orig * K * sizeof(ST_double));
    means = (ST_double *)malloc(K * sizeof(ST_double));
    stdevs = (ST_double *)malloc(K * sizeof(ST_double));
    tss = (ST_double *)malloc(K * sizeof(ST_double));

    /* Allocate weight array if using weights */
    ST_double *weights = NULL;
    if (has_weights) {
        weights = (ST_double *)malloc(N_orig * sizeof(ST_double));
    }

    if (!factors || !hash_tables || !mask || !data || !means || !stdevs || !tss ||
        (has_weights && !weights)) {
        if (factors) free(factors);
        if (hash_tables) free(hash_tables);
        if (mask) free(mask);
        if (data) free(data);
        if (means) free(means);
        if (stdevs) free(stdevs);
        if (tss) free(tss);
        if (weights) free(weights);
        SF_error("creghdfe: memory allocation failed\n");
        return 1;
    }

    /* Initialize */
    for (i = 0; i < N_orig; i++) mask[i] = 1;
    for (k = 0; k < K; k++) means[k] = 0.0;

    /* Allocate per-factor arrays */
    for (g = 0; g < G; g++) {
        factors[g].levels = (ST_int *)malloc(N_orig * sizeof(ST_int));
        factors[g].num_obs = N_orig;
        hash_tables[g] = hash_create(N_orig);

        if (!factors[g].levels || !hash_tables[g]) {
            for (i = 0; i <= g; i++) {
                if (factors[i].levels) free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
                if (hash_tables[i]) hash_destroy(hash_tables[i]);
            }
            free(factors); free(hash_tables); free(mask);
            free(data); free(means); free(stdevs); free(tss);
            return 1;
        }
    }

    /* SINGLE PASS: Read all K data variables AND G FE variables */
    idx = 0;
    for (obs = in1; obs <= in2; obs++) {
        if (!SF_ifobs(obs)) continue;

        /* Read K data variables (depvar + indepvars) */
        for (k = 0; k < K; k++) {
            if (SF_vdata(k + 1, obs, &val)) {
                for (g = 0; g < G; g++) {
                    if (factors[g].levels) free(factors[g].levels);
                    if (factors[g].counts) free(factors[g].counts);
                    if (hash_tables[g]) hash_destroy(hash_tables[g]);
                }
                free(factors); free(hash_tables); free(mask);
                free(data); free(means); free(stdevs); free(tss);
                if (weights) free(weights);
                return 198;
            }
            data[k * N_orig + idx] = val;
            means[k] += val;
        }

        /* Read G FE variables (at positions K+1 to K+G) */
        for (g = 0; g < G; g++) {
            if (SF_vdata(K + g + 1, obs, &val)) {
                for (i = 0; i < G; i++) {
                    if (factors[i].levels) free(factors[i].levels);
                    if (factors[i].counts) free(factors[i].counts);
                    if (hash_tables[i]) hash_destroy(hash_tables[i]);
                }
                free(factors); free(hash_tables); free(mask);
                free(data); free(means); free(stdevs); free(tss);
                if (weights) free(weights);
                return 198;
            }
            factors[g].levels[idx] = hash_insert_or_get(hash_tables[g], (ST_int)val);
        }

        /* Read weight variable if using weights */
        /* Weight var position: K + G + (1 if cluster) + 1 */
        if (has_weights) {
            ST_int weight_var_idx = K + G + (vcetype == 2 ? 1 : 0) + 1;
            if (SF_vdata(weight_var_idx, obs, &val)) {
                for (i = 0; i < G; i++) {
                    if (factors[i].levels) free(factors[i].levels);
                    if (factors[i].counts) free(factors[i].counts);
                    if (hash_tables[i]) hash_destroy(hash_tables[i]);
                }
                free(factors); free(hash_tables); free(mask);
                free(data); free(means); free(stdevs); free(tss);
                free(weights);
                return 198;
            }
            weights[idx] = val;
        }

        idx++;
    }

    /* Compute means */
    for (k = 0; k < K; k++) {
        means[k] /= N_orig;
    }

    /* Finalize hash tables and allocate counts */
    for (g = 0; g < G; g++) {
        factors[g].num_levels = hash_tables[g]->size;
        factors[g].counts = (ST_int *)calloc(factors[g].num_levels, sizeof(ST_int));
        if (!factors[g].counts) {
            for (i = 0; i < G; i++) {
                if (factors[i].levels) free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
                if (hash_tables[i]) hash_destroy(hash_tables[i]);
            }
            free(factors); free(hash_tables); free(mask);
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
                    if (hash_tables[i]) hash_destroy(hash_tables[i]);
                }
                free(factors); free(hash_tables); free(mask);
                free(data); free(means); free(stdevs); free(tss);
                free(weights);
                return 1;
            }
            for (i = 0; i < N_orig; i++) {
                weighted_counts_orig[g][factors[g].levels[i] - 1] += weights[i];
            }
        }
    }

    t_read = get_time_sec();

    /* ================================================================
     * STEP 2: Iteratively drop singletons
     * ================================================================ */
    num_singletons = 0;
    N = N_orig;

    if (drop_singletons) {
        iter = 0;
        do {
            num_singletons_iter = 0;

            for (g = 0; g < G; g++) {
                ST_int found = find_singletons(&factors[g], mask);
                num_singletons_iter += found;
            }

            if (num_singletons_iter > 0) {
                num_singletons += num_singletons_iter;
                N -= num_singletons_iter;

                for (g = 0; g < G; g++) {
                    update_factor_counts(&factors[g], mask, N_orig);
                }
            }

            iter++;
        } while (num_singletons_iter > 0 && iter < max_iter_singleton && N > 1);

        if (verbose >= 1 && num_singletons > 0) {
            sprintf(msg, "{txt}(dropped %d singleton observations)\n", num_singletons);
            SF_display(msg);
        }
    }

    t_singleton = get_time_sec();

    /* Recount levels after singleton removal */
    ST_int *levels_remaining = NULL;
    for (g = 0; g < G; g++) {
        levels_remaining = (ST_int *)calloc(factors[g].num_levels, sizeof(ST_int));
        if (levels_remaining) {
            ST_int num_levels_after = 0;
            for (i = 0; i < N_orig; i++) {
                if (mask[i]) {
                    ST_int level = factors[g].levels[i] - 1;
                    if (level >= 0 && level < factors[g].num_levels && !levels_remaining[level]) {
                        levels_remaining[level] = 1;
                        num_levels_after++;
                    }
                }
            }
            factors[g].num_levels = num_levels_after;
            free(levels_remaining);
        }
    }

    /* ================================================================
     * STEP 3: Compute DOF and mobility groups
     * ================================================================ */
    df_a = 0;
    for (g = 0; g < G; g++) {
        df_a += factors[g].num_levels;
    }

    if (compute_dof && G >= 2) {
        /* Build compact level arrays for remaining observations only */
        ST_int *fe1_compact = (ST_int *)malloc(N * sizeof(ST_int));
        ST_int *fe2_compact = (ST_int *)malloc(N * sizeof(ST_int));

        if (fe1_compact && fe2_compact) {
            ST_int *remap1 = (ST_int *)calloc(hash_tables[0]->size + 1, sizeof(ST_int));
            ST_int *remap2 = (ST_int *)calloc(hash_tables[1]->size + 1, sizeof(ST_int));

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

                mobility_groups = count_connected_components(
                    fe1_compact, fe2_compact, N,
                    factors[0].num_levels, factors[1].num_levels
                );

                if (mobility_groups < 0) mobility_groups = 1;

                free(remap1);
                free(remap2);
            }
            free(fe1_compact);
            free(fe2_compact);
        }

        df_a -= mobility_groups;
    } else if (G == 1) {
        mobility_groups = 0;
    }

    t_dof = get_time_sec();

    if (verbose >= 2) {
        sprintf(msg, "{txt}   DOF: df_a=%d, mobility_groups=%d\n", df_a, mobility_groups);
        SF_display(msg);
    }

    /* ================================================================
     * STEP 4: Set up global state for CG solver (compacted data)
     * ================================================================ */
    cleanup_state();
    g_state = (HDFE_State *)calloc(1, sizeof(HDFE_State));

    if (!g_state) {
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
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
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    for (g = 0; g < G; g++) {
        g_state->factors[g].num_levels = factors[g].num_levels;
        g_state->factors[g].max_level = factors[g].num_levels - 1;  /* Remapped levels are 0 to n-1 */
        g_state->factors[g].has_intercept = 1;
        g_state->factors[g].levels = (ST_int *)malloc(N * sizeof(ST_int));
        g_state->factors[g].counts = (ST_double *)calloc(factors[g].num_levels, sizeof(ST_double));
        g_state->factors[g].means = (ST_double *)malloc(factors[g].num_levels * sizeof(ST_double));
        g_state->factors[g].weighted_counts = NULL;

        /* Allocate weighted_counts if using weights */
        if (has_weights) {
            g_state->factors[g].weighted_counts = (ST_double *)calloc(factors[g].num_levels, sizeof(ST_double));
        }

        if (!g_state->factors[g].levels || !g_state->factors[g].counts || !g_state->factors[g].means ||
            (has_weights && !g_state->factors[g].weighted_counts)) {
            cleanup_state();
            for (i = 0; i < G; i++) {
                if (factors[i].levels) free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
                if (hash_tables[i]) hash_destroy(hash_tables[i]);
            }
            free(factors); free(hash_tables); free(mask);
            free(data); free(means); free(stdevs); free(tss);
            if (weights) free(weights);
            for (i = 0; i < G; i++) {
                if (weighted_counts_orig[i]) free(weighted_counts_orig[i]);
            }
            return 1;
        }

        /* Remap to contiguous levels and copy */
        ST_int *remap = (ST_int *)calloc(hash_tables[g]->size + 1, sizeof(ST_int));
        if (remap) {
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
        }
    }

    /* Allocate thread buffers */
    g_state->thread_cg_r = (ST_double **)calloc(num_threads, sizeof(ST_double *));
    g_state->thread_cg_u = (ST_double **)calloc(num_threads, sizeof(ST_double *));
    g_state->thread_cg_v = (ST_double **)calloc(num_threads, sizeof(ST_double *));
    g_state->thread_proj = (ST_double **)calloc(num_threads, sizeof(ST_double *));
    g_state->thread_fe_means = (ST_double **)calloc(num_threads * G, sizeof(ST_double *));

    for (t = 0; t < num_threads; t++) {
        g_state->thread_cg_r[t] = (ST_double *)malloc(N * sizeof(ST_double));
        g_state->thread_cg_u[t] = (ST_double *)malloc(N * sizeof(ST_double));
        g_state->thread_cg_v[t] = (ST_double *)malloc(N * sizeof(ST_double));
        g_state->thread_proj[t] = (ST_double *)malloc(N * sizeof(ST_double));
        for (g = 0; g < G; g++) {
            g_state->thread_fe_means[t * G + g] = (ST_double *)malloc(
                g_state->factors[g].num_levels * sizeof(ST_double));
        }
    }

    /* ================================================================
     * STEP 5: Compact data (remove singleton rows)
     * ================================================================ */
    ST_double *data_compact = (ST_double *)malloc(N * K * sizeof(ST_double));
    ST_double *means_compact = (ST_double *)malloc(K * sizeof(ST_double));

    if (!data_compact || !means_compact) {
        if (data_compact) free(data_compact);
        if (means_compact) free(means_compact);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    /* Compact data and recompute means (weighted if using weights) */
    for (k = 0; k < K; k++) means_compact[k] = 0.0;
    ST_double sum_w_compact = 0.0;

    idx = 0;
    for (i = 0; i < N_orig; i++) {
        if (mask[i]) {
            ST_double w = (has_weights) ? weights[i] : 1.0;
            for (k = 0; k < K; k++) {
                data_compact[k * N + idx] = data[k * N_orig + i];
                means_compact[k] += w * data[k * N_orig + i];
            }
            sum_w_compact += w;
            idx++;
        }
    }

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
                if (hash_tables[g]) hash_destroy(hash_tables[g]);
            }
            free(factors); free(hash_tables); free(mask);
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
            /* Also normalize the weighted_counts in factors */
            for (g = 0; g < G; g++) {
                if (g_state->factors[g].weighted_counts) {
                    for (i = 0; i < g_state->factors[g].num_levels; i++) {
                        g_state->factors[g].weighted_counts[i] *= scale;
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
        ST_double ss = 0.0;
        ST_double *col = &data[k * N];
        ST_double mean_k = means[k];

        if (has_weights && g_state->weights != NULL) {
            /* Weighted TSS: sum(w * (x - mean)^2) */
            for (idx = 0; idx < N; idx++) {
                ST_double dev = col[idx] - mean_k;
                ss += g_state->weights[idx] * dev * dev;
            }
            /* For normalized weights (aw/pw), sum(w) = N, so df is still N-1 */
            /* For fweight, sum(w) = sum_weights, so df is sum_weights - 1 */
            ST_double df_stdev = (weight_type == 2) ? (g_state->sum_weights - 1.0) : (ST_double)(N - 1);
            tss[k] = ss;
            stdevs[k] = sqrt(ss / df_stdev);
        } else {
            /* Unweighted TSS */
            #pragma omp simd reduction(+:ss)
            for (idx = 0; idx < N; idx++) {
                ST_double dev = col[idx] - mean_k;
                ss += dev * dev;
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

    /* Partial out via CG solver */
    iter_results = (ST_int *)malloc(K * sizeof(ST_int));
    if (!iter_results) {
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    if (verbose >= 1) {
        SF_display("{txt}   Partialling out and solving OLS...\n");
    }

#ifdef _OPENMP
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
#endif
    for (k = 0; k < K; k++) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        iter_results[k] = cg_solve_column_threaded(g_state, &data[k * N], tid);
    }

    /* Find max iterations */
    ST_int max_iters = 0;
    for (k = 0; k < K; k++) {
        ST_int iters = iter_results[k] < 0 ? -iter_results[k] : iter_results[k];
        if (iters > max_iters) max_iters = iters;
    }

    /* Save iteration count to Stata scalar */
    SF_scal_save("__creghdfe_iterations", (ST_double)max_iters);

    if (verbose >= 1) {
        sprintf(msg, "{txt}   Converged in %d iterations\n", max_iters);
        SF_display(msg);
    }

    free(iter_results);

    t_partial = get_time_sec();

    /* ================================================================
     * STEP 8: OLS with collinearity detection
     * ================================================================ */
    K_x = K - 1;  /* Excluding y */
    /* For fweight, use sum(weights) as effective N in df calculation (reghdfe.mata line 3594) */
    ST_int N_eff_df = (weight_type == 2) ? (ST_int)g_state->sum_weights : N;
    df_r = N_eff_df - df_a;

    /* Allocate for collinearity detection */
    xtx = (ST_double *)malloc(K_x * K_x * sizeof(ST_double));
    is_collinear = (ST_int *)malloc(K_x * sizeof(ST_int));

    if (!xtx || !is_collinear) {
        if (xtx) free(xtx);
        if (is_collinear) free(is_collinear);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    memset(is_collinear, 0, K_x * sizeof(ST_int));

    /* Check for FE-absorbed collinearity */
    ST_double collinear_tol = 1e-9;
    num_collinear = 0;
    for (k = 0; k < K_x; k++) {
        ST_double xx_partial = fast_dot(&data[(k+1) * N], &data[(k+1) * N], N);
        ST_double xx_orig = tss[k + 1];
        if (xx_orig > 0 && (xx_partial / xx_orig) <= collinear_tol) {
            is_collinear[k] = 1;
            num_collinear++;
            sprintf(scalar_name, "__creghdfe_collinear_varnum_%d", k + 1);
            SF_scal_save(scalar_name, 1.0);
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
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    /* Recount collinear */
    num_collinear = 0;
    for (k = 0; k < K_x; k++) {
        if (is_collinear[k]) num_collinear++;
    }

    K_keep = K_x - num_collinear;

    if (verbose >= 1 && num_collinear > 0) {
        sprintf(msg, "{txt}   Dropped %d collinear variable%s\n",
                num_collinear, num_collinear == 1 ? "" : "s");
        SF_display(msg);
    }

    /* Store collinearity flags */
    SF_scal_save("__creghdfe_num_collinear", (ST_double)num_collinear);
    for (k = 0; k < K_x; k++) {
        sprintf(scalar_name, "__creghdfe_collinear_%d", k + 1);
        SF_scal_save(scalar_name, (ST_double)is_collinear[k]);
    }

    if (K_keep == 0) {
        SF_scal_save("__creghdfe_K_keep", 0.0);
        SF_scal_save("__creghdfe_ols_N", (ST_double)N);
        free(xtx); free(is_collinear);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
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
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    idx = 0;
    for (k = 0; k < K_x; k++) {
        if (!is_collinear[k]) {
            keep_idx[idx++] = k + 1;  /* +1 to skip y */
        }
    }

    /* Allocate for OLS */
    ST_int K_with_cons = K_keep + 1;
    data_keep = (ST_double *)malloc(N * (K_keep + 1) * sizeof(ST_double));
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
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
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
        /* Weighted TSS_within = sum(w_i * y_i^2) - with normalized weights for aw/pw */
        tss_within = 0.0;
        ST_double w_scale = 1.0;
        if (weight_type == 1 || weight_type == 3) {
            w_scale = (ST_double)N / g_state->sum_weights;
        }
        for (idx = 0; idx < N; idx++) {
            ST_double w = g_state->weights[idx] * w_scale;
            tss_within += w * data_keep[idx] * data_keep[idx];
        }
    } else {
        compute_xtx_xty(data_keep, N, K_keep + 1, xtx_keep, xty_keep);
        tss_within = fast_dot(data_keep, data_keep, N);
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
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
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
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
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
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
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
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
        free(data); free(means); free(stdevs); free(tss);
        return 1;
    }

    for (j = 0; j < K_keep; j++) {
        side[j] = 0.0;
        for (i = 0; i < K_keep; i++) {
            side[j] -= means_x[i] * inv_xx_x[i * K_keep + j];
        }
    }

    ST_double corner = 1.0 / N;
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

    /* Build data with constant for VCE */
    ST_double *data_with_cons = (ST_double *)malloc(N * (K_with_cons + 1) * sizeof(ST_double));
    if (!data_with_cons) {
        free(data_keep); free(xtx_keep); free(xty_keep);
        free(beta_keep); free(inv_xx_keep); free(V_keep); free(means_x);
        free(keep_idx); free(xtx); free(is_collinear);
        cleanup_state();
        for (g = 0; g < G; g++) {
            if (factors[g].levels) free(factors[g].levels);
            if (factors[g].counts) free(factors[g].counts);
            if (hash_tables[g]) hash_destroy(hash_tables[g]);
        }
        free(factors); free(hash_tables); free(mask);
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
                if (hash_tables[g]) hash_destroy(hash_tables[g]);
            }
            free(factors); free(hash_tables); free(mask);
            free(data); free(means); free(stdevs); free(tss);
            return 1;
        }

        /* Cluster variable is at position K+G+1 */
        ST_int cluster_var_idx = K + G + 1;
        IntHashTable *cluster_ht = hash_create(N);

        if (!cluster_ht) {
            free(cluster_ids);
            free(data_with_cons);
            free(data_keep); free(xtx_keep); free(xty_keep);
            free(beta_keep); free(inv_xx_keep); free(V_keep); free(means_x);
            free(keep_idx); free(xtx); free(is_collinear);
            cleanup_state();
            for (g = 0; g < G; g++) {
                if (factors[g].levels) free(factors[g].levels);
                if (factors[g].counts) free(factors[g].counts);
                if (hash_tables[g]) hash_destroy(hash_tables[g]);
            }
            free(factors); free(hash_tables); free(mask);
            free(data); free(means); free(stdevs); free(tss);
            return 1;
        }

        /* Read cluster IDs and check for FE variables nested within cluster */
        ST_int *cluster_raw = (ST_int *)malloc(N * sizeof(ST_int));
        if (!cluster_raw) {
            hash_destroy(cluster_ht);
            free(cluster_ids);
            free(data_with_cons);
            free(data_keep); free(xtx_keep); free(xty_keep);
            free(beta_keep); free(inv_xx_keep); free(V_keep); free(means_x);
            free(keep_idx); free(xtx); free(is_collinear);
            cleanup_state();
            for (g = 0; g < G; g++) {
                if (factors[g].levels) free(factors[g].levels);
                if (factors[g].counts) free(factors[g].counts);
                if (hash_tables[g]) hash_destroy(hash_tables[g]);
            }
            free(factors); free(hash_tables); free(mask);
            free(data); free(means); free(stdevs); free(tss);
            return 1;
        }

        idx = 0;
        i = 0;
        for (obs = in1; obs <= in2; obs++) {
            if (!SF_ifobs(obs)) continue;

            if (mask[i]) {  /* Not a singleton */
                if (SF_vdata(cluster_var_idx, obs, &val)) {
                    free(cluster_raw);
                    hash_destroy(cluster_ht);
                    free(cluster_ids);
                    free(data_with_cons);
                    free(data_keep); free(xtx_keep); free(xty_keep);
                    free(beta_keep); free(inv_xx_keep); free(V_keep); free(means_x);
                    free(keep_idx); free(xtx); free(is_collinear);
                    cleanup_state();
                    for (g = 0; g < G; g++) {
                        if (factors[g].levels) free(factors[g].levels);
                        if (factors[g].counts) free(factors[g].counts);
                        if (hash_tables[g]) hash_destroy(hash_tables[g]);
                    }
                    free(factors); free(hash_tables);
                    free(data); free(means); free(stdevs); free(tss);
                    return 198;
                }
                cluster_raw[idx] = (ST_int)val;
                cluster_ids[idx] = hash_insert_or_get(cluster_ht, (ST_int)val) - 1;  /* 0-based */
                idx++;
            }
            i++;
        }

        num_clusters = cluster_ht->size;
        hash_destroy(cluster_ht);

        /* Check if any FE variable is identical to the cluster variable */
        for (g = 0; g < G; g++) {
            ST_int *fe_raw = (ST_int *)malloc(N * sizeof(ST_int));
            if (!fe_raw) continue;

            idx = 0;
            i = 0;
            ST_int fe_var_idx = K + g + 1;
            ST_int read_ok = 1;
            for (obs = in1; obs <= in2; obs++) {
                if (!SF_ifobs(obs)) continue;
                if (mask[i]) {
                    if (SF_vdata(fe_var_idx, obs, &val)) {
                        read_ok = 0;
                        break;
                    }
                    fe_raw[idx] = (ST_int)val;
                    idx++;
                }
                i++;
            }

            if (read_ok) {
                ST_int all_match = 1;
                for (idx = 0; idx < N; idx++) {
                    if (fe_raw[idx] != cluster_raw[idx]) {
                        all_match = 0;
                        break;
                    }
                }

                if (all_match) {
                    df_a_nested_computed += g_state->factors[g].num_levels;
                    if (verbose >= 1) {
                        sprintf(msg, "{txt}   FE #%d is nested within cluster (df_a_nested += %d)\n",
                                g + 1, g_state->factors[g].num_levels);
                        SF_display(msg);
                    }
                }
            }

            free(fe_raw);
        }

        free(cluster_raw);
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
                y_hat += data_keep[(k + 1) * N + idx] * beta_keep[k];
            }
            resid[idx] = data_keep[idx] - y_hat;
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

        rss = rss * stdev_y * stdev_y;
        tss_within = tss_within * stdev_y * stdev_y;

        /* Destandardize residuals */
        if (resid) {
            for (idx = 0; idx < N; idx++) {
                resid[idx] *= stdev_y;
            }
        }

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

    /* Store residuals back to Stata if requested */
    if (compute_resid && resid_var_idx > 0 && resid) {
        ST_int nvars_plugin = SF_nvars();
        if (verbose >= 1) {
            sprintf(msg, "{txt}   Storing residuals to var index %d of %d (N=%d, N_orig=%d, in1=%d, in2=%d)\n",
                    resid_var_idx, nvars_plugin, N, N_orig, in1, in2);
            SF_display(msg);
            sprintf(msg, "{txt}   First 3 residuals: %g, %g, %g\n", resid[0], resid[1], resid[2]);
            SF_display(msg);
        }

        /* Debug: Try to read the current value of the resid variable for obs 1 */
        ST_double test_val_before;
        SF_vdata(resid_var_idx, 1, &test_val_before);
        if (verbose >= 1) {
            sprintf(msg, "{txt}   Before store: resid var[1] = %g (miss=%g)\n", test_val_before, SV_missval);
            SF_display(msg);
        }
        /* Debug: Try a direct store to obs 1 and read back */
        SF_vstore(resid_var_idx, 1, 123.456);
        ST_double test_val_direct;
        SF_vdata(resid_var_idx, 1, &test_val_direct);
        if (verbose >= 1) {
            sprintf(msg, "{txt}   After direct store to var %d: val = %g\n", resid_var_idx, test_val_direct);
            SF_display(msg);
        }
        /* Try storing to var 1 (price) to see if that works */
        ST_double price_before;
        SF_vdata(1, 1, &price_before);
        SF_vstore(1, 1, 99999.0);
        ST_double price_after;
        SF_vdata(1, 1, &price_after);
        if (verbose >= 1) {
            sprintf(msg, "{txt}   Price[1] before=%g, after store 99999=%g\n", price_before, price_after);
            SF_display(msg);
        }
        /* Restore original price */
        SF_vstore(1, 1, price_before);

        idx = 0;
        i = 0;
        ST_int stored_count = 0;
        ST_int skipped_ifobs = 0;
        ST_int skipped_mask = 0;
        for (obs = in1; obs <= in2; obs++) {
            if (!SF_ifobs(obs)) {
                skipped_ifobs++;
                continue;
            }
            if (mask[i]) {
                ST_retcode rc = SF_vstore(resid_var_idx, obs, resid[idx]);
                if (rc != 0) {
                    if (verbose >= 1 && stored_count < 3) {
                        sprintf(msg, "{err}   SF_vstore failed with rc=%d for obs=%d, idx=%d, val=%g\n",
                                rc, obs, idx, resid[idx]);
                        SF_display(msg);
                    }
                } else {
                    stored_count++;
                }
                idx++;
            } else {
                skipped_mask++;
            }
            i++;
        }

        /* Debug: Check if the value was stored */
        ST_double test_val_after;
        SF_vdata(resid_var_idx, 1, &test_val_after);
        if (verbose >= 1) {
            sprintf(msg, "{txt}   After store: resid var[1] = %g\n", test_val_after);
            SF_display(msg);
            sprintf(msg, "{txt}   Stored %d residuals (skipped: %d ifobs, %d mask)\n",
                    stored_count, skipped_ifobs, skipped_mask);
            SF_display(msg);
        }
    }

    /* Free residuals */
    if (resid) free(resid);

    t_vce = get_time_sec();

    /* ================================================================
     * STEP 10: Store all results
     * ================================================================ */

    /* HDFE init results */
    /* For fweight, report N as sum of weights (like reghdfe) */
    if (weight_type == 2 && has_weights) {
        SF_scal_save("__creghdfe_N", g_state->sum_weights);
    } else {
        SF_scal_save("__creghdfe_N", (ST_double)N);
    }
    SF_scal_save("__creghdfe_num_singletons", (ST_double)num_singletons);
    for (g = 0; g < G; g++) {
        sprintf(scalar_name, "__creghdfe_num_levels_%d", g + 1);
        SF_scal_save(scalar_name, (ST_double)factors[g].num_levels);
    }
    SF_scal_save("__creghdfe_df_a", (ST_double)df_a);
    SF_scal_save("__creghdfe_df_a_nested_computed", (ST_double)df_a_nested);
    SF_scal_save("__creghdfe_mobility_groups", (ST_double)mobility_groups);

    /* OLS results */
    /* For fweight, report N as sum of weights (like reghdfe) */
    if (weight_type == 2 && has_weights) {
        SF_scal_save("__creghdfe_ols_N", g_state->sum_weights);
    } else {
        SF_scal_save("__creghdfe_ols_N", (ST_double)N);
    }
    SF_scal_save("__creghdfe_K_keep", (ST_double)K_keep);
    SF_scal_save("__creghdfe_has_cons", 1.0);
    SF_scal_save("__creghdfe_rss", rss);
    SF_scal_save("__creghdfe_tss_within", tss_within);
    SF_scal_save("__creghdfe_tss", tss[0]);

    /* Store betas */
    for (k = 0; k < K_keep; k++) {
        sprintf(scalar_name, "__creghdfe_beta_%d", k + 1);
        SF_scal_save(scalar_name, beta_keep[k]);
    }
    SF_scal_save("__creghdfe_cons", beta_keep[K_keep]);

    /* Store VCE matrix directly */
    for (i = 0; i < K_with_cons; i++) {
        for (j = 0; j < K_with_cons; j++) {
            SF_mat_store("__creghdfe_V", i + 1, j + 1, V_keep[i * K_with_cons + j]);
        }
    }

    /* Store number of clusters */
    if (vcetype == 2 && num_clusters > 0) {
        SF_scal_save("__creghdfe_N_clust", (ST_double)num_clusters);
    }

    if (verbose >= 2) {
        sprintf(msg, "{txt}   Timing: read=%.3fs singleton=%.3fs dof=%.3fs partial=%.3fs ols=%.3fs vce=%.3fs total=%.3fs\n",
                t_read - t_start, t_singleton - t_read, t_dof - t_singleton,
                t_partial - t_dof, t_ols - t_partial, t_vce - t_ols, t_vce - t_start);
        SF_display(msg);
    }

    /* ================================================================
     * Cleanup
     * ================================================================ */
    free(means_x);
    free(data_keep); free(xtx_keep); free(xty_keep); free(beta_keep);
    free(inv_xx_keep); free(V_keep); free(keep_idx); free(xtx); free(is_collinear);
    if (cluster_ids) free(cluster_ids);

    /* Clean up local factor data */
    for (g = 0; g < G; g++) {
        if (factors[g].levels) free(factors[g].levels);
        if (factors[g].counts) free(factors[g].counts);
        if (hash_tables[g]) hash_destroy(hash_tables[g]);
    }
    free(factors); free(hash_tables); free(mask);
    free(data); free(means); free(stdevs); free(tss);

    /* Clean up global state to prevent memory leaks between calls */
    cleanup_state();

    return 0;
}

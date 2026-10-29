/*
 * creghdfe_regress.c
 *
 * Full regression command for creghdfe
 * Orchestrates HDFE init, partial out, and OLS
 * Part of the ctools Stata plugin suite
 */

/* Resolve output indices in the plugin varlist, including hidden Stata variables. */
#define SD_SAFEMODE
#include "creghdfe_regress.h"
#include "creghdfe_result.h"
#include "creghdfe_input.h"
#include "hdfe/hdfe_remap.h"
#include "hdfe/hdfe_sample.h"
#include "creghdfe_hdfe.h"
#include "hdfe/hdfe_solver.h"
#include "../ctools_ols.h"
#include "../ctools_matrix.h"
#include "creghdfe_vce.h"
#include "../ctools_hdfe_utils.h"
#include "../ctools_config.h"
#include "../ctools_types.h"  /* For ctools_data_load */
#include "../ctools_spi.h"  /* Error-checking SPI wrappers */

/* Recover the additive FE component y - Xb - residual, not the already
 * orthogonal regression residual. Backfitting starts at zero and preserves the
 * weighted zero-mean normalization for each intercept effect. */
static ST_retcode store_fixed_effects(HDFE_State *S, const perm_idx_t *obs_map,
    const ST_int *mask, ST_int N_orig, ST_int K_keep, const ST_int *keep_idx,
    const ST_double *beta, const ST_double *resid, ST_int first_var)
{
    ST_int N = S->N, G = S->G, idx = 0;
    ST_retcode rc = 0;
    ST_double *error = malloc((size_t)N * sizeof(*error));
    ST_double **alpha = calloc((size_t)G, sizeof(*alpha));
    ST_double **update = calloc((size_t)G, sizeof(*update));
    ST_double scale = 1.0;
    if (!error || !alpha || !update) { rc = 920; goto cleanup; }
    for (ST_int g = 0; g < G; g++) {
        alpha[g] = calloc((size_t)S->factors[g].num_levels, sizeof(**alpha));
        update[g] = calloc((size_t)S->factors[g].num_levels, sizeof(**update));
        if (!alpha[g] || !update[g]) { rc = 920; goto cleanup; }
    }
    for (ST_int row = 0; row < N_orig; row++) {
        if (!mask[row]) continue;
        ST_double y, x;
        rc = SF_vdata(1, (ST_int)obs_map[row], &y);
        if (rc) goto cleanup;
        error[idx] = y - beta[K_keep] - resid[idx];
        for (ST_int k = 0; k < K_keep; k++) {
            rc = SF_vdata(keep_idx[k] + 1, (ST_int)obs_map[row], &x);
            if (rc) goto cleanup;
            error[idx] -= beta[k] * x;
        }
        if (fabs(error[idx]) > scale) scale = fabs(error[idx]);
        idx++;
    }
    ST_int converged = 0;
    for (ST_int iteration = 0; iteration < S->maxiter; iteration++) {
        ST_double largest_update = 0.0;
        for (ST_int g = 0; g < G; g++) {
            FE_Factor *f = &S->factors[g];
            memset(update[g], 0, (size_t)f->num_levels * sizeof(**update));
            for (ST_int i = 0; i < N; i++) {
                ST_double weight = S->weights ? S->weights[i] : 1.0;
                update[g][f->levels[i] - 1] += weight * error[i];
            }
            for (ST_int level = 0; level < f->num_levels; level++) {
                ST_double count = S->weights ? f->weighted_counts[level] : f->counts[level];
                update[g][level] = count > 0 ? update[g][level] / count : 0;
                alpha[g][level] += update[g][level];
                if (fabs(update[g][level]) > largest_update) largest_update = fabs(update[g][level]);
            }
            for (ST_int i = 0; i < N; i++) error[i] -= update[g][f->levels[i] - 1];
        }
        if (largest_update <= 1e-14 * scale) { converged = 1; break; }
    }
    if (!converged) { rc = 430; goto cleanup; }
    for (ST_int g = 0; g < G; g++) {
        idx = 0;
        for (ST_int row = 0; row < N_orig; row++) {
            if (!mask[row]) continue;
            rc = SF_vstore(first_var + g, (ST_int)obs_map[row], alpha[g][S->factors[g].levels[idx++] - 1]);
            if (rc) goto cleanup;
        }
    }
cleanup:
    for (ST_int g = 0; g < G; g++) {
        if (alpha) free(alpha[g]);
        if (update) free(update[g]);
    }
    free(alpha); free(update); free(error);
    return rc;
}

/* A mobility group's smallest value of the first absorbed effect */
typedef struct {
    ST_double value;
    ST_int root;
} creghdfe_group_min;

static int creghdfe_group_min_cmp(const void *a, const void *b)
{
    ST_double va = ((const creghdfe_group_min *)a)->value;
    ST_double vb = ((const creghdfe_group_min *)b)->value;
    return (va > vb) - (va < vb);
}

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
 *   __creghdfe_fe_redundant_1, ...: redundant coefficients per FE
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
    ST_int mobility_groups = 1;
    ST_double df_a = 0;
    ST_int compute_resid = 0;
    ST_int resid_var_idx = 0;
    ST_retcode output_rc = 0;

    /* Data arrays */
    ST_double *data = NULL;  /* N x K matrix (column-major) */
    ST_double *xtx = NULL, *xtx_keep = NULL, *xty_keep = NULL;
    ST_double *beta_keep = NULL, *inv_xx_keep = NULL, *V_keep = NULL;
    ST_double *data_keep = NULL;
    ST_int *is_collinear = NULL, *keep_idx = NULL;
    ST_int num_collinear, K_keep, K_x, vcetype;
    ST_double df_a_nested;
    ST_double df_r;
    ST_double tss_within, rss;
    ST_double *means = NULL, *stdevs = NULL, *tss = NULL;
    ST_int maxiter, standardize;
    ST_int num_threads;
    ST_int *cluster_ids = NULL;
    ST_int num_clusters = 0;
    perm_idx_t *obs_map = NULL;  /* Maps filtered index to 1-based Stata obs */

    ST_int *root_to_group = NULL;
    ST_double *weights = NULL;
    double *cluster_raw_values = NULL;
    ST_int *group_assignments = NULL;
    ST_double *data_compact = NULL;
    ST_double *means_compact = NULL;
    ST_int *compact_map = NULL;
    ST_double *weights_compact = NULL;
    ST_double *means_x = NULL;
    ST_double *inv_xx_x = NULL;
    ST_double *data_with_cons = NULL;
    double *cluster_values_compact = NULL;
    ST_double *resid = NULL;

    (void)argc;  /* Unused */
    (void)argv;  /* Unused */

    t_start = get_time_sec();

    in1 = SF_in1();
    in2 = SF_in2();

    creghdfe_input input;
    ST_retcode input_rc = creghdfe_read_input(&input);
    if (input_rc) return input_rc;
    K = input.K;
    G = input.G;
    drop_singletons = input.drop_singletons;
    verbose = input.verbose;
    maxiter = input.maxiter;
    standardize = input.standardize;
    vcetype = input.vcetype;
    ST_int has_weights = input.has_weights;
    ST_int weight_type = input.weight_type;
    compute_dof = input.compute_dof;
    ST_int dof_adjust_type = input.dof_adjust_type;
    ST_int dof_clusters = input.dof_clusters;
    ST_int savefe = input.savefe;
    ST_int savefe_var_idx = input.savefe_var_idx;
    ST_int compute_groupvar = input.compute_groupvar;
    ST_int groupvar_var_idx = input.groupvar_var_idx;
    compute_resid = input.compute_resid;
    resid_var_idx = input.resid_var_idx;
    ST_int use_quad = input.use_quad;
    ST_double tolerance = input.tolerance;
    df_a_nested = input.df_a_nested;

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

    /* Every variable position below follows the wrapper's layout: K data
     * vars, G FE vars, [cluster], [weight], [resid], [groupvar], [G savefe
     * vars], final-sample flag, [sample marker]. Refuse any other length
     * rather than read or write the wrong variables. */
    {
        ST_int expected_nvars = K + G + (vcetype == 2) + (has_weights != 0) +
                                (compute_resid != 0) + (compute_groupvar != 0) +
                                (savefe ? G : 0) + 1 + input.touse;
        if (SF_nvars() != expected_nvars) {
            char msg[160];
            snprintf(msg, sizeof(msg),
                     "creghdfe: plugin received %d variables, expected %d\n",
                     (int)SF_nvars(), (int)expected_nvars);
            SF_error(msg);
            return 198;
        }
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
        return 920;
    }

    /* Variables 1..K are numeric data, K+1..K+G are FE vars */
    for (i = 0; i < total_vars; i++) {
        var_indices[i] = i + 1;
    }

    /* Load all variables in parallel with if/in filtering */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);

    int load_flags = CTOOLS_LOAD_NO_SORT_ORDER | (input.touse ? CTOOLS_LOAD_MARKER : CTOOLS_LOAD_CHECK_IF);
    stata_retcode load_rc = ctools_data_load(&filtered, var_indices, total_vars, 0, 0, load_flags);
    free(var_indices);

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        SF_error("creghdfe: parallel data load failed\n");
        return ctools_stata_rc(load_rc);
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

    filtered.obs_map = NULL;  /* obs_map is now owned by this command. */
    /* Allocate factor structures (no hash tables needed with sort-based remapping) */
    factors = (FactorData *)calloc(G, sizeof(FactorData));
    mask = (ST_int *)ctools_safe_malloc2((size_t)N_orig, sizeof(ST_int));

    /* The data matrix (column-major) is built in STEP 5 from the loaded
     * columns, for the retained rows only */
    means = (ST_double *)malloc(K * sizeof(ST_double));
    stdevs = (ST_double *)malloc(K * sizeof(ST_double));
    tss = (ST_double *)malloc(K * sizeof(ST_double));

    /* Allocate weight array if using weights */

    if (has_weights) {
        weights = (ST_double *)ctools_safe_malloc2((size_t)N_orig, sizeof(ST_double));
    }

    if (!factors || !mask || !means || !stdevs || !tss ||
        (has_weights && !weights)) {

        SF_error("creghdfe: memory allocation failed\n");
        output_rc = 920;
        goto cleanup;
    }

    /* Initialize */
    for (i = 0; i < N_orig; i++) mask[i] = 1;
    if (has_weights) {
        ST_int weight_var_pos = K + G + (vcetype == 2 ? 1 : 0);
        memcpy(weights, filtered.data.vars[weight_var_pos].data.dbl, N_orig * sizeof(ST_double));
    }

    t_copy = get_time_sec();

    /* FE variables to 1-based levels with counts, one factor per thread */
    ST_retcode remap_rc;
    {
        const ST_double *fe_values[10];
        for (g = 0; g < G; g++) fe_values[g] = filtered.data.vars[K + g].data.dbl;
        remap_rc = hdfe_load_factors(factors, G, fe_values, N_orig);
    }

    /* Save cluster variable raw values BEFORE freeing filtered data (if clustering).
     * Also check if cluster variable matches any FE variable (common optimization). */

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
                if (cluster_data[0] != fe_data[0] ||
                    cluster_data[N_orig/2] != fe_data[N_orig/2] ||
                    cluster_data[N_orig-1] != fe_data[N_orig-1]) {
                    continue;
                }

                /* Full comparison - data is already filtered, compare directly */
                ST_int all_match = 1;
                for (idx = 0; idx < N_orig && all_match; idx++) {
                    if (cluster_data[idx] != fe_data[idx]) {
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

    /* filtered.data stays loaded until STEP 5 has compacted its columns;
     * obs_map is kept for write-back of residuals/groupvar/savefe */

    if (remap_rc) {
        SF_error("creghdfe: FE remapping failed\n");
        output_rc = remap_rc;
        goto cleanup;
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
        num_singletons = -1;
        if (fe_levels) {
            for (g = 0; g < G; g++) {
                fe_levels[g] = (ST_int *)factors[g].levels;
            }

            /* As in reghdfe, fweights drop only levels whose total weight is 1;
             * aweights and pweights use the unweighted rule. */
            num_singletons = ctools_remove_singletons(fe_levels, G, N_orig, mask,
                (has_weights && weight_type == 2) ? weights : NULL, (verbose >= 1));
            free(fe_levels);
        }

        /* Count remaining observations (a failure leaves N = 0, handled below) */
        N = 0;
        if (num_singletons >= 0) {
            for (i = 0; i < N_orig; i++) {
                if (mask[i]) N++;
            }
        }
    }

    /* Check if all observations were dropped as singletons */
    if (N == 0) {
        if (num_singletons >= 0) {
            if ((output_rc = ctools_scal_save("__creghdfe_N", 0.0))) goto cleanup;
            if ((output_rc = ctools_scal_save("__creghdfe_num_singletons", (ST_double)num_singletons))) goto cleanup;
            if ((output_rc = ctools_scal_save("__creghdfe_K_keep", 0.0))) goto cleanup;
            SF_error("creghdfe: all observations are singletons\n");
        } else {
            SF_error("creghdfe: singleton removal failed (out of memory)\n");
        }

        output_rc = num_singletons >= 0 ? 2001 : 920;
        goto cleanup;  /* 2001: all singletons */
    }

    t_singleton = get_time_sec();

    if (compute_groupvar) {
        group_assignments = (ST_int *)calloc(N, sizeof(ST_int));
        if (!group_assignments) { output_rc = 920; goto cleanup; }
    }

    /* ================================================================
     * STEP 3-4: Set up global state for CG solver on the retained rows.
     * Redundant coefficients and mobility groups are counted below from the
     * compacted factors.
     * ================================================================ */
    cleanup_state();
    g_state = (HDFE_State *)calloc(1, sizeof(HDFE_State));

    if (!g_state) {
        output_rc = 920;
        goto cleanup;
    }

    g_state->K = K;
    g_state->in1 = in1;
    g_state->in2 = in2;
    g_state->maxiter = maxiter;
    g_state->tolerance = tolerance;
    g_state->verbose = verbose;
    g_state->num_threads = num_threads;
    g_state->factors_initialized = 1;
    g_state->has_weights = has_weights;
    g_state->weight_type = weight_type;
    g_state->weights = NULL;  /* Will be set after compaction */
    g_state->sum_weights = 0.0;

    /* Factors: surviving levels renumbered 1..num_levels, with counts and the
     * (not yet normalized) weighted counts */
    {
        ST_int *level_ptrs[10];
        ST_int level_max[10];
        for (g = 0; g < G; g++) {
            level_ptrs[g] = factors[g].levels;
            level_max[g] = factors[g].num_levels;
        }
        ST_retcode factor_rc = hdfe_build_factors(g_state, G, level_ptrs, level_max, N_orig, mask, N,
                                                  has_weights ? weights : NULL, 0);
        if (factor_rc) {
            output_rc = factor_rc;
            goto cleanup;
        }
    }

    /* Levels left after singleton removal */
    df_a = 0;
    for (g = 0; g < G; g++) {
        factors[g].num_levels = g_state->factors[g].num_levels;
        df_a += factors[g].num_levels;
    }

    /* Redundant coefficients of each effect, counted as reghdfe counts them
     * (dofadjustments: 0=all, 1=none, 2=firstpair, 3=pairwise). Cluster
     * nesting is not known yet; STEP 9 recounts when it finds a nested effect. */
    ST_int redundant[10] = {0};
    ST_int dof_pairs = dof_adjust_type == 1 ? 0 :
        dof_adjust_type == 2 ? CTOOLS_DOF_FIRSTPAIR : CTOOLS_DOF_PAIRWISE;
    if (compute_dof || !dof_pairs) {
        if (ctools_hdfe_redundant_intercepts(g_state->factors, G, N, NULL, dof_pairs, redundant)) {
            output_rc = 920;
            goto cleanup;
        }
        for (g = 0; g < G; g++) df_a -= redundant[g];
    }

    /* Mobility groups of the first two effects, for e(N_hdfe_extended) and
     * groupvar(); one more is reported for each later effect, as before. */
    if (dof_pairs && compute_dof && G >= 2) {
        mobility_groups = redundant[1];
        if (G > 2 && dof_pairs == CTOOLS_DOF_PAIRWISE) mobility_groups += G - 2;

        if (group_assignments) {
            const ST_int *fe1_compact = g_state->factors[0].levels;
            const ST_int *fe2_compact = g_state->factors[1].levels;
            /* Union-find over the levels of both effects. reghdfe's zigzag
             * search starts each new subgraph at the next unvisited level of
             * the first effect in sorted order, so groups are numbered by the
             * smallest value of the first effect they contain (the compact
             * level numbering follows first appearance instead). */
            if (factors[0].num_levels > INT_MAX - factors[1].num_levels) {
                output_rc = 920; goto cleanup;
            }
            ST_int total_nodes = factors[0].num_levels + factors[1].num_levels;
            ctools_UnionFind *uf = ctools_uf_create(total_nodes);
            root_to_group = (ST_int *)calloc(total_nodes, sizeof(ST_int));
            if (!uf || !root_to_group) {
                ctools_uf_destroy(uf);
                output_rc = 920; goto cleanup;
            }
            for (idx = 0; idx < N; idx++) {
                ctools_uf_union(uf, fe1_compact[idx] - 1,
                                factors[0].num_levels + fe2_compact[idx] - 1);
            }
            const ST_double *fe1_values = filtered.data.vars[K].data.dbl;
            creghdfe_group_min *group_min = (creghdfe_group_min *)malloc(
                (size_t)total_nodes * sizeof(creghdfe_group_min));
            if (!group_min) {
                ctools_uf_destroy(uf);
                output_rc = 920; goto cleanup;
            }
            /* root_to_group temporarily holds 1 + the root's group_min slot */
            ST_int num_roots = 0;
            idx = 0;
            for (ST_int row = 0; row < N_orig; row++) {
                if (!mask[row]) continue;
                ST_int root = ctools_uf_find(uf, fe1_compact[idx++] - 1);
                if (root_to_group[root] == 0) {
                    group_min[num_roots].value = fe1_values[row];
                    group_min[num_roots].root = root;
                    root_to_group[root] = ++num_roots;
                } else if (fe1_values[row] < group_min[root_to_group[root] - 1].value) {
                    group_min[root_to_group[root] - 1].value = fe1_values[row];
                }
            }
            qsort(group_min, (size_t)num_roots, sizeof(creghdfe_group_min), creghdfe_group_min_cmp);
            for (ST_int a = 0; a < num_roots; a++) root_to_group[group_min[a].root] = a + 1;
            free(group_min);
            for (idx = 0; idx < N; idx++) {
                group_assignments[idx] = root_to_group[ctools_uf_find(uf, fe1_compact[idx] - 1)];
            }
            ctools_uf_destroy(uf);
            free(root_to_group); root_to_group = NULL;
        }
    } else if (!dof_pairs || G == 1) {
        if (dof_pairs) mobility_groups = 0;
        /* All observations in group 1 if groupvar requested */
        if (group_assignments) {
            for (idx = 0; idx < N; idx++) {
                group_assignments[idx] = 1;
            }
        }
    }
    for (g = 0; g < G; g++) {
        snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_redundant_%d", g + 1);
        if ((output_rc = ctools_scal_save(scalar_name, (ST_double)redundant[g]))) goto cleanup;
    }

    g_state->df_a = df_a;
    g_state->mobility_groups = mobility_groups;

    t_dof = get_time_sec();

    /* Allocate inv_counts, inv_weighted_counts, and thread buffers */
    g_state->num_threads = num_threads;
    if (ctools_hdfe_alloc_buffers(g_state, 0, K) != 0) {
        output_rc = 920;
        goto cleanup;
    }

    /* ================================================================
     * STEP 5: Compact data (remove singleton rows)
     * ================================================================ */
    data_compact = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K, sizeof(ST_double));
    means_compact = (ST_double *)ctools_safe_malloc2((size_t)K, sizeof(ST_double));

    if (!data_compact || !means_compact) {
        output_rc = 920;
        goto cleanup;
    }

    /* Compact data and recompute means (weighted if using weights).
     * Build compact index first, then parallelize column copy + means. */
    for (k = 0; k < K; k++) means_compact[k] = 0.0;

    /* Build compact index map: compact_idx[j] = destination index for j-th valid obs */
    compact_map = (ST_int *)malloc(N * sizeof(ST_int));
    if (!compact_map) {
        output_rc = 920;
        goto cleanup;
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
        const ST_double *src_col = filtered.data.vars[k].data.dbl;
        ST_double *dst_col = data_compact + (size_t)k * N;
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
    free(compact_map); compact_map = NULL;

    /* Divide by sum of weights (or N if unweighted) */
    ST_double mean_divisor = (has_weights) ? sum_w_compact : (ST_double)N;
    if (!isfinite(mean_divisor) || mean_divisor <= 0.0) {
        output_rc = 498; goto cleanup;
    }
    for (k = 0; k < K; k++) {
        means_compact[k] /= mean_divisor;
    }

    /* The loaded columns are no longer needed; use the compacted matrix */
    stata_data_free(&filtered.data);
    data = data_compact; data_compact = NULL;
    free(means);
    means = means_compact; means_compact = NULL;

    /* Compact weights and assign to g_state */
    if (has_weights) {
        weights_compact = (ST_double *)malloc(N * sizeof(ST_double));
        if (!weights_compact) {
            output_rc = 920;
            goto cleanup;
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

        free(weights); weights = NULL;
        g_state->weights = weights_compact; weights_compact = NULL;
        g_state->sum_weights = sum_w;

    }

    /* ================================================================
     * STEP 6: Compute TSS and stdevs on compacted data (weighted if using weights)
     * Note: for aweight/pweight, weights have been normalized so sum(w) = N
     * ================================================================ */
    for (k = 0; k < K; k++) {
        ST_double *col = &data[(size_t)k * N];
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
            ST_double *col = &data[(size_t)k * N];
            ST_double inv_stdev = 1.0 / stdevs[k];
            #pragma omp simd
            for (idx = 0; idx < N; idx++) {
                col[idx] *= inv_stdev;
            }
        }
    }

    /* Partial out via CG solver using shared helper */
    HDFE_SolveResult projection = partial_out_columns(g_state, data, N, K, num_threads);
    if (projection.status) {
        SF_error("creghdfe: fixed-effect projection did not converge\n");

        output_rc = projection.status;
        goto cleanup;
    }
    ST_int max_iters = projection.iterations;

    /* Save iteration count to Stata scalar */
    if ((output_rc = ctools_scal_save("__creghdfe_iterations", (ST_double)max_iters))) goto cleanup;

    t_partial = get_time_sec();

    /* ================================================================
     * STEP 8: OLS with collinearity detection
     * ================================================================ */
    K_x = K - 1;  /* Excluding y */
    /* For fweight, use sum(weights) as effective N in df calculation (reghdfe.mata line 3594) */
    ST_double N_eff_df = (weight_type == 2) ? g_state->sum_weights : N;
    df_r = N_eff_df - df_a;

    /* Allocate for collinearity detection - cast to size_t to prevent 32-bit overflow */
    xtx = (ST_double *)malloc((size_t)K_x * K_x * sizeof(ST_double));
    is_collinear = (ST_int *)malloc(K_x * sizeof(ST_int));

    if (!xtx || !is_collinear) {
        output_rc = 920;
        goto cleanup;
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
        ST_double xx_partial = fast_dot(&data[(size_t)(k+1) * N], &data[(size_t)(k+1) * N], N);
        ST_double xx_orig = tss[k + 1];
        /* Only mark collinear if essentially zero (both relative and absolute) */
        if (xx_orig > 0 && (xx_partial / xx_orig) <= collinear_tol && xx_partial < 1e-30) {
            is_collinear[k] = 1;
            num_collinear++;
            snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_collinear_varnum_%d", k + 1);
            if ((output_rc = ctools_scal_save(scalar_name, 1.0))) goto cleanup;
        }
    }

    /* Compute X'X */
    ctools_matmul_ata(data + N, N, K_x, xtx);

    /* Detect numerical collinearity via Cholesky */
    ST_int num_numerical_collinear = detect_collinearity(xtx, K_x, is_collinear, verbose);
    if (num_numerical_collinear < 0) {
        output_rc = 920;
        goto cleanup;
    }

    /* Recount collinear */
    num_collinear = 0;
    for (k = 0; k < K_x; k++) {
        if (is_collinear[k]) num_collinear++;
    }

    K_keep = K_x - num_collinear;

    /* Store collinearity flags */
    if ((output_rc = ctools_scal_save("__creghdfe_num_collinear", (ST_double)num_collinear))) goto cleanup;
    for (k = 0; k < K_x; k++) {
        snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_collinear_%d", k + 1);
        if ((output_rc = ctools_scal_save(scalar_name, (ST_double)is_collinear[k]))) goto cleanup;
    }

    ST_int sample_var_idx = 0;
    if (SF_scal_use("__creghdfe_sample_idx", &val) == 0)
        sample_var_idx = (ST_int)val;
    if (sample_var_idx > 0) {
        for (i = 0; i < N_orig && !output_rc; i++) {
            if (mask[i]) output_rc = SF_vstore(sample_var_idx, (ST_int)obs_map[i], 1.0);
        }
    }

    if (output_rc) goto cleanup;

    /* K_keep == 0 (every regressor omitted) continues through the same path,
     * as reghdfe does: the model is then the fixed effects and the constant.
     * K_alloc keeps the K_keep-sized allocations nonempty in that case. */
    const size_t K_alloc = (K_keep > 0) ? (size_t)K_keep : 1;

    /* Build index of non-collinear variables */
    keep_idx = (ST_int *)malloc(K_alloc * sizeof(ST_int));
    if (!keep_idx) {
        output_rc = 920;
        goto cleanup;
    }

    idx = 0;
    for (k = 0; k < K_x; k++) {
        if (!is_collinear[k]) {
            keep_idx[idx++] = k + 1;  /* +1 to skip y */
        }
    }

    /* Allocate for OLS - cast to size_t to prevent 32-bit overflow */
    ST_int K_with_cons = K_keep + 1;
    /* y and the kept regressors. When nothing is omitted this is the data
     * matrix itself and no copy is made. */
    const int keep_is_data = (K_keep == K_x);
    data_keep = keep_is_data ? data
                             : (ST_double *)malloc((size_t)N * (K_keep + 1) * sizeof(ST_double));
    xtx_keep = (ST_double *)malloc(K_alloc * K_alloc * sizeof(ST_double));
    xty_keep = (ST_double *)malloc(K_alloc * sizeof(ST_double));
    beta_keep = (ST_double *)malloc(K_with_cons * sizeof(ST_double));
    inv_xx_keep = (ST_double *)malloc(K_with_cons * K_with_cons * sizeof(ST_double));
    V_keep = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));
    means_x = (ST_double *)malloc(K_alloc * sizeof(ST_double));

    if (!data_keep || !xtx_keep || !xty_keep || !beta_keep || !inv_xx_keep || !V_keep || !means_x) {
        output_rc = 920;
        goto cleanup;
    }

    /* Copy y and non-collinear X */
    if (!keep_is_data) {
        memcpy(data_keep, data, (size_t)N * sizeof(ST_double));
        for (k = 0; k < K_keep; k++) {
            memcpy(data_keep + (size_t)(k + 1) * N, data + (size_t)keep_idx[k] * N,
                   (size_t)N * sizeof(ST_double));
        }
    }
    for (k = 0; k < K_keep; k++) {
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
    inv_xx_x = (ST_double *)malloc(K_alloc * K_alloc * sizeof(ST_double));
    if (!inv_xx_x) {
        output_rc = 920;
        goto cleanup;
    }

    if (K_keep > 0) {
        memcpy(inv_xx_x, xtx_keep, (size_t)K_keep * K_keep * sizeof(ST_double));
        if (ctools_cholesky(inv_xx_x, K_keep) != 0) {
            SF_error("creghdfe: X'X not positive definite\n");

            output_rc = 198;
            goto cleanup;
        }

        stata_retcode inverse_status_0 = ctools_invert_from_cholesky(inv_xx_x, K_keep, inv_xx_x);
        if (inverse_status_0 != STATA_OK) {
            output_rc = ctools_stata_rc(inverse_status_0);
            goto cleanup;
        }
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

    /* Extend inv_xx with the constant (block partition formula). For
     * fweight the 1'W1 term is sum(weights). */
    ST_double N_corner = (weight_type == 2 && has_weights) ? g_state->sum_weights : (ST_double)N;
    ctools_inverse_add_constant(inv_xx_x, K_keep, means_x, N_corner, inv_xx_keep);

    /* Compute constant */
    ST_double y_mean = means[0];
    ST_double xb_mean = 0.0;
    for (k = 0; k < K_keep; k++) {
        xb_mean += means_x[k] * beta_keep[k];
    }
    beta_keep[K_keep] = y_mean - xb_mean;

    free(inv_xx_x); inv_xx_x = NULL;

    /* Build data with constant for the robust and cluster VCEs (the
     * unadjusted VCE does not use it) - cast to size_t to prevent 32-bit overflow */
    if (vcetype == 1 || vcetype == 2) {
        data_with_cons = (ST_double *)malloc((size_t)N * (K_with_cons + 1) * sizeof(ST_double));
        if (!data_with_cons) {
            output_rc = 920;
            goto cleanup;
        }

        /* y (partialled), X with the means added back, and the constant column */
        #pragma omp parallel for schedule(static) if((size_t)N * (size_t)K_with_cons > 200000)
        for (k = 0; k <= K_with_cons; k++) {
            ST_double *dst = data_with_cons + (size_t)k * N;
            if (k == 0) {
                memcpy(dst, data_keep, (size_t)N * sizeof(ST_double));
            } else if (k == K_with_cons) {
                for (ST_int row = 0; row < N; row++) dst[row] = 1.0;
            } else {
                const ST_double *src = data_keep + (size_t)k * N;
                const ST_double mean_k = means_x[k - 1];
                for (ST_int row = 0; row < N; row++) dst[row] = src[row] + mean_k;
            }
        }
    }

    /* ================================================================
     * STEP 9: Compute VCE
     * ================================================================ */

    /* Read cluster variable if clustering */
    ST_double df_a_nested_computed = 0;  /* Track FEs nested within cluster */
    ST_int nested[10] = {0};

    if (vcetype == 2) {
        cluster_ids = (ST_int *)malloc(N * sizeof(ST_int));
        if (!cluster_ids) {
            output_rc = 920;
            goto cleanup;
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
                if (!dof_clusters) {
                    is_nested = 0;
                } else if (g == cluster_matches_fe) {
                    is_nested = 1;  /* FE is the cluster variable - trivially nested */
                } else {
                    is_nested = ctools_fe_nested_in_cluster(
                        g_state->factors[g].levels, g_state->factors[g].num_levels,
                        cluster_ids, N);
                    if (is_nested < 0) { output_rc = 920; goto cleanup; }
                }
                nested[g] = is_nested;
                if (is_nested) df_a_nested_computed += g_state->factors[g].num_levels;
                snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_nested_%d", g + 1);
                if ((output_rc = ctools_scal_save(scalar_name, (ST_double)is_nested))) goto cleanup;
            }
        } else {
            /* Cluster variable is different from all FE variables.
             * Use saved cluster_raw_values and remap with sort (still much faster than SF_vdata + hash). */
            if (!cluster_raw_values) {
                output_rc = 920;
                goto cleanup;
            }

            /* Extract cluster values for non-singleton observations and remap using sort */
            cluster_values_compact = (double *)malloc(N * sizeof(double));

            if (!cluster_values_compact) {
                output_rc = 920;
                goto cleanup;
            }

            /* Compact cluster values (remove singletons) */
            idx = 0;
            for (i = 0; i < N_orig; i++) {
                if (mask[i]) {
                    cluster_values_compact[idx] = cluster_raw_values[i];
                    idx++;
                }
            }

            /* Dense 0-based ids in sorted value order */
            if (ctools_numeric_to_cluster_ids(cluster_values_compact, N, cluster_ids, &num_clusters) != 0) {
                output_rc = 920;
                goto cleanup;
            }

            free(cluster_values_compact); cluster_values_compact = NULL;

            /* Check if any FE is nested within cluster */
            for (g = 0; g < G; g++) {
                ST_int is_nested = !dof_clusters ? 0 : ctools_fe_nested_in_cluster(
                    g_state->factors[g].levels, g_state->factors[g].num_levels,
                    cluster_ids, N);
                if (is_nested < 0) { output_rc = 920; goto cleanup; }
                nested[g] = is_nested;
                if (is_nested) df_a_nested_computed += g_state->factors[g].num_levels;
                snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_nested_%d", g + 1);
                if ((output_rc = ctools_scal_save(scalar_name, (ST_double)is_nested))) goto cleanup;
            }
        }

        /* Free saved cluster raw values */
        free(cluster_raw_values);
        cluster_raw_values = NULL;
    } else {
        /* No clustering - save 0 for all FE nested status */
        for (g = 0; g < G; g++) {
            snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_nested_%d", g + 1);
            if ((output_rc = ctools_scal_save(scalar_name, 0.0))) goto cleanup;
        }
    }

    /* Use computed df_a_nested if we found nested FEs, otherwise use passed value */
    if (df_a_nested_computed > 0) {
        df_a_nested = df_a_nested_computed;

        /* reghdfe sets nested effects aside before it counts the redundant
         * coefficients of the others, so the count depends on the nesting.
         * df_a keeps the nested levels; they are removed through df_a_nested. */
        if (compute_dof || !dof_pairs) {
            if (ctools_hdfe_redundant_intercepts(g_state->factors, G, N, nested, dof_pairs, redundant)) {
                output_rc = 920;
                goto cleanup;
            }
            df_a = 0;
            for (g = 0; g < G; g++) {
                df_a += g_state->factors[g].num_levels - (nested[g] ? 0 : redundant[g]);
                snprintf(scalar_name, sizeof(scalar_name), "__creghdfe_fe_redundant_%d", g + 1);
                if ((output_rc = ctools_scal_save(scalar_name, (ST_double)redundant[g]))) goto cleanup;
            }
            df_r = N_eff_df - df_a;
            g_state->df_a = df_a;
        }
    }

    df_r -= K_keep;

    /* Compute residuals (needed for robust/cluster VCE and optionally stored) */
    resid = (ST_double *)malloc(N * sizeof(ST_double));
    if (!resid) {
        SF_error("creghdfe: memory allocation failed for residuals\n");
        output_rc = 920;
        free(data_with_cons); data_with_cons = NULL;
        goto cleanup;
    }
    if (resid) {
        /* y - X*beta: each fitted value accumulates the products over the
         * regressors in order, in plain double arithmetic (no FMA, as Mata's
         * matrix multiply) */
        memset(resid, 0, (size_t)N * sizeof(ST_double));
        ctools_matvec_add(data_keep + N, N, K_keep, beta_keep, resid);
        #pragma omp parallel for schedule(static) if(N > 200000)
        for (idx = 0; idx < N; idx++) {
            resid[idx] = data_keep[idx] - resid[idx];
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
        output_rc = compute_vce_unadjusted(inv_xx_keep, rss, df_r, K_with_cons, V_keep);
    } else if (vcetype == 1 || vcetype == 2) {
        if (resid) {
            /* Pass weights for weighted VCE (NULL if no weights) */
            ST_double *vce_weights = has_weights ? g_state->weights : NULL;
            /* N_eff = sum(weights) for fweight, else N */
            ST_double N_eff = (weight_type == 2) ? g_state->sum_weights : N;
            if (vcetype == 1) {
                output_rc = compute_vce_robust(data_with_cons, resid, inv_xx_keep, vce_weights, weight_type, N, N_eff, K_with_cons, df_a, V_keep);
            } else {
                ST_int df_m_cluster = K_keep;
                output_rc = compute_vce_cluster(data_with_cons, resid, inv_xx_keep, vce_weights, weight_type, cluster_ids, N, N_eff, K_with_cons, num_clusters, V_keep, df_m_cluster, df_a, df_a_nested);
            }
        }
    }

    free(data_with_cons); data_with_cons = NULL;
    if (output_rc) {
        SF_error("creghdfe: covariance calculation failed\n");
        free(resid); resid = NULL;
        goto cleanup;
    }

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
                if (!output_rc) output_rc = SF_vstore(resid_var_idx, (ST_int)obs_map[i], resid[idx]);
                idx++;
            }
        }
    }

    /* Store groupvar (mobility group assignments) back to Stata if requested */
    if (compute_groupvar && groupvar_var_idx > 0 && group_assignments) {
        idx = 0;
        for (i = 0; i < N_orig; i++) {
            if (mask[i]) {
                if (!output_rc) output_rc = SF_vstore(groupvar_var_idx, (ST_int)obs_map[i], (ST_double)group_assignments[idx]);
                idx++;
            }
        }
    }

    if (savefe && savefe_var_idx > 0 && resid && !output_rc) {
        output_rc = store_fixed_effects(g_state, obs_map, mask, N_orig,
            K_keep, keep_idx, beta_keep, resid, savefe_var_idx);
    }

    /* Free residuals and group assignments */
    free(resid);
    resid = NULL;
    free(group_assignments);
    group_assignments = NULL;

    if (output_rc) goto cleanup;
    t_vce = get_time_sec();

    /* ================================================================
     * STEP 10: Store all results
     * ================================================================ */

    const creghdfe_result result = {
        .weight_type = weight_type,
        .has_weights = has_weights,
        .N = N,
        .num_singletons = num_singletons,
        .G = G,
        .mobility_groups = mobility_groups,
        .K_keep = K_keep,
        .K_with_cons = K_with_cons,
        .vcetype = vcetype,
        .num_clusters = num_clusters,
        .df_a = df_a,
        .df_a_nested = df_a_nested,
        .rss = rss,
        .tss_within = tss_within,
        .tss = tss,
        .beta_keep = beta_keep,
        .V_keep = V_keep,
        .factors = factors,
        .g_state = g_state,
    };
    output_rc = creghdfe_post_results(&result);
    if (output_rc) goto cleanup;

    /* Save timing results to Stata scalars */
    if ((output_rc = ctools_scal_save("_creghdfe_time_load", t_load - t_start))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_copy", t_copy - t_load))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_remap", t_remap - t_copy))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_singleton", t_singleton - t_remap))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_dof", t_dof - t_singleton))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_partial", t_partial - t_dof))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_ols", t_ols - t_partial))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_vce", t_vce - t_ols))) goto cleanup;
    if ((output_rc = ctools_scal_save("_creghdfe_time_total", t_vce - t_start))) goto cleanup;
    CTOOLS_SAVE_THREAD_INFO("_creghdfe");

    /* ================================================================
     * Cleanup
     * ================================================================ */
cleanup:
    ctools_filtered_data_free(&filtered);
    ctools_aligned_free(obs_map);
    for (g = 0; g < G; g++) {
        if (factors) {
            free(factors[g].levels);
            free(factors[g].counts);
        }
    }
    free(factors); free(mask);
    if (data_keep == data) data_keep = NULL;  /* shared when nothing was omitted */
    free(data); free(means); free(stdevs); free(tss);
    free(data_keep); free(xtx_keep); free(xty_keep); free(beta_keep);
    free(inv_xx_keep); free(V_keep); free(keep_idx); free(xtx); free(is_collinear);
    free(cluster_ids);
    free(weights);
    free(cluster_raw_values);
    free(group_assignments);
    free(data_compact);
    free(means_compact);
    free(compact_map);
    free(weights_compact);
    free(means_x);
    free(inv_xx_x);
    free(data_with_cons);
    free(cluster_values_compact);
    free(root_to_group);
    free(resid);
    cleanup_state();
    return output_rc;
}

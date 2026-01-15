/*
 * creghdfe_hdfe.c
 *
 * HDFE initialization, singleton detection, factor management
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_hdfe.h"
#include "creghdfe_utils.h"
#include "../ctools_hdfe_utils.h"

/* Define the global state pointer */
HDFE_State *g_state = NULL;

/*
 * Clean up global state
 */
void cleanup_state(void)
{
    ST_int g, t;

    if (g_state == NULL) return;

    if (g_state->factors != NULL) {
        for (g = 0; g < g_state->G; g++) {
            if (g_state->factors[g].levels != NULL) free(g_state->factors[g].levels);
            if (g_state->factors[g].counts != NULL) free(g_state->factors[g].counts);
            if (g_state->factors[g].weighted_counts != NULL) free(g_state->factors[g].weighted_counts);
            if (g_state->factors[g].means != NULL) free(g_state->factors[g].means);
        }
        free(g_state->factors);
    }

    if (g_state->weights != NULL) free(g_state->weights);

    /* Free per-thread buffers */
    if (g_state->thread_cg_r != NULL) {
        for (t = 0; t < g_state->num_threads; t++) {
            if (g_state->thread_cg_r[t]) free(g_state->thread_cg_r[t]);
        }
        free(g_state->thread_cg_r);
    }
    if (g_state->thread_cg_u != NULL) {
        for (t = 0; t < g_state->num_threads; t++) {
            if (g_state->thread_cg_u[t]) free(g_state->thread_cg_u[t]);
        }
        free(g_state->thread_cg_u);
    }
    if (g_state->thread_cg_v != NULL) {
        for (t = 0; t < g_state->num_threads; t++) {
            if (g_state->thread_cg_v[t]) free(g_state->thread_cg_v[t]);
        }
        free(g_state->thread_cg_v);
    }
    if (g_state->thread_proj != NULL) {
        for (t = 0; t < g_state->num_threads; t++) {
            if (g_state->thread_proj[t]) free(g_state->thread_proj[t]);
        }
        free(g_state->thread_proj);
    }
    if (g_state->thread_fe_means != NULL) {
        for (t = 0; t < g_state->num_threads * g_state->G; t++) {
            if (g_state->thread_fe_means[t]) free(g_state->thread_fe_means[t]);
        }
        free(g_state->thread_fe_means);
    }

    free(g_state);
    g_state = NULL;
}

/* ========================================================================
 * HDFE Initialization command
 * ======================================================================== */

/*
 * Full HDFE initialization (OPTIMIZED VERSION)
 *
 * Varlist: fe1_var fe2_var ... feG_var levels1_var levels2_var ... levelsG_var touse_var
 * First G variables are input FE variables
 * Next G variables are output level variables
 * Last variable is touse
 *
 * Scalars:
 *   __creghdfe_G: number of FE groups
 *   __creghdfe_drop_singletons: whether to drop singletons
 *   __creghdfe_verbose: verbosity level
 *   __creghdfe_compute_dof: whether to compute DOF/mobility groups (optimization)
 *
 * Returns (via scalars):
 *   __creghdfe_N: final number of observations
 *   __creghdfe_num_singletons: total singletons dropped
 *   __creghdfe_num_levels_1, __creghdfe_num_levels_2, ...: levels per FE
 *   __creghdfe_df_a: degrees of freedom absorbed (if compute_dof=1)
 *   __creghdfe_mobility_groups: mobility groups (if compute_dof=1 and G>=2)
 */
ST_retcode do_hdfe_init(int argc, char *argv[])
{
    ST_int G, N, N_orig, in1, in2;
    ST_int g, i, obs, idx;
    ST_double val;
    ST_int drop_singletons, verbose, compute_dof;
    ST_int num_singletons;
    ST_int *mask = NULL;  /* 1 = keep, 0 = drop */
    FactorData *factors = NULL;
    IntHashTable **hash_tables = NULL;
    ST_int *raw_values = NULL;  /* Temporary buffer for single-pass reading */
    char msg[256];
    char scalar_name[64];
    double t_start, t_read, t_singleton, t_dof, t_write;
    ST_int max_iter = 100;  /* Safety limit for singleton iterations */
    ST_int mobility_groups = 1;  /* Default: 1 mobility group */
    ST_int df_a = 0;

    (void)argc;  /* Unused */
    (void)argv;  /* Unused */

    in1 = SF_in1();
    in2 = SF_in2();

    /* Read parameters */
    SF_scal_use("__creghdfe_G", &val); G = (ST_int)val;
    SF_scal_use("__creghdfe_drop_singletons", &val); drop_singletons = (ST_int)val;
    SF_scal_use("__creghdfe_verbose", &val); verbose = (ST_int)val;

    /* Check if DOF computation is requested (optional parameter) */
    compute_dof = 0;
    if (SF_scal_use("__creghdfe_compute_dof", &val) == 0) {
        compute_dof = (ST_int)val;
    }

    if (G < 1 || G > 10) {
        SF_error("creghdfe: invalid number of FE groups\n");
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

    if (verbose >= 2) {
        sprintf(msg, "{txt}   C plugin: HDFE init (N=%d, G=%d, compute_dof=%d)\n", N_orig, G, compute_dof);
        SF_display(msg);
    }

    t_start = get_time_sec();

    /* Allocate structures */
    factors = (FactorData *)calloc(G, sizeof(FactorData));
    hash_tables = (IntHashTable **)calloc(G, sizeof(IntHashTable *));
    mask = (ST_int *)malloc(N_orig * sizeof(ST_int));
    raw_values = (ST_int *)malloc(G * sizeof(ST_int));  /* One value per FE for current obs */

    if (!factors || !hash_tables || !mask || !raw_values) {
        if (factors) free(factors);
        if (hash_tables) free(hash_tables);
        if (mask) free(mask);
        if (raw_values) free(raw_values);
        SF_error("creghdfe: memory allocation failed\n");
        return 1;
    }

    /* Initialize mask (all observations initially kept) */
    for (i = 0; i < N_orig; i++) {
        mask[i] = 1;
    }

    /* Allocate per-factor arrays */
    for (g = 0; g < G; g++) {
        factors[g].levels = (ST_int *)malloc(N_orig * sizeof(ST_int));
        factors[g].num_obs = N_orig;
        hash_tables[g] = hash_create(N_orig);

        if (!factors[g].levels || !hash_tables[g]) {
            /* Cleanup on error */
            for (i = 0; i <= g; i++) {
                if (factors[i].levels) free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
                if (hash_tables[i]) hash_destroy(hash_tables[i]);
            }
            free(factors);
            free(hash_tables);
            free(mask);
            free(raw_values);
            return 1;
        }
    }

    /* STEP 1: Read ALL FE variables in SINGLE PASS (cache-friendly optimization) */
    idx = 0;
    for (obs = in1; obs <= in2; obs++) {
        if (!SF_ifobs(obs)) continue;

        /* Read all G FE values for this observation */
        for (g = 0; g < G; g++) {
            if (SF_vdata(g + 1, obs, &val)) {
                /* Cleanup on error */
                for (i = 0; i < G; i++) {
                    if (factors[i].levels) free(factors[i].levels);
                    if (factors[i].counts) free(factors[i].counts);
                    if (hash_tables[i]) hash_destroy(hash_tables[i]);
                }
                free(factors);
                free(hash_tables);
                free(mask);
                free(raw_values);
                return 198;
            }
            /* Hash and assign level immediately */
            factors[g].levels[idx] = hash_insert_or_get(hash_tables[g], (ST_int)val);
        }
        idx++;
    }

    /* Finalize hash tables and allocate counts */
    for (g = 0; g < G; g++) {
        factors[g].num_levels = hash_tables[g]->size;

        /* Allocate and compute counts */
        factors[g].counts = (ST_int *)calloc(factors[g].num_levels, sizeof(ST_int));
        if (!factors[g].counts) {
            for (i = 0; i < G; i++) {
                if (factors[i].levels) free(factors[i].levels);
                if (factors[i].counts) free(factors[i].counts);
                if (hash_tables[i]) hash_destroy(hash_tables[i]);
            }
            free(factors);
            free(hash_tables);
            free(mask);
            free(raw_values);
            return 1;
        }
        for (i = 0; i < N_orig; i++) {
            factors[g].counts[factors[g].levels[i] - 1]++;
        }
    }

    free(raw_values);
    raw_values = NULL;

    t_read = get_time_sec();

    /* STEP 2: Iteratively drop singletons using shared utility */
    num_singletons = 0;
    N = N_orig;

    if (drop_singletons) {
        /* Build array of level pointers for shared utility */
        ST_int **fe_levels = (ST_int **)malloc(G * sizeof(ST_int *));
        if (fe_levels) {
            for (g = 0; g < G; g++) {
                fe_levels[g] = factors[g].levels;
            }

            num_singletons = ctools_remove_singletons(fe_levels, G, N_orig, mask, max_iter, (verbose >= 2));
            free(fe_levels);

            /* Count remaining observations */
            N = 0;
            for (i = 0; i < N_orig; i++) {
                if (mask[i]) N++;
            }
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

    /* STEP 3: Compute DOF and mobility groups if requested */
    if (compute_dof) {
        /* Compute df_a = sum of all num_levels */
        df_a = 0;
        for (g = 0; g < G; g++) {
            df_a += factors[g].num_levels;
        }

        /* Compute mobility groups for first two FEs (if G >= 2) */
        if (G >= 2) {
            /* Build compact level arrays for remaining observations only */
            ST_int *fe1_compact = (ST_int *)malloc(N * sizeof(ST_int));
            ST_int *fe2_compact = (ST_int *)malloc(N * sizeof(ST_int));

            if (fe1_compact && fe2_compact) {
                /* Remap levels to be contiguous after singleton removal */
                ST_int *remap1 = (ST_int *)calloc(factors[0].num_levels + hash_tables[0]->size + 1, sizeof(ST_int));
                ST_int *remap2 = (ST_int *)calloc(factors[1].num_levels + hash_tables[1]->size + 1, sizeof(ST_int));

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

                    /* Now compute connected components */
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

            /* Subtract mobility groups from df_a */
            df_a -= mobility_groups;
        } else {
            /* Single FE: no redundancy */
            mobility_groups = 0;
        }

        if (verbose >= 2) {
            sprintf(msg, "{txt}   DOF: df_a=%d, mobility_groups=%d\n", df_a, mobility_groups);
            SF_display(msg);
        }
    }

    t_dof = get_time_sec();

    /* STEP 4: Write results back to Stata */
    ST_int touse_var = 2 * G + 1;

    for (g = 0; g < G; g++) {
        idx = 0;
        for (obs = in1; obs <= in2; obs++) {
            if (!SF_ifobs(obs)) continue;

            if (mask[idx]) {
                if (SF_vstore(G + g + 1, obs, (ST_double)factors[g].levels[idx])) {
                    for (i = 0; i < G; i++) {
                        if (factors[i].levels) free(factors[i].levels);
                        if (factors[i].counts) free(factors[i].counts);
                        if (hash_tables[i]) hash_destroy(hash_tables[i]);
                    }
                    free(factors);
                    free(hash_tables);
                    free(mask);
                    return 198;
                }
            } else {
                SF_vstore(G + g + 1, obs, SV_missval);
            }
            idx++;
        }
    }

    /* Update touse variable */
    idx = 0;
    for (obs = in1; obs <= in2; obs++) {
        if (!SF_ifobs(obs)) continue;

        if (!mask[idx]) {
            SF_vstore(touse_var, obs, 0.0);
        }
        idx++;
    }

    t_write = get_time_sec();

    /* STEP 5: Store factors in global state for reuse */
    cleanup_state();
    g_state = (HDFE_State *)calloc(1, sizeof(HDFE_State));

    if (g_state) {
        g_state->G = G;
        g_state->N = N;
        g_state->in1 = in1;
        g_state->in2 = in2;
        g_state->factors_initialized = 1;
        g_state->df_a = df_a;
        g_state->mobility_groups = mobility_groups;

        /* Allocate and copy factors to global state */
        g_state->factors = (FE_Factor *)calloc(G, sizeof(FE_Factor));
        if (g_state->factors) {
            for (g = 0; g < G; g++) {
                g_state->factors[g].num_levels = factors[g].num_levels;
                g_state->factors[g].has_intercept = 1;

                /* Allocate and copy levels (compacted for remaining obs only) */
                g_state->factors[g].levels = (ST_int *)malloc(N * sizeof(ST_int));
                g_state->factors[g].counts = (ST_double *)calloc(factors[g].num_levels, sizeof(ST_double));

                if (g_state->factors[g].levels && g_state->factors[g].counts) {
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
                                idx++;
                            }
                        }
                        free(remap);
                    }
                }

                /* Allocate means buffer */
                g_state->factors[g].means = (ST_double *)malloc(
                    factors[g].num_levels * sizeof(ST_double));
            }
        }
    }

    /* Store results in scalars */
    SF_scal_save("__creghdfe_N", (ST_double)N);
    SF_scal_save("__creghdfe_num_singletons", (ST_double)num_singletons);

    for (g = 0; g < G; g++) {
        sprintf(scalar_name, "__creghdfe_num_levels_%d", g + 1);
        SF_scal_save(scalar_name, (ST_double)factors[g].num_levels);
    }

    if (compute_dof) {
        SF_scal_save("__creghdfe_df_a", (ST_double)df_a);
        SF_scal_save("__creghdfe_mobility_groups", (ST_double)mobility_groups);
    }

    if (verbose >= 2) {
        sprintf(msg, "{txt}   HDFE init: read=%.3fs singleton=%.3fs dof=%.3fs write=%.3fs total=%.3fs\n",
                t_read - t_start, t_singleton - t_read, t_dof - t_singleton, t_write - t_dof, t_write - t_start);
        SF_display(msg);
    }

    /* Cleanup local copies */
    for (g = 0; g < G; g++) {
        if (factors[g].levels) free(factors[g].levels);
        if (factors[g].counts) free(factors[g].counts);
        if (hash_tables[g]) hash_destroy(hash_tables[g]);
    }
    free(factors);
    free(hash_tables);
    free(mask);

    /* Clean up global state - no longer needed after returning results */
    cleanup_state();

    return 0;
}

/* ========================================================================
 * Mobility groups computation command
 * ======================================================================== */

/*
 * Compute mobility groups (connected components) between two FE variables
 */
ST_retcode do_mobility_groups(int argc, char *argv[])
{
    ST_int N, in1, in2;
    ST_int obs, idx;
    ST_double val;
    ST_int *fe1_levels = NULL;
    ST_int *fe2_levels = NULL;
    ST_int num_levels1, num_levels2;
    ST_int verbose;
    ST_int mobility_groups;
    char msg[256];
    double t_start, t_read, t_compute;

    (void)argc;  /* Unused */
    (void)argv;  /* Unused */

    in1 = SF_in1();
    in2 = SF_in2();

    /* Count valid observations */
    N = 0;
    for (obs = in1; obs <= in2; obs++) {
        if (SF_ifobs(obs)) N++;
    }
    if (N <= 0) {
        SF_error("creghdfe: no observations\n");
        return 198;
    }

    /* Read parameters */
    SF_scal_use("__creghdfe_num_levels1", &val); num_levels1 = (ST_int)val;
    SF_scal_use("__creghdfe_num_levels2", &val); num_levels2 = (ST_int)val;
    SF_scal_use("__creghdfe_verbose", &val); verbose = (ST_int)val;

    if (verbose >= 2) {
        sprintf(msg, "{txt}   C plugin: computing mobility groups (N=%d, L1=%d, L2=%d)\n",
                N, num_levels1, num_levels2);
        SF_display(msg);
    }

    t_start = get_time_sec();

    /* Allocate level arrays and hash tables */
    fe1_levels = (ST_int *)malloc(N * sizeof(ST_int));
    fe2_levels = (ST_int *)malloc(N * sizeof(ST_int));
    IntHashTable *ht1 = hash_create(N);
    IntHashTable *ht2 = hash_create(N);

    if (!fe1_levels || !fe2_levels || !ht1 || !ht2) {
        if (fe1_levels) free(fe1_levels);
        if (fe2_levels) free(fe2_levels);
        if (ht1) hash_destroy(ht1);
        if (ht2) hash_destroy(ht2);
        return 1;
    }

    /* Read FE values and hash to level indices */
    idx = 0;
    for (obs = in1; obs <= in2; obs++) {
        if (!SF_ifobs(obs)) continue;

        if (SF_vdata(1, obs, &val)) {
            free(fe1_levels); free(fe2_levels);
            hash_destroy(ht1); hash_destroy(ht2);
            return 198;
        }
        fe1_levels[idx] = hash_insert_or_get(ht1, (ST_int)val);

        if (SF_vdata(2, obs, &val)) {
            free(fe1_levels); free(fe2_levels);
            hash_destroy(ht1); hash_destroy(ht2);
            return 198;
        }
        fe2_levels[idx] = hash_insert_or_get(ht2, (ST_int)val);

        idx++;
    }

    /* Use actual number of unique levels found */
    num_levels1 = ht1->size;
    num_levels2 = ht2->size;

    hash_destroy(ht1);
    hash_destroy(ht2);

    t_read = get_time_sec();

    /* Compute connected components */
    mobility_groups = count_connected_components(
        fe1_levels, fe2_levels, N, num_levels1, num_levels2
    );

    t_compute = get_time_sec();

    if (mobility_groups < 0) {
        free(fe1_levels); free(fe2_levels);
        SF_error("creghdfe: failed to compute mobility groups\n");
        return 1;
    }

    if (verbose >= 2) {
        sprintf(msg, "{txt}   Mobility groups: %d (read=%.3fs, compute=%.3fs)\n",
                mobility_groups, t_read - t_start, t_compute - t_read);
        SF_display(msg);
    }

    /* Store result */
    SF_scal_save("__creghdfe_mobility_groups", (ST_double)mobility_groups);

    free(fe1_levels);
    free(fe2_levels);

    return 0;
}

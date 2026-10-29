/*
 * cppmlhdfe_irls.c
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

/* Resolve output indices in the plugin varlist, including hidden Stata variables. */
#define SD_SAFEMODE
#include <float.h>
#include <limits.h>
#include "cppmlhdfe_irls.h"
#include "cppmlhdfe_result.h"
#include "cppmlhdfe_input.h"
#include "cppmlhdfe_internal.h"
#include "cppmlhdfe_separation.h"
#include "hdfe/hdfe_remap.h"
#include "hdfe/hdfe_solver.h"
#include "hdfe/hdfe_sample.h"
#include "../ctools_matrix.h"
#include "../ctools_hdfe_utils.h"
#include "../ctools_config.h"
#include "../ctools_types.h"
#include "../ctools_spi.h"
#include "../ctools_ols.h"

/* Forward declarations for helpers defined at bottom of file */

/* Reference CG stopping controls are essential in the separation problem:
 * a near-singular projection can otherwise change the retained sample. */


/* Rank-revealing weighted QR, with a second orthogonalization pass. */

/* Recover raw per-level intercept/slope coefficients from the saved sum.
 * All factors participate even when only some outputs were requested. */
static int ppml_save_effects(HDFE_State *S,double *d,const perm_idx_t *obs_map)
{
    int targets[10],any=0,rc=0;
    double *alpha[10]={0};
    for(int g=0;g<S->G;g++) {
        char key[32];snprintf(key,sizeof(key),"save_fe_%d",g+1);
        targets[g]=(int)ppml_option(key,0);any|=targets[g]>0;
    }
    if(!any)return 0;
    for(int g=0;g<S->G;g++) {
        alpha[g]=calloc((size_t)S->factors[g].num_levels,sizeof(double));
        if(!alpha[g]){rc=920;goto done;}
    }
    double scale=1;
    for(int i=0;i<S->N;i++)scale=fmax(scale,fabs(d[i]));
    int iteration;
    for(iteration=0;iteration<S->maxiter;iteration++) {
        for(int g=0;g<S->G;g++) {
            FE_Factor *f=&S->factors[g];int joint=ppml_joint_slope(S,g);
            ppml_project(S,d,g);
            double *a=S->thread_fe_means[g];
            for(int l=0;l<f->num_levels;l++)alpha[g][l]+=a[l];
            if(joint) {
                FE_Factor *h=&S->factors[g+1];double *b=S->thread_fe_means[g+1];
                for(int l=0;l<f->num_levels;l++){alpha[g][l]-=h->slope_center[l]*b[l];alpha[g+1][l]+=b[l];}
                g++;
            }
            else if(f->slope_center)for(int h=0;h<S->G;h++) {
                FE_Factor *intercept=&S->factors[h];
                if(!intercept->slope && intercept->num_levels==f->num_levels &&
                    !memcmp(intercept->levels,f->levels,(size_t)S->N*sizeof(int))) {
                    for(int l=0;l<f->num_levels;l++)alpha[h][l]-=f->slope_center[l]*a[l];
                    break;
                }
            }
        }
        double error=0;for(int i=0;i<S->N;i++)error=fmax(error,fabs(d[i]));
        if(error<=fmax(S->tolerance,1e-10)*scale)break;
    }
    if(iteration==S->maxiter){rc=430;goto done;}
    for(int g=0;g<S->G && !rc;g++)if(targets[g]) {
        for(int i=0;i<S->N && !rc;i++)rc=SF_vstore(targets[g],(int)obs_map[i],alpha[g][S->factors[g].levels[i]-1]);
    }
done:
    for(int g=0;g<S->G;g++)free(alpha[g]);
    return rc;
}

/* X'WX (K x K, symmetric) and X'Wz in double-double precision, matching
 * ppmlhdfe/reghdfe's quadcross(): small differences in each IRLS iteration
 * accumulate in the mu trajectory and are amplified through cond(X'WX) in the
 * final VCE. data is [z, X] column-major with N rows; cols lists the K
 * regressors used (NULL = the first K). Each element is an independent
 * sequential sum, so the rows are computed in parallel and the sums of a
 * row are interleaved. xty_dd, when not NULL, also receives X'Wz as
 * unrounded double-doubles. Returns 0, or 920 on allocation failure. */
static ST_retcode ppml_normal_equations(const ST_double *data, ST_int N, const ST_int *cols,
                                        ST_int K, const ST_double *w, ST_double *xtx,
                                        ST_double *xty, dd_real *xty_dd)
{
    ST_int a;
    int failed = 0;

    #pragma omp parallel if(K > 1 && (size_t)N * (size_t)K > 100000)
    {
        /* per-thread list of right-hand vectors: z, then x_0 .. x_a */
        const ST_double **targets = (const ST_double **)malloc(((size_t)K + 1) * sizeof(*targets));
        dd_real *sums = (dd_real *)malloc(((size_t)K + 1) * sizeof(*sums));
        if (!targets || !sums) {
            #pragma omp atomic write
            failed = 1;
        }
        #pragma omp for schedule(dynamic)
        for (a = 0; a < K; a++) {
            if (!targets || !sums) continue;
            const ST_double *xa = data + (size_t)((cols ? cols[a] : a) + 1) * N;
            targets[0] = data;
            for (ST_int b = 0; b <= a; b++)
                targets[b + 1] = data + (size_t)((cols ? cols[b] : b) + 1) * N;
            dd_dot_weighted_multi(xa, targets, a + 2, w, N, sums);
            if (xty_dd) xty_dd[a] = sums[0];
            xty[a] = sums[0].hi + sums[0].lo;
            for (ST_int b = 0; b <= a; b++) {
                xtx[a * K + b] = sums[b + 1].hi + sums[b + 1].lo;
                xtx[b * K + a] = sums[b + 1].hi + sums[b + 1].lo;
            }
        }
        free(targets);
        free(sums);
    }
    return failed ? 920 : 0;
}

/* Global state pointer for CG solver (same pattern as creghdfe) */
static HDFE_State *g_ppml_state = NULL;

static void cleanup_ppml_state(void)
{
    if (g_ppml_state) {
        free(g_ppml_state->weights);
        g_ppml_state->weights = NULL;
        ctools_hdfe_state_cleanup(g_ppml_state);
        free(g_ppml_state);
        g_ppml_state = NULL;
    }
}

void cppmlhdfe_cleanup_state(void) { cleanup_ppml_state(); }


/* ========================================================================
 * Main PPML regression
 * ======================================================================== */

ST_retcode do_ppml_regression(int argc, char *argv[])
{
    ST_int K, G, N_orig, N, in1, in2;
    ST_retcode output_rc = 0;
    ST_int sample_var_idx = 0;
    ST_int k, g, i, j, idx;
    ST_double val;
    ST_int verbose;
    ST_int num_singletons;
    ST_int *mask = NULL;
    PPMLFactorData *factors = NULL;
    char scalar_name[64];
    double t_start, t_load, t_remap, t_singleton, t_dof, t_irls, t_vce;
    ST_int mobility_groups = 1;
    ST_int df_a = 0;
    ST_int num_threads;

    /* Data arrays */
    ST_double *data = NULL;     /* N x K matrix (column-major) */
    ST_int *cluster_ids = NULL;
    ST_int num_clusters = 0;
    perm_idx_t *obs_map = NULL;

    /* Heap state released at ppml_cleanup; every error path leaves through it */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    ST_double *weights = NULL, *offset_arr = NULL, *data_compact = NULL;
    ST_double *w_user_compact = NULL, *offset_compact = NULL;
    double *cluster_raw_values = NULL;
    ST_double *mu = NULL, *eta = NULL, *irls_w = NULL, *z_orig = NULL, *aug_copy = NULL;
    ST_double *xtx = NULL, *xty = NULL, *irls_beta = NULL, *inv_xx = NULL;
    ST_int *sep_m = NULL;
    ST_double *xb_work = NULL;
    dd_real *xty_dd = NULL;
    ST_double *stdev_x = NULL;
    ST_int *irls_keep_idx = NULL;
    ST_double *xtx_k = NULL;
    ST_double *xty_k = NULL;
    ST_double *beta_k = NULL;
    ST_double *inv_xx_k = NULL;
    ST_int *is_collinear = NULL;
    ST_int *keep_idx = NULL;
    ST_double *beta_final = NULL;
    ST_double *inv_xx_final = NULL;
    ST_double *V_keep = NULL;
    ST_double *xtx_keep = NULL;
    ST_double *xty_keep = NULL;
    ST_double *w_reg = NULL;
    ST_double *means_x = NULL;
    ST_double *inv_xx_ext = NULL;
    ST_double *X_eff = NULL;
    ST_double *vce_resid = NULL;
    ST_double *V_ext = NULL;
    ST_double *separation_certificate = NULL;
    ST_double *previous_working = NULL;

    (void)argc;
    (void)argv;

    PPML_Options configuration;
    ppml_read_options(&configuration);
    const PPML_Options *options = &configuration;
    output_rc = ctools_scal_save("__cppmlhdfe_api", 3);
    if (output_rc) return output_rc;
    t_start = get_time_sec();

    in1 = SF_in1();
    in2 = SF_in2();

    cppmlhdfe_input input;
    ST_retcode input_rc = cppmlhdfe_read_input(&input);
    if (input_rc) return input_rc;
    K = input.K;
    G = input.G;
    verbose = input.verbose;
    ST_int maxiter = input.maxiter;
    ST_int vcetype = input.vcetype;
    ST_int cluster_terms = input.cluster_terms;
    ST_int has_weights = input.has_weights;
    ST_int weight_type = input.weight_type;
    ST_int has_offset = input.has_offset;
    ST_int irls_maxiter = input.irls_maxiter;
    ST_int dof_adjust_type = input.dof_adjust_type;
    ST_int dof_flags = input.dof_flags;
    ST_int keep_singletons = input.keep_singletons;
    ST_int check_fe = input.check_fe;
    ST_int check_relu = input.check_relu;
    ST_double tolerance = input.tolerance;
    ST_double irls_tol = input.irls_tol;
    ST_int slope_only[10], slope_idx[10];
    memcpy(slope_only, input.slope_only, sizeof(slope_only));
    memcpy(slope_idx, input.slope_idx, sizeof(slope_idx));

    /* Validation */
    if (G < 1 || G > 10) {
        SF_error("cppmlhdfe: invalid number of FE groups (must be 1-10)\n");
        return 198;
    }
    if (K < 2) {
        SF_error("cppmlhdfe: need at least depvar and one indepvar (K < 2)\n");
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
    ST_int total_vars = K + G + cluster_terms + (has_weights ? 1 : 0) + (has_offset ? 1 : 0);

    int *var_indices = (int *)malloc(total_vars * sizeof(int));
    if (!var_indices) {
        SF_error("cppmlhdfe: memory allocation failed\n");
        return 920;
    }
    for (i = 0; i < total_vars; i++) {
        var_indices[i] = i + 1;
    }


    int load_flags = CTOOLS_LOAD_NO_SORT_ORDER | (input.touse ? CTOOLS_LOAD_MARKER : CTOOLS_LOAD_CHECK_IF);
    stata_retcode load_rc = ctools_data_load(&filtered, var_indices, total_vars, 0, 0, load_flags);
    free(var_indices);

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        SF_error("cppmlhdfe: parallel data load failed\n");
        return ctools_stata_rc(load_rc);
    }

    N_orig = (ST_int)filtered.data.nobs;
    obs_map = filtered.obs_map;
    filtered.obs_map = NULL;  /* obs_map is now owned by this command */

    if (N_orig <= 0) {
        SF_error("cppmlhdfe: no observations\n");
        output_rc = 2000;
        goto ppml_cleanup;
    }

    t_load = get_time_sec();

    /* Allocate working arrays */
    factors = (PPMLFactorData *)calloc(G, sizeof(PPMLFactorData));
    mask = (ST_int *)malloc(N_orig * sizeof(ST_int));
    data = (ST_double *)malloc((size_t)N_orig * K * sizeof(ST_double));

    if (has_weights) {
        weights = (ST_double *)malloc(N_orig * sizeof(ST_double));
    }

    if (has_offset) {
        offset_arr = (ST_double *)malloc(N_orig * sizeof(ST_double));
    }

    if (!factors || !mask || !data || (has_weights && !weights) || (has_offset && !offset_arr)) {
        SF_error("cppmlhdfe: memory allocation failed\n");
        output_rc = 920;
        goto ppml_cleanup;
    }

    for (i = 0; i < N_orig; i++) mask[i] = 1;

    /* Copy numeric data to column-major */
    for (k = 0; k < K; k++) {
        double *src = filtered.data.vars[k].data.dbl;
        double *dst = &data[(size_t)k * N_orig];
        memcpy(dst, src, N_orig * sizeof(double));
    }

    /* Copy weights */
    if (has_weights) {
        ST_int weight_pos = K + G + cluster_terms;
        memcpy(weights, filtered.data.vars[weight_pos].data.dbl, N_orig * sizeof(ST_double));
    }

    /* Copy offset */
    if (has_offset) {
        ST_int offset_pos = K + G + cluster_terms + (has_weights ? 1 : 0);
        memcpy(offset_arr, filtered.data.vars[offset_pos].data.dbl, N_orig * sizeof(ST_double));
    }

    /* Validate y >= 0 */
    {
        ST_double *y = data;  /* First column */
        for (i = 0; i < N_orig; i++) {
            if (y[i] < 0.0) {
                SF_error("cppmlhdfe: dependent variable must be non-negative\n");
                output_rc = 198;
                goto ppml_cleanup;
            }
        }
    }

    /* FE variables to 1-based levels with counts, one factor per thread */
    ST_retcode remap_rc;
    {
        const ST_double *fe_values[10];
        for (g = 0; g < G; g++) fe_values[g] = filtered.data.vars[K + g].data.dbl;
        remap_rc = hdfe_load_factors(factors, G, fe_values, N_orig);
    }

    /* Save cluster raw values before freeing filtered data */
    if (vcetype == 2) {
        cluster_raw_values = malloc((size_t)N_orig * cluster_terms * sizeof(double));
        if (cluster_raw_values) {
            for (ST_int term = 0; term < cluster_terms; term++)
                memcpy(cluster_raw_values + (size_t)term * N_orig,
                       filtered.data.vars[K + G + term].data.dbl,
                       N_orig * sizeof(double));
        }
    }

    /* Free filtered data, keep obs_map */
    stata_data_free(&filtered.data);

    if (remap_rc) {
        SF_error("cppmlhdfe: FE remapping failed\n");
        output_rc = remap_rc;
        goto ppml_cleanup;
    }

    t_remap = get_time_sec();

    /* ================================================================
     * STEP 2: Singleton removal
     * ================================================================ */
    ST_int num_separated = 0;
    ST_int *selection_levels[10];
    ST_int selection_nlevels[10];
    for (g = 0; g < G; g++) {
        selection_levels[g] = factors[g].levels;
        selection_nlevels[g] = factors[g].num_levels;
    }
    ST_int selection_rc = ppml_select_fe_sample_ex(data, selection_levels,
        selection_nlevels, G, N_orig, mask, &num_singletons, &num_separated,
        slope_only, keep_singletons, check_fe && !keep_singletons);
    /* ppmlhdfe counts the initial FE screen as singleton removal. Explicit
     * regressor/rectifier separation is reported separately. */
    num_singletons += num_separated;
    num_separated = 0;
    N = 0;
    for (i = 0; i < N_orig; i++) if (mask[i]) N++;
    if (selection_rc || N == 0 || (vcetype == 2 && !cluster_raw_values)) {
        if (selection_rc || N > 0) {
            SF_error("cppmlhdfe: sample selection allocation failed\n");
            output_rc = 920;
        } else {
            SF_error("cppmlhdfe: no observations remain after singleton/separation removal\n");
            output_rc = 2001;
        }
        goto ppml_cleanup;
    }
    if (verbose && (num_singletons || num_separated)) {
        char msg[160];
        snprintf(msg, sizeof(msg), "cppmlhdfe: dropped %d singletons and %d separated observations\n",
                 num_singletons, num_separated);
        SF_display(msg);
    }

    t_singleton = get_time_sec();

    /* ================================================================
     * STEP 3-4: Build HDFE state on the retained rows, then count the
     * absorbed degrees of freedom from its factors
     * ================================================================ */
    cleanup_ppml_state();
    g_ppml_state = (HDFE_State *)calloc(1, sizeof(HDFE_State));
    if (!g_ppml_state) {
        output_rc = 920;
        goto ppml_cleanup;
    }

    g_ppml_state->K = K;
    g_ppml_state->in1 = in1;
    g_ppml_state->in2 = in2;
    g_ppml_state->maxiter = maxiter;
    g_ppml_state->tolerance = tolerance;
    g_ppml_state->verbose = verbose;
    g_ppml_state->num_threads = num_threads;
    g_ppml_state->factors_initialized = 1;
    g_ppml_state->has_weights = 1;  /* IRLS always uses weights */
    g_ppml_state->weight_type = 1;  /* Treat as aweight for CG solver */
    g_ppml_state->weights = NULL;
    g_ppml_state->sum_weights = 0.0;

    /* Factors on the retained rows: levels renumbered 1..num_levels, with
     * zeroed weighted counts for the IRLS weights set later */
    ST_retcode factor_rc = hdfe_build_factors(g_ppml_state, G, selection_levels, selection_nlevels,
                                              N_orig, mask, N, NULL, 1);
    if (factor_rc) {
        output_rc = factor_rc;
        goto ppml_cleanup;
    }

    /* Levels left after selection; mobility groups of the first pair */
    for (g = 0; g < G; g++)
        factors[g].num_levels = g_ppml_state->factors[g].num_levels;
    ctools_compute_hdfe_dof(g_ppml_state->factors, G, N, dof_adjust_type, &df_a, &mobility_groups);
    g_ppml_state->df_a = df_a;
    g_ppml_state->mobility_groups = mobility_groups;

    t_dof = get_time_sec();

    /* ================================================================
     * STEP 5: Compact data (remove singleton rows)
     * ================================================================ */
    data_compact = (ST_double *)malloc((size_t)N * K * sizeof(ST_double));

    if (has_weights) w_user_compact = (ST_double *)malloc(N * sizeof(ST_double));
    if (has_offset) offset_compact = (ST_double *)malloc(N * sizeof(ST_double));

    if (!data_compact || (has_weights && !w_user_compact) || (has_offset && !offset_compact)) {
        output_rc = 920;
        goto ppml_cleanup;
    }

    /* Compact data columns */
    hdfe_compact_columns(data, N_orig, K, mask, data_compact, N);

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
    data_compact = NULL;
    if (weights) free(weights);
    weights = NULL;
    if (offset_arr) free(offset_arr);
    offset_arr = NULL;

    /* Compact all row identities with exactly the same original-row mask.
     * After this point every array, including clusters and obs_map, has N rows. */
    idx = 0;
    for (i = 0; i < N_orig; i++) {
        if (mask[i]) {
            obs_map[idx] = obs_map[i];
            if (cluster_raw_values) {
                for (ST_int term = 0; term < cluster_terms; term++)
                    cluster_raw_values[(size_t)term * N_orig + idx] =
                        cluster_raw_values[(size_t)term * N_orig + i];
            }
            idx++;
        }
    }
    free(mask);
    mask = NULL;

    /* ================================================================
     * STEP 7: Allocate IRLS arrays and HDFE buffers
     * ================================================================ */
    mu = (ST_double *)malloc(N * sizeof(ST_double));
    eta = (ST_double *)malloc(N * sizeof(ST_double));
    irls_w = (ST_double *)malloc(N * sizeof(ST_double));
    g_ppml_state->weights = irls_w;  /* owned by the HDFE state from here on */
    /* z_orig: saves pre-partialled working depvar for FE reconstruction.
     * aug_copy: [z, X1, ..., Xk] for CG partialling (modified in place). */
    ST_int K_x = K - 1;
    ST_int K_aug = K_x + 1;  /* working depvar + X vars */
    z_orig = (ST_double *)malloc(N * sizeof(ST_double));
    aug_copy = (ST_double *)malloc((size_t)N * K_aug * sizeof(ST_double));

    /* Pre-allocate IRLS working buffers (reused across iterations) */
    xtx = (ST_double *)malloc(K_x * K_x * sizeof(ST_double));
    xty = (ST_double *)malloc(K_x * sizeof(ST_double));
    irls_beta = (ST_double *)malloc(K_x * sizeof(ST_double));
    inv_xx = (ST_double *)malloc(K_x * K_x * sizeof(ST_double));
    sep_m = (ST_int *)malloc(N * sizeof(ST_int));
    xb_work = (ST_double *)malloc(N * sizeof(ST_double));

    if (!mu || !eta || !irls_w || !z_orig || !aug_copy ||
        !xtx || !xty || !irls_beta || !inv_xx || !sep_m || !xb_work) {
        output_rc = 920;
        goto ppml_cleanup;
    }

    /* Continuous FE loadings follow the same compact sample as y and X. */
    for (g = 0; g < G; g++) {
        if (!slope_only[g]) continue;
        FE_Factor *f = &g_ppml_state->factors[g];
        f->slope = malloc((size_t)N*sizeof(*f->slope));
        if (!f->slope) {
            /* Allocation failure uses the common path once IRLS owns weights. */
            output_rc = 920;
            break;
        }
        f->has_intercept = 0;
        f->num_slopes = 1;
        for (j = 0; j < G; j++) {
            FE_Factor *intercept = &g_ppml_state->factors[j];
            if (!slope_only[j] && f->num_levels == intercept->num_levels &&
                !memcmp(f->levels, intercept->levels, (size_t)N*sizeof(ST_int))) {
                f->slope_center = calloc((size_t)f->num_levels, sizeof(ST_double));
                if (!f->slope_center) output_rc = 920;
                break;
            }
        }
        if (output_rc) break;
        for (i = 0; i < N; i++) {
            ST_retcode read_rc = SF_vdata(slope_idx[g], (ST_int)obs_map[i], &f->slope[i]);
            if (read_rc) { output_rc = read_rc; break; }
        }
        if (output_rc) break;
    }

    /* Allocate HDFE buffers (alloc_proj=0, max_columns=K_aug) */
    if (ctools_hdfe_alloc_buffers(g_ppml_state, 0, K_aug) != 0) {
        output_rc = 920;
        goto ppml_cleanup;
    }

    if (output_rc) goto ppml_cleanup;
    /* ================================================================
     * STEP 8: Initialize IRLS
     * Initialize eta from OLS of log(max(y,1)) on X with unit weights
     * ================================================================ */
    ST_double *y = data;  /* First column of data */
    ST_double *X = data + N;  /* Columns 1..K-1 */

    /* Standardize by sample stdev before IRLS. Retain the reference's 1e-3
     * regressor scale floor, but keep response scaling invariant even for
     * very small outcomes. */
    ST_double stdev_y = 1.0;
    stdev_x = (ST_double *)malloc(K_x * sizeof(ST_double));
    if (!stdev_x) { output_rc = 920; goto ppml_cleanup; }
    {
        /* Center before squaring: E[x^2]-E[x]^2 loses all precision when
         * the mean is large relative to the within-sample variation. */
        for (k = -1; k < K_x; k++) {
            ST_double *column = k < 0 ? y : X + (size_t)k * N;
            dd_real total = {0.0, 0.0};
            for (i = 0; i < N; i++) total = dd_add_d(total, column[i] / N);
            ST_double mean = total.hi + total.lo;
            ST_double max_delta = 0.0;
            for (i = 0; i < N; i++) max_delta = fmax(max_delta, fabs(column[i] - mean));
            /* A constant column has no variation. The rounded sum of
             * column[i]/N need not reproduce the constant exactly, and a
             * roundoff-level max_delta would otherwise give an outcome scale
             * of order 1e-16 (amplifying every later rounding error). */
            for (i = 1; i < N && column[i] == column[0]; i++) {}
            if (N > 0 && i == N) { mean = column[0]; max_delta = 0.0; }
            dd_real squares = {0.0, 0.0};
            if (max_delta > 0.0) {
                for (i = 0; i < N; i++) {
                    ST_double delta = (column[i] - mean) / max_delta;
                    squares = dd_add_d(squares, delta * delta);
                }
            }
            ST_double sd = max_delta * sqrt((squares.hi + squares.lo) / (N - 1));
            /* The response scale must not change the fit or trigger spurious
             * separation for small positive outcomes. */
            if (k < 0) {
                stdev_y = sd > 0.0 ? sd : fmax(fabs(mean), 1e-3);
                sd = stdev_y;
            } else {
                sd = fmax(sd, 1e-3);
                stdev_x[k] = sd;
            }
            if(!options->standardize) {
                sd=1;
                if(k<0)stdev_y=1;else stdev_x[k]=1;
            }
            ST_double fixed_scales = 0;
            SF_scal_use("__cppmlhdfe_scale_fixed", &fixed_scales);
            if (fixed_scales) {
                SF_mat_el("__cppmlhdfe_input_scales", 1, k+2, &sd);
                if (k < 0) stdev_y = sd; else stdev_x[k] = sd;
            } else if ((output_rc = ctools_mat_store("__cppmlhdfe_input_scales", 1, k+2, sd)))
                goto ppml_cleanup;
            for (i = 0; i < N; i++) column[i] /= sd;
        }
    }

    separation_certificate=calloc((size_t)N,sizeof(double));
    previous_working=calloc((size_t)N,sizeof(double));
    if(!separation_certificate || !previous_working) {output_rc=920;goto ppml_cleanup;}
    if(options->check_simplex) {
        int nsep=ppml_simplex(options, g_ppml_state,y,X,w_user_compact,weight_type,N,K_x,sep_m);
        if(nsep<0){output_rc=-nsep;goto ppml_cleanup;}
        if(nsep) {
            sample_var_idx=(int)ppml_option("sample_idx",0);
            for(i=0;i<N && !output_rc;i++)if(!sep_m[i])output_rc=SF_vstore(sample_var_idx,(int)obs_map[i],1);
            if(!output_rc)output_rc=ctools_scal_save("__cppmlhdfe_refit",1);
            if(!output_rc)output_rc=ctools_scal_save("__cppmlhdfe_num_singletons",num_singletons);
            if(!output_rc)output_rc=ctools_scal_save("__cppmlhdfe_num_separated",num_separated+nsep);
            goto ppml_cleanup;
        }
    }
    ST_int relu_separated = check_relu ? ppml_relu_ex(options, g_ppml_state, y, X, N, K_x, sep_m, num_threads,separation_certificate) : 0;
    if (relu_separated < 0) { output_rc = -relu_separated; goto ppml_cleanup; }
    /* ReLU diagnostics refer to its input sample, excluding rows already
     * removed by FE screening, singleton removal, or simplex. */
    if(check_relu) for(i=0;i<N && !output_rc;i++) {
        if(options->tag_idx)output_rc=SF_vstore(options->tag_idx,(int)obs_map[i],sep_m[i]);
        if(!output_rc && options->certificate_idx)output_rc=SF_vstore(options->certificate_idx,(int)obs_map[i],separation_certificate[i]);
    }
    if(output_rc || ppml_option("diagnostic",0))goto ppml_cleanup;
    if (relu_separated) {
        if (verbose) {
            char message[128];
            snprintf(message, sizeof(message), "cppmlhdfe: ReLU identified %d separated observations\n", relu_separated);
            SF_display(message);
        }
        if (SF_scal_use("__cppmlhdfe_sample_idx", &val) == 0) sample_var_idx = (ST_int)val;
        for (i = 0; i < N && !output_rc; i++)
            if (!sep_m[i]) output_rc = SF_vstore(sample_var_idx, (ST_int)obs_map[i], 1.0);
        if (!output_rc) output_rc = ctools_scal_save("__cppmlhdfe_refit", 2);
        if (!output_rc) output_rc = ctools_scal_save("__cppmlhdfe_num_singletons", num_singletons);
        if (!output_rc) output_rc = ctools_scal_save("__cppmlhdfe_num_separated", num_separated + relu_separated);
        goto ppml_cleanup;
    }

    /* Initial means: simple, OLS, or a user-supplied mean variable. */
    {
        double ym=0,ws=0;
        for(i=0;i<N;i++){double w=w_user_compact?w_user_compact[i]:1;ym+=w*y[i];ws+=w;irls_w[i]=w;}
        ym/=ws;
        if(options->guess==1) {
            for(i=0;i<N;i++)aug_copy[i]=log(y[i]+ym/100);
            memcpy(aug_copy+N,X,(size_t)N*K_x*sizeof(double));
            ppml_update_fe_weights(g_ppml_state,irls_w,N);
            g_ppml_state->tolerance=fmax(options->start_tol,tolerance);
            HDFE_SolveResult pr=ppml_partial(options, g_ppml_state,aug_copy,N,K_aug,num_threads);
            if(pr.status){output_rc=pr.status;goto ppml_cleanup;}
            output_rc=ppml_qrsolve(aug_copy+N,aug_copy,irls_w,N,K_x,NULL,mu);
            if(output_rc)goto ppml_cleanup;
            for(i=0;i<N;i++)mu[i]=exp(log1p(y[i])-mu[i]);
        }
        else if(options->guess==2) {
            double average=0,sd=0;
            for(i=0;i<N;i++) {
                output_rc=SF_vdata(options->guess_idx,(int)obs_map[i],&mu[i]);
                if(output_rc || !isfinite(mu[i]) || mu[i]>=SV_missval){output_rc=198;goto ppml_cleanup;}
                average+=mu[i]/N;
            }
            if(options->standardize) {
                for(i=0;i<N;i++)sd+=(mu[i]-average)*(mu[i]-average);
                sd=fmax(sqrt(sd/fmax(N-1,1)),1e-3);
                for(i=0;i<N;i++)mu[i]/=sd;
            }
        }
        else for(i=0;i<N;i++)mu[i]=.5*(y[i]+ym);
        for(i=0;i<N;i++) {
            mu[i]=fmax(mu[i],fmax(.05*y[i],1e-3));
            eta[i]=log(mu[i])-(offset_compact?offset_compact[i]:0);
        }
    }

    ST_double user_weight_scale = 1.0;
    if (w_user_compact) {
        user_weight_scale = 0.0;
        for (i = 0; i < N; i++) user_weight_scale += w_user_compact[i] / N;
    }

    /* A common scale cancels from weighted projections and the sandwich. */
    for (i = 0; i < N; i++) {
        ST_double w_u = (w_user_compact != NULL) ? w_user_compact[i] : 1.0;
        irls_w[i] = mu[i] * (w_u / user_weight_scale);
        if (irls_w[i] < DBL_MIN) irls_w[i] = DBL_MIN;  /* Floor to prevent zero weights */
    }

    /* Set weighted_counts for initial weights */
    ppml_update_fe_weights(g_ppml_state, irls_w, N);

    /* ================================================================
     * STEP 8b: Pre-IRLS collinearity detection
     * Partial out FE from X with initial weights, then detect collinearity
     * in X'X to identify variables collinear with the absorbed FEs.
     * ================================================================ */

    ST_int K_irls = K_x;  /* number of non-collinear X vars for IRLS */
    if(options->remove_collinear) {
        /* Set up [z_dummy, X] for partialling - z doesn't matter, just X */
        for (i = 0; i < N; i++)
            aug_copy[i] = 0.0;  /* dummy z */
        for (k = 0; k < K_x; k++)
            memcpy(&aug_copy[(size_t)(k + 1) * N], &X[(size_t)k * N], N * sizeof(ST_double));

        HDFE_SolveResult projection = ppml_partial(options, g_ppml_state, aug_copy, N, K_aug, num_threads);
        if (projection.status) {
            SF_error("cppmlhdfe: fixed-effect projection did not converge\n");
            output_rc = projection.status;
            goto ppml_cleanup;
        }

        /* Build X'X from partialled X columns */
        ST_double *xtx_pre = (ST_double *)calloc(K_x * K_x, sizeof(ST_double));
        if (!xtx_pre) { output_rc = 920; goto ppml_cleanup; }
        if (xtx_pre) {
            ctools_matmul_ata(aug_copy + N, N, K_x, xtx_pre);
            ST_int *pre_collinear = (ST_int *)calloc(K_x, sizeof(ST_int));
            if (!pre_collinear) { free(xtx_pre); output_rc = 920; goto ppml_cleanup; }
            if (pre_collinear) {
                detect_collinearity(xtx_pre, K_x, pre_collinear, verbose);
                ST_int n_coll = 0;
                for (k = 0; k < K_x; k++)
                    if (pre_collinear[k]) n_coll++;
                K_irls = K_x - n_coll;
                if (n_coll > 0 && K_irls > 0) {
                    irls_keep_idx = (ST_int *)malloc(K_irls * sizeof(ST_int));
                    if (!irls_keep_idx) { free(pre_collinear); free(xtx_pre); output_rc = 920; goto ppml_cleanup; }
                    if (irls_keep_idx) {
                        idx = 0;
                        for (k = 0; k < K_x; k++)
                            if (!pre_collinear[k])
                                irls_keep_idx[idx++] = k;
                    }
                    if (verbose) {
                        char msg[128];
                        snprintf(msg, sizeof(msg),
                            "(cppmlhdfe: %d variable(s) collinear with FE, excluded from IRLS)\n", n_coll);
                        SF_display(msg);
                    }
                }
                free(pre_collinear);
            }
            free(xtx_pre);
        }
    }

    /* Allocate compact IRLS arrays if collinearity detected */

    if (irls_keep_idx && K_irls < K_x) {
        xtx_k = (ST_double *)malloc(K_irls * K_irls * sizeof(ST_double));
        xty_k = (ST_double *)malloc(K_irls * sizeof(ST_double));
        beta_k = (ST_double *)malloc(K_irls * sizeof(ST_double));
        inv_xx_k = (ST_double *)malloc(K_irls * K_irls * sizeof(ST_double));
        if (!xtx_k || !xty_k || !beta_k || !inv_xx_k) { output_rc = 920; goto ppml_cleanup; }
    }

    /* ================================================================
     * STEP 9: IRLS LOOP
     * ================================================================ */
    ST_double deviance = 1e30, deviance_old;
    ST_int irls_iter;
    ST_int irls_converged = 0, ok_iterations=0, halving_count=0;
    memset(sep_m,0,(size_t)N*sizeof(int));
    ST_double log_eps_history[3] = {0.0, 0.0, 0.0};
    ST_int eps_history_count = 0;

    /* Adaptive CG tolerance matching ppmlhdfe:
     * Start at 1e-4 (fast early iterations), tighten toward 1e-9
     * (precise final iterations) based on relative deviance change. */
    ST_double target_inner_tol = fmax(1e-12, fmin(1e-9, 0.1 * tolerance));
    if (SF_scal_use("__cppmlhdfe_inner_tol", &val) == 0 && val > 0) target_inner_tol = val;
    ST_double start_inner_tol = options->start_tol;
    g_ppml_state->tolerance = options->heuristic ? fmax(start_inner_tol,target_inner_tol) : target_inner_tol;
    /* The reference may only stop after an iteration solved exactly, and it
     * decides that from its own lagged inner-tolerance schedule. The native
     * schedule tightens faster, so tracking only it stopped one iteration
     * earlier than the reference at default tolerances (with weights from a
     * less converged iterate). Shadow the reference schedule for the stopping
     * decision only; the native solves keep their own tolerance. Unlike the
     * reference, an iteration at the inner-tolerance floor counts as exact:
     * without the heuristic and with tol < 11*floor the reference never
     * leaves its fast solver and cannot stop. */
    ST_double ref_inner_tol = fmax(start_inner_tol, irls_tol);
    ST_double ref_alt_tol = DBL_MAX;

    for (irls_iter = 0; irls_iter < irls_maxiter; irls_iter++) {
        deviance_old = deviance;

        /* (a) Compute working depvar: z = eta + (y - mu)/mu
         * Save in z_orig for FE reconstruction; copy into aug_copy for CG. */
        #pragma omp parallel for schedule(static) if (N > 100000)
        for (i = 0; i < N; i++) {
            ST_double z = eta[i] + (y[i] - mu[i]) / mu[i];
            z_orig[i] = z;
            aug_copy[i] = options->exact_partial || !irls_iter ? z : aug_copy[i]+z-previous_working[i];
            previous_working[i]=z;
        }

        /* Reuse the previous residuals when inexact partialling is requested. */
        if (options->exact_partial || !irls_iter) for (k = 0; k < K_x; k++) {
            memcpy(&aug_copy[(size_t)(k + 1) * N], &X[(size_t)k * N], N * sizeof(ST_double));
        }

        /* (d) Partial out via CG solver using current IRLS weights */
        HDFE_SolveResult projection = ppml_partial(options, g_ppml_state, aug_copy, N, K_aug, num_threads);
        if (projection.status) {
            SF_error("cppmlhdfe: fixed-effect projection did not converge\n");
            output_rc = projection.status;
            goto ppml_cleanup;
        }

        /* (e) Solve WLS: beta = (X~'W X~)^{-1} X~'W z~ using partialled data */
        memset(irls_beta, 0, K_x * sizeof(ST_double));

        if (options->exact_solver || !options->remove_collinear) {
            output_rc=ppml_qrsolve(aug_copy+N,aug_copy,irls_w,N,K_x,irls_beta,NULL);
            if(output_rc)goto ppml_cleanup;
        } else if (irls_keep_idx && K_irls < K_x && xtx_k && xty_k && beta_k && inv_xx_k) {
            /* Compact solve: only non-collinear columns.
             * Use dd_real (double-double) precision for X'WX and X'Wz to match
             * ppmlhdfe/reghdfe's quadcross(). Small precision differences in
             * each IRLS iteration accumulate in the mu trajectory, which then
             * gets amplified through cond(X'WX) in the final VCE. */
            ST_int ki, kj;
            output_rc = ppml_normal_equations(aug_copy, N, irls_keep_idx, K_irls, irls_w, xtx_k, xty_k, NULL);
            if (output_rc) goto ppml_cleanup;

            memcpy(inv_xx_k, xtx_k, K_irls * K_irls * sizeof(ST_double));
            if (ctools_cholesky(inv_xx_k, K_irls) != 0) {
                break;
            }
            output_rc = ctools_stata_rc(ctools_invert_from_cholesky(inv_xx_k, K_irls, inv_xx_k));
            if (output_rc) goto ppml_cleanup;

            for (ki = 0; ki < K_irls; ki++) {
                beta_k[ki] = 0.0;
                for (kj = 0; kj < K_irls; kj++)
                    beta_k[ki] += inv_xx_k[ki * K_irls + kj] * xty_k[kj];
            }
            /* Expand to full beta */
            for (ki = 0; ki < K_irls; ki++)
                irls_beta[irls_keep_idx[ki]] = beta_k[ki];
        } else if (K_irls > 0) {
            /* Full solve: all columns (no collinearity detected).
             * Use dd_real for X'WX and X'Wz to match quadcross() precision. */
            output_rc = ppml_normal_equations(aug_copy, N, NULL, K_x, irls_w, xtx, xty, NULL);
            if (output_rc) goto ppml_cleanup;

            memcpy(inv_xx, xtx, K_x * K_x * sizeof(ST_double));
            if (ctools_cholesky(inv_xx, K_x) != 0) break;
            output_rc = ctools_stata_rc(ctools_invert_from_cholesky(inv_xx, K_x, inv_xx));
            if (output_rc) goto ppml_cleanup;

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
        ST_double max_eta_change = 0.0;
        /* X~ * beta, each row accumulated over the regressors in order */
        memset(xb_work, 0, (size_t)N * sizeof(ST_double));
        ctools_matvec_add(aug_copy + N, N, K_x, irls_beta, xb_work);
        int eta_not_finite = 0;
        #pragma omp parallel for schedule(static) reduction(max : max_eta_change) reduction(| : eta_not_finite) if (N > 100000)
        for (i = 0; i < N; i++) {
            /* FE component: original z minus partialled z */
            ST_double fe_part = z_orig[i] - aug_copy[i];
            ST_double z_hat = xb_work[i] + fe_part;
            ST_double off = (offset_compact != NULL) ? offset_compact[i] : 0.0;
            /* eta = X*beta + FE (without offset); z didn't include offset */
            if (!isfinite(z_hat)) { eta_not_finite |= 1; continue; }
            /* Zero outcomes with negligible means may have no finite FE
             * coefficient after a tolerance-limited separation check. They
             * must not prevent convergence of the identified parameters. */
            if (y[i] > 0 || mu[i] > 1e-12)
                max_eta_change = fmax(max_eta_change, fabs(z_hat - eta[i]));
            eta[i] = z_hat;
            mu[i] = exp(eta[i] + off);
            /* Clamp mu to prevent overflow/underflow */
            if (mu[i] > 1e18) mu[i] = 1e18;
            if (mu[i] < 1e-18) mu[i] = 1e-18;
        }
        if (eta_not_finite) output_rc = 430;

        if (output_rc) break;
        if(options->check_mu) {
            double min_eta=DBL_MAX;
            for(i=0;i<N;i++)if(y[i]>0)min_eta=fmin(min_eta,eta[i]+(offset_compact?offset_compact[i]:0));
            double cutoff=log(options->mu_tol)+fmin(min_eta+5,0);
            for(i=0;i<N;i++) {
                sep_m[i]|=y[i]==0 && eta[i]+(offset_compact?offset_compact[i]:0)<=cutoff;
                if(sep_m[i])mu[i]=1.4210854715202004e-14;
            }
        }

        /* (g) Compute deviance */
        deviance = ppml_compute_deviance(y, mu, w_user_compact, N) / user_weight_scale;

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
            /* Match the reference's decision to finish with an exact solve:
             * extrapolate log(eps) with a cubic through iterations 0..3,
             * where the intercept is zero. One-way demeaning is exact; for
             * multiple FE dimensions also require the target inner tolerance. */
            ST_double predicted_eps = DBL_MAX;
            if (eps_history_count == 3) {
                predicted_eps = exp(4.0 * log_eps_history[0] -
                                    6.0 * log_eps_history[1] +
                                    4.0 * log_eps_history[2]);
            }
            int final_solve = options->exact_solver || K_irls == 0 || halving_count > 0 ||
                              ref_inner_tol <= 11.0 * irls_tol ||
                              ref_inner_tol <= 1.1 * target_inner_tol ||
                              (options->heuristic && predicted_eps <= irls_tol);
            /* Deviance is locally quadratic; covariance depends on weights
             * to first order. Also require stable fitted log-means so an
             * early deviance plateau cannot leave inaccurate standard errors
             * in explicitly tight fits. Preserve default ppmlhdfe stopping
             * behavior; it can intentionally return less accurate weights. */
            if (rel_change < irls_tol && final_solve &&
                (irls_tol >= 1e-10 || max_eta_change <= 1e-8) &&
                (G == 1 || (g_ppml_state->tolerance <= target_inner_tol &&
                            ref_inner_tol <= 1.1 * target_inner_tol))) {
                ok_iterations++;
                if(ok_iterations>=options->min_ok) {
                    irls_converged = 1;
                    irls_iter++;
                    break;
                }
            }

            if(options->step_halving && deviance>deviance_old && halving_count<options->max_halving) {
                for(i=0;i<N;i++) {
                    double wu=w_user_compact?w_user_compact[i]:1;
                    double old=log(irls_w[i]*user_weight_scale/wu)-(offset_compact?offset_compact[i]:0);
                    eta[i]=options->halving_memory*old+(1-options->halving_memory)*eta[i];
                    if(halving_count && eta[i]<-10)eta[i]=-10;
                    mu[i]=exp(eta[i]+(offset_compact?offset_compact[i]:0));
                }
                halving_count++;ok_iterations=0;deviance=deviance_old;
            } else {
                halving_count=0;
                /* Reference schedule: tighten at least 10x once eps falls
                 * below the inner tolerance, using the previous iteration's
                 * eps for the larger jump; never below the target. */
                if (rel_change < ref_inner_tol) {
                    ref_inner_tol = fmax(fmin(0.1 * ref_inner_tol, ref_alt_tol), target_inner_tol);
                    ref_alt_tol = pow(10.0, -ceil(log10(1.0 / fmax(0.1 * rel_change, DBL_EPSILON))));
                }
                /* The reference partials inexactly by default; native exact
                 * partialling is a default, not a request, so only an explicit
                 * request takes the reference's exact-partial clamp. */
                if (options->ref_exact_partial && ref_inner_tol > irls_tol) ref_inner_tol = 0.1 * irls_tol;
                if (options->heuristic && predicted_eps <= irls_tol && irls_tol <= rel_change &&
                    ref_inner_tol >= 1.1 * target_inner_tol && ref_inner_tol <= 10.1 * target_inner_tol)
                    ref_inner_tol *= 0.1;
            }

            if (rel_change > 0.0 && isfinite(rel_change)) {
                if (eps_history_count < 3) log_eps_history[eps_history_count++] = log(rel_change);
                else {
                    log_eps_history[0] = log_eps_history[1];
                    log_eps_history[1] = log_eps_history[2];
                    log_eps_history[2] = log(rel_change);
                }
            }

            /* Adaptive CG tolerance: tighten as IRLS converges.
             * Matches ppmlhdfe's scheme: new_tol = max(target, eps * 0.01) */
            if (options->heuristic && rel_change < 0.1) {
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
        #pragma omp parallel for schedule(static) if (N > 100000)
        for (i = 0; i < N; i++) {
            ST_double w_u = (w_user_compact != NULL) ? w_user_compact[i] : 1.0;
            irls_w[i] = mu[i] * (w_u / user_weight_scale);
            if (irls_w[i] < DBL_MIN) irls_w[i] = DBL_MIN;
        }
        ppml_update_fe_weights(g_ppml_state, irls_w, N);


    }

    if (!irls_converged && !output_rc) {
        SF_error("cppmlhdfe: IRLS did not converge\n");
        output_rc = 430;
    }
    if (output_rc) goto ppml_cleanup;

    if(options->check_mu) {
        int nsep=0;for(i=0;i<N;i++)nsep+=sep_m[i];
        if(nsep) {
            sample_var_idx=(int)ppml_option("sample_idx",0);
            for(i=0;i<N && !output_rc;i++) {
                if(!sep_m[i])output_rc=SF_vstore(sample_var_idx,(int)obs_map[i],1);
            }
            if(!output_rc)output_rc=ctools_scal_save("__cppmlhdfe_refit",3);
            if(!output_rc)output_rc=ctools_scal_save("__cppmlhdfe_num_singletons",num_singletons);
            if(!output_rc)output_rc=ctools_scal_save("__cppmlhdfe_num_separated",num_separated+nsep);
            goto ppml_cleanup;
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
        if(options->exact_partial || !irls_iter) {
            for (k = 0; k < K_x; k++)
                memcpy(&aug_copy[(size_t)(k + 1) * N], &X[(size_t)k * N], N * sizeof(ST_double));
        }
        HDFE_SolveResult projection = ppml_partial(options, g_ppml_state, aug_copy, N, K_aug, num_threads);
        if (projection.status) {
            SF_error("cppmlhdfe: fixed-effect projection did not converge\n");
            output_rc = projection.status;
            goto ppml_cleanup;
        }
    }

    /* Collinearity detection */
    is_collinear = (ST_int *)calloc(K_x, sizeof(ST_int));
    if (!is_collinear) { output_rc = 920; goto ppml_cleanup; }
    ST_int num_collinear = 0;
    ST_int K_keep;

    if (is_collinear) {
        ST_double *xtx_final = (ST_double *)malloc(K_x * K_x * sizeof(ST_double));
        if (!xtx_final) { output_rc = 920; goto ppml_cleanup; }
        if (xtx_final) {
            ctools_matmul_ata(aug_copy + N, N, K_x, xtx_final);
            detect_collinearity(xtx_final, K_x, is_collinear, verbose);
            free(xtx_final);
        }
        for (k = 0; k < K_x; k++)
            if (is_collinear[k]) num_collinear++;
    }
    K_keep = K_x - num_collinear;

    /* Store collinearity flags */
    if ((output_rc = ctools_scal_save("__cppmlhdfe_num_collinear", (ST_double)num_collinear)))
        goto ppml_cleanup;
    for (k = 0; k < K_x; k++) {
        snprintf(scalar_name, sizeof(scalar_name), "__cppmlhdfe_collinear_%d", k + 1);
        if ((output_rc = ctools_scal_save(scalar_name, (ST_double)(is_collinear ? is_collinear[k] : 0))))
            goto ppml_cleanup;
    }

    /* Build non-collinear index */
    keep_idx = (ST_int *)malloc(K_keep * sizeof(ST_int));
    if (!keep_idx && K_keep > 0) { output_rc = 920; goto ppml_cleanup; }
    if (keep_idx && is_collinear) {
        idx = 0;
        for (k = 0; k < K_x; k++) {
            if (!is_collinear[k])
                keep_idx[idx++] = k;
        }
    }

    /* Final solve: beta and inv(X'WX)
     * Use converged IRLS weights for the final normal equations so stored
     * coefficients match the terminal PPML optimum. */
    beta_final = (ST_double *)calloc(K_keep, sizeof(ST_double));
    if (!beta_final && K_keep > 0) { output_rc = 920; goto ppml_cleanup; }

    ST_double N_ref = (ST_double)N;
    if (has_weights && weight_type == 2 && w_user_compact) {
        N_ref = 0.0;
        for (i = 0; i < N; i++) N_ref += w_user_compact[i];
    }
    ST_double sum_irls_w = 0.0;
    for (i = 0; i < N; i++) sum_irls_w += irls_w[i];
    ST_double w_scale = (sum_irls_w > 0.0) ? (N_ref / sum_irls_w) : 1.0;
    w_reg = (ST_double *)malloc(N * sizeof(ST_double));
    if (!w_reg) { output_rc = 920; goto ppml_cleanup; }
    if (w_reg) {
        for (i = 0; i < N; i++) w_reg[i] = irls_w[i] * w_scale;
    }

    if (K_keep > 0 && keep_idx && beta_final) {
        xtx_keep = (ST_double *)calloc(K_keep * K_keep, sizeof(ST_double));
        xty_keep = (ST_double *)calloc(K_keep, sizeof(ST_double));
        inv_xx_final = (ST_double *)malloc(K_keep * K_keep * sizeof(ST_double));
        V_keep = (ST_double *)calloc(K_keep * K_keep, sizeof(ST_double));

        if (!xtx_keep || !xty_keep || !inv_xx_final || !V_keep) { output_rc = 920; goto ppml_cleanup; }
        if (xtx_keep && xty_keep && inv_xx_final && V_keep) {

            /* X'WX and X'Wz in double-double precision */
            /* the unrounded X'Wz sums are kept for the refinement below */
            xty_dd = (dd_real *)malloc((size_t)K_keep * sizeof(dd_real));
            if (!xty_dd) { output_rc = 920; goto ppml_cleanup; }
            output_rc = ppml_normal_equations(aug_copy, N, keep_idx, K_keep, w_reg, xtx_keep, xty_keep, xty_dd);
            if (output_rc) goto ppml_cleanup;

            /* Solve with iterative refinement for maximum beta precision. */
            memcpy(inv_xx_final, xtx_keep, K_keep * K_keep * sizeof(ST_double));
            if (ctools_cholesky(inv_xx_final, K_keep) != 0) { output_rc = 498; goto ppml_cleanup; }
            {
                output_rc = ctools_stata_rc(ctools_invert_from_cholesky(inv_xx_final, K_keep, inv_xx_final));
                if (output_rc) goto ppml_cleanup;
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
                        /* X'Wz row i, as accumulated for the normal equations */
                        dd_real acc = xty_dd[i];
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

    /* Compute IRLS-weighted means of original X columns (before partialling).
     * These means are added back to the partialled X for VCE computation.
     * Uses IRLS weights = mean(x, HDFE.weight) matching ppmlhdfe. */
    ST_int K_with_cons = K_keep + 1;
    means_x = (ST_double *)calloc(K_keep, sizeof(ST_double));
    if (!means_x && K_keep > 0) { output_rc = 920; goto ppml_cleanup; }
    if (means_x && K_keep > 0) {
        ST_double sum_w_eff = 0.0;
        for (i = 0; i < N; i++) sum_w_eff += irls_w[i];
        for (k = 0; k < K_keep; k++) {
            ST_double wmean = 0.0;
            const ST_double *xk = &X[(size_t)keep_idx[k] * N];
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
    inv_xx_ext = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));
    if (!inv_xx_ext) { output_rc = 920; goto ppml_cleanup; }
    /* N_ref = 1'W1 for the rescaled converged IRLS weights */
    if (K_keep == 0 || (inv_xx_final && means_x))
        ctools_inverse_add_constant(inv_xx_final, K_keep, means_x, N_ref, inv_xx_ext);

    /* Build X_eff with means added back + constant column (K_with_cons columns) */

    {
        X_eff = (ST_double *)malloc((size_t)N * K_with_cons * sizeof(ST_double));
        if (!X_eff) { output_rc = 920; goto ppml_cleanup; }
        if (X_eff) {
            for (k = 0; k < K_keep; k++) {
                const ST_double *xp = &aug_copy[(size_t)(keep_idx[k] + 1) * N];
                ST_double mk = means_x ? means_x[k] : 0.0;
                for (i = 0; i < N; i++) {
                    X_eff[(size_t)k * N + i] = xp[i] + mk;
                }
            }
            /* Constant column */
            for (i = 0; i < N; i++) {
                X_eff[(size_t)K_keep * N + i] = 1.0;
            }
        }
    }

    /* Compute final OLS residual from partialled data. */
    vce_resid = (ST_double *)malloc(N * sizeof(ST_double));
    if (!vce_resid) { output_rc = 920; goto ppml_cleanup; }
    if (vce_resid) {
        if (beta_final && K_keep > 0) {
            for (i = 0; i < N; i++) {
                ST_double xb = 0.0;
                for (k = 0; k < K_keep; k++) {
                    xb += aug_copy[(size_t)(keep_idx[k] + 1) * N + i] * beta_final[k];
                }
                vce_resid[i] = aug_copy[i] - xb;
            }
        }
        else {
            for (i = 0; i < N; i++)
                vce_resid[i] = aug_copy[i];
        }
    }

    /* Count each marginal cluster and detect FE nesting on the retained sample.
     * Subset indices are bitmasks; powers of two are marginal dimensions. */
    ST_int nested[10] = {0};
    if (vcetype == 2) {
        cluster_ids = malloc(N * sizeof(ST_int));
        if (!cluster_ids) { output_rc = 920; goto ppml_cleanup; }
        num_clusters = N;
        ST_int dim = 0;
        for (ST_int subset = 1; subset <= cluster_terms; subset <<= 1) {
            ST_int ncl;
            if (ctools_numeric_to_cluster_ids(
                    cluster_raw_values + (size_t)(subset - 1) * N_orig,
                    N, cluster_ids, &ncl)) {
                output_rc = 920; goto ppml_cleanup;
            }
            if (ncl < num_clusters) num_clusters = ncl;
            snprintf(scalar_name, sizeof(scalar_name), "__cppmlhdfe_N_clust_%d", ++dim);
            if ((output_rc = ctools_scal_save(scalar_name, ncl))) goto ppml_cleanup;
            for (g = 0; g < G; g++) {
                if (slope_only[g] || !(dof_flags & 4)) continue;
                ST_int is_nested = ctools_fe_nested_in_cluster(
                    g_ppml_state->factors[g].levels,
                    g_ppml_state->factors[g].num_levels, cluster_ids, N);
                if (is_nested < 0) { output_rc = 920; goto ppml_cleanup; }
                nested[g] |= is_nested;
            }
        }
        if (num_clusters < 2) {
            SF_error("cppmlhdfe: each cluster dimension requires at least two retained clusters\n");
            output_rc = 459; goto ppml_cleanup;
        }
    }
    for (g = 0; g < G; g++) {
        snprintf(scalar_name, sizeof(scalar_name), "__cppmlhdfe_fe_nested_%d", g + 1);
        if ((output_rc = ctools_scal_save(scalar_name, nested[g]))) goto ppml_cleanup;
    }

    /* Match reghdfe's pairwise mobility and cluster-nesting accounting.
     * These degrees of freedom describe the model; PPML's asymptotic
     * sandwich uses only N/(N-1) or G/(G-1), not residual-DF scaling. */
    ST_int redundant[10] = {0}, exact[10] = {0};
    ST_int df_initial = 0, df_nested = 0, df_redundant = 0, first_intercept = -1, pairs = 0;
    for (g = 0; g < G; g++) if (nested[g]) df_nested += factors[g].num_levels;
    for (g = 0; g < G; g++) {
        FE_Factor *fg = &g_ppml_state->factors[g];
        df_initial += fg->num_levels;
        if (nested[g]) { redundant[g] = fg->num_levels; exact[g] = 1; continue; }
        if (slope_only[g]) {
            if (dof_flags & 8) {
                ST_double *lo = malloc((size_t)fg->num_levels*sizeof(*lo));
                ST_double *hi = malloc((size_t)fg->num_levels*sizeof(*hi));
                if (!lo || !hi) { free(lo); free(hi); output_rc = 920; goto ppml_cleanup; }
                for (i = 0; i < fg->num_levels; i++) { lo[i] = DBL_MAX; hi[i] = -DBL_MAX; }
                for (i = 0; i < N; i++) {
                    ST_int l = fg->levels[i]-1;
                    lo[l] = fmin(lo[l], fg->slope[i]); hi[l] = fmax(hi[l], fg->slope[i]);
                }
                ST_int has_int = 0;
                for (j = 0; j < G; j++) {
                    FE_Factor *fj = &g_ppml_state->factors[j];
                    if (!slope_only[j] && fg->num_levels == fj->num_levels &&
                        !memcmp(fg->levels, fj->levels, (size_t)N*sizeof(ST_int))) has_int = 1;
                }
                for (i = 0; i < fg->num_levels; i++)
                    redundant[g] += has_int ? hi[i]-lo[i] <= 1e-6 : fabs(lo[i])+fabs(hi[i]) <= 1e-6;
                free(lo); free(hi);
            }
            continue;
        }
        redundant[g] = first_intercept >= 0 || df_nested > 0;
        if (first_intercept < 0) { first_intercept = g; exact[g] = 1; }
        if (!(dof_flags & 3)) continue;
        for (j = 0; j < g; j++) {
            if (slope_only[j] || nested[j]) continue;
            if ((dof_flags & 2) && pairs) break;
            ST_int groups = count_connected_components(g_ppml_state->factors[j].levels,
                fg->levels, N, g_ppml_state->factors[j].num_levels, fg->num_levels);
            if (groups < 0) { output_rc = 920; goto ppml_cleanup; }
            if (!pairs) exact[g] = 1;
            pairs++;
            if (groups > redundant[g]) redundant[g] = groups;
        }
    }
    for (g = 0; g < G; g++) {
        df_redundant += redundant[g];
        snprintf(scalar_name, sizeof(scalar_name), "__cppmlhdfe_fe_redundant_%d", g+1);
        if ((output_rc = ctools_scal_save(scalar_name, redundant[g]))) goto ppml_cleanup;
        snprintf(scalar_name, sizeof(scalar_name), "__cppmlhdfe_fe_exact_%d", g+1);
        if ((output_rc = ctools_scal_save(scalar_name, exact[g]))) goto ppml_cleanup;
    }
    df_a = df_initial - df_redundant;
    if ((output_rc = ctools_scal_save("__cppmlhdfe_df_initial", df_initial))) goto ppml_cleanup;
    if ((output_rc = ctools_scal_save("__cppmlhdfe_df_nested", df_nested))) goto ppml_cleanup;
    if ((output_rc = ctools_scal_save("__cppmlhdfe_df_redundant", df_redundant))) goto ppml_cleanup;

    /* Compute VCE: PPML sandwich estimator matching creghdfe/reghdfe.
     * Uses extended (K_keep+1)x(K_keep+1) system with constant column,
     * then extracts K_keep×K_keep submatrix for reporting.
     *
     * For fweight: D = extended H^{-1} (raw), resid = (y-mu), weights = fw
     * For others: D extended & normalized, resid = OLS_resid, weights = irls_w as aweight */
  /* (K_with_cons)x(K_with_cons) full VCE */
    if (vce_resid && X_eff && inv_xx_ext) {
        V_ext = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));
        if (!V_ext) { output_rc = 920; goto ppml_cleanup; }

        if (vcetype == 0) {
            /* Unadjusted VCE: V = H^{-1} (Poisson dispersion = 1) */
            for (i = 0; i < K_keep; i++)
                for (j = 0; j < K_keep; j++)
                    V_keep[i * K_keep + j] = inv_xx_final[i * K_keep + j];

        } else if (V_ext && ((vcetype == 1) || (vcetype == 2 && cluster_ids))) {

            ctools_vce_data d;
            /* For a single sandwich, the centered slope system is equivalent
             * and avoids cancellation from large regressor means. Multiway
             * PSD correction needs the full system including the constant. */
            ST_int vce_K = K_with_cons;
            if (cluster_terms <= 1) {
                for (k = 0; k < K_keep; k++)
                    memcpy(X_eff + (size_t)k * N,
                           aug_copy + (size_t)(keep_idx[k] + 1) * N,
                           N * sizeof(ST_double));
            }
            if (cluster_terms <= 1) {
                for (i = 0; i < K_keep; i++) {
                    inv_xx_ext[i*K_with_cons+K_keep] = 0;
                    inv_xx_ext[K_keep*K_with_cons+i] = 0;
                }
                inv_xx_ext[K_keep*K_with_cons+K_keep] = 1.0/N_ref;
            }
            d.X_eff = X_eff;
            d.D = inv_xx_ext;
            ST_double *vce_resid_work = vce_resid;
            if (has_weights && weight_type == 2 && w_user_compact) {
                vce_resid_work = (ST_double *)malloc(N * sizeof(ST_double));
                if (!vce_resid_work) { output_rc = 920; goto ppml_cleanup; }
                memcpy(vce_resid_work, vce_resid, N * sizeof(ST_double));
                for (i = 0; i < N; i++) {
                    ST_double fw = w_user_compact[i];
                    ST_double wr = w_reg ? w_reg[i] : irls_w[i];
                    if (fw > 0.0 && vce_resid_work) {
                        vce_resid_work[i] *= wr / fw;
                    }
                }
                d.weights = w_user_compact;
                d.weight_type = 2;
                d.normalize_weights = 0;
            }
            else {
                d.weights = w_reg ? w_reg : irls_w;
                d.weight_type = 1;
                d.normalize_weights = 1;
            }
            d.resid = vce_resid_work ? vce_resid_work : vce_resid;
            d.N = N;
            d.K = vce_K;

            if (vcetype == 1) {
                ST_double dof_adj = N_ref / (N_ref - 1.0);
                output_rc = ctools_vce_robust(&d, dof_adj, V_ext);
            } else {
                ST_double dof_adj = (ST_double)num_clusters / (num_clusters - 1);
                ST_double *V_term = calloc((size_t)K_with_cons * K_with_cons, sizeof(ST_double));
                if (!V_term) output_rc = 920;
                for (ST_int subset = 1; subset <= cluster_terms && !output_rc; subset++) {
                    ST_int ncl, bits = 0;
                    for (ST_int bit = subset; bit; bit >>= 1) bits += bit & 1;
                    if (ctools_numeric_to_cluster_ids(
                            cluster_raw_values + (size_t)(subset - 1) * N_orig,
                            N, cluster_ids, &ncl)) {
                        output_rc = 920; break;
                    }
                    output_rc = ctools_vce_cluster(&d, cluster_ids, ncl, dof_adj, V_term);
                    ST_double sign = bits % 2 ? 1.0 : -1.0;
                    for (i = 0; i < K_with_cons * K_with_cons && !output_rc; i++)
                        V_ext[i] += sign * V_term[i];
                }
                free(V_term);
            }
            if (vce_resid_work && vce_resid_work != vce_resid) free(vce_resid_work);

            if (output_rc) {
                SF_error("cppmlhdfe: covariance calculation failed\n");
                goto ppml_cleanup;
            }

            /* Centered-system sandwich is numerically stable even when a
             * regressor has a huge mean. Transform only the constant row and
             * column back to the reported parameterization. */
            if (cluster_terms <= 1) {
                ST_double corner = V_ext[K_keep*K_with_cons+K_keep];
                for (i = 0; i < K_keep; i++) {
                    corner -= 2*means_x[i]*V_ext[i*K_with_cons+K_keep];
                    for (j = 0; j < K_keep; j++)
                        corner += means_x[i]*means_x[j]*V_ext[i*K_with_cons+j];
                }
                for (i = 0; i < K_keep; i++) {
                    ST_double side = V_ext[i*K_with_cons+K_keep];
                    for (j = 0; j < K_keep; j++) side -= V_ext[i*K_with_cons+j]*means_x[j];
                    V_ext[i*K_with_cons+K_keep] = V_ext[K_keep*K_with_cons+i] = side;
                }
                V_ext[K_keep*K_with_cons+K_keep] = corner;
            }
            for (i = 0; i < K_with_cons; i++) {
                if ((output_rc = ctools_mat_store("__cppmlhdfe_scales", 1, i+1,
                                                  i < K_keep ? stdev_x[keep_idx[i]] : 1.0)))
                    goto ppml_cleanup;
                for (j = 0; j < K_with_cons; j++)
                    if ((output_rc = ctools_mat_store("__cppmlhdfe_V_full", i+1, j+1, V_ext[i*K_with_cons+j])))
                        goto ppml_cleanup;
            }

            /* Extract K_keep×K_keep submatrix (dropping constant row/col) */
            for (i = 0; i < K_keep; i++)
                for (j = 0; j < K_keep; j++)
                    V_keep[i * K_keep + j] = V_ext[i * vce_K + j];
        }
    }
    free(V_ext); V_ext = NULL;

    t_vce = get_time_sec();

    /* ================================================================
     * STEP 12: Un-standardize and store results
     * ================================================================ */

    /* Weighted means as sum(w*v)/sum(w) with double-double sums, as the
     * reference's quad-precision mean(v, w). Summing w/sum_w*v instead
     * does not reproduce a constant v, so d was not centered exactly. The
     * constant follows the reference order: mean(log mu) - xbar*b
     * - mean(offset) + log(stdev_y). */
    ST_double constant, mean_d;
    {
        dd_real sw = {0.0, 0.0}, slog = {0.0, 0.0}, soff = {0.0, 0.0}, sd = {0.0, 0.0};
        for (i = 0; i < N; i++) {
            ST_double w = irls_w[i];
            ST_double d = z_orig[i] - vce_resid[i];
            for (k = 0; k < K_keep; k++) d -= X[(size_t)keep_idx[k]*N+i]*beta_final[k];
            sw = dd_add_d(sw, w);
            slog = dd_add_d(slog, w * log(mu[i]));
            if (offset_compact) soff = dd_add_d(soff, w * offset_compact[i]);
            sd = dd_add_d(sd, w * d);
        }
        ST_double wsum = sw.hi + sw.lo;
        constant = (slog.hi + slog.lo) / wsum;
        mean_d = (sd.hi + sd.lo) / wsum;
        for (k = 0; k < K_keep; k++) constant -= means_x[k]*beta_final[k];
        if (offset_compact) constant -= (soff.hi + soff.lo) / wsum;
        constant += log(stdev_y);
    }
    if ((output_rc = ctools_scal_save("__cppmlhdfe_cons", constant))) goto ppml_cleanup;
    ST_int d_var_idx = 0;
    if (SF_scal_use("__cppmlhdfe_d_idx", &val) == 0) d_var_idx = (ST_int)val;
    if (d_var_idx > 0) {
        for (i = 0; i < N && !output_rc; i++) {
            ST_double d = z_orig[i] - vce_resid[i] - mean_d;
            for (k = 0; k < K_keep; k++) d -= X[(size_t)keep_idx[k]*N+i]*beta_final[k];
            output_rc = SF_vstore(d_var_idx, (ST_int)obs_map[i], d);
        }
    }

    for(i=0;i<N;i++) {
        separation_certificate[i]=z_orig[i]-vce_resid[i]-mean_d;
        for(k=0;k<K_keep;k++)separation_certificate[i]-=X[(size_t)keep_idx[k]*N+i]*beta_final[k];
    }
    /* A failed d() store above must not be replaced by this call's status */
    if(!output_rc)output_rc=ppml_save_effects(g_ppml_state,separation_certificate,obs_map);
    if(output_rc)goto ppml_cleanup;

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

    if (SF_scal_use("__cppmlhdfe_sample_idx", &val) == 0)
        sample_var_idx = (ST_int)val;
    if (sample_var_idx > 0) {
        for (i = 0; i < N && !output_rc; i++)
            output_rc = SF_vstore(sample_var_idx, (ST_int)obs_map[i], 1.0);
    }

    const cppmlhdfe_result result = {
        .N = N,
        .weight_type = weight_type,
        .has_weights = has_weights,
        .num_singletons = num_singletons,
        .num_separated = num_separated,
        .K_keep = K_keep,
        .df_a = df_a,
        .mobility_groups = mobility_groups,
        .irls_iter = irls_iter,
        .irls_converged = irls_converged,
        .G = G,
        .cluster_terms = cluster_terms,
        .vcetype = vcetype,
        .num_clusters = num_clusters,
        .deviance = deviance,
        .stdev_y = stdev_y,
        .user_weight_scale = user_weight_scale,
        .ll = ll,
        .ll_0 = ll_0,
        .w_user_compact = w_user_compact,
        .beta_final = beta_final,
        .V_keep = V_keep,
        .factors = factors,
    };
    if (!output_rc) output_rc = cppmlhdfe_post_results(&result);
    if (output_rc) goto ppml_cleanup;

    /* Timing */
    if ((output_rc = ctools_scal_save("_cppmlhdfe_time_load", t_load - t_start))) goto ppml_cleanup;
    if ((output_rc = ctools_scal_save("_cppmlhdfe_time_remap", t_remap - t_load))) goto ppml_cleanup;
    if ((output_rc = ctools_scal_save("_cppmlhdfe_time_singleton", t_singleton - t_remap))) goto ppml_cleanup;
    if ((output_rc = ctools_scal_save("_cppmlhdfe_time_dof", t_dof - t_singleton))) goto ppml_cleanup;
    if ((output_rc = ctools_scal_save("_cppmlhdfe_time_irls", t_irls - t_dof))) goto ppml_cleanup;
    if ((output_rc = ctools_scal_save("_cppmlhdfe_time_vce", t_vce - t_irls))) goto ppml_cleanup;
    if ((output_rc = ctools_scal_save("_cppmlhdfe_time_total", t_vce - t_start))) goto ppml_cleanup;
    CTOOLS_SAVE_THREAD_INFO("_cppmlhdfe");

    /* ================================================================
     * Cleanup
     * ================================================================ */
ppml_cleanup:
    free(V_ext);
    free(separation_certificate);
    free(previous_working);
    free(cluster_raw_values);
    free(mu); free(eta);
    free(z_orig); free(aug_copy);
    free(xtx); free(xty); free(irls_beta); free(inv_xx); free(sep_m); free(xb_work);
    free(xty_dd);
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
    if (w_reg) free(w_reg);
    if (is_collinear) free(is_collinear);
    if (keep_idx) free(keep_idx);
    if (X_eff) free(X_eff);
    if (means_x) free(means_x);
    if (inv_xx_ext) free(inv_xx_ext);
    if (cluster_ids) free(cluster_ids);
    if (stdev_x) free(stdev_x);
    if (w_user_compact) free(w_user_compact);
    if (offset_compact) free(offset_compact);

    if (factors) {
        for (g = 0; g < G; g++) {
            free(factors[g].levels);
            free(factors[g].counts);
        }
        free(factors);
    }
    free(mask);
    free(data); free(data_compact);
    free(weights); free(offset_arr);
    ctools_aligned_free(obs_map);
    ctools_filtered_data_free(&filtered);

    /* State owns irls_w; cleanup releases it together with the FE buffers. */
    cleanup_ppml_state();

    return output_rc;
}

/* ========================================================================
 * Helper: Update FE weighted counts from IRLS weights
 * ======================================================================== */


/* ========================================================================
 * Helper: Compute Poisson deviance
 * D = 2 * sum[y*log(y/mu) - (y - mu)]  (y*log(y/mu) = 0 when y = 0)
 * ======================================================================== */


/* ========================================================================
 * Helper: Compute Poisson log-likelihood
 * ll = sum[y*log(mu) - mu - log(y!)]
 * ======================================================================== */

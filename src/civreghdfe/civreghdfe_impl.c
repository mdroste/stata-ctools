/*
    civreghdfe_impl.c
    Instrumental Variables Regression with High-Dimensional Fixed Effects

    Implements 2SLS/IV regression with HDFE absorption. Reuses shared
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
#include "civreghdfe_result.h"
#include "civreghdfe_input.h"

#include "../ctools_config.h"
#include "../ctools_runtime.h"
#include "../ctools_ols.h"
#include "../ctools_matrix.h"
#include "../ctools_types.h"  /* For ctools_data_load */
#include "../ctools_spi.h"  /* Error-checking SPI wrappers */
#include "civreghdfe_matrix.h"
#include "civreghdfe_estimate.h"

#include "civreghdfe_vce.h"
#include "civreghdfe_tests.h"
#include "../ctools_hdfe_utils.h"
#include "../hdfe/hdfe_sample.h"
#include "../ctools_unroll.h"

static HDFE_State *iv_state = NULL;

/* Weighted sum of squares sum(w x^2); w may be NULL */
static ST_double civ_wss(const ST_double *x, const ST_double *w, ST_int N)
{
    ST_double s = 0.0;
    if (w) { for (ST_int i = 0; i < N; i++) s += w[i] * x[i] * x[i]; }
    else   { for (ST_int i = 0; i < N; i++) s += x[i] * x[i]; }
    return s;
}

/* Weighted sum of squared deviations from the weighted mean; w may be NULL */
static ST_double civ_wss_dev(const ST_double *x, const ST_double *w, ST_int N)
{
    ST_double sw = 0.0, m = 0.0, s = 0.0;
    for (ST_int i = 0; i < N; i++) {
        ST_double wi = w ? w[i] : 1.0;
        sw += wi;
        m += wi * x[i];
    }
    if (sw <= 0.0) return 0.0;
    m /= sw;
    for (ST_int i = 0; i < N; i++) {
        ST_double d = x[i] - m;
        s += (w ? w[i] : 1.0) * d * d;
    }
    return s;
}

/* Clear keep[i] for every row whose value in v is missing (any missing code,
   or NaN). Same rows as SF_is_missing() for values read from Stata. */
static void civ_drop_missing(const ST_double *v, ST_int i0, ST_int i1, ST_double missval,
                             ST_int *keep)
{
    for (ST_int i = i0; i < i1; i++) {
        if (!(v[i] < missval)) keep[i] = 0;
    }
}

/* dst[j] = src[rows[j]], j < N: the retained rows of one column, in order */
static void civ_gather(ST_double *dst, const ST_double *src, const ST_int *rows, ST_int N)
{
    for (ST_int j = 0; j < N; j++) dst[j] = src[rows[j]];
}

/*
    Read the 1 x n matrix of 0/1 column flags that the ado passes for
    partial(), orthog(), endogtest() and redundant() (one flag per expanded
    plugin column). Returns the number of flagged columns, or -1 if the
    matrix is missing or has the wrong shape.
*/
static ST_int read_flag_matrix(const char *name, ST_int n, ST_int *flags)
{
    char mname[64];
    snprintf(mname, sizeof(mname), "%s", name);
    if (n <= 0 || SF_row(mname) != 1 || SF_col(mname) != n) return -1;
    ST_int count = 0;
    for (ST_int k = 0; k < n; k++) {
        ST_double v = 0.0;
        if (SF_mat_el(mname, 1, k + 1, &v) != 0) return -1;
        flags[k] = (v != 0.0);
        count += flags[k];
    }
    return count;
}

/*
    Map the ado's flags for orthog()/endogtest()/redundant() (matrix `name`,
    one flag per column of the original list) onto the columns that survived
    the collinearity checks (dropped[j] != 0 marks a removed column; NULL =
    none removed). Writes the distinct 1-based positions among the surviving
    columns to out[0..n_flag-1]. A flagged column that was dropped cannot be
    tested (ivreg2 rejects it as well): its option code and position are
    saved for the ado's message and 198 is returned.
*/
static ST_retcode map_test_cols(const char *name, const char *opt, ST_int opt_code,
                                ST_int n_flag, ST_int n_orig, const ST_int *dropped,
                                ST_int *out)
{
    char msg[200];
    if (n_flag <= 0) return 0;
    ST_int *flags = (ST_int *)calloc((size_t)(n_orig > 0 ? n_orig : 1), sizeof(ST_int));
    if (!flags) {
        SF_error("civreghdfe: Memory allocation failed for test flags\n");
        return 920;
    }
    if (read_flag_matrix(name, n_orig, flags) != n_flag) {
        snprintf(msg, sizeof(msg), "civreghdfe: %s() column flags do not match the variable list\n", opt);
        SF_error(msg);
        free(flags);
        return 198;
    }
    ST_int pos = 0, n = 0;
    for (ST_int j = 0; j < n_orig; j++) {
        ST_int keep = dropped ? (dropped[j] == 0) : 1;
        pos += keep;
        if (!flags[j]) continue;
        if (!keep) {
            snprintf(msg, sizeof(msg),
                     "civreghdfe: a variable in %s() was dropped as collinear or absorbed by the fixed effects\n", opt);
            SF_error(msg);
            ST_retcode save_rc = ctools_scal_save("__civreghdfe_dropped_opt", (ST_double)opt_code);
            if (!save_rc) save_rc = ctools_scal_save("__civreghdfe_dropped_col", (ST_double)(j + 1));
            free(flags);
            return save_rc ? save_rc : 198;
        }
        out[n++] = pos;
    }
    free(flags);
    return 0;
}

/*
    FWL step for partial(). data (N x total_cols, column-major) is already
    FE-demeaned; is_partial[k] flags X_exog column exog_start + k, whose copy
    in Z is zexog_start + k. Every other column is residualized on the
    partial columns P by W-weighted least squares, (P'WP)^-1 P'Wx with the
    estimator's weights: ivreg2's s_partial after ivreghdfe's HDFE
    partial_out. Partial columns absorbed by the fixed effects (weighted SS
    <= absorb_tol times their pre-partialling SS in ss_exog) or collinear
    with earlier partial columns are left out of P and marked 2, as ivreg2
    drops collinear partial() variables; the rest keep flag 1.
*/
static ST_retcode fwl_partial_out(ST_double *data, ST_int N, ST_int total_cols,
                                  ST_int exog_start, ST_int zexog_start, ST_int K_exog,
                                  ST_int *is_partial, const ST_double *w,
                                  const ST_double *ss_exog, ST_double absorb_tol)
{
    ST_int n_cand = 0, n_keep = 0;
    for (ST_int k = 0; k < K_exog; k++) if (is_partial[k]) n_cand++;
    if (n_cand == 0) return 0;

    ST_int *cols = (ST_int *)malloc((size_t)n_cand * sizeof(ST_int));
    ST_int *coll = (ST_int *)calloc((size_t)n_cand, sizeof(ST_int));
    ST_double *P = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)n_cand, sizeof(ST_double));
    ST_double *PtP = (ST_double *)ctools_safe_calloc3((size_t)n_cand, (size_t)n_cand, sizeof(ST_double));
    ST_double *PtP_inv = (ST_double *)ctools_safe_calloc3((size_t)n_cand, (size_t)n_cand, sizeof(ST_double));
    ST_double *PtWX = (ST_double *)ctools_safe_calloc3((size_t)n_cand, (size_t)total_cols, sizeof(ST_double));
    ST_double *coef = (ST_double *)ctools_safe_calloc3((size_t)n_cand, (size_t)total_cols, sizeof(ST_double));
    ST_retcode rc = 0;

    if (!cols || !coll || !P || !PtP || !PtP_inv || !PtWX || !coef) {
        SF_error("civreghdfe: Memory allocation failed for FWL workspace\n");
        rc = 920;
        goto fwl_done;
    }

    /* P = partial columns that were not absorbed by the fixed effects */
    for (ST_int k = 0; k < K_exog; k++) {
        if (!is_partial[k]) continue;
        const ST_double *x = data + (size_t)(exog_start + k) * N;
        ST_double ss = civ_wss(x, w, N);
        if (ss < 1e-30 || ss <= absorb_tol * ss_exog[k]) {
            is_partial[k] = 2;
            continue;
        }
        memcpy(P + (size_t)n_keep * N, x, (size_t)N * sizeof(ST_double));
        cols[n_keep++] = k;
    }

    /* Drop partial columns collinear with earlier ones (column order, as
       ivreg2's _rmcollright), using the equilibrated P'WP */
    if (n_keep > 1) {
        ctools_matmul_atdb(P, P, w, N, n_keep, n_keep, PtP);
        if (detect_collinearity_scaled(PtP, n_keep, coll, 0) < 0) {
            SF_error("civreghdfe: Collinearity detection failed for partial()\n");
            rc = 920;
            goto fwl_done;
        }
        ST_int m = 0;
        for (ST_int a = 0; a < n_keep; a++) {
            if (coll[a]) {
                is_partial[cols[a]] = 2;
                continue;
            }
            if (m != a) {
                memcpy(P + (size_t)m * N, P + (size_t)a * N, (size_t)N * sizeof(ST_double));
                cols[m] = cols[a];
            }
            m++;
        }
        n_keep = m;
    }
    if (n_keep == 0) goto fwl_done;

    /* coef = (P'WP)^-1 P'W data for every column */
    ctools_matmul_atdb(P, P, w, N, n_keep, n_keep, PtP);
    memcpy(PtP_inv, PtP, (size_t)n_keep * n_keep * sizeof(ST_double));
    if (ctools_cholesky(PtP_inv, n_keep) != 0 ||
        ctools_invert_from_cholesky(PtP_inv, n_keep, PtP_inv) != 0) {
        SF_error("civreghdfe: cannot partial out the partial() variables (P'WP is singular)\n");
        rc = 498;
        goto fwl_done;
    }
    ctools_matmul_atdb(P, data, w, N, n_keep, total_cols, PtWX);
    ctools_matmul_ab(PtP_inv, PtWX, n_keep, n_keep, total_cols, coef);

    /* Residualize all columns except the partial columns and their Z copies */
    #pragma omp parallel for schedule(dynamic) if ((size_t)N * total_cols > 100000)
    for (ST_int c = 0; c < total_cols; c++) {
        if (c >= exog_start && c < exog_start + K_exog && is_partial[c - exog_start]) continue;
        if (c >= zexog_start && c < zexog_start + K_exog && is_partial[c - zexog_start]) continue;
        ST_double *x = data + (size_t)c * N;
        const ST_double *b = coef + (size_t)c * n_keep;
        for (ST_int i = 0; i < N; i++) {
            ST_double fitted = 0.0;
            for (ST_int p = 0; p < n_keep; p++) fitted += P[(size_t)p * N + i] * b[p];
            x[i] -= fitted;
        }
    }

fwl_done:
    free(cols); free(coll); free(P); free(PtP); free(PtP_inv); free(PtWX); free(coef);
    return rc;
}


/* 1 if x is constant within the groups of one absorbed factor, which the
   fixed effects then remove exactly (conservatively 1 on allocation failure) */
static ST_int civ_const_in_factor(const ST_double *x, ST_int N, const HDFE_State *state)
{
    for (ST_int g = 0; g < state->G; g++) {
        const FE_Factor *f = &state->factors[g];
        ST_int maxl = (f->max_level > f->num_levels) ? f->max_level : f->num_levels;
        ST_double *first = (ST_double *)malloc(((size_t)maxl + 1) * sizeof(ST_double));
        unsigned char *seen = (unsigned char *)calloc((size_t)maxl + 1, 1);
        if (!first || !seen) { free(first); free(seen); return 1; }
        ST_int constant = 1;
        for (ST_int i = 0; i < N && constant; i++) {
            ST_int l = f->levels[i];
            if (l < 0 || l > maxl) constant = 0;
            else if (!seen[l]) { seen[l] = 1; first[l] = x[i]; }
            else if (x[i] != first[l]) constant = 0;
        }
        free(first); free(seen);
        if (constant) return 1;
    }
    return 0;
}

/*
    ivreghdfe checks collinearity before partialling, so a column absorbed by
    the fixed effects reaches its estimation as partialling noise, and it
    counts in ivreghdfe's ranks (degrees of freedom, J df) unless Stata's
    invsym() drops it, which it does when the noise has a sum of squares of
    at most 1e-19. That happens when the fixed effects remove the column
    exactly: it is constant within the groups of one absorbed factor, or the
    projection is exact (one factor, balanced panels), where the partialled
    column here is exactly zero too. Otherwise (a column absorbed only
    jointly by unbalanced fixed effects, e.g. age = year - birth year with
    absorb(id year)) reghdfe's iterative solver leaves noise above that. So
    an absorbed column counts when its partialled SS is not zero and it is
    not constant within any single factor (columns absorbed by partial()
    alone, without fixed effects, never count). It is still omitted here.
*/
static ST_int civ_absorbed_counts(ST_double ss_partial, const ST_double *orig, ST_int N,
                                  const HDFE_State *state)
{
    return state->G > 0 && ss_partial > 0.0 && !civ_const_in_factor(orig, N, state);
}

/*
    Full IV regression with HDFE.

    Variable layout from Stata:
    [y, X_endog_1, ..., X_endog_Ke, X_exog_1, ..., X_exog_Kx, Z_1, ..., Z_Kz, FE_1, ..., FE_G, cluster, weight,
     tsset panel, tsset time]

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
    /* Heap state released at iv_cleanup; every error path leaves through it */
    ST_retcode iv_rc = 0;
    ST_int *cluster_ids = NULL, *cluster2_ids = NULL, *valid_mask = NULL, *singleton_mask = NULL;
    ST_int **fe_levels_c = NULL;
    ST_int *level_max = NULL, *is_partial = NULL, *is_collinear_x = NULL, *test_idx = NULL, *ts_int = NULL;
    ST_double *y_c = NULL, *X_endog_c = NULL, *X_exog_c = NULL, *Z_c = NULL, *weights_c = NULL;
    ST_int *cluster_ids_c = NULL, *cluster2_ids_c = NULL;
    ST_double *all_data = NULL, *beta = NULL, *V = NULL, *first_stage_F = NULL;
    perm_idx_t *est_obs = NULL;
    HDFE_State *state = NULL;

    /* Start total timing */
    double t_total_start = ctools_timer_seconds();
    double t_load = 0, t_extract = 0, t_singleton = 0, t_hdfe_setup = 0, t_fwl = 0, t_partial_out = 0;
    double t_postproc = 0, t_estimate = 0, t_store = 0;

    ST_int in1 = SF_in1();
    ST_int in2 = SF_in2();
    ST_int N_total = in2 - in1 + 1;

    civreghdfe_input input;
    ST_retcode input_rc = civreghdfe_read_input(&input);
    if (input_rc) return input_rc;
    ST_int K_endog = input.K_endog;
    ST_int K_exog = input.K_exog;
    ST_int K_iv = input.K_iv;
    ST_int G = input.G;
    ST_int has_cluster = input.has_cluster;
    ST_int has_cluster2 = input.has_cluster2;
    ST_int has_weights = input.has_weights;
    ST_int weight_type = input.weight_type;
    ST_int vce_type = input.vce_type;
    ST_int maxiter = input.maxiter;
    ST_int verbose = input.verbose;
    ST_int est_method = input.est_method;
    ST_int kernel_type = input.kernel_type;
    ST_int bw = input.bw;
    ST_int kiefer = input.kiefer;
    ST_int ts_mode = input.ts_mode;
    ST_int ts_cols = input.ts_cols;
    ST_int dofminus = input.dofminus;
    ST_int sdofminus_opt = input.sdofminus_opt;
    ST_int nopartialsmall = input.nopartialsmall;
    ST_int center = input.center;
    ST_double tolerance = input.tolerance;
    ST_double kclass_user = input.kclass_user;
    ST_double fuller_alpha = input.fuller_alpha;
    ST_double tdelta = input.tdelta;

    if (center) {
        SF_error("civreghdfe: center is not supported\n");
        return 198;
    }
    /* The IV code indexes its K x K workspaces in int: j * K + i overflows
       once K exceeds 46340 (cqreg has the same cap) */
    if (K_iv > 46340 || K_exog + K_endog > 46340) {
        SF_error("civreghdfe: at most 46340 regressors and instruments are supported\n");
        return 198;
    }

    /* Calculate variable positions */
    /* Layout: [y, X_endog (Ke), X_exog (Kx), excluded instruments (Kz - Kx), FE (G),
                cluster?, cluster2?, weight?, tsset panel?, tsset time?] */
    ST_int var_y_idx = 0;  /* Index in loaded_data.vars[] */
    ST_int var_endog_start_idx = 1;
    ST_int var_exog_start_idx = var_endog_start_idx + K_endog;
    /* K_iv counts the exogenous regressors (the first columns of Z) and the
       excluded instruments; only the latter are sent a second time */
    const ST_int K_excl_in = K_iv - K_exog;
    ST_int var_excl_start_idx = var_exog_start_idx + K_exog;
    ST_int var_fe_start_idx = var_excl_start_idx + K_excl_in;
    ST_int var_cluster_idx = has_cluster ? (var_fe_start_idx + G) : -1;
    ST_int var_cluster2_idx = has_cluster2 ? (var_fe_start_idx + G + (has_cluster ? 1 : 0)) : -1;
    ST_int var_weight_idx = has_weights ? (var_fe_start_idx + G + (has_cluster ? 1 : 0) + (has_cluster2 ? 1 : 0)) : -1;
    /* tsset columns (panel, time) follow the weight */
    ST_int var_ts_start_idx = var_fe_start_idx + G + (has_cluster ? 1 : 0) + (has_cluster2 ? 1 : 0) + (has_weights ? 1 : 0);

    ST_int K_total = K_exog + K_endog;  /* Total regressors */

    /* Start data load timing */
    double t_load_start = ctools_timer_seconds();

    /* ================================================================
     * PARALLEL DATA LOADING using ctools_data_load
     * This handles if/in filtering at load time, loading only filtered observations.
     * ================================================================ */

    /* Calculate total variables to load */
    if (K_excl_in < 0) {
        SF_error("civreghdfe: fewer instruments than exogenous regressors\n");
        return 198;
    }
    ST_int total_vars = 1 + K_endog + K_exog + K_excl_in + G +
                        (has_cluster ? 1 : 0) + (has_cluster2 ? 1 : 0) + (has_weights ? 1 : 0) + ts_cols;

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

    stata_retcode load_rc = ctools_data_load(&filtered, var_indices, total_vars, 0, 0, CTOOLS_LOAD_CHECK_IF | CTOOLS_LOAD_NO_SORT_ORDER);
    free(var_indices);

    /* Capture pure SPI data load time */
    t_load = ctools_timer_seconds() - t_load_start;

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        SF_error("civreghdfe: Parallel data load failed\n");
        return ctools_stata_rc(load_rc);
    }

    /* Start extraction timer */
    double t_extract_start = ctools_timer_seconds();

    /* Get filtered observation count - data is already filtered */
    N_total = (ST_int)filtered.data.nobs;

    if (N_total == 0) {
        SF_error("civreghdfe: No observations selected\n");
        iv_rc = 2001;
        goto iv_cleanup;
    }

    /* Allocate cluster ID arrays (needed for string conversion + sentinel marking).
     * All other data (y, X, Z, FE, weights) is read directly from filtered.data
     * during missing check and compaction, avoiding redundant N_total-sized copies. */
    cluster_ids = has_cluster ? (ST_int *)ctools_safe_malloc2((size_t)N_total, sizeof(ST_int)) : NULL;
    cluster2_ids = has_cluster2 ? (ST_int *)ctools_safe_malloc2((size_t)N_total, sizeof(ST_int)) : NULL;

    if ((has_cluster && !cluster_ids) || (has_cluster2 && !cluster2_ids)) {
        SF_error("civreghdfe: Memory allocation failed\n");
        iv_rc = 920;
        goto iv_cleanup;
    }

    /* Extract cluster IDs (needs early processing for string variables) */
    if (has_cluster) {
        if (filtered.data.vars[var_cluster_idx].type == STATA_TYPE_STRING) {
            ST_int n_groups = 0;
            if (ctools_strings_to_cluster_ids(filtered.data.vars[var_cluster_idx].data.str,
                                               N_total, cluster_ids, &n_groups) != 0) {
                SF_error("civreghdfe: Failed to convert string cluster variable to group IDs\n");
                iv_rc = 920;
                goto iv_cleanup;
            }
        } else {
            double *src = filtered.data.vars[var_cluster_idx].data.dbl;
            ST_int n_groups = 0;
            if (ctools_numeric_to_cluster_ids(src, N_total, cluster_ids, &n_groups) != 0) {
                iv_rc = 920;
                goto iv_cleanup;
            }
        }
    }
    if (has_cluster2) {
        if (filtered.data.vars[var_cluster2_idx].type == STATA_TYPE_STRING) {
            ST_int n_groups = 0;
            if (ctools_strings_to_cluster_ids(filtered.data.vars[var_cluster2_idx].data.str,
                                               N_total, cluster2_ids, &n_groups) != 0) {
                SF_error("civreghdfe: Failed to convert string cluster2 variable to group IDs\n");
                iv_rc = 920;
                goto iv_cleanup;
            }
        } else {
            double *src = filtered.data.vars[var_cluster2_idx].data.dbl;
            ST_int n_groups = 0;
            if (ctools_numeric_to_cluster_ids(src, N_total, cluster2_ids, &n_groups) != 0) {
                iv_rc = 920;
                goto iv_cleanup;
            }
        }
    }

    /* End extraction, start missing value check timing.
     * Data arrays (y, X, Z, FE, weights) are NOT copied — we read directly
     * from filtered.data during missing value check and compaction. */
    t_extract = ctools_timer_seconds() - t_extract_start;
    double t_missing_start = ctools_timer_seconds();

    /* Drop observations with missing values. Columns are scanned one at a
     * time over blocks of rows (sequential reads; blocks in parallel),
     * directly from filtered.data. */
    valid_mask = (ST_int *)malloc((size_t)N_total * sizeof(ST_int));
    if (!valid_mask) {
        SF_error("civreghdfe: Memory allocation failed for valid_mask\n");
        iv_rc = 920;
        goto iv_cleanup;
    }
    ST_int N_valid = 0;
    {
        const ST_double missval = SV_missval;
        /* y, X_endog, X_exog, excluded instruments and FE variables are the
           first var_fe_start_idx + G columns */
        const ST_int n_lead = var_fe_start_idx + G;
        const ST_int block = 16384;
        /* N_total + block - 1 would overflow int near INT_MAX rows */
        const ST_int num_blocks = N_total / block + (N_total % block != 0);
        ST_int blk;
        #pragma omp parallel for schedule(static) reduction(+ : N_valid) if (N_total > 100000)
        for (blk = 0; blk < num_blocks; blk++) {
            const ST_int i0 = blk * block;
            const ST_int i1 = (N_total - i0 < block) ? N_total : i0 + block;
            for (ST_int ii = i0; ii < i1; ii++) valid_mask[ii] = 1;
            for (ST_int c = 0; c < n_lead; c++)
                civ_drop_missing(filtered.data.vars[c].data.dbl, i0, i1, missval, valid_mask);
            for (ST_int k = 0; k < ts_cols; k++)
                civ_drop_missing(filtered.data.vars[var_ts_start_idx + k].data.dbl, i0, i1, missval, valid_mask);
            if (has_weights) {
                const ST_double *w = filtered.data.vars[var_weight_idx].data.dbl;
                for (ST_int ii = i0; ii < i1; ii++) {
                    if (!(w[ii] < missval) || w[ii] <= 0) valid_mask[ii] = 0;
                }
            }
            /* missing cluster values were marked -1 during extraction */
            if (has_cluster) {
                for (ST_int ii = i0; ii < i1; ii++) if (cluster_ids[ii] < 0) valid_mask[ii] = 0;
            }
            if (has_cluster2) {
                for (ST_int ii = i0; ii < i1; ii++) if (cluster2_ids[ii] < 0) valid_mask[ii] = 0;
            }
            for (ST_int ii = i0; ii < i1; ii++) N_valid += valid_mask[ii];
        }
    }

    if (N_valid <= K_total || N_valid < K_iv) {
        SF_error("civreghdfe: Insufficient observations\n");
        iv_rc = 2001;
        goto iv_cleanup;
    }

    /* Compact FE levels from raw doubles with O(N) counting-based remap.
     * This replaces the old approach of: remap_values_sorted (O(N log N)) +
     * separate integer compaction. Now we compact and remap in one step. */
    fe_levels_c = (ST_int **)calloc((size_t)G, sizeof(ST_int *));
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
        iv_rc = 920;
        goto iv_cleanup;
    }
    {
        /* Compact FE doubles to temp array, then O(N) remap to contiguous integers */
        double *fe_compact = (double *)malloc((size_t)N_valid * sizeof(double));
        if (!fe_compact) {
            iv_rc = 920;
            goto iv_cleanup;
        }

        for (ST_int g = 0; g < G; g++) {
            double *fe_src = filtered.data.vars[var_fe_start_idx + g].data.dbl;
            ST_int idx = 0;
            for (ST_int ii = 0; ii < N_total; ii++) {
                if (!valid_mask[ii]) continue;
                fe_compact[idx++] = fe_src[ii];
            }

            ST_int num_levels_g = 0;
            ST_int *counts_g = NULL;
            if (remap_and_count(fe_compact, N_valid, fe_levels_c[g], &num_levels_g,
                                &counts_g, NULL, NULL) != 0) {
                SF_error("civreghdfe: FE level remapping failed\n");
                free(fe_compact);
                iv_rc = 920;
                goto iv_cleanup;
            }
            free(counts_g);  /* Counts recomputed later in HDFE setup */
        }
        free(fe_compact);
    }

    /* End missing value check/compact timing, start singleton timing */
    double t_missing = ctools_timer_seconds() - t_missing_start;
    double t_singleton_start = ctools_timer_seconds();

    /* Step 2: Singleton detection on compacted FE levels */
    singleton_mask = (ST_int *)malloc((size_t)N_valid * sizeof(ST_int));
    if (!singleton_mask) {
        SF_error("civreghdfe: Memory allocation failed for singleton mask\n");
        iv_rc = 920;
        goto iv_cleanup;
    }
    /* ivreghdfe selects the sample with reghdfe: fweights drop only levels
     * whose total fweight is 1; aweights and pweights use the unweighted rule. */
    ST_double *fw_c = NULL;
    if (has_weights && weight_type == 2) {
        fw_c = (ST_double *)malloc((size_t)N_valid * sizeof(ST_double));
        if (fw_c) {
            ST_int idx = 0;
            for (ST_int ii = 0; ii < N_total; ii++)
                if (valid_mask[ii]) fw_c[idx++] = filtered.data.vars[var_weight_idx].data.dbl[ii];
        }
    }
    ST_int num_singletons_total = (has_weights && weight_type == 2 && !fw_c) ? -1 :
        ctools_remove_singletons(fe_levels_c, G, N_valid, singleton_mask, fw_c, verbose);
    free(fw_c);
    if (num_singletons_total < 0) {
        SF_error("civreghdfe: Memory allocation failed in singleton detection\n");
        iv_rc = 920;
        goto iv_cleanup;
    }

    /* Step 3: Build combined index mapping original obs → final position.
     * Combine valid_mask and singleton_mask into one compaction pass. */
    ST_int N = N_valid - num_singletons_total;

    /* Allocate final arrays */
    y_c = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double));
    X_endog_c = (K_endog > 0) ? (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_endog, sizeof(ST_double)) : NULL;
    X_exog_c = (K_exog > 0) ? (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_exog, sizeof(ST_double)) : NULL;
    /* Z_c holds the K_iv instruments followed by the ts_cols tsset columns
       (panel, time), which share its lifetime */
    Z_c = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)(K_iv + ts_cols), sizeof(ST_double));
    ST_double *ts_c = Z_c ? Z_c + (size_t)N * K_iv : NULL;
    weights_c = has_weights ? (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double)) : NULL;
    cluster_ids_c = has_cluster ? (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int)) : NULL;
    cluster2_ids_c = has_cluster2 ? (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int)) : NULL;

    if (!y_c || !Z_c || (K_endog > 0 && !X_endog_c) || (K_exog > 0 && !X_exog_c) ||
        (has_weights && !weights_c) || (has_cluster && !cluster_ids_c) ||
        (has_cluster2 && !cluster2_ids_c)) {
        SF_error("civreghdfe: Memory allocation failed for compacted arrays\n");
        iv_rc = 920;
        goto iv_cleanup;
    }

    /* Compaction: the retained rows (valid and not singleton) of every
     * column, read directly from filtered.data, one column per thread. */
    {
        ST_int *rows = (ST_int *)ctools_safe_malloc2((size_t)(N > 0 ? N : 1), sizeof(ST_int));
        if (!rows) {
            SF_error("civreghdfe: Memory allocation failed for compacted arrays\n");
            iv_rc = 920;
            goto iv_cleanup;
        }
        ST_int valid_idx = 0;  /* Index into singleton_mask (N_valid-sized) */
        ST_int out_idx = 0;    /* Index into final arrays (N-sized) */
        for (ST_int ii = 0; ii < N_total; ii++) {
            if (!valid_mask[ii]) continue;
            if (singleton_mask[valid_idx++]) rows[out_idx++] = ii;
        }
        for (ST_int j = 0; j < N; j++) {
            filtered.obs_map[j] = filtered.obs_map[rows[j]];
            if (has_weights) weights_c[j] = filtered.data.vars[var_weight_idx].data.dbl[rows[j]];
            if (has_cluster) cluster_ids_c[j] = cluster_ids[rows[j]];
            if (has_cluster2) cluster2_ids_c[j] = cluster2_ids[rows[j]];
        }
        /* y, X_endog, X_exog, excluded instruments, tsset columns */
        const ST_int n_gather = 1 + K_endog + K_exog + K_excl_in + ts_cols;
        ST_int c;
        #pragma omp parallel for schedule(dynamic) if ((size_t)N * (size_t)n_gather > 200000)
        for (c = 0; c < n_gather; c++) {
            const ST_double *src;
            ST_double *dst;
            if (c == 0) {
                src = filtered.data.vars[var_y_idx].data.dbl;
                dst = y_c;
            } else if (c < 1 + K_endog) {
                src = filtered.data.vars[var_endog_start_idx + (c - 1)].data.dbl;
                dst = X_endog_c + (size_t)(c - 1) * N;
            } else if (c < 1 + K_endog + K_exog) {
                src = filtered.data.vars[var_exog_start_idx + (c - 1 - K_endog)].data.dbl;
                dst = X_exog_c + (size_t)(c - 1 - K_endog) * N;
            } else if (c < 1 + K_endog + K_exog + K_excl_in) {
                ST_int k = c - 1 - K_endog - K_exog;
                src = filtered.data.vars[var_excl_start_idx + k].data.dbl;
                dst = Z_c + (size_t)(K_exog + k) * N;
            } else {
                ST_int k = c - 1 - K_endog - K_exog - K_excl_in;
                src = filtered.data.vars[var_ts_start_idx + k].data.dbl;
                dst = ts_c + (size_t)k * N;
            }
            civ_gather(dst, src, rows, N);
        }
        free(rows);
        /* Z = [X_exog, excluded instruments] */
        if (K_exog > 0) memcpy(Z_c, X_exog_c, (size_t)N * K_exog * sizeof(ST_double));
    }

    /* Compact FE levels in place using singleton_mask. The write index never
     * exceeds the read index, so this cannot fail after the data were compacted. */
    if (num_singletons_total > 0) {
        for (ST_int g = 0; g < G; g++) {
            ctools_compact_array_int(fe_levels_c[g], fe_levels_c[g], singleton_mask, N_valid, N);
        }
    }

    free(singleton_mask);
    singleton_mask = NULL;

    /* Weight normalization:
       - aw (type=1) and pw (type=3): normalize weights so sum(w) = N
       - fw (type=2): compute N_eff = sum(fw) for reporting/DOF */
    ST_double N_eff = (ST_double)N;
    if (has_weights && weights_c) {
        if (weight_type == 1 || weight_type == 3) {
            /* aweight or pweight: normalize so sum(w) = N */
            ST_double sum_w_raw = 0.0;
            for (ST_int i = 0; i < N; i++) sum_w_raw += weights_c[i];
            if (sum_w_raw > 0.0) {
                ST_double scale = (ST_double)N / sum_w_raw;
                for (ST_int i = 0; i < N; i++) weights_c[i] *= scale;
            }
        } else if (weight_type == 2) {
            /* fweight: N_eff = sum(fw) */
            ST_double sum_fw = 0.0;
            for (ST_int i = 0; i < N; i++) sum_fw += weights_c[i];
            N_eff = sum_fw;
        }
    }

    /* Retain the compacted observation mapping for residual output. */
    est_obs = filtered.obs_map;
    filtered.obs_map = NULL;

    /* Free original arrays and filtered data (kept alive for direct reads during compaction) */
    free(cluster_ids); free(cluster2_ids); free(valid_mask);
    cluster_ids = cluster2_ids = valid_mask = NULL;
    ctools_filtered_data_free(&filtered);

    /* Check we still have enough observations */
    if (N <= K_total || N < K_iv) {
        SF_error("civreghdfe: Insufficient observations after singleton removal\n");
        iv_rc = 2001;
        goto iv_cleanup;
    }

    /* End singleton timing, start HDFE setup timing */
    t_singleton = ctools_timer_seconds() - t_singleton_start;
    double t_hdfe_setup_start = ctools_timer_seconds();

    /* Initialize the HDFE state on the estimation sample */
    state = (HDFE_State *)calloc(1, sizeof(HDFE_State));
    level_max = (ST_int *)calloc((size_t)G + 1, sizeof(ST_int));
    if (!state || !level_max) {
        SF_error("civreghdfe: Memory allocation failed for HDFE state\n");
        iv_rc = 920;
        goto iv_cleanup;
    }
    state->K = 1 + K_endog + K_exog + K_iv;  /* Total columns to demean */
    state->in1 = 1;
    state->in2 = N;
    state->has_weights = has_weights;
    state->weight_type = weight_type;
    state->weights = weights_c;
    state->maxiter = maxiter;
    state->tolerance = tolerance;
    state->verbose = verbose;

    /* Factors: the CG solver expects levels 1..num_levels, so the levels
       left after singleton removal are renumbered and counted */
    for (ST_int g = 0; g < G; g++) {
        for (ST_int i = 0; i < N; i++) {
            if (fe_levels_c[g][i] > level_max[g]) level_max[g] = fe_levels_c[g][i];
        }
    }
    ST_retcode factor_rc = hdfe_build_factors(state, G, fe_levels_c, level_max, N, NULL, N,
                                              has_weights ? weights_c : NULL, 0);
    free(level_max);
    level_max = NULL;
    for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
    free(fe_levels_c);
    fe_levels_c = NULL;
    if (factor_rc) {
        SF_error("civreghdfe: Memory allocation failed for FE factors\n");
        iv_rc = factor_rc;
        goto iv_cleanup;
    }

    /* Allocate inv_counts, inv_weighted_counts, and CG solver buffers */
    ST_int num_threads = ctools_get_max_threads();
    if (num_threads < 1) num_threads = 1;
    state->num_threads = num_threads;

    /* Columns partialled out: y, endogenous, exogenous and excluded
       instruments (the buffer count is capped at this number) */
    ST_int max_partial_cols = 1 + K_endog + K_iv;
    if (ctools_hdfe_alloc_buffers(state, 0, max_partial_cols) != 0) {
        SF_error("civreghdfe: Memory allocation failed for HDFE buffers\n");
        iv_rc = 920;
        goto iv_cleanup;
    }

    state->factors_initialized = 1;

    /* Record this estimator's active state for lifecycle tracking */
    iv_state = state;

    /* End HDFE setup timing, start FWL timing */
    t_hdfe_setup = ctools_timer_seconds() - t_hdfe_setup_start;
    double t_fwl_start = ctools_timer_seconds();

    /* FWL partialling: partial() arrives as one 0/1 flag per X_exog column
       (matrix __civreghdfe_partial_flags, aligned with the expanded exogenous
       regressors). is_partial: 1 = partialled out, 2 = dropped as redundant. */
    ST_double dval_partial;
    ST_int n_partial = 0;

    if (SF_scal_use("__civreghdfe_n_partial", &dval_partial) == 0) {
        n_partial = (ST_int)dval_partial;
    }

    /* exogenous regressors as sent (the columns of X_exog_c and the first
       ones of Z_c), before partial() removes some */
    const ST_int K_exog_in = K_exog;

    if (n_partial > 0 && K_exog > 0) {
        is_partial = (ST_int *)calloc(K_exog, sizeof(ST_int));
        ST_int n_flagged = is_partial ? read_flag_matrix("__civreghdfe_partial_flags", K_exog, is_partial) : -1;
        if (n_flagged != n_partial) {
            ST_retcode flag_rc = is_partial ? 198 : 920;
            SF_error(is_partial ? "civreghdfe: partial() column flags do not match the exogenous regressors\n"
                                : "civreghdfe: Memory allocation failed for partial indices\n");
            iv_rc = flag_rc;
            goto iv_cleanup;
        }
    }

    /* End FWL timing, start partial out timing */
    t_fwl = ctools_timer_seconds() - t_fwl_start;
    double t_partial_start = ctools_timer_seconds();

    /* Partial out FEs from all variables using CG solver */
    ST_int total_cols = 1 + K_endog + K_exog + K_iv;

    /* Combine all data into one array for parallel processing. One spare row
       (total_cols values after the N x total_cols block) holds each column's
       pre-partialling sum of squares; it shares all_data's lifetime. */
    all_data = (ST_double *)ctools_safe_malloc3((size_t)N + 1, (size_t)total_cols, sizeof(ST_double));
    if (!all_data) {
        SF_error("civreghdfe: Memory allocation failed for demeaning buffer\n");
        iv_rc = 920;
        goto iv_cleanup;
    }

    /* Column blocks of all_data: Z = [copy of X_exog | excluded instruments],
       then y, X_endog and X_exog. The fixed effects are partialled out of
       the columns from the excluded instruments on; Z's copy of X_exog is
       filled from the partialled X_exog afterwards, which gives exactly what
       partialling it separately would. */
    const ST_int col_excl = K_exog;
    const ST_int col_y = K_iv;
    const ST_int col_endog = K_iv + 1;
    const ST_int col_exog = K_iv + 1 + K_endog;
    if (K_excl_in > 0) {
        memcpy(all_data + (size_t)N * col_excl, Z_c + (size_t)N * K_exog, (size_t)N * K_excl_in * sizeof(ST_double));
    }
    memcpy(all_data + (size_t)N * col_y, y_c, N * sizeof(ST_double));
    if (K_endog > 0) {
        memcpy(all_data + (size_t)N * col_endog, X_endog_c, (size_t)N * K_endog * sizeof(ST_double));
    }
    if (K_exog > 0) {
        memcpy(all_data + (size_t)N * col_exog, X_exog_c, (size_t)N * K_exog * sizeof(ST_double));
    }

    /* Pre-partialling (weighted) sum of squared deviations of every column.
       A column whose partialled SS is at most absorb_tol times this was
       absorbed by the fixed effects (or by partial()) and is omitted, as
       ivreghdfe does; absorb_tol is reghdfe's collinear-with-FE tolerance. */
    const ST_double absorb_tol = (tolerance / 10.0 < 1e-6) ? tolerance / 10.0 : 1e-6;
    ST_double *ss_pre = all_data + (size_t)N * total_cols;
    #pragma omp parallel for schedule(dynamic) if ((size_t)N * total_cols > 100000)
    for (ST_int c = col_excl; c < total_cols; c++) {
        ss_pre[c] = civ_wss_dev(all_data + (size_t)c * N, weights_c, N);
    }
    for (ST_int k = 0; k < K_exog; k++) ss_pre[k] = ss_pre[col_exog + k];
    const ST_double *ss_endog = ss_pre + col_endog;
    ST_double *ss_exog = ss_pre + col_exog;  /* compacted with X_exog below */
    const ST_double *ss_excl = ss_pre + col_excl;

    /* Partial out the FEs from every column in one pass, then the partial()
       variables by weighted FWL on the FE-demeaned data: the residual of the
       W-weighted projection on [FE, P] is M_P~ M_FE x with P~ = M_FE P, the
       way ivreghdfe runs HDFE partial_out and then s_partial. */
    HDFE_SolveResult projection = partial_out_columns(state, all_data + (size_t)N * col_excl, N,
                                                      total_cols - col_excl, num_threads);
    ST_retcode fwl_rc = 0;
    if (!projection.status && K_exog > 0) {
        memcpy(all_data, all_data + (size_t)N * col_exog, (size_t)N * K_exog * sizeof(ST_double));
    }
    if (!projection.status && is_partial) {
        fwl_rc = fwl_partial_out(all_data, N, total_cols, col_exog, 0, K_exog,
                                 is_partial, weights_c, ss_exog, absorb_tol);
    }

    if (projection.status || fwl_rc) {
        if (projection.status) SF_error("civreghdfe: fixed-effect projection did not converge\n");
        iv_rc = projection.status ? projection.status : fwl_rc;
        goto iv_cleanup;
    }

    /* End partial out timing, start post-processing timing */
    t_partial_out = ctools_timer_seconds() - t_partial_start;
    double t_postproc_start = ctools_timer_seconds();

    /* Extract demeaned data */
    ST_double *y_dem = all_data + (size_t)N * col_y;
    ST_double *X_endog_dem = (K_endog > 0) ? all_data + (size_t)N * col_endog : NULL;
    ST_double *X_exog_dem = (K_exog > 0) ? all_data + (size_t)N * col_exog : NULL;
    ST_double *Z_dem = all_data;

    /* Remove partial variables from X_exog AND from Z (instruments) after convergence.
       Z contains [exog_vars, excluded_instruments], so partial exog vars must be removed
       from both X_exog and the first K_exog columns of Z. Partial columns dropped as
       redundant (flag 2) are removed as well, and reported back for the ado's note. */
    if (is_partial) {
        ST_int K_exog_new = K_exog - n_partial;
        ST_int K_excl = K_iv - K_exog;  /* Number of excluded instruments */

        for (ST_int k = 0; k < K_exog; k++) {
            if (is_partial[k] == 2 &&
                (iv_rc = ctools_mat_store("__civreghdfe_partial_flags", 1, k + 1, 2.0)))
                goto iv_cleanup;
        }

        /* Remove partial columns from X_exog (and their pre-partialling SS) */
        if (K_exog_new > 0) {
            ST_int new_idx = 0;
            for (ST_int k = 0; k < K_exog; k++) {
                if (!is_partial[k]) {
                    if (new_idx != k) {
                        memcpy(X_exog_dem + (size_t)new_idx * N, X_exog_dem + (size_t)k * N, N * sizeof(ST_double));
                        ss_exog[new_idx] = ss_exog[k];
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
                    memcpy(Z_dem + (size_t)new_idx * N, Z_dem + (size_t)k * N, N * sizeof(ST_double));
                }
                new_idx++;
            }
        }
        /* Shift excluded instruments down */
        for (ST_int k = 0; k < K_excl; k++) {
            memcpy(Z_dem + (size_t)new_idx * N, Z_dem + (size_t)(K_exog + k) * N, N * sizeof(ST_double));
            new_idx++;
        }

        K_exog = K_exog_new;
        K_iv = K_exog_new + K_excl;
        X_exog_dem = (K_exog > 0) ? all_data + (size_t)N * col_exog : NULL;
        /* is_partial stays until the absorbed columns are counted below: it
           maps the remaining X_exog columns back to X_exog_c */
    }

    /* Compute df_a (absorbed degrees of freedom) and mobility groups */
    ST_int df_a = 0;
    ST_int df_a_nested = 0;  /* Levels from FE nested within cluster */
    ST_int mobility_groups = 0;

    ctools_compute_hdfe_dof(state->factors, G, N, 0, &df_a, &mobility_groups);

    for (ST_int g = 0; g < G; g++) {
        /* Save per-FE num_levels to Stata scalars for absorbed DOF table */
        char scalar_name[64];
        snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_num_levels_%d", (int)(g + 1));
        if ((iv_rc = ctools_scal_save(scalar_name, (ST_double)state->factors[g].num_levels)))
            goto iv_cleanup;
    }

    /* Save mobility groups for absorbed DOF table display */
    if ((iv_rc = ctools_scal_save("__civreghdfe_mobility_groups", (ST_double)mobility_groups)))
        goto iv_cleanup;

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
            iv_rc = 920;
            goto iv_cleanup;
        }
    }

    /* Remap cluster2 IDs for two-way clustering */
    ST_int num_clusters2 = 0;
    if (has_cluster2) {
        if (ctools_remap_cluster_ids(cluster2_ids_c, N, &num_clusters2) != 0) {
            SF_error("civreghdfe: Failed to remap cluster2 IDs (memory allocation error)\n");
            iv_rc = 920;
            goto iv_cleanup;
        }
    }

    if (!kiefer && (vce_type == 2 || vce_type == 3) &&
        (num_clusters < 2 || (has_cluster2 && num_clusters2 < 2))) {
        SF_error("civreghdfe: clustered VCE requires at least two retained clusters in each dimension\n");
        iv_rc = 459;
        goto iv_cleanup;
    }

    /* Detect FEs nested within cluster variable using data-based check.
     * This works for both numeric and string cluster variables. */
    /* A failed allocation must not silently skip the detection: it changes
     * the cluster dof */
    if (has_cluster) {
        df_a_nested = 0;
        for (ST_int g = 0; g < G; g++) {
            ST_int is_nested = ctools_fe_nested_in_cluster(
                state->factors[g].levels, state->factors[g].num_levels,
                cluster_ids_c, N);
            if (is_nested < 0) {
                SF_error("civreghdfe: Memory allocation failed for the nesting check\n");
                iv_rc = 920;
                goto iv_cleanup;
            }
            if (is_nested) df_a_nested += state->factors[g].num_levels;

            char scalar_name[64];
            snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_fe_nested_%d", (int)(g + 1));
            if ((iv_rc = ctools_scal_save(scalar_name, (ST_double)is_nested)))
                goto iv_cleanup;
        }
        /* Recompute VCE DOF with data-based nesting results */
        df_a_for_vce = df_a - df_a_nested;
        nested_adj = (df_a_nested > 0) ? 1 : 0;
    }

    /* For Kiefer, the panel variable is set as "cluster" for HAC lag structure,
       not true clustering. The FE nesting detection incorrectly sets df_a_for_vce=0
       because the FE and "cluster" are the same variable. Override: use full df_a
       to match ivreghdfe's sdofminus = absorb_ct. */
    if (kiefer) {
        df_a_for_vce = df_a;
        nested_adj = 0;
    }

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
    is_collinear_x = (ST_int *)calloc(K_total, sizeof(ST_int));
    ST_int num_collinear_x = 0;

    if (!is_collinear_x) {
        SF_error("civreghdfe: Memory allocation failed for collinearity arrays\n");
        iv_rc = 920;
        goto iv_cleanup;
    }

    /* Stage 1: columns absorbed by the fixed effects (or by partial()): the
       partialled (weighted) SS is zero or negligible relative to the
       pre-partialling SS (reghdfe's collinear-with-FE rule). What remains is
       numerical noise; it is omitted and must not reach stage 2, whose
       equilibration would rescale it to unit variance. Those that ivreghdfe
       keeps in its ranks still count in the degrees of freedom
       (civ_absorbed_counts). Every absorbed exogenous regressor counts in the
       Cragg-Donald and KP Wald F dof, whose instrument count ivreghdfe takes
       from its variable lists (ranktest partials those regressors out, so
       the statistics themselves do not change). An absorbed excluded
       instrument makes ivreghdfe's statistics missing (exact zero) or
       noise-based, so there the model without it is the reference. */
    ST_int n_abs_x_exog = 0, n_abs_x_endog = 0, n_abs_z_excl = 0, n_abs_exog_all = 0;
    {
        ST_int j_in = 0;   /* X_exog_c column of X_exog column k (partial() removed) */
        for (ST_int k = 0; k < K_exog; k++, j_in++) {
            while (is_partial && j_in < K_exog_in && is_partial[j_in]) j_in++;
            ST_double xx_partial = civ_wss(X_exog_dem + (size_t)k * N, weights_c, N);
            if (xx_partial < 1e-30 || xx_partial <= absorb_tol * ss_exog[k]) {
                is_collinear_x[k] = 1;
                num_collinear_x++;
                n_abs_exog_all++;
                if (j_in < K_exog_in &&
                    civ_absorbed_counts(xx_partial, X_exog_c + (size_t)j_in * N, N, state))
                    n_abs_x_exog++;
            }
        }
        for (ST_int k = 0; k < K_endog; k++) {
            ST_double xx_partial = civ_wss(X_endog_dem + (size_t)k * N, weights_c, N);
            if (xx_partial < 1e-30 || xx_partial <= absorb_tol * ss_endog[k]) {
                is_collinear_x[K_exog + k] = 1;
                num_collinear_x++;
                if (civ_absorbed_counts(xx_partial, X_endog_c + (size_t)k * N, N, state))
                    n_abs_x_endog++;
            }
        }
    }
    free(is_partial);
    is_partial = NULL;

    /* Stage 2: Numerical collinearity via Cholesky on X'X
       Build concatenated X = [X_exog, X_endog] in column-major order */
    if (K_total > 1) {
        ST_double *XtX = (ST_double *)ctools_safe_calloc3((size_t)K_total, (size_t)K_total, sizeof(ST_double));
        if (!XtX) {
            SF_error("civreghdfe: Memory allocation failed for XtX\n");
            iv_rc = 920;
            goto iv_cleanup;
        }

        /* Compute X'X where X = [X_exog_dem, X_endog_dem] using block matmul.
         * XtX (K_total x K_total) has blocks:
         *   [Xe'Xe  Xe'Xn]
         *   [Xn'Xe  Xn'Xn]
         * Use ctools_matmul_atb for each block (OpenMP + SIMD) instead of
         * K_total² individual fast_dot calls. */
        if (K_exog > 0 && K_endog > 0) {
            /* Both exog and endog present — compute 3 blocks (4th by symmetry) */
            ST_double *blk_ee = (ST_double *)ctools_safe_malloc3((size_t)K_exog, (size_t)K_exog, sizeof(ST_double));
            ST_double *blk_en = (ST_double *)ctools_safe_malloc3((size_t)K_exog, (size_t)K_endog, sizeof(ST_double));
            ST_double *blk_nn = (ST_double *)ctools_safe_malloc3((size_t)K_endog, (size_t)K_endog, sizeof(ST_double));
            if (blk_ee && blk_en && blk_nn) {
                ctools_matmul_ata(X_exog_dem, N, K_exog, blk_ee);
                ctools_matmul_atb(X_exog_dem, X_endog_dem, N, K_exog, K_endog, blk_en);
                ctools_matmul_ata(X_endog_dem, N, K_endog, blk_nn);

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
                    const ST_double *xi = (ii < K_exog) ? X_exog_dem + (size_t)ii * N
                                                        : X_endog_dem + (size_t)(ii - K_exog) * N;
                    for (ST_int jj = ii; jj < K_total; jj++) {
                        const ST_double *xj = (jj < K_exog) ? X_exog_dem + (size_t)jj * N
                                                            : X_endog_dem + (size_t)(jj - K_exog) * N;
                        ST_double val = fast_dot(xi, xj, N);
                        XtX[jj * K_total + ii] = val;
                        XtX[ii * K_total + jj] = val;
                    }
                }
            }
            free(blk_ee); free(blk_en); free(blk_nn);
        } else if (K_exog > 0) {
            ctools_matmul_ata(X_exog_dem, N, K_exog, XtX);
        } else {
            ctools_matmul_ata(X_endog_dem, N, K_endog, XtX);
        }

        /* Stage-1 columns are excluded from the factorization and stay
         * flagged; the rest are screened on the equilibrated X'X, which keeps
         * mixed-scale variables (e.g., x_tiny / x_huge) from producing false
         * positives. */
        ST_int num_chol_collinear = detect_collinearity_scaled(XtX, K_total, is_collinear_x, verbose);
        free(XtX);

        if (num_chol_collinear < 0) {
            SF_error("civreghdfe: Collinearity detection failed\n");
            iv_rc = 920;
            goto iv_cleanup;
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
        if ((iv_rc = ctools_scal_save("__civreghdfe_num_collinear", (ST_double)num_collinear_x)))
            goto iv_cleanup;
        for (ST_int k = 0; k < K_total_orig; k++) {
            snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_collinear_%d", (int)(k + 1));
            if ((iv_rc = ctools_scal_save(scalar_name, (ST_double)is_collinear_x[k])))
                goto iv_cleanup;
        }
    }

    /* Compact X arrays in-place if collinear columns found */
    if (num_collinear_x > 0) {
        /* Compact X_exog_dem */
        ST_int new_K_exog = 0;
        for (ST_int k = 0; k < K_exog; k++) {
            if (!is_collinear_x[k]) {
                if (new_K_exog != k) {
                    memcpy(X_exog_dem + (size_t)new_K_exog * N, X_exog_dem + (size_t)k * N, N * sizeof(ST_double));
                }
                new_K_exog++;
            }
        }

        /* Compact X_endog_dem */
        ST_int new_K_endog = 0;
        for (ST_int k = 0; k < K_endog; k++) {
            if (!is_collinear_x[K_exog + k]) {
                if (new_K_endog != k) {
                    memcpy(X_endog_dem + (size_t)new_K_endog * N, X_endog_dem + (size_t)k * N, N * sizeof(ST_double));
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
                    memcpy(Z_dem + (size_t)new_z_idx * N, Z_dem + (size_t)k * N, N * sizeof(ST_double));
                }
                new_z_idx++;
            }
        }
        /* Shift excluded instruments down after compacted exog portion */
        for (ST_int k = 0; k < K_excl; k++) {
            if (new_z_idx != K_exog + k) {
                memcpy(Z_dem + (size_t)new_z_idx * N, Z_dem + (size_t)(K_exog + k) * N, N * sizeof(ST_double));
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

    /* Z collinearity detection. Excluded instruments absorbed by the fixed
       effects (the stage-1 rule above) get zero rows and columns of Z'Z
       before the equilibration, so they are dropped rather than kept as
       rescaled noise. The exogenous block of Z was checked with X. */
    ST_int K_excl_z = K_iv - K_exog;  /* excluded instruments as sent by the ado */
    ST_int *is_collinear_z = (ST_int *)calloc((size_t)(K_iv > 0 ? K_iv : 1), sizeof(ST_int));
    ST_int num_collinear_z = 0;

    if (is_collinear_z) {
        for (ST_int j = 0; j < K_excl_z; j++) {
            ST_double zz_partial = civ_wss(Z_dem + (size_t)(K_exog + j) * N, weights_c, N);
            if (zz_partial < 1e-30 || zz_partial <= absorb_tol * ss_excl[j]) {
                is_collinear_z[K_exog + j] = 1;
                if (civ_absorbed_counts(zz_partial, Z_c + (size_t)(K_exog_in + j) * N, N, state))
                    n_abs_z_excl++;
            }
        }

        ST_double *ZtZ = (K_iv > 1) ? (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_iv, sizeof(ST_double)) : NULL;
        if (ZtZ) {
            ctools_matmul_ata(Z_dem, N, K_iv, ZtZ);
            detect_collinearity_scaled(ZtZ, K_iv, is_collinear_z, verbose);
        }
        free(ZtZ);

        /* Compact Z_dem in-place */
        ST_int new_K_iv = 0;
        for (ST_int k = 0; k < K_iv; k++) {
            if (!is_collinear_z[k]) {
                if (new_K_iv != k) {
                    memcpy(Z_dem + (size_t)new_K_iv * N, Z_dem + (size_t)k * N, N * sizeof(ST_double));
                }
                new_K_iv++;
            }
        }
        num_collinear_z = K_iv - new_K_iv;
        K_iv = new_K_iv;

        if (verbose && num_collinear_z > 0) {
            char msg[128];
            snprintf(msg, sizeof(msg), "civreghdfe: Dropped %d collinear instrument(s)\n",
                     (int)num_collinear_z);
            SF_display(msg);
        }
    }

    /* orthog()/endogtest()/redundant(): map the ado's per-column flags onto
       the excluded instruments and endogenous regressors that survived the
       checks above; test_idx holds the three position lists back to back. */
    ST_int n_orthog_t = 0, n_endogtest_t = 0, n_redund_t = 0;
    {
        ST_double dv;
        if (SF_scal_use("__civreghdfe_n_orthog", &dv) == 0 && dv > 0) n_orthog_t = (ST_int)dv;
        if (SF_scal_use("__civreghdfe_n_endogtest", &dv) == 0 && dv > 0) n_endogtest_t = (ST_int)dv;
        if (SF_scal_use("__civreghdfe_n_redundant", &dv) == 0 && dv > 0) n_redund_t = (ST_int)dv;
    }
    test_idx = (ST_int *)malloc(((size_t)n_orthog_t + n_endogtest_t + n_redund_t + 1) * sizeof(ST_int));
    ST_retcode test_rc = test_idx ? 0 : 920;
    if (test_rc) SF_error("civreghdfe: Memory allocation failed for test columns\n");
    if (!test_rc)
        test_rc = map_test_cols("__civreghdfe_orthog_flags", "orthog", 1, n_orthog_t, K_excl_z,
                                is_collinear_z ? is_collinear_z + K_exog : NULL, test_idx);
    if (!test_rc)
        test_rc = map_test_cols("__civreghdfe_endogtest_flags", "endogtest", 2, n_endogtest_t, K_endog_orig,
                                is_collinear_x + K_exog_orig, test_idx + n_orthog_t);
    if (!test_rc)
        test_rc = map_test_cols("__civreghdfe_redundant_flags", "redundant", 3, n_redund_t, K_excl_z,
                                is_collinear_z ? is_collinear_z + K_exog : NULL,
                                test_idx + n_orthog_t + n_endogtest_t);
    free(is_collinear_z);

    /* Post-compaction identification check */
    if (K_iv < K_total || test_rc) {
        if (!test_rc) SF_error("civreghdfe: Model is underidentified after removing collinear variables\n");
        iv_rc = test_rc ? test_rc : 481;
        goto iv_cleanup;
    }

    /* Store compacted dimensions */
    if ((iv_rc = ctools_scal_save("__civreghdfe_K_keep", (ST_double)K_total)))
        goto iv_cleanup;

    /* Allocate output arrays */
    beta = (ST_double *)calloc(K_total, sizeof(ST_double));
    V = (ST_double *)ctools_safe_calloc3((size_t)K_total, (size_t)K_total, sizeof(ST_double));
    first_stage_F = (ST_double *)calloc(K_endog, sizeof(ST_double));

    /* Check allocations - critical for preventing NULL pointer dereference */
    if (!beta || !V || !first_stage_F) {
        SF_error("civreghdfe: Memory allocation failed for output arrays\n");
        iv_rc = 920;
        goto iv_cleanup;
    }

    /* End post-processing timing, start estimation timing */
    t_postproc = ctools_timer_seconds() - t_postproc_start;
    double t_estimate_start = ctools_timer_seconds();

    /* Compute k-class IV estimation (2SLS, LIML, Fuller, etc.) */
    /* For cluster VCE, pass df_a_for_vce (excluding nested FE levels) and nested_adj */
    ST_retcode rc = 0;

    /* tsset structure for the kernel estimators: panel ids and the time
       index (t - tmin)/tdelta, so that kernel lags are time differences
       within a panel (ivreg2's L. operators) and Driscoll-Kraay periods
       follow calendar order */
    civreghdfe_ts ts_info;
    memset(&ts_info, 0, sizeof(ts_info));
    ts_int = NULL;  /* panel ids (N), then time index (N) */
    if (ts_mode > 0) {
        ts_int = (ST_int *)ctools_safe_malloc3((size_t)N, 2, sizeof(ST_int));
        if (!ts_int) {
            SF_error("civreghdfe: Memory allocation failed for the time structure\n");
            rc = 920;
        } else {
            const ST_double *tv = ts_c + (size_t)(ts_cols - 1) * N;
            ST_int *tidx = ts_int + N;
            ST_double tmin = tv[0], tmax = tv[0];
            for (ST_int i = 1; i < N; i++) {
                if (tv[i] < tmin) tmin = tv[i];
                if (tv[i] > tmax) tmax = tv[i];
            }
            /* The time index is an int: refuse spans of 2^31-2 or more delta
               units, which would wrap it and break the kernel lag pairing
               (e.g. a %tc clock with delta 1 ms over more than 24.8 days) */
            if (!((tmax - tmin) / tdelta < 2147483646.0)) {
                SF_error("civreghdfe: the tsset time span is too long for kernel estimators "
                         "(2^31-2 or more periods of delta); tsset with a coarser delta()\n");
                rc = 198;
            }
            ST_int tmax_idx = 0;
            for (ST_int i = 0; i < N && !rc; i++) {
                tidx[i] = (ST_int)llround((tv[i] - tmin) / tdelta);
                if (tidx[i] > tmax_idx) tmax_idx = tidx[i];
            }
            /* Time span of the estimation sample: the ado reports T = span + 1
               as Kiefer's bandwidth; as ivreg2, a bandwidth may not exceed
               span / delta */
            if (!rc) rc = ctools_scal_save("__civreghdfe_tspan", tmax - tmin);
            if (!rc && kernel_type > 0 && !kiefer && (ST_double)bw > (tmax - tmin) / tdelta) {
                SF_error("invalid bandwidth in option bw() - cannot exceed timespan of data\n");
                rc = 198;
            }
            ST_int num_panels = 1;
            if (ts_mode == 2 && !rc) {
                if (ctools_numeric_to_cluster_ids(ts_c, N, ts_int, &num_panels) != 0) {
                    SF_error("civreghdfe: Memory allocation failed for the panel structure\n");
                    rc = 920;
                }
                for (ST_int i = 0; i < N && !rc; i++) ts_int[i] += 1;
            }
            ts_info.panel = (ts_mode == 2) ? ts_int : NULL;
            ts_info.num_panels = num_panels;
            ts_info.tidx = tidx;
            ts_info.tmax_idx = tmax_idx;
            ts_info.tdelta = tdelta;
        }
    }

    civreghdfe_test_cols test_cols = {
        test_idx, n_orthog_t,
        test_idx + n_orthog_t, n_endogtest_t,
        test_idx + n_orthog_t + n_endogtest_t, n_redund_t
    };

    /* ivreghdfe with absorb() implies small and partialsmall: its small-sample
       sdofminus is sdofminus() + the partialled-out regressors + absorb_ct
       (max(1, df_a), the nested-in-cluster convention), and nopartialsmall
       is ignored. dofminus() is the large-sample adjustment. */
    ST_int sdof_extra = sdofminus_opt + ((G > 0 || !nopartialsmall) ? n_partial : 0);
    ST_double small_value = 0, has_cons_value = 0;
    SF_scal_use("__civreghdfe_small", &small_value);
    SF_scal_use("__civreghdfe_has_cons", &has_cons_value);
    ST_int small = G > 0 || small_value != 0;
    ST_int dofminus_vce = dofminus;
    ST_int df_a_small = df_a_for_vce + sdof_extra;
    if (sdof_extra > 0 && df_a_for_vce == 0 && nested_adj == 1) df_a_small = 1 + sdof_extra;

    ST_double lambda = 1.0;
    if (!rc) rc = ivest_compute_2sls(
        y_dem, X_exog_dem, X_endog_dem, Z_dem,
        weights_c, weight_type,
        N, N_eff, K_exog, K_endog, K_iv,
        beta, V, first_stage_F,
        vce_type, cluster_ids_c, num_clusters,
        cluster2_ids_c, num_clusters2,
        df_a_small, nested_adj, verbose,
        est_method, kclass_user, fuller_alpha, &lambda,
        kernel_type, bw, kiefer,
        ts_mode > 0 ? &ts_info : NULL,
        (G > 0) ? (df_a_for_vce > 0 ? df_a_for_vce : 1) + sdof_extra : sdofminus_opt, center,
        &test_cols,
        n_abs_x_exog + n_abs_x_endog, n_abs_x_exog + n_abs_z_excl, n_abs_exog_all,
        dofminus_vce, small
    );
    free(test_idx);
    free(ts_int);
    test_idx = ts_int = NULL;

    if (rc != STATA_OK) {
        iv_rc = rc;
        goto iv_cleanup;
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
    y_mean = (G > 0 || has_cons_value != 0) ? y_mean / sum_w : 0.0;

    /* Compute TSS (centered) = sum(w * (y - ybar)^2) */
    for (ST_int i = 0; i < N; i++) {
        ST_double w = (weights_c && weight_type != 0) ? weights_c[i] : 1.0;
        ST_double dev = y_dem[i] - y_mean;
        tss += w * dev * dev;
    }

    ST_double resid_idx_value = 0;
    SF_scal_use("__civreghdfe_resid_idx", &resid_idx_value);
    ST_int resid_idx = (ST_int)resid_idx_value;
    ST_retcode residual_rc = 0;
    /* Residuals y - (X_exog * beta[0:K_exog-1] + X_endog * beta[K_exog:K_total-1]),
       each fitted value accumulated in that column order, then
       RSS = sum(w * resid^2) in row order */
    ST_double *fit = (ST_double *)ctools_safe_calloc2((size_t)N, sizeof(ST_double));
    if (!fit) {
        SF_error("civreghdfe: Memory allocation failed for residuals\n");
        iv_rc = 920;
        goto iv_cleanup;
    }
    if (K_exog > 0) ctools_matvec_add(X_exog_dem, N, K_exog, beta, fit);
    if (K_endog > 0) ctools_matvec_add(X_endog_dem, N, K_endog, beta + K_exog, fit);
    for (ST_int i = 0; i < N; i++) {
        ST_double resid = y_dem[i] - fit[i];
        if (resid_idx > 0 && !residual_rc)
            residual_rc = (_stata_)->safestore(resid_idx, (ST_int)est_obs[i], resid);
        ST_double w = (weights_c && weight_type != 0) ? weights_c[i] : 1.0;
        rss += w * resid * resid;
    }
    free(fit);

    ST_double r2 = (tss > 0) ? 1.0 - rss / tss : 0.0;

    /* Compute df_r:
       - For clustered VCE (vce_type == 2, 3) and Driscoll-Kraay (4, clusters
         = the time periods of the estimation sample): df_r = num_clusters - 1
       - Otherwise: df_r = N - K - df_a - dofminus
       - nopartialsmall: exclude partialled variables from K for DOF calc */
    /* FE-absorbed columns kept as noise by ivreghdfe (stage 1 above) count */
    ST_int n_abs_x = n_abs_x_exog + n_abs_x_endog;
    ST_int K_for_dof = (G > 0) ? K_total + n_abs_x + sdof_extra
                               : K_total + n_abs_x + sdof_extra;
    /* In double: N_eff (the sum of fweights) may exceed INT_MAX */
    ST_double df_r_val;
    if (kiefer) {
        /* Kiefer VCE: df_r = N_eff - K - df_a - dofminus (not cluster-based) */
        df_r_val = N_eff - K_for_dof - df_a - dofminus;
    } else if (has_cluster2) {
        /* Two-way clustering: df_r = min(G1, G2) - 1 */
        ST_int min_clust = (num_clusters < num_clusters2) ? num_clusters : num_clusters2;
        df_r_val = min_clust - 1;
    } else if (has_cluster) {
        df_r_val = num_clusters - 1;
    } else {
        df_r_val = N_eff - K_for_dof - df_a - dofminus;
    }
    if (df_r_val <= 0) df_r_val = 1;
    /* For rmse, use sdofminus = max(1, df_a_for_vce) when nested, df_a otherwise
       (matches ivreghdfe line 670: if HDFE.df_a=0, force absorb_ct to 1) */
    ST_int sdofminus = (has_cluster && df_a_nested > 0) ?
                       (df_a_for_vce > 0 ? df_a_for_vce : 1) : df_a;
    ST_double rmse_den = (ST_double)N_eff - K_total - n_abs_x - sdofminus
                         - sdof_extra - dofminus_vce;
    if (!small) rmse_den = N_eff - dofminus;
    ST_double rmse = sqrt(rss / (rmse_den > 0 ? rmse_den : 1));

    /* End stats timing, start store timing */
    double t_stats = ctools_timer_seconds() - t_stats_start;
    double t_store_start = ctools_timer_seconds();

    const civreghdfe_result result = {
        .df_r_val = df_r_val,
        .df_a = df_a,
        .df_a_for_vce = df_a_for_vce,
        .K_total = K_total,
        .n_abs_x = n_abs_x,
        .has_cluster2 = has_cluster2,
        .num_clusters = num_clusters,
        .num_clusters2 = num_clusters2,
        .has_cluster = has_cluster,
        .est_method = est_method,
        .K_exog_orig = K_exog_orig,
        .K_endog_orig = K_endog_orig,
        .K_endog = K_endog,
        .N_eff = N_eff,
        .rss = rss,
        .tss = tss,
        .r2 = r2,
        .rmse = rmse,
        .lambda = lambda,
        .beta = beta,
        .V = V,
        .first_stage_F = first_stage_F,
        .is_collinear_x = is_collinear_x,
        .residual_rc = residual_rc,
    };
    iv_rc = civreghdfe_post_results(&result);
    if (iv_rc) goto iv_cleanup;

    /* End store timing, compute total */
    t_store = ctools_timer_seconds() - t_store_start;
    double t_total = ctools_timer_seconds() - t_total_start;

    /* Save timing scalars (checked, as every other result) */
    {
        const char *const time_names[] = {
            "_civreghdfe_time_load", "_civreghdfe_time_extract", "_civreghdfe_time_missing",
            "_civreghdfe_time_singleton", "_civreghdfe_time_remap", "_civreghdfe_time_fwl",
            "_civreghdfe_time_partial", "_civreghdfe_time_dof", "_civreghdfe_time_estimate",
            "_civreghdfe_time_stats", "_civreghdfe_time_store", "_civreghdfe_time_total"};
        const double time_values[] = {t_load, t_extract, t_missing, t_singleton,
                                      t_hdfe_setup, t_fwl, t_partial_out, t_postproc,
                                      t_estimate, t_stats, t_store, t_total};
        for (size_t t = 0; t < sizeof(time_values) / sizeof(time_values[0]); t++) {
            if ((iv_rc = ctools_scal_save(time_names[t], time_values[t])))
                goto iv_cleanup;
        }
    }
    CTOOLS_SAVE_THREAD_INFO("_civreghdfe");

iv_cleanup:
    /* Single exit: releases whatever the stage reached had allocated,
       including after a failed result publication */
    free(cluster_ids); free(cluster2_ids); free(valid_mask); free(singleton_mask);
    if (fe_levels_c) {
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
    }
    free(level_max); free(is_partial); free(is_collinear_x); free(test_idx); free(ts_int);
    free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
    free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
    free(beta); free(V); free(first_stage_F);
    if (state) {
        ctools_hdfe_state_cleanup(state);
        free(state);
    }
    iv_state = NULL;
    ctools_aligned_free(est_obs);
    ctools_filtered_data_free(&filtered);

    return iv_rc;
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
 * Cleanup function for civreghdfe persistent state.
 * Frees the global HDFE state if allocated.
 * Safe to call multiple times (idempotent).
 */
void civreghdfe_cleanup_state(void)
{
    /* civreghdfe cleans up iv_state at the end of do_iv_regression(),
       so this is just a safety measure for interrupted execution */
    if (iv_state != NULL) {
        /* Note: Full cleanup would require knowing G and num_threads,
           which we don't have here. The main function handles full cleanup.
           This just nulls the pointer to prevent double-free issues. */
        iv_state = NULL;
    }
}

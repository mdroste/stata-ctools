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


/* Expand omitted columns and reorder [exog,endog] to Stata's [endog,exog].
 * A failed map allocation must not fall back to indexing the compact V with
 * full-model indices. Checked SPI writes also protect undersized host matrices. */
static ST_retcode store_iv_matrices(const ST_double *beta, const ST_double *V,
                                   ST_int K_exog_orig, ST_int K_endog_orig,
                                   ST_int K_total, const ST_int *is_collinear)
{
    ST_int full_K = K_exog_orig + K_endog_orig;
    ST_int *map = (ST_int *)malloc((size_t)full_K * sizeof(ST_int));
    if (!map) return 920;
    ST_int compact = 0;
    for (ST_int k = 0; k < full_K; k++)
        map[k] = is_collinear[k] ? -1 : compact++;
    if (compact != K_total) {
        free(map);
        return 503;
    }

    ST_retcode rc = 0;
    for (ST_int i = 0; i < full_K && !rc; i++) {
        ST_int fi = i < K_endog_orig ? K_exog_orig + i : i - K_endog_orig;
        ST_int ci = map[fi];
        rc = ctools_mat_store("__civreghdfe_b", 1, i + 1, ci < 0 ? 0.0 : beta[ci]);
        for (ST_int j = 0; j < full_K && !rc; j++) {
            ST_int fj = j < K_endog_orig ? K_exog_orig + j : j - K_endog_orig;
            ST_int cj = map[fj];
            ST_double value = ci < 0 || cj < 0 ? 0.0 : V[(size_t)cj * K_total + ci];
            rc = ctools_mat_store("__civreghdfe_V", i + 1, j + 1, value);
        }
    }
    free(map);
    return rc;
}

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
            ctools_scal_save("__civreghdfe_dropped_opt", (ST_double)opt_code);
            ctools_scal_save("__civreghdfe_dropped_col", (ST_double)(j + 1));
            free(flags);
            return 198;
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
        for (ST_int a = 0; a < n_keep; a++) {
            ST_double d = PtP[a * n_keep + a];
            if (d > 0.0) {
                ST_double s = 1.0 / sqrt(d);
                for (ST_int b = 0; b < n_keep; b++) {
                    PtP[b * n_keep + a] *= s;
                    PtP[a * n_keep + b] *= s;
                }
            }
        }
        if (detect_collinearity(PtP, n_keep, coll, 0) < 0) {
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
    /* nested_fe_index from .ado is no longer used — data-based detection in C handles it */

    /* New scalars for estimation method */
    ST_int est_method = 0;
    ST_double kclass_user = 0, fuller_alpha = 0;
    SF_scal_use("__civreghdfe_est_method", &dval); est_method = (ST_int)dval;
    SF_scal_use("__civreghdfe_kclass", &dval); kclass_user = dval;
    SF_scal_use("__civreghdfe_fuller", &dval); fuller_alpha = dval;

    /* HAC parameters */
    ST_int kernel_type = 0, bw = 0, kiefer = 0;
    SF_scal_use("__civreghdfe_kernel", &dval); kernel_type = (ST_int)dval;
    SF_scal_use("__civreghdfe_bw", &dval); bw = (ST_int)dval;
    SF_scal_use("__civreghdfe_kiefer", &dval); kiefer = (ST_int)dval;

    /* tsset structure for the kernel estimators (HAC, AC, Driscoll-Kraay,
       Kiefer): ts_mode 1 = the time variable, 2 = the panel and time
       variables, passed after the weight; tdelta = tsset delta */
    ST_int ts_mode = 0;
    ST_double tdelta = 1.0;
    if (SF_scal_use("__civreghdfe_ts_mode", &dval) == 0 && dval > 0) ts_mode = (dval >= 2) ? 2 : 1;
    if (SF_scal_use("__civreghdfe_tdelta", &dval) == 0 && dval > 0 && !SF_is_missing(dval)) tdelta = dval;
    ST_int ts_cols = (ts_mode == 2) ? 2 : ts_mode;

    /* DOF adjustment parameters */
    ST_int dofminus = 0, sdofminus_opt = 0, nopartialsmall = 0, center = 0;
    SF_scal_use("__civreghdfe_dofminus", &dval); dofminus = (ST_int)dval;
    SF_scal_use("__civreghdfe_sdofminus", &dval); sdofminus_opt = (ST_int)dval;
    SF_scal_use("__civreghdfe_nopartialsmall", &dval); nopartialsmall = (ST_int)dval;
    SF_scal_use("__civreghdfe_center", &dval); center = (ST_int)dval;

    if (center) {
        SF_error("civreghdfe: center is not supported\n");
        return 198;
    }

    /* Calculate variable positions */
    /* Layout: [y, X_endog (Ke), X_exog (Kx), Z (Kz), FE (G), cluster?, cluster2?, weight?,
                tsset panel?, tsset time?] */
    ST_int var_y_idx = 0;  /* Index in loaded_data.vars[] */
    ST_int var_endog_start_idx = 1;
    ST_int var_exog_start_idx = var_endog_start_idx + K_endog;
    ST_int var_iv_start_idx = var_exog_start_idx + K_exog;
    ST_int var_fe_start_idx = var_iv_start_idx + K_iv;
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
    ST_int total_vars = 1 + K_endog + K_exog + K_iv + G +
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

    /* Allocate cluster ID arrays (needed for string conversion + sentinel marking).
     * All other data (y, X, Z, FE, weights) is read directly from filtered.data
     * during missing check and compaction, avoiding redundant N_total-sized copies. */
    ST_int *cluster_ids = has_cluster ? (ST_int *)ctools_safe_malloc2((size_t)N_total, sizeof(ST_int)) : NULL;
    ST_int *cluster2_ids = has_cluster2 ? (ST_int *)ctools_safe_malloc2((size_t)N_total, sizeof(ST_int)) : NULL;

    if ((has_cluster && !cluster_ids) || (has_cluster2 && !cluster2_ids)) {
        SF_error("civreghdfe: Memory allocation failed\n");
        free(cluster_ids); free(cluster2_ids);
        ctools_filtered_data_free(&filtered);
        return 920;
    }

    /* Extract cluster IDs (needs early processing for string variables) */
    ST_int i;
    if (has_cluster) {
        if (filtered.data.vars[var_cluster_idx].type == STATA_TYPE_STRING) {
            ST_int n_groups = 0;
            if (ctools_strings_to_cluster_ids(filtered.data.vars[var_cluster_idx].data.str,
                                               N_total, cluster_ids, &n_groups) != 0) {
                SF_error("civreghdfe: Failed to convert string cluster variable to group IDs\n");
                free(cluster_ids); free(cluster2_ids);
                ctools_filtered_data_free(&filtered);
                return 920;
            }
        } else {
            double *src = filtered.data.vars[var_cluster_idx].data.dbl;
            ST_int n_groups = 0;
            if (ctools_numeric_to_cluster_ids(src, N_total, cluster_ids, &n_groups) != 0) {
                free(cluster_ids); free(cluster2_ids);
                ctools_filtered_data_free(&filtered);
                return 920;
            }
        }
    }
    if (has_cluster2) {
        if (filtered.data.vars[var_cluster2_idx].type == STATA_TYPE_STRING) {
            ST_int n_groups = 0;
            if (ctools_strings_to_cluster_ids(filtered.data.vars[var_cluster2_idx].data.str,
                                               N_total, cluster2_ids, &n_groups) != 0) {
                SF_error("civreghdfe: Failed to convert string cluster2 variable to group IDs\n");
                free(cluster_ids); free(cluster2_ids);
                ctools_filtered_data_free(&filtered);
                return 920;
            }
        } else {
            double *src = filtered.data.vars[var_cluster2_idx].data.dbl;
            ST_int n_groups = 0;
            if (ctools_numeric_to_cluster_ids(src, N_total, cluster2_ids, &n_groups) != 0) {
                free(cluster_ids); free(cluster2_ids);
                ctools_filtered_data_free(&filtered);
                return 920;
            }
        }
    }

    /* End extraction, start missing value check timing.
     * Data arrays (y, X, Z, FE, weights) are NOT copied — we read directly
     * from filtered.data during missing value check and compaction. */
    t_extract = ctools_timer_seconds() - t_extract_start;
    double t_missing_start = ctools_timer_seconds();

    /* Drop observations with missing values.
     * Read directly from filtered.data to avoid redundant N_total-sized copies. */
    ST_int *valid_mask = (ST_int *)calloc(N_total, sizeof(ST_int));
    if (!valid_mask) {
        SF_error("civreghdfe: Memory allocation failed for valid_mask\n");
        free(cluster_ids); free(cluster2_ids);
        ctools_filtered_data_free(&filtered);
        return 920;
    }
    ST_int N_valid = 0;

    for (i = 0; i < N_total; i++) {
        int is_valid = 1;

        /* Check y */
        if (SF_is_missing(filtered.data.vars[var_y_idx].data.dbl[i])) is_valid = 0;

        /* Check X_endog */
        for (ST_int k = 0; k < K_endog && is_valid; k++) {
            if (SF_is_missing(filtered.data.vars[var_endog_start_idx + k].data.dbl[i])) is_valid = 0;
        }

        /* Check X_exog */
        for (ST_int k = 0; k < K_exog && is_valid; k++) {
            if (SF_is_missing(filtered.data.vars[var_exog_start_idx + k].data.dbl[i])) is_valid = 0;
        }

        /* Check Z */
        for (ST_int k = 0; k < K_iv && is_valid; k++) {
            if (SF_is_missing(filtered.data.vars[var_iv_start_idx + k].data.dbl[i])) is_valid = 0;
        }

        /* Check FE variables for missing values */
        for (ST_int g = 0; g < G && is_valid; g++) {
            if (SF_is_missing(filtered.data.vars[var_fe_start_idx + g].data.dbl[i])) is_valid = 0;
        }

        /* Check weights */
        if (has_weights && is_valid) {
            double w = filtered.data.vars[var_weight_idx].data.dbl[i];
            if (SF_is_missing(w) || w <= 0) is_valid = 0;
        }

        /* Check cluster variables for missing values (marked as -1 during extraction) */
        if (has_cluster && is_valid) {
            if (cluster_ids[i] < 0) is_valid = 0;
        }
        if (has_cluster2 && is_valid) {
            if (cluster2_ids[i] < 0) is_valid = 0;
        }

        /* Check the tsset panel and time variables */
        for (ST_int k = 0; k < ts_cols && is_valid; k++) {
            if (SF_is_missing(filtered.data.vars[var_ts_start_idx + k].data.dbl[i])) is_valid = 0;
        }

        valid_mask[i] = is_valid;
        if (is_valid) N_valid++;
    }

    if (N_valid < K_total + K_iv + 1) {
        SF_error("civreghdfe: Insufficient observations\n");
        free(cluster_ids); free(cluster2_ids); free(valid_mask);
        ctools_filtered_data_free(&filtered);
        return 2001;
    }

    /* Compact FE levels from raw doubles with O(N) counting-based remap.
     * This replaces the old approach of: remap_values_sorted (O(N log N)) +
     * separate integer compaction. Now we compact and remap in one step. */
    ST_int **fe_levels_c = (ST_int **)calloc((size_t)G, sizeof(ST_int *));
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
        free(cluster_ids); free(cluster2_ids); free(valid_mask);
        ctools_filtered_data_free(&filtered);
        return 920;
    }
    {
        /* Compact FE doubles to temp array, then O(N) remap to contiguous integers */
        double *fe_compact = (double *)malloc((size_t)N_valid * sizeof(double));
        if (!fe_compact) {
            for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
            free(fe_levels_c);
            free(cluster_ids); free(cluster2_ids); free(valid_mask);
            ctools_filtered_data_free(&filtered);
            return 920;
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
                for (ST_int fg = 0; fg < G; fg++) free(fe_levels_c[fg]);
                free(fe_levels_c);
                free(cluster_ids); free(cluster2_ids); free(valid_mask);
                ctools_filtered_data_free(&filtered);
                return 920;
            }
            free(counts_g);  /* Counts recomputed later in HDFE setup */
        }
        free(fe_compact);
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
        free(cluster_ids); free(cluster2_ids); free(valid_mask);
        ctools_filtered_data_free(&filtered);
        return 920;
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
        free(singleton_mask);
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
        free(cluster_ids); free(cluster2_ids); free(valid_mask);
        ctools_filtered_data_free(&filtered);
        return 920;
    }

    /* Step 3: Build combined index mapping original obs → final position.
     * Combine valid_mask and singleton_mask into one compaction pass. */
    ST_int N = N_valid - num_singletons_total;

    /* Allocate final arrays */
    ST_double *y_c = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double));
    ST_double *X_endog_c = (K_endog > 0) ? (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_endog, sizeof(ST_double)) : NULL;
    ST_double *X_exog_c = (K_exog > 0) ? (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_exog, sizeof(ST_double)) : NULL;
    /* Z_c holds the K_iv instruments followed by the ts_cols tsset columns
       (panel, time), which share its lifetime */
    ST_double *Z_c = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)(K_iv + ts_cols), sizeof(ST_double));
    ST_double *ts_c = Z_c ? Z_c + (size_t)N * K_iv : NULL;
    ST_double *weights_c = has_weights ? (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double)) : NULL;
    ST_int *cluster_ids_c = has_cluster ? (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int)) : NULL;
    ST_int *cluster2_ids_c = has_cluster2 ? (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int)) : NULL;

    if (!y_c || !Z_c || (K_endog > 0 && !X_endog_c) || (K_exog > 0 && !X_exog_c) ||
        (has_weights && !weights_c) || (has_cluster && !cluster_ids_c) ||
        (has_cluster2 && !cluster2_ids_c)) {
        SF_error("civreghdfe: Memory allocation failed for compacted arrays\n");
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        free(singleton_mask);
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
        free(cluster_ids); free(cluster2_ids); free(valid_mask);
        ctools_filtered_data_free(&filtered);
        return 920;
    }

    /* Single fused compaction pass: skip invalid AND singleton observations.
     * Read directly from filtered.data — no intermediate N_total-sized copies. */
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
            filtered.obs_map[out_idx] = filtered.obs_map[ii];
            y_c[out_idx] = filtered.data.vars[var_y_idx].data.dbl[ii];
            for (ST_int k = 0; k < K_endog; k++)
                X_endog_c[(size_t)k * N + out_idx] = filtered.data.vars[var_endog_start_idx + k].data.dbl[ii];
            for (ST_int k = 0; k < K_exog; k++)
                X_exog_c[(size_t)k * N + out_idx] = filtered.data.vars[var_exog_start_idx + k].data.dbl[ii];
            for (ST_int k = 0; k < K_iv; k++)
                Z_c[(size_t)k * N + out_idx] = filtered.data.vars[var_iv_start_idx + k].data.dbl[ii];
            for (ST_int k = 0; k < ts_cols; k++)
                ts_c[(size_t)k * N + out_idx] = filtered.data.vars[var_ts_start_idx + k].data.dbl[ii];
            if (has_weights) weights_c[out_idx] = filtered.data.vars[var_weight_idx].data.dbl[ii];
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
    perm_idx_t *est_obs = filtered.obs_map;
    filtered.obs_map = NULL;

    /* Free original arrays and filtered data (kept alive for direct reads during compaction) */
    free(cluster_ids); free(cluster2_ids); free(valid_mask);
    ctools_filtered_data_free(&filtered);

    /* Check we still have enough observations */
    if (N < K_total + K_iv + 1) {
        SF_error("civreghdfe: Insufficient observations after singleton removal\n");
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
        ctools_aligned_free(est_obs);
        return 2001;
    }

    /* End singleton timing, start HDFE setup timing */
    t_singleton = ctools_timer_seconds() - t_singleton_start;
    double t_hdfe_setup_start = ctools_timer_seconds();

    /* Initialize HDFE state (reuse creghdfe infrastructure) */
    /* This requires setting up FE_Factor structures */
    HDFE_State *state = (HDFE_State *)calloc(1, sizeof(HDFE_State));
    if (!state) {
        SF_error("civreghdfe: Memory allocation failed for HDFE state\n");
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
        ctools_aligned_free(est_obs);
        return 920;
    }
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
    if (!state->factors) {
        SF_error("civreghdfe: Memory allocation failed for FE factors\n");
        free(state);
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
        ctools_aligned_free(est_obs);
        return 920;
    }
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
            ctools_aligned_free(est_obs);
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
    }

    /* Allocate inv_counts, inv_weighted_counts, and CG solver buffers */
    ST_int num_threads = ctools_get_max_threads();
    if (num_threads > 8) num_threads = 8;
    state->num_threads = num_threads;

    /* Max columns partialled out: y + endogenous + exogenous + instruments */
    ST_int max_partial_cols = 1 + K_endog + K_exog + K_iv;
    if (ctools_hdfe_alloc_buffers(state, 0, max_partial_cols) != 0) {
        SF_error("civreghdfe: Memory allocation failed for HDFE buffers\n");
        ctools_hdfe_state_cleanup(state);
        free(state);
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        free(fe_levels_c);
        ctools_aligned_free(est_obs);
        return 920;
    }

    state->factors_initialized = 1;

    /* Set global state for creghdfe solver functions */
    g_state = state;

    /* End HDFE setup timing, start FWL timing */
    t_hdfe_setup = ctools_timer_seconds() - t_hdfe_setup_start;
    double t_fwl_start = ctools_timer_seconds();

    /* FWL partialling: partial() arrives as one 0/1 flag per X_exog column
       (matrix __civreghdfe_partial_flags, aligned with the expanded exogenous
       regressors). is_partial: 1 = partialled out, 2 = dropped as redundant. */
    ST_double dval_partial;
    ST_int n_partial = 0;
    ST_int *is_partial = NULL;

    if (SF_scal_use("__civreghdfe_n_partial", &dval_partial) == 0) {
        n_partial = (ST_int)dval_partial;
    }

    if (n_partial > 0 && K_exog > 0) {
        is_partial = (ST_int *)calloc(K_exog, sizeof(ST_int));
        ST_int n_flagged = is_partial ? read_flag_matrix("__civreghdfe_partial_flags", K_exog, is_partial) : -1;
        if (n_flagged != n_partial) {
            ST_retcode flag_rc = is_partial ? 198 : 920;
            SF_error(is_partial ? "civreghdfe: partial() column flags do not match the exogenous regressors\n"
                                : "civreghdfe: Memory allocation failed for partial indices\n");
            free(is_partial);
            free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
            free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
            free(fe_levels_c);
            ctools_hdfe_state_cleanup(state);
            free(state);
            g_state = NULL;
            ctools_aligned_free(est_obs);
            return flag_rc;
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
    ST_double *all_data = (ST_double *)ctools_safe_malloc3((size_t)N + 1, (size_t)total_cols, sizeof(ST_double));
    if (!all_data) {
        SF_error("civreghdfe: Memory allocation failed for demeaning buffer\n");
        free(is_partial);
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        free(fe_levels_c);  /* levels freed by state cleanup */
        ctools_hdfe_state_cleanup(state);
        free(state);
        g_state = NULL;
        ctools_aligned_free(est_obs);
        return 920;
    }

    /* Copy y */
    memcpy(all_data, y_c, N * sizeof(ST_double));

    /* Copy X_endog */
    if (K_endog > 0) {
        memcpy(all_data + N, X_endog_c, (size_t)N * K_endog * sizeof(ST_double));
    }

    /* Copy X_exog */
    if (K_exog > 0) {
        memcpy(all_data + (size_t)N * (1 + K_endog), X_exog_c, (size_t)N * K_exog * sizeof(ST_double));
    }

    /* Copy Z */
    memcpy(all_data + (size_t)N * (1 + K_endog + K_exog), Z_c, (size_t)N * K_iv * sizeof(ST_double));

    /* Pre-partialling (weighted) sum of squared deviations of every column.
       A column whose partialled SS is at most absorb_tol times this was
       absorbed by the fixed effects (or by partial()) and is omitted, as
       ivreghdfe does; absorb_tol is reghdfe's collinear-with-FE tolerance. */
    const ST_double absorb_tol = (tolerance / 10.0 < 1e-6) ? tolerance / 10.0 : 1e-6;
    ST_double *ss_pre = all_data + (size_t)N * total_cols;
    #pragma omp parallel for schedule(dynamic) if ((size_t)N * total_cols > 100000)
    for (ST_int c = 0; c < total_cols; c++) {
        ss_pre[c] = civ_wss_dev(all_data + (size_t)c * N, weights_c, N);
    }
    const ST_double *ss_endog = ss_pre + 1;
    ST_double *ss_exog = ss_pre + 1 + K_endog;  /* compacted with X_exog below */
    const ST_double *ss_excl = ss_pre + 1 + K_endog + 2 * K_exog;  /* after Z's copy of X_exog */

    /* Partial out the FEs from every column in one pass, then the partial()
       variables by weighted FWL on the FE-demeaned data: the residual of the
       W-weighted projection on [FE, P] is M_P~ M_FE x with P~ = M_FE P, the
       way ivreghdfe runs HDFE partial_out and then s_partial. */
    HDFE_SolveResult projection = partial_out_columns(state, all_data, N, total_cols, num_threads);
    ST_retcode fwl_rc = 0;
    if (!projection.status && is_partial) {
        fwl_rc = fwl_partial_out(all_data, N, total_cols, 1 + K_endog, 1 + K_endog + K_exog, K_exog,
                                 is_partial, weights_c, ss_exog, absorb_tol);
    }

    if (projection.status || fwl_rc) {
        if (projection.status) SF_error("civreghdfe: fixed-effect projection did not converge\n");
        free(all_data); free(is_partial);
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        free(fe_levels_c); ctools_aligned_free(est_obs);
        ctools_hdfe_state_cleanup(state); free(state); g_state = NULL;
        return projection.status ? projection.status : fwl_rc;
    }

    /* End partial out timing, start post-processing timing */
    t_partial_out = ctools_timer_seconds() - t_partial_start;
    double t_postproc_start = ctools_timer_seconds();

    /* Extract demeaned data */
    ST_double *y_dem = all_data;
    ST_double *X_endog_dem = (K_endog > 0) ? all_data + N : NULL;
    ST_double *X_exog_dem = (K_exog > 0) ? all_data + (size_t)N * (1 + K_endog) : NULL;
    ST_double *Z_dem = all_data + (size_t)N * (1 + K_endog + K_exog);

    /* Remove partial variables from X_exog AND from Z (instruments) after convergence.
       Z contains [exog_vars, excluded_instruments], so partial exog vars must be removed
       from both X_exog and the first K_exog columns of Z. Partial columns dropped as
       redundant (flag 2) are removed as well, and reported back for the ado's note. */
    if (is_partial) {
        ST_int K_exog_new = K_exog - n_partial;
        ST_int K_excl = K_iv - K_exog;  /* Number of excluded instruments */

        for (ST_int k = 0; k < K_exog; k++) {
            if (is_partial[k] == 2)
                ctools_mat_store("__civreghdfe_partial_flags", 1, k + 1, 2.0);
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
        X_exog_dem = (K_exog > 0) ? all_data + (size_t)N * (1 + K_endog) : NULL;

        free(is_partial);
        is_partial = NULL;
    }

    /* Compute df_a (absorbed degrees of freedom) and mobility groups */
    ST_int df_a = 0;
    ST_int df_a_nested = 0;  /* Levels from FE nested within cluster */
    ST_int mobility_groups = 0;

    ctools_compute_hdfe_dof(state->factors, G, N, &df_a, &mobility_groups);

    for (ST_int g = 0; g < G; g++) {
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
            ctools_aligned_free(est_obs);
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
            ctools_aligned_free(est_obs);
            return 920;
        }
    }

    if (!kiefer && (vce_type == 2 || vce_type == 3) &&
        (num_clusters < 2 || (has_cluster2 && num_clusters2 < 2))) {
        SF_error("civreghdfe: clustered VCE requires at least two retained clusters in each dimension\n");
        free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        free(fe_levels_c);
        ctools_hdfe_state_cleanup(state);
        free(state); g_state = NULL;
        ctools_aligned_free(est_obs);
        return 459;
    }

    /* Detect FEs nested within cluster variable using data-based check.
     * This works for both numeric and string cluster variables. */
    ST_int *fe_nested = (ST_int *)calloc(G, sizeof(ST_int));
    if (has_cluster && fe_nested) {
        df_a_nested = 0;
        for (ST_int g = 0; g < G; g++) {
            ST_int is_nested = ctools_fe_nested_in_cluster(
                state->factors[g].levels, state->factors[g].num_levels,
                cluster_ids_c, N);
            if (is_nested < 0) is_nested = 0;
            fe_nested[g] = is_nested;
            if (is_nested) df_a_nested += state->factors[g].num_levels;

            char scalar_name[64];
            snprintf(scalar_name, sizeof(scalar_name), "__civreghdfe_fe_nested_%d", (int)(g + 1));
            ctools_scal_save(scalar_name, (ST_double)is_nested);
        }
        /* Recompute VCE DOF with data-based nesting results */
        df_a_for_vce = df_a - df_a_nested;
        nested_adj = (df_a_nested > 0) ? 1 : 0;
    }
    if (fe_nested) free(fe_nested);

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
    ST_int *is_collinear_x = (ST_int *)calloc(K_total, sizeof(ST_int));
    ST_int num_collinear_x = 0;

    if (!is_collinear_x) {
        SF_error("civreghdfe: Memory allocation failed for collinearity arrays\n");
        free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        ctools_hdfe_state_cleanup(state);
        free(state);
        g_state = NULL;
        ctools_aligned_free(est_obs);
        return 920;
    }

    /* Stage 1: columns absorbed by the fixed effects (or by partial()): the
       partialled (weighted) SS is zero or negligible relative to the
       pre-partialling SS (reghdfe's collinear-with-FE rule). What remains is
       numerical noise, which ivreghdfe omits; it must not reach stage 2,
       whose equilibration would rescale it to unit variance. */
    {
        for (ST_int k = 0; k < K_exog; k++) {
            ST_double xx_partial = civ_wss(X_exog_dem + (size_t)k * N, weights_c, N);
            if (xx_partial < 1e-30 || xx_partial <= absorb_tol * ss_exog[k]) {
                is_collinear_x[k] = 1;
                num_collinear_x++;
            }
        }
        for (ST_int k = 0; k < K_endog; k++) {
            ST_double xx_partial = civ_wss(X_endog_dem + (size_t)k * N, weights_c, N);
            if (xx_partial < 1e-30 || xx_partial <= absorb_tol * ss_endog[k]) {
                is_collinear_x[K_exog + k] = 1;
                num_collinear_x++;
            }
        }
    }

    /* Stage 2: Numerical collinearity via Cholesky on X'X
       Build concatenated X = [X_exog, X_endog] in column-major order */
    if (K_total > 1) {
        ST_double *XtX = (ST_double *)ctools_safe_calloc3((size_t)K_total, (size_t)K_total, sizeof(ST_double));
        if (!XtX) {
            free(is_collinear_x);
            SF_error("civreghdfe: Memory allocation failed for XtX\n");
            free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
            free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
            ctools_hdfe_state_cleanup(state);
            free(state);
            g_state = NULL;
            ctools_aligned_free(est_obs);
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
            ST_double *blk_ee = (ST_double *)ctools_safe_malloc3((size_t)K_exog, (size_t)K_exog, sizeof(ST_double));
            ST_double *blk_en = (ST_double *)ctools_safe_malloc3((size_t)K_exog, (size_t)K_endog, sizeof(ST_double));
            ST_double *blk_nn = (ST_double *)ctools_safe_malloc3((size_t)K_endog, (size_t)K_endog, sizeof(ST_double));
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
            ctools_matmul_atb(X_exog_dem, X_exog_dem, N, K_exog, K_exog, XtX);
        } else {
            ctools_matmul_atb(X_endog_dem, X_endog_dem, N, K_endog, K_endog, XtX);
        }

        /* Stage-1 columns get zero rows and columns, so the equilibration
         * leaves them at zero and the Cholesky below flags them without
         * involving the other columns. detect_collinearity() resets the flag
         * array; the stage-1 flags are restored after it. */
        ST_int *stage1_x = (ST_int *)malloc((size_t)K_total * sizeof(ST_int));
        if (stage1_x) memcpy(stage1_x, is_collinear_x, (size_t)K_total * sizeof(ST_int));
        for (ST_int ii = 0; ii < K_total; ii++) {
            if (!is_collinear_x[ii]) continue;
            for (ST_int jj = 0; jj < K_total; jj++) {
                XtX[jj * K_total + ii] = 0.0;
                XtX[ii * K_total + jj] = 0.0;
            }
        }

        /* Equilibrate XtX to correlation matrix before collinearity detection.
         * This makes the Cholesky-based detection scale-invariant, preventing
         * false positives from mixed-scale variables (e.g., x_tiny / x_huge). */
        {
            for (ST_int ii = 0; ii < K_total; ii++) {
                ST_double d = XtX[ii * K_total + ii];
                if (d > 0.0) {
                    ST_double s = 1.0 / sqrt(d);
                    for (ST_int jj = 0; jj < K_total; jj++) {
                        XtX[jj * K_total + ii] *= s;
                        XtX[ii * K_total + jj] *= s;
                    }
                }
            }
        }

        ST_int num_chol_collinear = detect_collinearity(XtX, K_total, is_collinear_x, verbose);
        free(XtX);
        if (stage1_x) {
            for (ST_int k = 0; k < K_total; k++) if (stage1_x[k]) is_collinear_x[k] = 1;
            free(stage1_x);
        }

        if (num_chol_collinear < 0) {
            free(is_collinear_x);
            SF_error("civreghdfe: Collinearity detection failed\n");
            free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
            free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
            ctools_hdfe_state_cleanup(state);
            free(state);
            g_state = NULL;
            ctools_aligned_free(est_obs);
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
            }
        }

        ST_double *ZtZ = (K_iv > 1) ? (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_iv, sizeof(ST_double)) : NULL;
        ST_int *stage1_z = (K_iv > 1) ? (ST_int *)malloc((size_t)K_iv * sizeof(ST_int)) : NULL;
        if (ZtZ && stage1_z) {
            memcpy(stage1_z, is_collinear_z, (size_t)K_iv * sizeof(ST_int));
            ctools_matmul_atb(Z_dem, Z_dem, N, K_iv, K_iv, ZtZ);
            for (ST_int ii = 0; ii < K_iv; ii++) {
                if (!stage1_z[ii]) continue;
                for (ST_int jj = 0; jj < K_iv; jj++) {
                    ZtZ[jj * K_iv + ii] = 0.0;
                    ZtZ[ii * K_iv + jj] = 0.0;
                }
            }

            /* Equilibrate ZtZ for scale-invariant collinearity detection */
            for (ST_int ii = 0; ii < K_iv; ii++) {
                ST_double d = ZtZ[ii * K_iv + ii];
                if (d > 0.0) {
                    ST_double s = 1.0 / sqrt(d);
                    for (ST_int jj = 0; jj < K_iv; jj++) {
                        ZtZ[jj * K_iv + ii] *= s;
                        ZtZ[ii * K_iv + jj] *= s;
                    }
                }
            }

            detect_collinearity(ZtZ, K_iv, is_collinear_z, verbose);
            for (ST_int k = 0; k < K_iv; k++) if (stage1_z[k]) is_collinear_z[k] = 1;
        }
        free(ZtZ);
        free(stage1_z);

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
    ST_int *test_idx = (ST_int *)malloc(((size_t)n_orthog_t + n_endogtest_t + n_redund_t + 1) * sizeof(ST_int));
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
        free(is_collinear_x);
        free(test_idx);
        if (!test_rc) SF_error("civreghdfe: Model is underidentified after removing collinear variables\n");
        free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        ctools_hdfe_state_cleanup(state);
        free(state);
        g_state = NULL;
        ctools_aligned_free(est_obs);
        return test_rc ? test_rc : 481;
    }

    /* Store compacted dimensions */
    ctools_scal_save("__civreghdfe_K_keep", (ST_double)K_total);

    /* Allocate output arrays */
    ST_double *beta = (ST_double *)calloc(K_total, sizeof(ST_double));
    ST_double *V = (ST_double *)ctools_safe_calloc3((size_t)K_total, (size_t)K_total, sizeof(ST_double));
    ST_double *first_stage_F = (ST_double *)calloc(K_endog, sizeof(ST_double));

    /* Check allocations - critical for preventing NULL pointer dereference */
    if (!beta || !V || !first_stage_F) {
        SF_error("civreghdfe: Memory allocation failed for output arrays\n");
        free(beta); free(V); free(first_stage_F);
        free(is_collinear_x); free(test_idx);
        free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        ctools_hdfe_state_cleanup(state);
        free(state);
        g_state = NULL;
        ctools_aligned_free(est_obs);
        return 920;  /* Memory allocation error */
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
    ST_int *ts_int = NULL;  /* panel ids (N), then time index (N) */
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
            ST_int tmax_idx = 0;
            for (ST_int i = 0; i < N; i++) {
                tidx[i] = (ST_int)llround((tv[i] - tmin) / tdelta);
                if (tidx[i] > tmax_idx) tmax_idx = tidx[i];
            }
            /* Time span of the estimation sample: the ado reports T = span + 1
               as Kiefer's bandwidth; as ivreg2, a bandwidth may not exceed
               span / delta */
            ctools_scal_save("__civreghdfe_tspan", tmax - tmin);
            if (kernel_type > 0 && !kiefer && (ST_double)bw > (tmax - tmin) / tdelta) {
                SF_error("invalid bandwidth in option bw() - cannot exceed timespan of data\n");
                rc = 198;
            }
            ST_int num_panels = 1;
            if (ts_mode == 2) {
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

    ST_double lambda = 1.0;
    if (!rc) rc = ivest_compute_2sls(
        y_dem, X_exog_dem, X_endog_dem, Z_dem,
        weights_c, weight_type,
        N, (ST_int)N_eff, K_exog, K_endog, K_iv,
        beta, V, first_stage_F,
        vce_type, cluster_ids_c, num_clusters,
        cluster2_ids_c, num_clusters2,
        df_a_for_vce, nested_adj, verbose,
        est_method, kclass_user, fuller_alpha, &lambda,
        kernel_type, bw, kiefer,
        ts_mode > 0 ? &ts_info : NULL,
        (G > 0) ? (df_a_for_vce > 0 ? df_a_for_vce : 1) : sdofminus_opt, center,
        &test_cols
    );
    free(test_idx);
    free(ts_int);

    if (rc != STATA_OK) {
        /* Cleanup and return */
        free(is_collinear_x);
        free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
        free(beta); free(V); free(first_stage_F);
        ctools_hdfe_state_cleanup(state);
        free(state);
        g_state = NULL;
        ctools_aligned_free(est_obs);
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

    ST_double resid_idx_value = 0;
    SF_scal_use("__civreghdfe_resid_idx", &resid_idx_value);
    ST_int resid_idx = (ST_int)resid_idx_value;
    ST_retcode residual_rc = 0;
    /* Compute RSS = sum(w * resid^2) */
    /* Compute fitted values: yhat = X_exog * beta[0:K_exog-1] + X_endog * beta[K_exog:K_total-1] */
    for (ST_int i = 0; i < N; i++) {
        ST_double fitted = 0.0;
        for (ST_int j = 0; j < K_exog; j++) {
            fitted += X_exog_dem[(size_t)j * N + i] * beta[j];
        }
        for (ST_int j = 0; j < K_endog; j++) {
            fitted += X_endog_dem[(size_t)j * N + i] * beta[K_exog + j];
        }
        ST_double resid = y_dem[i] - fitted;
        if (resid_idx > 0 && !residual_rc)
            residual_rc = (_stata_)->safestore(resid_idx, (ST_int)est_obs[i], resid);
        ST_double w = (weights_c && weight_type != 0) ? weights_c[i] : 1.0;
        rss += w * resid * resid;
    }

    ST_double r2 = (tss > 0) ? 1.0 - rss / tss : 0.0;

    /* Compute df_r:
       - For clustered VCE (vce_type == 2, 3) and Driscoll-Kraay (4, clusters
         = the time periods of the estimation sample): df_r = num_clusters - 1
       - Otherwise: df_r = N - K - df_a - dofminus
       - nopartialsmall: exclude partialled variables from K for DOF calc */
    ST_int K_for_dof = nopartialsmall ? (K_total - n_partial) : K_total;
    ST_int df_r_val;
    if (kiefer) {
        /* Kiefer VCE: df_r = N_eff - K - df_a - dofminus (not cluster-based) */
        df_r_val = (ST_int)(N_eff - K_for_dof - df_a - dofminus);
    } else if (has_cluster2) {
        /* Two-way clustering: df_r = min(G1, G2) - 1 */
        ST_int min_clust = (num_clusters < num_clusters2) ? num_clusters : num_clusters2;
        df_r_val = min_clust - 1;
    } else if (has_cluster) {
        df_r_val = num_clusters - 1;
    } else {
        df_r_val = (ST_int)(N_eff - K_for_dof - df_a - dofminus);
    }
    if (df_r_val <= 0) df_r_val = 1;
    /* For rmse, use sdofminus = max(1, df_a_for_vce) when nested, df_a otherwise
       (matches ivreghdfe line 670: if HDFE.df_a=0, force absorb_ct to 1) */
    ST_int sdofminus = (has_cluster && df_a_nested > 0) ?
                       (df_a_for_vce > 0 ? df_a_for_vce : 1) : df_a;
    ST_double rmse = sqrt(rss / ((ST_double)N_eff - K_total - sdofminus > 0 ? (ST_double)N_eff - K_total - sdofminus : 1));

    /* End stats timing, start store timing */
    double t_stats = ctools_timer_seconds() - t_stats_start;
    double t_store_start = ctools_timer_seconds();

    /* Store results to Stata using checked wrappers */
    /* Scalars - use error-checking wrappers to prevent silent failures */
    ctools_scal_save("__civreghdfe_N", N_eff);
    ctools_scal_save("__civreghdfe_df_r", (ST_double)df_r_val);
    ctools_scal_save("__civreghdfe_df_a", (ST_double)df_a);
    ctools_scal_save("__civreghdfe_df_a_for_vce", (ST_double)df_a_for_vce);
    ctools_scal_save("__civreghdfe_K", (ST_double)K_total);
    ctools_scal_save("__civreghdfe_rss", rss);
    ctools_scal_save("__civreghdfe_tss", tss);
    ctools_scal_save("__civreghdfe_r2", r2);
    ctools_scal_save("__civreghdfe_rmse", rmse);

    if (has_cluster2) {
        /* Two-way clustering: N_clust = min(G1, G2), store both counts separately */
        ST_int min_clust = (num_clusters < num_clusters2) ? num_clusters : num_clusters2;
        ctools_scal_save("__civreghdfe_N_clust", (ST_double)min_clust);
        ctools_scal_save("__civreghdfe_N_clust1", (ST_double)num_clusters);
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

    ST_retcode matrix_rc = store_iv_matrices(beta, V, K_exog_orig,
                                             K_endog_orig, K_total, is_collinear_x);
    if (!residual_rc) residual_rc = matrix_rc;
    if (residual_rc) goto result_cleanup;

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

result_cleanup:
    /* Cleanup, including failed result publication. */
    free(is_collinear_x);
    free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
    free(weights_c); free(cluster_ids_c); free(cluster2_ids_c);
    free(beta); free(V); free(first_stage_F);

    /* Cleanup state */
    ctools_hdfe_state_cleanup(state);
    free(fe_levels_c);
    free(state);
    g_state = NULL;

    ctools_aligned_free(est_obs);

    return residual_rc;
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
    /* civreghdfe cleans up g_state at the end of do_iv_regression(),
       so this is just a safety measure for interrupted execution */
    if (g_state != NULL) {
        /* Note: Full cleanup would require knowing G and num_threads,
           which we don't have here. The main function handles full cleanup.
           This just nulls the pointer to prevent double-free issues. */
        g_state = NULL;
    }
}

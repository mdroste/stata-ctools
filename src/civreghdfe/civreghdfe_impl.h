/*
    civreghdfe_impl.h
    Instrumental Variables Regression with High-Dimensional Fixed Effects

    This module implements 2SLS/IV regression with HDFE absorption,
    reusing the creghdfe infrastructure for HDFE demeaning.
*/

#ifndef CIVREGHDFE_IMPL_H
#define CIVREGHDFE_IMPL_H

#include "../stplugin.h"
#include "../ctools_types.h"
#include "../creghdfe/creghdfe_types.h"
#include "../creghdfe/creghdfe_hdfe.h"
#include "../creghdfe/creghdfe_solver.h"
#include "../creghdfe/creghdfe_ols.h"
#include "../creghdfe/creghdfe_vce.h"
#include "../creghdfe/creghdfe_utils.h"

/*
    Main entry point for civreghdfe plugin calls.

    Subcommands:
    - "iv_regression" - Full IV/2SLS regression with HDFE

    Returns STATA_OK on success, error code on failure.
*/
ST_retcode civreghdfe_main(const char *args);

/*
    Perform k-class IV estimation on demeaned data (includes 2SLS, LIML, Fuller).

    After HDFE partialling out, this computes:
    1. First stage: X_endog = Z * (Z'Z)^-1 * Z'X_endog + residual
    2. For LIML/Fuller: compute lambda (min eigenvalue for k-class)
    3. Compute k-class estimator: beta = ((1-k)X'X + kX'P_Z X)^-1 ((1-k)X'y + kX'P_Z y)
    4. Corrected VCE using original X, not X_hat

    Parameters:
    - y: demeaned dependent variable (N x 1)
    - X_exog: demeaned exogenous regressors (N x K_exog), may be NULL if K_exog=0
    - X_endog: demeaned endogenous regressors (N x K_endog)
    - Z: demeaned instruments (N x K_iv), must include exogenous vars
    - weights: optional weights array
    - weight_type: 0=none, 1=aweight, 2=fweight, 3=pweight
    - N: number of observations
    - K_exog: number of exogenous regressors (excluding constant)
    - K_endog: number of endogenous regressors
    - K_iv: total number of instruments (must be >= K_endog + K_exog)
    - beta: output coefficient vector (K_exog + K_endog)
    - V: output variance-covariance matrix ((K+1) x (K+1) for constant)
    - first_stage_F: output array of first-stage F-statistics (K_endog)
    - vce_type: 0=unadjusted, 1=robust, 2=cluster
    - cluster_ids: cluster IDs (may be NULL)
    - num_clusters: number of clusters
    - df_a: absorbed degrees of freedom
    - nested_adj: 1 if FE nested in cluster
    - verbose: 1 for verbose output
    - est_method: 0=2SLS, 1=LIML, 2=Fuller, 3=kclass, 4=GMM2S, 5=CUE
    - kclass_user: user-specified k for kclass method
    - fuller_alpha: Fuller modification parameter
    - lambda_out: output for LIML lambda value (may be NULL)

    Returns STATA_OK on success.
*/
ST_retcode compute_2sls(
    const ST_double *y,
    const ST_double *X_exog,
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_double *beta,
    ST_double *V,
    ST_double *first_stage_F,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_int df_a,
    ST_int nested_adj,
    ST_int verbose,
    ST_int est_method,
    ST_double kclass_user,
    ST_double fuller_alpha,
    ST_double *lambda_out
);

#endif /* CIVREGHDFE_IMPL_H */

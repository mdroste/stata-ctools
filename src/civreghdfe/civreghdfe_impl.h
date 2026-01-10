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
    Perform 2SLS estimation on demeaned data.

    After HDFE partialling out, this computes:
    1. First stage: X_endog = Z * (Z'Z)^-1 * Z'X_endog + residual
    2. Fitted values: X_hat = X_exog + X_endog_fitted
    3. Second stage: y = X_hat * beta + error
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
    ST_int verbose
);

#endif /* CIVREGHDFE_IMPL_H */

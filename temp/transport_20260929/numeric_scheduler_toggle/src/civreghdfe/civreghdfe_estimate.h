/*
    civreghdfe_estimate.h
    Core IV Estimation: 2SLS, LIML, Fuller, GMM2S, CUE

    This module implements the k-class family of IV estimators.
    VCE computation is handled separately in civreghdfe_vce.c
*/

#ifndef CIVREGHDFE_ESTIMATE_H
#define CIVREGHDFE_ESTIMATE_H

#include "../stplugin.h"
#include "../ctools_types.h"
#include "civreghdfe_omega.h"

/*
    Estimation method constants
*/
#define CIVREGHDFE_EST_2SLS   0
#define CIVREGHDFE_EST_LIML   1
#define CIVREGHDFE_EST_FULLER 2
#define CIVREGHDFE_EST_KCLASS 3
#define CIVREGHDFE_EST_GMM2S  4
#define CIVREGHDFE_EST_CUE    5

/*
    IV Estimation context - holds pre-computed matrices for estimation
*/
typedef struct {
    /* Input dimensions */
    ST_int N;           /* Number of observations */
    ST_int K_exog;      /* Number of exogenous regressors */
    ST_int K_endog;     /* Number of endogenous regressors */
    ST_int K_iv;        /* Number of instruments (must >= K_exog + K_endog) */
    ST_int K_total;     /* K_exog + K_endog */

    /* Weights */
    const ST_double *weights;
    ST_int weight_type; /* 0=none, 1=aweight, 2=fweight, 3=pweight */

    /* Original data (not owned, do not free) */
    const ST_double *Z;     /* Instruments (N x K_iv) - needed for GMM2S/CUE */
    const ST_double *y;     /* Dependent variable (N x 1) */

    /* Pre-computed matrices (owned by context, freed in ivest_free_context) */
    ST_double *ZtZ;         /* Z'Z (K_iv x K_iv) */
    ST_double *ZtZ_inv;     /* (Z'Z)^-1 (K_iv x K_iv) */
    ST_double *ZtX;         /* Z'X (K_iv x K_total) */
    ST_double *Zty;         /* Z'y (K_iv x 1) */
    ST_double *XtPzX;       /* X'P_Z X (K_total x K_total) */
    ST_double *XtPzy;       /* X'P_Z y (K_total x 1) */
    ST_double *temp_kiv_ktotal; /* Temp storage: (Z'Z)^-1 Z'X (K_iv x K_total) */

    /* Combined data arrays (owned by context) */
    ST_double *X_all;       /* [X_exog, X_endog] (N x K_total) */

    /* Estimation method parameters */
    ST_int est_method;      /* CIVREGHDFE_EST_* constant */
    ST_double kclass_user;  /* User-specified k for kclass */
    ST_double fuller_alpha; /* Fuller modification parameter */

    /* Moment covariance of the chosen VCE (CUE weighting matrix) */
    const civ_omega *om;

    /* Verbose output */
    ST_int verbose;
} IVEstContext;

/*
    Compute GMM2S (two-step efficient GMM) estimator.

    Uses the optimal weighting matrix W = S^-1, where S is the moment
    covariance of the chosen VCE built from the first-step (2SLS) residuals
    (civ_omega_build: robust, cluster, two-way, HAC, AC, Driscoll-Kraay or
    Kiefer), as ivreg2's s_egmm.

    Parameters:
    - ctx: Initialized estimation context
    - y: Dependent variable
    - S: First-step moment covariance (K_iv x K_iv, raw)
    - beta: Output coefficients
    - resid: Output residuals
    - XZWZX_inv_out: Output GMM Hessian inverse (K_total x K_total), may be NULL

    Returns STATA_OK on success.
*/
ST_retcode ivest_compute_gmm2s(
    IVEstContext *ctx,
    const ST_double *y,
    const ST_double *S,
    ST_double *beta,
    ST_double *resid,
    ST_double *XZWZX_inv_out
);

/*
    Compute CUE (Continuously Updated Estimator).

    Minimizes J(b) = g(b)' S(b)^-1 g(b), where S(b) is the moment covariance
    of the chosen VCE (ctx->om) evaluated at the residuals of b.

    Parameters:
    - ctx: Initialized estimation context (ctx->om set)
    - y: Dependent variable
    - initial_beta: Starting point from 2SLS/GMM2S
    - beta: Output coefficients
    - resid: Output residuals
    - max_iter: Maximum iterations
    - tol: Convergence tolerance
    - XZWZX_inv_out: Output final CUE Hessian inverse (K_total x K_total), may be NULL

    Returns STATA_OK on success.
*/
ST_retcode ivest_compute_cue(
    IVEstContext *ctx,
    const ST_double *y,
    const ST_double *initial_beta,
    ST_double *beta,
    ST_double *resid,
    ST_int max_iter,
    ST_double tol,
    ST_double *XZWZX_inv_out
);

/*
    Columns selected by orthog(), endogtest() and redundant(): distinct 1-based
    positions among the final (post-collinearity) excluded instruments
    (orthog, redundant) or endogenous regressors (endogtest). Built in
    civreghdfe_impl.c from the ado's per-column flags.
*/
typedef struct {
    const ST_int *orthog;     ST_int n_orthog;
    const ST_int *endogtest;  ST_int n_endogtest;
    const ST_int *redundant;  ST_int n_redundant;
} civreghdfe_test_cols;

/*
    Full k-class IV estimation with VCE and diagnostics.

    This is the main 2SLS/LIML/Fuller/GMM2S/CUE estimation function.
    Computes coefficients, VCE, first-stage F, and diagnostic tests.

    Parameters:
    - y: Dependent variable (N x 1)
    - X_exog: Exogenous regressors (N x K_exog), may be NULL
    - X_endog: Endogenous regressors (N x K_endog)
    - Z: Instruments (N x K_iv)
    - weights, weight_type: Optional weighting
    - N, K_exog, K_endog, K_iv: Dimensions
    - beta: Output coefficients (K_total x 1)
    - V: Output VCE matrix (K_total x K_total)
    - first_stage_F: Output first-stage F-stats (K_endog x 1), may be NULL
    - vce_type: 0=unadjusted, 1=robust, 2=cluster, 3=hac, 4=cluster2
    - cluster_ids, num_clusters: Cluster structure for clustered VCE
    - cluster2_ids, num_clusters2: Second cluster for two-way clustering
    - df_a: Absorbed degrees of freedom
    - nested_adj: Adjustment for nested clusters
    - verbose: Print debug output
    - est_method: 0=2SLS, 1=LIML, 2=Fuller, 3=kclass, 4=GMM2S, 5=CUE
    - kclass_user: User-specified k for kclass estimator
    - fuller_alpha: Fuller modification parameter
    - lambda_out: Output LIML lambda (may be NULL)
    - kernel_type, bw: kernel parameters (HAC with vce_type 1, AC with 0,
      Driscoll-Kraay with 4)
    - kiefer: Use Kiefer (1980) homoskedastic within-panel VCE
    - ts: tsset panel/time structure (required for kernels, Driscoll-Kraay
      and Kiefer; NULL otherwise)
    - test_cols: Columns for the orthog/endogtest/redundant tests (NULL = none)
    - n_absorbed_x, n_absorbed_iv: regressors and instruments absorbed by the
      fixed effects (omitted) that still count in the degrees of freedom, as
      in ivreghdfe's rank counts
    - n_absorbed_iv_all: all absorbed exogenous regressors (included
      instruments), which count in the Cragg-Donald / KP Wald F dof as in
      ivreghdfe's variable-list instrument count

    Returns STATA_OK on success.
*/
ST_retcode ivest_compute_2sls(
    const ST_double *y,
    const ST_double *X_exog,
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_double *beta,
    ST_double *V,
    ST_double *first_stage_F,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    const ST_int *cluster2_ids,
    ST_int num_clusters2,
    ST_int df_a,
    ST_int nested_adj,
    ST_int verbose,
    ST_int est_method,
    ST_double kclass_user,
    ST_double fuller_alpha,
    ST_double *lambda_out,
    ST_int kernel_type,
    ST_int bw,
    ST_int kiefer,
    const civreghdfe_ts *ts,
    ST_int sdofminus,
    ST_int center,
    const civreghdfe_test_cols *test_cols,
    ST_int n_absorbed_x,
    ST_int n_absorbed_iv,
    ST_int n_absorbed_iv_all
);

#endif /* CIVREGHDFE_ESTIMATE_H */

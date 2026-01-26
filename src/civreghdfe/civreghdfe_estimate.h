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

    /* Verbose output */
    ST_int verbose;
} IVEstContext;

/*
    Initialize IV estimation context and compute basic matrices.

    Computes:
    - Z'Z and (Z'Z)^-1
    - Z'X and Z'y
    - X'P_Z X and X'P_Z y

    Parameters:
    - ctx: Context to initialize (caller provides storage)
    - y: Dependent variable (N x 1)
    - X_exog: Exogenous regressors (N x K_exog), may be NULL
    - X_endog: Endogenous regressors (N x K_endog)
    - Z: Instruments (N x K_iv)
    - weights, weight_type: Optional weighting
    - N, K_exog, K_endog, K_iv: Dimensions
    - verbose: Print debug info

    Returns STATA_OK on success.
*/
ST_retcode ivest_init_context(
    IVEstContext *ctx,
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
    ST_int verbose
);

/*
    Free resources allocated by ivest_init_context.
*/
void ivest_free_context(IVEstContext *ctx);

/*
    Compute k-class estimator coefficients.

    The k-class estimator is:
    beta = ((1-k)*X'X + k*X'P_Z X)^-1 ((1-k)*X'y + k*X'P_Z y)

    When k=1: 2SLS
    When k=lambda (min eigenvalue): LIML
    When k=lambda - alpha/(N-K_iv): Fuller

    Parameters:
    - ctx: Initialized estimation context
    - y: Dependent variable (N x 1)
    - kclass: The k parameter
    - beta: Output coefficients (K_total x 1)
    - resid: Output residuals (N x 1), may be NULL

    Returns STATA_OK on success.
*/
ST_retcode ivest_compute_kclass(
    IVEstContext *ctx,
    const ST_double *y,
    ST_double kclass,
    ST_double *beta,
    ST_double *resid
);

/*
    Compute GMM2S (two-step efficient GMM) estimator.

    Uses optimal weighting matrix W = (Z'ΩZ)^-1 where Ω = diag(e²)
    or cluster-robust covariance.

    Parameters:
    - ctx: Initialized estimation context
    - y: Dependent variable
    - initial_resid: Initial residuals from 2SLS (N x 1)
    - cluster_ids: Cluster IDs (NULL for heteroskedastic)
    - num_clusters: Number of clusters
    - beta: Output coefficients
    - resid: Output residuals
    - XZWZX_inv_out: Output GMM Hessian inverse (K_total x K_total), may be NULL

    Returns STATA_OK on success.
*/
ST_retcode ivest_compute_gmm2s(
    IVEstContext *ctx,
    const ST_double *y,
    const ST_double *initial_resid,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_double *beta,
    ST_double *resid,
    ST_double *XZWZX_inv_out
);

/*
    Compute CUE (Continuously Updated Estimator).

    Iteratively re-weights until convergence.

    Parameters:
    - ctx: Initialized estimation context
    - y: Dependent variable
    - initial_beta: Starting point from 2SLS/GMM2S
    - vce_type: 0=homoskedastic, 1=robust, 2=cluster, 3=twoway cluster
    - cluster_ids: Cluster IDs (NULL for non-cluster)
    - num_clusters: Number of clusters
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
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_double *beta,
    ST_double *resid,
    ST_int max_iter,
    ST_double tol,
    ST_double *XZWZX_inv_out
);

/*
    Compute first-stage F-statistics for each endogenous variable.

    Parameters:
    - ctx: Initialized estimation context
    - X_endog: Endogenous regressors (N x K_endog)
    - df_a: Absorbed degrees of freedom
    - first_stage_F: Output F-statistics (K_endog x 1)

    Returns STATA_OK on success.
*/
ST_retcode ivest_first_stage_f(
    IVEstContext *ctx,
    const ST_double *X_endog,
    ST_int df_a,
    ST_double *first_stage_F
);

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
    - kernel_type, bw: HAC kernel parameters
    - kiefer: Use Kiefer (1980) homoskedastic within-panel VCE
    - hac_panel_ids, num_hac_panels: Panel IDs for panel-aware HAC (NULL if not used)

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
    const ST_int *hac_panel_ids,
    ST_int num_hac_panels
);

#endif /* CIVREGHDFE_ESTIMATE_H */

/*
    civreghdfe_vce.h
    Variance-Covariance Estimation for IV Regression

    Implements unadjusted, robust (HC), clustered, and HAC VCE
    for k-class IV estimators.
*/

#ifndef CIVREGHDFE_VCE_H
#define CIVREGHDFE_VCE_H

#include "../stplugin.h"
#include "../ctools_types.h"
#include "civreghdfe_estimate.h"

/*
    VCE type constants
*/
#define CIVREGHDFE_VCE_UNADJUSTED 0
#define CIVREGHDFE_VCE_ROBUST     1
#define CIVREGHDFE_VCE_CLUSTER    2

/*
    HAC kernel types
*/
#define CIVREGHDFE_KERNEL_NONE      0
#define CIVREGHDFE_KERNEL_BARTLETT  1
#define CIVREGHDFE_KERNEL_PARZEN    2
#define CIVREGHDFE_KERNEL_QS        3
#define CIVREGHDFE_KERNEL_TRUNCATED 4
#define CIVREGHDFE_KERNEL_TUKEY     5

/*
    Compute unadjusted VCE for k-class estimator.

    V = sigma^2 * (XkX)^-1
    where sigma^2 = RSS / df_r

    Parameters:
    - XkX_inv: (XkX)^-1, the k-class matrix inverse (K_total x K_total)
    - rss: Residual sum of squares
    - df_r: Residual degrees of freedom (N - K_total - df_a)
    - K_total: Number of regressors
    - V: Output VCE matrix (K_total x K_total)
*/
void ivvce_compute_unadjusted(
    const ST_double *XkX_inv,
    ST_double rss,
    ST_int df_r,
    ST_int K_total,
    ST_double *V
);

/*
    Compute robust (HC1) VCE for k-class estimator.

    V = (XkX)^-1 * Meat * (XkX)^-1 * dof_adj
    where Meat = sum_i (PzX_i * e_i)^2

    Parameters:
    - ctx: IV estimation context
    - resid: Residuals (N x 1)
    - XkX_inv: Inverse of k-class matrix (K_total x K_total)
    - df_a: Absorbed degrees of freedom
    - V: Output VCE matrix (K_total x K_total)
*/
void ivvce_compute_robust(
    IVEstContext *ctx,
    const ST_double *resid,
    const ST_double *XkX_inv,
    ST_int df_a,
    ST_double *V
);

/*
    Compute clustered VCE for k-class estimator.

    V = (XkX)^-1 * Meat * (XkX)^-1 * dof_adj
    where Meat = sum_c (sum_i in c (PzX_i * e_i))^2

    Parameters:
    - ctx: IV estimation context
    - resid: Residuals (N x 1)
    - XkX_inv: Inverse of k-class matrix
    - cluster_ids: Cluster IDs (1-indexed)
    - num_clusters: Number of clusters
    - df_a: Absorbed degrees of freedom
    - nested_adj: 1 if FE nested in cluster, else 0
    - V: Output VCE matrix
*/
void ivvce_compute_cluster(
    IVEstContext *ctx,
    const ST_double *resid,
    const ST_double *XkX_inv,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_int df_a,
    ST_int nested_adj,
    ST_double *V
);

/*
    Compute HAC (Heteroskedasticity and Autocorrelation Consistent) VCE.

    Uses kernel weighting for cross-lag products.

    Parameters:
    - ctx: IV estimation context
    - resid: Residuals (N x 1)
    - XkX_inv: Inverse of k-class matrix
    - kernel_type: CIVREGHDFE_KERNEL_* constant
    - bw: Bandwidth
    - df_a: Absorbed degrees of freedom
    - V: Output VCE matrix
*/
void ivvce_compute_hac(
    IVEstContext *ctx,
    const ST_double *resid,
    const ST_double *XkX_inv,
    ST_int kernel_type,
    ST_int bw,
    ST_int df_a,
    ST_double *V
);

/*
    Compute efficient GMM2S VCE.

    When using optimal weights W = (Z'ΩZ)^-1, the VCE is:
    V = (X'ZWZ'X)^-1

    Parameters:
    - XZWZX_inv: Inverse of X'ZWZ'X (K_total x K_total)
    - K_total: Number of regressors
    - V: Output VCE matrix
*/
void ivvce_compute_gmm2s(
    const ST_double *XZWZX_inv,
    ST_int K_total,
    ST_double *V
);

/*
    Compute Z'ΩZ matrix for heteroskedastic case.

    Ω = diag(e²)
    Returns Z' * diag(e²) * Z

    Parameters:
    - Z: Instruments (N x K_iv)
    - resid: Residuals (N x 1)
    - weights: Optional weights
    - weight_type: Weight type
    - N: Number of observations
    - K_iv: Number of instruments
    - ZOmegaZ: Output matrix (K_iv x K_iv)
*/
void ivvce_compute_ZOmegaZ_robust(
    const ST_double *Z,
    const ST_double *resid,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int K_iv,
    ST_double *ZOmegaZ
);

/*
    Compute Z'ΩZ matrix for cluster case.

    Returns sum_c (Z_c'e_c)(e_c'Z_c)

    Parameters:
    - Z: Instruments (N x K_iv)
    - resid: Residuals (N x 1)
    - weights: Optional weights
    - weight_type: Weight type
    - cluster_ids: Cluster IDs (1-indexed)
    - N: Number of observations
    - K_iv: Number of instruments
    - num_clusters: Number of clusters
    - ZOmegaZ: Output matrix (K_iv x K_iv)
*/
void ivvce_compute_ZOmegaZ_cluster(
    const ST_double *Z,
    const ST_double *resid,
    const ST_double *weights,
    ST_int weight_type,
    const ST_int *cluster_ids,
    ST_int N,
    ST_int K_iv,
    ST_int num_clusters,
    ST_double *ZOmegaZ
);

/*
    Full VCE computation with P_Z X calculation.

    This is the main entry point for VCE computation when Z is available.
    Handles unadjusted, robust (HC), clustered, and HAC VCE types.

    Parameters:
    - Z: Instruments (N x K_iv)
    - resid: Residuals (N x 1)
    - temp_kiv_ktotal: (Z'Z)^-1 Z'X (K_iv x K_total)
    - XkX_inv: Inverse of k-class matrix (K_total x K_total)
    - weights, weight_type: Weighting
    - N, K_total, K_iv: Dimensions
    - vce_type: CIVREGHDFE_VCE_* constant
    - cluster_ids, num_clusters: Clustering info
    - df_a, nested_adj: DOF adjustments
    - kernel_type, bw: HAC parameters
    - V: Output VCE matrix (K_total x K_total)
*/
void ivvce_compute_full(
    const ST_double *Z,
    const ST_double *resid,
    const ST_double *temp_kiv_ktotal,
    const ST_double *XkX_inv,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int K_total,
    ST_int K_iv,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_int df_a,
    ST_int nested_adj,
    ST_int kernel_type,
    ST_int bw,
    ST_double *V
);

#endif /* CIVREGHDFE_VCE_H */

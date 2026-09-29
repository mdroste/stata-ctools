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
#define CIVREGHDFE_VCE_CLUSTER2   3  /* Two-way clustering */
#define CIVREGHDFE_VCE_DKRAAY     4  /* Driscoll-Kraay (clustered by time with HAC) */

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
    Compute two-way clustered VCE using Cameron-Gelbach-Miller (2011) formula.

    V_twoway = V_cluster1 + V_cluster2 - V_intersection
    where V_intersection is the VCE clustering on the intersection of the two variables.

    Parameters:
    - Z: Instruments (N x K_iv)
    - resid: Residuals (N x 1)
    - temp_kiv_ktotal: (Z'Z)^-1 Z'X (K_iv x K_total)
    - XkX_inv: Inverse of k-class matrix (K_total x K_total)
    - weights, weight_type: Weighting
    - N, K_total, K_iv: Dimensions
    - cluster1_ids: First cluster IDs (1-indexed)
    - num_clusters1: Number of first-dimension clusters
    - cluster2_ids: Second cluster IDs (1-indexed)
    - num_clusters2: Number of second-dimension clusters
    - df_a: Absorbed degrees of freedom
    - V: Output VCE matrix (K_total x K_total)
*/
void ivvce_compute_twoway(
    const ST_double *Z,
    const ST_double *resid,
    const ST_double *temp_kiv_ktotal,
    const ST_double *XkX_inv,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
    ST_int K_total,
    ST_int K_iv,
    const ST_int *cluster1_ids,
    ST_int num_clusters1,
    const ST_int *cluster2_ids,
    ST_int num_clusters2,
    ST_int df_a,
    ST_double *V
);

/*
    Sandwich VCE from a Z-space moment covariance (the kernel VCEs: HAC,
    AC, Driscoll-Kraay, Kiefer).

    V = XkX_inv * A'SA * XkX_inv * dof_adj, A = (Z'Z)^-1 Z'X

    Parameters:
    - temp_kiv_ktotal: A = (Z'Z)^-1 Z'X (K_iv x K_total)
    - S: Moment covariance (K_iv x K_iv, raw, from civ_omega_build)
    - XkX_inv: Inverse of k-class matrix (K_total x K_total)
    - K_total, K_iv: Dimensions
    - dof_adj: Small-sample factor
    - V: Output VCE matrix (K_total x K_total)

    Returns 0 or 920.
*/
ST_retcode ivvce_compute_sandwich(
    const ST_double *temp_kiv_ktotal,
    const ST_double *S,
    const ST_double *XkX_inv,
    ST_int K_total,
    ST_int K_iv,
    ST_double dof_adj,
    ST_double *V
);

/*
    Full VCE computation with P_Z X calculation.

    This is the main entry point for VCE computation when Z is available.
    Handles unadjusted, robust (HC) and one-way clustered VCE types.

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
    - V: Output VCE matrix (K_total x K_total)
*/
ST_retcode ivvce_compute_full(
    const ST_double *Z,
    const ST_double *resid,
    const ST_double *temp_kiv_ktotal,
    const ST_double *XkX_inv,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
    ST_int K_total,
    ST_int K_iv,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_int df_a,
    ST_int nested_adj,
    ST_double *V
);

#endif /* CIVREGHDFE_VCE_H */

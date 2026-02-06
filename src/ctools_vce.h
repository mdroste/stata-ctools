/*
 * ctools_vce.h
 *
 * Shared sandwich VCE computation for creghdfe and civreghdfe.
 * Computes V = D * meat * D * dof_adj for both OLS and IV estimators.
 * Part of the ctools Stata plugin suite.
 */

#ifndef CTOOLS_VCE_H
#define CTOOLS_VCE_H

#include "stplugin.h"

/*
 * Data needed for sandwich VCE computation.
 *
 * The sandwich formula is V = D * meat * D * dof_adj, where:
 * - For OLS:  X_eff = X,     D = (X'X)^-1
 * - For IV:   X_eff = P_Z*X, D = (XkX)^-1
 */
typedef struct {
    const ST_double *X_eff;   /* Effective regressors: N x K, column-major */
    const ST_double *D;       /* Bread: (X'X)^-1 or (XkX)^-1, K x K */
    const ST_double *resid;   /* Residuals: N x 1 */
    const ST_double *weights; /* Weights: N x 1, NULL if unweighted */
    ST_int weight_type;       /* 0=none, 1=aweight, 2=fweight, 3=pweight */
    ST_int N;                 /* Number of observations */
    ST_int K;                 /* Number of effective regressors */
} ctools_vce_data;

/*
 * Compute robust (HC1) sandwich VCE.
 *
 * V = D * meat * D * dof_adj
 * where meat = sum_i w_i * (X_eff_i)(X_eff_i)'
 * with w_i = e_i^2 (unweighted), e_i^2 * fw_i (fweight), or (e_i * w_norm_i)^2 (aw/pw)
 *
 * Parameters:
 * - d: VCE data (X_eff, D, residuals, weights)
 * - dof_adj: Degrees of freedom adjustment scalar
 * - V: Output K x K VCE matrix
 */
void ctools_vce_robust(const ctools_vce_data *d, ST_double dof_adj, ST_double *V);

/*
 * Compute clustered sandwich VCE.
 *
 * V = D * meat * D * dof_adj
 * where meat = sum_c (sum_{i in c} X_eff_i * e_w_i)(sum_{i in c} X_eff_i * e_w_i)'
 *
 * Parameters:
 * - d: VCE data (X_eff, D, residuals, weights)
 * - cluster_ids: Cluster IDs (0-indexed), N x 1
 * - num_clusters: Number of unique clusters
 * - dof_adj: Degrees of freedom adjustment scalar
 * - V: Output K x K VCE matrix
 */
void ctools_vce_cluster(const ctools_vce_data *d,
                        const ST_int *cluster_ids, ST_int num_clusters,
                        ST_double dof_adj, ST_double *V);

#endif /* CTOOLS_VCE_H */

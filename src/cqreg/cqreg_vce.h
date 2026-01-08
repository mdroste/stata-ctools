/*
 * cqreg_vce.h
 *
 * Variance-covariance estimation for quantile regression.
 * Implements IID, robust (sandwich), and clustered VCE.
 * Part of the ctools suite.
 */

#ifndef CQREG_VCE_H
#define CQREG_VCE_H

#include "cqreg_types.h"

/* ============================================================================
 * Main Interface
 * ============================================================================ */

/*
 * Compute variance-covariance matrix based on VCE type.
 * Dispatches to the appropriate estimator.
 *
 * Parameters:
 *   state - cqreg state with regression results
 *   X     - Design matrix (N x K, column-major)
 *
 * Returns:
 *   0 on success, error code otherwise
 */
ST_int cqreg_compute_vce(cqreg_state *state, const ST_double *X);

/* ============================================================================
 * IID VCE
 * ============================================================================ */

/*
 * Compute IID variance: V = (1/n) * sparsity^2 * q*(1-q) * (X'X)^{-1}
 *
 * This assumes errors are i.i.d. and uses kernel density estimation
 * for the sparsity function.
 *
 * Parameters:
 *   V         - Output: variance-covariance matrix (K x K)
 *   X         - Design matrix (N x K, column-major)
 *   residuals - Regression residuals (N)
 *   N, K      - Dimensions
 *   q         - Quantile
 *   sparsity  - Estimated sparsity (1/f(0))
 *
 * Returns:
 *   0 on success, -1 on failure
 */
ST_int cqreg_vce_iid(ST_double *V,
                     const ST_double *X,
                     const ST_double *residuals,
                     ST_int N, ST_int K,
                     ST_double q,
                     ST_double sparsity);

/* ============================================================================
 * Robust (Sandwich) VCE
 * ============================================================================ */

/*
 * Compute robust variance using Powell sandwich estimator.
 * V = J^{-1} * I * J^{-1} where:
 *   J = (1/n) * sum_i f_i(0) * X_i X_i' (density-weighted Hessian)
 *   I = tau*(1-tau) * (1/n) * X'X       (score variance)
 *
 * This allows for heteroskedasticity in the error distribution.
 *
 * Parameters:
 *   V         - Output: variance-covariance matrix (K x K)
 *   X         - Design matrix (N x K, column-major)
 *   residuals - Regression residuals (N)
 *   N, K      - Dimensions
 *   q         - Quantile
 *   bandwidth - Kernel bandwidth (probability scale)
 *   sparsity  - Estimated sparsity (for bandwidth conversion to residual scale)
 *
 * Returns:
 *   0 on success, -1 on failure
 */
ST_int cqreg_vce_robust(ST_double *V,
                        const ST_double *X,
                        const ST_double *residuals,
                        ST_int N, ST_int K,
                        ST_double q,
                        ST_double bandwidth,
                        ST_double sparsity);

/* ============================================================================
 * Cluster-Robust VCE
 * ============================================================================ */

/*
 * Compute cluster-robust variance.
 * V = D * M_cluster * D with small-sample adjustment.
 *
 * Parameters:
 *   V           - Output: variance-covariance matrix (K x K)
 *   X           - Design matrix (N x K, column-major)
 *   residuals   - Regression residuals (N)
 *   cluster_ids - Cluster identifiers for each observation (N)
 *   num_clusters- Number of unique clusters
 *   N, K        - Dimensions
 *   q           - Quantile
 *   bandwidth   - Kernel bandwidth
 *
 * Returns:
 *   0 on success, -1 on failure
 */
ST_int cqreg_vce_cluster(ST_double *V,
                         const ST_double *X,
                         const ST_double *residuals,
                         const ST_int *cluster_ids,
                         ST_int num_clusters,
                         ST_int N, ST_int K,
                         ST_double q,
                         ST_double bandwidth);

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

/*
 * Compute (X'X)^{-1} using Cholesky decomposition.
 *
 * Parameters:
 *   XtX_inv - Output: inverse matrix (K x K)
 *   X       - Design matrix (N x K, column-major)
 *   N, K    - Dimensions
 *
 * Returns:
 *   0 on success, -1 if matrix is singular
 */
ST_int cqreg_compute_xtx_inv(ST_double *XtX_inv,
                             const ST_double *X,
                             ST_int N, ST_int K);

/*
 * Compute sandwich product: V = A * B * A'
 * All matrices are K x K, A is typically (X'X)^{-1}.
 *
 * Parameters:
 *   V - Output: sandwich product (K x K)
 *   A - Left matrix (K x K)
 *   B - Middle matrix (K x K)
 *   K - Dimension
 */
void cqreg_sandwich_product(ST_double *V,
                            const ST_double *A,
                            const ST_double *B,
                            ST_int K);

/*
 * Count unique clusters and create mapping.
 *
 * Parameters:
 *   cluster_ids  - Raw cluster identifiers (N)
 *   N            - Number of observations
 *   cluster_map  - Output: cluster index for each observation (N)
 *   num_clusters - Output: number of unique clusters
 *
 * Returns:
 *   0 on success, -1 on failure
 */
ST_int cqreg_map_clusters(const ST_int *cluster_ids,
                          ST_int N,
                          ST_int *cluster_map,
                          ST_int *num_clusters);

#endif /* CQREG_VCE_H */

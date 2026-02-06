/*
 * ctools_matrix.h
 *
 * Shared matrix operations and sandwich VCE for creghdfe and civreghdfe.
 * Provides optimized matrix multiplication with OpenMP and SIMD,
 * plus robust and clustered variance-covariance estimation.
 * Part of the ctools Stata plugin suite.
 */

#ifndef CTOOLS_MATRIX_H
#define CTOOLS_MATRIX_H

#include "stplugin.h"

/* ============================================================================
 * Matrix Multiplication
 * ============================================================================ */

/*
 * Matrix multiply C = A' * B
 * A is N x K1, B is N x K2, Result C is K1 x K2
 * All matrices stored in column-major order.
 *
 * OPTIMIZED: Uses K-way unrolled dot products and OpenMP parallelization
 */
void ctools_matmul_atb(const ST_double * restrict A, const ST_double * restrict B,
                       ST_int N, ST_int K1, ST_int K2,
                       ST_double * restrict C);

/*
 * Matrix multiply C = A * B
 * A is K1 x K2, B is K2 x K3, Result C is K1 x K3
 * All matrices stored in column-major order.
 *
 * OPTIMIZED: Uses cache-friendly loop order and OpenMP parallelization
 */
void ctools_matmul_ab(const ST_double * restrict A, const ST_double * restrict B,
                      ST_int K1, ST_int K2, ST_int K3,
                      ST_double * restrict C);

/*
 * Weighted matrix multiply C = A' * diag(w) * B
 * A is N x K1, B is N x K2, w is N x 1 (weights, may be NULL)
 * Result C is K1 x K2
 *
 * OPTIMIZED: Uses SIMD-accelerated weighted dot products and OpenMP parallelization
 */
void ctools_matmul_atdb(const ST_double * restrict A, const ST_double * restrict B,
                        const ST_double * restrict w, ST_int N, ST_int K1, ST_int K2,
                        ST_double * restrict C);

/* ============================================================================
 * Sandwich VCE
 *
 * Computes V = D * meat * D * dof_adj for both OLS and IV estimators.
 * ============================================================================ */

#ifndef CTOOLS_VCE_DATA_DEFINED
#define CTOOLS_VCE_DATA_DEFINED

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
    ST_int normalize_weights; /* 1=normalize aw/pw weights (OLS), 0=use raw (IV) */
} ctools_vce_data;

#endif /* CTOOLS_VCE_DATA_DEFINED */

void ctools_vce_robust(const ctools_vce_data *d, ST_double dof_adj, ST_double *V);

void ctools_vce_cluster(const ctools_vce_data *d,
                        const ST_int *cluster_ids, ST_int num_clusters,
                        ST_double dof_adj, ST_double *V);

#endif /* CTOOLS_MATRIX_H */

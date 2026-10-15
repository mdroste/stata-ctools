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
 * Gram matrix C = A' * A
 * A is N x K, Result C is K x K (symmetric), column-major.
 *
 * Each element is the same unrolled dot product ctools_matmul_atb(A, A, ...)
 * computes, so the two give identical results; only the upper triangle is
 * computed (in parallel) and then mirrored.
 */
void ctools_matmul_ata(const ST_double *A, ST_int N, ST_int K, ST_double *C);

/*
 * y += X * b, with X an N x K column-major matrix and b its K coefficients.
 * Every y[i] accumulates X[i,k] * b[k] in column order k = 0..K-1 with
 * ordinary double arithmetic, exactly as a per-row loop over the columns
 * would; rows are processed in cache-sized blocks, in parallel.
 */
void ctools_matvec_add(const ST_double *X, ST_int N, ST_int K, const ST_double *b, ST_double *y);

/*
 * C = X * B, with X N x K, B K x M and C N x M, all column-major.
 * Each C[i,m] is 0 + X[i,0]*B[0,m] + ... + X[i,K-1]*B[K-1,m], summed in that
 * order (the accumulation order of ctools_matvec_add()).
 */
void ctools_matmul_xb(const ST_double *X, ST_int N, ST_int K, const ST_double *B, ST_int M,
                      ST_double *C);

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

/* Return 0 on success, 920 on allocation failure, 498 on invalid inputs.
 * Callers must not publish a covariance when the return is nonzero. */
ST_retcode ctools_vce_robust(const ctools_vce_data *d, ST_double dof_adj, ST_double *V);

ST_retcode ctools_vce_cluster(const ctools_vce_data *d,
                        const ST_int *cluster_ids, ST_int num_clusters,
                        ST_double dof_adj, ST_double *V);

#endif /* CTOOLS_MATRIX_H */

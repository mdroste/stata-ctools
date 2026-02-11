/*
 * cqreg_linalg.h
 *
 * Linear algebra operations for quantile regression IPM solver.
 * Includes optimized dot products, Cholesky decomposition, and matrix operations.
 * Part of the ctools suite.
 */

#ifndef CQREG_LINALG_H
#define CQREG_LINALG_H

#include "cqreg_types.h"

/* ============================================================================
 * Dot Products (with various optimizations)
 * ============================================================================ */

/*
 * Fast dot product with 8-way loop unrolling.
 * Returns: sum_i(x[i] * y[i])
 */
ST_double cqreg_dot(const ST_double * CQREG_RESTRICT x,
                    const ST_double * CQREG_RESTRICT y,
                    ST_int N);


/*
 * Dot product of a vector with itself (squared norm).
 * Returns: sum_i(x[i]^2)
 */
ST_double cqreg_dot_self(const ST_double * CQREG_RESTRICT x, ST_int N);

/* ============================================================================
 * Vector Operations
 * ============================================================================ */

/*
 * Vector copy: dst = src
 */
void cqreg_vcopy(ST_double * CQREG_RESTRICT dst,
                 const ST_double * CQREG_RESTRICT src,
                 ST_int N);


/* ============================================================================
 * Matrix Operations
 * ============================================================================ */


/*
 * Column-major matrix-vector multiply: y = A * x
 * A is M x N (column-major: A[j*M + i] = A[i,j]), x is N, y is M
 */
void cqreg_matvec_col(ST_double * CQREG_RESTRICT y,
                      const ST_double * CQREG_RESTRICT A,
                      const ST_double * CQREG_RESTRICT x,
                      ST_int M, ST_int N);


/*
 * Compute X' * D * X where D is diagonal.
 * X is N x K (column-major), D is N, result is K x K (symmetric, stored full).
 * Parallel implementation with OpenMP.
 */
void cqreg_xtdx(ST_double * CQREG_RESTRICT XDX,
                const ST_double * CQREG_RESTRICT X,
                const ST_double * CQREG_RESTRICT D,
                ST_int N, ST_int K);

/*
 * Compute X' * v for column-major X.
 * X is N x K (column-major), v is N, result is K.
 */
void cqreg_xtv(ST_double * CQREG_RESTRICT result,
               const ST_double * CQREG_RESTRICT X,
               const ST_double * CQREG_RESTRICT v,
               ST_int N, ST_int K);

/* ============================================================================
 * Cholesky Decomposition and Solve
 * ============================================================================ */

/*
 * Cholesky decomposition: A = L * L'
 * A is K x K symmetric positive definite (stored full, uses lower triangle).
 * L is stored in the lower triangle of A (in-place).
 *
 * Returns: 0 on success, -1 if matrix is not positive definite (collinear)
 */
ST_int cqreg_cholesky(ST_double * CQREG_RESTRICT A, ST_int K);


/*
 * Solve A * x = b where A = L * L' has been Cholesky factored.
 * L is K x K lower triangular, b is K.
 * Solution stored in b (in-place).
 */
void cqreg_solve_cholesky(const ST_double * CQREG_RESTRICT L,
                          ST_double * CQREG_RESTRICT b,
                          ST_int K);

/*
 * Compute inverse of A from its Cholesky factor L.
 * L is K x K lower triangular.
 * Result is stored in Ainv (K x K, full symmetric matrix).
 */
void cqreg_invert_cholesky(ST_double * CQREG_RESTRICT Ainv,
                           const ST_double * CQREG_RESTRICT L,
                           ST_int K);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/*
 * Add regularization to diagonal: A[i,i] += lambda
 * Used when Cholesky fails due to near-singularity.
 */
void cqreg_add_regularization(ST_double * CQREG_RESTRICT A, ST_int K, ST_double lambda);


/* ============================================================================
 * Quantile and Statistics Functions
 * ============================================================================ */

/*
 * Compute sample quantile by sorting and selecting (O(N log N)).
 * Creates a copy of y internally - does not modify input.
 * tau: quantile (0 < tau < 1), e.g., 0.5 for median
 * Returns: the tau-th sample quantile
 */
ST_double cqreg_compute_quantile(const ST_double *y, ST_int N, ST_double tau);

/*
 * Compute raw sum of deviations (check function sum against a fixed value).
 * Returns: tau * sum(y[i] - q | y[i] > q) + (1-tau) * sum(q - y[i] | y[i] <= q)
 */
ST_double cqreg_sum_raw_deviations(const ST_double *y, ST_int N, ST_double q, ST_double tau);

#endif /* CQREG_LINALG_H */

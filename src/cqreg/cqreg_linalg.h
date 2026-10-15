/*
 * cqreg_linalg.h
 *
 * Linear algebra operations for the quantile regression solver: dot products,
 * Cholesky decomposition and selection.
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
 * Cholesky Decomposition and Solve
 * ============================================================================ */

/*
 * Cholesky decomposition: A = L * L'
 * A is K x K symmetric positive definite (stored full, uses lower triangle).
 * L is stored in the lower triangle of A (in-place).
 *
 * Returns: 0 on success, -1 if matrix is not positive definite (collinear)
 */
/* Dense row-major K x K with leading dimension K. In-place ordinary-double
 * lower factorization; upper triangle is retained. Returns 0 or -1 for a
 * nonpositive pivot; on failure the lower factor can be partial. No allocation.
 * The solve below consumes an existing factor, overwriting b in place without
 * scratch or refactorization. Keep this distinct from ctools_solve_cholesky. */
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
 * Select order statistic k in place, with expected O(N) work and a bounded
 * partition budget. Requires finite values, N > 0 and 0 <= k < N.
 */
ST_double cqreg_select(ST_double *values, ST_int N, ST_int k);

/*
 * Compute sample quantile by selection (expected O(N)).
 * Creates a copy of y internally - does not modify input.
 * tau: quantile (0 < tau < 1), e.g., 0.5 for median
 * Returns: the tau-th sample quantile, or NaN when the working copy cannot
 * be allocated (a quantile of finite values is never NaN)
 */
ST_double cqreg_compute_quantile(const ST_double *y, ST_int N, ST_double tau);

/*
 * Compute raw sum of deviations (check function sum against a fixed value).
 * Returns: tau * sum(y[i] - q | y[i] > q) + (1-tau) * sum(q - y[i] | y[i] <= q)
 */
ST_double cqreg_sum_raw_deviations(const ST_double *y, ST_int N, ST_double q, ST_double tau);

#endif /* CQREG_LINALG_H */

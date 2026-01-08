/*
 * cqreg_fn.h
 *
 * Frisch-Newton exact LP solver for quantile regression.
 * Based on Portnoy & Koenker (1997) and Chernozhukov et al. (2020).
 *
 * This solver gives LP-exact solutions (not smoothed approximations),
 * enabling preprocessing for 10x+ speedup.
 */

#ifndef CQREG_FN_H
#define CQREG_FN_H

#include "cqreg_types.h"

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/* Tolerance for identifying active set (observations on hyperplane) */
#define FN_ACTIVE_TOL       1e-6

/* Tolerance for dual feasibility (X'g = 0) - relative to sqrt(N) */
#define FN_DUAL_TOL         5e-2

/* Maximum outer iterations (active set refinement) */
#define FN_MAX_OUTER_ITER   200

/* Maximum inner iterations per active set */
#define FN_MAX_INNER_ITER   20

/* Step-back factor from boundary in ratio test */
#define FN_STEP_BACK        0.9995

/* Regularization for Cholesky */
#define FN_REG              1e-12

/* ============================================================================
 * Main Solver Interface
 * ============================================================================ */

/*
 * Solve quantile regression using Frisch-Newton exact LP method.
 *
 * This is an exact LP solver (not smoothed), enabling preprocessing.
 * The algorithm:
 * 1. Initialize β via OLS
 * 2. Identify active set (observations with |residual| < tolerance)
 * 3. Solve weighted normal equations on active set
 * 4. Update β with ratio test for step length
 * 5. Check optimality: X'g = 0 where g = subgradient
 * 6. If not optimal, expand active set and repeat
 *
 * Parameters:
 *   ipm     - Pre-allocated IPM state (reuses arrays)
 *   y       - Dependent variable (N)
 *   X       - Design matrix (N x K, column-major)
 *   q       - Quantile (0 < q < 1)
 *   beta    - Output: coefficient estimates (K)
 *
 * Returns:
 *   Number of iterations on success (positive), negative on failure.
 */
ST_int cqreg_fn_solve(cqreg_ipm_state *ipm,
                       const ST_double *y,
                       const ST_double *X,
                       ST_double q,
                       ST_double *beta);

/*
 * Solve quantile regression with warm-start initialization.
 *
 * Same as cqreg_fn_solve but uses initial_beta as starting point instead
 * of computing OLS. This is much faster when the initial estimate is close
 * to the true solution (e.g., for auxiliary solves at τ±h in Siddiqui method).
 *
 * Parameters:
 *   ipm          - Pre-allocated IPM state
 *   y            - Dependent variable (N)
 *   X            - Design matrix (N x K, column-major)
 *   q            - Quantile (0 < q < 1)
 *   initial_beta - Initial coefficient estimate for warm-start (K)
 *   beta         - Output: coefficient estimates (K)
 *
 * Returns:
 *   Number of iterations on success (positive), negative on failure.
 */
ST_int cqreg_fn_solve_warmstart(cqreg_ipm_state *ipm,
                                 const ST_double *y,
                                 const ST_double *X,
                                 ST_double q,
                                 const ST_double *initial_beta,
                                 ST_double *beta);

/*
 * Solve QR with preprocessing (Chernozhukov et al. algorithm).
 *
 * Uses initial subsample of size m = (N * (log K + 1))^{2/3}
 * then iteratively expands until optimal on full sample.
 *
 * Parameters: Same as cqreg_fn_solve()
 * Returns: Total iterations across all subsamples.
 */
ST_int cqreg_fn_preprocess_solve(cqreg_ipm_state *ipm,
                                  const ST_double *y,
                                  const ST_double *X,
                                  ST_double q,
                                  ST_double *beta);

/* ============================================================================
 * Internal Functions (exposed for testing)
 * ============================================================================ */

/*
 * Identify active set: observations with |residual| < tolerance.
 *
 * Parameters:
 *   r           - Residuals (N)
 *   N           - Number of observations
 *   tol         - Tolerance for active identification
 *   active_idx  - Output: indices of active observations (size N buffer)
 *
 * Returns:
 *   Number of active observations (m).
 */
ST_int fn_identify_active_set(const ST_double *r, ST_int N, ST_double tol,
                               ST_int *active_idx);

/*
 * Form reduced normal equations on active set.
 *
 * Computes X_A'X_A and X_A'g for observations in active set.
 *
 * Parameters:
 *   X           - Full design matrix (N x K, column-major)
 *   g           - Subgradient vector (N)
 *   N, K        - Dimensions
 *   active_idx  - Indices of active observations
 *   m           - Number of active observations
 *   XAX         - Output: X_A'X_A (K x K)
 *   XAg         - Output: X_A'g (K)
 */
void fn_form_reduced_system(const ST_double *X, const ST_double *g,
                             ST_int N, ST_int K,
                             const ST_int *active_idx, ST_int m,
                             ST_double *XAX, ST_double *XAg);

/*
 * Compute step length via ratio test.
 *
 * Finds maximum α such that u + α*Δu ≥ 0 and v + α*Δv ≥ 0.
 *
 * Parameters:
 *   r           - Current residuals (N)
 *   delta_r     - Change in residuals from step (N)
 *   N           - Number of observations
 *   q           - Quantile
 *
 * Returns:
 *   Step length α in (0, 1].
 */
ST_double fn_ratio_test(const ST_double *r, const ST_double *delta_r,
                         ST_int N, ST_double q);

/*
 * Check optimality on full sample.
 *
 * Verifies X'g = 0 where g_i = q for r_i > 0, g_i = q-1 for r_i < 0.
 *
 * Parameters:
 *   X           - Design matrix (N x K, column-major)
 *   r           - Residuals (N)
 *   N, K        - Dimensions
 *   q           - Quantile
 *   Xg          - Output: X'g (K) - caller provides buffer
 *
 * Returns:
 *   1 if optimal (||X'g|| < tolerance), 0 otherwise.
 */
ST_int fn_check_optimality(const ST_double *X, const ST_double *r,
                            ST_int N, ST_int K, ST_double q,
                            ST_double *Xg);

#endif /* CQREG_FN_H */

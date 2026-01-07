/*
 * cqreg_ipm.h
 *
 * Interior Point Method solver for quantile regression.
 * Implements primal-dual IPM with Mehrotra predictor-corrector.
 * Part of the ctools suite.
 */

#ifndef CQREG_IPM_H
#define CQREG_IPM_H

#include "cqreg_types.h"

/* ============================================================================
 * Main Solver Interface
 * ============================================================================ */

/*
 * Solve quantile regression problem using Interior Point Method.
 *
 * Minimizes: sum_i [ q * u_i + (1-q) * v_i ]
 * Subject to: y - X * beta = u - v, u >= 0, v >= 0
 *
 * Parameters:
 *   ipm     - Pre-allocated IPM state (from cqreg_ipm_create)
 *   y       - Dependent variable (N)
 *   X       - Design matrix (N x K, column-major)
 *   q       - Quantile (0 < q < 1)
 *   beta    - Output: coefficient estimates (K)
 *
 * Returns:
 *   Number of iterations on success, negative on failure.
 */
ST_int cqreg_ipm_solve(cqreg_ipm_state *ipm,
                       const ST_double *y,
                       const ST_double *X,
                       ST_double q,
                       ST_double *beta);

/*
 * Get the objective value (sum of weighted absolute deviations).
 * Call after cqreg_ipm_solve().
 */
ST_double cqreg_ipm_get_objective(const cqreg_ipm_state *ipm, ST_double q);

/*
 * Get residuals: r = y - X * beta = u - v
 * Call after cqreg_ipm_solve().
 *
 * Parameters:
 *   ipm      - Solved IPM state
 *   residuals - Output array (N)
 */
void cqreg_ipm_get_residuals(const cqreg_ipm_state *ipm, ST_double *residuals);

/* ============================================================================
 * Internal Functions (exposed for testing)
 * ============================================================================ */

/*
 * Initialize primal and dual variables for IPM.
 * Uses a simple heuristic based on least squares solution.
 */
void cqreg_ipm_initialize(cqreg_ipm_state *ipm,
                          const ST_double *y,
                          const ST_double *X,
                          ST_double q);

/*
 * Compute primal and dual residuals.
 */
void cqreg_ipm_compute_residuals(cqreg_ipm_state *ipm,
                                 const ST_double *y,
                                 const ST_double *X,
                                 ST_double q);

/*
 * Check convergence based on residuals and duality gap.
 * Returns: 1 if converged, 0 otherwise.
 */
ST_int cqreg_ipm_check_convergence(cqreg_ipm_state *ipm);

/*
 * Compute diagonal scaling matrix D.
 * D[i] = 1 / (lambda_u[i]/u[i] + lambda_v[i]/v[i])
 */
void cqreg_ipm_compute_scaling(cqreg_ipm_state *ipm);

/*
 * Form and factorize normal equations.
 * Computes X' * D * X and its Cholesky factorization.
 * Returns: 0 on success, -1 if Cholesky fails.
 */
ST_int cqreg_ipm_form_normal_equations(cqreg_ipm_state *ipm, const ST_double *X);

/*
 * Compute affine scaling direction (predictor step).
 */
void cqreg_ipm_affine_direction(cqreg_ipm_state *ipm,
                                const ST_double *X,
                                const ST_double *y,
                                ST_double q);

/*
 * Compute combined (predictor-corrector) direction.
 */
void cqreg_ipm_combined_direction(cqreg_ipm_state *ipm,
                                  const ST_double *X,
                                  const ST_double *y,
                                  ST_double q,
                                  ST_double sigma);

/*
 * Compute step length maintaining positivity of u, v, lambda_u, lambda_v.
 * Returns: maximum step length in (0, 1].
 */
ST_double cqreg_ipm_step_length(cqreg_ipm_state *ipm, ST_int affine);

/*
 * Update variables: x = x + alpha * delta_x.
 */
void cqreg_ipm_update_variables(cqreg_ipm_state *ipm, ST_double alpha);

/*
 * Compute complementarity: sum(u .* lambda_u + v .* lambda_v)
 */
ST_double cqreg_ipm_complementarity(const cqreg_ipm_state *ipm);

/* ============================================================================
 * IRLS Alternative Solver
 * ============================================================================ */

/*
 * Solve quantile regression using Iteratively Reweighted Least Squares.
 * Simpler and more robust than IPM for small/medium datasets.
 *
 * Parameters: Same as cqreg_ipm_solve()
 * Returns: Number of iterations on success, negative on failure.
 */
ST_int cqreg_irls_solve(cqreg_ipm_state *ipm,
                        const ST_double *y,
                        const ST_double *X,
                        ST_double q,
                        ST_double *beta);

/* ============================================================================
 * Preprocessing Solver (Chernozhukov et al. 2020)
 * ============================================================================ */

/*
 * Solve quantile regression using preprocessing algorithm.
 *
 * This is the exact algorithm from:
 *   Chernozhukov, Fernandez-Val, Melly (2020)
 *   "Fast Algorithms for the Quantile Regression Process"
 *
 * Algorithm:
 * 1. Get initial beta via OLS
 * 2. Predict residual signs
 * 3. Select reduced subsample of size m = O((n*log(k))^{2/3})
 * 4. Solve QR on subsample
 * 5. Check optimality on full sample
 * 6. Add violated observations, repeat until optimal
 *
 * This is EXACT (not approximate) - numerically identical to full solve.
 *
 * Parameters:
 *   ipm     - Pre-allocated IPM state (from cqreg_ipm_create)
 *   y       - Dependent variable (N)
 *   X       - Design matrix (N x K, column-major)
 *   q       - Quantile (0 < q < 1)
 *   beta    - Output: coefficient estimates (K)
 *
 * Returns:
 *   Number of total iterations on success, negative on failure.
 */
ST_int cqreg_preprocess_solve(cqreg_ipm_state *ipm,
                               const ST_double *y,
                               const ST_double *X,
                               ST_double q,
                               ST_double *beta);

#endif /* CQREG_IPM_H */

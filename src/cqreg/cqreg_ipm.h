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

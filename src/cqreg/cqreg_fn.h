/*
 * cqreg_fn.h
 *
 * Frisch-Newton interior point method solver for quantile regression.
 * Based on the primal-dual IPM with Mehrotra predictor-corrector.
 */

#ifndef CQREG_FN_H
#define CQREG_FN_H

#include "cqreg_types.h"

/* ============================================================================
 * Main Solver Interface
 * ============================================================================ */

/*
 * Solve quantile regression using Frisch-Newton interior point method.
 *
 * The algorithm uses primal-dual IPM with Mehrotra predictor-corrector:
 * 1. Initialize from OLS solution
 * 2. Iterate with predictor-corrector steps
 * 3. Converge when duality gap < tolerance
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

#endif /* CQREG_FN_H */

/*
 * cbinscatter_resid.h
 *
 * Residualization functions for cbinscatter
 * Supports OLS residualization (controls) and HDFE (absorb)
 * Part of the ctools Stata plugin suite
 */

#ifndef CBINSCATTER_RESID_H
#define CBINSCATTER_RESID_H

#include "cbinscatter_types.h"

/* Rank-aware normal-equation solve; redundant coefficients are set to zero.
 * Nonzero return means failure, never a request to skip adjustment. */
ST_int cbinscatter_solve_with_collinearity(const ST_double *XtX,
    const ST_double *rhs, ST_int K, ST_double *beta);

/* ========================================================================
 * OLS Residualization (Controls Only)
 * ======================================================================== */

/*
 * Residualize y and x on control variables using OLS
 * y_resid = y - X*(X'X)^-1*X'y
 * x_resid = x - X*(X'X)^-1*X'x
 *
 * Parameters:
 *   y           - y values (N), modified in-place
 *   x           - x values (N), modified in-place
 *   controls    - control variables (N x K), column-major
 *   N           - number of observations
 *   K           - number of control variables
 *   weights     - observation weights (NULL if unweighted)
 *   weight_type - 0=none, 1=aweight, 2=fweight, 3=pweight
 *
 * Returns:
 *   0 on success, error code on failure
 */
ST_retcode ols_residualize(
    ST_double *y,
    ST_double *x,
    const ST_double *controls,
    ST_int N,
    ST_int K,
    const ST_double *weights,
    ST_int weight_type
);

/* ========================================================================
 * HDFE Residualization (Absorb Option)
 * ======================================================================== */

/*
 * Residualize y and x using high-dimensional fixed effects
 * One factor is a single exact projection; several use the shared CG solver
 * (src/hdfe) with symmetric Kaczmarz transformations
 *
 * Parameters:
 *   y           - y values (N), modified in-place
 *   x           - x values (N), modified in-place
 *   fe_vars     - fixed effect variables (N x G), column-major, integer levels
 *   N           - number of observations
 *   G           - number of FE groups
 *   weights     - observation weights (NULL if unweighted)
 *   weight_type - 0=none, 1=aweight, 2=fweight, 3=pweight
 *   maxiter     - maximum CG iterations
 *   tolerance   - convergence tolerance
 *   verbose     - print iteration info
 *   dropped     - output: number of singleton observations dropped
 *
 * Returns:
 *   0 on success, CBINSCATTER_ERR_HDFE when the projections have not
 *   converged after maxiter iterations, another error code on failure
 */
ST_retcode hdfe_residualize(
    ST_double *y,
    ST_double *x,
    const ST_int *fe_vars,
    ST_int N,
    ST_int G,
    const ST_double *weights,
    ST_int weight_type,
    ST_int maxiter,
    ST_double tolerance,
    ST_int verbose,
    ST_int *dropped
);

/* ========================================================================
 * Combined Residualization (Controls + Absorb)
 * ======================================================================== */

/*
 * Residualize y and x on both controls and FE
 * Strategy: First HDFE residualize all (y, x, controls),
 *           then OLS residualize partialled y,x on partialled controls
 *
 * Parameters:
 *   y           - y values (N), modified in-place
 *   x           - x values (N), modified in-place
 *   controls    - control variables (N x K_ctrl), column-major
 *   fe_vars     - fixed effect variables (N x G), integer levels
 *   N           - number of observations
 *   K_ctrl      - number of control variables
 *   G           - number of FE groups
 *   weights     - observation weights (NULL if unweighted)
 *   weight_type - 0=none, 1=aweight, 2=fweight, 3=pweight
 *   maxiter     - maximum CG iterations
 *   tolerance   - convergence tolerance
 *   verbose     - print iteration info
 *   dropped     - output: number of singleton observations dropped
 *
 * Returns:
 *   0 on success, CBINSCATTER_ERR_HDFE when the projections have not
 *   converged after maxiter iterations, another error code on failure
 */
ST_retcode combined_residualize(
    ST_double *y,
    ST_double *x,
    ST_double *controls,
    const ST_int *fe_vars,
    ST_int N,
    ST_int K_ctrl,
    ST_int G,
    const ST_double *weights,
    ST_int weight_type,
    ST_int maxiter,
    ST_double tolerance,
    ST_int verbose,
    ST_int *dropped
);

/* ========================================================================
 * Batch HDFE Residualization
 * ======================================================================== */

/*
 * Residualize multiple variables with the same absorbed effects. One factor
 * is a single exact projection; several use the shared conjugate-gradient
 * solver (src/hdfe), each column converging to tolerance as in creghdfe.
 * Returns CBINSCATTER_ERR_HDFE if a column does not converge within maxiter
 * iterations. Failure leaves working columns modified; callers must not
 * publish them.
 *
 * Parameters:
 *   vars        - K_vars pointers to N-length arrays, each modified in-place
 *   K_vars      - number of variables to residualize
 *   fe_vars     - fixed effect variables (N x G), integer levels
 *   N           - number of observations
 *   G           - number of FE groups
 *   weights     - observation weights (NULL if unweighted)
 *   weight_type - 0=none, 1=aweight, 2=fweight, 3=pweight
 *   maxiter     - maximum iterations
 *   tolerance   - convergence tolerance
 *   verbose     - print iteration info
 *   dropped     - output: number of singleton observations dropped
 */
ST_retcode hdfe_residualize_batch(
    ST_double **vars,
    ST_int K_vars,
    const ST_int *fe_vars,
    ST_int N,
    ST_int G,
    const ST_double *weights,
    ST_int weight_type,
    ST_int maxiter,
    ST_double tolerance,
    ST_int verbose,
    ST_int *dropped
);

/*
 * Weighted normal equations for R regressors given as column pointers:
 * XtX[j * R + k] (k >= j, upper triangle only) = sum_i (w_i x_ij) x_ik and
 * Xty[j] = sum_i (w_i x_ij) y_i, each summed over the observations in order
 * (weights NULL = 1). When bin_ids is not NULL, rows whose bin id is outside
 * 1..num_bins are skipped. The rows of the triangle are independent sums and
 * are computed in parallel. XtX and Xty are overwritten.
 */
void cbinscatter_normal_equations(
    const ST_double *const *cols,
    ST_int R,
    const ST_double *y,
    const ST_double *weights,
    const ST_int *bin_ids,
    ST_int num_bins,
    ST_int N,
    ST_double *XtX,
    ST_double *Xty
);

#endif /* CBINSCATTER_RESID_H */

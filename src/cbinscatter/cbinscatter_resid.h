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
 * Uses CG solver with symmetric Kaczmarz transformations
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
 *   0 on success, error code on failure
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
 *   0 on success, error code on failure
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
 * Residualize multiple variables simultaneously using HDFE.
 * More efficient than calling hdfe_residualize_y_only in a loop because
 * FE level arrays are traversed once per sweep instead of once per variable.
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

#endif /* CBINSCATTER_RESID_H */

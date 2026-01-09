/*
 * cbinscatter_fit.h
 *
 * Line fitting functions for cbinscatter
 * Part of the ctools Stata plugin suite
 */

#ifndef CBINSCATTER_FIT_H
#define CBINSCATTER_FIT_H

#include "cbinscatter_types.h"

/* ========================================================================
 * Line Fitting
 * ======================================================================== */

/*
 * Fit polynomial to scatter data (not bin means)
 *
 * Parameters:
 *   y           - y values
 *   x           - x values
 *   weights     - observation weights (NULL if unweighted)
 *   N           - number of observations
 *   order       - polynomial order (1=linear, 2=quadratic, 3=cubic)
 *   coefs       - output: coefficients (order+1 values: const, x, x^2, ...)
 *   r2          - output: R-squared
 *
 * Returns:
 *   0 on success, error code on failure
 */
ST_retcode fit_polynomial(
    const ST_double *y,
    const ST_double *x,
    const ST_double *weights,
    ST_int N,
    ST_int order,
    ST_double *coefs,
    ST_double *r2
);

/*
 * Fit lines for all by-groups in results
 *
 * Parameters:
 *   y           - y values (full dataset)
 *   x           - x values (full dataset)
 *   by_groups   - by-group assignments (NULL if no by)
 *   weights     - weights (NULL if unweighted)
 *   N           - total observations
 *   config      - binscatter configuration
 *   results     - results structure (groups will be updated with fit coefs)
 *
 * Returns:
 *   0 on success, error code on failure
 */
ST_retcode fit_all_groups(
    const ST_double *y,
    const ST_double *x,
    const ST_int *by_groups,
    const ST_double *weights,
    ST_int N,
    const BinscatterConfig *config,
    BinscatterResults *results
);

/*
 * Evaluate polynomial at given x values
 *
 * Parameters:
 *   x           - x values to evaluate at
 *   N           - number of x values
 *   coefs       - polynomial coefficients
 *   order       - polynomial order
 *   fitted      - output: fitted values
 */
void eval_polynomial(
    const ST_double *x,
    ST_int N,
    const ST_double *coefs,
    ST_int order,
    ST_double *fitted
);

#endif /* CBINSCATTER_FIT_H */

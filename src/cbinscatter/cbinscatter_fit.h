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

#endif /* CBINSCATTER_FIT_H */

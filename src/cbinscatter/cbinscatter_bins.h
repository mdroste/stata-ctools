/*
 * cbinscatter_bins.h
 *
 * Bin computation functions for cbinscatter
 * Part of the ctools Stata plugin suite
 */

#ifndef CBINSCATTER_BINS_H
#define CBINSCATTER_BINS_H

#include "cbinscatter_types.h"

/* ========================================================================
 * Quantile Bin Computation
 * ======================================================================== */

/*
 * Compute quantile cutpoints for binning
 *
 * Parameters:
 *   x_sorted    - x values sorted in ascending order
 *   N           - number of observations
 *   nquantiles  - number of bins
 *   weights     - observation weights (NULL if unweighted)
 *   w_sorted    - weights in same order as x_sorted (NULL if unweighted)
 *   cutpoints   - output array of size (nquantiles + 1)
 *
 * Returns:
 *   0 on success, error code on failure
 */
ST_retcode compute_quantile_cutpoints(
    const ST_double *x_sorted,
    ST_int N,
    ST_int nquantiles,
    const ST_double *w_sorted,
    ST_double *cutpoints
);

/*
 * Assign bin IDs based on cutpoints
 * Ties go to lower bin (consistent with Stata binscatter)
 *
 * Parameters:
 *   x           - x values (not sorted)
 *   N           - number of observations
 *   cutpoints   - quantile cutpoints (size nquantiles + 1)
 *   nquantiles  - number of bins
 *   bin_ids     - output: bin assignment for each observation (1-based)
 */
void assign_bins(
    const ST_double *x,
    ST_int N,
    const ST_double *cutpoints,
    ST_int nquantiles,
    ST_int *bin_ids
);

/* ========================================================================
 * Bin Statistics Computation
 * ======================================================================== */

/*
 * Compute bin means (and optionally SE)
 *
 * Parameters:
 *   y           - y values
 *   x           - x values
 *   bin_ids     - bin assignments (1-based)
 *   weights     - observation weights (NULL if unweighted)
 *   N           - number of observations
 *   nquantiles  - number of bins
 *   compute_se  - whether to compute standard errors
 *   result      - output: ByGroupResult to populate
 *
 * Returns:
 *   0 on success, error code on failure
 */
ST_retcode compute_bin_statistics(
    const ST_double *y,
    const ST_double *x,
    const ST_int *bin_ids,
    const ST_double *weights,
    ST_int N,
    ST_int nquantiles,
    ST_int compute_se,
    ByGroupResult *result
);

/* ========================================================================
 * Single Group Processing
 * ======================================================================== */

/*
 * Compute bins for a single by-group (or entire dataset if no by())
 *
 * Parameters:
 *   y           - y values for this group
 *   x           - x values for this group
 *   weights     - weights for this group (NULL if unweighted)
 *   N           - number of observations in group
 *   config      - binscatter configuration
 *   result      - output: ByGroupResult to populate
 *
 * Returns:
 *   0 on success, error code on failure
 */
ST_retcode compute_bins_single_group(
    ST_double *y,
    ST_double *x,
    ST_double *weights,
    ST_int N,
    const BinscatterConfig *config,
    ByGroupResult *result
);

/* ========================================================================
 * Discrete Mode
 * ======================================================================== */

/*
 * Compute bins in discrete mode (one bin per unique x value)
 *
 * Parameters:
 *   y           - y values
 *   x           - x values
 *   weights     - weights (NULL if unweighted)
 *   N           - number of observations
 *   compute_se  - whether to compute standard errors
 *   result      - output: ByGroupResult to populate
 *
 * Returns:
 *   0 on success, error code on failure
 */
ST_retcode compute_bins_discrete(
    const ST_double *y,
    const ST_double *x,
    const ST_double *weights,
    ST_int N,
    ST_int compute_se,
    ByGroupResult *result
);

#endif /* CBINSCATTER_BINS_H */

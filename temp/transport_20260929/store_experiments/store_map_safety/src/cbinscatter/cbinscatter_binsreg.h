/*
 * cbinscatter_binsreg.h
 *
 * Binsreg method implementation for cbinscatter
 * Implements the Cattaneo et al. "On Binscatter" approach
 * Part of the ctools Stata plugin suite
 */

#ifndef CBINSCATTER_BINSREG_H
#define CBINSCATTER_BINSREG_H

#include "cbinscatter_types.h"

/*
 * Adjust bin Y-means using binsreg regression method
 *
 * This function implements the correct binsreg approach when controls are present:
 * Instead of: mean(Y_resid | bin) where Y_resid = Y - W*gamma_pooled
 * We compute: beta_j from reg Y = sum(beta_j * I(bin=j)) + W*gamma_partial + e
 *
 * Algorithm (Frisch-Waugh-Lovell for categorical regressors):
 * 1. Compute within-bin means of Y and each control: E[Y|bin], E[W_k|bin]
 * 2. Create demeaned versions: Y_dm = Y - E[Y|bin], W_dm = W - E[W|bin]
 * 3. Run OLS: Y_dm ~ W_dm to get gamma_partial
 * 4. Update bin y_mean: beta_j = E[Y|bin=j] - sum_k(E[W_k|bin=j] * gamma_k)
 *
 * Parameters:
 *   y           - y values (NOT residualized)
 *   controls    - control variables (column-major: K x N)
 *   bin_ids     - bin assignments (1-based)
 *   weights     - observation weights (NULL if unweighted)
 *   N           - number of observations
 *   K           - number of control variables
 *   result      - ByGroupResult with bins to update (bins[].y_mean will be modified)
 *
 * Returns:
 *   0 on success, error code on failure
 */
ST_retcode adjust_bins_binsreg(
    const ST_double *y,
    const ST_double *controls,
    const ST_int *bin_ids,
    const ST_double *weights,
    ST_int N,
    ST_int K,
    ByGroupResult *result
);

/*
 * Adjust bin Y-means using binsreg method with fixed effects
 *
 * For absorb(), this first HDFE-demeans Y and controls, then applies
 * the binsreg adjustment.
 *
 * Parameters:
 *   y           - y values (NOT residualized)
 *   controls    - control variables (column-major: K_ctrl x N), can be NULL
 *   fe_vars     - fixed effect group IDs (column-major: K_fe x N)
 *   bin_ids     - bin assignments (1-based)
 *   weights     - observation weights (NULL if unweighted)
 *   N           - number of observations
 *   K_ctrl      - number of control variables (can be 0)
 *   K_fe        - number of fixed effect variables
 *   weight_type - type of weights (0=none, 1=aweight, 2=fweight, 3=pweight)
 *   maxiter     - maximum iterations for HDFE convergence
 *   tolerance   - convergence tolerance for HDFE
 *   result      - ByGroupResult with bins to update
 *
 * Returns:
 *   0 on success, error code on failure
 */
ST_retcode adjust_bins_binsreg_hdfe(
    ST_double *y,
    ST_double *controls,
    ST_int *fe_vars,
    const ST_int *bin_ids,
    const ST_double *weights,
    ST_int N,
    ST_int K_ctrl,
    ST_int K_fe,
    ST_int weight_type,
    ST_int maxiter,
    ST_double tolerance,
    ByGroupResult *result
);

#endif /* CBINSCATTER_BINSREG_H */

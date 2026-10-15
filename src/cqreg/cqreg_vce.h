/*
 * cqreg_vce.h
 *
 * Variance-covariance estimation for quantile regression.
 * Implements IID, robust (sandwich), and clustered VCE.
 * Part of the ctools suite.
 */

#ifndef CQREG_VCE_H
#define CQREG_VCE_H

#include "cqreg_types.h"

/* ============================================================================
 * Main Interface
 * ============================================================================ */

/*
 * Compute variance-covariance matrix based on VCE type.
 * Dispatches to the appropriate estimator.
 *
 * Parameters:
 *   state - cqreg state with regression results
 *   X     - Design matrix (N x K, column-major)
 *
 * Returns:
 *   0 on success, 920 when a sandwich temporary cannot be allocated,
 *   another nonzero code otherwise
 */
ST_int cqreg_compute_vce(cqreg_state *state, const ST_double *X);

#endif /* CQREG_VCE_H */

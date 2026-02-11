/*
 * cqreg_sparsity.h
 *
 * Sparsity estimation for quantile regression variance computation.
 * Implements bandwidth selection and kernel density estimation.
 * Part of the ctools suite.
 */

#ifndef CQREG_SPARSITY_H
#define CQREG_SPARSITY_H

#include "cqreg_types.h"

/* ============================================================================
 * Main Interface
 * ============================================================================ */

/*
 * Estimate sparsity (1/f(F^{-1}(q))) from regression residuals.
 *
 * Parameters:
 *   sp        - Pre-allocated sparsity state
 *   residuals - Regression residuals (N)
 *
 * Returns:
 *   Estimated sparsity value (stored in sp->sparsity)
 */
ST_double cqreg_estimate_sparsity(cqreg_sparsity_state *sp,
                                  const ST_double *residuals);

/* ============================================================================
 * Bandwidth Selection
 * ============================================================================ */

/*
 * Compute bandwidth using the specified method.
 */
ST_double cqreg_compute_bandwidth(ST_int N, ST_double q, cqreg_bw_method method);


#endif /* CQREG_SPARSITY_H */

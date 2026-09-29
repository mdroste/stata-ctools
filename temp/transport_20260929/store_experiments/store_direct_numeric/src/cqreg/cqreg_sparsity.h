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
 * Estimate sparsity (1/f(F^{-1}(q))) from regression residuals by qreg's
 * residual method: the difference quotient of the residual percentiles after
 * dropping the sp->n_basis observations listed in sp->basis (the LP basis).
 *
 * Parameters:
 *   sp        - Pre-allocated sparsity state (basis, n_basis and, if
 *               positive, full_bandwidth set by the caller)
 *   residuals - Regression residuals (N)
 *
 * Returns:
 *   Estimated sparsity value (stored in sp->sparsity). sp->status is 0, or the
 *   Stata error code (498, 2000) when qreg would stop with an error; the
 *   message has then been displayed.
 */
ST_double cqreg_estimate_sparsity(cqreg_sparsity_state *sp,
                                  const ST_double *residuals);

/*
 * qreg's kernel density method: kernel bandwidth sp->kbwidth from the
 * residuals and sp->kfactor, kernel weights kval_i = K(r_i/kbwidth) for the
 * robust Hessian, and the iid sparsity kbwidth / mean(kval), which is
 * returned (infinite if every weight is zero). sp->status as above.
 */
ST_double cqreg_estimate_kernel_sparsity(cqreg_sparsity_state *sp,
                                         const ST_double *residuals,
                                         ST_double *kval);

/* ============================================================================
 * Bandwidth Selection
 * ============================================================================ */

/*
 * Compute bandwidth using the specified method.
 */
ST_double cqreg_compute_bandwidth(ST_int N, ST_double q, cqreg_bw_method method);


#endif /* CQREG_SPARSITY_H */

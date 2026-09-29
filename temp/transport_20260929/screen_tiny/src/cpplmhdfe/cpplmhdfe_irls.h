/*
 * cpplmhdfe_irls.h
 *
 * IRLS (Iteratively Reweighted Least Squares) loop for PPML estimation
 * Part of the ctools Stata plugin suite
 */

#ifndef CPPLMHDFE_IRLS_H
#define CPPLMHDFE_IRLS_H

#include "cpplmhdfe_types.h"

/*
 * Main PPML regression function.
 *
 * Performs:
 *   1. Data loading and validation (y >= 0)
 *   2. FE setup: remap, singleton removal, DOF
 *   3. Separation detection
 *   4. IRLS loop with CG solver for FE demeaning
 *   5. Collinearity detection
 *   6. VCE computation (robust/cluster)
 *   7. Result storage
 *
 * Expected scalars from Stata:
 *   __cpplmhdfe_K, __cpplmhdfe_G, vcetype, weights, etc.
 *
 * Returns:
 *   0 on success, Stata error code on failure
 */
ST_retcode do_ppml_regression(int argc, char *argv[]);

#endif /* CPPLMHDFE_IRLS_H */

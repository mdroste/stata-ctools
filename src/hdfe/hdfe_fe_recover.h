/*
 * hdfe_fe_recover.h
 *
 * Fixed-effect coefficients from the fixed-effect component of a fit,
 * without the Stata SPI (callers load and store the variables).
 * Part of the ctools Stata plugin suite
 */

#ifndef CTOOLS_HDFE_FE_RECOVER_H
#define CTOOLS_HDFE_FE_RECOVER_H

#include "hdfe_types.h"

/*
 * Recover the additive fixed effects alpha from d = y - Xb - _cons - resid
 * (the N rows of S, in S's row order), as reghdfe's store_alphas does for its
 * default solver: forward Kaczmarz sweeps over the effects accelerated by
 * steepest descent, with each sweep's level means accumulated into the
 * coefficients. Starting from zero, every effect after the first then has a
 * zero (weighted) sum within each connected component it forms with the
 * effect before it, which fixes the normalization reghdfe reports. The
 * iteration stops when the weighted mean of |step| / (|d| + scale) is at
 * most S->tolerance, with scale = max(1, max |d|), and fails after
 * S->maxiter steps.
 *
 * error: the N values of d, overwritten. alpha: G arrays of
 * S->factors[g].num_levels zeros, filled with the coefficients; levels with
 * no weight stay zero. Returns 0, 1 (cancelled), 920 (allocation), or
 * 430 (no convergence).
 */
ST_retcode hdfe_fe_recover(const HDFE_State *S, ST_double *error, ST_double **alpha);

#endif /* CTOOLS_HDFE_FE_RECOVER_H */

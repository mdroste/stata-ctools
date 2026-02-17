/*
 * cpplmhdfe_separation.h
 *
 * Separation detection for PPML estimation
 * Part of the ctools Stata plugin suite
 */

#ifndef CPPLMHDFE_SEPARATION_H
#define CPPLMHDFE_SEPARATION_H

#include "cpplmhdfe_types.h"

/*
 * Detect separation via FE-level screening.
 * For each FE level, if all y_i = 0 for observations in that level,
 * mark those observations as separated.
 *
 * Parameters:
 *   y:          Dependent variable array (N)
 *   factors:    Array of FE_Factor structs (G factors)
 *   G:          Number of FE factors
 *   N:          Number of observations
 *   sep_mask:   Output mask (1 = separated, 0 = not). Must be pre-allocated to N.
 *
 * Returns: Number of separated observations found
 */
ST_int ppml_detect_separation_fe(
    const ST_double *y,
    const FE_Factor *factors,
    ST_int G,
    ST_int N,
    ST_int *sep_mask
);

/*
 * Detect separation via mu-based screening (post-iteration).
 * For observations where y=0 and mu < sep_tol, mark as separated.
 *
 * Parameters:
 *   y:          Dependent variable (N)
 *   mu:         Fitted values (N)
 *   sep_tol:    Separation tolerance
 *   N:          Number of observations
 *   sep_mask:   Output mask (1 = separated, 0 = not). Must be pre-allocated to N.
 *
 * Returns: Number of newly separated observations found
 */
ST_int ppml_detect_separation_mu(
    const ST_double *y,
    const ST_double *mu,
    ST_double sep_tol,
    ST_int N,
    ST_int *sep_mask
);

#endif /* CPPLMHDFE_SEPARATION_H */

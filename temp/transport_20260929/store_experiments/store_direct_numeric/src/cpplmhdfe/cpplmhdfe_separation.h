/*
 * cpplmhdfe_separation.h
 *
 * Separation detection for PPML estimation
 * Part of the ctools Stata plugin suite
 */

#ifndef CPPLMHDFE_SEPARATION_H
#define CPPLMHDFE_SEPARATION_H

#include "cpplmhdfe_types.h"

/* Select on original rows until singleton and all-zero FE removals stabilize.
 * levels are 1-based, mask is initialized by the caller; -1 means allocation
 * or invalid-level failure. No compact/original index spaces are mixed. */
ST_int ppml_select_fe_sample(const ST_double *y, ST_int *const *levels,
    const ST_int *num_levels, ST_int G, ST_int N, ST_int *mask,
    ST_int *num_singletons, ST_int *num_separated);

ST_int ppml_select_fe_sample_ex(const ST_double *y, ST_int *const *levels,
    const ST_int *num_levels, ST_int G, ST_int N, ST_int *mask,
    ST_int *num_singletons, ST_int *num_separated,
    const ST_int *slope_only, ST_int keep_singletons, ST_int check_separation);

#endif /* CPPLMHDFE_SEPARATION_H */

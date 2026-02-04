/*
 * cpsmatch_impl.h
 *
 * C-accelerated propensity score matching for Stata
 * Part of the ctools suite
 *
 * Implements fast matching algorithms similar to psmatch2:
 *   - Nearest neighbor matching (with/without replacement)
 *   - Radius/caliper matching
 *   - Kernel matching
 *   - Mahalanobis distance matching
 */

#ifndef CPSMATCH_IMPL_H
#define CPSMATCH_IMPL_H

#include "stplugin.h"

/*
 * Main entry point for cpsmatch command.
 *
 * @param args  Command arguments as space-separated string
 *              Format: "treat_idx pscore_idx outcome_idx out_weight_idx out_match_idx
 *                       nobs method neighbor caliper common kernel_type bwidth
 *                       with_replace ties descending"
 *
 * @return      0 on success, Stata error code on failure
 */
ST_retcode cpsmatch_main(const char *args);

#endif /* CPSMATCH_IMPL_H */

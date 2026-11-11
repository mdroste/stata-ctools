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
 *   - with ate, ai() and psmatch2's _id/_n# (Mahalanobis matching runs in Mata)
 */

#ifndef CPSMATCH_IMPL_H
#define CPSMATCH_IMPL_H

#include "stplugin.h"

/*
 * Main entry point for cpsmatch command.
 *
 * @param args  Command arguments as space-separated string
 *              Format: "treat_idx pscore_idx outcome_idx out_weight_idx id_idx
 *                       support_idx nobs method neighbor caliper common kernel_type
 *                       bwidth with_replace ties descending", then optional flags:
 *              presup= arm= nnvar= ids= pdif= matched= self= (plugin variable
 *              positions), ai=#, ate, altvar, w2float, stats, touse
 *
 * @return      0 on success, Stata error code on failure
 */
ST_retcode cpsmatch_main(const char *args);

#endif /* CPSMATCH_IMPL_H */

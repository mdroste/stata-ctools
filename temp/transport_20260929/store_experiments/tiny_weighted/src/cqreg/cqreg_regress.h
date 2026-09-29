/*
 * cqreg_regress.h
 *
 * Full quantile regression command orchestrator.
 * Coordinates data loading, IPM solver, VCE, and result storage.
 * Part of the ctools suite.
 */

#ifndef CQREG_REGRESS_H
#define CQREG_REGRESS_H

#include "stplugin.h"

/*
 * Execute full quantile regression.
 *
 * This is the main entry point called from cqreg_impl.c.
 * It orchestrates the complete regression:
 *   1. Parse arguments from Stata scalars
 *   2. Load data (y, X, optionally weights and cluster IDs)
 *   3. Run IPM solver
 *   4. Estimate sparsity and compute VCE
 *   5. Store results back to Stata scalars/matrices
 *
 * Arguments from Stata scalars:
 *   __cqreg_quantile     - Quantile (0 < q < 1)
 *   __cqreg_K            - Number of regressors (including constant)
 *   __cqreg_G            - Must be 0; absorb() is rejected
 *   __cqreg_vce_type     - 0=iid, 1=robust, 2=cluster
 *   __cqreg_bw_method    - 0=hsheather, 1=bofinger, 2=chamberlain
 *   __cqreg_verbose      - Verbosity level
 *   __cqreg_tolerance    - IPM tolerance
 *   __cqreg_maxiter      - IPM max iterations
 *
 * Variables passed via plugin call:
 *   First K vars: depvar, indepvars (constant added internally)
 *   Next 1 var: cluster variable (if vce_type == 2)
 *
 * Returns:
 *   0 on success, Stata error code otherwise
 */
ST_retcode cqreg_full_regression(const char *args);

#endif /* CQREG_REGRESS_H */

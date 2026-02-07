/*
 * creghdfe_ols.h
 *
 * Backward-compatibility header.
 * All OLS utilities have been moved to the shared ctools_ols.h.
 * Part of the ctools Stata plugin suite.
 */

#ifndef CREGHDFE_OLS_H
#define CREGHDFE_OLS_H

#include "../ctools_ols.h"

/* Backward compatibility macros - use shared implementations */
#define cholesky ctools_cholesky
#define invert_from_cholesky ctools_invert_from_cholesky

#endif /* CREGHDFE_OLS_H */

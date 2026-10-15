/*
    civreghdfe_impl.h
    Instrumental Variables Regression with High-Dimensional Fixed Effects

    This module implements 2SLS/IV regression with HDFE absorption,
    reusing the creghdfe infrastructure for HDFE demeaning.
*/

#ifndef CIVREGHDFE_IMPL_H
#define CIVREGHDFE_IMPL_H

#include "../stplugin.h"
#include "../ctools_types.h"
#include "hdfe/hdfe_types.h"
#include "hdfe/hdfe_solver.h"
#include "../ctools_ols.h"
#include "hdfe/hdfe_remap.h"

/*
    Main entry point for civreghdfe plugin calls.

    Subcommands:
    - "iv_regression" - Full IV/2SLS regression with HDFE

    Returns STATA_OK on success, error code on failure.
*/
ST_retcode civreghdfe_main(const char *args);

/*
 * Cleanup function for civreghdfe persistent state.
 * Frees the global HDFE state if allocated.
 * Safe to call multiple times (idempotent).
 * Note: civreghdfe normally cleans up after itself, but this
 * handles cases where execution was interrupted.
 */
void civreghdfe_cleanup_state(void);

#endif /* CIVREGHDFE_IMPL_H */

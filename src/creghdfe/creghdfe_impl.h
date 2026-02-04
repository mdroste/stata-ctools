/*
 * creghdfe_impl.h
 *
 * Main entry point for creghdfe command (C-accelerated reghdfe)
 * Part of the ctools Stata plugin suite
 */

#ifndef CREGHDFE_IMPL_H
#define CREGHDFE_IMPL_H

#include "stplugin.h"
#include "creghdfe_types.h"

/*
 * Main entry point for creghdfe command.
 *
 * Dispatches to appropriate subcommand based on args:
 *   - "hdfe_init": Initialize HDFE factors, detect singletons
 *   - "mobility_groups": Compute connected components between FEs
 *   - "full_regression": Complete regression (init + partial out + OLS + VCE)
 *
 * Parameters:
 *   args - Command arguments string (subcommand name)
 *
 * Returns:
 *   0 on success, Stata error code on failure
 */
ST_retcode creghdfe_main(const char *args);

/*
 * Cleanup function for creghdfe persistent state.
 * Frees the global HDFE state if allocated.
 * Safe to call multiple times (idempotent).
 * Note: creghdfe normally cleans up after itself, but this
 * handles cases where execution was interrupted.
 */
void creghdfe_cleanup_state(void);

#endif /* CREGHDFE_IMPL_H */

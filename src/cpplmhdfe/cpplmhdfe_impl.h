/*
 * cpplmhdfe_impl.h
 *
 * Main entry point for cpplmhdfe command (C-accelerated ppmlhdfe)
 * Part of the ctools Stata plugin suite
 */

#ifndef CPPLMHDFE_IMPL_H
#define CPPLMHDFE_IMPL_H

#include "stplugin.h"
#include "cpplmhdfe_types.h"

/*
 * Main entry point for cpplmhdfe command.
 *
 * Dispatches to appropriate subcommand based on args:
 *   - "full_regression": Complete PPML regression (IRLS + VCE)
 *
 * Parameters:
 *   args - Command arguments string (subcommand name)
 *
 * Returns:
 *   0 on success, Stata error code on failure
 */
ST_retcode cpplmhdfe_main(const char *args);

/*
 * Cleanup function for cpplmhdfe persistent state.
 * Safe to call multiple times (idempotent).
 */
void cpplmhdfe_cleanup_state(void);

#endif /* CPPLMHDFE_IMPL_H */

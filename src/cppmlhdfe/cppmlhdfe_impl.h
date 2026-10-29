/*
 * cppmlhdfe_impl.h
 *
 * Main entry point for cppmlhdfe command (C-accelerated ppmlhdfe)
 * Part of the ctools Stata plugin suite
 */

#ifndef CPPMLHDFE_IMPL_H
#define CPPMLHDFE_IMPL_H

#include "stplugin.h"
#include "cppmlhdfe_types.h"

/*
 * Main entry point for cppmlhdfe command.
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
ST_retcode cppmlhdfe_main(const char *args);

/*
 * Cleanup function for cppmlhdfe persistent state.
 * Safe to call multiple times (idempotent).
 */
void cppmlhdfe_cleanup_state(void);

#endif /* CPPMLHDFE_IMPL_H */

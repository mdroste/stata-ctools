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

#endif /* CREGHDFE_IMPL_H */

/*
 * cqreg_impl.h
 *
 * Main entry point for C-accelerated quantile regression.
 * Part of the ctools suite.
 */

#ifndef CQREG_IMPL_H
#define CQREG_IMPL_H

#include "stplugin.h"

/*
 * Main entry point for cqreg command.
 * Called from ctools_plugin.c dispatcher.
 *
 * Arguments:
 *   args - Command-specific arguments as space-separated string
 *          Subcommands:
 *            "full_regression" - Complete QR with VCE
 *            "ipm_only"        - Just run IPM solver (for testing)
 *
 * Returns:
 *   0 on success, Stata error code on failure
 */
ST_retcode cqreg_main(const char *args);

#endif /* CQREG_IMPL_H */

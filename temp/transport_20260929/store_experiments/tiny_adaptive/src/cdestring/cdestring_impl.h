/*
 * cdestring_impl.h
 *
 * High-performance string-to-numeric conversion for Stata
 * Part of the ctools suite
 *
 * Replaces Stata's built-in destring command with a parallelized
 * C implementation for better performance on large datasets.
 *
 * Features:
 *   - Parallel conversion across observations using OpenMP
 *   - Fast double parsing using ctools_parse_double_fast
 *   - Support for all destring options: ignore(), force, percent, dpcomma, float
 *   - Handles multiple variables in a single call
 */

#ifndef CDESTRING_IMPL_H
#define CDESTRING_IMPL_H

#include "stplugin.h"

/*
 * Main entry point for cdestring command.
 *
 * Called from ctools_plugin.c dispatcher.
 *
 * Arguments format (space-separated):
 *   var_indices: space-separated 1-based variable indices (string vars to convert)
 *   gen_indices: space-separated 1-based variable indices (destination numeric vars)
 *   Options (keyword-based):
 *     ignore=<chars>  - Characters to strip before parsing
 *     force           - Also store variables that contain nonnumeric values
 *     percent         - Remove %; divide a variable by 100 if any value had %
 *     dpcomma         - First comma is the decimal point; "." is nonnumeric
 *     float           - Outputs are float (range check, percent rounding)
 *     nvars=<n>       - Number of variables being processed
 *
 * Sets local cdestring_failed to the per-variable counts of nonnumeric
 * values; without force such variables are not stored (the ado leaves them
 * untouched, as native destring does).
 *
 * @param args  Command arguments as space-separated string
 * @return      0 on success, Stata error code on failure
 */
ST_retcode cdestring_main(const char *args);

#endif /* CDESTRING_IMPL_H */

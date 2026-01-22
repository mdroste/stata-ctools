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
 *     force           - Convert non-numeric to missing instead of error
 *     percent         - Remove % and divide by 100
 *     dpcomma         - Use comma as decimal separator
 *     float           - Generate as float (not used in C, handled in ado)
 *     nvars=<n>       - Number of variables being processed
 *
 * @param args  Command arguments as space-separated string
 * @return      0 on success, Stata error code on failure
 */
ST_retcode cdestring_main(const char *args);

#endif /* CDESTRING_IMPL_H */

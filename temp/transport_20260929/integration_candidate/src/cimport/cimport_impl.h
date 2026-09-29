/*
 * cimport_impl.h
 *
 * High-performance CSV import for Stata
 * Part of the ctools suite
 *
 * Replaces "import delimited" with a C-accelerated implementation
 * featuring multi-threaded parsing, SIMD scanning, and direct DTA writing.
 */

#ifndef CIMPORT_IMPL_H
#define CIMPORT_IMPL_H

#include "stplugin.h"

/*
 * Main entry point for cimport command.
 *
 * Arguments passed via args string:
 *   "scan filename delimiter [noheader] [verbose]"
 *   "load filename delimiter [noheader] [verbose]"
 *   "writedta filename delimiter output=path [noheader] [verbose]"
 *
 * Returns:
 *   0 on success, Stata error code on failure
 */
ST_retcode cimport_main(const char *args);

/*
 * Cleanup function for cimport persistent state.
 * Frees the global parsed CSV context cache if loaded.
 * Safe to call multiple times (idempotent).
 */
void cimport_cleanup_cache(void);

#endif /* CIMPORT_IMPL_H */

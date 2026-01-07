/*
 * ctools_error.h - Unified error reporting utilities for ctools
 *
 * Provides consistent error message formatting across all ctools modules.
 * All functions output via Stata's SF_error() and SF_display().
 *
 * Part of the ctools suite for Stata.
 */

#ifndef CTOOLS_ERROR_H
#define CTOOLS_ERROR_H

#include <stdarg.h>

/*
 * Display an informational message to Stata.
 * Automatically adds newline if not present.
 *
 * Example: ctools_msg("csort", "Loaded %zu observations", nobs);
 */
void ctools_msg(const char *module, const char *fmt, ...);

/*
 * Display an error message to Stata.
 * Automatically adds newline if not present.
 *
 * Example: ctools_error("csort", "Failed to allocate memory");
 */
void ctools_error(const char *module, const char *fmt, ...);

/*
 * Display a standard "memory allocation failed" error.
 *
 * Example: ctools_error_alloc("csort");
 * Output:  "csort: memory allocation failed"
 */
void ctools_error_alloc(const char *module);

/*
 * Display a standard "out of memory" error with context.
 *
 * Example: ctools_error_alloc_ctx("csort", "thread arguments");
 * Output:  "csort: failed to allocate thread arguments"
 */
void ctools_error_alloc_ctx(const char *module, const char *context);

/*
 * Display a verbose/debug message (only if verbose flag is set).
 * Use for timing and progress information.
 *
 * Example: ctools_verbose("csort", verbose, "Sorted in %.1f ms", elapsed);
 */
void ctools_verbose(const char *module, int verbose, const char *fmt, ...);

/*
 * Convenience macros for common patterns
 */

/* Check allocation and return error code if NULL */
#define CTOOLS_CHECK_ALLOC(ptr, module, retval) \
    do { \
        if ((ptr) == NULL) { \
            ctools_error_alloc(module); \
            return (retval); \
        } \
    } while (0)

/* Check allocation with cleanup before return */
#define CTOOLS_CHECK_ALLOC_CLEANUP(ptr, module, cleanup, retval) \
    do { \
        if ((ptr) == NULL) { \
            ctools_error_alloc(module); \
            cleanup; \
            return (retval); \
        } \
    } while (0)

#endif /* CTOOLS_ERROR_H */

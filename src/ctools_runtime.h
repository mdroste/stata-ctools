/*
 * ctools_runtime.h - Runtime infrastructure for ctools
 *
 * Consolidates error reporting, timing, and lifecycle cleanup utilities.
 *
 * Part of the ctools suite for Stata.
 */

#ifndef CTOOLS_RUNTIME_H
#define CTOOLS_RUNTIME_H

#include <stdarg.h>
#include <stdint.h>
#include "stplugin.h"

/* ============================================================================
 * Error Reporting
 *
 * Provides consistent error message formatting across all ctools modules.
 * All functions output via Stata's SF_error() and SF_display().
 * ============================================================================ */

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

/* Convenience macros for common patterns */

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

/* ============================================================================
 * High-Resolution Timer
 *
 * Cross-platform timing utilities:
 * - macOS (mach_absolute_time)
 * - Linux (clock_gettime with CLOCK_MONOTONIC)
 * - Windows (QueryPerformanceCounter)
 * - Fallback (gettimeofday)
 * ============================================================================ */

/*
 * Initialize the timer subsystem.
 * Call once at program start. Safe to call multiple times.
 */
void ctools_timer_init(void);

/*
 * Get current time in seconds (double precision).
 * Resolution: ~1 microsecond on most platforms.
 * Returns: seconds since an arbitrary epoch (use for elapsed time only)
 */
double ctools_timer_seconds(void);

/*
 * Get current time in milliseconds (double precision).
 * Resolution: ~1 microsecond on most platforms.
 * Returns: milliseconds since an arbitrary epoch (use for elapsed time only)
 */
double ctools_timer_ms(void);

/*
 * Convenience macro for timing a code block.
 * Usage:
 *   double elapsed;
 *   CTOOLS_TIME_BLOCK(elapsed) {
 *       // code to time
 *   }
 *   printf("Elapsed: %.3f ms\n", elapsed);
 */
#define CTOOLS_TIME_BLOCK(elapsed_ms) \
    for (double _t_start = ctools_timer_ms(), _t_done = 0; \
         !_t_done; \
         elapsed_ms = ctools_timer_ms() - _t_start, _t_done = 1)

/*
 * Simple timer start/end macros for phase timing.
 * Usage:
 *   CTOOLS_TIMER_START(load);
 *   // ... loading code ...
 *   CTOOLS_TIMER_END(load);
 *   printf("Load took %.3f sec\n", _timer_load_elapsed);
 */
#define CTOOLS_TIMER_START(name) \
    double _timer_##name##_start = ctools_timer_seconds()

#define CTOOLS_TIMER_END(name) \
    double _timer_##name##_elapsed = ctools_timer_seconds() - _timer_##name##_start

/*
 * Timer macros that store elapsed time in a variable.
 * Usage:
 *   double t_load;
 *   CTOOLS_TIMER_BEGIN(load);
 *   // ... loading code ...
 *   CTOOLS_TIMER_STORE(load, t_load);
 */
#define CTOOLS_TIMER_BEGIN(name) \
    double _timer_##name##_begin = ctools_timer_seconds()

#define CTOOLS_TIMER_STORE(name, var) \
    (var) = ctools_timer_seconds() - _timer_##name##_begin

/* ============================================================================
 * Lifecycle Cleanup
 *
 * Centralized cleanup for persistent global state. Called automatically at
 * plugin exit to prevent memory leaks from interrupted operations.
 *
 * Design principles:
 * - All cleanup functions are idempotent (safe to call multiple times)
 * - Cleanup functions handle NULL/uninitialized state gracefully
 * - State required for multi-phase operations (cmerge, cimport) is preserved
 *   when calling the same command, but cleaned when switching commands
 * ============================================================================ */

/*
 * Command identifiers for tracking which command is running.
 * Used to preserve state for multi-phase operations within the same command.
 */
typedef enum {
    CTOOLS_CMD_NONE = 0,
    CTOOLS_CMD_CSORT,
    CTOOLS_CMD_CMERGE,
    CTOOLS_CMD_CIMPORT,
    CTOOLS_CMD_CEXPORT,
    CTOOLS_CMD_CREGHDFE,
    CTOOLS_CMD_CIVREGHDFE,
    CTOOLS_CMD_CQREG,
    CTOOLS_CMD_CBINSCATTER,
    CTOOLS_CMD_CENCODE,
    CTOOLS_CMD_CDECODE,
    CTOOLS_CMD_CWINSOR,
    CTOOLS_CMD_CDESTRING,
    CTOOLS_CMD_CSAMPLE,
    CTOOLS_CMD_CBSAMPLE,
    CTOOLS_CMD_CRANGESTAT,
    CTOOLS_CMD_CPSMATCH,
    CTOOLS_CMD_OTHER
} ctools_command_t;

/*
 * Set the current command being executed.
 * Call at the start of each command dispatch.
 */
void ctools_set_current_command(ctools_command_t cmd);

/*
 * Get the current command being executed.
 */
ctools_command_t ctools_get_current_command(void);

/*
 * Clean up stale state from other commands.
 * Called at the start of each plugin invocation to free state from
 * previously interrupted commands, while preserving state for the
 * current command's multi-phase operations.
 *
 * @param current_cmd  The command about to execute (its state is preserved)
 */
void ctools_cleanup_stale_state(ctools_command_t current_cmd);

/*
 * Force cleanup of all persistent state.
 * Called on explicit user request or when Stata is exiting.
 * Cleans up ALL state regardless of which command is running.
 */
void ctools_cleanup_all(void);

/*
 * Individual command cleanup functions.
 * These are called by ctools_cleanup_stale_state() and ctools_cleanup_all().
 * Each function is idempotent and handles uninitialized state gracefully.
 */

/* cmerge: Clears the using dataset cache (g_using_cache) */
void cmerge_cleanup_cache(void);

/* cimport: Clears the parsed CSV context cache (g_cimport_ctx) */
void cimport_cleanup_cache(void);

/* creghdfe: Clears the HDFE state (g_state) - usually already cleaned */
void creghdfe_cleanup_state(void);

/* civreghdfe: Clears the IV-HDFE state - usually already cleaned */
void civreghdfe_cleanup_state(void);

/* cexport: Clears any export context - usually stack-allocated */
void cexport_cleanup_state(void);

#endif /* CTOOLS_RUNTIME_H */

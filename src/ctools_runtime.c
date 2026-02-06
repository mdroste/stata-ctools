/*
 * ctools_runtime.c - Runtime infrastructure for ctools
 *
 * Consolidates error reporting, timing, and lifecycle cleanup.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "stplugin.h"
#include "ctools_runtime.h"

/* ============================================================================
 * Error Reporting
 * ============================================================================ */

/* Maximum message buffer size */
#define CTOOLS_MSG_BUFSIZE 1024

/*
 * Display an informational message to Stata.
 */
void ctools_msg(const char *module, const char *fmt, ...)
{
    char buf[CTOOLS_MSG_BUFSIZE];
    va_list args;
    int len;

    /* Format: "module: message\n" */
    if (module && module[0]) {
        len = snprintf(buf, sizeof(buf), "%s: ", module);
    } else {
        len = 0;
    }

    va_start(args, fmt);
    len += vsnprintf(buf + len, sizeof(buf) - len, fmt, args);
    va_end(args);

    /* Add newline if not present */
    if (len > 0 && len < (int)sizeof(buf) - 1 && buf[len - 1] != '\n') {
        buf[len] = '\n';
        buf[len + 1] = '\0';
    }

    SF_display(buf);
}

/*
 * Display an error message to Stata.
 */
void ctools_error(const char *module, const char *fmt, ...)
{
    char buf[CTOOLS_MSG_BUFSIZE];
    va_list args;
    int len;

    /* Format: "module: message\n" */
    if (module && module[0]) {
        len = snprintf(buf, sizeof(buf), "%s: ", module);
    } else {
        len = 0;
    }

    va_start(args, fmt);
    len += vsnprintf(buf + len, sizeof(buf) - len, fmt, args);
    va_end(args);

    /* Add newline if not present */
    if (len > 0 && len < (int)sizeof(buf) - 1 && buf[len - 1] != '\n') {
        buf[len] = '\n';
        buf[len + 1] = '\0';
    }

    SF_error(buf);
}

/*
 * Display a standard "memory allocation failed" error.
 */
void ctools_error_alloc(const char *module)
{
    ctools_error(module, "memory allocation failed");
}

/*
 * Display a standard "out of memory" error with context.
 */
void ctools_error_alloc_ctx(const char *module, const char *context)
{
    ctools_error(module, "failed to allocate %s", context);
}

/*
 * Display a verbose/debug message (only if verbose flag is set).
 */
void ctools_verbose(const char *module, int verbose, const char *fmt, ...)
{
    char buf[CTOOLS_MSG_BUFSIZE];
    va_list args;
    int len;

    if (!verbose) return;

    /* Format: "module: message\n" or just "message\n" for indented output */
    if (module && module[0]) {
        len = snprintf(buf, sizeof(buf), "%s: ", module);
    } else {
        len = 0;
    }

    va_start(args, fmt);
    len += vsnprintf(buf + len, sizeof(buf) - len, fmt, args);
    va_end(args);

    /* Add newline if not present */
    if (len > 0 && len < (int)sizeof(buf) - 1 && buf[len - 1] != '\n') {
        buf[len] = '\n';
        buf[len + 1] = '\0';
    }

    SF_display(buf);
}

/* ============================================================================
 * High-Resolution Timer
 * ============================================================================ */

#include "ctools_runtime.h"

/* Platform detection */
#if defined(__APPLE__) && defined(__MACH__)
    #define CTOOLS_TIMER_MACH 1
    #include <mach/mach_time.h>
#elif defined(_WIN32) || defined(_WIN64)
    #define CTOOLS_TIMER_WINDOWS 1
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
#elif defined(__linux__) || defined(__unix__) || defined(_POSIX_VERSION)
    #define CTOOLS_TIMER_POSIX 1
    #include <time.h>
#else
    #define CTOOLS_TIMER_FALLBACK 1
    #include <sys/time.h>
#endif

/* Platform-specific state - use simple volatile flags, no fancy initializers */
#if defined(CTOOLS_TIMER_MACH)
    static mach_timebase_info_data_t g_timebase_info;
    static volatile int g_timer_initialized = 0;

#elif defined(CTOOLS_TIMER_WINDOWS)
    static LARGE_INTEGER g_frequency;
    static volatile int g_timer_initialized = 0;

#else
    /* POSIX and fallback don't need initialization */
    static volatile int g_timer_initialized = 1;
#endif

/*
 * Initialize the timer subsystem.
 * Simple double-checked pattern - safe enough for our use case since
 * worst case is redundant initialization, not corruption.
 */
void ctools_timer_init(void) {
    if (g_timer_initialized) {
        return;
    }

#if defined(CTOOLS_TIMER_MACH)
    mach_timebase_info(&g_timebase_info);
    g_timer_initialized = 1;

#elif defined(CTOOLS_TIMER_WINDOWS)
    QueryPerformanceFrequency(&g_frequency);
    g_timer_initialized = 1;

#else
    /* Nothing to initialize */
    g_timer_initialized = 1;
#endif
}

/*
 * Get current time in seconds.
 */
double ctools_timer_seconds(void) {
    /* Auto-initialize on first call */
    if (!g_timer_initialized) {
        ctools_timer_init();
    }

#if defined(CTOOLS_TIMER_MACH)
    uint64_t t = mach_absolute_time();
    /* Convert to nanoseconds, then to seconds */
    double nanos = (double)t * g_timebase_info.numer / g_timebase_info.denom;
    return nanos / 1e9;

#elif defined(CTOOLS_TIMER_WINDOWS)
    LARGE_INTEGER counter;
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart / (double)g_frequency.QuadPart;

#elif defined(CTOOLS_TIMER_POSIX)
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;

#else /* Fallback: gettimeofday */
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
#endif
}

/*
 * Get current time in milliseconds.
 */
double ctools_timer_ms(void) {
    return ctools_timer_seconds() * 1000.0;
}

/* ============================================================================
 * Lifecycle Cleanup
 * ============================================================================ */

/* Track the current command for multi-phase operation support */
static ctools_command_t g_current_command = CTOOLS_CMD_NONE;
static ctools_command_t g_previous_command = CTOOLS_CMD_NONE;

void ctools_set_current_command(ctools_command_t cmd)
{
    g_previous_command = g_current_command;
    g_current_command = cmd;
}

ctools_command_t ctools_get_current_command(void)
{
    return g_current_command;
}

/*
 * Clean up stale state from other commands.
 *
 * Logic:
 * - If the current command is the same as the previous command, preserve state
 *   (supports multi-phase operations like cmerge load_using -> execute)
 * - If switching to a different command, clean up the previous command's state
 * - Special case: CTOOLS_CMD_NONE cleans everything (used at plugin exit)
 */
void ctools_cleanup_stale_state(ctools_command_t current_cmd)
{
    /* If same command as before, preserve state for multi-phase operations */
    if (current_cmd != CTOOLS_CMD_NONE && current_cmd == g_previous_command) {
        return;
    }

    /* Clean up state from other commands */

    /* cmerge cache - only clean if NOT running cmerge */
    if (current_cmd != CTOOLS_CMD_CMERGE) {
        cmerge_cleanup_cache();
    }

    /* cimport cache - only clean if NOT running cimport */
    if (current_cmd != CTOOLS_CMD_CIMPORT) {
        cimport_cleanup_cache();
    }

    /* creghdfe and civreghdfe clean themselves up, but call for safety */
    if (current_cmd != CTOOLS_CMD_CREGHDFE) {
        creghdfe_cleanup_state();
    }
    if (current_cmd != CTOOLS_CMD_CIVREGHDFE) {
        civreghdfe_cleanup_state();
    }

    /* cexport usually has no persistent state, but call for completeness */
    if (current_cmd != CTOOLS_CMD_CEXPORT) {
        cexport_cleanup_state();
    }
}

void ctools_cleanup_all(void)
{
    /* Force cleanup of everything */
    cmerge_cleanup_cache();
    cimport_cleanup_cache();
    creghdfe_cleanup_state();
    civreghdfe_cleanup_state();
    cexport_cleanup_state();

    /* Reset command tracking */
    g_current_command = CTOOLS_CMD_NONE;
    g_previous_command = CTOOLS_CMD_NONE;
}

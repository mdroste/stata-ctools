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

static const ctools_command_descriptor *active;
static int phase, next_phase;

void ctools_command_reset(void)
{
    if (active && active->cleanup) active->cleanup();
    active = NULL;
    phase = next_phase = 0;
}

int ctools_command_begin(const ctools_command_descriptor *cmd, const char *args)
{
    char op[32] = {0};
    int required = 0, preserve = 0;
    next_phase = 0;
    if (args) sscanf(args, "%31s", op);
    if (!cmd) { ctools_command_reset(); return 198; }
    switch (cmd->owner) {
    case CTOOLS_CMD_CMERGE:
        if (!strcmp(op, "load_using")) next_phase = 1;
        else if (!strcmp(op, "execute")) required = 1;
        break;
    case CTOOLS_CMD_CSPLIT:
        if (!strcmp(op, "scan")) next_phase = 1;
        else if (!strcmp(op, "write")) required = 1;
        break;
    case CTOOLS_CMD_CRANGEJOIN:
        if (!strcmp(op, "using")) next_phase = 1;
        else if (!strcmp(op, "prepare")) { required = 1; next_phase = 2; }
        else if (!strcmp(op, "write")) required = 2;
        break;
    case CTOOLS_CMD_CIMPORT:
        if (!strcmp(op, "scan")) next_phase = 1;
        else if (!strcmp(op, "blob")) { required = 1; next_phase = 1; }
        /* A direct load may parse its file without a preceding scan. */
        else if (!strcmp(op, "load")) preserve = 1;
        break;
    case CTOOLS_CMD_CIO:
        if (!strcmp(op, "scan")) next_phase = 1;
        else if (!strcmp(op, "column") || !strcmp(op, "label") ||
                 !strcmp(op, "blob") || !strcmp(op, "textnumbers") ||
                 !strcmp(op, "load")) { required = 1; next_phase = 1; }
        break;
    default: break;
    }
    int target_phase = next_phase;
    int same = active && active->owner == cmd->owner;
    if (required && (!same || phase != required)) {
        ctools_command_reset();
        return 198;
    }
    if (!same || (!required && !preserve)) ctools_command_reset();
    active = cmd;
    next_phase = target_phase;
    return 0;
}

void ctools_command_finish(int rc)
{
    if (rc || !next_phase) ctools_command_reset();
    else phase = next_phase;
}

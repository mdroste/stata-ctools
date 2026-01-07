/*
 * ctools_error.c - Unified error reporting utilities for ctools
 *
 * Provides consistent error message formatting across all ctools modules.
 */

#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "stplugin.h"
#include "ctools_error.h"

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

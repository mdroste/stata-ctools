/*
 * cexport_format.h
 * Formatting utilities for cexport
 *
 * Contains functions for:
 * - Double-to-string conversion (optimized)
 * - String quoting and escaping
 * - Row and header formatting
 */

#ifndef CEXPORT_FORMAT_H
#define CEXPORT_FORMAT_H

#include <stddef.h>
#include <stdbool.h>

#include "stplugin.h"
#include "cexport_context.h"

/* ========================================================================
   Row Formatting
   ======================================================================== */

/* Row formatters return the row length as int. A row longer than INT_MAX
   bytes returns this code instead; a larger buffer cannot help. */
#define CEXPORT_ROW_TOO_LONG (-3)

/*
    Format a single row of ALL NUMERIC variables into the output buffer.

    FAST PATH: No type checking, no string handling.
    Used when ctx->all_numeric is true for maximum throughput.

    @param ctx       Export context
    @param row_idx   Row index (0-based)
    @param buf       Output buffer
    @param buf_size  Size of output buffer
    @return          Number of bytes written, -1 on buffer overflow, or
                     CEXPORT_ROW_TOO_LONG
*/
int cexport_format_row_numeric(const cexport_context *ctx, size_t row_idx,
                               char *buf, size_t buf_size);

/*
    Format a single row into the output buffer.

    Handles both numeric and string variables.

    @param ctx       Export context
    @param row_idx   Row index (0-based)
    @param buf       Output buffer
    @param buf_size  Size of output buffer
    @return          Number of bytes written, -1 on buffer overflow, or
                     CEXPORT_ROW_TOO_LONG
*/
int cexport_format_row(const cexport_context *ctx, size_t row_idx,
                       char *buf, size_t buf_size);

/*
    Format one observation read directly from Stata (no loaded copy).
    Requires ctx->is_string; strL variables are not supported.

    @param ctx       Export context
    @param obs       1-based Stata observation
    @param buf       Output buffer
    @param buf_size  Size of output buffer
    @return          Bytes written, -1 when the buffer is too small (retry with
                     a larger buffer), -2 when Stata data cannot be read, or
                     CEXPORT_ROW_TOO_LONG
*/
int cexport_format_row_stream(const cexport_context *ctx, ST_int obs,
                              char *buf, size_t buf_size);

/*
    Format the header row into a buffer.

    @param ctx       Export context
    @param buf       Output buffer
    @param buf_size  Size of output buffer
    @return          Number of bytes written, or -1 on error
*/
int cexport_format_header(const cexport_context *ctx, char *buf, size_t buf_size);

/* ========================================================================
   Date and Time Display Formats
   ======================================================================== */

/*
    Compile a %t display format (daily, clock, weekly, monthly, quarterly,
    half-yearly or yearly). Returns NULL for a format this writer does not
    display itself, or when memory runs out; free() the result.
*/
cexport_timefmt *cexport_timefmt_compile(const char *format);

/* Bytes a cell in the format can take. */
size_t cexport_timefmt_reserve(const cexport_timefmt *format);

/* Bytes any numeric cell can take when written as a plain number. */
#define CEXPORT_NUMBER_RESERVE 48

#endif /* CEXPORT_FORMAT_H */

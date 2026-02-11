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

#include "cexport_context.h"

/* ========================================================================
   Row Formatting
   ======================================================================== */

/*
    Format a single row of ALL NUMERIC variables into the output buffer.

    FAST PATH: No type checking, no string handling.
    Used when ctx->all_numeric is true for maximum throughput.

    @param ctx       Export context
    @param row_idx   Row index (0-based)
    @param buf       Output buffer
    @param buf_size  Size of output buffer
    @return          Number of bytes written, or -1 on buffer overflow
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
    @return          Number of bytes written, or -1 on buffer overflow
*/
int cexport_format_row(const cexport_context *ctx, size_t row_idx,
                       char *buf, size_t buf_size);

/*
    Format the header row into a buffer.

    @param ctx       Export context
    @param buf       Output buffer
    @param buf_size  Size of output buffer
    @return          Number of bytes written, or -1 on error
*/
int cexport_format_header(const cexport_context *ctx, char *buf, size_t buf_size);

#endif /* CEXPORT_FORMAT_H */

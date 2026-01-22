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
   Double to String Conversion
   ======================================================================== */

/*
    Fast double to string - optimized version.

    Key optimizations:
    1. Fast path for integers (handles ~80% of typical data)
    2. Fast path for "simple decimals" with few decimal places
    3. Direct formatting without leading zero (Stata style)
    4. Single-pass trailing zero removal
    5. Fallback to snprintf only for complex cases (scientific notation)

    @param val            Value to convert
    @param buf            Output buffer
    @param buf_size       Size of output buffer
    @param missing_as_dot If true, Stata missing values are written as "."
    @param vtype          Variable storage type (for precision)
    @return               Number of characters written (not including null terminator)
*/
int cexport_double_to_str(double val, char *buf, int buf_size,
                          bool missing_as_dot, vartype_t vtype);

/* ========================================================================
   String Quoting and Escaping
   ======================================================================== */

/*
    Check if a string needs quoting.

    A string needs quoting if it contains:
    - The delimiter character
    - A double quote
    - A newline or carriage return

    @param str        String to check (NULL returns false)
    @param delimiter  Delimiter character
    @return           true if quoting is needed
*/
bool cexport_string_needs_quoting(const char *str, char delimiter);

/*
    Write a quoted string to a buffer.

    Doubles any internal quotes (CSV escaping convention).

    @param str       String to quote (NULL produces empty quoted string)
    @param buf       Output buffer
    @param buf_size  Size of output buffer
    @return          Number of characters written
*/
int cexport_write_quoted_string(const char *str, char *buf, int buf_size);

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

/*
 * cexport_format.c
 * Formatting utilities for cexport
 *
 * Contains optimized double-to-string conversion, string quoting,
 * and row/header formatting functions.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

#include "stplugin.h"
#include "ctools_types.h"
#include "cexport_format.h"
#include "cexport_context.h"

/* ========================================================================
   Configuration Constants
   ======================================================================== */

#define CEXPORT_DBL_BUFFER_SIZE  32  /* Buffer for double formatting */

/* ========================================================================
   Fast Double to String Conversion

   Custom implementation that avoids sprintf overhead.
   Handles integers, decimals, scientific notation, and Stata missing values.
   ======================================================================== */

/* Check if value is a Stata missing value */
static inline bool is_stata_missing(double val)
{
    return SF_is_missing(val);
}

/*
    Powers of 10 for fast decimal formatting.
    POW10[i] = 10^i for i in [0, 16]
*/
static const uint64_t POW10[17] = {
    1ULL,
    10ULL,
    100ULL,
    1000ULL,
    10000ULL,
    100000ULL,
    1000000ULL,
    10000000ULL,
    100000000ULL,
    1000000000ULL,
    10000000000ULL,
    100000000000ULL,
    1000000000000ULL,
    10000000000000ULL,
    100000000000000ULL,
    1000000000000000ULL,
    10000000000000000ULL
};

/*
    Format unsigned integer with leading zeros to exactly 'width' digits.
    Used for fractional part formatting.
    Returns number of characters written.
*/
static inline int format_frac_digits(uint64_t val, char *buf, int width)
{
    /* Format right-to-left */
    char temp[20];
    int pos = width;
    temp[pos] = '\0';

    while (pos > 0) {
        pos--;
        temp[pos] = '0' + (val % 10);
        val /= 10;
    }

    /* Copy and strip trailing zeros */
    int len = width;
    while (len > 1 && temp[len - 1] == '0') {
        len--;
    }

    memcpy(buf, temp, len);
    return len;
}

int cexport_double_to_str(double val, char *buf, int buf_size,
                          bool missing_as_dot, vartype_t vtype)
{
    (void)buf_size;  /* Used for safety in fallback path */

    /* Handle Stata missing values */
    if (is_stata_missing(val)) {
        if (missing_as_dot) {
            buf[0] = '.';
            buf[1] = '\0';
            return 1;
        } else {
            buf[0] = '\0';
            return 0;
        }
    }

    /* For float storage type, truncate to single precision before formatting.
     * This ensures we match Stata's native export delimited output exactly.
     * Without this, the float-to-double conversion introduces spurious precision
     * (e.g., 0.1f becomes 0.10000000149011612 as double). */
    if (vtype == VARTYPE_FLOAT) {
        val = (double)(float)val;
    }

    /* Handle special floating point values */
    if (isnan(val)) {
        buf[0] = '.';
        buf[1] = '\0';
        return 1;
    }
    if (isinf(val)) {
        if (val > 0) {
            memcpy(buf, "inf", 4);
            return 3;
        } else {
            memcpy(buf, "-inf", 5);
            return 4;
        }
    }

    /* FAST PATH 1: Exact integers
     * Handles ~80% of typical Stata data (IDs, counts, years, etc.)
     */
    if (val == (double)(int64_t)val && val >= -9007199254740992.0 && val <= 9007199254740992.0) {
        return ctools_int64_to_str((int64_t)val, buf);
    }

    /* Determine precision based on storage type */
    int max_frac_digits = (vtype == VARTYPE_FLOAT) ? 8 : 16;

    /* Handle sign */
    int pos = 0;
    double abs_val = val;
    if (val < 0) {
        buf[pos++] = '-';
        abs_val = -val;
    }

    /* FAST PATH 2: Simple decimals (most common in real data)
     * Values that can be represented exactly with limited decimal places.
     * Handles values like 1.5, 0.25, 123.456, etc.
     *
     * Upper bound: 1e9 to avoid precision issues with large floats.
     * Values >= 1e9 may not convert back exactly when reconstructed.
     */
    if (abs_val < 1e9 && abs_val > 1e-10) {
        /* Extract integer and fractional parts */
        double int_part_d = floor(abs_val);
        double frac_part = abs_val - int_part_d;

        /* Check if fractional part is "simple" (few decimal places needed) */
        if (frac_part > 0) {
            /* Try progressively more decimal places until we get exact match */
            for (int decimals = 1; decimals <= max_frac_digits && decimals <= 15; decimals++) {
                double scale = (double)POW10[decimals];
                double scaled_frac = frac_part * scale;
                uint64_t frac_int = (uint64_t)(scaled_frac + 0.5);

                /* Check for rounding to next integer */
                if (frac_int >= POW10[decimals]) {
                    int_part_d += 1.0;
                    frac_int = 0;
                    frac_part = 0;
                    break;
                }

                /* Check if this representation is exact enough */
                double reconstructed = int_part_d + (double)frac_int / scale;
                double rel_error = fabs(reconstructed - abs_val) / abs_val;

                if (rel_error < 1e-15 || (vtype == VARTYPE_FLOAT && rel_error < 1e-8)) {
                    /* Found exact representation! Format directly. */
                    uint64_t int_part = (uint64_t)int_part_d;

                    if (int_part == 0) {
                        /* Value between -1 and 1: Stata omits leading zero */
                        buf[pos++] = '.';
                    } else {
                        /* Format integer part */
                        pos += ctools_uint64_to_str_fast(int_part, buf + pos);
                        buf[pos++] = '.';
                    }

                    /* Format fractional part with trailing zeros stripped */
                    if (frac_int == 0) {
                        /* No fractional part after rounding - remove decimal point */
                        pos--;
                    } else {
                        pos += format_frac_digits(frac_int, buf + pos, decimals);
                    }

                    buf[pos] = '\0';
                    return pos;
                }
            }
        } else {
            /* No fractional part - should have been caught by integer path,
             * but handle it anyway */
            uint64_t int_part = (uint64_t)int_part_d;
            pos += ctools_uint64_to_str_fast(int_part, buf + pos);
            buf[pos] = '\0';
            return pos;
        }
    }

    /* FALLBACK: Use snprintf for complex cases
     * - Very small numbers (scientific notation)
     * - Very large numbers
     * - Numbers requiring many decimal places
     *
     * Stata uses scientific notation for very small numbers (< 0.0001)
     * with format like "1.00000000000e-10" (11 decimal places in mantissa)
     */
    char temp[48];
    int len;

    /* Check if this will use scientific notation (very small numbers) */
    if (abs_val > 0 && abs_val < 0.0001) {
        /* Use Stata-compatible scientific notation format
         * Stata uses 11 decimal places for doubles, 7 for floats */
        if (vtype == VARTYPE_FLOAT) {
            len = snprintf(temp, sizeof(temp), "%.7e", val);
        } else {
            len = snprintf(temp, sizeof(temp), "%.11e", val);
        }
        /* Copy directly - no post-processing needed for scientific notation */
        memcpy(buf + pos, temp + (temp[0] == '-' ? 1 : 0), len - (temp[0] == '-' ? 1 : 0));
        pos += len - (temp[0] == '-' ? 1 : 0);
        buf[pos] = '\0';
        return pos;
    }

    /* For other cases, use %g format */
    if (vtype == VARTYPE_FLOAT) {
        len = snprintf(temp, sizeof(temp), "%.8g", val);
    } else {
        len = snprintf(temp, sizeof(temp), "%.16g", val);
    }

    /* Process snprintf output in single pass:
     * - Skip leading zero for values between -1 and 1
     * - Strip trailing zeros after decimal point
     */
    int src = 0;
    int has_sign = (temp[0] == '-');
    if (has_sign) {
        src++;
    }

    /* Skip leading zero: "0." -> "." */
    int skip_zero = (temp[src] == '0' && temp[src + 1] == '.');
    if (skip_zero) {
        src++;
    }

    /* Find decimal point and exponent positions */
    int dot_pos = -1;
    int exp_pos = -1;
    for (int i = src; i < len; i++) {
        if (temp[i] == '.') dot_pos = i;
        else if (temp[i] == 'e' || temp[i] == 'E') {
            exp_pos = i;
            break;
        }
    }

    /* Copy to output, stripping trailing zeros */
    pos = 0;
    if (has_sign) {
        buf[pos++] = '-';
    }

    if (dot_pos >= 0 && exp_pos < 0) {
        /* Has decimal, no exponent - can strip trailing zeros */
        int last_nonzero = len - 1;
        while (last_nonzero > dot_pos && temp[last_nonzero] == '0') {
            last_nonzero--;
        }
        /* If only decimal point left, remove it too */
        if (last_nonzero == dot_pos) {
            last_nonzero--;
        }
        for (int i = src; i <= last_nonzero; i++) {
            buf[pos++] = temp[i];
        }
    } else if (dot_pos >= 0 && exp_pos > 0) {
        /* Has decimal and exponent - copy as-is (Stata keeps trailing zeros in sci notation) */
        for (int i = src; i < len; i++) {
            buf[pos++] = temp[i];
        }
    } else {
        /* No decimal point - copy as-is */
        for (int i = src; i < len; i++) {
            buf[pos++] = temp[i];
        }
    }

    buf[pos] = '\0';
    return pos;
}

/* ========================================================================
   String Quoting and Escaping
   ======================================================================== */

bool cexport_string_needs_quoting(const char *str, char delimiter)
{
    if (str == NULL) return false;

    for (const char *p = str; *p != '\0'; p++) {
        if (*p == delimiter || *p == '"' || *p == '\n' || *p == '\r') {
            return true;
        }
    }
    return false;
}

int cexport_write_quoted_string(const char *str, char *buf, int buf_size)
{
    /* Handle edge cases: need at least 3 bytes for empty quoted string '""' + null */
    if (buf_size < 1) {
        return 0;
    }
    if (buf_size < 3) {
        buf[0] = '\0';
        return 0;
    }
    if (str == NULL) {
        buf[0] = '"';
        buf[1] = '"';
        buf[2] = '\0';
        return 2;
    }

    int pos = 0;
    buf[pos++] = '"';

    for (const char *p = str; *p != '\0' && pos < buf_size - 2; p++) {
        if (*p == '"') {
            /* Escape quote by doubling it */
            if (pos < buf_size - 3) {
                buf[pos++] = '"';
                buf[pos++] = '"';
            }
        } else {
            buf[pos++] = *p;
        }
    }

    buf[pos++] = '"';
    buf[pos] = '\0';
    return pos;
}

/* ========================================================================
   Row Formatting
   ======================================================================== */

int cexport_format_row_numeric(const cexport_context *ctx, size_t row_idx,
                               char *buf, size_t buf_size)
{
    size_t pos = 0;
    size_t nvars = ctx->data.nvars;
    char delimiter = ctx->delimiter;
    const char *line_end = ctx->line_ending;
    int line_end_len = ctx->line_ending_len;

    for (size_t j = 0; j < nvars; j++) {
        double val = ctx->data.vars[j].data.dbl[row_idx];
        vartype_t vtype = ctx->vartypes[j];

        /* Check buffer space (assume max 32 chars for number + delimiter/newline) */
        if (pos + 34 + (size_t)line_end_len > buf_size) {
            return -1;
        }

        /* Write directly to output buffer */
        int field_len = cexport_double_to_str(val, buf + pos, 32, false, vtype);
        pos += field_len;

        /* Add delimiter or line ending */
        if (j < nvars - 1) {
            buf[pos++] = delimiter;
        } else {
            /* Copy line ending (1 or 2 bytes) */
            memcpy(buf + pos, line_end, line_end_len);
            pos += line_end_len;
        }
    }

    return (int)pos;
}

int cexport_format_row(const cexport_context *ctx, size_t row_idx,
                       char *buf, size_t buf_size)
{
    char field_buf[256];  /* Small buffer - most fields are <50 chars */
    size_t pos = 0;
    size_t nvars = ctx->data.nvars;
    char delimiter = ctx->delimiter;
    const char *line_end = ctx->line_ending;
    int line_end_len = ctx->line_ending_len;

    for (size_t j = 0; j < nvars; j++) {
        const stata_variable *var = &ctx->data.vars[j];
        int field_len;

        if (var->type == STATA_TYPE_DOUBLE) {
            /* FAST PATH: Numeric variable - inline for speed */
            double val = var->data.dbl[row_idx];
            vartype_t vtype = ctx->vartypes[j];

            /* Check buffer space (assume max 32 chars for number + line ending) */
            if (pos + 34 + (size_t)line_end_len > buf_size) {
                return -1;
            }

            /* Write directly to output buffer */
            field_len = cexport_double_to_str(val, buf + pos, sizeof(field_buf), false, vtype);
            pos += field_len;
        } else {
            /* String variable */
            const char *str = var->data.str[row_idx];
            if (str == NULL) str = "";

            bool need_quote = ctx->quote_strings ||
                             (ctx->quote_if_needed && cexport_string_needs_quoting(str, delimiter));

            if (need_quote) {
                field_len = cexport_write_quoted_string(str, field_buf, sizeof(field_buf));
                if (pos + field_len + 2 + (size_t)line_end_len > buf_size) return -1;
                memcpy(buf + pos, field_buf, field_len);
            } else {
                field_len = (int)strlen(str);
                if (pos + field_len + 2 + (size_t)line_end_len > buf_size) return -1;
                memcpy(buf + pos, str, field_len);
            }
            pos += field_len;
        }

        /* Add delimiter or line ending */
        if (j < nvars - 1) {
            buf[pos++] = delimiter;
        } else {
            /* Copy line ending (1 or 2 bytes) */
            memcpy(buf + pos, line_end, line_end_len);
            pos += line_end_len;
        }
    }

    return (int)pos;
}

int cexport_format_header(const cexport_context *ctx, char *buf, size_t buf_size)
{
    size_t pos = 0;
    const char *line_end = ctx->line_ending;
    int line_end_len = ctx->line_ending_len;

    for (size_t j = 0; j < ctx->nvars; j++) {
        const char *name = ctx->varnames[j];
        size_t name_len = strlen(name);

        /* Check if name needs quoting */
        bool need_quote = ctx->quote_strings ||
                         (ctx->quote_if_needed && cexport_string_needs_quoting(name, ctx->delimiter));

        if (need_quote) {
            if (pos + name_len * 2 + 4 + (size_t)line_end_len > buf_size) {
                return -1;
            }
            int quoted_len = cexport_write_quoted_string(name, buf + pos, buf_size - pos);
            pos += quoted_len;
        } else {
            if (pos + name_len + 2 + (size_t)line_end_len > buf_size) {
                return -1;
            }
            memcpy(buf + pos, name, name_len);
            pos += name_len;
        }

        if (j < ctx->nvars - 1) {
            buf[pos++] = ctx->delimiter;
        } else {
            /* Copy line ending (1 or 2 bytes) */
            memcpy(buf + pos, line_end, line_end_len);
            pos += line_end_len;
        }
    }

    return (int)pos;
}

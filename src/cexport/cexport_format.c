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
#include "ctools_simd.h"
#include "cexport_format.h"
#include "cexport_context.h"

/* ========================================================================
   Configuration Constants
   ======================================================================== */

#define CEXPORT_DBL_BUFFER_SIZE  48  /* Buffer for double formatting */

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

static int cexport_double_to_str(double val, char *buf, int buf_size,
                                 bool missing_as_dot, vartype_t vtype)
{
    /* buf_size is used for bounds checking in snprintf fallback */

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

    /* Handle sign */
    int pos = 0;
    double abs_val = val;
    if (val < 0) {
        buf[pos++] = '-';
        abs_val = -val;
    }

    /* Stata uses scientific notation for very small numbers (< 0.0001)
     * with format like "1.00000000000e-10" (11 decimal places in mantissa) */
    if (abs_val > 0 && abs_val < 0.0001) {
        char temp[48];
        int len;
        if (vtype == VARTYPE_FLOAT) {
            len = snprintf(temp, sizeof(temp), "%.7e", val);
        } else {
            len = snprintf(temp, sizeof(temp), "%.11e", val);
        }
        memcpy(buf + pos, temp + (temp[0] == '-' ? 1 : 0), len - (temp[0] == '-' ? 1 : 0));
        pos += len - (temp[0] == '-' ? 1 : 0);
        buf[pos] = '\0';
        return pos;
    }

    /* FAST PATH: Custom %g-equivalent for normal decimals.
     * Avoids snprintf overhead (~150ns per call) by using integer arithmetic
     * with the pow10 lookup table. Handles values in [0.0001, 1e16). */
    int sig_digits = (vtype == VARTYPE_FLOAT) ? 8 : 16;

    /* Find decimal exponent: abs_val is in [10^exp10, 10^(exp10+1)) */
    int exp10;
    if (abs_val >= 10.0) {
        if (abs_val >= 1e8) {
            if (abs_val >= 1e12) {
                if (abs_val >= 1e14) { exp10 = (abs_val >= 1e15) ? 15 : 14; }
                else { exp10 = (abs_val >= 1e13) ? 13 : 12; }
            } else {
                if (abs_val >= 1e10) { exp10 = (abs_val >= 1e11) ? 11 : 10; }
                else { exp10 = (abs_val >= 1e9) ? 9 : 8; }
            }
        } else {
            if (abs_val >= 1e4) {
                if (abs_val >= 1e6) { exp10 = (abs_val >= 1e7) ? 7 : 6; }
                else { exp10 = (abs_val >= 1e5) ? 5 : 4; }
            } else {
                if (abs_val >= 1e2) { exp10 = (abs_val >= 1e3) ? 3 : 2; }
                else { exp10 = 1; }
            }
        }
    } else if (abs_val >= 1.0) {
        exp10 = 0;
    } else if (abs_val >= 0.1) {
        exp10 = -1;
    } else if (abs_val >= 0.01) {
        exp10 = -2;
    } else if (abs_val >= 0.001) {
        exp10 = -3;
    } else {
        exp10 = -4;  /* abs_val >= 0.0001 guaranteed by check above */
    }

    /* If %g would use scientific notation (exponent >= sig_digits), use snprintf */
    if (exp10 >= sig_digits) {
        goto snprintf_fallback;
    }

    {
        /* Scale to extract sig_digits significant digits as an integer */
        int scale = sig_digits - 1 - exp10;
        if (scale < 0 || scale > 22) goto snprintf_fallback;

        uint64_t sig = (uint64_t)llround(abs_val * ctools_pow10_table[scale]);

        /* Handle rounding overflow (e.g., 9.999...9 rounds up) */
        uint64_t max_sig = (sig_digits == 8) ? 100000000ULL : 10000000000000000ULL;
        if (sig >= max_sig) {
            exp10++;
            if (exp10 >= sig_digits) goto snprintf_fallback;
            scale = sig_digits - 1 - exp10;
            if (scale < 0 || scale > 22) goto snprintf_fallback;
            sig = (uint64_t)llround(abs_val * ctools_pow10_table[scale]);
        }

        /* Convert significant digits to string */
        char digit_buf[24];
        int digit_len = ctools_uint64_to_str(sig, digit_buf);

        /* Pad with leading zeros if digit_len < sig_digits
         * (can happen when leading significant digits round to fewer digits) */
        while (digit_len < sig_digits) {
            memmove(digit_buf + 1, digit_buf, digit_len + 1);
            digit_buf[0] = '0';
            digit_len++;
        }

        if (exp10 >= 0) {
            /* Value >= 1.0: integer part is first (exp10+1) digits */
            int int_part_len = exp10 + 1;
            memcpy(buf + pos, digit_buf, int_part_len);
            pos += int_part_len;

            /* Find last non-zero fractional digit */
            int frac_end = digit_len;
            while (frac_end > int_part_len && digit_buf[frac_end - 1] == '0') {
                frac_end--;
            }

            if (frac_end > int_part_len) {
                buf[pos++] = '.';
                memcpy(buf + pos, digit_buf + int_part_len, frac_end - int_part_len);
                pos += frac_end - int_part_len;
            }
        } else {
            /* Value < 1.0: Stata strips leading zero → ".xxxx" */
            buf[pos++] = '.';

            /* Leading zeros after decimal point */
            int leading_zeros = -exp10 - 1;
            for (int z = 0; z < leading_zeros; z++) {
                buf[pos++] = '0';
            }

            /* Strip trailing zeros from significant digits */
            int frac_end = digit_len;
            while (frac_end > 0 && digit_buf[frac_end - 1] == '0') {
                frac_end--;
            }

            memcpy(buf + pos, digit_buf, frac_end);
            pos += frac_end;
        }

        buf[pos] = '\0';
        return pos;
    }

snprintf_fallback:
    {
        /* Rare fallback for very large numbers needing scientific notation */
        char temp[48];
        int len;
        if (vtype == VARTYPE_FLOAT) {
            len = snprintf(temp, sizeof(temp), "%.8g", val);
        } else {
            len = snprintf(temp, sizeof(temp), "%.16g", val);
        }

        int src = 0;
        int has_sign = (temp[0] == '-');
        if (has_sign) src++;

        /* Skip leading zero: "0." -> "." */
        if (temp[src] == '0' && temp[src + 1] == '.') src++;

        /* Find decimal point and exponent positions */
        int dot_pos = -1;
        int exp_found = -1;
        for (int i = src; i < len; i++) {
            if (temp[i] == '.') dot_pos = i;
            else if (temp[i] == 'e' || temp[i] == 'E') { exp_found = i; break; }
        }

        pos = 0;
        if (has_sign) buf[pos++] = '-';

        if (dot_pos >= 0 && exp_found < 0) {
            int last_nonzero = len - 1;
            while (last_nonzero > dot_pos && temp[last_nonzero] == '0') last_nonzero--;
            if (last_nonzero == dot_pos) last_nonzero--;
            for (int i = src; i <= last_nonzero; i++) buf[pos++] = temp[i];
        } else {
            for (int i = src; i < len; i++) buf[pos++] = temp[i];
        }

        /* Safety: ensure we don't exceed buf_size */
        if (pos >= buf_size) pos = buf_size - 1;
        buf[pos] = '\0';
        return pos;
    }
}

/* ========================================================================
   String Quoting and Escaping
   ======================================================================== */

static bool cexport_string_needs_quoting(const char *str, char delimiter)
{
    if (str == NULL) return false;

    size_t len = strlen(str);
    if (len == 0) return false;

    const char *p = str;

#if CTOOLS_HAS_AVX2
    {
        const __m256i v_delim = _mm256_set1_epi8(delimiter);
        const __m256i v_quote = _mm256_set1_epi8('"');
        const __m256i v_nl    = _mm256_set1_epi8('\n');
        const __m256i v_cr    = _mm256_set1_epi8('\r');

        while (len >= 32) {
            __m256i chunk = _mm256_loadu_si256((const __m256i *)p);
            __m256i cmp = _mm256_or_si256(
                _mm256_or_si256(_mm256_cmpeq_epi8(chunk, v_delim), _mm256_cmpeq_epi8(chunk, v_quote)),
                _mm256_or_si256(_mm256_cmpeq_epi8(chunk, v_nl), _mm256_cmpeq_epi8(chunk, v_cr)));
            if (_mm256_movemask_epi8(cmp)) return true;
            p += 32;
            len -= 32;
        }
    }
    /* SSE2 tail for 16-31 bytes */
    if (len >= 16) {
        const __m128i v_delim = _mm_set1_epi8(delimiter);
        const __m128i v_quote = _mm_set1_epi8('"');
        const __m128i v_nl    = _mm_set1_epi8('\n');
        const __m128i v_cr    = _mm_set1_epi8('\r');

        __m128i chunk = _mm_loadu_si128((const __m128i *)p);
        __m128i cmp = _mm_or_si128(
            _mm_or_si128(_mm_cmpeq_epi8(chunk, v_delim), _mm_cmpeq_epi8(chunk, v_quote)),
            _mm_or_si128(_mm_cmpeq_epi8(chunk, v_nl), _mm_cmpeq_epi8(chunk, v_cr)));
        if (_mm_movemask_epi8(cmp)) return true;
        p += 16;
        len -= 16;
    }
#elif CTOOLS_HAS_SSE2
    {
        const __m128i v_delim = _mm_set1_epi8(delimiter);
        const __m128i v_quote = _mm_set1_epi8('"');
        const __m128i v_nl    = _mm_set1_epi8('\n');
        const __m128i v_cr    = _mm_set1_epi8('\r');

        while (len >= 16) {
            __m128i chunk = _mm_loadu_si128((const __m128i *)p);
            __m128i cmp = _mm_or_si128(
                _mm_or_si128(_mm_cmpeq_epi8(chunk, v_delim), _mm_cmpeq_epi8(chunk, v_quote)),
                _mm_or_si128(_mm_cmpeq_epi8(chunk, v_nl), _mm_cmpeq_epi8(chunk, v_cr)));
            if (_mm_movemask_epi8(cmp)) return true;
            p += 16;
            len -= 16;
        }
    }
#elif CTOOLS_HAS_NEON
    {
        const uint8x16_t v_delim = vdupq_n_u8((uint8_t)delimiter);
        const uint8x16_t v_quote = vdupq_n_u8('"');
        const uint8x16_t v_nl    = vdupq_n_u8('\n');
        const uint8x16_t v_cr    = vdupq_n_u8('\r');

        while (len >= 16) {
            uint8x16_t chunk = vld1q_u8((const uint8_t *)p);
            uint8x16_t cmp = vorrq_u8(
                vorrq_u8(vceqq_u8(chunk, v_delim), vceqq_u8(chunk, v_quote)),
                vorrq_u8(vceqq_u8(chunk, v_nl), vceqq_u8(chunk, v_cr)));
            uint64x2_t cmp64 = vreinterpretq_u64_u8(cmp);
            if (vgetq_lane_u64(cmp64, 0) || vgetq_lane_u64(cmp64, 1)) return true;
            p += 16;
            len -= 16;
        }
    }
#endif

    /* Scalar tail */
    while (len > 0) {
        char c = *p;
        if (c == delimiter || c == '"' || c == '\n' || c == '\r') return true;
        p++;
        len--;
    }
    return false;
}

static int cexport_write_quoted_string(const char *str, char *buf, int buf_size)
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
    size_t nvars = ctx->filtered.data.nvars;
    char delimiter = ctx->delimiter;
    const char *line_end = ctx->line_ending;
    int line_end_len = ctx->line_ending_len;

    for (size_t j = 0; j < nvars; j++) {
        double val = ctx->filtered.data.vars[j].data.dbl[row_idx];
        vartype_t vtype = ctx->vartypes[j];

        /* Check buffer space for number + delimiter/newline */
        if (pos + CEXPORT_DBL_BUFFER_SIZE + 2 + (size_t)line_end_len > buf_size) {
            return -1;
        }

        /* Write directly to output buffer */
        int field_len = cexport_double_to_str(val, buf + pos, CEXPORT_DBL_BUFFER_SIZE, false, vtype);
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
    size_t nvars = ctx->filtered.data.nvars;
    char delimiter = ctx->delimiter;
    const char *line_end = ctx->line_ending;
    int line_end_len = ctx->line_ending_len;

    for (size_t j = 0; j < nvars; j++) {
        const stata_variable *var = &ctx->filtered.data.vars[j];
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

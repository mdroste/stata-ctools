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

/* 5^k and 10^k for the exact decimal scaling below */
static const uint64_t cexport_pow5[28] = {
    1ULL, 5ULL, 25ULL, 125ULL, 625ULL, 3125ULL, 15625ULL, 78125ULL, 390625ULL,
    1953125ULL, 9765625ULL, 48828125ULL, 244140625ULL, 1220703125ULL,
    6103515625ULL, 30517578125ULL, 152587890625ULL, 762939453125ULL,
    3814697265625ULL, 19073486328125ULL, 95367431640625ULL, 476837158203125ULL,
    2384185791015625ULL, 11920928955078125ULL, 59604644775390625ULL,
    298023223876953125ULL, 1490116119384765625ULL, 7450580596923828125ULL
};
static const uint64_t cexport_pow10[20] = {
    1ULL, 10ULL, 100ULL, 1000ULL, 10000ULL, 100000ULL, 1000000ULL, 10000000ULL,
    100000000ULL, 1000000000ULL, 10000000000ULL, 100000000000ULL,
    1000000000000ULL, 10000000000000ULL, 100000000000000ULL,
    1000000000000000ULL, 10000000000000000ULL, 100000000000000000ULL,
    1000000000000000000ULL, 10000000000000000000ULL
};

/* Decimal digits in v > 0 */
static inline int cexport_digit_count(uint64_t v)
{
    int t = ((64 - __builtin_clzll(v)) * 1233) >> 12;
    return t + (v >= cexport_pow10[t]);
}

/* Write v right to left ending before end, exactly width digits (zero padded
 * when width > 0, otherwise all digits of v > 0); returns the first digit. */
static inline char *cexport_write_digits(uint64_t v, char *end, int width)
{
    char *p = end;
    if (width > 0) {
        for (; width >= 2; width -= 2) {
            const char *pair = CTOOLS_DIGIT_PAIRS + (v % 100) * 2;
            v /= 100;
            *--p = pair[1];
            *--p = pair[0];
        }
        if (width) *--p = (char)('0' + v % 10);
        return p;
    }
    while (v >= 100) {
        const char *pair = CTOOLS_DIGIT_PAIRS + (v % 100) * 2;
        v /= 100;
        *--p = pair[1];
        *--p = pair[0];
    }
    if (v >= 10) {
        *--p = CTOOLS_DIGIT_PAIRS[v * 2 + 1];
        *--p = CTOOLS_DIGIT_PAIRS[v * 2];
    } else {
        *--p = (char)('0' + v);
    }
    return p;
}

/* All Stata missing codes compare >= SV_missval, so only those values (and
 * NaN) need Stata's classification callback. */
static inline bool is_stata_missing(double val)
{
    return !(val < SV_missval) && SF_is_missing(val);
}

/* Scale a positive binary64 by 10^places and round half away from zero.
 * Multiplying in binary64 first can change the last exported decimal digit.
 * The exact product fits in 128 bits for the maximum 23 decimal places
 * used in the 19-significant-digit intermediate below. */
static uint64_t cexport_decimal_integer(double value, int places)
{
    uint64_t bits;
    memcpy(&bits, &value, sizeof(bits));
    uint64_t mantissa = (bits & UINT64_C(0x000fffffffffffff)) | (UINT64_C(1) << 52);
    int shift = (int)((bits >> 52) & 0x7ff) - 1023 - 52 + places;
    __uint128_t product = mantissa;
    if (places >= 0 && places < 28) product *= cexport_pow5[places];
    else for (int i = 0; i < places; i++) product *= 5;
    if (shift >= 0) return (uint64_t)(product << shift);
    int right = -shift;
    return (uint64_t)((product + ((__uint128_t)1 << (right - 1))) >> right);
}

static int cexport_double_to_str(double val, char *buf, int buf_size,
                                 bool missing_as_dot, vartype_t vtype)
{
    /* Handle Stata missing values */
    if (is_stata_missing(val)) {
        /* Extended missings are values, not empty cells. */
        uint64_t bits, base;
        double missing = SV_missval;
        memcpy(&bits, &val, sizeof(bits));
        memcpy(&base, &missing, sizeof(base));
        uint64_t code = (bits - base) >> 40;
        if (bits > base && code >= 1 && code <= 26) {
            buf[0] = '.';
            buf[1] = (char)('a' + code - 1);
            buf[2] = '\0';
            return 2;
        }
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
    if (val > -1e16 && val < 1e16 && val == (double)(int64_t)val) {
        return ctools_int64_to_str((int64_t)val, buf);
    }

    /* Stata's default general format uses 16 significant digits for doubles
     * and 8 for floats, choosing fixed notation when it fits the field. The
     * scientific mantissa has fixed precision, reduced for 3-digit exponents. */
    const int digits = vtype == VARTYPE_FLOAT ? 8 : 16;
    const int width = vtype == VARTYPE_FLOAT ? 14 : 18;
    double magnitude = fabs(val);
    int exponent = (int)floor(log10(magnitude));
    if (magnitude >= 1e-5 && exponent < digits) {
        int decimals = digits - 1 - exponent;
        if (decimals > 16) decimals = 16;
        if (decimals < 0) decimals = 0;
        /* Stata rounds decimal halfway cases away from zero. printf's
         * ties-to-even conversion differs for values such as 2384959.25f. */
        /* Retain a 19-digit decimal intermediate before display rounding.
         * Rounding the binary value directly gives different results near
         * decimal halfway cases (covered by the native comparison suite). */
        int precise_places = 18 - exponent;
        uint64_t precise = cexport_decimal_integer(magnitude, precise_places);
        int drop = precise_places - decimals;
        uint64_t divisor = 1;
        if (drop >= 0 && drop < 20) divisor = cexport_pow10[drop];
        else for (int i = decimals; i < precise_places; i++) divisor *= 10;
        uint64_t scaled = (precise + divisor / 2) / divisor;
        /* Stata omits trailing fractional zeros and the leading integer zero
         * (".05"). Strip zeros arithmetically, then write in place. */
        int dec = decimals;
        if (scaled == 0) {
            int len = (val < 0) + (dec == 0);
            if (len < buf_size) {
                char *p = buf;
                if (val < 0) *p++ = '-';
                if (dec == 0) *p++ = '0';
                *p = '\0';
                return len;
            }
        } else {
            while (dec >= 8 && scaled % 100000000 == 0) { scaled /= 100000000; dec -= 8; }
            if (dec >= 4 && scaled % 10000 == 0) { scaled /= 10000; dec -= 4; }
            if (dec >= 2 && scaled % 100 == 0) { scaled /= 100; dec -= 2; }
            if (dec >= 1 && scaled % 10 == 0) { scaled /= 10; dec -= 1; }
            uint64_t whole = scaled / cexport_pow10[dec];
            uint64_t frac = scaled - whole * cexport_pow10[dec];
            int len = (val < 0) + (whole ? cexport_digit_count(whole) : 0) + (dec ? dec + 1 : 0);
            if (len < buf_size) {
                char *p = buf + len;
                *p = '\0';
                if (dec) {
                    p = cexport_write_digits(frac, p, dec);
                    *--p = '.';
                }
                if (whole) p = cexport_write_digits(whole, p, 0);
                if (val < 0) *--p = '-';
                return len;
            }
        }
    }
    int decimals = width - 7 - (abs(exponent) >= 100);
    int len = snprintf(buf, (size_t)buf_size, "%.*e", decimals, val);
    return len >= buf_size ? -1 : len;
}

/* ========================================================================
   String Quoting and Escaping
   ======================================================================== */

/* Native export quotes a field only when it contains the delimiter or a
 * double quote; embedded line breaks and padding are written as they are. */
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

        while (len >= 32) {
            __m256i chunk = _mm256_loadu_si256((const __m256i *)p);
            __m256i cmp = _mm256_or_si256(_mm256_cmpeq_epi8(chunk, v_delim), _mm256_cmpeq_epi8(chunk, v_quote));
            if (_mm256_movemask_epi8(cmp)) return true;
            p += 32;
            len -= 32;
        }
    }
    /* SSE2 tail for 16-31 bytes */
    if (len >= 16) {
        const __m128i v_delim = _mm_set1_epi8(delimiter);
        const __m128i v_quote = _mm_set1_epi8('"');

        __m128i chunk = _mm_loadu_si128((const __m128i *)p);
        __m128i cmp = _mm_or_si128(_mm_cmpeq_epi8(chunk, v_delim), _mm_cmpeq_epi8(chunk, v_quote));
        if (_mm_movemask_epi8(cmp)) return true;
        p += 16;
        len -= 16;
    }
#elif CTOOLS_HAS_SSE2
    {
        const __m128i v_delim = _mm_set1_epi8(delimiter);
        const __m128i v_quote = _mm_set1_epi8('"');

        while (len >= 16) {
            __m128i chunk = _mm_loadu_si128((const __m128i *)p);
            __m128i cmp = _mm_or_si128(_mm_cmpeq_epi8(chunk, v_delim), _mm_cmpeq_epi8(chunk, v_quote));
            if (_mm_movemask_epi8(cmp)) return true;
            p += 16;
            len -= 16;
        }
    }
#elif CTOOLS_HAS_NEON
    {
        const uint8x16_t v_delim = vdupq_n_u8((uint8_t)delimiter);
        const uint8x16_t v_quote = vdupq_n_u8('"');

        while (len >= 16) {
            uint8x16_t chunk = vld1q_u8((const uint8_t *)p);
            uint8x16_t cmp = vorrq_u8(vceqq_u8(chunk, v_delim), vceqq_u8(chunk, v_quote));
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
        if (c == delimiter || c == '"') return true;
        p++;
        len--;
    }
    return false;
}

static int cexport_write_quoted_string(const char *str, char *buf, size_t buf_size)
{
    if (!str) str = "";
    size_t need = 3; /* outer quotes and terminator */
    for (const char *p = str; *p; p++) {
        size_t add = *p == '"' ? 2 : 1;
        if (need > buf_size || add > buf_size - need) return -1;
        need += add;
    }
    if (need > buf_size) return -1;
    size_t pos = 0;
    buf[pos++] = '"';
    for (const char *p = str; *p; p++) {
        if (*p == '"') buf[pos++] = '"';
        buf[pos++] = *p;
    }
    buf[pos++] = '"';
    buf[pos] = '\0';
    return (int)pos;
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
            if (pos + CEXPORT_DBL_BUFFER_SIZE + 2 + (size_t)line_end_len > buf_size) {
                return -1;
            }

            /* Write directly to output buffer */
            field_len = cexport_double_to_str(val, buf + pos, CEXPORT_DBL_BUFFER_SIZE, false, vtype);
            pos += field_len;
        } else {
            /* String variable */
            const char *str = var->data.str[row_idx];
            if (str == NULL) str = "";

            /* Display-formatted numbers and dates are never quoted, even when
             * they contain the delimiter (native export behaves the same). */
            bool need_quote = ctx->vartypes[j] != VARTYPE_FORMATTED &&
                              (ctx->quote_strings ||
                               (ctx->quote_if_needed && cexport_string_needs_quoting(str, delimiter)));

            if (need_quote) {
                field_len = cexport_write_quoted_string(str, buf + pos, buf_size - pos);
                if (field_len < 0 || pos + (size_t)field_len + 2 + (size_t)line_end_len > buf_size) return -1;
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

/* Stata fixed strings hold at most 2045 bytes; SF_sdata adds a terminator. */
#define CEXPORT_STR_READ_SIZE 2046

int cexport_format_row_stream(const cexport_context *ctx, ST_int obs,
                              char *buf, size_t buf_size)
{
    size_t pos = 0;
    const size_t nvars = ctx->nvars;
    const char delimiter = ctx->delimiter;
    const int line_end_len = ctx->line_ending_len;
    ST_IIIDp read_number = (_stata_)->safevdata;
    ST_IIIS read_string = (_stata_)->sdata;

    for (size_t j = 0; j < nvars; j++) {
        ST_int var = (ST_int)(j + 1);
        int field_len;

        if (!ctx->is_string[j]) {
            double val;
            if (pos + CEXPORT_DBL_BUFFER_SIZE + 2 + (size_t)line_end_len > buf_size) return -1;
            if (read_number(var, obs, &val)) return -2;
            field_len = cexport_double_to_str(val, buf + pos, CEXPORT_DBL_BUFFER_SIZE, false, ctx->vartypes[j]);
        } else {
            /* Read directly into the output; quote in place only when needed. */
            if (pos + CEXPORT_STR_READ_SIZE + 2 + (size_t)line_end_len > buf_size) return -1;
            char *dst = buf + pos;
            if (read_string(var, obs, dst)) return -2;
            field_len = (int)strlen(dst);
            bool need_quote = ctx->vartypes[j] != VARTYPE_FORMATTED &&
                              (ctx->quote_strings ||
                               (ctx->quote_if_needed && cexport_string_needs_quoting(dst, delimiter)));
            if (need_quote) {
                char value[CEXPORT_STR_READ_SIZE];
                memcpy(value, dst, (size_t)field_len + 1);
                field_len = cexport_write_quoted_string(value, dst, buf_size - pos);
                if (field_len < 0 || pos + (size_t)field_len + 2 + (size_t)line_end_len > buf_size) return -1;
            }
        }
        pos += field_len;

        if (j < nvars - 1) {
            buf[pos++] = delimiter;
        } else {
            memcpy(buf + pos, ctx->line_ending, line_end_len);
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
        bool need_quote = (ctx->quote_if_needed && cexport_string_needs_quoting(name, ctx->delimiter));

        if (need_quote) {
            if (pos + name_len * 2 + 4 + (size_t)line_end_len > buf_size) {
                return -1;
            }
            int quoted_len = cexport_write_quoted_string(name, buf + pos, buf_size - pos);
            if (quoted_len < 0) return -1;
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

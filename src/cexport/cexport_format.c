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
#include <limits.h>
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
    /* log10 rounds to the nearest double, so just below a power of ten it
     * can return that power: 99.99999999999999 gave 2 and was written "100"
     * with a digit lost (native writes it in full). Powers of ten through
     * 1e19 are exact, so the decimal exponent is checked against them. */
    if (exponent >= 1 && exponent <= 19 && magnitude < (double)cexport_pow10[exponent])
        exponent--;
    else if (exponent >= -1 && exponent <= 18 && magnitude >= (double)cexport_pow10[exponent + 1])
        exponent++;
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

/* Returns the quoted length, -1 when buf is too small, or CEXPORT_ROW_TOO_LONG
 * when the quoted field alone exceeds INT_MAX bytes. */
static int cexport_write_quoted_string(const char *str, char *buf, size_t buf_size)
{
    if (!str) str = "";
    size_t need = 3; /* outer quotes and terminator */
    for (const char *p = str; *p; p++) {
        need += *p == '"' ? 2 : 1;
        if (need - 1 > (size_t)INT_MAX) return CEXPORT_ROW_TOO_LONG;
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
   Date and Time Display Formats

   export delimited writes a variable with a %t format as
   strtrim(string(value, format)). The formats compiled here are the daily,
   clock, weekly, monthly, quarterly, half-yearly and yearly ones built from
   the documented display codes below. Any other format is refused and stays
   with Stata, as do values outside the years 100-9999 (which Stata shows as
   numbers): the ado checks the variable's range before it sends a format.
   ======================================================================== */

#define CEXPORT_TIMEFMT_MAX 64

enum {
    TF_LITERAL, TF_CC, TF_cc, TF_YY, TF_yy, TF_JJJ, TF_jjj,
    TF_Mon, TF_Month, TF_mon, TF_month, TF_NN, TF_nn, TF_DD, TF_dd,
    TF_Dayname, TF_Day, TF_Da, TF_day, TF_da, TF_h, TF_q, TF_WW, TF_ww,
    TF_HH, TF_Hh, TF_hH, TF_hh, TF_MM, TF_mm, TF_SS, TF_ss,
    TF_am, TF_a_m, TF_AM, TF_A_M
};

struct cexport_timefmt {
    char unit;                                /* d c w m q h y */
    int count;
    unsigned char code[CEXPORT_TIMEFMT_MAX];
    char literal[CEXPORT_TIMEFMT_MAX];        /* for TF_LITERAL */
};

/* Display codes, a longer code before any code it begins with. clock codes
 * belong to %tc only. */
static const struct {
    const char *text;
    unsigned char code;
    bool clock;
} tf_codes[] = {
    {"Dayname", TF_Dayname, false}, {"Month", TF_Month, false}, {"month", TF_month, false},
    {"A.M.", TF_A_M, true}, {"a.m.", TF_a_m, true},
    {"JJJ", TF_JJJ, false}, {"jjj", TF_jjj, false},
    {"Mon", TF_Mon, false}, {"mon", TF_mon, false},
    {"Day", TF_Day, false}, {"day", TF_day, false},
    {"CC", TF_CC, false}, {"cc", TF_cc, false}, {"YY", TF_YY, false}, {"yy", TF_yy, false},
    {"NN", TF_NN, false}, {"nn", TF_nn, false}, {"DD", TF_DD, false}, {"dd", TF_dd, false},
    {"Da", TF_Da, false}, {"da", TF_da, false}, {"WW", TF_WW, false}, {"ww", TF_ww, false},
    {"HH", TF_HH, true}, {"Hh", TF_Hh, true}, {"hH", TF_hH, true}, {"hh", TF_hh, true},
    {"MM", TF_MM, true}, {"mm", TF_mm, true}, {"SS", TF_SS, true}, {"ss", TF_ss, true},
    {"AM", TF_AM, true}, {"am", TF_am, true},
    {"h", TF_h, false}, {"q", TF_q, false}
};

cexport_timefmt *cexport_timefmt_compile(const char *format)
{
    const char *p = format;
    if (!p || *p++ != '%') return NULL;
    if (*p == '-') p++;
    if (*p++ != 't') return NULL;
    const char unit = *p;
    if (!unit || !strchr("dcwmqhy", unit)) return NULL;
    p++;
    /* Stata's default for each unit */
    const char *detail = *p ? p
        : unit == 'd' ? "DDmonCCYY" : unit == 'c' ? "DDmonCCYY_HH:MM:SS"
        : unit == 'w' ? "CCYY!www" : unit == 'm' ? "CCYY!mnn"
        : unit == 'q' ? "CCYY!qq" : unit == 'h' ? "CCYY!hh" : "CCYY";

    cexport_timefmt *f = (cexport_timefmt *)calloc(1, sizeof(*f));
    if (!f) return NULL;
    f->unit = unit;
    for (p = detail; *p;) {
        if (*p == '+') { /* separates codes, displays nothing */
            p++;
            continue;
        }
        if (f->count >= CEXPORT_TIMEFMT_MAX) goto unsupported;
        char literal = 0;
        if (*p == '!') {
            /* !c displays the character c; a multibyte c stays with Stata */
            if (!p[1] || (unsigned char)p[1] >= 128) goto unsupported;
            literal = p[1];
            p += 2;
        } else if (*p == '_') {
            literal = ' ';
            p++;
        } else if (*p == '.' && p[1] == 's') {
            goto unsupported; /* fractions of a second */
        } else if (strchr(".,:-/\\", *p)) {
            literal = *p++;
        }
        if (literal) {
            f->code[f->count] = TF_LITERAL;
            f->literal[f->count++] = literal;
            continue;
        }
        size_t k, n = sizeof(tf_codes) / sizeof(tf_codes[0]);
        for (k = 0; k < n; k++)
            if (!strncmp(p, tf_codes[k].text, strlen(tf_codes[k].text))) break;
        if (k == n || (tf_codes[k].clock && unit != 'c')) goto unsupported;
        f->code[f->count++] = tf_codes[k].code;
        p += strlen(tf_codes[k].text);
    }
    return f;

unsupported:
    free(f);
    return NULL;
}

/* Bytes a cell in this format can take: no code displays more than 9. */
size_t cexport_timefmt_reserve(const cexport_timefmt *f)
{
    return 9 * (size_t)f->count + 2;
}

static inline long long tf_floor_div(long long a, long long b)
{
    long long q = a / b;
    return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
}

/* Days from 01jan1960 to a civil date and back (proleptic Gregorian). */
static long long tf_days(long long y, int m, int d)
{
    y -= m <= 2;
    const long long era = tf_floor_div(y, 400);
    const long long yoe = y - era * 400;
    const long long doy = (153 * (m > 2 ? m - 3 : m + 9) + 2) / 5 + d - 1;
    const long long doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + doe - 715815;
}

static void tf_civil(long long days, long long *y, int *m, int *d)
{
    const long long z = days + 715815;
    const long long era = tf_floor_div(z, 146097);
    const long long doe = z - era * 146097;
    const long long yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    const long long doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    const long long mp = (5 * doy + 2) / 153;
    *d = (int)(doy - (153 * mp + 2) / 5 + 1);
    *m = (int)(mp < 10 ? mp + 3 : mp - 9);
    *y = yoe + era * 400 + (*m <= 2);
}

static inline char *tf_put_text(char *p, const char *text, size_t n, bool lower)
{
    for (size_t i = 0; i < n; i++) {
        char c = text[i];
        *p++ = (lower && c >= 'A' && c <= 'Z') ? (char)(c - 'A' + 'a') : c;
    }
    return p;
}

static inline char *tf_put_number(char *p, long long v, int width)
{
    char digits[24];
    int n = 0;
    if (v < 0) {
        *p++ = '-';
        v = -v;
    }
    do {
        digits[n++] = (char)('0' + v % 10);
        v /= 10;
    } while (v);
    for (int i = n; i < width; i++) *p++ = '0';
    while (n) *p++ = digits[--n];
    return p;
}

/* The display text of value, without leading and trailing blanks. Returns
 * its length, or -1 when the value is not one this code displays (the
 * caller then writes the number). out holds cexport_timefmt_reserve(). */
static int cexport_timefmt_render(const cexport_timefmt *f, double value, char *out)
{
    static const char *const months[12] = {"January", "February", "March", "April", "May", "June", "July",
                                           "August", "September", "October", "November", "December"};
    static const char *const days[7] = {"Sunday", "Monday", "Tuesday", "Wednesday",
                                        "Thursday", "Friday", "Saturday"};
    if (!(value > -1e15 && value < 1e15)) return -1;
    const long long whole = (long long)floor(value);
    long long day = 0, year;
    int month, mday, hour = 0, minute = 0, second = 0;
    switch (f->unit) {
    case 'd':
        day = whole;
        break;
    case 'c': {
        day = tf_floor_div(whole, 86400000LL);
        const long long rest = whole - day * 86400000LL;
        hour = (int)(rest / 3600000);
        minute = (int)(rest / 60000 % 60);
        second = (int)(rest / 1000 % 60);
        break;
    }
    case 'w': {
        const long long yr = tf_floor_div(whole, 52);
        day = tf_days(1960 + yr, 1, 1) + 7 * (whole - yr * 52);
        break;
    }
    case 'm': {
        const long long yr = tf_floor_div(whole, 12);
        day = tf_days(1960 + yr, (int)(whole - yr * 12) + 1, 1);
        break;
    }
    case 'q': {
        const long long yr = tf_floor_div(whole, 4);
        day = tf_days(1960 + yr, 3 * (int)(whole - yr * 4) + 1, 1);
        break;
    }
    case 'h': {
        const long long yr = tf_floor_div(whole, 2);
        day = tf_days(1960 + yr, 6 * (int)(whole - yr * 2) + 1, 1);
        break;
    }
    default: /* 'y' */
        day = tf_days(whole, 1, 1);
        break;
    }
    tf_civil(day, &year, &month, &mday);
    if (year < 100 || year > 9999) return -1;
    const int doy = (int)(day - tf_days(year, 1, 1)) + 1;
    const int dow = (int)(day - tf_floor_div(day + 5, 7) * 7 + 5); /* 01jan1960 was a Friday */
    const int week = (doy - 1) / 7 + 1 > 52 ? 52 : (doy - 1) / 7 + 1;
    const int hour12 = hour % 12 ? hour % 12 : 12;

    char *p = out;
    for (int i = 0; i < f->count; i++) {
        switch (f->code[i]) {
        case TF_LITERAL: *p++ = f->literal[i]; break;
        case TF_CC: p = tf_put_number(p, year / 100, 2); break;
        case TF_cc: p = tf_put_number(p, year / 100, 0); break;
        case TF_YY: p = tf_put_number(p, year % 100, 2); break;
        case TF_yy: p = tf_put_number(p, year % 100, 0); break;
        case TF_JJJ: p = tf_put_number(p, doy, 3); break;
        case TF_jjj: p = tf_put_number(p, doy, 0); break;
        case TF_Mon: p = tf_put_text(p, months[month - 1], 3, false); break;
        case TF_mon: p = tf_put_text(p, months[month - 1], 3, true); break;
        case TF_Month: p = tf_put_text(p, months[month - 1], strlen(months[month - 1]), false); break;
        case TF_month: p = tf_put_text(p, months[month - 1], strlen(months[month - 1]), true); break;
        case TF_NN: p = tf_put_number(p, month, 2); break;
        case TF_nn: p = tf_put_number(p, month, 0); break;
        case TF_DD: p = tf_put_number(p, mday, 2); break;
        case TF_dd: p = tf_put_number(p, mday, 0); break;
        case TF_Dayname: p = tf_put_text(p, days[dow], strlen(days[dow]), false); break;
        case TF_Day: p = tf_put_text(p, days[dow], 3, false); break;
        case TF_Da: p = tf_put_text(p, days[dow], 2, false); break;
        case TF_day: p = tf_put_text(p, days[dow], 3, true); break;
        case TF_da: p = tf_put_text(p, days[dow], 2, true); break;
        case TF_h: p = tf_put_number(p, (month - 1) / 6 + 1, 0); break;
        case TF_q: p = tf_put_number(p, (month - 1) / 3 + 1, 0); break;
        case TF_WW: p = tf_put_number(p, week, 2); break;
        case TF_ww: p = tf_put_number(p, week, 0); break;
        case TF_HH: p = tf_put_number(p, hour, 2); break;
        case TF_hH: p = tf_put_number(p, hour, 0); break;
        case TF_Hh: p = tf_put_number(p, hour12, 2); break;
        case TF_hh: p = tf_put_number(p, hour12, 0); break;
        case TF_MM: p = tf_put_number(p, minute, 2); break;
        case TF_mm: p = tf_put_number(p, minute, 0); break;
        case TF_SS: p = tf_put_number(p, second, 2); break;
        case TF_ss: p = tf_put_number(p, second, 0); break;
        case TF_am: p = tf_put_text(p, hour < 12 ? "am" : "pm", 2, false); break;
        case TF_AM: p = tf_put_text(p, hour < 12 ? "AM" : "PM", 2, false); break;
        case TF_a_m: p = tf_put_text(p, hour < 12 ? "a.m." : "p.m.", 4, false); break;
        default: p = tf_put_text(p, hour < 12 ? "A.M." : "P.M.", 4, false); break;
        }
    }
    /* strtrim() */
    char *start = out;
    while (start < p && *start == ' ') start++;
    while (p > start && p[-1] == ' ') p--;
    const size_t len = (size_t)(p - start);
    if (start != out) memmove(out, start, len);
    out[len] = '\0';
    return (int)len;
}

/* The cell of a numeric variable written in its %t display format or as its
 * value label's text. buf holds at least ctx->cell_reserve[j] bytes. Returns
 * the length, or a negative row-formatter code. */
static int cexport_write_cell(const cexport_context *ctx, size_t j, double val, char *buf, size_t cap)
{
    const vartype_t vtype = ctx->vartypes[j];
    if (ctx->cell_kind[j] == CEXPORT_CELL_TIME) {
        /* A missing value is written as in a plain column: nothing for .,
         * the code for .a-.z. A display-formatted cell is never quoted. */
        if (!is_stata_missing(val)) {
            int len = cexport_timefmt_render(ctx->timefmt[j], val, buf);
            if (len >= 0) return len;
        }
        return cexport_double_to_str(val, buf, CEXPORT_DBL_BUFFER_SIZE, false, vtype);
    }

    /* Label text is a string cell, also when the text is empty (export
     * delimited then writes an empty cell, not the number). */
    const char *text = ctx->label_text(ctx->labelmap[j], val);
    if (text) {
        if (ctx->quote_strings ||
            (ctx->quote_if_needed && cexport_string_needs_quoting(text, ctx->delimiter)))
            return cexport_write_quoted_string(text, buf, cap);
        const size_t len = strlen(text);
        if (len + 1 > cap) return -1;
        memcpy(buf, text, len + 1);
        return (int)len;
    }
    /* A value without a label is the number of an unlabeled column. Under
     * quote it is written bare, as export delimited writes it; otherwise it
     * shares the label text's cell rules and is quoted when it holds the
     * delimiter. */
    char number[CEXPORT_DBL_BUFFER_SIZE];
    const int len = cexport_double_to_str(val, number, CEXPORT_DBL_BUFFER_SIZE, false, vtype);
    if (len < 0) return -1;
    if (!ctx->quote_strings && ctx->quote_if_needed &&
        cexport_string_needs_quoting(number, ctx->delimiter))
        return cexport_write_quoted_string(number, buf, cap);
    if ((size_t)len + 1 > cap) return -1;
    memcpy(buf, number, (size_t)len + 1);
    return len;
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
        int field_len;

        if (ctx->cell_kind && ctx->cell_kind[j]) {
            const size_t reserve = ctx->cell_reserve[j];
            if (pos + reserve + 2 + (size_t)line_end_len > buf_size) return -1;
            field_len = cexport_write_cell(ctx, j, val, buf + pos, reserve);
            if (field_len < 0) return field_len;
        } else {
            /* Check buffer space for number + delimiter/newline */
            if (pos + CEXPORT_DBL_BUFFER_SIZE + 2 + (size_t)line_end_len > buf_size) {
                return -1;
            }

            /* Write directly to output buffer */
            field_len = cexport_double_to_str(val, buf + pos, CEXPORT_DBL_BUFFER_SIZE, false, vtype);
        }
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

    return pos > (size_t)INT_MAX ? CEXPORT_ROW_TOO_LONG : (int)pos;
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

            if (ctx->cell_kind && ctx->cell_kind[j]) {
                const size_t reserve = ctx->cell_reserve[j];
                if (pos + reserve + 2 + (size_t)line_end_len > buf_size) return -1;
                field_len = cexport_write_cell(ctx, j, val, buf + pos, reserve);
                if (field_len < 0) return field_len;
            } else {
                /* Check buffer space (assume max 32 chars for number + line ending) */
                if (pos + CEXPORT_DBL_BUFFER_SIZE + 2 + (size_t)line_end_len > buf_size) {
                    return -1;
                }

                /* Write directly to output buffer */
                field_len = cexport_double_to_str(val, buf + pos, CEXPORT_DBL_BUFFER_SIZE, false, vtype);
            }
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
                if (field_len == CEXPORT_ROW_TOO_LONG) return field_len;
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

    /* Loaded strL cells can make a row longer than the int return reports. */
    return pos > (size_t)INT_MAX ? CEXPORT_ROW_TOO_LONG : (int)pos;
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
            const bool special = ctx->cell_kind && ctx->cell_kind[j];
            const size_t reserve = special ? ctx->cell_reserve[j] : CEXPORT_DBL_BUFFER_SIZE;
            if (pos + reserve + 2 + (size_t)line_end_len > buf_size) return -1;
            if (read_number(var, obs, &val)) return -2;
            if (special) {
                field_len = cexport_write_cell(ctx, j, val, buf + pos, reserve);
                if (field_len < 0) return field_len;
            } else {
                field_len = cexport_double_to_str(val, buf + pos, CEXPORT_DBL_BUFFER_SIZE, false, ctx->vartypes[j]);
            }
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

    return pos > (size_t)INT_MAX ? CEXPORT_ROW_TOO_LONG : (int)pos;
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

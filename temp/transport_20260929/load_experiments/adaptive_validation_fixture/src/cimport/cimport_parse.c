/*
 * cimport_parse.c
 * CSV parsing functions for cimport
 */

#include <string.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include "cimport_parse.h"
#include "../ctools_types.h"
#include "stplugin.h"

#include "../ctools_simd.h"

#include "cimport_locale.inc"

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

bool cimport_is_whitespace(char c)
{
    return c == ' ' || c == '\t' || c == '\r';
}

/* ============================================================================
 * SIMD-Accelerated Scanning
 * ============================================================================ */

#if CTOOLS_HAS_AVX2
static const char *cimport_find_delim_or_newline_simd(const char *ptr, const char *end, char delim)
{
    const __m256i v_newline = _mm256_set1_epi8('\n');
    const __m256i v_delim = _mm256_set1_epi8(delim);
    const __m256i v_quote = _mm256_set1_epi8('"');

    /* AVX2: Process 32 bytes at a time */
    while (ptr + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256((const __m256i *)ptr);
        __m256i cmp_nl = _mm256_cmpeq_epi8(chunk, v_newline);
        __m256i cmp_delim = _mm256_cmpeq_epi8(chunk, v_delim);
        __m256i cmp_quote = _mm256_cmpeq_epi8(chunk, v_quote);
        __m256i cmp = _mm256_or_si256(_mm256_or_si256(cmp_nl, cmp_delim), cmp_quote);
        int mask = _mm256_movemask_epi8(cmp);

        if (mask) {
            int pos = __builtin_ctz(mask);
            return ptr + pos;
        }
        ptr += 32;
    }

    /* SSE2: Handle 16-31 byte remainder */
    if (ptr + 16 <= end) {
        __m128i v_newline_128 = _mm_set1_epi8('\n');
        __m128i v_delim_128 = _mm_set1_epi8(delim);
        __m128i v_quote_128 = _mm_set1_epi8('"');

        __m128i chunk = _mm_loadu_si128((const __m128i *)ptr);
        __m128i cmp_nl = _mm_cmpeq_epi8(chunk, v_newline_128);
        __m128i cmp_delim = _mm_cmpeq_epi8(chunk, v_delim_128);
        __m128i cmp_quote = _mm_cmpeq_epi8(chunk, v_quote_128);
        __m128i cmp = _mm_or_si128(_mm_or_si128(cmp_nl, cmp_delim), cmp_quote);
        int mask = _mm_movemask_epi8(cmp);

        if (mask) {
            int pos = __builtin_ctz(mask);
            return ptr + pos;
        }
        ptr += 16;
    }

    /* Scalar tail for <16 bytes */
    while (ptr < end) {
        char c = *ptr;
        if (c == '\n' || c == delim || c == '"') return ptr;
        ptr++;
    }
    return end;
}

#elif CTOOLS_HAS_SSE2
static const char *cimport_find_delim_or_newline_simd(const char *ptr, const char *end, char delim)
{
    const __m128i v_newline = _mm_set1_epi8('\n');
    const __m128i v_delim = _mm_set1_epi8(delim);
    const __m128i v_quote = _mm_set1_epi8('"');

    while (ptr + 16 <= end) {
        __m128i chunk = _mm_loadu_si128((const __m128i *)ptr);
        __m128i cmp_nl = _mm_cmpeq_epi8(chunk, v_newline);
        __m128i cmp_delim = _mm_cmpeq_epi8(chunk, v_delim);
        __m128i cmp_quote = _mm_cmpeq_epi8(chunk, v_quote);
        __m128i cmp = _mm_or_si128(_mm_or_si128(cmp_nl, cmp_delim), cmp_quote);
        int mask = _mm_movemask_epi8(cmp);

        if (mask) {
            int pos = __builtin_ctz(mask);
            return ptr + pos;
        }
        ptr += 16;
    }

    while (ptr < end) {
        char c = *ptr;
        if (c == '\n' || c == delim || c == '"') return ptr;
        ptr++;
    }
    return end;
}

#elif CTOOLS_HAS_NEON
static const char *cimport_find_delim_or_newline_simd(const char *ptr, const char *end, char delim)
{
    const uint8x16_t v_newline = vdupq_n_u8('\n');
    const uint8x16_t v_delim = vdupq_n_u8(delim);
    const uint8x16_t v_quote = vdupq_n_u8('"');

    while (ptr + 16 <= end) {
        uint8x16_t chunk = vld1q_u8((const uint8_t *)ptr);
        uint8x16_t cmp_nl = vceqq_u8(chunk, v_newline);
        uint8x16_t cmp_delim = vceqq_u8(chunk, v_delim);
        uint8x16_t cmp_quote = vceqq_u8(chunk, v_quote);
        uint8x16_t cmp = vorrq_u8(vorrq_u8(cmp_nl, cmp_delim), cmp_quote);

        /* Use ctzll to find first matching byte position directly,
         * matching the technique used in find_next_row_loose NEON */
        uint64x2_t cmp64 = vreinterpretq_u64_u8(cmp);
        uint64_t lo = vgetq_lane_u64(cmp64, 0);
        uint64_t hi = vgetq_lane_u64(cmp64, 1);
        if (lo) {
            return ptr + (__builtin_ctzll(lo) >> 3);
        }
        if (hi) {
            return ptr + 8 + (__builtin_ctzll(hi) >> 3);
        }
        ptr += 16;
    }

    while (ptr < end) {
        char c = *ptr;
        if (c == '\n' || c == delim || c == '"') return ptr;
        ptr++;
    }
    return end;
}

#else
/* Scalar fallback */
static const char *cimport_find_delim_or_newline_simd(const char *ptr, const char *end, char delim)
{
    while (ptr < end) {
        char c = *ptr;
        if (c == '\n' || c == delim || c == '"') return ptr;
        ptr++;
    }
    return end;
}
#endif

/* ============================================================================
 * Row Boundary Detection
 * ============================================================================ */

static const char *cimport_find_next_row_loose(const char *ptr, const char *end)
{
#if CTOOLS_HAS_AVX2
    {
        const __m256i v_newline = _mm256_set1_epi8('\n');
        while (ptr + 32 <= end) {
            __m256i chunk = _mm256_loadu_si256((const __m256i *)ptr);
            __m256i cmp = _mm256_cmpeq_epi8(chunk, v_newline);
            int mask = _mm256_movemask_epi8(cmp);
            if (mask) {
                return ptr + __builtin_ctz(mask) + 1;
            }
            ptr += 32;
        }
        /* SSE2 tail for 16-31 byte remainder */
        if (ptr + 16 <= end) {
            __m128i v_nl = _mm_set1_epi8('\n');
            __m128i chunk = _mm_loadu_si128((const __m128i *)ptr);
            __m128i cmp = _mm_cmpeq_epi8(chunk, v_nl);
            int mask = _mm_movemask_epi8(cmp);
            if (mask) {
                return ptr + __builtin_ctz(mask) + 1;
            }
            ptr += 16;
        }
    }
#elif CTOOLS_HAS_SSE2
    {
        const __m128i v_newline = _mm_set1_epi8('\n');
        while (ptr + 16 <= end) {
            __m128i chunk = _mm_loadu_si128((const __m128i *)ptr);
            __m128i cmp = _mm_cmpeq_epi8(chunk, v_newline);
            int mask = _mm_movemask_epi8(cmp);
            if (mask) {
                return ptr + __builtin_ctz(mask) + 1;
            }
            ptr += 16;
        }
    }
#elif CTOOLS_HAS_NEON
    {
        const uint8x16_t v_newline = vdupq_n_u8('\n');
        while (ptr + 16 <= end) {
            uint8x16_t chunk = vld1q_u8((const uint8_t *)ptr);
            uint8x16_t cmp = vceqq_u8(chunk, v_newline);
            uint64x2_t cmp64 = vreinterpretq_u64_u8(cmp);
            uint64_t lo = vgetq_lane_u64(cmp64, 0);
            uint64_t hi = vgetq_lane_u64(cmp64, 1);
            if (lo) {
                return ptr + (__builtin_ctzll(lo) >> 3) + 1;
            }
            if (hi) {
                return ptr + 8 + (__builtin_ctzll(hi) >> 3) + 1;
            }
            ptr += 16;
        }
    }
#endif
    /* Scalar tail */
    while (ptr < end) {
        if (*ptr == '\n') {
            return ptr + 1;
        }
        ptr++;
    }
    return end;
}

/*
 * Strict-mode row boundary: the first newline outside quotes. The quote state
 * toggles on every quote byte, exactly as cimport_parse_row_fast does (an
 * escaped "" inside a quoted field toggles twice), so in_quotes gives the state
 * at ptr. Only boundary searches, header detection and skipped rows use this,
 * so a memchr scan is fast enough; block-wise quote counting must not skip a
 * newline between a closing quote and a later opening quote.
 */
const char *cimport_find_next_row_strict_from(const char *ptr, const char *end,
                                              char quote, bool in_quotes)
{
    const char *nl = NULL;   /* next newline at or after ptr; end if none */
    while (ptr < end) {
        if (in_quotes) {
            const char *q = memchr(ptr, quote, (size_t)(end - ptr));
            if (!q) return end;
            ptr = q + 1;
            in_quotes = false;
            continue;
        }
        if (nl == NULL || nl < ptr) {
            nl = memchr(ptr, '\n', (size_t)(end - ptr));
            if (!nl) nl = end;
        }
        const char *q = memchr(ptr, quote, (size_t)(nl - ptr));
        if (q) {
            ptr = q + 1;
            in_quotes = true;
            continue;
        }
        return nl < end ? nl + 1 : end;
    }
    return end;
}

const char *cimport_find_next_row_strict(const char *ptr, const char *end, char quote)
{
    return cimport_find_next_row_strict_from(ptr, end, quote, false);
}

const char *cimport_find_next_row(const char *ptr, const char *end, char quote, CImportBindQuotesMode bindquotes)
{
    if (bindquotes == CIMPORT_BINDQUOTES_STRICT) {
        return cimport_find_next_row_strict(ptr, end, quote);
    } else {
        return cimport_find_next_row_loose(ptr, end);
    }
}

/* ============================================================================
 * Field Parsing
 * ============================================================================ */

int cimport_parse_row_fast(const char *start, const char *end, char delim, char quote,
                           CImportFieldRef *fields, int max_fields, const char *file_base,
                           CImportBindQuotesMode bindquotes,
                           const char **next_row_out, bool *had_unmatched_quote)
{
    int field_count = 0;
    const char *ptr = start;
    const char *field_start = start;
    bool in_quotes = false;
    bool current_field_has_quote = false;

    while (ptr < end && field_count < max_fields) {
        if (!in_quotes) {
            const char *found = cimport_find_delim_or_newline_simd(ptr, end, delim);

            if (found < end) {
                char c = *found;

                if (c == '"' && bindquotes == CIMPORT_BINDQUOTES_NOBIND) {
                    current_field_has_quote = true;
                    ptr = found + 1;
                    continue;
                }
                if (c == '"') {
                    in_quotes = true;
                    current_field_has_quote = true;
                    ptr = found + 1;
                    continue;
                }

                fields[field_count].offset = (uint64_t)(field_start - file_base);
                fields[field_count].length = (uint32_t)(found - field_start)
                    | (current_field_has_quote ? CIMPORT_FIELD_QUOTED_FLAG : 0);
                field_count++;
                current_field_has_quote = false;

                if (c == '\n') {
                    /* Row complete — newline found outside quotes */
                    if (next_row_out) *next_row_out = found + 1;
                    if (had_unmatched_quote) *had_unmatched_quote = false;
                    return field_count;
                }

                field_start = found + 1;
                ptr = field_start;
                continue;
            }
            ptr = end;
        } else {
            char c = *ptr;
            if (c == '\n' && bindquotes != CIMPORT_BINDQUOTES_STRICT) {
                /* LOOSE mode: newline always terminates row, even inside quotes.
                 * Record current field up to the newline. */
                fields[field_count].offset = (uint64_t)(field_start - file_base);
                fields[field_count].length = (uint32_t)(ptr - field_start)
                    | (current_field_has_quote ? CIMPORT_FIELD_QUOTED_FLAG : 0);
                field_count++;
                if (next_row_out) *next_row_out = ptr + 1;
                if (had_unmatched_quote) *had_unmatched_quote = true;
                return field_count;
            }
            if (c == quote) {
                if (ptr + 1 < end && *(ptr + 1) == quote) {
                    ptr += 2;
                    continue;
                }
                in_quotes = false;
            }
            ptr++;
        }
    }

    /* Trailing field: reached end of buffer or max_fields */
    if (field_start < end && field_count < max_fields) {
        const char *field_end = end;
        while (field_end > field_start && (field_end[-1] == '\r' || field_end[-1] == '\n')) {
            field_end--;
        }
        if (field_end > field_start) {
            fields[field_count].offset = (uint64_t)(field_start - file_base);
            fields[field_count].length = (uint32_t)(field_end - field_start)
                | (current_field_has_quote ? CIMPORT_FIELD_QUOTED_FLAG : 0);
            field_count++;
        }
    }

    /* Set fused-mode outputs */
    if (next_row_out) {
        if (field_count >= max_fields && ptr < end) {
            /* Stopped due to max_fields (outside quotes): skip the rest of
             * the row. Strict mode keeps quoted newlines inside the row, as
             * the chunk-boundary search assumes. */
            if (bindquotes == CIMPORT_BINDQUOTES_STRICT) {
                *next_row_out = cimport_find_next_row_strict_from(ptr, end, quote, false);
            } else {
                while (ptr < end && *ptr != '\n') ptr++;
                *next_row_out = (ptr < end) ? ptr + 1 : end;
            }
        } else {
            *next_row_out = end;
        }
    }
    if (had_unmatched_quote) *had_unmatched_quote = in_quotes;

    return field_count;
}

int cimport_extract_field_fast(const char *file_base, CImportFieldRef *field,
                                char *output, int max_len, char quote)
{
    const char *src = file_base + field->offset;
    int src_len = CIMPORT_FIELD_LENGTH(*field);
    int out_len = 0;
    bool strip_all = quote == '\1';
    /* output == NULL counts the decoded length (no buffer bound). */
    char discard;
    const bool count_only = output == NULL;
    if (count_only) { output = &discard; max_len = INT_MAX; }
#define CIMPORT_PUT(c) do { if (!count_only) output[out_len] = (c); out_len++; } while (0)
    if (strip_all) quote = '"';

    /* Only strip trailing CR/LF (not spaces/tabs - those are preserved like Stata) */
    while (src_len > 0 && (src[src_len-1] == '\r' || src[src_len-1] == '\n')) {
        src_len--;
    }
    if(strip_all || quote=='\2') {
        for(int i=0;i<src_len && out_len<max_len-1;i++) {
            if(strip_all && src[i]=='"')continue;
            CIMPORT_PUT(src[i]);
            if(quote=='\2' && src[i]=='"' && i+1<src_len && src[i+1]=='"')i++;
        }
        if (!count_only) output[out_len]=0;
        return out_len;
    }

    /* Check for quotes - need to look past leading whitespace to find them */
    const char *trimmed_src = src;
    int trimmed_len = src_len;
    while (trimmed_len > 0 && (*trimmed_src == ' ' || *trimmed_src == '\t')) {
        trimmed_src++;
        trimmed_len--;
    }
    while (trimmed_len > 0 && (trimmed_src[trimmed_len-1] == ' ' || trimmed_src[trimmed_len-1] == '\t')) {
        trimmed_len--;
    }

    bool starts_with_quote = (trimmed_len >= 1 && trimmed_src[0] == quote);
    bool ends_with_quote = (trimmed_len >= 1 && trimmed_src[trimmed_len-1] == quote);
    bool is_quoted = (trimmed_len >= 2 && starts_with_quote && ends_with_quote);

    /* Special case: field is just a single quote - treat as empty (matches Stata) */
    if (trimmed_len == 1 && trimmed_src[0] == quote) {
        if (!count_only) output[0] = '\0';
        return 0;
    }

    if (is_quoted) {
        /* Properly quoted field - strip quotes and handle escaped quotes */
        trimmed_src++;
        trimmed_len -= 2;

        for (int i = 0; i < trimmed_len && out_len < max_len - 1; i++) {
            if (trimmed_src[i] == quote && i + 1 < trimmed_len && trimmed_src[i + 1] == quote) {
                CIMPORT_PUT(quote);
                i++;
            } else {
                CIMPORT_PUT(trimmed_src[i]);
            }
        }
    } else if (starts_with_quote && !ends_with_quote) {
        /* Orphan leading quote - strip it */
        trimmed_src++;
        trimmed_len--;
        for (int i = 0; i < trimmed_len && out_len < max_len - 1; i++) {
            if (trimmed_src[i] == quote && i + 1 < trimmed_len && trimmed_src[i + 1] == quote) {
                CIMPORT_PUT(quote);
                i++;
            } else {
                CIMPORT_PUT(trimmed_src[i]);
            }
        }
    } else if (strip_all && !starts_with_quote && ends_with_quote) {
        /* Orphan trailing quote - strip it */
        trimmed_len--;
        for (int i = 0; i < trimmed_len && out_len < max_len - 1; i++) {
            if (trimmed_src[i] == quote && i + 1 < trimmed_len && trimmed_src[i + 1] == quote) {
                CIMPORT_PUT(quote);
                i++;
            } else {
                CIMPORT_PUT(trimmed_src[i]);
            }
        }
    } else {
        /* No quotes - copy as-is from ORIGINAL src (preserving whitespace like Stata) */
        int copy_len = (src_len < max_len - 1) ? src_len : max_len - 1;
        if (!count_only) memcpy(output, src, copy_len);
        out_len = copy_len;
    }

    if (!count_only) output[out_len] = '\0';
    return out_len;
#undef CIMPORT_PUT
}

/* ============================================================================
 * Type Detection
 * ============================================================================ */

bool cimport_field_looks_numeric_sep(const char *src, int len, char dec_sep, char grp_sep)
{
    while (len > 0 && (*src == ' ' || *src == '\t')) { src++; len--; }
    if (len >= 2 && src[0] == '"' && src[len-1] == '"') { src++; len-=2; }
    return cimport_field_looks_unquoted_numeric_sep(src,len,dec_sep,grp_sep);
}

bool cimport_field_looks_unquoted_numeric_sep(const char *src, int len, char dec_sep, char grp_sep)
{
    /* Line-ending bytes are not field content; an empty field is missing.
     * Native import reads a field of only spaces/tabs as text. */
    while (len > 0 && (src[len-1] == '\r' || src[len-1] == '\n')) len--;
    if (len == 0) return true;
    while (len > 0 && (*src == ' ' || *src == '\t')) { src++; len--; }
    if (len == 0) return false;

    /* Skip trailing whitespace */
    while (len > 0 && (src[len-1] == ' ' || src[len-1] == '\t' ||
                       src[len-1] == '\r' || src[len-1] == '\n')) { len--; }
    if (len == 0) return true;

    /* Stata extended missing values remain numeric; NA/NaN are strings. */
    if (len == 1 && *src == '.') return true;
    if (len == 2 && src[0] == '.' && src[1] >= 'a' && src[1] <= 'z') return true;

    if(cimport_locale_enabled) {
        double value;return cimport_locale_number(src,len,&value,SV_missval);
    }

    /* Scan the ENTIRE field to validate numeric format */
    const char *p = src;
    const char *end = src + len;
    bool has_decimal = false;
    bool has_digits = false;

    /* Optional leading sign */
    if (p < end && (*p == '-' || *p == '+')) p++;

    /* Must have at least one character left */
    if (p >= end) return false;

    /* Scan digits, optional decimal, and optional group separators */
    while (p < end) {
        char c = *p;
        if (c >= '0' && c <= '9') {
            has_digits = true;
            p++;
        } else if (c == dec_sep && !has_decimal) {
            has_decimal = true;
            p++;
        } else if (grp_sep != '\0' && c == grp_sep) {
            /* Group separator (e.g., thousand separator) - skip it */
            p++;
        } else if (c == 'e' || c == 'E') {
            /* Scientific notation - validate exponent */
            if (!has_digits) return false;  /* Need digits before E */
            p++;
            if (p < end && (*p == '-' || *p == '+')) p++;
            if (p >= end || *p < '0' || *p > '9') return false;  /* Need exponent digits */
            while (p < end && *p >= '0' && *p <= '9') p++;
            break;  /* End of number */
        } else {
            /* Invalid character - not numeric */
            return false;
        }
    }

    /* Must have consumed all input and seen at least one digit */
    return (p == end) && has_digits;
}

bool cimport_parse_number(const char *src, int len, double *value, double missing,
                          char decimal, char group)
{
    while (len && (*src == ' ' || *src == '\t')) { src++; len--; }
    while (len && (src[len-1] == ' ' || src[len-1] == '\t' || src[len-1] == '\r' || src[len-1] == '\n')) len--;
    if (len >= 2 && src[0] == '"' && src[len-1] == '"') { src++; len -= 2; }
    return cimport_parse_unquoted_number(src,len,value,missing,decimal,group);
}

bool cimport_parse_unquoted_number(const char *src, int len, double *value, double missing,
                                   char decimal, char group)
{
    while (len && (src[len-1] == '\r' || src[len-1] == '\n')) len--;
    int unpadded = len;
    while (len && (*src == ' ' || *src == '\t')) { src++; len--; }
    while (len && (src[len-1] == ' ' || src[len-1] == '\t')) len--;
    if (!len || (len == 1 && src[0] == '.')) { *value = missing; return true; }
    /* Native keeps .a-.z only when unpadded; " .a" is system missing. */
    if (len == 2 && src[0] == '.' && src[1] >= 'a' && src[1] <= 'z' && len != unpadded) {
        *value = missing;
        return true;
    }
    if (len == 2 && src[0] == '.' && src[1] >= 'a' && src[1] <= 'z') {
        uint64_t bits;
        memcpy(&bits, &missing, sizeof(bits));
        bits += (uint64_t)(src[1] - 'a' + 1) << 40;
        memcpy(value, &bits, sizeof(bits));
        return true;
    }
    if(cimport_locale_enabled)return cimport_locale_number(src,len,value,missing);
    return ctools_parse_double_with_separators(src, len, value, missing, decimal, group);
}

bool cimport_analyze_numeric_with_sep(const char *file_base, CImportFieldRef *field, char quote,
                                       char dec_sep, char grp_sep,
                                       double *out_value, bool *out_is_integer)
{
    const char *src = file_base + field->offset;
    int len = CIMPORT_FIELD_LENGTH(*field);

    while (len > 0 && (*src == ' ' || *src == '\t' || *src == quote)) { src++; len--; }
    while (len > 0 && (src[len-1] == ' ' || src[len-1] == '\t' || src[len-1] == quote ||
                       src[len-1] == '\r' || src[len-1] == '\n')) { len--; }

    if (len == 0) return false;

    /* Single decimal separator = missing */
    if (len == 1 && *src == dec_sep) return false;
    if (len == 2 && (src[0] == 'N' || src[0] == 'n') && (src[1] == 'A' || src[1] == 'a')) return false;
    if (len == 3 && (src[0] == 'N' || src[0] == 'n') && (src[1] == 'a' || src[1] == 'A') &&
        (src[2] == 'N' || src[2] == 'n')) return false;

    double val;
    if (!cimport_parse_number(src, len, &val, SV_missval, dec_sep, grp_sep)) return false;

    if (!isfinite(val) || val >= SV_missval) return false;
    *out_value = val;
    *out_is_integer = (floor(val) == val);

    return true;
}

bool cimport_locale_active(void)
{
    return cimport_locale_enabled;
}

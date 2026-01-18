/*
 * cimport_parse.c
 * CSV parsing functions for cimport
 */

#include <string.h>
#include "cimport_parse.h"
#include "../ctools_types.h"
#include "stplugin.h"

/* SIMD headers */
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    #include <arm_neon.h>
    #define CIMPORT_USE_NEON 1
#elif defined(__SSE2__)
    #include <emmintrin.h>
    #define CIMPORT_USE_SSE2 1
#endif

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

bool cimport_is_whitespace(char c)
{
    return c == ' ' || c == '\t' || c == '\r';
}

bool cimport_row_has_unmatched_quote(const char *start, const char *end, char quote)
{
    int quote_count = 0;
    for (const char *p = start; p < end; p++) {
        if (*p == quote) quote_count++;
    }
    return (quote_count % 2) == 1;
}

/* ============================================================================
 * SIMD-Accelerated Scanning
 * ============================================================================ */

#if CIMPORT_USE_NEON
const char *cimport_find_delim_or_newline_simd(const char *ptr, const char *end, char delim)
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

        uint64x2_t cmp64 = vreinterpretq_u64_u8(cmp);
        if (vgetq_lane_u64(cmp64, 0) || vgetq_lane_u64(cmp64, 1)) {
            for (int i = 0; i < 16 && ptr + i < end; i++) {
                char c = ptr[i];
                if (c == '\n' || c == delim || c == '"') return ptr + i;
            }
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

#elif CIMPORT_USE_SSE2
const char *cimport_find_delim_or_newline_simd(const char *ptr, const char *end, char delim)
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

#else
/* Scalar fallback */
const char *cimport_find_delim_or_newline_simd(const char *ptr, const char *end, char delim)
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

const char *cimport_find_next_row_loose(const char *ptr, const char *end)
{
    while (ptr < end) {
        if (*ptr == '\n') {
            return ptr + 1;
        }
        ptr++;
    }
    return end;
}

const char *cimport_find_next_row_strict(const char *ptr, const char *end, char quote)
{
    bool in_quotes = false;

    while (ptr < end) {
        char c = *ptr;

        if (c == quote) {
            in_quotes = !in_quotes;
        } else if (!in_quotes && c == '\n') {
            return ptr + 1;
        }
        ptr++;
    }
    return end;
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
                           CImportFieldRef *fields, int max_fields, const char *file_base)
{
    int field_count = 0;
    const char *ptr = start;
    const char *field_start = start;
    bool in_quotes = false;

    while (ptr < end && field_count < max_fields) {
        if (!in_quotes) {
            const char *found = cimport_find_delim_or_newline_simd(ptr, end, delim);

            if (found < end) {
                char c = *found;

                if (c == '"') {
                    in_quotes = true;
                    ptr = found + 1;
                    continue;
                }

                fields[field_count].offset = (uint64_t)(field_start - file_base);
                fields[field_count].length = (uint32_t)(found - field_start);
                field_count++;

                if (c == '\r' || c == '\n') {
                    return field_count;
                }

                field_start = found + 1;
                ptr = field_start;
                continue;
            }
            ptr = end;
        } else {
            char c = *ptr;
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

    if (field_start < end && field_count < max_fields) {
        const char *field_end = end;
        while (field_end > field_start && (field_end[-1] == '\r' || field_end[-1] == '\n')) {
            field_end--;
        }
        if (field_end > field_start) {
            fields[field_count].offset = (uint64_t)(field_start - file_base);
            fields[field_count].length = (uint32_t)(field_end - field_start);
            field_count++;
        }
    }

    return field_count;
}

bool cimport_field_contains_quote(const char *src, int len, char quote)
{
    /* Unrolled 8-byte check for performance */
    while (len >= 8) {
        if (src[0] == quote || src[1] == quote || src[2] == quote || src[3] == quote ||
            src[4] == quote || src[5] == quote || src[6] == quote || src[7] == quote) {
            return true;
        }
        src += 8;
        len -= 8;
    }
    while (len > 0) {
        if (*src == quote) return true;
        src++;
        len--;
    }
    return false;
}

int cimport_extract_field_fast(const char *file_base, CImportFieldRef *field,
                                char *output, int max_len, char quote)
{
    const char *src = file_base + field->offset;
    int src_len = field->length;
    int out_len = 0;

    /* Only strip trailing CR/LF (not spaces/tabs - those are preserved like Stata) */
    while (src_len > 0 && (src[src_len-1] == '\r' || src[src_len-1] == '\n')) {
        src_len--;
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
        output[0] = '\0';
        return 0;
    }

    if (is_quoted) {
        /* Properly quoted field - strip quotes and handle escaped quotes */
        trimmed_src++;
        trimmed_len -= 2;

        for (int i = 0; i < trimmed_len && out_len < max_len - 1; i++) {
            if (trimmed_src[i] == quote && i + 1 < trimmed_len && trimmed_src[i + 1] == quote) {
                output[out_len++] = quote;
                i++;
            } else {
                output[out_len++] = trimmed_src[i];
            }
        }
    } else if (starts_with_quote && !ends_with_quote) {
        /* Orphan leading quote - strip it */
        trimmed_src++;
        trimmed_len--;
        for (int i = 0; i < trimmed_len && out_len < max_len - 1; i++) {
            if (trimmed_src[i] == quote && i + 1 < trimmed_len && trimmed_src[i + 1] == quote) {
                output[out_len++] = quote;
                i++;
            } else {
                output[out_len++] = trimmed_src[i];
            }
        }
    } else if (!starts_with_quote && ends_with_quote) {
        /* Orphan trailing quote - strip it */
        trimmed_len--;
        for (int i = 0; i < trimmed_len && out_len < max_len - 1; i++) {
            if (trimmed_src[i] == quote && i + 1 < trimmed_len && trimmed_src[i + 1] == quote) {
                output[out_len++] = quote;
                i++;
            } else {
                output[out_len++] = trimmed_src[i];
            }
        }
    } else {
        /* No quotes - copy as-is from ORIGINAL src (preserving whitespace like Stata) */
        int copy_len = (src_len < max_len - 1) ? src_len : max_len - 1;
        memcpy(output, src, copy_len);
        out_len = copy_len;
    }

    output[out_len] = '\0';
    return out_len;
}

int cimport_extract_field_unquoted(const char *file_base, CImportFieldRef *field,
                                    char *output, int max_len)
{
    const char *src = file_base + field->offset;
    int src_len = field->length;

    /* Only strip trailing CR/LF - preserve spaces/tabs like Stata */
    while (src_len > 0 && (src[src_len-1] == '\r' || src[src_len-1] == '\n')) {
        src_len--;
    }

    int copy_len = (src_len < max_len - 1) ? src_len : max_len - 1;
    memcpy(output, src, copy_len);
    output[copy_len] = '\0';
    return copy_len;
}

/* ============================================================================
 * Type Detection
 * ============================================================================ */

bool cimport_field_looks_numeric(const char *src, int len)
{
    while (len > 0 && (*src == ' ' || *src == '\t')) { src++; len--; }
    if (len == 0) return true;
    if (*src == '"' && len >= 2) { src++; len -= 2; }
    if (len == 0) return true;

    char c = *src;
    if (c >= '0' && c <= '9') return true;
    if (c == '-' || c == '+' || c == '.') return true;

    if (len >= 2 && (c == 'N' || c == 'n')) {
        char c2 = src[1];
        if (c2 == 'A' || c2 == 'a') return true;
        if (len >= 3 && (c2 == 'a' || c2 == 'A')) {
            char c3 = src[2];
            if (c3 == 'N' || c3 == 'n') return true;
        }
    }

    return false;
}

bool cimport_analyze_numeric_fast(const char *file_base, CImportFieldRef *field, char quote,
                                   double *out_value, bool *out_is_integer)
{
    const char *src = file_base + field->offset;
    int len = field->length;

    while (len > 0 && (*src == ' ' || *src == '\t' || *src == quote)) { src++; len--; }
    while (len > 0 && (src[len-1] == ' ' || src[len-1] == '\t' || src[len-1] == quote ||
                       src[len-1] == '\r' || src[len-1] == '\n')) { len--; }

    if (len == 0) return false;

    if (len == 1 && *src == '.') return false;
    if (len == 2 && (src[0] == 'N' || src[0] == 'n') && (src[1] == 'A' || src[1] == 'a')) return false;
    if (len == 3 && (src[0] == 'N' || src[0] == 'n') && (src[1] == 'a' || src[1] == 'A') &&
        (src[2] == 'N' || src[2] == 'n')) return false;

    double val;
    if (!ctools_parse_double_fast(src, len, &val, SV_missval)) return false;

    *out_value = val;

    *out_is_integer = true;
    for (int i = 0; i < len; i++) {
        if (src[i] == '.') {
            *out_is_integer = false;
            break;
        }
    }

    return true;
}

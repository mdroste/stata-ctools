/*
 * cimport_parse.c
 * CSV parsing functions for cimport
 */

#include <string.h>
#include "cimport_parse.h"
#include "../ctools_types.h"
#include "stplugin.h"

#include "../ctools_simd.h"

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

#if CTOOLS_HAS_AVX2
const char *cimport_find_delim_or_newline_simd(const char *ptr, const char *end, char delim)
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

#elif CTOOLS_HAS_NEON
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

const char *cimport_find_next_row_strict(const char *ptr, const char *end, char quote)
{
    bool in_quotes = false;

#if CTOOLS_HAS_AVX2
    {
        const __m256i v_newline = _mm256_set1_epi8('\n');
        const __m256i v_quote = _mm256_set1_epi8(quote);

        while (ptr + 32 <= end) {
            __m256i chunk = _mm256_loadu_si256((const __m256i *)ptr);
            __m256i cmp_q = _mm256_cmpeq_epi8(chunk, v_quote);
            int q_mask = _mm256_movemask_epi8(cmp_q);

            if (!in_quotes) {
                /* Not in quotes: scan for newline or quote */
                __m256i cmp_nl = _mm256_cmpeq_epi8(chunk, v_newline);
                int nl_mask = _mm256_movemask_epi8(cmp_nl);
                int combined = q_mask | nl_mask;
                if (!combined) {
                    ptr += 32;
                    continue;
                }
                /* Process first interesting character */
                int pos = __builtin_ctz(combined);
                if (nl_mask && (!q_mask || (__builtin_ctz(nl_mask) < __builtin_ctz(q_mask)))) {
                    /* Newline comes first */
                    return ptr + __builtin_ctz(nl_mask) + 1;
                }
                /* Quote comes first — enter quoted mode */
                ptr += pos + 1;
                in_quotes = true;
            } else {
                /* In quotes: SIMD bulk quote counting */
                if (!q_mask) {
                    /* No quotes in this chunk — skip it entirely */
                    ptr += 32;
                    continue;
                }
                int quote_count = __builtin_popcount(q_mask);
                if ((quote_count & 1) == 0) {
                    /* Even number of quotes — state unchanged, skip chunk */
                    ptr += 32;
                    continue;
                }
                /* Odd quotes — state flips. Find the last quote. */
                int last_quote_pos = 31 - __builtin_clz(q_mask);
                ptr += last_quote_pos + 1;
                in_quotes = false;
                /* Check for newline between last quote+1 and end of chunk */
                /* (will be caught by the next iteration or scalar tail) */
            }
        }
    }
#elif CTOOLS_HAS_SSE2
    {
        const __m128i v_newline = _mm_set1_epi8('\n');
        const __m128i v_quote = _mm_set1_epi8(quote);

        while (ptr + 16 <= end) {
            __m128i chunk = _mm_loadu_si128((const __m128i *)ptr);
            __m128i cmp_q = _mm_cmpeq_epi8(chunk, v_quote);
            int q_mask = _mm_movemask_epi8(cmp_q);

            if (!in_quotes) {
                __m128i cmp_nl = _mm_cmpeq_epi8(chunk, v_newline);
                int nl_mask = _mm_movemask_epi8(cmp_nl);
                int combined = q_mask | nl_mask;
                if (!combined) {
                    ptr += 16;
                    continue;
                }
                int pos = __builtin_ctz(combined);
                if (nl_mask && (!q_mask || (__builtin_ctz(nl_mask) < __builtin_ctz(q_mask)))) {
                    return ptr + __builtin_ctz(nl_mask) + 1;
                }
                ptr += pos + 1;
                in_quotes = true;
            } else {
                if (!q_mask) {
                    ptr += 16;
                    continue;
                }
                int quote_count = __builtin_popcount(q_mask);
                if ((quote_count & 1) == 0) {
                    ptr += 16;
                    continue;
                }
                int last_quote_pos = 15 - __builtin_clz(q_mask << 16);
                ptr += last_quote_pos + 1;
                in_quotes = false;
            }
        }
    }
#elif CTOOLS_HAS_NEON
    {
        const uint8x16_t v_newline = vdupq_n_u8('\n');
        const uint8x16_t v_quote = vdupq_n_u8(quote);

        while (ptr + 16 <= end) {
            uint8x16_t chunk = vld1q_u8((const uint8_t *)ptr);
            uint8x16_t cmp_q = vceqq_u8(chunk, v_quote);

            /* Convert NEON comparison to a bitmask (1 bit per byte) */
            /* Use shift + narrow + reinterpret approach for byte mask */
            uint64x2_t q64 = vreinterpretq_u64_u8(cmp_q);
            uint64_t q_lo = vgetq_lane_u64(q64, 0);
            uint64_t q_hi = vgetq_lane_u64(q64, 1);
            /* Extract bit-per-byte mask: each matching byte is 0xFF,
             * so test high bit of each byte */
            uint64_t q_mask_lo = q_lo & 0x8080808080808080ULL;
            uint64_t q_mask_hi = q_hi & 0x8080808080808080ULL;
            bool has_quotes = (q_mask_lo | q_mask_hi) != 0;

            if (!in_quotes) {
                uint8x16_t cmp_nl = vceqq_u8(chunk, v_newline);
                uint8x16_t cmp = vorrq_u8(cmp_nl, cmp_q);
                uint64x2_t cmp64 = vreinterpretq_u64_u8(cmp);
                uint64_t lo = vgetq_lane_u64(cmp64, 0);
                uint64_t hi = vgetq_lane_u64(cmp64, 1);
                if (!lo && !hi) {
                    ptr += 16;
                    continue;
                }
                /* Find first match position */
                int pos;
                if (lo) {
                    pos = __builtin_ctzll(lo) >> 3;
                } else {
                    pos = 8 + (__builtin_ctzll(hi) >> 3);
                }
                /* Check if it's a newline */
                uint64x2_t nl64 = vreinterpretq_u64_u8(cmp_nl);
                uint64_t nl_lo = vgetq_lane_u64(nl64, 0);
                uint64_t nl_hi = vgetq_lane_u64(nl64, 1);
                int nl_pos = 16;
                if (nl_lo) nl_pos = __builtin_ctzll(nl_lo) >> 3;
                else if (nl_hi) nl_pos = 8 + (__builtin_ctzll(nl_hi) >> 3);

                if (nl_pos <= pos) {
                    return ptr + nl_pos + 1;
                }
                /* Quote comes first */
                ptr += pos + 1;
                in_quotes = true;
            } else {
                /* In quotes: SIMD bulk quote counting */
                if (!has_quotes) {
                    ptr += 16;
                    continue;
                }
                /* Count quote bytes using popcount on the high-bit mask */
                int quote_count = __builtin_popcountll(q_mask_lo) + __builtin_popcountll(q_mask_hi);
                if ((quote_count & 1) == 0) {
                    ptr += 16;
                    continue;
                }
                /* Odd quotes — find last quote position */
                int last_pos;
                if (q_mask_hi) {
                    last_pos = 8 + (63 - __builtin_clzll(q_mask_hi)) / 8;
                } else {
                    last_pos = (63 - __builtin_clzll(q_mask_lo)) / 8;
                }
                ptr += last_pos + 1;
                in_quotes = false;
            }
        }
    }
#endif

    /* Scalar tail */
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
            /* Stopped due to max_fields: scan forward for newline */
            while (ptr < end && *ptr != '\n') ptr++;
            *next_row_out = (ptr < end) ? ptr + 1 : end;
        } else {
            *next_row_out = end;
        }
    }
    if (had_unmatched_quote) *had_unmatched_quote = in_quotes;

    return field_count;
}

bool cimport_field_contains_quote(const char *src, int len, char quote)
{
#if CTOOLS_HAS_AVX2
    /* AVX2: Scan 32 bytes at a time */
    const __m256i v_quote = _mm256_set1_epi8(quote);

    while (len >= 32) {
        __m256i chunk = _mm256_loadu_si256((const __m256i *)src);
        __m256i cmp = _mm256_cmpeq_epi8(chunk, v_quote);
        if (_mm256_movemask_epi8(cmp)) {
            return true;
        }
        src += 32;
        len -= 32;
    }

    /* SSE2: Handle 16-31 byte remainder */
    if (len >= 16) {
        __m128i v_quote_128 = _mm_set1_epi8(quote);
        __m128i chunk = _mm_loadu_si128((const __m128i *)src);
        __m128i cmp = _mm_cmpeq_epi8(chunk, v_quote_128);
        if (_mm_movemask_epi8(cmp)) {
            return true;
        }
        src += 16;
        len -= 16;
    }
#elif CTOOLS_HAS_SSE2
    /* SSE2: Scan 16 bytes at a time */
    const __m128i v_quote = _mm_set1_epi8(quote);

    while (len >= 16) {
        __m128i chunk = _mm_loadu_si128((const __m128i *)src);
        __m128i cmp = _mm_cmpeq_epi8(chunk, v_quote);
        if (_mm_movemask_epi8(cmp)) {
            return true;
        }
        src += 16;
        len -= 16;
    }
#elif CTOOLS_HAS_NEON
    /* NEON: Scan 16 bytes at a time */
    const uint8x16_t v_quote = vdupq_n_u8(quote);

    while (len >= 16) {
        uint8x16_t chunk = vld1q_u8((const uint8_t *)src);
        uint8x16_t cmp = vceqq_u8(chunk, v_quote);
        uint64x2_t cmp64 = vreinterpretq_u64_u8(cmp);
        if (vgetq_lane_u64(cmp64, 0) || vgetq_lane_u64(cmp64, 1)) {
            return true;
        }
        src += 16;
        len -= 16;
    }
#else
    /* Scalar: Unrolled 8-byte check for performance */
    while (len >= 8) {
        if (src[0] == quote || src[1] == quote || src[2] == quote || src[3] == quote ||
            src[4] == quote || src[5] == quote || src[6] == quote || src[7] == quote) {
            return true;
        }
        src += 8;
        len -= 8;
    }
#endif

    /* Scalar tail */
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
    int src_len = CIMPORT_FIELD_LENGTH(*field);
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
    int src_len = CIMPORT_FIELD_LENGTH(*field);

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

bool cimport_field_looks_numeric_sep(const char *src, int len, char dec_sep, char grp_sep)
{
    /* Skip leading whitespace */
    while (len > 0 && (*src == ' ' || *src == '\t')) { src++; len--; }
    if (len == 0) return true;  /* Empty -> treat as missing numeric */

    /* Skip surrounding quotes */
    if (*src == '"' && len >= 2) { src++; len -= 2; }
    if (len == 0) return true;

    /* Skip trailing whitespace */
    while (len > 0 && (src[len-1] == ' ' || src[len-1] == '\t' ||
                       src[len-1] == '\r' || src[len-1] == '\n')) { len--; }
    if (len == 0) return true;

    /* Check for NA/NaN */
    if (len == 2 && (src[0] == 'N' || src[0] == 'n') && (src[1] == 'A' || src[1] == 'a')) return true;
    if (len == 3 && (src[0] == 'N' || src[0] == 'n') && (src[1] == 'a' || src[1] == 'A') &&
        (src[2] == 'N' || src[2] == 'n')) return true;

    /* Check for single decimal separator (Stata missing) */
    if (len == 1 && *src == dec_sep) return true;

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
    if (!ctools_parse_double_with_separators(src, len, &val, SV_missval, dec_sep, grp_sep)) return false;

    *out_value = val;

    /* Check if integer by looking for decimal separator */
    *out_is_integer = true;
    for (int i = 0; i < len; i++) {
        if (src[i] == dec_sep) {
            *out_is_integer = false;
            break;
        }
    }

    return true;
}

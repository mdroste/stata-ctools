/*
 * cimport_parse.h
 * CSV parsing functions for cimport
 *
 * Provides SIMD-accelerated parsing for CSV files including:
 * - Row boundary detection (loose and strict quote modes)
 * - Field extraction with quote handling
 * - Type detection heuristics
 */

#ifndef CIMPORT_PARSE_H
#define CIMPORT_PARSE_H

#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include "../ctools_eisel_lemire.h"
#include "../ctools_config.h"

/* Bindquotes mode for row boundary detection */
typedef enum {
    CIMPORT_BINDQUOTES_LOOSE = 0,   /* Each line is a row, ignore quotes */
    CIMPORT_BINDQUOTES_NOBIND = 2,
    CIMPORT_BINDQUOTES_STRICT = 1   /* Respect quotes: fields can span lines */
} CImportBindQuotesMode;

/* Field reference - offset and length into memory-mapped file.
 * High bit of length serves as quote flag to avoid redundant SIMD scans. */
typedef struct {
    uint64_t offset;
    uint32_t length;   /* Bit 31 = has_quote flag, bits 0-30 = actual length */
} CImportFieldRef;

#define CIMPORT_FIELD_QUOTED_FLAG 0x80000000u
#define CIMPORT_FIELD_LENGTH(f) ((f).length & 0x7FFFFFFFu)

/* Parsed row. Fields are contiguous in the file (one delimiter apart), so
 * the row keeps the first field's offset and each field's length (with the
 * quote flag); offsets are recovered by summing lengths. A uint64_t array of
 * cached numeric values follows the lengths (see cimport_row_values). */
typedef struct {
    uint64_t offset;
    uint16_t num_fields;
    uint32_t lengths[];
} CImportParsedRow;

/* Cached value slot for fields parsed from text at store time (a NaN
 * payload never produced by the number parsers). */
#define CIMPORT_NO_VALUE UINT64_C(0xFFF8C0DE00000001)

static inline size_t cimport_row_values_offset(int num_fields)
{
    return (offsetof(CImportParsedRow, lengths) + sizeof(uint32_t) * (size_t)num_fields + 7) & ~(size_t)7;
}

static inline size_t cimport_row_size(int num_fields)
{
    return cimport_row_values_offset(num_fields) + sizeof(uint64_t) * (size_t)num_fields;
}

static inline uint64_t *cimport_row_values(const CImportParsedRow *row)
{
    return (uint64_t *)((char *)row + cimport_row_values_offset(row->num_fields));
}

/* Copy parser output into a row; false if fields are not contiguous. */
static inline bool cimport_fill_row(CImportParsedRow *row, const CImportFieldRef *fields, int num_fields)
{
    uint64_t next = num_fields ? fields[0].offset : 0;
    row->offset = next;
    row->num_fields = (uint16_t)num_fields;
    for (int i = 0; i < num_fields; i++) {
        if (fields[i].offset != next) return false;
        row->lengths[i] = fields[i].length;
        next += CIMPORT_FIELD_LENGTH(fields[i]) + 1;
    }
    return true;
}

/* Field f of a row (f < num_fields); O(f) because offsets are implicit. */
static inline CImportFieldRef cimport_row_field(const CImportParsedRow *row, int f)
{
    CImportFieldRef field = {row->offset, 0};
    for (int i = 0; i < f; i++) field.offset += (row->lengths[i] & 0x7FFFFFFFu) + 1;
    field.length = row->lengths[f];
    return field;
}

/* ============================================================================
 * SIMD-Accelerated Scanning
 * ============================================================================ */

/* ============================================================================
 * Row Boundary Detection
 * ============================================================================ */

/* Find next row - STRICT mode: respects quotes (fields can span lines) */
const char *cimport_find_next_row_strict(const char *ptr, const char *end, char quote);
/* Same, starting inside a quoted field when in_quotes is true. */
const char *cimport_find_next_row_strict_from(const char *ptr, const char *end,
                                              char quote, bool in_quotes);

/* Wrapper that chooses mode based on bindquotes setting */
const char *cimport_find_next_row(const char *ptr, const char *end, char quote, CImportBindQuotesMode bindquotes);

/* ============================================================================
 * Field Parsing
 * ============================================================================ */

/*
 * Parse a row into field references (fast, zero-copy).
 *
 * @param start              Start of row (or full buffer for fused mode)
 * @param end                End of buffer
 * @param delim              Delimiter character
 * @param quote              Quote character
 * @param fields             Output array for field references
 * @param max_fields         Maximum fields to parse
 * @param file_base          Base pointer for offset calculation
 * @param bindquotes         Quote binding mode (LOOSE: \n always terminates row;
 *                           STRICT: \n inside quotes is part of field)
 * @param next_row_out       If non-NULL, set to start of next row (past \n) or end
 * @param had_unmatched_quote If non-NULL, set to true if row ended inside quotes
 * @return                   Number of fields parsed
 */
int cimport_parse_row_fast(const char *start, const char *end, char delim, char quote,
                           CImportFieldRef *fields, int max_fields, const char *file_base,
                           CImportBindQuotesMode bindquotes,
                           const char **next_row_out, bool *had_unmatched_quote);

/*
 * Extract field value with quote handling.
 * Handles: quoted fields, escaped quotes (""), orphan quotes.
 *
 * @param file_base  Base pointer of memory-mapped file
 * @param field      Field reference
 * @param output     Output buffer
 * @param max_len    Maximum output length
 * @param quote      Quote character
 * @return           Length of extracted string
 */
int cimport_extract_field_fast(const char *file_base, CImportFieldRef *field,
                                char *output, int max_len, char quote);
/* With output == NULL the decoded length is counted without writing. */

/* ============================================================================
 * Type Detection
 * ============================================================================ */

/* Quick check if field looks like a number (for type inference) */
int cimport_configure_locale(const char *name,const char *decimal,const char *group);
/* True while a parselocale() profile replaces the default number grammar. */
bool cimport_locale_active(void);

bool cimport_parse_unquoted_number(const char *src, int len, double *value, double missing,
                                   char decimal, char group);

bool cimport_parse_number(const char *src, int len, double *value, double missing,
                          char decimal, char group);

bool cimport_field_looks_unquoted_numeric_sep(const char *src, int len, char dec_sep, char grp_sep);

bool cimport_field_looks_numeric_sep(const char *src, int len, char dec_sep, char grp_sep);

/*
 * Single-pass parser for plain fields: optional '-', at most 19 digits with at
 * most one '.', optional trailing CR. Valid only for the default grammar
 * ('.' decimals, no group separator, no locale). On success the value is
 * bit-identical to cimport_parse_unquoted_number (both are correctly rounded)
 * and the field passes the numeric test; otherwise callers use those.
 */
static inline bool cimport_parse_plain_decimal(const char *src, int len, double *value)
{
    if (len > 0 && src[len - 1] == '\r') len--;
    const char *p = src, *end = src + len;
    bool negative = p < end && *p == '-';
    if (negative) p++;
    uint64_t mantissa = 0;
    int digits = 0, frac = -1;
    for (; p < end; p++) {
        unsigned d = (unsigned)(unsigned char)*p - '0';
        if (d < 10) {
            mantissa = mantissa * 10 + d;
            digits++;
            if (frac >= 0) frac++;
        } else if (*p == '.' && frac < 0) {
            frac = 0;
        } else {
            return false;
        }
    }
    if (digits == 0 || digits > 19) return false;
    return eisel_lemire_compute(mantissa, frac > 0 ? -frac : 0, negative, value);
}

/*
 * Analyze numeric field and extract value.
 *
 * @param file_base    Base pointer of memory-mapped file
 * @param field        Field reference
 * @param quote        Quote character
 * @param out_value    Output: parsed numeric value
 * @param out_is_integer Output: true if value is integer
 * @return             true if successfully parsed as number
 */
/*
 * Analyze numeric field with custom separators.
 *
 * @param file_base    Base pointer of memory-mapped file
 * @param field        Field reference
 * @param quote        Quote character
 * @param dec_sep      Decimal separator ('.' or ',')
 * @param grp_sep      Group/thousands separator ('\0' = none)
 * @param out_value    Output: parsed numeric value
 * @param out_is_integer Output: true if value is integer
 * @return             true if successfully parsed as number
 */
bool cimport_analyze_numeric_with_sep(const char *file_base, CImportFieldRef *field, char quote,
                                       char dec_sep, char grp_sep,
                                       double *out_value, bool *out_is_integer);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/* Check if character is whitespace (space or tab, not newline) */
bool cimport_is_whitespace(char c);

#endif /* CIMPORT_PARSE_H */

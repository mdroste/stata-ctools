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
#include "../ctools_config.h"

/* Bindquotes mode for row boundary detection */
typedef enum {
    CIMPORT_BINDQUOTES_LOOSE = 0,   /* Each line is a row, ignore quotes */
    CIMPORT_BINDQUOTES_STRICT = 1   /* Respect quotes: fields can span lines */
} CImportBindQuotesMode;

/* Field reference - offset and length into memory-mapped file */
typedef struct {
    uint64_t offset;
    uint32_t length;
} CImportFieldRef;

/* Parsed row - number of fields plus flexible array of field refs */
typedef struct {
    uint16_t num_fields;
    CImportFieldRef fields[];
} CImportParsedRow;

/* ============================================================================
 * SIMD-Accelerated Scanning
 * ============================================================================ */

/*
 * Find next delimiter, newline, or quote character using SIMD.
 * Falls back to scalar loop on platforms without SIMD.
 *
 * @param ptr    Current position in buffer
 * @param end    End of buffer
 * @param delim  Delimiter character (e.g., ',')
 * @return       Pointer to found character or end if not found
 */
const char *cimport_find_delim_or_newline_simd(const char *ptr, const char *end, char delim);

/* ============================================================================
 * Row Boundary Detection
 * ============================================================================ */

/* Find next row - LOOSE mode: each line is a row (ignores quotes) */
const char *cimport_find_next_row_loose(const char *ptr, const char *end);

/* Find next row - STRICT mode: respects quotes (fields can span lines) */
const char *cimport_find_next_row_strict(const char *ptr, const char *end, char quote);

/* Wrapper that chooses mode based on bindquotes setting */
const char *cimport_find_next_row(const char *ptr, const char *end, char quote, CImportBindQuotesMode bindquotes);

/* ============================================================================
 * Field Parsing
 * ============================================================================ */

/*
 * Parse a row into field references (fast, zero-copy).
 *
 * @param start      Start of row
 * @param end        End of row (typically next newline)
 * @param delim      Delimiter character
 * @param quote      Quote character
 * @param fields     Output array for field references
 * @param max_fields Maximum fields to parse
 * @param file_base  Base pointer for offset calculation
 * @return           Number of fields parsed
 */
int cimport_parse_row_fast(const char *start, const char *end, char delim, char quote,
                           CImportFieldRef *fields, int max_fields, const char *file_base);

/* Check if field contains quote character (fast 8-byte unrolled check) */
bool cimport_field_contains_quote(const char *src, int len, char quote);

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

/* Extract unquoted field (fast path, no quote handling) */
int cimport_extract_field_unquoted(const char *file_base, CImportFieldRef *field,
                                    char *output, int max_len);

/* ============================================================================
 * Type Detection
 * ============================================================================ */

/* Quick check if field looks like a number (for type inference) */
bool cimport_field_looks_numeric_sep(const char *src, int len, char dec_sep, char grp_sep);

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

/* Check if row has unmatched quote (odd number of quotes) */
bool cimport_row_has_unmatched_quote(const char *start, const char *end, char quote);

#endif /* CIMPORT_PARSE_H */

/*
 * cimport_xlsx_xml.h
 * Lightweight SAX-style XML parser for XLSX files
 *
 * Memory-efficient parser designed specifically for XLSX worksheet processing.
 * Supports only the subset of XML needed for XLSX parsing.
 */

#ifndef CIMPORT_XLSX_XML_H
#define CIMPORT_XLSX_XML_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

/* Maximum sizes for XML parsing */
#define XLSX_XML_MAX_TAG_LEN      64
#define XLSX_XML_MAX_ATTR_NAME    32
#define XLSX_XML_MAX_ATTR_VALUE   256
#define XLSX_XML_MAX_ATTRIBUTES   16

/* ============================================================================
 * XML Attribute Structure
 * ============================================================================ */

typedef struct {
    char name[XLSX_XML_MAX_ATTR_NAME];
    char value[XLSX_XML_MAX_ATTR_VALUE];
} xlsx_xml_attr;

/* ============================================================================
 * XML Parser Events
 * ============================================================================ */

typedef enum {
    XLSX_XML_START_ELEMENT,   /* Opening tag: <tag ...> */
    XLSX_XML_END_ELEMENT,     /* Closing tag: </tag> or self-closing /> */
    XLSX_XML_TEXT,            /* Text content between tags */
    XLSX_XML_ERROR            /* Parse error */
} xlsx_xml_event_type;

/* Event data passed to callback */
typedef struct {
    xlsx_xml_event_type type;

    /* For START_ELEMENT events */
    const char *tag_name;           /* Tag name (without namespace prefix) */
    const char *tag_name_full;      /* Full tag name (with prefix if any) */
    const xlsx_xml_attr *attrs;     /* Array of attributes */
    int num_attrs;                  /* Number of attributes */
    bool is_self_closing;           /* True if self-closing tag <.../> */

    /* For END_ELEMENT events */
    /* tag_name and tag_name_full are set */

    /* For TEXT events */
    const char *text;               /* Text content (not null-terminated) */
    size_t text_len;                /* Length of text */

    /* For ERROR events */
    const char *error_msg;          /* Error description */
    size_t error_pos;               /* Position in input where error occurred */
} xlsx_xml_event;

/* ============================================================================
 * Callback Function Type
 *
 * Return true to continue parsing, false to stop.
 * ============================================================================ */

typedef bool (*xlsx_xml_callback)(const xlsx_xml_event *event, void *user_data);

/* ============================================================================
 * Parser Context
 * ============================================================================ */

typedef struct xlsx_xml_parser xlsx_xml_parser;

/* Create a new parser with the given callback */
xlsx_xml_parser *xlsx_xml_parser_create(xlsx_xml_callback callback, void *user_data);

/* Destroy a parser */
void xlsx_xml_parser_destroy(xlsx_xml_parser *parser);

/* Parse XML data. Can be called multiple times with chunks of data.
 * Set is_final to true on the last chunk.
 * Returns true on success (or if callback returned false to stop),
 * false on parse error. */
bool xlsx_xml_parse(xlsx_xml_parser *parser, const char *data, size_t len, bool is_final);

/* Get the last error message (valid after parse returns false) */
const char *xlsx_xml_get_error(xlsx_xml_parser *parser);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/* Get attribute value by name from an event, returns NULL if not found */
const char *xlsx_xml_get_attr(const xlsx_xml_event *event, const char *name);

/* Decode XML entities in a string (in-place modification)
 * Handles: &amp; &lt; &gt; &quot; &apos; and numeric refs &#N; &#xN;
 * Returns new length */
size_t xlsx_xml_decode_entities(char *str, size_t len);

/* Parse a cell reference like "A1" into column (0-based) and row (1-based)
 * Returns true on success */
bool xlsx_xml_parse_cell_ref(const char *ref, int *col, int *row);

/* Convert column letters to 0-based index (A=0, B=1, ..., AA=26, etc.) */
int xlsx_xml_col_to_index(const char *col_str);

/* Convert 0-based column index to letters (0=A, 1=B, ..., 26=AA, etc.)
 * Buffer must be at least 4 bytes */
void xlsx_xml_index_to_col(int index, char *buf);

#endif /* CIMPORT_XLSX_XML_H */

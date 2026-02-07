/*
    cexport_xlsx.c
    Excel (.xlsx) export implementation

    Creates XLSX files (ZIP archives containing XML) from Stata data.
    Part of the ctools suite.
*/

#include "cexport_xlsx.h"
#include "../cimport/miniz/miniz.h"
#include "../ctools_threads.h"
#include "../ctools_runtime.h"
#include "../ctools_arena.h"
#include "../ctools_types.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

/* Maximum buffer sizes */
#define MAX_PATH_LEN 4096
#define MAX_SHEET_NAME 31
#define MAX_CELL_CONTENT 32767
#define INITIAL_XML_SIZE (1024 * 1024)  /* 1MB initial XML buffer */

/* Context for XLSX export */
typedef struct {
    char filename[MAX_PATH_LEN];
    char sheet_name[MAX_SHEET_NAME + 1];

    /* Options */
    bool firstrow;
    bool replace;
    bool nolabel;
    bool verbose;
    bool keepcellfmt;

    /* Starting cell (default A1) */
    int start_col;       /* 1-based column number */
    int start_row;       /* 1-based row number */

    /* Missing value replacement (NULL = empty cell) */
    char *missing_value;

    /* Data info */
    ST_int nobs;
    ST_int nvars;
    char **varnames;
    int *vartypes;  /* 0=string, 1-5=numeric types */

    /* Shared strings for deduplication */
    char **shared_strings;
    size_t num_shared_strings;
    size_t shared_strings_capacity;

    /* Preserved styles from existing file (for keepcellfmt) */
    char *preserved_styles;
    size_t preserved_styles_len;

    /* Error message buffer */
    char error_message[512];
} XLSXExportContext;

/* Forward declarations */
static void xlsx_context_init(XLSXExportContext *ctx);
static void xlsx_context_free(XLSXExportContext *ctx);
static ST_retcode xlsx_parse_args(XLSXExportContext *ctx, const char *args);
static ST_retcode xlsx_load_var_metadata(XLSXExportContext *ctx);
static ST_retcode xlsx_write_file(XLSXExportContext *ctx, ctools_filtered_data *filtered);
static void xlsx_escape_xml(const char *src, char *dst, size_t dst_size);
static void xlsx_col_to_letters(int col, char *letters);

/* ============================================================================
 * XML Content Generators
 * ============================================================================ */

static const char *CONTENT_TYPES_XML =
    "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
    "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">\n"
    "  <Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>\n"
    "  <Default Extension=\"xml\" ContentType=\"application/xml\"/>\n"
    "  <Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>\n"
    "  <Override PartName=\"/xl/worksheets/sheet1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>\n"
    "  <Override PartName=\"/xl/styles.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml\"/>\n"
    "  <Override PartName=\"/xl/sharedStrings.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml\"/>\n"
    "</Types>";

static const char *RELS_XML =
    "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
    "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\n"
    "  <Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"xl/workbook.xml\"/>\n"
    "</Relationships>";

static const char *STYLES_XML =
    "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
    "<styleSheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">\n"
    "  <fonts count=\"1\"><font><sz val=\"11\"/><name val=\"Calibri\"/></font></fonts>\n"
    "  <fills count=\"2\"><fill><patternFill patternType=\"none\"/></fill><fill><patternFill patternType=\"gray125\"/></fill></fills>\n"
    "  <borders count=\"1\"><border><left/><right/><top/><bottom/><diagonal/></border></borders>\n"
    "  <cellStyleXfs count=\"1\"><xf numFmtId=\"0\" fontId=\"0\" fillId=\"0\" borderId=\"0\"/></cellStyleXfs>\n"
    "  <cellXfs count=\"1\"><xf numFmtId=\"0\" fontId=\"0\" fillId=\"0\" borderId=\"0\" xfId=\"0\"/></cellXfs>\n"
    "</styleSheet>";

static const char *WORKBOOK_RELS_XML =
    "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
    "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\n"
    "  <Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet1.xml\"/>\n"
    "  <Relationship Id=\"rId2\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles\" Target=\"styles.xml\"/>\n"
    "  <Relationship Id=\"rId3\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/sharedStrings\" Target=\"sharedStrings.xml\"/>\n"
    "</Relationships>";

/* Generate workbook.xml with sheet name */
static char *xlsx_gen_workbook_xml(const char *sheet_name) {
    char *xml = malloc(2048);
    if (!xml) return NULL;

    char escaped_name[MAX_SHEET_NAME * 6 + 1];
    xlsx_escape_xml(sheet_name, escaped_name, sizeof(escaped_name));

    snprintf(xml, 2048,
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
        "<workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
        "xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\">\n"
        "  <sheets>\n"
        "    <sheet name=\"%s\" sheetId=\"1\" r:id=\"rId1\"/>\n"
        "  </sheets>\n"
        "</workbook>",
        escaped_name);

    return xml;
}

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

/* Escape special XML characters */
static void xlsx_escape_xml(const char *src, char *dst, size_t dst_size) {
    size_t di = 0;
    for (size_t si = 0; src[si] && di < dst_size - 6; si++) {
        switch (src[si]) {
            case '&':
                memcpy(dst + di, "&amp;", 5);
                di += 5;
                break;
            case '<':
                memcpy(dst + di, "&lt;", 4);
                di += 4;
                break;
            case '>':
                memcpy(dst + di, "&gt;", 4);
                di += 4;
                break;
            case '"':
                memcpy(dst + di, "&quot;", 6);
                di += 6;
                break;
            case '\'':
                memcpy(dst + di, "&apos;", 6);
                di += 6;
                break;
            default:
                /* Filter control characters except tab, newline, carriage return */
                if ((unsigned char)src[si] >= 32 || src[si] == '\t' ||
                    src[si] == '\n' || src[si] == '\r') {
                    dst[di++] = src[si];
                }
                break;
        }
    }
    dst[di] = '\0';
}

/* Convert column number (1-based) to Excel letters (A, B, ..., Z, AA, AB, ...) */
static void xlsx_col_to_letters(int col, char *letters) {
    char buf[8];
    int i = 0;

    while (col > 0) {
        int rem = (col - 1) % 26;
        buf[i++] = 'A' + rem;
        col = (col - 1) / 26;
    }

    /* Reverse */
    for (int j = 0; j < i; j++) {
        letters[j] = buf[i - 1 - j];
    }
    letters[i] = '\0';
}

/* Parse cell reference (e.g., "B5") into column and row numbers (1-based) */
static void xlsx_parse_cell_ref(const char *cell_ref, int *col, int *row) {
    const char *p = cell_ref;
    *col = 0;
    *row = 0;

    /* Parse letters (column) */
    while (*p && isalpha((unsigned char)*p)) {
        *col = (*col) * 26 + (toupper((unsigned char)*p) - 'A' + 1);
        p++;
    }

    /* Parse digits (row) */
    while (*p && isdigit((unsigned char)*p)) {
        *row = (*row) * 10 + (*p - '0');
        p++;
    }

    /* Default to 1 if not specified */
    if (*col == 0) *col = 1;
    if (*row == 0) *row = 1;
}

/* ============================================================================
 * Context Management
 * ============================================================================ */

static void xlsx_context_init(XLSXExportContext *ctx) {
    memset(ctx, 0, sizeof(*ctx));
    snprintf(ctx->sheet_name, sizeof(ctx->sheet_name), "Sheet1");
    ctx->firstrow = true;  /* Default to writing variable names */
    ctx->start_col = 1;    /* Default to column A */
    ctx->start_row = 1;    /* Default to row 1 */
    ctx->missing_value = NULL;  /* Default to empty cells for missing */
}

static void xlsx_context_free(XLSXExportContext *ctx) {
    if (ctx->varnames) {
        for (ST_int i = 0; i < ctx->nvars; i++) {
            free(ctx->varnames[i]);
        }
        free(ctx->varnames);
    }
    free(ctx->vartypes);
    free(ctx->missing_value);
    free(ctx->preserved_styles);

    if (ctx->shared_strings) {
        for (size_t i = 0; i < ctx->num_shared_strings; i++) {
            free(ctx->shared_strings[i]);
        }
        free(ctx->shared_strings);
    }
}

/* Extract styles.xml from existing XLSX file for keepcellfmt option */
static ST_retcode xlsx_load_preserved_styles(XLSXExportContext *ctx) {
    mz_zip_archive zip;
    memset(&zip, 0, sizeof(zip));

    /* Try to open existing file */
    if (!mz_zip_reader_init_file(&zip, ctx->filename, 0)) {
        /* File doesn't exist or can't be read - not an error, just no styles to preserve */
        return 0;
    }

    /* Find styles.xml in archive */
    int file_index = mz_zip_reader_locate_file(&zip, "xl/styles.xml", NULL, 0);
    if (file_index < 0) {
        mz_zip_reader_end(&zip);
        return 0;  /* No styles to preserve */
    }

    /* Get file info */
    mz_zip_archive_file_stat file_stat;
    if (!mz_zip_reader_file_stat(&zip, file_index, &file_stat)) {
        mz_zip_reader_end(&zip);
        return 0;
    }

    /* Allocate buffer for styles */
    ctx->preserved_styles = malloc(file_stat.m_uncomp_size + 1);
    if (!ctx->preserved_styles) {
        mz_zip_reader_end(&zip);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Memory allocation failed for preserved styles\n");
        return 920;
    }

    /* Extract styles.xml */
    if (!mz_zip_reader_extract_to_mem(&zip, file_index, ctx->preserved_styles,
                                       file_stat.m_uncomp_size, 0)) {
        free(ctx->preserved_styles);
        ctx->preserved_styles = NULL;
        mz_zip_reader_end(&zip);
        return 0;  /* Failed to extract - fall back to default styles */
    }

    ctx->preserved_styles[file_stat.m_uncomp_size] = '\0';
    ctx->preserved_styles_len = file_stat.m_uncomp_size;

    mz_zip_reader_end(&zip);

    if (ctx->verbose) {
        char msg[128];
        snprintf(msg, sizeof(msg), "Preserved styles from existing file (%zu bytes)\n",
                ctx->preserved_styles_len);
        SF_display(msg);
    }

    return 0;
}

/* ============================================================================
 * Argument Parsing
 * ============================================================================ */

static ST_retcode xlsx_parse_args(XLSXExportContext *ctx, const char *args) {
    if (!args || !*args) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "cexport xlsx: missing arguments\n");
        return 198;
    }

    char args_copy[4096];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    /* Parse filename (first token) */
    char *filename = strtok(args_copy, " \t");
    if (!filename) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "cexport xlsx: missing filename\n");
        return 198;
    }
    strncpy(ctx->filename, filename, sizeof(ctx->filename) - 1);

    /* Parse options */
    char *opt;
    while ((opt = strtok(NULL, " \t")) != NULL) {
        if (strncmp(opt, "sheet=", 6) == 0) {
            strncpy(ctx->sheet_name, opt + 6, MAX_SHEET_NAME);
            ctx->sheet_name[MAX_SHEET_NAME] = '\0';
        } else if (strcmp(opt, "firstrow") == 0) {
            ctx->firstrow = true;
        } else if (strcmp(opt, "nofirstrow") == 0) {
            ctx->firstrow = false;
        } else if (strcmp(opt, "replace") == 0) {
            ctx->replace = true;
        } else if (strcmp(opt, "nolabel") == 0) {
            ctx->nolabel = true;
        } else if (strcmp(opt, "verbose") == 0) {
            ctx->verbose = true;
        } else if (strncmp(opt, "cell=", 5) == 0) {
            xlsx_parse_cell_ref(opt + 5, &ctx->start_col, &ctx->start_row);
        } else if (strncmp(opt, "missing=", 8) == 0) {
            ctx->missing_value = strdup(opt + 8);
            if (ctx->missing_value == NULL) {
                return 920;  /* Memory allocation failed */
            }
        } else if (strcmp(opt, "keepcellfmt") == 0) {
            ctx->keepcellfmt = true;
        }
    }

    return 0;
}

/* ============================================================================
 * Metadata Loading
 * ============================================================================ */

static ST_retcode xlsx_load_var_metadata(XLSXExportContext *ctx) {
    /* Get number of observations and variables from Stata.
     * The ado passes `if `touse'' so touse is not in the varlist;
     * SF_nvars() returns the actual data variable count. */
    ctx->nobs = SF_nobs();
    ctx->nvars = SF_nvars();

    if (ctx->nvars <= 0) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "cexport xlsx: no variables to export\n");
        return 198;
    }

    /* Allocate arrays */
    ctx->varnames = calloc(ctx->nvars, sizeof(char *));
    ctx->vartypes = calloc(ctx->nvars, sizeof(int));

    if (!ctx->varnames || !ctx->vartypes) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "cexport xlsx: memory allocation failed\n");
        return 920;
    }

    /* Load variable names from global macro */
    char varnames_buf[32000];
    ST_retcode rc = SF_macro_use("CEXPORT_VARNAMES", varnames_buf, sizeof(varnames_buf));
    if (rc) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "cexport xlsx: failed to get variable names\n");
        return rc;
    }

    /* Parse variable names */
    char *name = strtok(varnames_buf, " ");
    for (ST_int i = 0; i < ctx->nvars && name; i++) {
        ctx->varnames[i] = strdup(name);
        if (ctx->varnames[i] == NULL) {
            snprintf(ctx->error_message, sizeof(ctx->error_message),
                     "cexport xlsx: memory allocation failed for variable name\n");
            return 920;
        }
        name = strtok(NULL, " ");
    }

    /* Load variable types from global macro */
    char vartypes_buf[8000];
    rc = SF_macro_use("CEXPORT_VARTYPES", vartypes_buf, sizeof(vartypes_buf));
    if (rc) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "cexport xlsx: failed to get variable types\n");
        return rc;
    }

    /* Parse variable types */
    char *type_str = strtok(vartypes_buf, " ");
    for (ST_int i = 0; i < ctx->nvars && type_str; i++) {
        ctx->vartypes[i] = atoi(type_str);
        type_str = strtok(NULL, " ");
    }

    return 0;
}

/* ============================================================================
 * Worksheet XML Generation
 * ============================================================================ */

/* Dynamic string buffer for XML generation */
typedef struct {
    char *data;
    size_t len;
    size_t capacity;
} StringBuffer;

static bool strbuf_init(StringBuffer *sb, size_t initial_cap) {
    sb->data = malloc(initial_cap);
    sb->len = 0;
    sb->capacity = sb->data ? initial_cap : 0;
    if (sb->data) sb->data[0] = '\0';
    return sb->data != NULL;
}

static void strbuf_append(StringBuffer *sb, const char *str) {
    if (!sb->data) return;
    size_t slen = strlen(str);
    if (sb->len + slen + 1 > sb->capacity) {
        size_t new_cap = sb->capacity * 2;
        if (new_cap < sb->len + slen + 1) new_cap = sb->len + slen + 1024;
        char *new_data = realloc(sb->data, new_cap);
        if (!new_data) return;
        sb->data = new_data;
        sb->capacity = new_cap;
    }
    memcpy(sb->data + sb->len, str, slen + 1);
    sb->len += slen;
}

static void strbuf_appendf(StringBuffer *sb, const char *fmt, ...) {
    char buf[4096];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    strbuf_append(sb, buf);
}

/* Ensure buffer has at least 'need' bytes available beyond current length */
static void strbuf_ensure(StringBuffer *sb, size_t need) {
    if (!sb->data) return;
    if (sb->len + need + 1 > sb->capacity) {
        size_t new_cap = sb->capacity * 2;
        if (new_cap < sb->len + need + 1) new_cap = sb->len + need + 1024;
        char *new_data = realloc(sb->data, new_cap);
        if (!new_data) return;
        sb->data = new_data;
        sb->capacity = new_cap;
    }
}

/* Append string with known length (avoids strlen) */
static void strbuf_append_n(StringBuffer *sb, const char *str, size_t slen) {
    if (!sb->data) return;
    strbuf_ensure(sb, slen);
    memcpy(sb->data + sb->len, str, slen);
    sb->len += slen;
    sb->data[sb->len] = '\0';
}

/* Write int64 directly into buffer (avoids snprintf) */
static void strbuf_append_int64(StringBuffer *sb, long long val) {
    char buf[24];
    int len = snprintf(buf, sizeof(buf), "%lld", val);
    if (len > 0) strbuf_append_n(sb, buf, (size_t)len);
}

/* Write double directly into buffer with specified format */
static void strbuf_append_double(StringBuffer *sb, double val, int precision) {
    char buf[64];
    int len;
    if (precision == 0) {
        len = snprintf(buf, sizeof(buf), "%.0f", val);
    } else {
        len = snprintf(buf, sizeof(buf), "%.15g", val);
    }
    if (len > 0) strbuf_append_n(sb, buf, (size_t)len);
}

/* Generate worksheet XML with data from bulk-loaded in-memory arrays.
 * filtered->data.vars[j] contains data for variable j (0-based, data vars only).
 * filtered->data.nobs is the number of filtered observations. */
static char *xlsx_gen_worksheet_xml(XLSXExportContext *ctx,
                                     ctools_filtered_data *filtered,
                                     size_t *out_len) {
    StringBuffer sb;
    if (!strbuf_init(&sb, INITIAL_XML_SIZE)) {
        *out_len = 0;
        return NULL;
    }

    size_t nobs_filtered = filtered->data.nobs;

    /* XML header */
    strbuf_append(&sb, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n");
    strbuf_append(&sb, "<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">\n");

    /* Dimension - use start cell offset */
    char start_col_letters[8], end_col_letters[8];
    xlsx_col_to_letters(ctx->start_col, start_col_letters);
    xlsx_col_to_letters(ctx->start_col + ctx->nvars - 1, end_col_letters);
    long long total_rows = (long long)nobs_filtered + (ctx->firstrow ? 1 : 0);
    long long end_row = (long long)ctx->start_row + total_rows - 1;
    strbuf_appendf(&sb, "  <dimension ref=\"%s%d:%s%lld\"/>\n",
                  start_col_letters, ctx->start_row, end_col_letters, end_row);

    /* Sheet data */
    strbuf_append(&sb, "  <sheetData>\n");

    long long row_num = ctx->start_row;  /* Start at specified row */

    /* Header row if firstrow option */
    if (ctx->firstrow) {
        strbuf_appendf(&sb, "    <row r=\"%lld\">\n", row_num);
        for (ST_int j = 0; j < ctx->nvars; j++) {
            char col_letters[8];
            xlsx_col_to_letters(ctx->start_col + j, col_letters);

            char escaped_name[256];
            xlsx_escape_xml(ctx->varnames[j], escaped_name, sizeof(escaped_name));

            strbuf_appendf(&sb, "      <c r=\"%s%lld\" t=\"inlineStr\"><is><t>%s</t></is></c>\n",
                          col_letters, row_num, escaped_name);
        }
        strbuf_append(&sb, "    </row>\n");
        row_num++;
    }

    /* Data rows - read from in-memory arrays instead of SPI calls */
    char escaped[MAX_CELL_CONTENT * 6 + 1];

    /* Prepare escaped missing value if specified */
    char escaped_missing[MAX_CELL_CONTENT * 6 + 1];
    if (ctx->missing_value) {
        xlsx_escape_xml(ctx->missing_value, escaped_missing, sizeof(escaped_missing));
    }

    /* Pre-compute column letters for all variables (avoids per-row recomputation) */
    char (*col_letters_arr)[8] = (char (*)[8])malloc(ctx->nvars * 8);
    size_t *col_letters_len = (size_t *)malloc(ctx->nvars * sizeof(size_t));
    if (col_letters_arr && col_letters_len) {
        for (ST_int j = 0; j < ctx->nvars; j++) {
            xlsx_col_to_letters(ctx->start_col + j, col_letters_arr[j]);
            col_letters_len[j] = strlen(col_letters_arr[j]);
        }
    }

    /* Pre-compute escaped missing length */
    size_t escaped_missing_len = ctx->missing_value ? strlen(escaped_missing) : 0;

    for (size_t i = 0; i < nobs_filtered; i++) {
        /* "    <row r=\"" + int + "\">\n" */
        strbuf_append_n(&sb, "    <row r=\"", 12);
        strbuf_append_int64(&sb, row_num);
        strbuf_append_n(&sb, "\">\n", 3);

        for (ST_int j = 0; j < ctx->nvars; j++) {
            const char *col_l = col_letters_arr ? col_letters_arr[j] : "A";
            size_t col_l_len = col_letters_arr ? col_letters_len[j] : 1;

            int vartype = ctx->vartypes[j];

            if (vartype == 0) {
                /* String variable */
                const char *str = filtered->data.vars[j].data.str[i];
                if (!str || str[0] == '\0') {
                    if (ctx->missing_value) {
                        /* "      <c r=\"" + col + row + "\" t=\"inlineStr\"><is><t>" + missing + "</t></is></c>\n" */
                        strbuf_append_n(&sb, "      <c r=\"", 12);
                        strbuf_append_n(&sb, col_l, col_l_len);
                        strbuf_append_int64(&sb, row_num);
                        strbuf_append_n(&sb, "\" t=\"inlineStr\"><is><t>", 23);
                        strbuf_append_n(&sb, escaped_missing, escaped_missing_len);
                        strbuf_append_n(&sb, "</t></is></c>\n", 14);
                    } else {
                        strbuf_append_n(&sb, "      <c r=\"", 12);
                        strbuf_append_n(&sb, col_l, col_l_len);
                        strbuf_append_int64(&sb, row_num);
                        strbuf_append_n(&sb, "\"/>\n", 4);
                    }
                } else {
                    xlsx_escape_xml(str, escaped, sizeof(escaped));
                    size_t esc_len = strlen(escaped);
                    strbuf_append_n(&sb, "      <c r=\"", 12);
                    strbuf_append_n(&sb, col_l, col_l_len);
                    strbuf_append_int64(&sb, row_num);
                    strbuf_append_n(&sb, "\" t=\"inlineStr\"><is><t>", 23);
                    strbuf_append_n(&sb, escaped, esc_len);
                    strbuf_append_n(&sb, "</t></is></c>\n", 14);
                }
            } else {
                /* Numeric variable */
                double val = filtered->data.vars[j].data.dbl[i];

                if (SF_is_missing(val)) {
                    if (ctx->missing_value) {
                        strbuf_append_n(&sb, "      <c r=\"", 12);
                        strbuf_append_n(&sb, col_l, col_l_len);
                        strbuf_append_int64(&sb, row_num);
                        strbuf_append_n(&sb, "\" t=\"inlineStr\"><is><t>", 23);
                        strbuf_append_n(&sb, escaped_missing, escaped_missing_len);
                        strbuf_append_n(&sb, "</t></is></c>\n", 14);
                    } else {
                        strbuf_append_n(&sb, "      <c r=\"", 12);
                        strbuf_append_n(&sb, col_l, col_l_len);
                        strbuf_append_int64(&sb, row_num);
                        strbuf_append_n(&sb, "\"/>\n", 4);
                    }
                } else {
                    strbuf_append_n(&sb, "      <c r=\"", 12);
                    strbuf_append_n(&sb, col_l, col_l_len);
                    strbuf_append_int64(&sb, row_num);
                    if (val == floor(val) && fabs(val) < 1e15) {
                        strbuf_append_n(&sb, "\"><v>", 5);
                        strbuf_append_double(&sb, val, 0);
                        strbuf_append_n(&sb, "</v></c>\n", 9);
                    } else {
                        strbuf_append_n(&sb, "\"><v>", 5);
                        strbuf_append_double(&sb, val, 15);
                        strbuf_append_n(&sb, "</v></c>\n", 9);
                    }
                }
            }
        }

        strbuf_append_n(&sb, "    </row>\n", 11);
        row_num++;
    }

    free(col_letters_arr);
    free(col_letters_len);

    strbuf_append(&sb, "  </sheetData>\n");
    strbuf_append(&sb, "</worksheet>");

    *out_len = sb.len;
    return sb.data;
}

/* Generate shared strings XML (empty for now, using inline strings) */
static const char *SHARED_STRINGS_XML =
    "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
    "<sst xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" count=\"0\" uniqueCount=\"0\">\n"
    "</sst>";

/* ============================================================================
 * XLSX File Writing
 * ============================================================================ */

static ST_retcode xlsx_write_file(XLSXExportContext *ctx, ctools_filtered_data *filtered) {
    mz_zip_archive zip;
    memset(&zip, 0, sizeof(zip));

    /* Initialize ZIP writer */
    if (!mz_zip_writer_init_file(&zip, ctx->filename, 0)) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to create XLSX file: %s\n", ctx->filename);
        return 603;
    }

    /* Add [Content_Types].xml */
    if (!mz_zip_writer_add_mem(&zip, "[Content_Types].xml",
                               CONTENT_TYPES_XML, strlen(CONTENT_TYPES_XML),
                               MZ_DEFAULT_COMPRESSION)) {
        mz_zip_writer_end(&zip);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to write [Content_Types].xml\n");
        return 603;
    }

    /* Add _rels/.rels */
    if (!mz_zip_writer_add_mem(&zip, "_rels/.rels",
                               RELS_XML, strlen(RELS_XML),
                               MZ_DEFAULT_COMPRESSION)) {
        mz_zip_writer_end(&zip);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to write _rels/.rels\n");
        return 603;
    }

    /* Add xl/workbook.xml */
    char *workbook_xml = xlsx_gen_workbook_xml(ctx->sheet_name);
    if (!workbook_xml) {
        mz_zip_writer_end(&zip);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Memory allocation failed for workbook.xml\n");
        return 920;
    }

    if (!mz_zip_writer_add_mem(&zip, "xl/workbook.xml",
                               workbook_xml, strlen(workbook_xml),
                               MZ_DEFAULT_COMPRESSION)) {
        free(workbook_xml);
        mz_zip_writer_end(&zip);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to write xl/workbook.xml\n");
        return 603;
    }
    free(workbook_xml);

    /* Add xl/_rels/workbook.xml.rels */
    if (!mz_zip_writer_add_mem(&zip, "xl/_rels/workbook.xml.rels",
                               WORKBOOK_RELS_XML, strlen(WORKBOOK_RELS_XML),
                               MZ_DEFAULT_COMPRESSION)) {
        mz_zip_writer_end(&zip);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to write xl/_rels/workbook.xml.rels\n");
        return 603;
    }

    /* Add xl/styles.xml - use preserved styles if keepcellfmt is enabled */
    const char *styles_to_write = STYLES_XML;
    size_t styles_len = strlen(STYLES_XML);
    if (ctx->keepcellfmt && ctx->preserved_styles) {
        styles_to_write = ctx->preserved_styles;
        styles_len = ctx->preserved_styles_len;
    }
    if (!mz_zip_writer_add_mem(&zip, "xl/styles.xml",
                               styles_to_write, styles_len,
                               MZ_DEFAULT_COMPRESSION)) {
        mz_zip_writer_end(&zip);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to write xl/styles.xml\n");
        return 603;
    }

    /* Add xl/sharedStrings.xml */
    if (!mz_zip_writer_add_mem(&zip, "xl/sharedStrings.xml",
                               SHARED_STRINGS_XML, strlen(SHARED_STRINGS_XML),
                               MZ_DEFAULT_COMPRESSION)) {
        mz_zip_writer_end(&zip);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to write xl/sharedStrings.xml\n");
        return 603;
    }

    /* Generate and add worksheet */
    size_t sheet_len;
    char *sheet_xml = xlsx_gen_worksheet_xml(ctx, filtered, &sheet_len);
    if (!sheet_xml) {
        mz_zip_writer_end(&zip);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Memory allocation failed for worksheet\n");
        return 920;
    }

    if (!mz_zip_writer_add_mem(&zip, "xl/worksheets/sheet1.xml",
                               sheet_xml, sheet_len,
                               MZ_DEFAULT_COMPRESSION)) {
        free(sheet_xml);
        mz_zip_writer_end(&zip);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to write xl/worksheets/sheet1.xml\n");
        return 603;
    }
    free(sheet_xml);

    /* Finalize ZIP */
    if (!mz_zip_writer_finalize_archive(&zip)) {
        mz_zip_writer_end(&zip);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to finalize XLSX archive\n");
        return 603;
    }

    mz_zip_writer_end(&zip);
    return 0;
}

/* ============================================================================
 * Main Entry Point
 * ============================================================================ */

ST_retcode cexport_xlsx_main(const char *args) {
    XLSXExportContext ctx;
    xlsx_context_init(&ctx);

    double t_start = ctools_timer_seconds();

    /* Parse arguments */
    ST_retcode rc = xlsx_parse_args(&ctx, args);
    if (rc) {
        SF_error(ctx.error_message);
        xlsx_context_free(&ctx);
        return rc;
    }

    /* Check if file exists (unless replace) */
    if (!ctx.replace) {
        FILE *f = fopen(ctx.filename, "r");
        if (f) {
            fclose(f);
            char err_buf[512];
            snprintf(err_buf, sizeof(err_buf),
                     "file %s already exists; use replace option\n", ctx.filename);
            SF_error(err_buf);
            xlsx_context_free(&ctx);
            return 602;
        }
    }

    /* Load preserved styles from existing file if keepcellfmt is enabled */
    if (ctx.keepcellfmt) {
        rc = xlsx_load_preserved_styles(&ctx);
        if (rc) {
            SF_error(ctx.error_message);
            xlsx_context_free(&ctx);
            return rc;
        }
    }

    /* Load variable metadata */
    rc = xlsx_load_var_metadata(&ctx);
    if (rc) {
        SF_error(ctx.error_message);
        xlsx_context_free(&ctx);
        return rc;
    }

    double t_meta = ctools_timer_seconds();

    /* Bulk load data from Stata into memory (parallel).
     * The ado passes `if `touse'' to the plugin call, so SF_ifobs()
     * filters correctly. We only need to load data variables (1..nvars). */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);

    {
        int *var_indices = (int *)malloc(ctx.nvars * sizeof(int));
        if (!var_indices) {
            SF_error("cexport xlsx: memory allocation failed for var indices\n");
            xlsx_context_free(&ctx);
            return 920;
        }
        for (ST_int j = 0; j < ctx.nvars; j++) {
            var_indices[j] = j + 1;
        }

        rc = ctools_data_load(&filtered, var_indices, ctx.nvars, 0, 0, CTOOLS_LOAD_CHECK_IF);
        free(var_indices);

        if (rc != 0) {
            SF_error("cexport xlsx: failed to bulk load data\n");
            ctools_filtered_data_free(&filtered);
            xlsx_context_free(&ctx);
            return 920;
        }
    }

    double t_load = ctools_timer_seconds();

    /* Write XLSX file using in-memory data */
    rc = xlsx_write_file(&ctx, &filtered);
    ctools_filtered_data_free(&filtered);
    if (rc) {
        SF_error(ctx.error_message);
        xlsx_context_free(&ctx);
        return rc;
    }

    double t_end = ctools_timer_seconds();

    /* Report timing if verbose */
    if (ctx.verbose) {
        char msg[256];
        snprintf(msg, sizeof(msg), "Metadata load: %.3f sec\n", t_meta - t_start);
        SF_display(msg);
        snprintf(msg, sizeof(msg), "Data load: %.3f sec\n", t_load - t_meta);
        SF_display(msg);
        snprintf(msg, sizeof(msg), "Write XLSX: %.3f sec\n", t_end - t_load);
        SF_display(msg);
        snprintf(msg, sizeof(msg), "Total: %.3f sec\n", t_end - t_start);
        SF_display(msg);
    }

    /* Save timing to macros for .ado */
    char macro_val[32];
    snprintf(macro_val, sizeof(macro_val), "%.6f", t_end - t_start);
    SF_macro_save("_cexport_xlsx_time", macro_val);

    xlsx_context_free(&ctx);
    return 0;
}

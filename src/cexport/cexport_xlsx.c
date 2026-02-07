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
#include "../ctools_config.h"

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
static size_t xlsx_escape_xml(const char *src, char *dst, size_t dst_size);
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

/* Escape special XML characters. Returns output length. */
static size_t xlsx_escape_xml(const char *src, char *dst, size_t dst_size) {
    /* Fast path: scan for characters that need escaping.
     * Most strings contain no special chars — do a single memcpy. */
    size_t src_len = strlen(src);
    if (src_len > 0 && src_len < dst_size) {
        bool needs_escape = false;
        for (size_t i = 0; i < src_len; i++) {
            unsigned char ch = (unsigned char)src[i];
            if (ch == '&' || ch == '<' || ch == '>' || ch == '"' || ch == '\'' || ch < 32) {
                needs_escape = true;
                break;
            }
        }
        if (!needs_escape) {
            memcpy(dst, src, src_len);
            dst[src_len] = '\0';
            return src_len;
        }
    }

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
    return di;
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

/* ============================================================================
 * Streaming ZIP Segments
 *
 * Instead of concatenating preamble + chunks + footer into one huge buffer,
 * we keep them as separate segments and stream through them via a read callback.
 * This eliminates ~500MB peak allocation for 10M×10 datasets.
 * ============================================================================ */

typedef struct {
    const char *data;
    size_t len;
} xlsx_segment;

typedef struct {
    xlsx_segment *segments;
    size_t num_segments;
    size_t total_size;
    /* Ownership tracking for cleanup */
    StringBuffer preamble_buf;
    StringBuffer *chunk_bufs;  /* Array of chunk StringBuffers (num_segments - 2) */
    size_t num_chunks;
    char (*col_letters_arr)[8];
    size_t *col_letters_len;
    char (*col_prefix)[24];
    size_t *col_prefix_len;
} xlsx_worksheet_segments;

static void xlsx_segments_free(xlsx_worksheet_segments *segs) {
    if (!segs) return;
    if (segs->chunk_bufs) {
        for (size_t i = 0; i < segs->num_chunks; i++) {
            free(segs->chunk_bufs[i].data);
        }
        free(segs->chunk_bufs);
    }
    free(segs->preamble_buf.data);
    free(segs->segments);
    free(segs->col_letters_arr);
    free(segs->col_letters_len);
    free(segs->col_prefix);
    free(segs->col_prefix_len);
}

/* Read callback for mz_zip_writer_add_read_buf_callback.
 * Maps file_ofs to the correct segment and copies data. */
static size_t xlsx_segment_read_callback(void *opaque, mz_uint64 file_ofs, void *pBuf, size_t n)
{
    xlsx_worksheet_segments *segs = (xlsx_worksheet_segments *)opaque;
    if (file_ofs >= segs->total_size) return 0;

    size_t remaining = n;
    size_t buf_pos = 0;
    char *out = (char *)pBuf;

    /* Find starting segment */
    mz_uint64 seg_start = 0;
    for (size_t i = 0; i < segs->num_segments; i++) {
        mz_uint64 seg_end = seg_start + segs->segments[i].len;
        if (file_ofs < seg_end) {
            /* Start reading from this segment */
            size_t offset_in_seg = (size_t)(file_ofs - seg_start);
            size_t avail = segs->segments[i].len - offset_in_seg;
            size_t to_copy = remaining < avail ? remaining : avail;
            memcpy(out + buf_pos, segs->segments[i].data + offset_in_seg, to_copy);
            buf_pos += to_copy;
            remaining -= to_copy;
            file_ofs += to_copy;

            /* Continue to subsequent segments if needed */
            for (size_t j = i + 1; j < segs->num_segments && remaining > 0; j++) {
                avail = segs->segments[j].len;
                to_copy = remaining < avail ? remaining : avail;
                memcpy(out + buf_pos, segs->segments[j].data, to_copy);
                buf_pos += to_copy;
                remaining -= to_copy;
            }
            break;
        }
        seg_start = seg_end;
    }

    return buf_pos;
}

/* ============================================================================
 * Fast Double-to-String Conversion
 * ============================================================================ */

/*
 * Fast double-to-string converter for XLSX numeric cells.
 * Produces output equivalent to snprintf("%.15g", val) but ~5-10x faster.
 * Returns number of characters written (not including null terminator).
 *
 * Strategy:
 * - Fixed-point path for values in [1e-4, 1e15)
 * - Scientific notation path for everything else
 * - Uses CTOOLS_DIGIT_PAIRS for fast 2-digit output
 * - Strips trailing zeros from fractional part
 *
 * Round-trip guarantee: strtod(output) == original value (15 significant digits).
 */
static int xlsx_fast_dtoa(double val, char *buf)
{
    char *start = buf;
    double original_val = val;

    /* Handle sign */
    if (val < 0) {
        *buf++ = '-';
        val = -val;
    }

    /* Handle zero */
    if (val == 0.0) {
        *buf++ = '0';
        *buf = '\0';
        return (int)(buf - start);
    }

    double abs_val = val;

    /* Choose path based on magnitude */
    if (abs_val >= 1e-4 && abs_val < 1e15) {
        /* Fixed-point path */
        int64_t int_part = (int64_t)abs_val;
        double frac = abs_val - (double)int_part;

        /* Write integer part */
        buf += ctools_int64_to_str(int_part, buf);

        /* Compute how many significant digits are left for the fractional part.
         * We want 15 total significant digits. */
        int int_digits = 0;
        {
            int64_t tmp = int_part;
            if (tmp == 0) {
                int_digits = 1;
            } else {
                while (tmp > 0) { int_digits++; tmp /= 10; }
            }
        }
        int frac_digits = 15 - int_digits;
        if (frac_digits < 0) frac_digits = 0;
        if (frac_digits > 15) frac_digits = 15;

        if (frac > 0.0 && frac_digits > 0) {
            /* Multiply fractional part to extract digits */
            *buf++ = '.';
            char *frac_start = buf;

            /* Scale frac to an integer with frac_digits digits */
            double scale = 1.0;
            for (int i = 0; i < frac_digits; i++) scale *= 10.0;
            int64_t frac_int = (int64_t)(frac * scale + 0.5);

            /* Handle rounding overflow */
            if (frac_int >= (int64_t)scale) {
                frac_int = (int64_t)scale - 1;
            }

            /* Write fractional digits with leading zeros */
            char frac_buf[20];
            int flen = 0;
            int64_t tmp = frac_int;
            if (tmp == 0) {
                frac_buf[flen++] = '0';
            } else {
                while (tmp > 0) {
                    frac_buf[flen++] = '0' + (char)(tmp % 10);
                    tmp /= 10;
                }
            }
            /* Pad leading zeros */
            for (int i = flen; i < frac_digits; i++) {
                *buf++ = '0';
            }
            /* Write digits in reverse */
            for (int i = flen - 1; i >= 0; i--) {
                *buf++ = frac_buf[i];
            }

            /* Strip trailing zeros */
            while (buf > frac_start && *(buf - 1) == '0') buf--;
            if (buf > frac_start && *(buf - 1) == '.') buf--;
        }

        *buf = '\0';

        /* Verify round-trip accuracy; fall back to snprintf if needed */
        {
            double check = strtod(start, NULL);
            if (check != original_val) {
                return snprintf(start, 64, "%.15g", original_val);
            }
        }

        return (int)(buf - start);
    } else {
        /* Scientific notation path — normalize to 1 <= m < 10 */
        int exponent = 0;
        double m = abs_val;

        if (m >= 10.0) {
            while (m >= 1e8) { m *= 1e-8; exponent += 8; }
            while (m >= 10.0) { m *= 0.1; exponent++; }
        } else if (m < 1.0 && m > 0.0) {
            while (m < 1e-7) { m *= 1e8; exponent -= 8; }
            while (m < 1.0) { m *= 10.0; exponent--; }
        }

        /* Extract 15 significant digits from mantissa */
        int64_t mantissa = (int64_t)(m * 1e14 + 0.5);
        if (mantissa >= (int64_t)1e15) {
            mantissa /= 10;
            exponent++;
        }

        /* Write mantissa digits */
        char mant_buf[20];
        int mlen = 0;
        {
            int64_t tmp = mantissa;
            while (tmp > 0) {
                mant_buf[mlen++] = '0' + (char)(tmp % 10);
                tmp /= 10;
            }
        }
        /* Pad to 15 digits */
        while (mlen < 15) mant_buf[mlen++] = '0';

        /* First digit */
        *buf++ = mant_buf[mlen - 1];

        /* Strip trailing zeros from remaining mantissa digits */
        int last_nonzero = 0;
        for (int i = 0; i < mlen - 1; i++) {
            if (mant_buf[i] != '0') last_nonzero = i;
        }

        if (last_nonzero > 0 || mant_buf[0] != '0') {
            *buf++ = '.';
            for (int i = mlen - 2; i >= (mlen - 1 - last_nonzero - 1) && i >= 0; i--) {
                *buf++ = mant_buf[i];
            }
            /* Trim trailing zeros after decimal point */
            while (*(buf - 1) == '0') buf--;
            if (*(buf - 1) == '.') buf--;
        }

        /* Write exponent */
        *buf++ = 'e';
        if (exponent >= 0) {
            *buf++ = '+';
        } else {
            *buf++ = '-';
            exponent = -exponent;
        }

        if (exponent >= 100) {
            *buf++ = '0' + (char)(exponent / 100);
            exponent %= 100;
            buf[0] = CTOOLS_DIGIT_PAIRS[exponent * 2];
            buf[1] = CTOOLS_DIGIT_PAIRS[exponent * 2 + 1];
            buf += 2;
        } else if (exponent >= 10) {
            buf[0] = CTOOLS_DIGIT_PAIRS[exponent * 2];
            buf[1] = CTOOLS_DIGIT_PAIRS[exponent * 2 + 1];
            buf += 2;
        } else {
            *buf++ = '0' + (char)exponent;
        }

        *buf = '\0';

        /* Verify round-trip */
        {
            double check = strtod(start, NULL);
            if (check != original_val) {
                return snprintf(start, 64, "%.15g", original_val);
            }
        }

        return (int)(buf - start);
    }
}

/* ============================================================================
 * Parallel Worksheet Generation
 * ============================================================================ */

/* Thread arguments for parallel XML chunk generation */
typedef struct {
    const XLSXExportContext *ctx;
    const ctools_filtered_data *filtered;
    size_t start_idx;           /* 0-based index into filtered data */
    size_t end_idx;             /* exclusive */
    long long start_row_num;    /* 1-based Excel row number */
    char (*col_letters_arr)[8]; /* shared pre-computed column letters */
    size_t *col_letters_len;    /* shared pre-computed column letter lengths */
    char (*col_prefix)[24];     /* pre-computed "      <c r=\"" + col letters */
    size_t *col_prefix_len;     /* length of each prefix */
    const char *escaped_missing;
    size_t escaped_missing_len;
    StringBuffer sb;            /* output buffer */
} xlsx_chunk_args_t;

/* Format a chunk of data rows into XML (called per thread) */
static void *xlsx_format_chunk(void *arg) {
    xlsx_chunk_args_t *a = (xlsx_chunk_args_t *)arg;
    StringBuffer *sb = &a->sb;
    const XLSXExportContext *ctx = a->ctx;
    const ctools_filtered_data *filtered = a->filtered;
    long long row_num = a->start_row_num;

    /* Per-thread escape buffer for strings (heap-allocated for thread safety) */
    char *escaped = malloc(MAX_CELL_CONTENT * 6 + 1);
    if (!escaped) return (void *)(intptr_t)-1;

    /* Per-cell buffer for building cell XML before single append */
    char cell_buf[128];

    /* Per-row ensure budget: row tags + per-cell numeric max */
    size_t row_budget = 64 + (size_t)ctx->nvars * 80;

    for (size_t i = a->start_idx; i < a->end_idx; i++) {
        /* Pre-ensure for this row's content (avoids per-cell capacity checks) */
        strbuf_ensure(sb, row_budget);

        /* Row open tag */
        int rlen = 0;
        memcpy(cell_buf, "    <row r=\"", 12); rlen = 12;
        rlen += ctools_int64_to_str((int64_t)row_num, cell_buf + rlen);
        memcpy(cell_buf + rlen, "\">\n", 3); rlen += 3;
        strbuf_append_n(sb, cell_buf, (size_t)rlen);

        for (ST_int j = 0; j < ctx->nvars; j++) {
            const char *prefix = a->col_prefix[j];
            size_t prefix_len = a->col_prefix_len[j];
            int vartype = ctx->vartypes[j];

            if (vartype == 0) {
                /* String variable */
                const char *str = filtered->data.vars[j].data.str[i];
                if (!str || str[0] == '\0') {
                    if (a->escaped_missing) {
                        int pos = 0;
                        memcpy(cell_buf, prefix, prefix_len); pos = (int)prefix_len;
                        pos += ctools_int64_to_str((int64_t)row_num, cell_buf + pos);
                        memcpy(cell_buf + pos, "\" t=\"inlineStr\"><is><t>", 23); pos += 23;
                        strbuf_append_n(sb, cell_buf, (size_t)pos);
                        strbuf_append_n(sb, a->escaped_missing, a->escaped_missing_len);
                        strbuf_append_n(sb, "</t></is></c>\n", 14);
                    } else {
                        int pos = 0;
                        memcpy(cell_buf, prefix, prefix_len); pos = (int)prefix_len;
                        pos += ctools_int64_to_str((int64_t)row_num, cell_buf + pos);
                        memcpy(cell_buf + pos, "\"/>\n", 4); pos += 4;
                        strbuf_append_n(sb, cell_buf, (size_t)pos);
                    }
                } else {
                    size_t esc_len = xlsx_escape_xml(str, escaped, MAX_CELL_CONTENT * 6 + 1);
                    int pos = 0;
                    memcpy(cell_buf, prefix, prefix_len); pos = (int)prefix_len;
                    pos += ctools_int64_to_str((int64_t)row_num, cell_buf + pos);
                    memcpy(cell_buf + pos, "\" t=\"inlineStr\"><is><t>", 23); pos += 23;
                    strbuf_append_n(sb, cell_buf, (size_t)pos);
                    strbuf_append_n(sb, escaped, esc_len);
                    strbuf_append_n(sb, "</t></is></c>\n", 14);
                }
            } else {
                /* Numeric variable */
                double val = filtered->data.vars[j].data.dbl[i];

                if (SF_is_missing(val)) {
                    if (a->escaped_missing) {
                        int pos = 0;
                        memcpy(cell_buf, prefix, prefix_len); pos = (int)prefix_len;
                        pos += ctools_int64_to_str((int64_t)row_num, cell_buf + pos);
                        memcpy(cell_buf + pos, "\" t=\"inlineStr\"><is><t>", 23); pos += 23;
                        strbuf_append_n(sb, cell_buf, (size_t)pos);
                        strbuf_append_n(sb, a->escaped_missing, a->escaped_missing_len);
                        strbuf_append_n(sb, "</t></is></c>\n", 14);
                    } else {
                        int pos = 0;
                        memcpy(cell_buf, prefix, prefix_len); pos = (int)prefix_len;
                        pos += ctools_int64_to_str((int64_t)row_num, cell_buf + pos);
                        memcpy(cell_buf + pos, "\"/>\n", 4); pos += 4;
                        strbuf_append_n(sb, cell_buf, (size_t)pos);
                    }
                } else {
                    /* Build entire numeric cell in cell_buf, single append */
                    int pos = 0;
                    memcpy(cell_buf, prefix, prefix_len); pos = (int)prefix_len;
                    pos += ctools_int64_to_str((int64_t)row_num, cell_buf + pos);
                    memcpy(cell_buf + pos, "\"><v>", 5); pos += 5;

                    if (val == floor(val) && fabs(val) < 1e15) {
                        pos += ctools_int64_to_str((int64_t)val, cell_buf + pos);
                    } else {
                        pos += xlsx_fast_dtoa(val, cell_buf + pos);
                    }

                    memcpy(cell_buf + pos, "</v></c>\n", 9); pos += 9;
                    strbuf_append_n(sb, cell_buf, (size_t)pos);
                }
            }
        }

        strbuf_append_n(sb, "    </row>\n", 11);
        row_num++;
    }

    free(escaped);
    return NULL;
}

/* Generate worksheet XML as segments (preamble + chunks + footer).
 * Returns true on success, false on failure.
 * Caller must call xlsx_segments_free(segs) when done. */
static bool xlsx_gen_worksheet_segments(XLSXExportContext *ctx,
                                         ctools_filtered_data *filtered,
                                         xlsx_worksheet_segments *segs) {
    memset(segs, 0, sizeof(*segs));
    size_t nobs_filtered = filtered->data.nobs;

    /* Pre-compute column letters for all variables */
    char (*col_letters_arr)[8] = (char (*)[8])malloc(ctx->nvars * 8);
    size_t *col_letters_len = (size_t *)malloc(ctx->nvars * sizeof(size_t));
    if (!col_letters_arr || !col_letters_len) {
        free(col_letters_arr);
        free(col_letters_len);
        return false;
    }
    for (ST_int j = 0; j < ctx->nvars; j++) {
        xlsx_col_to_letters(ctx->start_col + j, col_letters_arr[j]);
        col_letters_len[j] = strlen(col_letters_arr[j]);
    }

    /* Pre-compute cell prefixes: "      <c r=\"" + column letters */
    char (*col_prefix)[24] = (char (*)[24])malloc(ctx->nvars * 24);
    size_t *col_prefix_len = (size_t *)malloc(ctx->nvars * sizeof(size_t));
    if (!col_prefix || !col_prefix_len) {
        free(col_letters_arr);
        free(col_letters_len);
        free(col_prefix);
        free(col_prefix_len);
        return false;
    }
    for (ST_int j = 0; j < ctx->nvars; j++) {
        memcpy(col_prefix[j], "      <c r=\"", 12);
        memcpy(col_prefix[j] + 12, col_letters_arr[j], col_letters_len[j]);
        col_prefix_len[j] = 12 + col_letters_len[j];
    }

    /* Prepare escaped missing value */
    char escaped_missing[MAX_CELL_CONTENT * 6 + 1];
    size_t escaped_missing_len = 0;
    if (ctx->missing_value) {
        escaped_missing_len = xlsx_escape_xml(ctx->missing_value, escaped_missing, sizeof(escaped_missing));
    }

    /* Build XML preamble (header, dimension, optional header row) */
    StringBuffer preamble;
    if (!strbuf_init(&preamble, 4096)) {
        free(col_letters_arr);
        free(col_letters_len);
        free(col_prefix);
        free(col_prefix_len);
        return false;
    }

    strbuf_append(&preamble, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n");
    strbuf_append(&preamble, "<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">\n");

    char start_col_letters[8], end_col_letters[8];
    xlsx_col_to_letters(ctx->start_col, start_col_letters);
    xlsx_col_to_letters(ctx->start_col + ctx->nvars - 1, end_col_letters);
    long long total_rows = (long long)nobs_filtered + (ctx->firstrow ? 1 : 0);
    long long end_row = (long long)ctx->start_row + total_rows - 1;
    strbuf_appendf(&preamble, "  <dimension ref=\"%s%d:%s%lld\"/>\n",
                  start_col_letters, ctx->start_row, end_col_letters, end_row);

    strbuf_append(&preamble, "  <sheetData>\n");

    long long data_start_row = ctx->start_row;

    /* Header row if firstrow option */
    if (ctx->firstrow) {
        strbuf_appendf(&preamble, "    <row r=\"%lld\">\n", data_start_row);
        for (ST_int j = 0; j < ctx->nvars; j++) {
            char escaped_name[256];
            xlsx_escape_xml(ctx->varnames[j], escaped_name, sizeof(escaped_name));
            strbuf_appendf(&preamble, "      <c r=\"%s%lld\" t=\"inlineStr\"><is><t>%s</t></is></c>\n",
                          col_letters_arr[j], data_start_row, escaped_name);
        }
        strbuf_append(&preamble, "    </row>\n");
        data_start_row++;
    }

    /* Footer (static string, not freed with segments) */
    static const char footer[] = "  </sheetData>\n</worksheet>";
    size_t footer_len = sizeof(footer) - 1;

    /* Adaptive chunk sizing: target ~384KB per chunk (fits in L2 cache). */
    size_t est_bytes_per_row = (size_t)ctx->nvars * 50 + 50;
    size_t chunk_size;
    if (est_bytes_per_row > 0) {
        chunk_size = (384 * 1024) / est_bytes_per_row;
        if (chunk_size < 1000) chunk_size = 1000;
        if (chunk_size > 50000) chunk_size = 50000;
    } else {
        chunk_size = CTOOLS_EXPORT_CHUNK_SIZE;
    }
    int num_threads = ctools_get_max_threads();
    if (num_threads > 1 && nobs_filtered > 0) {
        size_t min_chunks = (size_t)num_threads * 4;
        size_t max_chunk_for_threads = (nobs_filtered + min_chunks - 1) / min_chunks;
        if (max_chunk_for_threads < chunk_size && max_chunk_for_threads >= 1000) {
            chunk_size = max_chunk_for_threads;
        }
    }
    size_t num_chunks = nobs_filtered > 0 ? (nobs_filtered + chunk_size - 1) / chunk_size : 1;

    /* Allocate chunk args */
    xlsx_chunk_args_t *chunks = calloc(num_chunks, sizeof(xlsx_chunk_args_t));
    if (!chunks) {
        free(preamble.data);
        free(col_letters_arr);
        free(col_letters_len);
        free(col_prefix);
        free(col_prefix_len);
        return false;
    }

    /* Initialize chunks */
    bool init_ok = true;
    for (size_t c = 0; c < num_chunks; c++) {
        chunks[c].ctx = ctx;
        chunks[c].filtered = filtered;
        chunks[c].start_idx = c * chunk_size;
        chunks[c].end_idx = (c + 1) * chunk_size;
        if (chunks[c].end_idx > nobs_filtered) chunks[c].end_idx = nobs_filtered;
        chunks[c].start_row_num = data_start_row + (long long)(c * chunk_size);
        chunks[c].col_letters_arr = col_letters_arr;
        chunks[c].col_letters_len = col_letters_len;
        chunks[c].col_prefix = col_prefix;
        chunks[c].col_prefix_len = col_prefix_len;
        chunks[c].escaped_missing = ctx->missing_value ? escaped_missing : NULL;
        chunks[c].escaped_missing_len = escaped_missing_len;

        size_t chunk_rows = chunks[c].end_idx - chunks[c].start_idx;
        size_t est_chunk_size = chunk_rows * est_bytes_per_row;
        if (est_chunk_size < INITIAL_XML_SIZE / 4) est_chunk_size = INITIAL_XML_SIZE / 4;

        if (!strbuf_init(&chunks[c].sb, est_chunk_size)) {
            init_ok = false;
            break;
        }
    }

    if (!init_ok) {
        for (size_t c = 0; c < num_chunks; c++) {
            free(chunks[c].sb.data);
        }
        free(chunks);
        free(preamble.data);
        free(col_letters_arr);
        free(col_letters_len);
        free(col_prefix);
        free(col_prefix_len);
        return false;
    }

    /* Execute chunks: parallel if pool available and enough work */
    bool use_parallel = (num_chunks > 1 && num_threads > 1);
    ctools_persistent_pool *pool = NULL;

    if (use_parallel) {
        pool = ctools_get_global_pool();
        if (!pool) use_parallel = false;
    }

    if (use_parallel) {
        for (size_t c = 0; c < num_chunks; c++) {
            ctools_persistent_pool_submit(pool, xlsx_format_chunk, &chunks[c]);
        }
        ctools_persistent_pool_wait(pool);
    } else {
        for (size_t c = 0; c < num_chunks; c++) {
            xlsx_format_chunk(&chunks[c]);
        }
    }

    /* Build segment array: preamble + num_chunks + footer */
    size_t total_segments = 1 + num_chunks + 1;
    segs->segments = (xlsx_segment *)malloc(total_segments * sizeof(xlsx_segment));
    if (!segs->segments) {
        for (size_t c = 0; c < num_chunks; c++) {
            free(chunks[c].sb.data);
        }
        free(chunks);
        free(preamble.data);
        free(col_letters_arr);
        free(col_letters_len);
        free(col_prefix);
        free(col_prefix_len);
        return false;
    }

    segs->num_segments = total_segments;
    segs->segments[0].data = preamble.data;
    segs->segments[0].len = preamble.len;

    /* Transfer chunk StringBuffers to segments */
    segs->chunk_bufs = (StringBuffer *)malloc(num_chunks * sizeof(StringBuffer));
    segs->num_chunks = num_chunks;
    size_t total_size = preamble.len + footer_len;

    for (size_t c = 0; c < num_chunks; c++) {
        segs->segments[1 + c].data = chunks[c].sb.data;
        segs->segments[1 + c].len = chunks[c].sb.len;
        if (segs->chunk_bufs) {
            segs->chunk_bufs[c] = chunks[c].sb;
        }
        total_size += chunks[c].sb.len;
    }

    segs->segments[total_segments - 1].data = footer;
    segs->segments[total_segments - 1].len = footer_len;

    segs->total_size = total_size;
    segs->preamble_buf = preamble;
    segs->col_letters_arr = col_letters_arr;
    segs->col_letters_len = col_letters_len;
    segs->col_prefix = col_prefix;
    segs->col_prefix_len = col_prefix_len;

    free(chunks);
    return true;
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

    /* Generate worksheet as segments and stream to ZIP.
     * Uses read callback API to avoid allocating a single contiguous buffer
     * for the full worksheet XML (saves ~500MB for 10M×10 datasets). */
    xlsx_worksheet_segments segs;
    if (!xlsx_gen_worksheet_segments(ctx, filtered, &segs)) {
        mz_zip_writer_end(&zip);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Memory allocation failed for worksheet\n");
        return 920;
    }

    if (!mz_zip_writer_add_read_buf_callback(&zip, "xl/worksheets/sheet1.xml",
            xlsx_segment_read_callback, &segs, segs.total_size,
            NULL, NULL, 0, MZ_BEST_SPEED, NULL, 0, NULL, 0)) {
        xlsx_segments_free(&segs);
        mz_zip_writer_end(&zip);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to write xl/worksheets/sheet1.xml\n");
        return 603;
    }
    xlsx_segments_free(&segs);

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

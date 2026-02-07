/*
 * cimport_xlsx.c
 * High-performance XLSX import for Stata
 *
 * Implements Excel file parsing using miniz for ZIP extraction
 * and custom XML parsing for worksheet data.
 */

#include "cimport_xlsx.h"
#include "cimport_xlsx_zip.h"
#include "cimport_xlsx_xml.h"
#include "../ctools_runtime.h"
#include "../ctools_arena.h"
#include "../ctools_threads.h"
#include "../ctools_types.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

/* Fast inline parser for non-negative integers (XLSX shared string indices, row numbers, etc.)
 * Returns -1 on empty/invalid input. */
static inline int xlsx_fast_atoi(const char *s)
{
    if (!s || !*s) return -1;
    int val = 0;
    while (*s >= '0' && *s <= '9') {
        val = val * 10 + (*s - '0');
        s++;
    }
    return val;
}

/* ============================================================================
 * Forward Declarations
 * ============================================================================ */

static bool workbook_callback(const xlsx_xml_event *event, void *user_data);
static bool shared_strings_callback(const xlsx_xml_event *event, void *user_data);
static bool worksheet_callback(const xlsx_xml_event *event, void *user_data);
static bool styles_callback(const xlsx_xml_event *event, void *user_data);
static void generate_varname(char *buf, int col, const char *header_value, int case_mode);
static double excel_date_to_stata(double excel_date);

/* Global cached context (survives between scan and load calls) */
static XLSXContext *g_xlsx_ctx = NULL;

void xlsx_clear_cached_context(void)
{
    if (g_xlsx_ctx) {
        xlsx_context_free(g_xlsx_ctx);
        free(g_xlsx_ctx);
        g_xlsx_ctx = NULL;
    }
}

/* Parser state structures */
typedef struct {
    XLSXContext *ctx;
    bool in_sheets;
} WorkbookParseState;

typedef struct {
    XLSXContext *ctx;
    bool in_si;
    bool in_t;
    char current_string[4096];
    size_t current_len;
} SharedStringsParseState;

typedef struct {
    XLSXContext *ctx;
    bool in_sheetdata;
    bool in_row;
    bool in_cell;
    bool in_value;
    bool in_inline_str;
    bool in_t;
    int current_row;
    int current_col;
    XLSXCellType current_type;
    int current_style;
    char value_buf[4096];
    size_t value_len;
    XLSXParsedRow *row_ptr;
} WorksheetParseState;

typedef struct {
    XLSXContext *ctx;
    bool in_numfmts;
    bool in_cellxfs;
    int current_numfmt_id;
    int current_xf_index;
    /* Track which numFmtIds are dates */
    bool is_date_numfmt[256];
} StylesParseState;

/* ============================================================================
 * Context Management
 * ============================================================================ */

void xlsx_context_init(XLSXContext *ctx)
{
    memset(ctx, 0, sizeof(XLSXContext));
    ctx->cell_range.start_col = -1;
    ctx->cell_range.start_row = -1;
    ctx->cell_range.end_col = -1;
    ctx->cell_range.end_row = -1;
    ctx->selected_sheet = 0;
    ctx->min_row = INT32_MAX;
    ctools_arena_init(&ctx->parse_arena, 0);
    ctools_arena_init(&ctx->cells_arena, 0);
}

void xlsx_context_free(XLSXContext *ctx)
{
    if (!ctx) return;

    /* Free zip archive */
    if (ctx->zip_archive) {
        xlsx_zip_close((xlsx_zip_archive *)ctx->zip_archive);
        ctx->zip_archive = NULL;
    }

    /* Free shared strings */
    if (ctx->shared_strings) {
        free(ctx->shared_strings);
        ctx->shared_strings = NULL;
    }
    if (ctx->shared_string_lengths) {
        free(ctx->shared_string_lengths);
        ctx->shared_string_lengths = NULL;
    }
    if (ctx->shared_strings_pool) {
        free(ctx->shared_strings_pool);
        ctx->shared_strings_pool = NULL;
    }

    /* Free date styles */
    if (ctx->date_styles) {
        free(ctx->date_styles);
        ctx->date_styles = NULL;
    }

    /* Free parsed rows (inline strings freed via parse_arena below) */
    if (ctx->rows) {
        for (int i = 0; i < ctx->num_rows; i++) {
            if (ctx->rows[i].cells && !ctx->rows[i].cells_from_arena) {
                free(ctx->rows[i].cells);
            }
        }
        free(ctx->rows);
        ctx->rows = NULL;
    }

    /* Free column-major arrays */
    if (ctx->cm_numeric) {
        for (int c = 0; c < ctx->cm_num_cols; c++) {
            free(ctx->cm_numeric[c]);
        }
        free(ctx->cm_numeric);
        ctx->cm_numeric = NULL;
    }
    if (ctx->cm_types) {
        for (int c = 0; c < ctx->cm_num_cols; c++) {
            free(ctx->cm_types[c]);
        }
        free(ctx->cm_types);
        ctx->cm_types = NULL;
    }
    if (ctx->cm_styles) {
        free(ctx->cm_styles);
        ctx->cm_styles = NULL;
    }

    /* Free arenas */
    ctools_arena_free(&ctx->cells_arena);
    ctools_arena_free(&ctx->parse_arena);

    /* Free col_stats */
    if (ctx->col_stats) {
        free(ctx->col_stats);
        ctx->col_stats = NULL;
    }

    /* Free columns */
    if (ctx->columns) {
        free(ctx->columns);
        ctx->columns = NULL;
    }

    /* Free column caches */
    if (ctx->col_cache) {
        for (int i = 0; i < ctx->num_columns; i++) {
            if (ctx->col_cache[i].string_data) {
                free(ctx->col_cache[i].string_data);
            }
            if (ctx->col_cache[i].numeric_data) {
                free(ctx->col_cache[i].numeric_data);
            }
            ctools_arena_free(&ctx->col_cache[i].string_arena);
        }
        free(ctx->col_cache);
        ctx->col_cache = NULL;
    }

    /* Free filename */
    if (ctx->filename) {
        free(ctx->filename);
        ctx->filename = NULL;
    }
}

/* ============================================================================
 * File Operations
 * ============================================================================ */

ST_retcode xlsx_open_file(XLSXContext *ctx, const char *filename)
{
    ctx->filename = strdup(filename);
    if (!ctx->filename) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to allocate memory for filename");
        return 920;
    }

    double t_start = ctools_timer_seconds();

    ctx->zip_archive = xlsx_zip_open_file(filename);
    if (!ctx->zip_archive) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to open XLSX file: %s", filename);
        return 601;
    }

    ctx->time_zip = ctools_timer_seconds() - t_start;

    /* Verify it's a valid XLSX (has required files) */
    xlsx_zip_archive *zip = (xlsx_zip_archive *)ctx->zip_archive;
    if (xlsx_zip_locate_file(zip, XLSX_PATH_WORKBOOK) == (size_t)-1) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Invalid XLSX file: missing workbook.xml");
        return 610;
    }

    return 0;
}

/* ============================================================================
 * Workbook Parsing (get sheet list)
 * ============================================================================ */

static bool workbook_callback(const xlsx_xml_event *event, void *user_data)
{
    WorkbookParseState *state = (WorkbookParseState *)user_data;
    XLSXContext *ctx = state->ctx;

    if (event->type == XLSX_XML_START_ELEMENT) {
        if (strcmp(event->tag_name, "sheets") == 0) {
            state->in_sheets = true;
        }
        else if (state->in_sheets && strcmp(event->tag_name, "sheet") == 0) {
            if (ctx->num_sheets >= XLSX_MAX_SHEETS) return true;

            XLSXSheetInfo *info = &ctx->sheets[ctx->num_sheets];

            const char *name = xlsx_xml_get_attr(event, "name");
            const char *sheet_id = xlsx_xml_get_attr(event, "sheetId");
            const char *rel_id = xlsx_xml_get_attr(event, "id");  /* r:id */

            if (name) {
                strncpy(info->name, name, XLSX_MAX_SHEET_NAME - 1);
                info->name[XLSX_MAX_SHEET_NAME - 1] = '\0';
            }
            if (sheet_id) {
                info->sheet_id = atoi(sheet_id);
            }
            if (rel_id) {
                strncpy(info->rel_id, rel_id, sizeof(info->rel_id) - 1);
            }

            /* Default sheet index is based on position */
            info->sheet_index = ctx->num_sheets + 1;

            ctx->num_sheets++;
        }
    }
    else if (event->type == XLSX_XML_END_ELEMENT) {
        if (strcmp(event->tag_name, "sheets") == 0) {
            state->in_sheets = false;
        }
    }

    return true;
}

ST_retcode xlsx_parse_workbook(XLSXContext *ctx)
{
    xlsx_zip_archive *zip = (xlsx_zip_archive *)ctx->zip_archive;

    size_t size;
    void *data = xlsx_zip_extract_file(zip, XLSX_PATH_WORKBOOK, &size);
    if (!data) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to extract workbook.xml");
        return 610;
    }

    WorkbookParseState state = {0};
    state.ctx = ctx;

    xlsx_xml_parser *parser = xlsx_xml_parser_create(workbook_callback, &state);
    if (!parser) {
        free(data);
        return 920;
    }

    bool success = xlsx_xml_parse(parser, (const char *)data, size, true);

    xlsx_xml_parser_destroy(parser);
    free(data);

    if (!success) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to parse workbook.xml");
        return 610;
    }

    return 0;
}

/* ============================================================================
 * Shared Strings Parsing
 * ============================================================================ */

static bool shared_strings_callback(const xlsx_xml_event *event, void *user_data)
{
    SharedStringsParseState *state = (SharedStringsParseState *)user_data;
    XLSXContext *ctx = state->ctx;

    if (event->type == XLSX_XML_START_ELEMENT) {
        if (strcmp(event->tag_name, "sst") == 0) {
            /* Pre-allocate from uniqueCount if available */
            const char *unique_count = xlsx_xml_get_attr(event, "uniqueCount");
            if (unique_count) {
                uint32_t count = (uint32_t)strtoul(unique_count, NULL, 10);
                if (count > 0 && count <= XLSX_MAX_SHARED_STRINGS) {
                    ctx->shared_strings = (char **)malloc(count * sizeof(char *));
                    if (ctx->shared_strings) {
                        ctx->shared_strings_capacity = count;
                    }
                    ctx->shared_string_lengths = (uint32_t *)malloc(count * sizeof(uint32_t));
                    /* Pre-allocate pool: estimate 32 bytes per string */
                    size_t pool_est = (size_t)count * 32;
                    ctx->shared_strings_pool = (char *)malloc(pool_est);
                    if (ctx->shared_strings_pool) {
                        ctx->shared_strings_pool_size = pool_est;
                    }
                }
            }
        }
        else if (strcmp(event->tag_name, "si") == 0) {
            state->in_si = true;
            state->current_len = 0;
            state->current_string[0] = '\0';
        }
        else if (state->in_si && strcmp(event->tag_name, "t") == 0) {
            state->in_t = true;
        }
    }
    else if (event->type == XLSX_XML_TEXT) {
        if (state->in_t) {
            /* Append text to current string */
            size_t copy_len = event->text_len;
            if (state->current_len + copy_len >= sizeof(state->current_string)) {
                copy_len = sizeof(state->current_string) - state->current_len - 1;
            }
            if (copy_len > 0) {
                memcpy(state->current_string + state->current_len,
                       event->text, copy_len);
                state->current_len += copy_len;
                state->current_string[state->current_len] = '\0';
            }
        }
    }
    else if (event->type == XLSX_XML_END_ELEMENT) {
        if (strcmp(event->tag_name, "t") == 0) {
            state->in_t = false;
        }
        else if (strcmp(event->tag_name, "si") == 0) {
            state->in_si = false;

            /* Add string to shared strings table */
            if (ctx->num_shared_strings >= ctx->shared_strings_capacity) {
                size_t new_cap = ctx->shared_strings_capacity == 0 ?
                                 1024 : ctx->shared_strings_capacity * 2;
                char **new_strings = (char **)realloc(ctx->shared_strings,
                                                       new_cap * sizeof(char *));
                if (!new_strings) return false;
                ctx->shared_strings = new_strings;
                uint32_t *new_lengths = (uint32_t *)realloc(ctx->shared_string_lengths,
                                                             new_cap * sizeof(uint32_t));
                if (!new_lengths) return false;
                ctx->shared_string_lengths = new_lengths;
                ctx->shared_strings_capacity = (uint32_t)new_cap;
            }

            /* Decode XML entities */
            size_t decoded_len = xlsx_xml_decode_entities(state->current_string,
                                                          state->current_len);

            /* Allocate from pool */
            size_t need = decoded_len + 1;
            if (ctx->shared_strings_pool_used + need > ctx->shared_strings_pool_size) {
                /* Expand pool */
                size_t new_size = ctx->shared_strings_pool_size == 0 ?
                                  (1024 * 1024) : ctx->shared_strings_pool_size * 2;
                while (ctx->shared_strings_pool_used + need > new_size) {
                    /* Check for overflow before doubling */
                    if (new_size > SIZE_MAX / 2) return false;
                    new_size *= 2;
                }
                char *new_pool = (char *)realloc(ctx->shared_strings_pool, new_size);
                if (!new_pool) return false;

                /* Update pointers if pool moved */
                if (ctx->shared_strings_pool && new_pool != ctx->shared_strings_pool) {
                    ptrdiff_t diff = new_pool - ctx->shared_strings_pool;
                    for (uint32_t i = 0; i < ctx->num_shared_strings; i++) {
                        ctx->shared_strings[i] += diff;
                    }
                }
                ctx->shared_strings_pool = new_pool;
                ctx->shared_strings_pool_size = new_size;
            }

            char *dest = ctx->shared_strings_pool + ctx->shared_strings_pool_used;
            memcpy(dest, state->current_string, decoded_len + 1);
            ctx->shared_strings[ctx->num_shared_strings] = dest;
            if (ctx->shared_string_lengths) {
                ctx->shared_string_lengths[ctx->num_shared_strings] = (uint32_t)decoded_len;
            }
            ctx->num_shared_strings++;
            ctx->shared_strings_pool_used += need;
        }
    }

    return true;
}

/* Inner: parse shared strings from pre-extracted buffer (does not free data) */
static ST_retcode xlsx_parse_shared_strings_buf(XLSXContext *ctx, const void *data, size_t size)
{
    SharedStringsParseState state = {0};
    state.ctx = ctx;

    xlsx_xml_parser *parser = xlsx_xml_parser_create(shared_strings_callback, &state);
    if (!parser) return 920;

    bool success = xlsx_xml_parse(parser, (const char *)data, size, true);
    xlsx_xml_parser_destroy(parser);

    if (!success) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to parse sharedStrings.xml");
        return 610;
    }
    return 0;
}

ST_retcode xlsx_parse_shared_strings(XLSXContext *ctx)
{
    xlsx_zip_archive *zip = (xlsx_zip_archive *)ctx->zip_archive;

    double t_start = ctools_timer_seconds();

    /* sharedStrings.xml may not exist if all cells are numbers */
    size_t idx = xlsx_zip_locate_file(zip, XLSX_PATH_SHARED_STRINGS);
    if (idx == (size_t)-1) {
        ctx->time_shared_strings = ctools_timer_seconds() - t_start;
        return 0;  /* Not an error */
    }

    size_t size;
    void *data = xlsx_zip_extract_to_heap(zip, idx, &size);
    if (!data) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to extract sharedStrings.xml");
        return 610;
    }

    ST_retcode rc = xlsx_parse_shared_strings_buf(ctx, data, size);
    free(data);

    ctx->time_shared_strings = ctools_timer_seconds() - t_start;
    return rc;
}

/* ============================================================================
 * Styles Parsing (date format detection)
 * ============================================================================ */

/* Built-in Excel number format IDs that are dates */
static bool is_builtin_date_format(int numFmtId)
{
    /* Standard date/time formats: 14-22, 27-36, 45-47, 50-58 */
    if (numFmtId >= 14 && numFmtId <= 22) return true;
    if (numFmtId >= 27 && numFmtId <= 36) return true;
    if (numFmtId >= 45 && numFmtId <= 47) return true;
    if (numFmtId >= 50 && numFmtId <= 58) return true;
    return false;
}

static bool styles_callback(const xlsx_xml_event *event, void *user_data)
{
    StylesParseState *state = (StylesParseState *)user_data;
    XLSXContext *ctx = state->ctx;

    if (event->type == XLSX_XML_START_ELEMENT) {
        if (strcmp(event->tag_name, "numFmts") == 0) {
            state->in_numfmts = true;
        }
        else if (state->in_numfmts && strcmp(event->tag_name, "numFmt") == 0) {
            const char *numFmtId = xlsx_xml_get_attr(event, "numFmtId");
            const char *formatCode = xlsx_xml_get_attr(event, "formatCode");

            if (numFmtId && formatCode) {
                int id = atoi(numFmtId);
                /* Check if format code looks like a date */
                bool is_date = false;
                const char *p = formatCode;
                while (*p) {
                    /* Look for date-related format characters */
                    if (*p == 'd' || *p == 'm' || *p == 'y' ||
                        *p == 'D' || *p == 'M' || *p == 'Y') {
                        is_date = true;
                        break;
                    }
                    p++;
                }
                if (id >= 0 && id < 256) {
                    state->is_date_numfmt[id] = is_date;
                }
            }
        }
        else if (strcmp(event->tag_name, "cellXfs") == 0) {
            state->in_cellxfs = true;
            state->current_xf_index = 0;
        }
        else if (state->in_cellxfs && strcmp(event->tag_name, "xf") == 0) {
            const char *numFmtId = xlsx_xml_get_attr(event, "numFmtId");
            bool is_date = false;

            if (numFmtId) {
                int id = atoi(numFmtId);
                is_date = is_builtin_date_format(id);
                if (!is_date && id >= 0 && id < 256) {
                    is_date = state->is_date_numfmt[id];
                }
            }

            /* Expand date_styles array if needed */
            if (state->current_xf_index >= ctx->num_styles) {
                int new_size = ctx->num_styles == 0 ? 64 : ctx->num_styles * 2;
                bool *new_styles = (bool *)realloc(ctx->date_styles,
                                                    new_size * sizeof(bool));
                if (!new_styles) {
                    return false;  /* Memory allocation failed */
                }
                memset(new_styles + ctx->num_styles, 0,
                       (new_size - ctx->num_styles) * sizeof(bool));
                ctx->date_styles = new_styles;
                ctx->num_styles = new_size;
            }

            if (state->current_xf_index < ctx->num_styles) {
                ctx->date_styles[state->current_xf_index] = is_date;
            }

            state->current_xf_index++;
        }
    }
    else if (event->type == XLSX_XML_END_ELEMENT) {
        if (strcmp(event->tag_name, "numFmts") == 0) {
            state->in_numfmts = false;
        }
        else if (strcmp(event->tag_name, "cellXfs") == 0) {
            state->in_cellxfs = false;
        }
    }

    return true;
}

/* Inner: parse styles from pre-extracted buffer (does not free data) */
static ST_retcode xlsx_parse_styles_buf(XLSXContext *ctx, const void *data, size_t size)
{
    StylesParseState state = {0};
    state.ctx = ctx;

    xlsx_xml_parser *parser = xlsx_xml_parser_create(styles_callback, &state);
    if (!parser) return 0;

    xlsx_xml_parse(parser, (const char *)data, size, true);
    xlsx_xml_parser_destroy(parser);

    return 0;
}

ST_retcode xlsx_parse_styles(XLSXContext *ctx)
{
    xlsx_zip_archive *zip = (xlsx_zip_archive *)ctx->zip_archive;

    size_t idx = xlsx_zip_locate_file(zip, XLSX_PATH_STYLES);
    if (idx == (size_t)-1) {
        return 0;  /* Styles are optional */
    }

    size_t size;
    void *data = xlsx_zip_extract_to_heap(zip, idx, &size);
    if (!data) {
        return 0;  /* Non-fatal */
    }

    ST_retcode rc = xlsx_parse_styles_buf(ctx, data, size);
    free(data);

    return rc;
}

/* ============================================================================
 * Worksheet Parsing
 * ============================================================================ */

/* Ensure col_stats covers at least col_idx+1 entries */
static CImportColumnInfo *ensure_col_stat(XLSXContext *ctx, int col_idx)
{
    if (col_idx < ctx->col_stats_capacity) {
        return &ctx->col_stats[col_idx];
    }
    int new_cap = ctx->col_stats_capacity == 0 ? 64 : ctx->col_stats_capacity;
    while (new_cap <= col_idx) new_cap *= 2;
    CImportColumnInfo *new_stats = (CImportColumnInfo *)realloc(
        ctx->col_stats, new_cap * sizeof(CImportColumnInfo));
    if (!new_stats) return NULL;
    /* Initialize new entries */
    for (int i = ctx->col_stats_capacity; i < new_cap; i++) {
        memset(&new_stats[i], 0, sizeof(CImportColumnInfo));
        new_stats[i].type = CIMPORT_COL_UNKNOWN;
        new_stats[i].num_subtype = CIMPORT_NUM_BYTE;
        new_stats[i].is_integer = true;
        new_stats[i].min_value = DBL_MAX;
        new_stats[i].max_value = -DBL_MAX;
    }
    ctx->col_stats = new_stats;
    ctx->col_stats_capacity = new_cap;
    return &ctx->col_stats[col_idx];
}

static XLSXParsedRow *ensure_row(XLSXContext *ctx, int row_num)
{
    /* Fast path: check last row (rows arrive in order from XML) */
    if (ctx->num_rows > 0 && ctx->rows[ctx->num_rows - 1].row_num == row_num) {
        return &ctx->rows[ctx->num_rows - 1];
    }

    /* Slow path: linear scan (for out-of-order or revisited rows) */
    for (int i = 0; i < ctx->num_rows; i++) {
        if (ctx->rows[i].row_num == row_num) {
            return &ctx->rows[i];
        }
    }

    /* Create new row */
    if (ctx->num_rows >= ctx->rows_capacity) {
        int new_cap = ctx->rows_capacity == 0 ? 1024 : ctx->rows_capacity * 2;
        XLSXParsedRow *new_rows = (XLSXParsedRow *)realloc(ctx->rows,
                                                           new_cap * sizeof(XLSXParsedRow));
        if (!new_rows) return NULL;
        ctx->rows = new_rows;
        ctx->rows_capacity = new_cap;
    }

    XLSXParsedRow *row = &ctx->rows[ctx->num_rows++];
    memset(row, 0, sizeof(XLSXParsedRow));
    row->row_num = row_num;

    /* Pre-allocate cells from arena if dimension is known */
    if (ctx->dimension_cols > 0) {
        int cell_cap = ctx->dimension_cols;
        XLSXCell *cells = (XLSXCell *)ctools_arena_alloc(
            &ctx->cells_arena, (size_t)cell_cap * sizeof(XLSXCell));
        if (cells) {
            memset(cells, 0, (size_t)cell_cap * sizeof(XLSXCell));
            row->cells = cells;
            row->capacity = cell_cap;
            row->cells_from_arena = true;
        }
    }

    if (row_num < ctx->min_row) ctx->min_row = row_num;
    if (row_num > ctx->max_row) ctx->max_row = row_num;

    return row;
}

static XLSXCell *add_cell(XLSXParsedRow *row, int col, XLSXContext *ctx)
{
    if (row->num_cells >= row->capacity) {
        int new_cap = row->capacity == 0 ? 16 : row->capacity * 2;
        if (row->cells_from_arena) {
            /* Arena cells can't be realloc'd — alloc new block and copy */
            XLSXCell *new_cells = (XLSXCell *)ctools_arena_alloc(
                &ctx->cells_arena, (size_t)new_cap * sizeof(XLSXCell));
            if (!new_cells) return NULL;
            memcpy(new_cells, row->cells, (size_t)row->num_cells * sizeof(XLSXCell));
            row->cells = new_cells;
            /* Old cells remain in arena — small waste, overflow is rare */
        } else {
            XLSXCell *new_cells = (XLSXCell *)realloc(row->cells,
                                                       new_cap * sizeof(XLSXCell));
            if (!new_cells) return NULL;
            row->cells = new_cells;
        }
        row->capacity = new_cap;
    }

    XLSXCell *cell = &row->cells[row->num_cells++];
    memset(cell, 0, sizeof(XLSXCell));
    cell->col = col;
    return cell;
}

/* ============================================================================
 * Fast Worksheet Scanner
 *
 * Replaces the SAX parser for worksheet parsing. Uses memchr to find tag
 * boundaries and processes complete cells as units, avoiding per-character
 * state machine overhead. Falls back to the SAX parser only when streaming
 * is unavailable.
 * ============================================================================ */

/* Find a short byte sequence in a buffer (like memmem but always available).
 * Optimized for short needles via memchr + memcmp. */
static inline const char *xlsx_memfind(const char *hay, size_t hlen,
                                       const char *needle, size_t nlen)
{
    if (nlen == 0) return hay;
    if (nlen > hlen) return NULL;
    const char *end = hay + hlen - nlen + 1;
    const char *p = hay;
    while (p < end) {
        p = (const char *)memchr(p, needle[0], (size_t)(end - p));
        if (!p) return NULL;
        if (nlen == 1 || memcmp(p + 1, needle + 1, nlen - 1) == 0) return p;
        p++;
    }
    return NULL;
}

/* Parse cell reference "A1", "BC99999" → (col 0-based, row 1-based).
 * Length-aware — does not require null termination. */
static inline bool scan_parse_cell_ref(const char *ref, size_t len,
                                       int *col, int *row)
{
    int c = 0, r = 0;
    size_t i = 0;
    while (i < len && ref[i] >= 'A' && ref[i] <= 'Z') {
        c = c * 26 + (ref[i] - 'A' + 1);
        i++;
    }
    if (c == 0) return false;
    c--;
    while (i < len && ref[i] >= '0' && ref[i] <= '9') {
        r = r * 10 + (ref[i] - '0');
        i++;
    }
    *col = c;
    *row = r;
    return (r > 0);
}

/* Parse non-negative integer from a non-null-terminated buffer */
static inline int scan_atoi(const char *s, size_t len)
{
    int val = 0;
    for (size_t i = 0; i < len && s[i] >= '0' && s[i] <= '9'; i++)
        val = val * 10 + (s[i] - '0');
    return val;
}

/* Find value of a single-character attribute within a tag's attribute region.
 * Looks for the pattern ` X="value"` where X is `attr`.
 * Returns pointer to value start, sets *vlen to value length. NULL if absent. */
static inline const char *scan_find_attr(const char *start, const char *end,
                                         char attr, size_t *vlen)
{
    for (const char *p = start; p + 3 < end; p++) {
        if (*p == ' ' && p[1] == attr && p[2] == '=' && p[3] == '"') {
            const char *val = p + 4;
            const char *q = (const char *)memchr(val, '"', (size_t)(end - val));
            if (q) { *vlen = (size_t)(q - val); return val; }
            *vlen = 0;
            return NULL;
        }
    }
    *vlen = 0;
    return NULL;
}

/* Handle <dimension ref="..."> — allocate cm arrays, rows, col_stats.
 * ref_str must be null-terminated. Shared by SAX callback and fast scanner. */
static void handle_dimension_ref(XLSXContext *ctx, const char *ref_str)
{
    const char *colon = strchr(ref_str, ':');
    if (!colon) return;

    int end_col, end_row;
    if (!xlsx_xml_parse_cell_ref(colon + 1, &end_col, &end_row)) return;

    ctx->dimension_cols = end_col + 1;
    ctx->dimension_rows = end_row;
    int ncols = ctx->dimension_cols;

    /* Pre-allocate rows array */
    if (ctx->rows_capacity == 0 && end_row > 0) {
        ctx->rows = (XLSXParsedRow *)calloc((size_t)end_row, sizeof(XLSXParsedRow));
        if (ctx->rows) ctx->rows_capacity = end_row;
    }

    /* Pre-allocate col_stats */
    if (ctx->col_stats_capacity == 0 && ncols > 0) {
        CImportColumnInfo *stats = (CImportColumnInfo *)calloc(
            (size_t)ncols, sizeof(CImportColumnInfo));
        if (stats) {
            for (int i = 0; i < ncols; i++) {
                stats[i].type = CIMPORT_COL_UNKNOWN;
                stats[i].num_subtype = CIMPORT_NUM_BYTE;
                stats[i].is_integer = true;
                stats[i].min_value = DBL_MAX;
                stats[i].max_value = -DBL_MAX;
            }
            ctx->col_stats = stats;
            ctx->col_stats_capacity = ncols;
        }
    }

    /* Allocate column-major arrays for direct write */
    if (ncols > 0 && end_row > 0 && !ctx->cm_active) {
        size_t nrows = (size_t)end_row;
        ctx->cm_numeric = (double **)calloc((size_t)ncols, sizeof(double *));
        ctx->cm_types = (uint8_t **)calloc((size_t)ncols, sizeof(uint8_t *));
        if (ctx->cm_numeric && ctx->cm_types) {
            bool alloc_ok = true;
            for (int c = 0; c < ncols && alloc_ok; c++) {
                ctx->cm_numeric[c] = (double *)malloc(nrows * sizeof(double));
                ctx->cm_types[c] = (uint8_t *)calloc(nrows, sizeof(uint8_t));
                if (!ctx->cm_numeric[c] || !ctx->cm_types[c]) {
                    alloc_ok = false;
                } else {
                    for (size_t r = 0; r < nrows; r++)
                        ctx->cm_numeric[c][r] = SV_missval;
                }
            }
            if (alloc_ok) {
                ctx->cm_num_rows = nrows;
                ctx->cm_num_cols = ncols;
                ctx->cm_active = true;
            } else {
                for (int c = 0; c < ncols; c++) {
                    free(ctx->cm_numeric[c]);
                    free(ctx->cm_types[c]);
                }
                free(ctx->cm_numeric);
                free(ctx->cm_types);
                ctx->cm_numeric = NULL;
                ctx->cm_types = NULL;
            }
        }
    }
}

/* Process one complete cell element from <c...> to </c>.
 * cell_start points to '<c', close_tag points to '</c>'.
 * Writes directly to cm arrays (or row-major fallback). */
static void xlsx_scan_cell(XLSXContext *ctx, const char *cell_start,
                           const char *close_tag, int fallback_row)
{
    /* Find end of opening <c ...> tag */
    const char *gt = (const char *)memchr(cell_start, '>',
                                          (size_t)(close_tag - cell_start));
    if (!gt) return;

    bool self_closing = (gt > cell_start && gt[-1] == '/');
    const char *attrs_start = cell_start + 2;  /* skip '<c' */
    const char *attrs_end = self_closing ? gt - 1 : gt;

    /* Extract r, t, s attributes */
    int col = -1, row = fallback_row;
    size_t vlen;

    const char *r_val = scan_find_attr(attrs_start, attrs_end, 'r', &vlen);
    if (r_val) scan_parse_cell_ref(r_val, vlen, &col, &row);
    if (col < 0) return;

    /* Apply cell range filter */
    XLSXCellRange *range = &ctx->cell_range;
    if ((range->start_col >= 0 && col < range->start_col) ||
        (range->end_col >= 0 && col > range->end_col))
        return;

    if (self_closing) return;  /* Empty cell <c .../> */

    const char *t_val = scan_find_attr(attrs_start, attrs_end, 't', &vlen);
    char type_first = (t_val && vlen > 0) ? t_val[0] : 0;
    size_t t_len = vlen;

    const char *s_val = scan_find_attr(attrs_start, attrs_end, 's', &vlen);
    int style = s_val ? scan_atoi(s_val, vlen) : -1;

    /* Determine cell type from t= attribute */
    XLSXCellType cell_type = XLSX_CELL_NUMBER;
    bool is_inline_str = false;
    if (type_first) {
        switch (type_first) {
        case 's':
            if (t_len == 1) cell_type = XLSX_CELL_SHARED_STRING;
            else if (t_len == 3 && t_val[1] == 't' && t_val[2] == 'r')
                cell_type = XLSX_CELL_STRING;
            break;
        case 'i':
            if (t_len >= 9) { cell_type = XLSX_CELL_STRING; is_inline_str = true; }
            break;
        case 'b':
            if (t_len == 1) cell_type = XLSX_CELL_BOOLEAN;
            break;
        case 'e':
            if (t_len == 1) cell_type = XLSX_CELL_ERROR;
            break;
        default: break;
        }
    }

    /* Extract value text from cell body.
     * Optimized: <v> content is always plain text (no nested elements), so
     * we use single memchr('<') to find both <v> open and value end. */
    const char *body = gt + 1;
    size_t body_len = (size_t)(close_tag - body);
    const char *val_text = NULL;
    size_t val_len = 0;

    if (is_inline_str) {
        /* Find <t>...</t> (or <t xml:space="preserve">...</t>) within <is> */
        const char *t_open = xlsx_memfind(body, body_len, "<t>", 3);
        if (!t_open) t_open = xlsx_memfind(body, body_len, "<t ", 3);
        if (t_open) {
            const char *t_gt = (const char *)memchr(t_open, '>',
                                                    (size_t)(close_tag - t_open));
            if (t_gt) {
                const char *ts = t_gt + 1;
                /* Text ends at next '<' (which is </t>) */
                const char *te = (const char *)memchr(ts, '<',
                                                      (size_t)(close_tag - ts));
                if (te) { val_text = ts; val_len = (size_t)(te - ts); }
            }
        }
    } else {
        /* Find first '<' in body — this is <v> (or <f> for formulas) */
        const char *first_lt = (const char *)memchr(body, '<', body_len);
        if (first_lt && first_lt + 2 < close_tag &&
            first_lt[1] == 'v' && first_lt[2] == '>') {
            /* <v> found — value ends at next '<' (which is </v>) */
            const char *vs = first_lt + 3;
            const char *ve = (const char *)memchr(vs, '<',
                                                  (size_t)(close_tag - vs));
            if (ve) { val_text = vs; val_len = (size_t)(ve - vs); }
        } else if (first_lt && first_lt + 1 < close_tag && first_lt[1] == 'f') {
            /* <f>formula</f><v>value</v> — skip formula, find <v> */
            const char *v_open = xlsx_memfind(first_lt, (size_t)(close_tag - first_lt),
                                              "<v>", 3);
            if (v_open) {
                const char *vs = v_open + 3;
                const char *ve = (const char *)memchr(vs, '<',
                                                      (size_t)(close_tag - vs));
                if (ve) { val_text = vs; val_len = (size_t)(ve - vs); }
            }
        }
    }

    if (!val_text || val_len == 0) return;

    /* Parse cell value */
    double numeric_val = SV_missval;
    int ss_idx = -1;
    const char *inline_str = NULL;

    switch (cell_type) {
    case XLSX_CELL_SHARED_STRING:
        ss_idx = scan_atoi(val_text, val_len);
        numeric_val = (double)ss_idx;
        break;
    case XLSX_CELL_NUMBER: {
        /* Need null-terminated string for parser */
        char num_buf[64];
        size_t copy = val_len < 63 ? val_len : 63;
        memcpy(num_buf, val_text, copy);
        num_buf[copy] = '\0';
        double parsed;
        if (ctools_parse_double_fast(num_buf, (int)copy, &parsed, SV_missval))
            numeric_val = parsed;
        else
            numeric_val = strtod(num_buf, NULL);
        if (style >= 0 && style < ctx->num_styles && ctx->date_styles[style])
            cell_type = XLSX_CELL_DATE;
        break;
    }
    case XLSX_CELL_BOOLEAN:
        numeric_val = (val_text[0] == '1') ? 1.0 : 0.0;
        break;
    case XLSX_CELL_STRING: {
        char *str_buf = (char *)ctools_arena_alloc(&ctx->parse_arena, val_len + 1);
        if (!str_buf) return;
        memcpy(str_buf, val_text, val_len);
        str_buf[val_len] = '\0';
        xlsx_xml_decode_entities(str_buf, val_len);
        inline_str = str_buf;
        break;
    }
    default:
        return;
    }

    /* Write to column-major arrays */
    if (ctx->cm_active && col >= 0 && col < ctx->cm_num_cols &&
        row > 0 && (size_t)row <= ctx->cm_num_rows) {
        size_t row_idx = (size_t)(row - 1);
        ctx->cm_numeric[col][row_idx] = numeric_val;
        ctx->cm_types[col][row_idx] = (uint8_t)cell_type;
    }

    /* Row-major fallback (when dimension unknown) */
    if (!ctx->cm_active) {
        XLSXParsedRow *row_ptr = ensure_row(ctx, row);
        if (row_ptr) {
            XLSXCell *cell = add_cell(row_ptr, col, ctx);
            if (cell) {
                cell->col = col;
                cell->row = row;
                cell->type = cell_type;
                cell->style_idx = style;
                switch (cell_type) {
                case XLSX_CELL_SHARED_STRING:
                    cell->value.shared_string_idx = ss_idx;
                    break;
                case XLSX_CELL_NUMBER:
                case XLSX_CELL_DATE:
                    cell->value.number = numeric_val;
                    break;
                case XLSX_CELL_BOOLEAN:
                    cell->value.boolean = (numeric_val == 1.0);
                    break;
                case XLSX_CELL_STRING:
                    cell->value.inline_string = inline_str;
                    break;
                default:
                    break;
                }
            }
        }
    }

    /* Track extent */
    if (col > ctx->max_col) ctx->max_col = col;

    /* Inline type inference */
    if (!ctx->allstring && !(ctx->firstrow && row == ctx->min_row)) {
        CImportColumnInfo *cs = ensure_col_stat(ctx, col);
        if (cs) {
            switch (cell_type) {
            case XLSX_CELL_SHARED_STRING:
            case XLSX_CELL_STRING: {
                cs->type = CIMPORT_COL_STRING;
                size_t slen = 0;
                if (cell_type == XLSX_CELL_SHARED_STRING && ss_idx >= 0 &&
                    (uint32_t)ss_idx < ctx->num_shared_strings) {
                    if (ctx->shared_string_lengths)
                        slen = ctx->shared_string_lengths[ss_idx];
                    else if (ctx->shared_strings[ss_idx])
                        slen = strlen(ctx->shared_strings[ss_idx]);
                } else if (inline_str) {
                    slen = strlen(inline_str);
                }
                if ((int)slen > cs->max_strlen) cs->max_strlen = (int)slen;
                break;
            }
            case XLSX_CELL_NUMBER:
            case XLSX_CELL_DATE:
                if (cs->type == CIMPORT_COL_UNKNOWN) cs->type = CIMPORT_COL_NUMERIC;
                if (cs->type == CIMPORT_COL_NUMERIC) {
                    double val = numeric_val;
                    if (cell_type == XLSX_CELL_DATE) val = excel_date_to_stata(numeric_val);
                    if (val < cs->min_value) cs->min_value = val;
                    if (val > cs->max_value) cs->max_value = val;
                    if (cs->is_integer && val != floor(val)) cs->is_integer = false;
                }
                break;
            case XLSX_CELL_BOOLEAN:
                if (cs->type == CIMPORT_COL_UNKNOWN) cs->type = CIMPORT_COL_NUMERIC;
                break;
            default:
                break;
            }
        }
    }
}

/* Scan a complete in-memory XML buffer for worksheet data.
 * buf must be NUL-terminated. */
static void xlsx_scan_buffer(XLSXContext *ctx, const char *buf, size_t len)
{
    const char *pos = buf;
    const char *end = buf + len;
    bool in_sheetdata = false;
    int current_row = 0;

    while (pos < end) {
        const char *lt = (const char *)memchr(pos, '<', (size_t)(end - pos));
        if (!lt || lt + 1 >= end) break;

        if (!in_sheetdata) {
            /* ---- Pre-sheetData preamble ---- */
            if (lt[1] == 'd' && lt + 10 <= end &&
                memcmp(lt + 1, "dimension", 9) == 0) {
                const char *gt = (const char *)memchr(lt, '>',
                                                      (size_t)(end - lt));
                if (!gt) break;
                const char *ra = xlsx_memfind(lt, (size_t)(gt - lt + 1),
                                              " ref=\"", 6);
                if (ra) {
                    const char *rv = ra + 6;
                    const char *rq = (const char *)memchr(rv, '"',
                                                          (size_t)(gt - rv));
                    if (rq) {
                        size_t rlen = (size_t)(rq - rv);
                        char ref_buf[64];
                        if (rlen < sizeof(ref_buf)) {
                            memcpy(ref_buf, rv, rlen);
                            ref_buf[rlen] = '\0';
                            handle_dimension_ref(ctx, ref_buf);
                        }
                    }
                }
                pos = gt + 1;
                continue;
            }
            if (lt[1] == 's' && lt + 10 <= end &&
                memcmp(lt + 1, "sheetData", 9) == 0) {
                const char *gt = (const char *)memchr(lt, '>',
                                                      (size_t)(end - lt));
                if (!gt) break;
                if (gt > lt && gt[-1] == '/') return;  /* Empty sheet */
                in_sheetdata = true;
                pos = gt + 1;
                continue;
            }
            const char *gt = (const char *)memchr(lt + 1, '>',
                                                  (size_t)(end - lt - 1));
            if (!gt) break;
            pos = gt + 1;
            continue;
        }

        /* ---- Inside <sheetData> ---- */
        char c1 = lt[1];

        if (c1 == 'c' && lt + 2 < end &&
            (lt[2] == ' ' || lt[2] == '>' || lt[2] == '/')) {
            const char *close_c = xlsx_memfind(lt + 2,
                                               (size_t)(end - lt - 2),
                                               "</c>", 4);
            if (!close_c) break;
            xlsx_scan_cell(ctx, lt, close_c, current_row);
            pos = close_c + 4;
        }
        else if (c1 == 'r' && lt + 4 <= end &&
                 lt[2] == 'o' && lt[3] == 'w' &&
                 (lt[4] == ' ' || lt[4] == '>' || lt[4] == '/')) {
            const char *gt = (const char *)memchr(lt, '>',
                                                  (size_t)(end - lt));
            if (!gt) break;

            size_t rlen;
            const char *r_val = scan_find_attr(lt + 4, gt, 'r', &rlen);
            current_row = r_val ? scan_atoi(r_val, rlen) : (current_row + 1);

            XLSXCellRange *range = &ctx->cell_range;
            if (range->end_row > 0 && current_row > range->end_row)
                return;
            if (range->start_row > 0 && current_row < range->start_row) {
                const char *row_close = xlsx_memfind(gt + 1,
                                                    (size_t)(end - gt - 1),
                                                    "</row>", 6);
                if (!row_close) break;
                pos = row_close + 6;
                continue;
            }

            if (current_row < ctx->min_row) ctx->min_row = current_row;
            if (current_row > ctx->max_row) ctx->max_row = current_row;
            if (!ctx->cm_active) ensure_row(ctx, current_row);

            pos = gt + 1;
        }
        else if (c1 == '/') {
            if (lt + 5 < end && lt[2] == 'r' && lt[3] == 'o' &&
                lt[4] == 'w' && lt[5] == '>') {
                pos = lt + 6;
            }
            else if (lt + 12 <= end &&
                     memcmp(lt + 2, "sheetData>", 10) == 0) {
                return;
            }
            else {
                const char *gt = (const char *)memchr(lt + 1, '>',
                                                      (size_t)(end - lt - 1));
                if (!gt) break;
                pos = gt + 1;
            }
        }
        else {
            const char *gt = (const char *)memchr(lt + 1, '>',
                                                  (size_t)(end - lt - 1));
            if (!gt) break;
            pos = gt + 1;
        }
    }
}

/* Main fast worksheet scanner.
 * Extracts worksheet via libdeflate (fast inflate) then scans in-memory.
 * Returns 0 on success, non-zero on error (caller falls back to SAX parser). */
static ST_retcode xlsx_scan_worksheet_fast(XLSXContext *ctx,
                                           xlsx_zip_archive *zip,
                                           size_t file_idx)
{
    size_t xml_size = 0;
    void *xml_data = xlsx_zip_extract_fast(zip, file_idx, &xml_size);
    if (!xml_data) return 610;

    /* NUL-terminate for safety */
    char *xml = (char *)xml_data;
    /* The extract_fast function allocates +1 byte, so this is safe */

    xlsx_scan_buffer(ctx, xml, xml_size);

    free(xml_data);
    return 0;
}

/* ---- SAX-based worksheet parser (fallback) ---- */

static bool worksheet_callback(const xlsx_xml_event *event, void *user_data)
{
    WorksheetParseState *state = (WorksheetParseState *)user_data;
    XLSXContext *ctx = state->ctx;
    const char *tag = event->tag_name;

    if (event->type == XLSX_XML_START_ELEMENT) {
        /* Tag dispatch via first character for common XLSX worksheet tags.
         * The worksheet tag set is small and fixed — this eliminates most
         * strcmp overhead on the ~40M events for a 10M cell dataset. */
        switch (tag[0]) {
        case 'c':
            if (tag[1] == '\0' && state->in_row) {
                /* <c> cell element */
                state->in_cell = true;
                state->value_len = 0;

                const char *r = xlsx_xml_get_attr(event, "r");
                const char *t = xlsx_xml_get_attr(event, "t");
                const char *s = xlsx_xml_get_attr(event, "s");

                if (r) {
                    int col, row;
                    if (xlsx_xml_parse_cell_ref(r, &col, &row)) {
                        state->current_col = col;
                        state->current_row = row;
                    }
                }

                XLSXCellRange *range = &ctx->cell_range;
                if (range->start_col >= 0 && state->current_col < range->start_col) {
                    state->in_cell = false;
                }
                if (range->end_col >= 0 && state->current_col > range->end_col) {
                    state->in_cell = false;
                }

                state->current_type = XLSX_CELL_NUMBER;
                if (t) {
                    switch (t[0]) {
                    case 's':
                        state->current_type = (t[1] == '\0') ? XLSX_CELL_SHARED_STRING :
                            (strcmp(t, "str") == 0) ? XLSX_CELL_STRING : XLSX_CELL_NUMBER;
                        break;
                    case 'i':
                        if (strcmp(t, "inlineStr") == 0) state->current_type = XLSX_CELL_STRING;
                        break;
                    case 'b':
                        if (t[1] == '\0') state->current_type = XLSX_CELL_BOOLEAN;
                        break;
                    case 'e':
                        if (t[1] == '\0') state->current_type = XLSX_CELL_ERROR;
                        break;
                    default:
                        break;
                    }
                }

                state->current_style = s ? xlsx_fast_atoi(s) : -1;
            }
            break;
        case 'v':
            if (tag[1] == '\0' && state->in_cell) {
                state->in_value = true;
                state->value_len = 0;
            }
            break;
        case 't':
            if (tag[1] == '\0' && state->in_inline_str) {
                state->in_t = true;
                state->value_len = 0;
            }
            break;
        case 'r':
            if (tag[1] == 'o' && tag[2] == 'w' && tag[3] == '\0' && state->in_sheetdata) {
                state->in_row = true;
                const char *r = xlsx_xml_get_attr(event, "r");
                state->current_row = r ? xlsx_fast_atoi(r) : (state->current_row + 1);

                XLSXCellRange *range = &ctx->cell_range;
                if (range->start_row > 0 && state->current_row < range->start_row) {
                    state->in_row = false;
                }
                if (range->end_row > 0 && state->current_row > range->end_row) {
                    return false;
                }

                if (state->in_row) {
                    /* Track min/max row for column-major path */
                    if (state->current_row < ctx->min_row) ctx->min_row = state->current_row;
                    if (state->current_row > ctx->max_row) ctx->max_row = state->current_row;

                    /* Only allocate row storage for row-major fallback path */
                    if (!ctx->cm_active) {
                        state->row_ptr = ensure_row(ctx, state->current_row);
                    }
                }
            }
            break;
        case 's':
            if (strcmp(tag, "sheetData") == 0) {
                state->in_sheetdata = true;
            }
            break;
        case 'i':
            if (tag[1] == 's' && tag[2] == '\0' && state->in_cell) {
                state->in_inline_str = true;
            }
            break;
        case 'd':
            if (strcmp(tag, "dimension") == 0) {
                const char *ref = xlsx_xml_get_attr(event, "ref");
                if (ref) handle_dimension_ref(ctx, ref);
            }
            break;
        default:
            break;
        }
    }
    else if (event->type == XLSX_XML_TEXT) {
        if (state->in_value || state->in_t) {
            size_t copy_len = event->text_len;
            if (state->value_len + copy_len >= sizeof(state->value_buf)) {
                copy_len = sizeof(state->value_buf) - state->value_len - 1;
            }
            if (copy_len > 0) {
                memcpy(state->value_buf + state->value_len, event->text, copy_len);
                state->value_len += copy_len;
            }
        }
    }
    else if (event->type == XLSX_XML_END_ELEMENT) {
        switch (tag[0]) {
        case 'c':
            if (tag[1] == '\0' && state->in_cell) {
                state->in_cell = false;
                state->value_buf[state->value_len] = '\0';

                if (state->value_len == 0) break;

                /* Parse cell value */
                XLSXCellType cell_type = state->current_type;
                double numeric_val = SV_missval;
                int ss_idx = -1;
                const char *inline_str = NULL;

                switch (cell_type) {
                case XLSX_CELL_SHARED_STRING:
                    ss_idx = xlsx_fast_atoi(state->value_buf);
                    numeric_val = (double)ss_idx;
                    break;
                case XLSX_CELL_NUMBER: {
                    double parsed_val;
                    if (ctools_parse_double_fast(state->value_buf, (int)state->value_len, &parsed_val, SV_missval)) {
                        numeric_val = parsed_val;
                    } else {
                        numeric_val = strtod(state->value_buf, NULL);
                    }
                    if (state->current_style >= 0 &&
                        state->current_style < ctx->num_styles &&
                        ctx->date_styles[state->current_style]) {
                        cell_type = XLSX_CELL_DATE;
                    }
                    break;
                }
                case XLSX_CELL_BOOLEAN:
                    numeric_val = (state->value_buf[0] == '1') ? 1.0 : 0.0;
                    break;
                case XLSX_CELL_STRING:
                    xlsx_xml_decode_entities(state->value_buf, state->value_len);
                    inline_str = ctools_arena_strdup(&ctx->parse_arena, state->value_buf);
                    if (!inline_str) return false;
                    break;
                default:
                    break;
                }

                int col = state->current_col;
                int row = state->current_row;

                /* Column-major direct write path */
                if (ctx->cm_active && col >= 0 && col < ctx->cm_num_cols &&
                    row > 0 && (size_t)row <= ctx->cm_num_rows) {
                    size_t row_idx = (size_t)(row - 1);
                    ctx->cm_numeric[col][row_idx] = numeric_val;
                    ctx->cm_types[col][row_idx] = (uint8_t)cell_type;
                }

                /* Row-major fallback path (when dimension unknown) */
                if (!ctx->cm_active && state->row_ptr) {
                    XLSXCell *cell = add_cell(state->row_ptr, col, ctx);
                    if (cell) {
                        cell->col = col;
                        cell->row = row;
                        cell->type = cell_type;
                        cell->style_idx = state->current_style;
                        switch (cell_type) {
                        case XLSX_CELL_SHARED_STRING:
                            cell->value.shared_string_idx = ss_idx;
                            break;
                        case XLSX_CELL_NUMBER:
                        case XLSX_CELL_DATE:
                            cell->value.number = numeric_val;
                            break;
                        case XLSX_CELL_BOOLEAN:
                            cell->value.boolean = (numeric_val == 1.0);
                            break;
                        case XLSX_CELL_STRING:
                            cell->value.inline_string = inline_str;
                            break;
                        default:
                            break;
                        }
                    }
                }

                /* Track extent */
                if (col > ctx->max_col) ctx->max_col = col;

                /* Inline type inference */
                if (!ctx->allstring && !(ctx->firstrow && row == ctx->min_row)) {
                    CImportColumnInfo *cs = ensure_col_stat(ctx, col);
                    if (cs) {
                        switch (cell_type) {
                        case XLSX_CELL_SHARED_STRING:
                        case XLSX_CELL_STRING:
                            cs->type = CIMPORT_COL_STRING;
                            {
                                size_t slen = 0;
                                if (cell_type == XLSX_CELL_SHARED_STRING && ss_idx >= 0 &&
                                    (uint32_t)ss_idx < ctx->num_shared_strings) {
                                    if (ctx->shared_string_lengths) {
                                        slen = ctx->shared_string_lengths[ss_idx];
                                    } else if (ctx->shared_strings[ss_idx]) {
                                        slen = strlen(ctx->shared_strings[ss_idx]);
                                    }
                                } else if (cell_type == XLSX_CELL_STRING && inline_str) {
                                    slen = strlen(inline_str);
                                }
                                if ((int)slen > cs->max_strlen) {
                                    cs->max_strlen = (int)slen;
                                }
                            }
                            break;
                        case XLSX_CELL_NUMBER:
                        case XLSX_CELL_DATE:
                            if (cs->type == CIMPORT_COL_UNKNOWN) {
                                cs->type = CIMPORT_COL_NUMERIC;
                            }
                            if (cs->type == CIMPORT_COL_NUMERIC) {
                                double val = numeric_val;
                                if (cell_type == XLSX_CELL_DATE) {
                                    val = excel_date_to_stata(numeric_val);
                                }
                                if (val < cs->min_value) cs->min_value = val;
                                if (val > cs->max_value) cs->max_value = val;
                                if (cs->is_integer && val != floor(val)) {
                                    cs->is_integer = false;
                                }
                            }
                            break;
                        case XLSX_CELL_BOOLEAN:
                            if (cs->type == CIMPORT_COL_UNKNOWN) {
                                cs->type = CIMPORT_COL_NUMERIC;
                            }
                            break;
                        default:
                            break;
                        }
                    }
                }
            }
            break;
        case 'v':
            if (tag[1] == '\0') state->in_value = false;
            break;
        case 't':
            if (tag[1] == '\0') state->in_t = false;
            break;
        case 'r':
            if (tag[1] == 'o' && tag[2] == 'w' && tag[3] == '\0') {
                state->in_row = false;
                state->row_ptr = NULL;
            }
            break;
        case 's':
            if (strcmp(tag, "sheetData") == 0) {
                state->in_sheetdata = false;
            }
            break;
        case 'i':
            if (tag[1] == 's' && tag[2] == '\0') {
                state->in_inline_str = false;
            }
            break;
        default:
            break;
        }
    }

    return true;
}

ST_retcode xlsx_parse_worksheet(XLSXContext *ctx)
{
    xlsx_zip_archive *zip = (xlsx_zip_archive *)ctx->zip_archive;

    double t_start = ctools_timer_seconds();

    /* Build worksheet path */
    char worksheet_path[128];
    int sheet_idx = ctx->sheets[ctx->selected_sheet].sheet_index;
    snprintf(worksheet_path, sizeof(worksheet_path), XLSX_WORKSHEET_FMT, sheet_idx);

    size_t file_idx = xlsx_zip_locate_file(zip, worksheet_path);
    if (file_idx == (size_t)-1) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to find worksheet: %s", worksheet_path);
        return 610;
    }

    /* Try fast memchr-based scanner first (streaming) */
    ST_retcode rc = xlsx_scan_worksheet_fast(ctx, zip, file_idx);
    if (rc == 0) {
        ctx->time_worksheet = ctools_timer_seconds() - t_start;
        return 0;
    }

    /* Fast scanner unavailable — fall back to SAX parser */
    WorksheetParseState state = {0};
    state.ctx = ctx;

    xlsx_xml_parser *parser = xlsx_xml_parser_create(worksheet_callback, &state);
    if (!parser) {
        return 920;
    }

    /* Full extraction fallback (streaming was unavailable) */
    size_t size;
    void *data = xlsx_zip_extract_to_heap(zip, file_idx, &size);
    if (!data) {
        xlsx_xml_parser_destroy(parser);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to extract worksheet: %s", worksheet_path);
        return 610;
    }
    xlsx_xml_parse(parser, (const char *)data, size, true);
    free(data);

    xlsx_xml_parser_destroy(parser);

    ctx->time_worksheet = ctools_timer_seconds() - t_start;

    return 0;
}

/* ============================================================================
 * Type Inference
 * ============================================================================ */

ST_retcode xlsx_infer_types(XLSXContext *ctx)
{
    double t_start = ctools_timer_seconds();

    /* Determine number of columns */
    ctx->num_columns = ctx->max_col + 1;

    /* Apply column range */
    int start_col = ctx->cell_range.start_col >= 0 ? ctx->cell_range.start_col : 0;
    int end_col = ctx->cell_range.end_col >= 0 ? ctx->cell_range.end_col : ctx->max_col;
    ctx->num_columns = end_col - start_col + 1;

    if (ctx->num_columns <= 0) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "No columns found in worksheet");
        return 2000;
    }

    /* Allocate columns */
    ctx->columns = (CImportColumnInfo *)calloc(ctx->num_columns, sizeof(CImportColumnInfo));
    if (!ctx->columns) {
        return 920;
    }

    /* Copy pre-computed stats from col_stats (populated during worksheet parse) */
    for (int c = 0; c < ctx->num_columns; c++) {
        int src_col = c + start_col;
        if (ctx->col_stats && src_col < ctx->col_stats_capacity) {
            ctx->columns[c] = ctx->col_stats[src_col];
        } else {
            ctx->columns[c].type = CIMPORT_COL_UNKNOWN;
            ctx->columns[c].num_subtype = CIMPORT_NUM_BYTE;
            ctx->columns[c].is_integer = true;
            ctx->columns[c].min_value = DBL_MAX;
            ctx->columns[c].max_value = -DBL_MAX;
            ctx->columns[c].max_strlen = 0;
        }
    }

    /* Handle allstring mode */
    if (ctx->allstring) {
        for (int c = 0; c < ctx->num_columns; c++) {
            ctx->columns[c].type = CIMPORT_COL_STRING;
        }
    }

    /* Determine header row */
    int header_row = -1;

    if (ctx->firstrow) {
        if (ctx->cm_active || ctx->num_rows > 0) {
            header_row = ctx->min_row;
        }
    }

    /* Extract column names from header row */
    if (header_row > 0) {
        if (ctx->cm_active && header_row >= 1 && (size_t)header_row <= ctx->cm_num_rows) {
            /* Column-major path: read header directly from cm arrays */
            size_t hdr_idx = (size_t)(header_row - 1);
            for (int c = 0; c < ctx->num_columns; c++) {
                int cm_col = c + start_col;
                if (cm_col < 0 || cm_col >= ctx->cm_num_cols) continue;
                uint8_t ctype = ctx->cm_types[cm_col][hdr_idx];
                const char *header_value = NULL;
                if (ctype == XLSX_CELL_SHARED_STRING) {
                    int ss_idx = (int)ctx->cm_numeric[cm_col][hdr_idx];
                    if (ss_idx >= 0 && (uint32_t)ss_idx < ctx->num_shared_strings) {
                        header_value = ctx->shared_strings[ss_idx];
                    }
                }
                /* Inline strings stored in parse_arena are not in cm_numeric;
                 * they're rare in headers so we skip them (fallback to default name) */
                generate_varname(ctx->columns[c].name, c, header_value, ctx->case_mode);
            }
        } else {
            /* Row-major fallback path */
            for (int i = 0; i < ctx->num_rows; i++) {
                if (ctx->rows[i].row_num == header_row) {
                    for (int j = 0; j < ctx->rows[i].num_cells; j++) {
                        XLSXCell *cell = &ctx->rows[i].cells[j];
                        int col_idx = cell->col - start_col;
                        if (col_idx < 0 || col_idx >= ctx->num_columns) continue;

                        const char *header_value = NULL;
                        if (cell->type == XLSX_CELL_SHARED_STRING &&
                            (uint32_t)cell->value.shared_string_idx < ctx->num_shared_strings) {
                            header_value = ctx->shared_strings[cell->value.shared_string_idx];
                        } else if (cell->type == XLSX_CELL_STRING) {
                            header_value = cell->value.inline_string;
                        }

                        generate_varname(ctx->columns[col_idx].name, col_idx,
                                         header_value, ctx->case_mode);
                    }
                    break;
                }
            }
        }
    }

    /* Generate default names for columns without headers */
    for (int c = 0; c < ctx->num_columns; c++) {
        if (ctx->columns[c].name[0] == '\0') {
            generate_varname(ctx->columns[c].name, c, NULL, ctx->case_mode);
        }
    }

    /* Finalize numeric subtypes (O(num_columns) — no cell iteration) */
    for (int c = 0; c < ctx->num_columns; c++) {
        CImportColumnInfo *col = &ctx->columns[c];

        if (col->type == CIMPORT_COL_UNKNOWN) {
            col->type = CIMPORT_COL_NUMERIC;  /* Empty column defaults to numeric */
        }

        if (col->type == CIMPORT_COL_STRING) {
            if (col->max_strlen < 1) col->max_strlen = 1;
            if (col->max_strlen > 2045) col->max_strlen = 2045;
        } else if (col->type == CIMPORT_COL_NUMERIC) {
            if (col->is_integer) {
                if (col->min_value >= -127 && col->max_value <= 100) {
                    col->num_subtype = CIMPORT_NUM_BYTE;
                } else if (col->min_value >= -32767 && col->max_value <= 32740) {
                    col->num_subtype = CIMPORT_NUM_INT;
                } else if (col->min_value >= -2147483647.0 && col->max_value <= 2147483620.0) {
                    col->num_subtype = CIMPORT_NUM_LONG;
                } else {
                    col->num_subtype = CIMPORT_NUM_DOUBLE;
                }
            } else {
                col->num_subtype = CIMPORT_NUM_DOUBLE;
            }
        }
    }

    ctx->time_type_infer = ctools_timer_seconds() - t_start;

    return 0;
}

/* ============================================================================
 * Cache Building
 * ============================================================================ */

/* Column-major fast path: build cache directly from cm_numeric/cm_types arrays.
 * For numeric columns, transfers ownership of the array (zero copy).
 * For string columns, resolves shared string indices in one pass. */
static ST_retcode xlsx_build_cache_columnar(XLSXContext *ctx)
{
    double t_start = ctools_timer_seconds();

    int start_col = ctx->cell_range.start_col >= 0 ? ctx->cell_range.start_col : 0;
    int data_start_row = ctx->min_row;
    if (ctx->firstrow) data_start_row = ctx->min_row + 1;

    int range_start_row = ctx->cell_range.start_row > 0 ? ctx->cell_range.start_row : data_start_row;
    int range_end_row = ctx->cell_range.end_row > 0 ? ctx->cell_range.end_row : ctx->max_row;

    int actual_data_start = data_start_row;
    if (actual_data_start < range_start_row) actual_data_start = range_start_row;

    /* cm arrays are 0-indexed (row 1 = index 0), convert to cm indices */
    size_t cm_start = (actual_data_start >= 1) ? (size_t)(actual_data_start - 1) : 0;
    size_t cm_end = (range_end_row >= 1) ? (size_t)(range_end_row - 1) : 0;
    if (cm_end >= ctx->cm_num_rows) cm_end = ctx->cm_num_rows - 1;

    size_t total_data_rows = (cm_end >= cm_start) ? (cm_end - cm_start + 1) : 0;
    if (total_data_rows == 0) {
        snprintf(ctx->error_message, sizeof(ctx->error_message), "No data rows found");
        return 2000;
    }

    ctx->col_cache = (CImportColumnCache *)calloc((size_t)ctx->num_columns,
                                                   sizeof(CImportColumnCache));
    if (!ctx->col_cache) return 920;

    for (int c = 0; c < ctx->num_columns; c++) {
        CImportColumnInfo *col = &ctx->columns[c];
        CImportColumnCache *cache = &ctx->col_cache[c];
        cache->count = total_data_rows;
        int cm_col = c + start_col;

        if (cm_col < 0 || cm_col >= ctx->cm_num_cols) {
            /* Column out of range — fill with missing */
            if (col->type == CIMPORT_COL_STRING) {
                cache->string_data = (char **)calloc(total_data_rows, sizeof(char *));
                if (!cache->string_data) return 920;
                ctools_arena_init(&cache->string_arena, 0);
            } else {
                cache->numeric_data = (double *)malloc(total_data_rows * sizeof(double));
                if (!cache->numeric_data) return 920;
                for (size_t i = 0; i < total_data_rows; i++)
                    cache->numeric_data[i] = SV_missval;
            }
            continue;
        }

        double *src_num = ctx->cm_numeric[cm_col];
        uint8_t *src_types = ctx->cm_types[cm_col];

        if (col->type == CIMPORT_COL_STRING) {
            /* String column: resolve shared string indices and inline strings */
            cache->string_data = (char **)calloc(total_data_rows, sizeof(char *));
            if (!cache->string_data) return 920;
            ctools_arena_init(&cache->string_arena, 0);

            for (size_t r = 0; r < total_data_rows; r++) {
                size_t cm_r = cm_start + r;
                uint8_t ct = src_types[cm_r];
                const char *str = NULL;
                size_t slen = 0;

                if (ct == XLSX_CELL_SHARED_STRING) {
                    int ss_idx = (int)src_num[cm_r];
                    if (ss_idx >= 0 && (uint32_t)ss_idx < ctx->num_shared_strings) {
                        str = ctx->shared_strings[ss_idx];
                        if (ctx->shared_string_lengths) {
                            slen = ctx->shared_string_lengths[ss_idx];
                        } else if (str) {
                            slen = strlen(str);
                        }
                    }
                } else if (ct == XLSX_CELL_NUMBER || ct == XLSX_CELL_DATE) {
                    /* Number in string column — convert */
                    char num_buf[64];
                    snprintf(num_buf, sizeof(num_buf), "%g", src_num[cm_r]);
                    str = num_buf;
                    slen = strlen(str);
                } else if (ct == XLSX_CELL_BOOLEAN) {
                    str = (src_num[cm_r] == 1.0) ? "1" : "0";
                    slen = 1;
                }
                /* Note: inline strings (XLSX_CELL_STRING) are stored in parse_arena
                 * but we can't easily recover them from cm_numeric. They're rare in
                 * practice. For now they show as empty in the string column. */

                if (str && slen > 0) {
                    char *copy = (char *)ctools_arena_alloc(&cache->string_arena, slen + 1);
                    if (copy) {
                        memcpy(copy, str, slen + 1);
                        cache->string_data[r] = copy;
                    }
                }
            }
        } else {
            /* Numeric column: check if we can do zero-copy transfer */
            if (cm_start == 0 && total_data_rows == ctx->cm_num_rows) {
                /* Perfect alignment: transfer ownership of the array */
                cache->numeric_data = src_num;
                ctx->cm_numeric[cm_col] = NULL; /* prevent double-free */

                /* Apply date conversion in-place */
                for (size_t r = 0; r < total_data_rows; r++) {
                    if (src_types[r] == XLSX_CELL_DATE) {
                        cache->numeric_data[r] = excel_date_to_stata(cache->numeric_data[r]);
                    }
                }
            } else {
                /* Subrange: copy the relevant slice */
                cache->numeric_data = (double *)malloc(total_data_rows * sizeof(double));
                if (!cache->numeric_data) return 920;
                for (size_t r = 0; r < total_data_rows; r++) {
                    size_t cm_r = cm_start + r;
                    uint8_t ct = src_types[cm_r];
                    if (ct == XLSX_CELL_NUMBER) {
                        cache->numeric_data[r] = src_num[cm_r];
                    } else if (ct == XLSX_CELL_DATE) {
                        cache->numeric_data[r] = excel_date_to_stata(src_num[cm_r]);
                    } else if (ct == XLSX_CELL_BOOLEAN) {
                        cache->numeric_data[r] = src_num[cm_r];
                    } else {
                        cache->numeric_data[r] = SV_missval;
                    }
                }
            }
        }
    }

    ctx->cache_ready = true;
    ctx->time_cache = ctools_timer_seconds() - t_start;
    return 0;
}

/* Row-major fallback: build cache from XLSXParsedRow arrays */
static ST_retcode xlsx_build_cache_rowmajor(XLSXContext *ctx)
{
    double t_start = ctools_timer_seconds();

    int start_col = ctx->cell_range.start_col >= 0 ? ctx->cell_range.start_col : 0;
    int data_start_row = ctx->min_row;
    if (ctx->firstrow) data_start_row = ctx->min_row + 1;

    int range_start_row = ctx->cell_range.start_row > 0 ? ctx->cell_range.start_row : data_start_row;
    int range_end_row = ctx->cell_range.end_row > 0 ? ctx->cell_range.end_row : ctx->max_row;

    /* Count empty rows within range */
    size_t total_data_rows = 0;
    if (range_end_row >= range_start_row) {
        total_data_rows = (size_t)(range_end_row - range_start_row + 1);
        if (ctx->firstrow && range_start_row == ctx->min_row) {
            total_data_rows--;
        }
    }

    if (total_data_rows == 0) {
        snprintf(ctx->error_message, sizeof(ctx->error_message), "No data rows found");
        return 2000;
    }

    ctx->col_cache = (CImportColumnCache *)calloc((size_t)ctx->num_columns,
                                                   sizeof(CImportColumnCache));
    if (!ctx->col_cache) return 920;

    for (int c = 0; c < ctx->num_columns; c++) {
        CImportColumnInfo *col = &ctx->columns[c];
        CImportColumnCache *cache = &ctx->col_cache[c];
        cache->count = total_data_rows;

        if (col->type == CIMPORT_COL_STRING) {
            cache->string_data = (char **)calloc(total_data_rows, sizeof(char *));
            if (!cache->string_data) return 920;
            ctools_arena_init(&cache->string_arena, 0);
        } else {
            cache->numeric_data = (double *)malloc(total_data_rows * sizeof(double));
            if (!cache->numeric_data) return 920;
            for (size_t i = 0; i < total_data_rows; i++) {
                cache->numeric_data[i] = SV_missval;
            }
        }
    }

    int actual_data_start = data_start_row;
    if (actual_data_start < range_start_row) actual_data_start = range_start_row;

    for (int i = 0; i < ctx->num_rows; i++) {
        XLSXParsedRow *row = &ctx->rows[i];
        if (row->row_num < actual_data_start) continue;
        if (row->row_num > range_end_row) continue;

        size_t cache_row = (size_t)(row->row_num - actual_data_start);
        if (cache_row >= total_data_rows) continue;

        for (int j = 0; j < row->num_cells; j++) {
            XLSXCell *cell = &row->cells[j];
            int col_idx = cell->col - start_col;
            if (col_idx < 0 || col_idx >= ctx->num_columns) continue;

            CImportColumnInfo *col = &ctx->columns[col_idx];
            CImportColumnCache *cache = &ctx->col_cache[col_idx];

            if (col->type == CIMPORT_COL_STRING) {
                const char *str = "";
                if (cell->type == XLSX_CELL_SHARED_STRING &&
                    (uint32_t)cell->value.shared_string_idx < ctx->num_shared_strings) {
                    str = ctx->shared_strings[cell->value.shared_string_idx];
                } else if (cell->type == XLSX_CELL_STRING) {
                    str = cell->value.inline_string ? cell->value.inline_string : "";
                } else if (cell->type == XLSX_CELL_NUMBER ||
                           cell->type == XLSX_CELL_DATE) {
                    char num_buf[64];
                    snprintf(num_buf, sizeof(num_buf), "%g", cell->value.number);
                    str = num_buf;
                }
                size_t len = strlen(str) + 1;
                char *copy = (char *)ctools_arena_alloc(&cache->string_arena, len);
                if (copy) {
                    memcpy(copy, str, len);
                    cache->string_data[cache_row] = copy;
                }
            } else {
                double val = SV_missval;
                switch (cell->type) {
                case XLSX_CELL_NUMBER:
                    val = cell->value.number;
                    break;
                case XLSX_CELL_DATE:
                    val = excel_date_to_stata(cell->value.number);
                    break;
                case XLSX_CELL_BOOLEAN:
                    val = cell->value.boolean ? 1.0 : 0.0;
                    break;
                default:
                    break;
                }
                cache->numeric_data[cache_row] = val;
            }
        }
    }

    ctx->cache_ready = true;
    ctx->time_cache = ctools_timer_seconds() - t_start;
    return 0;
}

ST_retcode xlsx_build_cache(XLSXContext *ctx)
{
    if (ctx->cm_active) {
        return xlsx_build_cache_columnar(ctx);
    }
    return xlsx_build_cache_rowmajor(ctx);
}

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static void generate_varname(char *buf, int col, const char *header_value, int case_mode)
{
    if (header_value && header_value[0]) {
        /* Use header value, sanitizing for Stata variable name rules */
        const char *src = header_value;
        int len = 0;

        /* Skip leading non-alpha characters */
        while (*src && !isalpha((unsigned char)*src) && *src != '_') src++;

        /* First char must be letter or underscore */
        if (!*src) {
            xlsx_xml_index_to_col(col, buf);
            return;
        }

        while (*src && len < CTOOLS_MAX_VARNAME_LEN - 1) {
            char c = *src++;
            if (isalnum((unsigned char)c) || c == '_') {
                if (case_mode == 1) c = (char)tolower((unsigned char)c);
                else if (case_mode == 2) c = (char)toupper((unsigned char)c);
                buf[len++] = c;
            }
        }
        buf[len] = '\0';

        if (len == 0) {
            xlsx_xml_index_to_col(col, buf);
        }
    } else {
        /* Generate column letter name (A, B, ..., AA, AB, etc.) */
        xlsx_xml_index_to_col(col, buf);
    }
}

static double excel_date_to_stata(double excel_date)
{
    /* Excel dates are days since 1899-12-30 (with 1900 bug)
     * Stata dates are days since 1960-01-01
     * Difference: 21916 days (1960-01-01 - 1899-12-30)
     * But Excel has the 1900 bug (treats 1900 as leap year), so subtract 1 more
     * for dates after Feb 28, 1900 */

    /* Check for Excel's 1900 date system bug */
    if (excel_date >= 60.0) {
        /* After Feb 28, 1900 - adjust for bug */
        excel_date -= 1.0;
    }

    /* Convert to Stata date */
    return excel_date - 21916.0;
}

int xlsx_select_sheet_by_name(XLSXContext *ctx, const char *name)
{
    if (!ctx || !name) return -1;

    for (int i = 0; i < ctx->num_sheets; i++) {
        if (strcmp(ctx->sheets[i].name, name) == 0) {
            ctx->selected_sheet = i;
            return i;
        }
    }
    return -1;
}

bool xlsx_parse_cellrange(const char *range_str, XLSXCellRange *range)
{
    if (!range_str || !range) return false;

    /* Initialize to -1 (meaning "all") */
    range->start_col = -1;
    range->start_row = -1;
    range->end_col = -1;
    range->end_row = -1;

    /* Find colon separator */
    const char *colon = strchr(range_str, ':');

    /* Parse start reference */
    const char *start = range_str;
    const char *end = colon ? colon : range_str + strlen(range_str);

    if (start < end) {
        char start_ref[32];
        size_t len = end - start;
        if (len >= sizeof(start_ref)) len = sizeof(start_ref) - 1;
        memcpy(start_ref, start, len);
        start_ref[len] = '\0';

        int col, row;
        if (xlsx_xml_parse_cell_ref(start_ref, &col, &row)) {
            range->start_col = col;
            range->start_row = row;
        }
    }

    /* Parse end reference */
    if (colon) {
        const char *end_ref_str = colon + 1;
        if (*end_ref_str) {
            int col, row;
            if (xlsx_xml_parse_cell_ref(end_ref_str, &col, &row)) {
                range->end_col = col;
                range->end_row = row;
            }
        }
    }

    return true;
}

/* ============================================================================
 * Parallel XLSX Store Worker
 * ============================================================================ */

typedef struct {
    CImportColumnInfo *col;
    CImportColumnCache *cache;
    ST_int var;
    ST_int stata_nobs;
} xlsx_store_task;

static void *xlsx_store_worker(void *arg)
{
    xlsx_store_task *task = (xlsx_store_task *)arg;
    if (task->col->type == CIMPORT_COL_STRING) {
        for (size_t r = 0; r < task->cache->count && (ST_int)(r + 1) <= task->stata_nobs; r++) {
            char *str = task->cache->string_data[r] ? task->cache->string_data[r] : (char *)"";
            SF_sstore(task->var, (ST_int)(r + 1), str);
        }
    } else {
        for (size_t r = 0; r < task->cache->count && (ST_int)(r + 1) <= task->stata_nobs; r++) {
            SF_vstore(task->var, (ST_int)(r + 1), task->cache->numeric_data[r]);
        }
    }
    return NULL;
}

/* ============================================================================
 * Main Entry Point
 * ============================================================================ */

/* Helper: parse arguments into an XLSXContext (options only, not mode/filename) */
static void xlsx_parse_options(XLSXContext *ctx, const char *args)
{
    char args_copy[4096];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    /* Skip mode and filename tokens */
    strtok(args_copy, " \t");
    strtok(NULL, " \t");

    char *opt;
    while ((opt = strtok(NULL, " \t")) != NULL) {
        if (strncmp(opt, "cellrange=", 10) == 0) {
            xlsx_parse_cellrange(opt + 10, &ctx->cell_range);
        } else if (strcmp(opt, "firstrow") == 0) {
            ctx->firstrow = true;
        } else if (strcmp(opt, "allstring") == 0) {
            ctx->allstring = true;
        } else if (strncmp(opt, "case=", 5) == 0) {
            if (strcmp(opt + 5, "lower") == 0) ctx->case_mode = 1;
            else if (strcmp(opt + 5, "upper") == 0) ctx->case_mode = 2;
        } else if (strcmp(opt, "verbose") == 0) {
            ctx->verbose = true;
        }
    }
}

/* Helper: extract sheet name from args (returns NULL if not specified) */
static char *xlsx_extract_sheet_name(const char *args, char *buf, size_t buf_size)
{
    char args_copy[4096];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    char *p = strstr(args_copy, "sheet=");
    if (!p) return NULL;

    const char *start = p + 6;
    const char *end = start;
    while (*end && !isspace((unsigned char)*end)) end++;

    size_t len = (size_t)(end - start);
    if (len >= buf_size) len = buf_size - 1;
    memcpy(buf, start, len);
    buf[len] = '\0';
    return buf;
}

/* Worker structs for parallel shared strings + styles parsing */
typedef struct {
    XLSXContext *ctx;
    const void *data;
    size_t size;
    ST_retcode rc;
} xlsx_ss_parse_task;

typedef struct {
    XLSXContext *ctx;
    const void *data;
    size_t size;
    ST_retcode rc;
} xlsx_styles_parse_task;

static void *xlsx_ss_parse_worker(void *arg)
{
    xlsx_ss_parse_task *task = (xlsx_ss_parse_task *)arg;
    task->rc = xlsx_parse_shared_strings_buf(task->ctx, task->data, task->size);
    return NULL;
}

static void *xlsx_styles_parse_worker(void *arg)
{
    xlsx_styles_parse_task *task = (xlsx_styles_parse_task *)arg;
    task->rc = xlsx_parse_styles_buf(task->ctx, task->data, task->size);
    return NULL;
}

/* Helper: full parse pipeline (open, workbook, shared strings, styles, worksheet, infer) */
static ST_retcode xlsx_full_parse(XLSXContext *ctx, const char *filename,
                                   const char *sheet_name)
{
    ST_retcode rc;

    rc = xlsx_open_file(ctx, filename);
    if (rc) { SF_error(ctx->error_message); return rc; }

    rc = xlsx_parse_workbook(ctx);
    if (rc) { SF_error(ctx->error_message); return rc; }

    if (ctx->num_sheets == 0) {
        SF_error("No sheets found in workbook\n");
        return 610;
    }

    if (sheet_name) {
        int idx = xlsx_select_sheet_by_name(ctx, sheet_name);
        if (idx < 0) {
            char err_buf[256];
            snprintf(err_buf, sizeof(err_buf), "Sheet not found: %s\n", sheet_name);
            SF_error(err_buf);
            return 601;
        }
    }

    /* Parse shared strings and styles in parallel if both exist.
     * ZIP extraction is NOT thread-safe, so extract both files to heap
     * sequentially first, then parse the XML in parallel. */
    xlsx_zip_archive *zip = (xlsx_zip_archive *)ctx->zip_archive;
    double t_ss_start = ctools_timer_seconds();

    size_t ss_idx = xlsx_zip_locate_file(zip, XLSX_PATH_SHARED_STRINGS);
    size_t st_idx = xlsx_zip_locate_file(zip, XLSX_PATH_STYLES);

    void *ss_data = NULL;
    size_t ss_size = 0;
    void *st_data = NULL;
    size_t st_size = 0;

    /* Sequential extraction */
    if (ss_idx != (size_t)-1) {
        ss_data = xlsx_zip_extract_to_heap(zip, ss_idx, &ss_size);
        if (!ss_data) {
            snprintf(ctx->error_message, sizeof(ctx->error_message),
                     "Failed to extract sharedStrings.xml");
            return 610;
        }
    }
    if (st_idx != (size_t)-1) {
        st_data = xlsx_zip_extract_to_heap(zip, st_idx, &st_size);
        /* Non-fatal if styles extraction fails */
    }

    /* Parallel parse if both buffers exist and thread pool available */
    ctools_persistent_pool *pool = ctools_get_global_pool();
    if (ss_data && st_data && pool) {
        xlsx_ss_parse_task ss_task = { .ctx = ctx, .data = ss_data, .size = ss_size, .rc = 0 };
        xlsx_styles_parse_task st_task = { .ctx = ctx, .data = st_data, .size = st_size, .rc = 0 };

        ctools_persistent_pool_submit(pool, xlsx_ss_parse_worker, &ss_task);
        ctools_persistent_pool_submit(pool, xlsx_styles_parse_worker, &st_task);
        ctools_persistent_pool_wait(pool);

        rc = ss_task.rc;
        free(ss_data);
        free(st_data);

        ctx->time_shared_strings = ctools_timer_seconds() - t_ss_start;

        if (rc) { SF_error(ctx->error_message); return rc; }
    } else {
        /* Sequential fallback */
        if (ss_data) {
            rc = xlsx_parse_shared_strings_buf(ctx, ss_data, ss_size);
            free(ss_data);
            ctx->time_shared_strings = ctools_timer_seconds() - t_ss_start;
            if (rc) { SF_error(ctx->error_message); return rc; }
        } else {
            ctx->time_shared_strings = ctools_timer_seconds() - t_ss_start;
        }
        if (st_data) {
            xlsx_parse_styles_buf(ctx, st_data, st_size);
            free(st_data);
        }
    }

    rc = xlsx_parse_worksheet(ctx);
    if (rc) { SF_error(ctx->error_message); return rc; }

    rc = xlsx_infer_types(ctx);
    if (rc) { SF_error(ctx->error_message); return rc; }

    return 0;
}

/* Helper: compute nobs from context */
static size_t xlsx_compute_nobs(XLSXContext *ctx)
{
    int data_start = ctx->min_row;
    if (ctx->firstrow) data_start++;
    int range_start = ctx->cell_range.start_row > 0 ? ctx->cell_range.start_row : data_start;
    int range_end = ctx->cell_range.end_row > 0 ? ctx->cell_range.end_row : ctx->max_row;
    if (range_start < data_start) range_start = data_start;

    if (range_end >= range_start)
        return (size_t)(range_end - range_start + 1);
    return 0;
}

/* Helper: save scan metadata to Stata macros */
static void xlsx_save_scan_macros(XLSXContext *ctx, size_t nobs)
{
    char buf[32];

    snprintf(buf, sizeof(buf), "%zu", nobs);
    SF_macro_save("_cimport_nobs", buf);

    snprintf(buf, sizeof(buf), "%d", ctx->num_columns);
    SF_macro_save("_cimport_nvar", buf);

    char varnames[8192] = "";
    char vartypes[1024] = "";
    char numtypes[1024] = "";
    char strlens[1024] = "";
    size_t varnames_pos = 0, vartypes_pos = 0, numtypes_pos = 0, strlens_pos = 0;
    int written;

    for (int c = 0; c < ctx->num_columns; c++) {
        CImportColumnInfo *col = &ctx->columns[c];
        const char *sep = (c > 0) ? " " : "";

        if (varnames_pos < sizeof(varnames) - 1) {
            written = snprintf(varnames + varnames_pos, sizeof(varnames) - varnames_pos,
                               "%s%s", sep, col->name);
            if (written > 0) varnames_pos += (size_t)written;
        }
        if (vartypes_pos < sizeof(vartypes) - 1) {
            written = snprintf(vartypes + vartypes_pos, sizeof(vartypes) - vartypes_pos,
                               "%s%d", sep, col->type == CIMPORT_COL_STRING ? 1 : 0);
            if (written > 0) vartypes_pos += (size_t)written;
        }
        if (numtypes_pos < sizeof(numtypes) - 1) {
            written = snprintf(numtypes + numtypes_pos, sizeof(numtypes) - numtypes_pos,
                               "%s%d", sep, col->num_subtype);
            if (written > 0) numtypes_pos += (size_t)written;
        }
        if (strlens_pos < sizeof(strlens) - 1) {
            written = snprintf(strlens + strlens_pos, sizeof(strlens) - strlens_pos,
                               "%s%d", sep, col->max_strlen);
            if (written > 0) strlens_pos += (size_t)written;
        }
    }

    SF_macro_save("_cimport_varnames", varnames);
    SF_macro_save("_cimport_vartypes", vartypes);
    SF_macro_save("_cimport_numtypes", numtypes);
    SF_macro_save("_cimport_strlens", strlens);
}

ST_retcode xlsx_import_main(const char *args)
{
    ST_retcode rc = 0;

    /* Parse arguments */
    if (!args || !*args) {
        SF_error("xlsx_import_main: no arguments\n");
        return 198;
    }

    char args_copy[4096];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    /* Parse mode (scan or load) */
    char *mode = strtok(args_copy, " \t");
    if (!mode) {
        SF_error("xlsx_import_main: missing mode\n");
        return 198;
    }

    bool is_scan = (strcmp(mode, "scan") == 0);
    bool is_load = (strcmp(mode, "load") == 0);
    if (!is_scan && !is_load) {
        char err_buf[256];
        snprintf(err_buf, sizeof(err_buf), "xlsx_import_main: invalid mode '%s'\n", mode);
        SF_error(err_buf);
        return 198;
    }

    /* Parse filename */
    char *filename_tok = strtok(NULL, " \t");
    if (!filename_tok) {
        SF_error("xlsx_import_main: missing filename\n");
        return 198;
    }

    char filename_buf[1024];
    strncpy(filename_buf, filename_tok, sizeof(filename_buf) - 1);
    filename_buf[sizeof(filename_buf) - 1] = '\0';
    char *filename = filename_buf;

    /* Extract sheet name */
    char sheet_name_buf[XLSX_MAX_SHEET_NAME];
    char *sheet_name = xlsx_extract_sheet_name(args, sheet_name_buf, sizeof(sheet_name_buf));

    if (is_scan) {
        /* Clear any previous cached context */
        xlsx_clear_cached_context();

        /* Allocate new context on heap for caching */
        g_xlsx_ctx = (XLSXContext *)calloc(1, sizeof(XLSXContext));
        if (!g_xlsx_ctx) {
            SF_error("xlsx_import_main: memory allocation failed\n");
            return 920;
        }
        xlsx_context_init(g_xlsx_ctx);
        xlsx_parse_options(g_xlsx_ctx, args);

        /* Full parse pipeline */
        rc = xlsx_full_parse(g_xlsx_ctx, filename, sheet_name);
        if (rc) {
            xlsx_clear_cached_context();
            return rc;
        }

        size_t nobs = xlsx_compute_nobs(g_xlsx_ctx);
        xlsx_save_scan_macros(g_xlsx_ctx, nobs);

        if (g_xlsx_ctx->verbose) {
            char msg[256];
            snprintf(msg, sizeof(msg), "XLSX scan: %d sheets, selected '%s'\n",
                     g_xlsx_ctx->num_sheets, g_xlsx_ctx->sheets[g_xlsx_ctx->selected_sheet].name);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Shared strings: %u\n", g_xlsx_ctx->num_shared_strings);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Rows: %d-%d, Columns: %d\n",
                     g_xlsx_ctx->min_row, g_xlsx_ctx->max_row, g_xlsx_ctx->num_columns);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Data rows (after range/firstrow): %zu\n", nobs);
            SF_display(msg);
        }

        /* Context remains cached for the subsequent load call */
    }
    else if (is_load) {
        XLSXContext *ctx;
        bool used_cache = false;

        /* Try to reuse cached context from scan */
        if (g_xlsx_ctx && g_xlsx_ctx->filename &&
            strcmp(g_xlsx_ctx->filename, filename) == 0) {
            ctx = g_xlsx_ctx;
            used_cache = true;
        } else {
            /* Cache miss — parse fresh */
            xlsx_clear_cached_context();

            g_xlsx_ctx = (XLSXContext *)calloc(1, sizeof(XLSXContext));
            if (!g_xlsx_ctx) {
                SF_error("xlsx_import_main: memory allocation failed\n");
                return 920;
            }
            xlsx_context_init(g_xlsx_ctx);
            xlsx_parse_options(g_xlsx_ctx, args);

            rc = xlsx_full_parse(g_xlsx_ctx, filename, sheet_name);
            if (rc) {
                xlsx_clear_cached_context();
                return rc;
            }
            ctx = g_xlsx_ctx;
        }

        /* Build cache */
        rc = xlsx_build_cache(ctx);
        if (rc) {
            SF_error(ctx->error_message);
            xlsx_clear_cached_context();
            return rc;
        }

        /* Store data to Stata (parallel by column) */
        ST_int nvar = SF_nvar();
        ST_int stata_nobs = SF_nobs();
        int store_cols = ctx->num_columns < nvar ? ctx->num_columns : nvar;

        ctools_persistent_pool *pool = ctools_get_global_pool();
        if (pool && store_cols > 1) {
            xlsx_store_task *tasks = (xlsx_store_task *)malloc(store_cols * sizeof(xlsx_store_task));
            if (tasks) {
                for (int c = 0; c < store_cols; c++) {
                    tasks[c].col = &ctx->columns[c];
                    tasks[c].cache = &ctx->col_cache[c];
                    tasks[c].var = c + 1;
                    tasks[c].stata_nobs = stata_nobs;
                }
                ctools_persistent_pool_submit_batch(pool, xlsx_store_worker,
                                                     tasks, store_cols, sizeof(xlsx_store_task));
                ctools_persistent_pool_wait(pool);
                free(tasks);
            } else {
                goto xlsx_store_sequential;
            }
        } else {
xlsx_store_sequential:;
            for (int c = 0; c < store_cols; c++) {
                CImportColumnInfo *col = &ctx->columns[c];
                CImportColumnCache *cache = &ctx->col_cache[c];
                ST_int var = c + 1;

                if (col->type == CIMPORT_COL_STRING) {
                    for (size_t r = 0; r < cache->count && (ST_int)(r + 1) <= stata_nobs; r++) {
                        char *str = cache->string_data[r] ? cache->string_data[r] : (char *)"";
                        SF_sstore(var, (ST_int)(r + 1), str);
                    }
                } else {
                    for (size_t r = 0; r < cache->count && (ST_int)(r + 1) <= stata_nobs; r++) {
                        SF_vstore(var, (ST_int)(r + 1), cache->numeric_data[r]);
                    }
                }
            }
        }

        if (ctx->verbose) {
            char msg[256];
            SF_display("XLSX load complete");
            if (used_cache) SF_display(" (cached from scan)");
            SF_display("\n");
            snprintf(msg, sizeof(msg), "  Time ZIP: %.3f s\n", ctx->time_zip);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Time shared strings: %.3f s\n", ctx->time_shared_strings);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Time worksheet: %.3f s\n", ctx->time_worksheet);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Time type inference: %.3f s\n", ctx->time_type_infer);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Time cache: %.3f s\n", ctx->time_cache);
            SF_display(msg);
        }

        /* Free cached context after load is complete */
        xlsx_clear_cached_context();
    }

    return 0;
}

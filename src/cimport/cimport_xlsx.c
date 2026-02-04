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
#include "../ctools_timer.h"
#include "../ctools_arena.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

/* ============================================================================
 * Forward Declarations
 * ============================================================================ */

static bool workbook_callback(const xlsx_xml_event *event, void *user_data);
static bool shared_strings_callback(const xlsx_xml_event *event, void *user_data);
static bool worksheet_callback(const xlsx_xml_event *event, void *user_data);
static bool styles_callback(const xlsx_xml_event *event, void *user_data);
static void generate_varname(char *buf, int col, const char *header_value, int case_mode);
static double excel_date_to_stata(double excel_date);

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
    if (ctx->shared_strings_pool) {
        free(ctx->shared_strings_pool);
        ctx->shared_strings_pool = NULL;
    }

    /* Free date styles */
    if (ctx->date_styles) {
        free(ctx->date_styles);
        ctx->date_styles = NULL;
    }

    /* Free parsed rows */
    if (ctx->rows) {
        for (int i = 0; i < ctx->num_rows; i++) {
            if (ctx->rows[i].cells) {
                /* Free any inline strings within cells */
                for (int j = 0; j < ctx->rows[i].num_cells; j++) {
                    if (ctx->rows[i].cells[j].type == XLSX_CELL_STRING &&
                        ctx->rows[i].cells[j].value.inline_string != NULL) {
                        free((void *)ctx->rows[i].cells[j].value.inline_string);
                    }
                }
                free(ctx->rows[i].cells);
            }
        }
        free(ctx->rows);
        ctx->rows = NULL;
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
        if (strcmp(event->tag_name, "si") == 0) {
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
                if (new_cap > XLSX_MAX_SHARED_STRINGS) {
                    new_cap = XLSX_MAX_SHARED_STRINGS;
                }
                char **new_strings = (char **)realloc(ctx->shared_strings,
                                                       new_cap * sizeof(char *));
                if (!new_strings) return false;
                ctx->shared_strings = new_strings;
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
            ctx->shared_strings[ctx->num_shared_strings++] = dest;
            ctx->shared_strings_pool_used += need;
        }
    }

    return true;
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

    SharedStringsParseState state = {0};
    state.ctx = ctx;

    xlsx_xml_parser *parser = xlsx_xml_parser_create(shared_strings_callback, &state);
    if (!parser) {
        free(data);
        return 920;
    }

    bool success = xlsx_xml_parse(parser, (const char *)data, size, true);

    xlsx_xml_parser_destroy(parser);
    free(data);

    ctx->time_shared_strings = ctools_timer_seconds() - t_start;

    if (!success) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to parse sharedStrings.xml");
        return 610;
    }

    return 0;
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

    StylesParseState state = {0};
    state.ctx = ctx;

    xlsx_xml_parser *parser = xlsx_xml_parser_create(styles_callback, &state);
    if (!parser) {
        free(data);
        return 0;
    }

    xlsx_xml_parse(parser, (const char *)data, size, true);

    xlsx_xml_parser_destroy(parser);
    free(data);

    return 0;
}

/* ============================================================================
 * Worksheet Parsing
 * ============================================================================ */

static XLSXParsedRow *ensure_row(XLSXContext *ctx, int row_num)
{
    /* Find or create row */
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

    if (row_num < ctx->min_row) ctx->min_row = row_num;
    if (row_num > ctx->max_row) ctx->max_row = row_num;

    return row;
}

static XLSXCell *add_cell(XLSXParsedRow *row, int col)
{
    if (row->num_cells >= row->capacity) {
        int new_cap = row->capacity == 0 ? 16 : row->capacity * 2;
        XLSXCell *new_cells = (XLSXCell *)realloc(row->cells,
                                                   new_cap * sizeof(XLSXCell));
        if (!new_cells) return NULL;
        row->cells = new_cells;
        row->capacity = new_cap;
    }

    XLSXCell *cell = &row->cells[row->num_cells++];
    memset(cell, 0, sizeof(XLSXCell));
    cell->col = col;
    return cell;
}

static bool worksheet_callback(const xlsx_xml_event *event, void *user_data)
{
    WorksheetParseState *state = (WorksheetParseState *)user_data;
    XLSXContext *ctx = state->ctx;

    if (event->type == XLSX_XML_START_ELEMENT) {
        if (strcmp(event->tag_name, "sheetData") == 0) {
            state->in_sheetdata = true;
        }
        else if (state->in_sheetdata && strcmp(event->tag_name, "row") == 0) {
            state->in_row = true;
            const char *r = xlsx_xml_get_attr(event, "r");
            state->current_row = r ? atoi(r) : (state->current_row + 1);

            /* Check cell range */
            XLSXCellRange *range = &ctx->cell_range;
            if (range->start_row > 0 && state->current_row < range->start_row) {
                state->in_row = false;  /* Skip this row */
            }
            if (range->end_row > 0 && state->current_row > range->end_row) {
                return false;  /* Stop parsing */
            }

            if (state->in_row) {
                state->row_ptr = ensure_row(ctx, state->current_row);
            }
        }
        else if (state->in_row && strcmp(event->tag_name, "c") == 0) {
            state->in_cell = true;
            state->value_len = 0;

            const char *r = xlsx_xml_get_attr(event, "r");
            const char *t = xlsx_xml_get_attr(event, "t");
            const char *s = xlsx_xml_get_attr(event, "s");

            /* Parse cell reference */
            if (r) {
                int col, row;
                if (xlsx_xml_parse_cell_ref(r, &col, &row)) {
                    state->current_col = col;
                    state->current_row = row;
                }
            }

            /* Check column range */
            XLSXCellRange *range = &ctx->cell_range;
            if (range->start_col >= 0 && state->current_col < range->start_col) {
                state->in_cell = false;  /* Skip */
            }
            if (range->end_col >= 0 && state->current_col > range->end_col) {
                state->in_cell = false;  /* Skip */
            }

            /* Cell type: s=shared string, str=inline, b=boolean, e=error, n=number (default) */
            state->current_type = XLSX_CELL_NUMBER;  /* Default */
            if (t) {
                if (strcmp(t, "s") == 0) {
                    state->current_type = XLSX_CELL_SHARED_STRING;
                } else if (strcmp(t, "str") == 0 || strcmp(t, "inlineStr") == 0) {
                    state->current_type = XLSX_CELL_STRING;
                } else if (strcmp(t, "b") == 0) {
                    state->current_type = XLSX_CELL_BOOLEAN;
                } else if (strcmp(t, "e") == 0) {
                    state->current_type = XLSX_CELL_ERROR;
                }
            }

            /* Style index */
            state->current_style = s ? atoi(s) : -1;

            if (state->current_col > ctx->max_col) {
                ctx->max_col = state->current_col;
            }
        }
        else if (state->in_cell && strcmp(event->tag_name, "v") == 0) {
            state->in_value = true;
            state->value_len = 0;
        }
        else if (state->in_cell && strcmp(event->tag_name, "is") == 0) {
            state->in_inline_str = true;
        }
        else if (state->in_inline_str && strcmp(event->tag_name, "t") == 0) {
            state->in_t = true;
            state->value_len = 0;
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
        if (strcmp(event->tag_name, "sheetData") == 0) {
            state->in_sheetdata = false;
        }
        else if (strcmp(event->tag_name, "row") == 0) {
            state->in_row = false;
            state->row_ptr = NULL;
        }
        else if (strcmp(event->tag_name, "v") == 0) {
            state->in_value = false;
        }
        else if (strcmp(event->tag_name, "t") == 0) {
            state->in_t = false;
        }
        else if (strcmp(event->tag_name, "is") == 0) {
            state->in_inline_str = false;
        }
        else if (strcmp(event->tag_name, "c") == 0 && state->in_cell) {
            state->in_cell = false;
            state->value_buf[state->value_len] = '\0';

            if (state->row_ptr && state->value_len > 0) {
                XLSXCell *cell = add_cell(state->row_ptr, state->current_col);
                if (cell) {
                    cell->col = state->current_col;
                    cell->row = state->current_row;
                    cell->type = state->current_type;
                    cell->style_idx = state->current_style;

                    switch (state->current_type) {
                    case XLSX_CELL_SHARED_STRING:
                        cell->value.shared_string_idx = atoi(state->value_buf);
                        break;
                    case XLSX_CELL_NUMBER:
                        cell->value.number = strtod(state->value_buf, NULL);
                        /* Check if this is a date based on style */
                        if (state->current_style >= 0 &&
                            state->current_style < ctx->num_styles &&
                            ctx->date_styles[state->current_style]) {
                            cell->type = XLSX_CELL_DATE;
                        }
                        break;
                    case XLSX_CELL_BOOLEAN:
                        cell->value.boolean = (state->value_buf[0] == '1');
                        break;
                    case XLSX_CELL_STRING:
                        /* Inline string - decode entities */
                        xlsx_xml_decode_entities(state->value_buf, state->value_len);
                        /* Store as shared string for uniformity */
                        cell->type = XLSX_CELL_STRING;
                        cell->value.inline_string = strdup(state->value_buf);
                        if (cell->value.inline_string == NULL) {
                            return false;  /* Memory allocation failed */
                        }
                        break;
                    default:
                        break;
                    }
                }
            }
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

    size_t size;
    void *data = xlsx_zip_extract_file(zip, worksheet_path, &size);
    if (!data) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to extract worksheet: %s", worksheet_path);
        return 610;
    }

    WorksheetParseState state = {0};
    state.ctx = ctx;

    xlsx_xml_parser *parser = xlsx_xml_parser_create(worksheet_callback, &state);
    if (!parser) {
        free(data);
        return 920;
    }

    xlsx_xml_parse(parser, (const char *)data, size, true);

    xlsx_xml_parser_destroy(parser);
    free(data);

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

    /* Initialize columns */
    for (int c = 0; c < ctx->num_columns; c++) {
        CImportColumnInfo *col = &ctx->columns[c];
        col->type = CIMPORT_COL_UNKNOWN;
        col->num_subtype = CIMPORT_NUM_BYTE;
        col->is_integer = true;
        col->min_value = DBL_MAX;
        col->max_value = -DBL_MAX;
        col->max_strlen = 0;
    }

    /* Determine header row */
    int data_start_row = ctx->min_row;
    int header_row = -1;

    if (ctx->firstrow && ctx->num_rows > 0) {
        header_row = ctx->min_row;
        data_start_row = ctx->min_row + 1;
    }

    /* Extract column names from header row */
    if (header_row > 0) {
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

    /* Generate default names for columns without headers */
    for (int c = 0; c < ctx->num_columns; c++) {
        if (ctx->columns[c].name[0] == '\0') {
            generate_varname(ctx->columns[c].name, c, NULL, ctx->case_mode);
        }
    }

    /* Scan data rows for type inference */
    for (int i = 0; i < ctx->num_rows; i++) {
        XLSXParsedRow *row = &ctx->rows[i];
        if (row->row_num < data_start_row) continue;

        for (int j = 0; j < row->num_cells; j++) {
            XLSXCell *cell = &row->cells[j];
            int col_idx = cell->col - start_col;
            if (col_idx < 0 || col_idx >= ctx->num_columns) continue;

            CImportColumnInfo *col = &ctx->columns[col_idx];

            if (ctx->allstring) {
                col->type = CIMPORT_COL_STRING;
            } else {
                switch (cell->type) {
                case XLSX_CELL_SHARED_STRING:
                case XLSX_CELL_STRING:
                    col->type = CIMPORT_COL_STRING;
                    {
                        const char *str = NULL;
                        if (cell->type == XLSX_CELL_SHARED_STRING &&
                            (uint32_t)cell->value.shared_string_idx < ctx->num_shared_strings) {
                            str = ctx->shared_strings[cell->value.shared_string_idx];
                        } else if (cell->type == XLSX_CELL_STRING) {
                            str = cell->value.inline_string;
                        }
                        if (str) {
                            int len = (int)strlen(str);
                            if (len > col->max_strlen) col->max_strlen = len;
                        }
                    }
                    break;

                case XLSX_CELL_NUMBER:
                case XLSX_CELL_DATE:
                    if (col->type == CIMPORT_COL_UNKNOWN) {
                        col->type = CIMPORT_COL_NUMERIC;
                    }
                    if (col->type == CIMPORT_COL_NUMERIC) {
                        double val = cell->value.number;
                        if (cell->type == XLSX_CELL_DATE) {
                            val = excel_date_to_stata(val);
                        }
                        if (val < col->min_value) col->min_value = val;
                        if (val > col->max_value) col->max_value = val;
                        if (col->is_integer && val != floor(val)) {
                            col->is_integer = false;
                        }
                    }
                    break;

                case XLSX_CELL_BOOLEAN:
                    if (col->type == CIMPORT_COL_UNKNOWN) {
                        col->type = CIMPORT_COL_NUMERIC;
                    }
                    break;

                default:
                    break;
                }
            }
        }
    }

    /* Finalize numeric subtypes */
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

ST_retcode xlsx_build_cache(XLSXContext *ctx)
{
    double t_start = ctools_timer_seconds();

    /* Determine data rows */
    int start_col = ctx->cell_range.start_col >= 0 ? ctx->cell_range.start_col : 0;
    int data_start_row = ctx->min_row;
    if (ctx->firstrow) data_start_row = ctx->min_row + 1;

    int range_start_row = ctx->cell_range.start_row > 0 ? ctx->cell_range.start_row : data_start_row;
    int range_end_row = ctx->cell_range.end_row > 0 ? ctx->cell_range.end_row : ctx->max_row;

    /* Count data rows */
    size_t total_data_rows = 0;
    for (int i = 0; i < ctx->num_rows; i++) {
        if (ctx->rows[i].row_num >= data_start_row &&
            ctx->rows[i].row_num >= range_start_row &&
            ctx->rows[i].row_num <= range_end_row) {
            total_data_rows++;
        }
    }

    /* Also count empty rows within range */
    if (range_end_row >= range_start_row) {
        total_data_rows = range_end_row - range_start_row + 1;
        if (ctx->firstrow && range_start_row == ctx->min_row) {
            total_data_rows--;
        }
    }

    if (total_data_rows == 0) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "No data rows found");
        return 2000;
    }

    /* Allocate caches */
    ctx->col_cache = (CImportColumnCache *)calloc(ctx->num_columns,
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
            /* Initialize with missing values */
            for (size_t i = 0; i < total_data_rows; i++) {
                cache->numeric_data[i] = SV_missval;
            }
        }
    }

    /* Build row index map */
    int actual_data_start = data_start_row;
    if (actual_data_start < range_start_row) actual_data_start = range_start_row;

    /* Fill caches */
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
                    /* Convert number to string */
                    static char num_buf[64];
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
 * Main Entry Point
 * ============================================================================ */

ST_retcode xlsx_import_main(const char *args)
{
    ST_retcode rc = 0;
    XLSXContext ctx;
    xlsx_context_init(&ctx);

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

    /* Save filename to separate buffer before we re-parse args_copy */
    char filename_buf[1024];
    strncpy(filename_buf, filename_tok, sizeof(filename_buf) - 1);
    filename_buf[sizeof(filename_buf) - 1] = '\0';
    char *filename = filename_buf;

    /* Parse options */
    char *opt;
    while ((opt = strtok(NULL, " \t")) != NULL) {
        if (strncmp(opt, "sheet=", 6) == 0) {
            /* Will be used after parsing workbook */
            /* Store for later */
        } else if (strncmp(opt, "cellrange=", 10) == 0) {
            xlsx_parse_cellrange(opt + 10, &ctx.cell_range);
        } else if (strcmp(opt, "firstrow") == 0) {
            ctx.firstrow = true;
        } else if (strcmp(opt, "allstring") == 0) {
            ctx.allstring = true;
        } else if (strncmp(opt, "case=", 5) == 0) {
            if (strcmp(opt + 5, "lower") == 0) ctx.case_mode = 1;
            else if (strcmp(opt + 5, "upper") == 0) ctx.case_mode = 2;
        } else if (strcmp(opt, "verbose") == 0) {
            ctx.verbose = true;
        }
    }

    /* Re-parse for sheet option since we need the context */
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    char *sheet_name = NULL;
    char *p = strstr(args_copy, "sheet=");
    if (p) {
        sheet_name = p + 6;
        char *end = sheet_name;
        while (*end && !isspace((unsigned char)*end)) end++;
        *end = '\0';
    }

    /* Open file */
    rc = xlsx_open_file(&ctx, filename);
    if (rc) {
        SF_error(ctx.error_message);
        xlsx_context_free(&ctx);
        return rc;
    }

    /* Parse workbook to get sheet list */
    rc = xlsx_parse_workbook(&ctx);
    if (rc) {
        SF_error(ctx.error_message);
        xlsx_context_free(&ctx);
        return rc;
    }

    if (ctx.num_sheets == 0) {
        SF_error("No sheets found in workbook\n");
        xlsx_context_free(&ctx);
        return 610;
    }

    /* Select sheet */
    if (sheet_name) {
        int idx = xlsx_select_sheet_by_name(&ctx, sheet_name);
        if (idx < 0) {
            char err_buf[256];
            snprintf(err_buf, sizeof(err_buf), "Sheet not found: %s\n", sheet_name);
            SF_error(err_buf);
            xlsx_context_free(&ctx);
            return 601;
        }
    }

    /* Parse shared strings */
    rc = xlsx_parse_shared_strings(&ctx);
    if (rc) {
        SF_error(ctx.error_message);
        xlsx_context_free(&ctx);
        return rc;
    }

    /* Parse styles for date detection */
    xlsx_parse_styles(&ctx);

    /* Parse worksheet */
    rc = xlsx_parse_worksheet(&ctx);
    if (rc) {
        SF_error(ctx.error_message);
        xlsx_context_free(&ctx);
        return rc;
    }

    /* Infer types */
    rc = xlsx_infer_types(&ctx);
    if (rc) {
        SF_error(ctx.error_message);
        xlsx_context_free(&ctx);
        return rc;
    }

    /* Calculate actual row count */
    int data_start = ctx.min_row;
    if (ctx.firstrow) data_start++;
    int range_start = ctx.cell_range.start_row > 0 ? ctx.cell_range.start_row : data_start;
    int range_end = ctx.cell_range.end_row > 0 ? ctx.cell_range.end_row : ctx.max_row;
    if (range_start < data_start) range_start = data_start;

    size_t nobs = 0;
    if (range_end >= range_start) {
        nobs = range_end - range_start + 1;
    }

    if (is_scan) {
        /* Save metadata to Stata macros */
        char buf[32];

        snprintf(buf, sizeof(buf), "%zu", nobs);
        SF_macro_save("_cimport_nobs", buf);

        snprintf(buf, sizeof(buf), "%d", ctx.num_columns);
        SF_macro_save("_cimport_nvar", buf);

        /* Build varnames string with bounds checking */
        char varnames[8192] = "";
        char vartypes[1024] = "";
        char numtypes[1024] = "";
        char strlens[1024] = "";
        size_t varnames_pos = 0, vartypes_pos = 0, numtypes_pos = 0, strlens_pos = 0;
        int written;

        for (int c = 0; c < ctx.num_columns; c++) {
            CImportColumnInfo *col = &ctx.columns[c];
            const char *sep = (c > 0) ? " " : "";

            /* Append to varnames with bounds check */
            if (varnames_pos < sizeof(varnames) - 1) {
                written = snprintf(varnames + varnames_pos, sizeof(varnames) - varnames_pos,
                                   "%s%s", sep, col->name);
                if (written > 0) varnames_pos += (size_t)written;
            }

            /* Append to vartypes */
            if (vartypes_pos < sizeof(vartypes) - 1) {
                written = snprintf(vartypes + vartypes_pos, sizeof(vartypes) - vartypes_pos,
                                   "%s%d", sep, col->type == CIMPORT_COL_STRING ? 1 : 0);
                if (written > 0) vartypes_pos += (size_t)written;
            }

            /* Append to numtypes */
            if (numtypes_pos < sizeof(numtypes) - 1) {
                written = snprintf(numtypes + numtypes_pos, sizeof(numtypes) - numtypes_pos,
                                   "%s%d", sep, col->num_subtype);
                if (written > 0) numtypes_pos += (size_t)written;
            }

            /* Append to strlens */
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

        if (ctx.verbose) {
            char msg[256];
            snprintf(msg, sizeof(msg), "XLSX scan: %d sheets, selected '%s'\n",
                     ctx.num_sheets, ctx.sheets[ctx.selected_sheet].name);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Shared strings: %u\n", ctx.num_shared_strings);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Rows: %d-%d, Columns: %d\n",
                     ctx.min_row, ctx.max_row, ctx.num_columns);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Data rows (after range/firstrow): %zu\n", nobs);
            SF_display(msg);
        }
    }
    else if (is_load) {
        /* Build cache */
        rc = xlsx_build_cache(&ctx);
        if (rc) {
            SF_error(ctx.error_message);
            xlsx_context_free(&ctx);
            return rc;
        }

        /* Store data to Stata */
        ST_int nvar = SF_nvar();
        ST_int stata_nobs = SF_nobs();

        for (int c = 0; c < ctx.num_columns && c < nvar; c++) {
            CImportColumnInfo *col = &ctx.columns[c];
            CImportColumnCache *cache = &ctx.col_cache[c];
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

        if (ctx.verbose) {
            char msg[256];
            SF_display("XLSX load complete\n");
            snprintf(msg, sizeof(msg), "  Time ZIP: %.3f s\n", ctx.time_zip);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Time shared strings: %.3f s\n", ctx.time_shared_strings);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Time worksheet: %.3f s\n", ctx.time_worksheet);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Time type inference: %.3f s\n", ctx.time_type_infer);
            SF_display(msg);
            snprintf(msg, sizeof(msg), "  Time cache: %.3f s\n", ctx.time_cache);
            SF_display(msg);
        }
    }

    xlsx_context_free(&ctx);
    return 0;
}

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
#include "../io/cio_excel_format.h"

/* ============================================================================
 * Forward Declarations
 * ============================================================================ */

static bool workbook_callback(const xlsx_xml_event *event, void *user_data);
static bool shared_strings_callback(const xlsx_xml_event *event, void *user_data);
static bool styles_callback(const xlsx_xml_event *event, void *user_data);

/* Tagged column-major slots store either a number or an arena-owned pointer.
 * memcpy avoids aliasing violations; supported plugin platforms use 64-bit pointers. */
_Static_assert(sizeof(const char *) <= sizeof(double), "XLSX slot must hold a pointer");

/* Parser state structures */
typedef struct {
    XLSXContext *ctx;
    bool in_sheets;
} WorkbookParseState;

/* Text of the current <si>: every <t> run except phonetic guides (<rPh>),
 * accumulated without a length limit and entity-decoded once complete. */
typedef struct {
    XLSXContext *ctx;
    bool in_si;
    bool in_t;
    int phonetic_depth;
    char *text;
    size_t text_len;
    size_t text_capacity;
    bool failed;
} SharedStringsParseState;

typedef struct {
    XLSXContext *ctx;
    bool in_numfmts;
    bool in_cellxfs;
    int current_numfmt_id;
    int current_xf_index;
    bool failed;
    /* Track which numFmtIds are dates */
    uint8_t native_numfmt[65536];
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
    snprintf(ctx->workbook_path, sizeof(ctx->workbook_path), "%s", XLSX_PATH_WORKBOOK);
    snprintf(ctx->shared_strings_path, sizeof(ctx->shared_strings_path), "%s",
             XLSX_PATH_SHARED_STRINGS);
    snprintf(ctx->styles_path, sizeof(ctx->styles_path), "%s", XLSX_PATH_STYLES);
    ctools_arena_init(&ctx->parse_arena, 0);
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

    free(ctx->excel_formats);ctx->excel_formats=NULL;
    /* Free date styles */
    if (ctx->date_styles) {
        free(ctx->date_styles);
        ctx->date_styles = NULL;
    }

    /* Free column-major arrays (inline strings freed via parse_arena below) */
    for (int c = 0; c < ctx->cm_num_cols; c++) {
        free(ctx->cm_numeric[c]);
        free(ctx->cm_types[c]);
        free(ctx->cm_formats[c]);
    }
    free(ctx->cm_numeric);
    free(ctx->cm_types);
    free(ctx->cm_formats);
    ctx->cm_numeric = NULL;
    ctx->cm_types = NULL;
    ctx->cm_formats = NULL;
    ctx->cm_num_cols = 0;
    ctx->cm_num_rows = 0;

    /* Free arenas */
    ctools_arena_free(&ctx->parse_arena);

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
 * Package Relationships
 *
 * Parts are located through relationship parts, not fixed names: _rels/.rels
 * names the workbook, and the workbook's .rels names each sheet (tab order is
 * not tied to sheet{N}.xml, and chart sheets live elsewhere), sharedStrings
 * and styles.
 * ============================================================================ */

typedef struct {
    XLSXContext *ctx;
    const char *dir;          /* directory of the source part ("" or ".../") */
    char *office_document;    /* non-NULL for _rels/.rels: the workbook part */
} RelsParseState;

/* Transitional and strict OOXML type URIs share their final segment. */
static bool xlsx_rel_type(const char *type, const char *name)
{
    size_t n = strlen(type), m = strlen(name);
    return n > m && type[n - m - 1] == '/' && !strcmp(type + n - m, name);
}

/* Resolve a relationship target against the directory of its source part;
 * a leading '/' starts at the package root. Handles "." and ".." segments. */
static bool xlsx_resolve_part(const char *dir, const char *target, char *out, size_t cap)
{
    char joined[2 * XLSX_MAX_PART_PATH];
    if (target[0] == '/') snprintf(joined, sizeof(joined), "%s", target + 1);
    else snprintf(joined, sizeof(joined), "%s%s", dir, target);
    size_t n = 0;
    for (const char *p = joined; *p;) {
        const char *slash = strchr(p, '/');
        size_t len = slash ? (size_t)(slash - p) : strlen(p);
        if (len == 2 && p[0] == '.' && p[1] == '.') {
            if (n) n--;
            while (n && out[n - 1] != '/') n--;
        } else if (len && !(len == 1 && p[0] == '.')) {
            if (n + len + 1 >= cap) return false;
            memcpy(out + n, p, len);
            n += len;
            out[n++] = '/';
        }
        p += len + (slash != NULL);
    }
    if (!n) return false;
    out[n - 1] = '\0';
    return true;
}

static bool rels_callback(const xlsx_xml_event *event, void *user_data)
{
    RelsParseState *state = (RelsParseState *)user_data;
    XLSXContext *ctx = state->ctx;
    if (event->type != XLSX_XML_START_ELEMENT || strcmp(event->tag_name, "Relationship"))
        return true;
    const char *id_attr = xlsx_xml_get_attr(event, "Id");
    const char *type = xlsx_xml_get_attr(event, "Type");
    const char *target_attr = xlsx_xml_get_attr(event, "Target");
    const char *mode = xlsx_xml_get_attr(event, "TargetMode");
    char path[XLSX_MAX_PART_PATH];
    if (!id_attr || !type || !target_attr || (mode && !strcmp(mode, "External")))
        return true;
    /* Attribute values are raw XML text (entities such as &amp;) */
    char id[256], target[XLSX_MAX_PART_PATH];
    snprintf(id, sizeof(id), "%s", id_attr);
    xlsx_xml_decode_entities(id, strlen(id));
    snprintf(target, sizeof(target), "%s", target_attr);
    xlsx_xml_decode_entities(target, strlen(target));
    if (!xlsx_resolve_part(state->dir, target, path, sizeof(path)))
        return true;
    if (state->office_document) {
        if (xlsx_rel_type(type, "officeDocument"))
            snprintf(state->office_document, XLSX_MAX_PART_PATH, "%s", path);
        return true;
    }
    if (xlsx_rel_type(type, "sharedStrings"))
        snprintf(ctx->shared_strings_path, sizeof(ctx->shared_strings_path), "%s", path);
    else if (xlsx_rel_type(type, "styles"))
        snprintf(ctx->styles_path, sizeof(ctx->styles_path), "%s", path);
    for (int i = 0; i < ctx->num_sheets; i++) {
        XLSXSheetInfo *info = &ctx->sheets[i];
        if (strcmp(info->rel_id, id)) continue;
        snprintf(info->path, sizeof(info->path), "%s", path);
        info->kind = xlsx_rel_type(type, "worksheet")     ? XLSX_SHEET_WORKSHEET
                   : xlsx_rel_type(type, "chartsheet")    ? XLSX_SHEET_CHART
                   : xlsx_rel_type(type, "dialogsheet")   ? XLSX_SHEET_DIALOG
                   : xlsx_rel_type(type, "xlMacrosheet") ||
                     xlsx_rel_type(type, "xlIntlMacrosheet") ? XLSX_SHEET_MACRO
                                                             : XLSX_SHEET_OTHER;
    }
    return true;
}

/* Apply the relationships of part (""  for the package). A missing .rels
 * part leaves the defaults in place. */
static ST_retcode xlsx_parse_rels(XLSXContext *ctx, const char *part, char *office_document)
{
    char dir[XLSX_MAX_PART_PATH], rels[2 * XLSX_MAX_PART_PATH];
    const char *slash = strrchr(part, '/');
    int dlen = slash ? (int)(slash - part + 1) : 0;
    snprintf(dir, sizeof(dir), "%.*s", dlen, part);
    snprintf(rels, sizeof(rels), "%s_rels/%s.rels", dir, part + dlen);
    size_t size;
    void *data = xlsx_zip_extract_file((xlsx_zip_archive *)ctx->zip_archive, rels, &size);
    if (!data) return 0;
    RelsParseState state = {ctx, dir, office_document};
    xlsx_xml_parser *parser = xlsx_xml_parser_create(rels_callback, &state);
    if (!parser) {
        free(data);
        return 920;
    }
    xlsx_xml_parse(parser, (const char *)data, size, true);
    xlsx_xml_parser_destroy(parser);
    free(data);
    return 0;
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

    /* Verify it's a valid XLSX: the package relationships name the workbook */
    xlsx_zip_archive *zip = (xlsx_zip_archive *)ctx->zip_archive;
    char workbook[XLSX_MAX_PART_PATH] = "";
    ST_retcode rc = xlsx_parse_rels(ctx, "", workbook);
    if (rc) return rc;
    if (*workbook && xlsx_zip_locate_file(zip, workbook) != (size_t)-1)
        snprintf(ctx->workbook_path, sizeof(ctx->workbook_path), "%s", workbook);
    if (xlsx_zip_locate_file(zip, ctx->workbook_path) == (size_t)-1) {
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
        if (strcmp(event->tag_name, "workbookPr") == 0) {
            const char *date1904 = xlsx_xml_get_attr(event, "date1904");
            ctx->date1904 = date1904 && (!strcmp(date1904, "1") || !strcmp(date1904, "true"));
        }
        else if (strcmp(event->tag_name, "sheets") == 0) {
            state->in_sheets = true;
        }
        else if (state->in_sheets && strcmp(event->tag_name, "sheet") == 0) {
            if (ctx->num_sheets >= XLSX_MAX_SHEETS) return true;

            XLSXSheetInfo *info = &ctx->sheets[ctx->num_sheets];

            const char *name = xlsx_xml_get_attr(event, "name");
            const char *sheet_id = xlsx_xml_get_attr(event, "sheetId");
            const char *rel_id = xlsx_xml_get_attr(event, "id");  /* r:id */

            /* Attribute values are raw XML text: "P&amp;L" names sheet P&L */
            if (name) {
                strncpy(info->name, name, XLSX_MAX_SHEET_NAME - 1);
                info->name[XLSX_MAX_SHEET_NAME - 1] = '\0';
                xlsx_xml_decode_entities(info->name, strlen(info->name));
            }
            if (sheet_id) {
                info->sheet_id = atoi(sheet_id);
            }
            if (rel_id) {
                strncpy(info->rel_id, rel_id, sizeof(info->rel_id) - 1);
                info->rel_id[sizeof(info->rel_id) - 1] = '\0';
                xlsx_xml_decode_entities(info->rel_id, strlen(info->rel_id));
            }

            /* Position names the part only if no relationship resolves it */
            info->sheet_index = ctx->num_sheets + 1;
            snprintf(info->path, sizeof(info->path), XLSX_WORKSHEET_FMT, info->sheet_index);
            info->kind = XLSX_SHEET_WORKSHEET;

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
    void *data = xlsx_zip_extract_file(zip, ctx->workbook_path, &size);
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

    ST_retcode rc = xlsx_parse_rels(ctx, ctx->workbook_path, NULL);
    if (rc) return rc;

    /* Without sheet(), import the first worksheet (tabs may be chart sheets) */
    ctx->selected_sheet = 0;
    for (int i = 0; i < ctx->num_sheets; i++) {
        if (ctx->sheets[i].kind == XLSX_SHEET_WORKSHEET) {
            ctx->selected_sheet = i;
            break;
        }
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
            state->in_t = false;
            state->phonetic_depth = 0;
            state->text_len = 0;
        }
        else if (state->in_si && strcmp(event->tag_name, "rPh") == 0) {
            state->phonetic_depth++;
        }
        else if (state->in_si && !state->phonetic_depth &&
                 strcmp(event->tag_name, "t") == 0) {
            state->in_t = true;
        }
    }
    else if (event->type == XLSX_XML_TEXT) {
        if (state->in_t) {
            /* Append text to current string, including whitespace-only runs */
            size_t need = state->text_len + event->text_len + 1;
            if (need > state->text_capacity) {
                size_t capacity = state->text_capacity ? state->text_capacity : 256;
                while (capacity < need) {
                    if (capacity > SIZE_MAX / 2) {
                        state->failed = true;
                        return false;
                    }
                    capacity *= 2;
                }
                char *text = (char *)realloc(state->text, capacity);
                if (!text) {
                    state->failed = true;
                    return false;
                }
                state->text = text;
                state->text_capacity = capacity;
            }
            memcpy(state->text + state->text_len, event->text, event->text_len);
            state->text_len += event->text_len;
        }
    }
    else if (event->type == XLSX_XML_END_ELEMENT) {
        if (strcmp(event->tag_name, "t") == 0) {
            state->in_t = false;
        }
        else if (strcmp(event->tag_name, "rPh") == 0) {
            if (state->phonetic_depth) state->phonetic_depth--;
        }
        else if (strcmp(event->tag_name, "si") == 0) {
            state->in_si = false;
            state->failed = true;  /* cleared once the string is stored */

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

            /* Decode XML entities once, over the whole string */
            const char *string = "";
            size_t decoded_len = 0;
            if (state->text_len) {
                state->text[state->text_len] = '\0';
                decoded_len = xlsx_xml_decode_entities(state->text, state->text_len);
                string = state->text;
            }

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
            memcpy(dest, string, decoded_len + 1);
            ctx->shared_strings[ctx->num_shared_strings] = dest;
            if (ctx->shared_string_lengths) {
                ctx->shared_string_lengths[ctx->num_shared_strings] = (uint32_t)decoded_len;
            }
            ctx->num_shared_strings++;
            ctx->shared_strings_pool_used += need;
            state->failed = false;
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
    free(state.text);

    /* A callback stops the parser only when memory runs out */
    if (state.failed) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Insufficient memory for sharedStrings.xml");
        return 909;
    }
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
    size_t idx = xlsx_zip_locate_file(zip, ctx->shared_strings_path);
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
                char decoded[4096];snprintf(decoded,sizeof(decoded),"%s",formatCode);
                xlsx_xml_decode_entities(decoded,strlen(decoded));
                if(id>=0 && id<65536)state->native_numfmt[id]=cio_excel_format_code(id,decoded);
            }
        }
        else if (strcmp(event->tag_name, "cellXfs") == 0) {
            state->in_cellxfs = true;
            state->current_xf_index = 0;
        }
        else if (state->in_cellxfs && strcmp(event->tag_name, "xf") == 0) {
            const char *numFmtId = xlsx_xml_get_attr(event, "numFmtId");
            uint8_t native_format=0;
            if(numFmtId) {
                int id=atoi(numFmtId);
                native_format=cio_excel_format_code(id,NULL);
                if(!native_format && id>=0 && id<65536)native_format=state->native_numfmt[id];
            }
            uint8_t is_date=cio_excel_date_kind(native_format);

            /* Expand date_styles array if needed */
            if (state->current_xf_index >= ctx->num_styles) {
                if (ctx->num_styles > INT_MAX / 2) { state->failed = true; return false; }
                int new_size = ctx->num_styles == 0 ? 64 : ctx->num_styles * 2;
                uint8_t *new_styles = (uint8_t *)realloc(ctx->date_styles,
                                                    new_size * sizeof(uint8_t));
                if (!new_styles) {
                    state->failed = true;
                    return false;  /* Memory allocation failed */
                }
                memset(new_styles + ctx->num_styles, 0,
                       (new_size - ctx->num_styles) * sizeof(uint8_t));
                ctx->date_styles = new_styles;
                uint8_t *new_formats=realloc(ctx->excel_formats,(size_t)new_size);
                if (!new_formats) { state->failed = true; return false; }
                memset(new_formats+ctx->num_styles,0,(size_t)(new_size-ctx->num_styles));
                ctx->excel_formats=new_formats;
                ctx->num_styles = new_size;
            }

            if (state->current_xf_index < ctx->num_styles) {
                ctx->date_styles[state->current_xf_index] = is_date;
                ctx->excel_formats[state->current_xf_index]=native_format;
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
    if (!parser) return 920;

    bool ok=xlsx_xml_parse(parser, (const char *)data, size, true);
    xlsx_xml_parser_destroy(parser);

    return state.failed ? 920 : (ok ? 0 : 610);
}

ST_retcode xlsx_parse_styles(XLSXContext *ctx)
{
    xlsx_zip_archive *zip = (xlsx_zip_archive *)ctx->zip_archive;

    size_t idx = xlsx_zip_locate_file(zip, ctx->styles_path);
    if (idx == (size_t)-1) {
        return 0;  /* Styles are optional */
    }

    size_t size;
    void *data = xlsx_zip_extract_to_heap(zip, idx, &size);
    if (!data) {
        return 610;  /* A present styles part must be readable. */
    }

    ST_retcode rc = xlsx_parse_styles_buf(ctx, data, size);
    free(data);

    return rc;
}

/* ============================================================================
 * Worksheet Scanner
 *
 * A memchr-driven pass over the in-memory worksheet part that processes each
 * <row> and <c> element as a unit. It accepts namespace prefixes (<x:c>),
 * either quote character, any XML whitespace, and rows or cells without r=
 * (positions then follow document order). <dimension> only seeds the first
 * allocation: storage grows to hold every cell that carries a value.
 * ============================================================================ */

#define XLSX_PARALLEL_MIN_BYTES (100 * 1024)  /* 100 KB minimum for parallel */
/* Larger <dimension> hints are not preallocated (a stray XFD1048576 would
 * request ~170 GB); the cells found size the storage instead. */
#define XLSX_DIMENSION_HINT_CELLS ((size_t)1 << 27)

static inline bool xlsx_space(char c)
{
    return c == ' ' || c == '\t' || c == '\n' || c == '\r';
}

/* Find a short byte sequence in a buffer (like memmem but always available). */
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

/* End of the markup at lt: the '>' of a tag (quoted attribute values may
 * contain '>'), or the last character of a comment, CDATA section or
 * processing instruction. NULL if the part is truncated. */
static const char *xlsx_markup_end(const char *lt, const char *end)
{
    const char *p = lt + 1;
    if (p < end && (*p == '!' || *p == '?')) {
        const char *close = *p == '?' ? "?>" : NULL;
        if (end - p >= 3 && !memcmp(p, "!--", 3)) close = "-->";
        else if (end - p >= 8 && !memcmp(p, "![CDATA[", 8)) close = "]]>";
        if (!close) return (const char *)memchr(p, '>', (size_t)(end - p));
        size_t n = strlen(close);
        const char *q = xlsx_memfind(p, (size_t)(end - p), close, n);
        return q ? q + n - 1 : NULL;
    }
    char quote = 0;
    for (; p < end; p++) {
        if (quote) {
            if (*p == quote) quote = 0;
        } else if (*p == '"' || *p == '\'') {
            quote = *p;
        } else if (*p == '>') {
            return p;
        }
    }
    return NULL;
}

/* Local name of the tag at lt, after "/" and any "prefix:". *attrs points
 * past the qualified name, where the attributes begin. */
static inline const char *xlsx_tag_name(const char *lt, const char *end, bool *closing,
                                        size_t *len, const char **attrs)
{
    const char *p = lt + 1;
    *closing = p < end && *p == '/';
    if (*closing) p++;
    const char *name = p;
    while (p < end && !xlsx_space(*p) && *p != '>' && *p != '/') {
        if (*p == ':') name = p + 1;
        p++;
    }
    *len = (size_t)(p - name);
    *attrs = p;
    return name;
}

#define XLSX_TAG_IS(name, len, literal) \
    ((len) == sizeof(literal) - 1 && !memcmp((name), (literal), sizeof(literal) - 1))

typedef struct {
    const char *value;
    size_t len;
} xlsx_span;

/* Values of unprefixed attributes of the start tag [p, gt), either quote. */
static void xlsx_attrs(const char *p, const char *gt, const char *const *names,
                       xlsx_span *values, int count)
{
    for (int i = 0; i < count; i++) values[i].value = NULL;
    while (p < gt) {
        while (p < gt && (xlsx_space(*p) || *p == '/')) p++;
        const char *name = p;
        while (p < gt && !xlsx_space(*p) && *p != '=') p++;
        size_t nlen = (size_t)(p - name);
        while (p < gt && xlsx_space(*p)) p++;
        if (p >= gt || *p != '=') return;
        p++;
        while (p < gt && xlsx_space(*p)) p++;
        if (p >= gt || (*p != '"' && *p != '\'')) return;
        char quote = *p++;
        const char *value = p;
        while (p < gt && *p != quote) p++;
        if (p >= gt) return;
        for (int i = 0; i < count; i++) {
            if (strlen(names[i]) == nlen && !memcmp(names[i], name, nlen)) {
                values[i].value = value;
                values[i].len = (size_t)(p - value);
            }
        }
        p++;
    }
}

/* "BC12" or "$BC$12" -> col 54 (0-based), row 12. False unless the whole
 * text is a reference within Excel's limits (XFD, row 1048576). */
static bool xlsx_cell_ref(const char *s, size_t n, int *col, int *row)
{
    size_t i = 0, first;
    int c = 0, r = 0;
    if (i < n && s[i] == '$') i++;
    for (first = i; i < n && i - first < 3; i++) {
        char ch = s[i];
        if (ch >= 'a' && ch <= 'z') ch = (char)(ch - 'a' + 'A');
        if (ch < 'A' || ch > 'Z') break;
        c = c * 26 + (ch - 'A' + 1);
    }
    if (i == first || c > XLSX_MAX_COLUMNS) return false;
    if (i < n && s[i] == '$') i++;
    for (first = i; i < n && i - first < 7 && s[i] >= '0' && s[i] <= '9'; i++)
        r = r * 10 + (s[i] - '0');
    if (i == first || i != n || r < 1 || r > XLSX_MAX_ROWS) return false;
    *col = c - 1;
    *row = r;
    return true;
}

/* Non-negative integer (row number, style or shared string index); -1 if
 * the text is not one. */
static int xlsx_index(const char *s, size_t n)
{
    size_t i = 0, first;
    long long v = 0;
    while (i < n && xlsx_space(s[i])) i++;
    for (first = i; i < n && s[i] >= '0' && s[i] <= '9' && v <= INT32_MAX; i++)
        v = v * 10 + (s[i] - '0');
    while (i < n && xlsx_space(s[i])) i++;
    return i == first || i != n || v > INT32_MAX ? -1 : (int)v;
}

/* Numeric <v> text. Like native, unparseable text keeps its numeric prefix
 * ("1,5" -> 1); text with no number (NaN, words) holds no value. */
static bool xlsx_number(const char *s, size_t n, double *out)
{
    double x;
    if (n > INT32_MAX || !ctools_parse_double_fast(s, (int)n, &x, SV_missval)) {
        char buffer[64];
        char *stop;
        if (n >= sizeof(buffer)) return false;
        memcpy(buffer, s, n);
        buffer[n] = '\0';
        x = strtod(buffer, &stop);
        if (stop == buffer) return false;
    }
    if (!isfinite(x) || x >= SV_missval) return false;
    *out = x;
    return true;
}

/* Grow column-major storage to at least rows x cols, keeping stored cells.
 * On failure the arrays keep their previous, consistent size. */
static bool xlsx_cm_reserve(XLSXContext *ctx, size_t rows, int cols)
{
    size_t old_rows = ctx->cm_num_rows;
    if (rows < old_rows) rows = old_rows;
    if (rows < 1) rows = 1;
    if (cols > XLSX_MAX_COLUMNS) return false;
    if (!ctx->cm_numeric) {
        /* Column pointers cover every Excel column; data columns are
         * allocated only as cells reach them. */
        ctx->cm_numeric = (double **)calloc(XLSX_MAX_COLUMNS, sizeof(double *));
        ctx->cm_types = (uint8_t **)calloc(XLSX_MAX_COLUMNS, sizeof(uint8_t *));
        ctx->cm_formats = (uint8_t **)calloc(XLSX_MAX_COLUMNS, sizeof(uint8_t *));
        if (!ctx->cm_numeric || !ctx->cm_types || !ctx->cm_formats) {
            free(ctx->cm_numeric);
            free(ctx->cm_types);
            free(ctx->cm_formats);
            ctx->cm_numeric = NULL;
            ctx->cm_types = NULL;
            ctx->cm_formats = NULL;
            return false;
        }
    }
    if (rows > old_rows) {
        for (int c = 0; c < ctx->cm_num_cols; c++) {
            double *numeric = (double *)realloc(ctx->cm_numeric[c], rows * sizeof(double));
            if (!numeric) return false;
            ctx->cm_numeric[c] = numeric;
            uint8_t *types = (uint8_t *)realloc(ctx->cm_types[c], rows);
            if (!types) return false;
            ctx->cm_types[c] = types;
            uint8_t *formats = (uint8_t *)realloc(ctx->cm_formats[c], rows);
            if (!formats) return false;
            ctx->cm_formats[c] = formats;
            for (size_t r = old_rows; r < rows; r++) numeric[r] = SV_missval;
            memset(types + old_rows, 0, rows - old_rows);
            memset(formats + old_rows, 0, rows - old_rows);
        }
        ctx->cm_num_rows = rows;
    }
    for (int c = ctx->cm_num_cols; c < cols; c++) {
        double *numeric = (double *)malloc(rows * sizeof(double));
        uint8_t *types = (uint8_t *)calloc(rows, 1);
        uint8_t *formats = (uint8_t *)calloc(rows, 1);
        if (!numeric || !types || !formats) {
            free(numeric);
            free(types);
            free(formats);
            return false;
        }
        for (size_t r = 0; r < rows; r++) numeric[r] = SV_missval;
        ctx->cm_numeric[c] = numeric;
        ctx->cm_types[c] = types;
        ctx->cm_formats[c] = formats;
        ctx->cm_num_cols = c + 1;
    }
    ctx->cm_num_rows = rows;
    ctx->cm_active = true;
    return true;
}

/* Clear stored cells and the used range before a rescan. */
static void xlsx_cm_clear(XLSXContext *ctx)
{
    for (int c = 0; c < ctx->cm_num_cols; c++) {
        for (size_t r = 0; r < ctx->cm_num_rows; r++) ctx->cm_numeric[c][r] = SV_missval;
        memset(ctx->cm_types[c], 0, ctx->cm_num_rows);
        memset(ctx->cm_formats[c], 0, ctx->cm_num_rows);
    }
    ctx->max_col = 0;
    ctx->min_row = INT32_MAX;
    ctx->max_row = 0;
}

/* Mutable scan state. Each parallel worker owns one; the serial scan sets
 * grow so that storage can be enlarged as cells arrive. */
typedef struct {
    XLSXContext *grow;          /* serial scan only */
    ctools_arena *arena;        /* inline strings (ctx->parse_arena or _arena_storage) */
    ctools_arena _arena_storage;
    int max_col, min_row, max_row;  /* used range */
    int cell_rows, cell_cols;       /* storage that the cells require */
    size_t seen, stored;
    int own_first, own_last;    /* parallel workers: rows this chunk may store
                                   (own_first 0 = no restriction) */
    bool overflow;              /* a cell fell outside fixed storage */
    bool implicit_rows;         /* a chunk opened with a row lacking r= */
    bool disorder;              /* a cell named a row owned by another chunk */
    bool failed;                /* out of memory */
} xlsx_scan_local;

static void xlsx_scan_local_init(xlsx_scan_local *loc)
{
    memset(loc, 0, sizeof(*loc));
    loc->max_col = -1;
    loc->min_row = INT32_MAX;
}

/* Fold a finished scan into the used range of ctx and the totals. */
static void xlsx_scan_merge(XLSXContext *ctx, xlsx_scan_local *total,
                            const xlsx_scan_local *loc)
{
    if (loc->max_col > ctx->max_col) ctx->max_col = loc->max_col;
    if (loc->min_row < ctx->min_row) ctx->min_row = loc->min_row;
    if (loc->max_row > ctx->max_row) ctx->max_row = loc->max_row;
    if (loc->cell_rows > total->cell_rows) total->cell_rows = loc->cell_rows;
    if (loc->cell_cols > total->cell_cols) total->cell_cols = loc->cell_cols;
    total->seen += loc->seen;
    total->stored += loc->stored;
    total->overflow |= loc->overflow;
    total->implicit_rows |= loc->implicit_rows;
    total->disorder |= loc->disorder;
    total->failed |= loc->failed;
}

/* Entity-decoded copy of cell text, owned by the scan arena. */
static const char *xlsx_arena_text(xlsx_scan_local *loc, const char *s, size_t n)
{
    char *copy = (char *)ctools_arena_alloc(loc->arena, n + 1);
    if (!copy) return NULL;
    memcpy(copy, s, n);
    copy[n] = '\0';
    xlsx_xml_decode_entities(copy, n);
    return copy;
}

/* Store one value. Fixed storage (parallel workers) records an overflow for
 * the rescan; the serial scan grows storage geometrically. */
static void xlsx_store_cell(const XLSXContext *ctx, xlsx_scan_local *loc, int row, int col,
                            XLSXCellType type, double number, const char *text,
                            uint8_t format)
{
    loc->seen++;
    if ((size_t)row > ctx->cm_num_rows || col >= ctx->cm_num_cols) {
        XLSXContext *grow = loc->grow;
        size_t rows = ctx->cm_num_rows;
        if (!grow) {
            loc->overflow = true;
            return;
        }
        if ((size_t)row > rows) {
            rows += rows / 2;
            if (rows < 1024) rows = 1024;
            if (rows > XLSX_MAX_ROWS) rows = XLSX_MAX_ROWS;
            if (rows < (size_t)row) rows = (size_t)row;
        }
        if (!xlsx_cm_reserve(grow, rows, col >= ctx->cm_num_cols ? col + 1 : 0)) {
            loc->failed = true;
            return;
        }
    }
    size_t i = (size_t)row - 1;
    if (type == XLSX_CELL_STRING)
        memcpy(&ctx->cm_numeric[col][i], &text, sizeof(text));
    else
        ctx->cm_numeric[col][i] = number;
    ctx->cm_types[col][i] = (uint8_t)type;
    ctx->cm_formats[col][i] = format;
    loc->stored++;
    if (col > loc->max_col) loc->max_col = col;
    if (row < loc->min_row) loc->min_row = row;
    if (row > loc->max_row) loc->max_row = row;
}

/* Parse one <c> element whose start tag runs from attrs to gt. The value is
 * the first <v>, or for inline strings the first <t> run outside phonetic
 * guides (<rPh>): native import excel also reads only the first run of rich
 * inline text. Returns the position after </c>, or NULL if the part is
 * truncated or memory runs out. */
static const char *xlsx_scan_cell(const XLSXContext *ctx, xlsx_scan_local *loc,
                                  const char *attrs, const char *gt, const char *end,
                                  int row, int *prev_col, bool skip)
{
    static const char *const names[] = {"r", "t", "s"};
    xlsx_span a[3];
    xlsx_attrs(attrs, gt, names, a, 3);
    int col, ref_row;
    if (a[0].value && xlsx_cell_ref(a[0].value, a[0].len, &col, &ref_row)) row = ref_row;
    else col = *prev_col + 1;  /* no usable r=: the cell follows its predecessor */
    *prev_col = col;
    if (gt[-1] == '/') return gt + 1;  /* <c/> holds no value */

    const XLSXCellRange *range = &ctx->cell_range;
    skip = skip || row < 1 || row > XLSX_MAX_ROWS || col >= XLSX_MAX_COLUMNS ||
           (range->start_col >= 0 && col < range->start_col) ||
           (range->end_col >= 0 && col > range->end_col);
    /* Parallel chunks store only their own rows, so no two workers write the
     * same cell; rows repeated or out of order across chunks (malformed
     * files) send the sheet to the serial rescan. */
    if (!skip && loc->own_first && (row < loc->own_first || row > loc->own_last)) {
        loc->disorder = true;
        skip = true;
    }
    /* Workers only measure cells beyond fixed storage; the rescan stores them. */
    bool outside = !skip && !loc->grow &&
                   ((size_t)row > ctx->cm_num_rows || col >= ctx->cm_num_cols);
    if (outside) {
        loc->overflow = true;
        if (row > loc->cell_rows) loc->cell_rows = row;
        if (col + 1 > loc->cell_cols) loc->cell_cols = col + 1;
    }

    const char *v = NULL, *run = NULL, *p = gt + 1;
    size_t vlen = 0, run_len = 0;
    bool in_is = false;
    int phonetic = 0;
    for (;;) {
        const char *lt = (const char *)memchr(p, '<', (size_t)(end - p));
        if (!lt) return NULL;
        bool closing;
        size_t len;
        const char *tag_attrs;
        const char *name = xlsx_tag_name(lt, end, &closing, &len, &tag_attrs);
        const char *tag_end = xlsx_markup_end(lt, end);
        if (!tag_end) return NULL;
        p = tag_end + 1;
        if (closing) {
            if (XLSX_TAG_IS(name, len, "c")) break;
            if (XLSX_TAG_IS(name, len, "is")) in_is = false;
            else if (XLSX_TAG_IS(name, len, "rPh") && phonetic) phonetic--;
            continue;
        }
        if (tag_end[-1] == '/') continue;  /* empty element */
        bool value = XLSX_TAG_IS(name, len, "v") && !v;
        bool first_run = XLSX_TAG_IS(name, len, "t") && in_is && !phonetic && !run;
        if (XLSX_TAG_IS(name, len, "is")) in_is = true;
        else if (XLSX_TAG_IS(name, len, "rPh")) phonetic++;
        if (value || first_run) {
            /* Character data never contains '<' */
            const char *text_end = (const char *)memchr(p, '<', (size_t)(end - p));
            if (!text_end) return NULL;
            if (value) {
                v = p;
                vlen = (size_t)(text_end - p);
            } else {
                run = p;
                run_len = (size_t)(text_end - p);
            }
            p = text_end;
        }
    }
    if (skip || outside) return p;

    const char *t = a[1].value;
    size_t tn = t ? a[1].len : 0;
    int style = a[2].value ? xlsx_index(a[2].value, a[2].len) : -1;
    uint8_t format = style >= 0 && style < ctx->num_styles ? ctx->excel_formats[style] : 0;
    if (tn >= 9 && t[0] == 'i') {  /* inlineStr */
        if (!run_len) return p;
        const char *text = xlsx_arena_text(loc, run, run_len);
        if (!text) {
            loc->failed = true;
            return NULL;
        }
        xlsx_store_cell(ctx, loc, row, col, XLSX_CELL_STRING, 0, text, format);
    } else if (!v || !vlen) {
        return p;
    } else if (tn == 1 && t[0] == 's') {
        int index = xlsx_index(v, vlen);
        if (index >= 0)
            xlsx_store_cell(ctx, loc, row, col, XLSX_CELL_SHARED_STRING, index, NULL, format);
    } else if (tn == 3 && !memcmp(t, "str", 3)) {  /* formula text result */
        const char *text = xlsx_arena_text(loc, v, vlen);
        if (!text) {
            loc->failed = true;
            return NULL;
        }
        xlsx_store_cell(ctx, loc, row, col, XLSX_CELL_STRING, 0, text, format);
    } else if (tn == 1 && t[0] == 'b') {
        xlsx_store_cell(ctx, loc, row, col, XLSX_CELL_BOOLEAN, v[0] == '1' ? 1.0 : 0.0,
                        NULL, format);
    } else if (!(tn == 1 && t[0] == 'e')) {  /* error values import as missing */
        double number;
        if (!xlsx_number(v, vlen, &number)) return p;
        XLSXCellType type = XLSX_CELL_NUMBER;
        if (style >= 0 && style < ctx->num_styles && ctx->date_styles[style])
            type = ctx->date_styles[style] == 2 ? XLSX_CELL_DATETIME : XLSX_CELL_DATE;
        xlsx_store_cell(ctx, loc, row, col, type, number, NULL, format);
    }
    return loc->failed ? NULL : p;
}

/* Scan the rows in [p, end) inside <sheetData>. Only the first region may
 * number rows that lack r= (implicit numbers depend on all earlier rows):
 * any other region must open with a numbered row, or it sets implicit_rows
 * so that the caller rescans serially. */
static void xlsx_scan_rows(const XLSXContext *ctx, xlsx_scan_local *loc,
                           const char *p, const char *end, bool first)
{
    static const char *const names[] = {"r"};
    const XLSXCellRange *range = &ctx->cell_range;
    int row = 0, prev_col = -1;
    bool numbered = first, skip = false;
    while (p < end) {
        const char *lt = (const char *)memchr(p, '<', (size_t)(end - p));
        if (!lt) return;
        bool closing;
        size_t len;
        const char *attrs;
        const char *name = xlsx_tag_name(lt, end, &closing, &len, &attrs);
        const char *gt = xlsx_markup_end(lt, end);
        if (!gt) return;
        p = gt + 1;
        if (closing) {
            if (XLSX_TAG_IS(name, len, "sheetData")) return;
        } else if (XLSX_TAG_IS(name, len, "c")) {
            p = xlsx_scan_cell(ctx, loc, attrs, gt, end, row, &prev_col, skip);
            if (!p) return;
        } else if (XLSX_TAG_IS(name, len, "row")) {
            xlsx_span r;
            xlsx_attrs(attrs, gt, names, &r, 1);
            int number = r.value ? xlsx_index(r.value, r.len) : -1;
            if (number >= 1 && number <= XLSX_MAX_ROWS) {
                row = number;
                numbered = true;
            } else if (numbered) {
                row++;
            } else {
                loc->implicit_rows = true;
                return;
            }
            prev_col = -1;
            if (range->end_row > 0 && row > range->end_row) return;
            skip = range->start_row > 0 && row < range->start_row;
            if (!skip) {
                if (row < loc->min_row) loc->min_row = row;
                if (row > loc->max_row) loc->max_row = row;
            }
        }
    }
}

/* ============================================================================
 * Parallel Worksheet Scanner
 *
 * Large <sheetData> regions are split at <row> boundaries and scanned on the
 * thread pool. Workers write to disjoint rows of fixed storage and report
 * the extents of cells that did not fit, so the caller can grow the storage
 * and scan again.
 * ============================================================================ */

typedef struct {
    const XLSXContext *ctx;       /* read-only context (cm arrays, shared strings, etc.) */
    const char *chunk_start;      /* start of this thread's region (at a <row boundary) */
    const char *chunk_end;        /* end of this thread's region */
    bool first;                   /* region starts at the beginning of <sheetData> */
    xlsx_scan_local local;        /* thread-local mutable state */
} xlsx_parallel_chunk;

static void *xlsx_parallel_scan_worker(void *arg)
{
    xlsx_parallel_chunk *chunk = (xlsx_parallel_chunk *)arg;
    xlsx_scan_rows(chunk->ctx, &chunk->local, chunk->chunk_start, chunk->chunk_end,
                   chunk->first);
    return NULL;
}

/* Next <row> start tag (any prefix) at or after p. */
static const char *xlsx_next_row(const char *p, const char *end)
{
    while (p < end && (p = (const char *)memchr(p, '<', (size_t)(end - p)))) {
        bool closing;
        size_t len;
        const char *attrs;
        const char *name = xlsx_tag_name(p, end, &closing, &len, &attrs);
        if (!closing && XLSX_TAG_IS(name, len, "row")) return p;
        p++;
    }
    return NULL;
}

/* Scan [begin, end) on the thread pool. Returns false, having scanned
 * nothing, when the region is too small or the pool is unavailable. */
static bool xlsx_scan_parallel(XLSXContext *ctx, const char *begin, const char *end,
                               xlsx_scan_local *total)
{
    ctools_persistent_pool *pool = ctools_get_global_pool();
    if (!pool) return false;
    size_t sheet_len = (size_t)(end - begin);
    if (sheet_len < XLSX_PARALLEL_MIN_BYTES) return false;
    int nthreads = (int)pool->num_workers;
    if (nthreads < 2) return false;
    if (nthreads > 16) nthreads = 16;

    const char *boundaries[17];
    int nchunks = 1;
    boundaries[0] = begin;
    size_t target_chunk_size = sheet_len / (size_t)nthreads;
    for (int i = 1; i < nthreads; i++) {
        const char *found = xlsx_next_row(begin + target_chunk_size * (size_t)i, end);
        if (found && found > boundaries[nchunks - 1]) boundaries[nchunks++] = found;
    }
    boundaries[nchunks] = end;
    if (nchunks < 2) return false;  /* Not enough rows for parallel */

    /* Chunk i owns the rows from its first row number up to the next chunk's
     * first row number - 1. A chunk whose first row lacks r= reports
     * implicit_rows itself. */
    int first_row[17];
    first_row[0] = 1;
    for (int i = 1; i < nchunks; i++) {
        static const char *const rnames[] = {"r"};
        bool closing;
        size_t len;
        const char *attrs;
        xlsx_tag_name(boundaries[i], end, &closing, &len, &attrs);
        const char *gt = xlsx_markup_end(boundaries[i], end);
        xlsx_span r = {0};
        if (gt) xlsx_attrs(attrs, gt, rnames, &r, 1);
        int number = r.value ? xlsx_index(r.value, r.len) : -1;
        first_row[i] = (number >= 1 && number <= XLSX_MAX_ROWS) ? number : 1;
    }

    xlsx_parallel_chunk *chunks = (xlsx_parallel_chunk *)calloc(
        (size_t)nchunks, sizeof(xlsx_parallel_chunk));
    if (!chunks) return false;
    for (int i = 0; i < nchunks; i++) {
        chunks[i].ctx = ctx;
        chunks[i].chunk_start = boundaries[i];
        chunks[i].chunk_end = boundaries[i + 1];
        chunks[i].first = i == 0;
        xlsx_scan_local_init(&chunks[i].local);
        chunks[i].local.own_first = first_row[i];
        chunks[i].local.own_last = i + 1 < nchunks ? first_row[i + 1] - 1 : XLSX_MAX_ROWS;
        ctools_arena_init(&chunks[i].local._arena_storage, 0);
        chunks[i].local.arena = &chunks[i].local._arena_storage;
    }

    bool ran = ctools_persistent_pool_submit_batch(pool, xlsx_parallel_scan_worker,
                                                   chunks, (size_t)nchunks,
                                                   sizeof(xlsx_parallel_chunk)) == 0;
    if (ran) ctools_persistent_pool_wait(pool);

    for (int i = 0; i < nchunks; i++) {
        if (ran) xlsx_scan_merge(ctx, total, &chunks[i].local);
        /* Retain inline strings until column caches have copied them. */
        ctools_arena *arena = chunks[i].local.arena;
        if (arena->first) {
            if (ctx->parse_arena.current) ctx->parse_arena.current->next = arena->first;
            else ctx->parse_arena.first = arena->first;
            ctx->parse_arena.current = arena->current;
            ctx->parse_arena.total_allocated += arena->total_allocated;
        }
    }
    free(chunks);
    return ran;
}

/* Scan the content of <sheetData>, [begin, end). Parallel workers write to
 * fixed storage; cells beyond it (no, degenerate or understated <dimension>)
 * trigger a second pass once storage covers the extents they reported. A
 * chunk that cannot number its rows falls back to the serial scan, which
 * grows storage as it goes. */
static ST_retcode xlsx_scan_sheet(XLSXContext *ctx, const char *begin, const char *end)
{
    if (ctx->dimension_rows > 0 && ctx->dimension_cols > 0 &&
        (size_t)ctx->dimension_rows * (size_t)ctx->dimension_cols <= XLSX_DIMENSION_HINT_CELLS)
        xlsx_cm_reserve(ctx, (size_t)ctx->dimension_rows, ctx->dimension_cols);

    xlsx_scan_local total;
    for (int pass = 0; pass < 2; pass++) {
        xlsx_scan_local_init(&total);
        if (!xlsx_scan_parallel(ctx, begin, end, &total)) break;
        if (total.failed) return 909;
        if (total.implicit_rows || total.disorder) {
            xlsx_cm_clear(ctx);
            break;
        }
        if (!total.overflow) {
            ctx->cells_seen = total.seen;
            ctx->cells_stored = total.stored;
            return 0;
        }
        if (!xlsx_cm_reserve(ctx, (size_t)total.cell_rows, total.cell_cols)) return 909;
    }

    xlsx_scan_local loc;
    xlsx_scan_local_init(&loc);
    loc.grow = ctx;
    loc.arena = &ctx->parse_arena;
    xlsx_scan_rows(ctx, &loc, begin, end, true);
    xlsx_scan_local_init(&total);
    xlsx_scan_merge(ctx, &total, &loc);
    if (loc.failed) return 909;
    ctx->cells_seen = loc.seen;
    ctx->cells_stored = loc.stored;
    return 0;
}

/* <dimension ref="A1:D20">, or one cell ("B3"): a hint for the first
 * allocation only. */
static void xlsx_dimension_hint(XLSXContext *ctx, const char *ref, size_t len)
{
    const char *colon = (const char *)memchr(ref, ':', len);
    const char *last = colon ? colon + 1 : ref;
    int col, row;
    if (!xlsx_cell_ref(last, (size_t)(ref + len - last), &col, &row)) return;
    ctx->dimension_cols = col + 1;
    ctx->dimension_rows = row;
}

/* </sheetData> sits near the end of the part: search backwards, stopping at
 * the first cell or row markup (a truncated part then ends the region). */
static const char *xlsx_sheetdata_end(const char *begin, const char *end)
{
    for (const char *p = end; p > begin;) {
        p--;
        if (*p != '<') continue;
        bool closing;
        size_t len;
        const char *attrs;
        const char *name = xlsx_tag_name(p, end, &closing, &len, &attrs);
        if (closing && XLSX_TAG_IS(name, len, "sheetData")) return p;
        if (XLSX_TAG_IS(name, len, "row") || XLSX_TAG_IS(name, len, "c") ||
            XLSX_TAG_IS(name, len, "v"))
            break;
    }
    return end;
}

/* Scan a complete worksheet part: the <dimension> hint, then <sheetData>. */
static ST_retcode xlsx_scan_buffer(XLSXContext *ctx, const char *buf, size_t len)
{
    static const char *const names[] = {"ref"};
    const char *p = buf, *end = buf + len;
    while (p < end) {
        const char *lt = (const char *)memchr(p, '<', (size_t)(end - p));
        if (!lt) return 0;
        bool closing;
        size_t nlen;
        const char *attrs;
        const char *name = xlsx_tag_name(lt, end, &closing, &nlen, &attrs);
        const char *gt = xlsx_markup_end(lt, end);
        if (!gt) return 0;
        p = gt + 1;
        if (closing) continue;
        if (XLSX_TAG_IS(name, nlen, "dimension")) {
            xlsx_span ref;
            xlsx_attrs(attrs, gt, names, &ref, 1);
            if (ref.value) xlsx_dimension_hint(ctx, ref.value, ref.len);
        } else if (XLSX_TAG_IS(name, nlen, "sheetData")) {
            if (gt[-1] == '/') return 0;  /* <sheetData/>: an empty sheet */
            return xlsx_scan_sheet(ctx, p, xlsx_sheetdata_end(p, end));
        }
    }
    return 0;
}

static const char *xlsx_sheet_kind_name(int kind)
{
    return kind == XLSX_SHEET_CHART    ? "a chart sheet"
         : kind == XLSX_SHEET_DIALOG   ? "a dialog sheet"
         : kind == XLSX_SHEET_MACRO    ? "a macro sheet"
                                       : "not a worksheet";
}

ST_retcode xlsx_parse_worksheet(XLSXContext *ctx)
{
    xlsx_zip_archive *zip = (xlsx_zip_archive *)ctx->zip_archive;

    double t_start = ctools_timer_seconds();

    if (ctx->selected_sheet < 0 || ctx->selected_sheet >= ctx->num_sheets) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "workbook contains no worksheets");
        return 610;
    }
    const XLSXSheetInfo *info = &ctx->sheets[ctx->selected_sheet];
    /* Native finds no worksheet of that name (r(601)); the caller treats a
     * default selection without any worksheet as an empty sheet. */
    if (info->kind != XLSX_SHEET_WORKSHEET) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "%s is %s; only worksheets can be imported", info->name,
                 xlsx_sheet_kind_name(info->kind));
        return 601;
    }

    size_t file_idx = xlsx_zip_locate_file(zip, info->path);
    if (file_idx == (size_t)-1) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to find worksheet: %s", info->path);
        return 610;
    }

    /* Inflate the whole part (libdeflate) and scan it in memory */
    size_t xml_size = 0;
    char *xml = (char *)xlsx_zip_extract_direct(zip, file_idx, &xml_size);
    if (!xml) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Failed to extract worksheet: %s", info->path);
        return 610;
    }
    ST_retcode rc = xlsx_scan_buffer(ctx, xml, xml_size);
    free(xml);

    /* Never report an empty or partial sheet as a successful import */
    if (rc == 909) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Insufficient memory for worksheet cells");
    } else if (!rc && ctx->cells_stored < ctx->cells_seen) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "%zu worksheet cells could not be stored",
                 ctx->cells_seen - ctx->cells_stored);
        rc = 610;
    }

    ctx->time_worksheet = ctools_timer_seconds() - t_start;

    return rc;
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
 * XLSX Column Store Worker
 *
 * No longer called by the plugin (Excel import goes through src/io/cio.c);
 * kept because validation/test_p2_native.py exercises its SPI failure paths.
 * ============================================================================ */

typedef struct {
    CImportColumnInfo *col;
    CImportColumnCache *cache;
    ST_int var;
    ST_int stata_nobs;
    ST_retcode error;
} xlsx_store_task;

__attribute__((unused))
static void *xlsx_store_worker(void *arg)
{
    xlsx_store_task *task = (xlsx_store_task *)arg;
    task->error = 0;
    int is_string = task->col->type == CIMPORT_COL_STRING;
    if (task->var < 1 || task->var > SF_nvars() ||
        task->cache->count > (size_t)task->stata_nobs ||
        SF_var_is_string(task->var) != is_string ||
        (is_string && SF_var_is_strl(task->var))) {
        task->error = 198;
        return (void *)1;
    }
    if (task->col->type == CIMPORT_COL_STRING) {
        for (size_t r = 0; r < task->cache->count && (ST_int)(r + 1) <= task->stata_nobs; r++) {
            char *str = task->cache->string_data[r] ? task->cache->string_data[r] : (char *)"";
            task->error = SF_sstore(task->var, (ST_int)(r + 1), str);
            if (task->error) return (void *)1;
        }
    } else {
        for (size_t r = 0; r < task->cache->count && (ST_int)(r + 1) <= task->stata_nobs; r++) {
            task->error = (_stata_)->safestore(task->var, (ST_int)(r + 1), task->cache->numeric_data[r]);
            if (task->error) return (void *)1;
        }
    }
    return NULL;
}

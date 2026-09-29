/*
 * cimport_impl.c - High-Performance CSV Import for Stata
 *
 * Part of the ctools suite. Replaces "import delimited" with a
 * C-accelerated implementation featuring:
 *   - Multi-threaded parallel parsing
 *   - SIMD-accelerated newline/delimiter scanning
 *   - Custom fast float parser
 *   - Arena allocator for memory efficiency
 *   - Column-major caching for fast SPI loading
 *
 * Based on fastimport by Claude (Anthropic)
 * License: MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdatomic.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <pthread.h>

/* Platform-specific includes */
#ifdef _WIN32
    #include <windows.h>
    #include <io.h>
    #define CIMPORT_WINDOWS 1
    #define strcasecmp _stricmp
#else
    #include <unistd.h>
    #include <strings.h>  /* For strcasecmp */
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "cimport_impl.h"
#include "../io/cio.h"

/* Helper modules (cimport_context.h includes arena, parse, and mmap headers) */
#include "cimport_context.h"
#include "cimport_mmap.h"
#include "ctools_threads.h"

/* High-resolution timing */
/* Timer: use ctools_timer_ms() directly from ctools_timer.h */

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/* Local: cimport may use more threads for parsing than general I/O */
#define CIMPORT_MAX_THREADS      32

/* Global cached context */
static CImportContext *g_cimport_ctx = NULL;
static ST_retcode g_cimport_parse_status = 601;

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

static void cimport_display_msg(const char *msg) {
    SF_display((char *)msg);
}

static void cimport_display_error(const char *msg) {
    SF_error((char *)msg);
}

static char cimport_extract_quote(const CImportContext *ctx) {
  return ctx->retain_quotes                             ? '\0'
         : ctx->strip_quotes_all                        ? '\1'
         : ctx->bindquotes == CIMPORT_BINDQUOTES_NOBIND ? '\2'
                                                        : ctx->quote_char;
}

/* Inference and cached values must see the same decoded quote policy. */
static bool cimport_context_number(CImportContext *ctx, CImportFieldRef *field,
                                   double *value) {
  int len = (int)CIMPORT_FIELD_LENGTH(*field);
  const char *src = ctx->file_data + field->offset;
  char stack[256], *decoded = NULL;
  if (field->length & CIMPORT_FIELD_QUOTED_FLAG) {
    decoded =
        (size_t)len + 1 <= sizeof(stack) ? stack : malloc((size_t)len + 1);
    if (!decoded) {
      atomic_store(&ctx->error_code, 920);
      return false;
    }
    len = cimport_extract_field_fast(ctx->file_data, field, decoded, len + 1,
                                     cimport_extract_quote(ctx));
    src = decoded;
  }
  bool valid = cimport_field_looks_unquoted_numeric_sep(
                   src, len, ctx->decimal_separator, ctx->group_separator) &&
               cimport_parse_unquoted_number(src, len, value, SV_missval,
                                             ctx->decimal_separator,
                                             ctx->group_separator);
  if (decoded && decoded != stack)
    free(decoded);
  return valid;
}
static bool cimport_context_numeric(CImportContext *ctx,
                                    CImportFieldRef *field) {
  double value;
  return cimport_context_number(ctx, field, &value);
}

/* Record a warning for unmatched quote (thread-safe) */
static void cimport_record_unmatched_quote(CImportContext *ctx,
                                           size_t row_num) {
  pthread_mutex_lock(&ctx->warning_mutex);
  if (ctx->num_unmatched_quote_warnings < CIMPORT_MAX_WARNINGS) {
    ctx->unmatched_quote_rows[ctx->num_unmatched_quote_warnings] = row_num;
  }
  ctx->num_unmatched_quote_warnings++;
  pthread_mutex_unlock(&ctx->warning_mutex);
}

/* Skip N rows from the start of data, returning pointer to start of row N+1 */
static const char *cimport_skip_rows(const char *data, const char *end, int n_rows, char quote_char, CImportBindQuotesMode bindquotes) {
    const char *ptr = data;
    for (int i = 0; i < n_rows && ptr < end; i++) {
        ptr = cimport_find_next_row(ptr, end, quote_char, bindquotes);
    }
    return ptr;
}

/* Number of `byte` bytes in [p, end), eight at a time. */
static size_t cimport_count_byte(const char *p, const char *end, char byte) {
    const uint64_t low7 = 0x7F7F7F7F7F7F7F7FULL;
    const uint64_t pattern = 0x0101010101010101ULL * (unsigned char)byte;
    size_t count = 0;
    while (end - p >= 8) {
        uint64_t w;
        memcpy(&w, p, 8);
        uint64_t x = w ^ pattern;                  /* zero byte where w == byte */
        uint64_t t = ~(((x & low7) + low7) | x | low7); /* 0x80 in each zero byte */
        count += (size_t)__builtin_popcountll(t);
        p += 8;
    }
    while (p < end) count += (*p++ == byte);
    return count;
}

/* Native rowrange() counts physical file lines: blank lines and newlines inside
 * strict quoted fields count too. The cursor holds the number of '\n' bytes
 * before pos and is only moved forward. */
typedef struct {
    const char *pos;
    size_t lines;
} CImportLineCursor;

/* 1-based line holding the byte at p (p at or after the cursor). */
static size_t cimport_line_at(CImportLineCursor *cursor, const char *p) {
    if (p > cursor->pos) {
        cursor->lines += cimport_count_byte(cursor->pos, p, '\n');
        cursor->pos = p;
    }
    return cursor->lines + 1;
}

/* Last line of the row [start, next): the line of its terminator. */
static size_t cimport_row_last_line(CImportLineCursor *cursor, const char *start, const char *next) {
    return cimport_line_at(cursor, next > start && next[-1] == '\n' ? next - 1 : next);
}

/* rowrange() bounds as file lines (0 when open); the ado passes them as given. */
static void cimport_row_range(size_t *first, size_t *last) {
    char text[32] = "";
    SF_macro_use("___cimport_startrow", text, sizeof(text) - 1);
    *first = (size_t)strtoull(text, NULL, 10);
    text[0] = '\0';
    SF_macro_use("___cimport_endrow", text, sizeof(text) - 1);
    *last = (size_t)strtoull(text, NULL, 10);
}

/* Native reads a row within rowrange() when the row ends on or after the first
 * line and the row read before it (or the header) ended before the last line. */
static bool cimport_line_in_range(size_t first, size_t last, size_t prev, size_t line) {
    return (!last || prev < last) && (!first || line >= first);
}

/* Empty rows: a blank line, or a CRLF blank line whose only field is "\r". */
static bool cimport_row_is_empty(const char *file_base, const CImportFieldRef *fields, int num_fields) {
    if (num_fields == 0) return true;
    if (num_fields > 1) return false;
    uint32_t len = CIMPORT_FIELD_LENGTH(fields[0]);
    return len == 0 || (len == 1 && file_base[fields[0].offset] == '\r');
}

/* Display warnings for unmatched quotes (Stata-compatible format) */
static void cimport_display_warnings(CImportContext *ctx) {
    if (ctx->num_unmatched_quote_warnings == 0) return;

    char msg[600];
    int to_display = (ctx->num_unmatched_quote_warnings < CIMPORT_MAX_WARNINGS)
                     ? ctx->num_unmatched_quote_warnings
                     : CIMPORT_MAX_WARNINGS;

    for (int i = 0; i < to_display; i++) {
        snprintf(msg, sizeof(msg),
            "Note: Unmatched quote while processing row %zu; this can be due to a formatting\n"
            "    problem in the file or because a quoted data element spans multiple\n"
            "    lines. You should carefully inspect your data after importing. Consider\n"
            "    using option bindquote(strict) if quoted data spans multiple lines or\n"
            "    option bindquote(nobind) if quotes are not used for binding data.\n",
            ctx->unmatched_quote_rows[i]);
        cimport_display_msg(msg);
    }

    if (ctx->num_unmatched_quote_warnings > CIMPORT_MAX_WARNINGS) {
        snprintf(msg, sizeof(msg),
            "(... %d additional unmatched quote warnings not shown)\n",
            ctx->num_unmatched_quote_warnings - CIMPORT_MAX_WARNINGS);
        cimport_display_msg(msg);
    }
}

/* ============================================================================
 * Variable Name Sanitization
 * ============================================================================ */

static void cimport_sanitize_varname(const char *input, int len, char *output) {
    int j = 0;

    while (len > 0 && (cimport_is_whitespace(*input) || *input == '"')) {
        input++;
        len--;
    }
    while (len > 0 && (cimport_is_whitespace(input[len-1]) || input[len-1] == '"')) {
        len--;
    }

    for (int i = 0; i < len && j < CTOOLS_MAX_VARNAME_LEN; i++) {
        char c = input[i];
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_' || (unsigned char)c >= 128) {
            output[j++] = c;
        } else if (c >= '0' && c <= '9' && j > 0) {
            output[j++] = c;
        }
    }

    if (j == 0) {
        output[0] = '\0';
    }

    output[j] = '\0';
}

/* ============================================================================
 * Type Inference
 * ============================================================================ */

static CImportNumericSubtype cimport_determine_numeric_subtype(CImportColumnInfo *col, bool use_double_for_decimals) {
    if (col->min_value < -2147483647.0 || col->max_value > 2147483620.0) {
        return CIMPORT_NUM_DOUBLE;
    }

    if (!col->is_integer) {
        /* The caller's numeric type preference applies to non-integers. */
        if (use_double_for_decimals) {
            return CIMPORT_NUM_DOUBLE;
        }
        return CIMPORT_NUM_FLOAT;
    }

    if (col->min_value >= -127 && col->max_value <= 100) return CIMPORT_NUM_BYTE;

    /* Every numeric value is checked before choosing integer storage. */
    if (col->min_value >= -32767 && col->max_value <= 32740) {
        return CIMPORT_NUM_INT;
    }

    return CIMPORT_NUM_LONG;
}

/* ============================================================================
 * Header and Type Inference
 * ============================================================================ */

/* The header is the first row after the varnames() skip that is not empty in
 * emptylines(skip) mode, as native reads it; the chunk parsers skip the same
 * empty rows before taking row 0. Returns its field count, or 0 when only
 * empty rows remain; *row_start and *next_row bound the row. */
static int cimport_parse_header_row(const CImportContext *ctx, CImportFieldRef *fields,
                                    const char **row_start, const char **next_row) {
    const char *end = ctx->file_data + ctx->file_size;
    const char *row = cimport_skip_rows(ctx->file_data, end, ctx->skip_rows, ctx->quote_char, ctx->bindquotes);
    const char *next = end;
    int count = 0;
    while (row < end) {
        count = cimport_parse_row_fast(row, end, ctx->delimiter, ctx->quote_char, fields, CTOOLS_MAX_COLUMNS,
                                       ctx->file_data, ctx->bindquotes, &next, NULL);
        if (ctx->emptylines_mode != CIMPORT_EMPTYLINES_SKIP ||
            !cimport_row_is_empty(ctx->file_data, fields, count)) break;
        row = next;
        count = 0;
    }
    if (row_start) *row_start = row;
    if (next_row) *next_row = row < end ? next : end;
    return count;
}

static int cimport_parse_header(CImportContext *ctx) {
    CImportFieldRef fields[CTOOLS_MAX_COLUMNS];
    const char *end = ctx->file_data + ctx->file_size;
    const char *header_start, *next;

    ctx->num_columns = cimport_parse_header_row(ctx, fields, &header_start, &next);
    /* Only empty rows: native imports nothing, like an empty file. */
    if (ctx->num_columns == 0) return 0;

    /* Auto headers are inferred from a text-to-number transition, as native
     * does: some first-row cell is text while its column's data are numbers.
     * The data are the rows native reads for inference: empty rows are skipped
     * and rowrange() lines apply after the header, so text with no such rows
     * is still a header. Merely seeing text in the first row would drop data
     * from all-string files. */
    char autoheader[8] = "";
    SF_macro_use("___cimport_autoheader", autoheader, sizeof(autoheader) - 1);
    if (strcmp(autoheader, "1") == 0) {
        ctx->has_header = false;
        CImportFieldRef next_fields[CTOOLS_MAX_COLUMNS];
        bool candidates[CTOOLS_MAX_COLUMNS];
        for (int i = 0; i < ctx->num_columns; i++) {
            int len = CIMPORT_FIELD_LENGTH(fields[i]);
            candidates[i] = len && !cimport_context_numeric(ctx,&fields[i]);
        }
        size_t first, last, prev = 0;
        cimport_row_range(&first, &last);
        CImportLineCursor cursor = {ctx->file_data, 0};
        if (first || last) prev = cimport_row_last_line(&cursor, header_start, next);
        for (int r = 0; r < 100 && next < end && !(last && prev >= last); ) {
            const char *row = next;
            if (first) {
                /* Rows ending before the first line are passed over unparsed. */
                const char *after = cimport_find_next_row(row, end, ctx->quote_char, ctx->bindquotes);
                if (cimport_row_last_line(&cursor, row, after) < first) {
                    next = after;
                    continue;
                }
            }
            int count = cimport_parse_row_fast(row, end, ctx->delimiter, ctx->quote_char,
                        next_fields, CTOOLS_MAX_COLUMNS, ctx->file_data, ctx->bindquotes, &next, NULL);
            if (ctx->emptylines_mode == CIMPORT_EMPTYLINES_SKIP &&
                cimport_row_is_empty(ctx->file_data, next_fields, count)) continue;
            if (first || last) prev = cimport_row_last_line(&cursor, row, next);
            r++;
            for (int i = 0; i < ctx->num_columns && i < count; i++) {
                if (!cimport_context_numeric(ctx,&next_fields[i])) candidates[i] = false;
            }
        }
        for (int i = 0; i < ctx->num_columns; i++) {
            if (candidates[i]) { ctx->has_header = true; break; }
        }
    }

    ctx->columns = ctools_safe_calloc2(ctx->num_columns, sizeof(CImportColumnInfo));
    if (!ctx->columns) return -1;

    char name_buf[CTOOLS_MAX_STRING_LEN];
    for (int i = 0; i < ctx->num_columns; i++) {
        if (ctx->has_header) {
            cimport_extract_field_fast(ctx->file_data, &fields[i], name_buf, sizeof(name_buf), cimport_extract_quote(ctx));
            cimport_sanitize_varname(name_buf, strlen(name_buf), ctx->columns[i].name);
            if (!ctx->columns[i].name[0])
                snprintf(ctx->columns[i].name, sizeof(ctx->columns[i].name), "v%d", i + 1);
        } else {
            snprintf(ctx->columns[i].name, CTOOLS_MAX_VARNAME_LEN, "v%d", i + 1);
        }
        ctx->columns[i].type = CIMPORT_COL_UNKNOWN;
        ctx->columns[i].max_strlen = 0;
    }

    return 0;
}

/* Helper to check if a column is in the force list */
static bool cimport_col_in_list(int col_1based, int *list, int count) {
    for (int i = 0; i < count; i++) {
        if (list[i] == 0 || list[i] == col_1based) return true;
    }
    return false;
}

/* String overrides take precedence over numeric overrides. */
static int cimport_column_force(const CImportContext *ctx, int col) {
    if (ctx->force_string_cols && cimport_col_in_list(col + 1, ctx->force_string_cols, ctx->num_force_string))
        return CIMPORT_FORCE_STRING;
    if (ctx->force_numeric_cols && cimport_col_in_list(col + 1, ctx->force_numeric_cols, ctx->num_force_numeric))
        return CIMPORT_FORCE_NUMERIC;
    return CIMPORT_FORCE_NONE;
}

/* Resolve overrides for columns [from, to) so parse workers index a table. */
static int cimport_resolve_col_force(CImportContext *ctx, int from, int to) {
    unsigned char *force = realloc(ctx->col_force, (size_t)(to > 0 ? to : 1));
    if (!force) return -1;
    ctx->col_force = force;
    for (int col = from; col < to; col++) force[col] = (unsigned char)cimport_column_force(ctx, col);
    return 0;
}

static void cimport_init_col_stats(CImportColumnParseStats *stats, int count) {
    memset(stats, 0, (size_t)count * sizeof(*stats));
    for (int i = 0; i < count; i++) {
        stats[i].num_min = DBL_MAX;
        stats[i].num_max = -DBL_MAX;
    }
}

/* Gather every per-field statistic type inference needs while the parse
 * workers already have the field in cache: the numeric test, decoded string
 * widths and bytes, and numeric extrema while the column can still be
 * numeric. Returns the field's value for the store when it is a plain
 * number (identical to cimport_context_number), else CIMPORT_NO_VALUE. */
static inline uint64_t cimport_accumulate_field(CImportContext *ctx, CImportColumnParseStats *stats,
                                                const CImportFieldRef *field, int force) {
    uint32_t flen = CIMPORT_FIELD_LENGTH(*field);
    if ((int)flen > stats->max_field_len) stats->max_field_len = flen;
    if (flen == 0) return CIMPORT_NO_VALUE;

    const char *src = ctx->file_data + field->offset;
    bool quoted = (field->length & CIMPORT_FIELD_QUOTED_FLAG) != 0;
    double val;
    bool have_value = false;
    int decoded;

    if (!quoted && ctx->plain_numbers && cimport_parse_plain_decimal(src, (int)flen, &val)) {
        /* Plain numbers pass the numeric test; one scan suffices. */
        have_value = true;
        decoded = (int)flen - (src[flen-1] == '\r');
    } else {
        if (quoted) stats->has_quotes = true;
        if (!stats->seen_string && !cimport_context_numeric(ctx, (CImportFieldRef *)field))
            stats->seen_string = true;
        if (ctx->retain_quotes || quoted) {
            decoded = cimport_extract_field_fast(ctx->file_data, (CImportFieldRef *)field, NULL, INT_MAX,
                                                 cimport_extract_quote(ctx));
        } else {
            decoded = (int)flen;
            while (decoded > 0 && (src[decoded-1] == '\r' || src[decoded-1] == '\n')) decoded--;
        }
    }
    if (decoded > stats->max_decoded_len) stats->max_decoded_len = decoded;
    stats->decoded_bytes += (uint64_t)decoded;

    uint64_t cached = CIMPORT_NO_VALUE;
    if (have_value) memcpy(&cached, &val, sizeof(cached));

    if (force == CIMPORT_FORCE_STRING || (force != CIMPORT_FORCE_NUMERIC && stats->seen_string)) return cached;
    if (!have_value && (!cimport_context_number(ctx, (CImportFieldRef *)field, &val) ||
                        !isfinite(val) || val >= SV_missval)) {
        return cached;
    }
    if (floor(val) != val) stats->num_nonint = true;
    if (val < stats->num_min) stats->num_min = val;
    if (val > stats->num_max) stats->num_max = val;
    return cached;
}

static void cimport_infer_column_types(CImportContext *ctx) {
    size_t total_data_rows = 0;
    for (int c = 0; c < ctx->num_chunks; c++) {
        size_t start_row = (c == 0 && ctx->has_header) ? 1 : 0;
        total_data_rows += ctx->chunks[c].num_rows - start_row;
    }
    char fixed_text[8] = "";
    SF_macro_use("___cimport_favorstrfixed", fixed_text, sizeof(fixed_text) - 1);

    /* Merge the per-chunk statistics gathered by the parse workers. Every data
     * row contributes: numeric extrema determine integer storage, and string
     * widths reflect decoded values rather than raw quoted fields. */
    for (int col = 0; col < ctx->num_columns; col++) {
        CImportColumnInfo *colinfo = &ctx->columns[col];
        int force = ctx->col_force[col];
        bool is_string = false, has_quotes = false;
        int decoded_len = 0, max_len = 0;
        uint64_t textbytes = 0;

        colinfo->is_integer = true;
        colinfo->min_value = DBL_MAX;
        colinfo->max_value = -DBL_MAX;
        colinfo->num_subtype = CIMPORT_NUM_DOUBLE;
        for (int c = 0; c < ctx->num_chunks; c++) {
            CImportParsedChunk *chunk = &ctx->chunks[c];
            if (!chunk->col_stats || col >= chunk->num_col_stats) continue;
            CImportColumnParseStats *stats = &chunk->col_stats[col];
            if (stats->seen_string) is_string = true;
            if (stats->has_quotes) has_quotes = true;
            if (stats->max_field_len > max_len) max_len = stats->max_field_len;
            if (stats->max_decoded_len > decoded_len) decoded_len = stats->max_decoded_len;
            textbytes += stats->decoded_bytes;
            if (stats->num_nonint) colinfo->is_integer = false;
            if (stats->num_min < colinfo->min_value) colinfo->min_value = stats->num_min;
            if (stats->num_max > colinfo->max_value) colinfo->max_value = stats->num_max;
        }

        colinfo->column_has_quotes = has_quotes;
        if (force == CIMPORT_FORCE_STRING || (force != CIMPORT_FORCE_NUMERIC && is_string)) {
            colinfo->type = CIMPORT_COL_STRING;
            if (max_len >= INT_MAX) { atomic_store(&ctx->error_code, 109); return; }
            int width = colinfo->inference_width > 1 ? colinfo->inference_width : 1;
            if (decoded_len > width) width = decoded_len;
            size_t average_count = total_data_rows;
            if (ctx->has_header && ctx->num_chunks && ctx->chunks[0].num_rows) {
                const CImportParsedRow *header = ctx->chunks[0].rows[0];
                if (col < header->num_fields) {
                    CImportFieldRef field = cimport_row_field(header, col);
                    textbytes += (uint64_t)cimport_extract_field_fast(ctx->file_data, &field, NULL, INT_MAX,
                                                                      cimport_extract_quote(ctx));
                    average_count++;
                }
            }
            colinfo->max_strlen = width;
            /* Native includes the heading in its integer mean and uses strL when
               the maximum width exceeds that mean by more than 74 bytes. */
            colinfo->use_strl = width > 2045 || (strcmp(fixed_text, "1") && total_data_rows &&
                width - (int)(textbytes / average_count) > 74);
            continue;
        }

        colinfo->type = CIMPORT_COL_NUMERIC;
        if (!colinfo->is_integer && ctx->numeric_type_mode == CIMPORT_NUMTYPE_FLOAT) {
            colinfo->num_subtype = CIMPORT_NUM_FLOAT;
        } else if (!colinfo->is_integer && ctx->numeric_type_mode == CIMPORT_NUMTYPE_DOUBLE) {
            colinfo->num_subtype = CIMPORT_NUM_DOUBLE;
        } else if (colinfo->min_value == DBL_MAX) {
            colinfo->num_subtype = force == CIMPORT_FORCE_NUMERIC ? CIMPORT_NUM_DOUBLE : CIMPORT_NUM_BYTE;
        } else {
            colinfo->num_subtype = cimport_determine_numeric_subtype(colinfo, false);
        }
    }
}

/* ============================================================================
 * Context Management
 * ============================================================================ */

static void cimport_free_context(CImportContext *ctx) {
    if (!ctx) return;

    if (ctx->chunks) {
        for (int c = 0; c < ctx->num_chunks; c++) {
            CImportParsedChunk *chunk = &ctx->chunks[c];
            free(chunk->rows);
            free(chunk->col_stats);
            ctools_arena_free(&chunk->arena);
        }
        free(ctx->chunks);
    }

    if (ctx->col_cache) {
        for (int i = 0; i < ctx->num_columns; i++) {
            free(ctx->col_cache[i].numeric_data);
            free(ctx->col_cache[i].string_data);
            ctools_arena_free(&ctx->col_cache[i].string_arena);
        }
        free(ctx->col_cache);
    }
    free(ctx->columns);
    free(ctx->filename);
    free(ctx->force_numeric_cols);
    free(ctx->force_string_cols);
    free(ctx->col_force);
    free(ctx->converted_data);  /* Free encoding-converted data if any */
    cimport_munmap_file(ctx);
    pthread_mutex_destroy(&ctx->warning_mutex);
    free(ctx);
}

static void cimport_clear_cached_context(void) {
    if (g_cimport_ctx) {
        cimport_free_context(g_cimport_ctx);
        g_cimport_ctx = NULL;
    }
}

/* Forward declarations */
static CImportContext *cimport_parse_csv(const char *filename, char delimiter, bool has_header, int skip_rows, bool verbose, CImportBindQuotesMode bindquotes,
                                          CImportNumericTypeMode numeric_type_mode, char decimal_sep, char group_sep,
                                          CImportEmptyLinesMode emptylines_mode, int max_quoted_rows,
                                          CImportEncoding requested_encoding,
                                          int *force_numeric_cols, int num_force_numeric,
                                          int *force_string_cols, int num_force_string);
static void cimport_build_column_cache(CImportContext *ctx);

#include "cimport_delimiters.inc"

/* ============================================================================
 * Parallel Chunk Parsing
 * ============================================================================ */

typedef struct {
    CImportContext *ctx;
    int chunk_id;
    const char *start;
    const char *end;
    CImportParsedChunk *chunk;
} CImportChunkParseTask;

/* Loose and nobind modes end a row at every newline (quotes never bind across
 * lines), so any line start is a safe boundary: use the last one at or before
 * target, or the next one when the data start has no newline before target. */
static size_t cimport_loose_boundary(const char *data, size_t file_size,
                                     size_t data_start, size_t target) {
    if (target >= file_size) return file_size;
    const char *start = data + data_start;
    const char *ptr = data + target;
    while (ptr > start && ptr[-1] != '\n') ptr--;
    if (ptr > start) return (size_t)(ptr - data);
    const char *nl = memchr(data + target, '\n', file_size - target);
    return nl ? (size_t)(nl + 1 - data) : file_size;
}

/* Parity of the number of quote bytes in [p, end). */
static int cimport_quote_parity(const char *p, const char *end, char quote) {
    return (int)(cimport_count_byte(p, end, quote) & 1);
}

typedef struct {
    const char *start;
    const char *end;
    char quote;
    int parity;
} CImportParityTask;

static void *cimport_parity_task(void *arg) {
    CImportParityTask *task = (CImportParityTask *)arg;
    task->parity = cimport_quote_parity(task->start, task->end, task->quote);
    return NULL;
}

/* Fill boundaries[1..num_chunks-1] with row starts near data_start + i*chunk_size.
 * Strict mode may only split at a newline outside quotes. The parser toggles its
 * quote state on every quote byte, so the state at each target is the parity
 * of all quote bytes since data_start (a row start); the segments are counted
 * in parallel. Boundaries are made monotonic, so a quoted field that runs past
 * the next target leaves an empty chunk instead of overlapping chunks. */
static void cimport_chunk_boundaries(CImportContext *ctx, size_t data_start,
                                     size_t chunk_size, int num_chunks,
                                     size_t *boundaries) {
    const char *data = ctx->file_data;
    size_t file_size = ctx->file_size;
    boundaries[0] = data_start;
    boundaries[num_chunks] = file_size;

    if (ctx->bindquotes == CIMPORT_BINDQUOTES_STRICT) {
        CImportParityTask ptasks[CIMPORT_MAX_THREADS];
        int nseg = num_chunks - 1;
        for (int i = 0; i < nseg; i++) {
            ptasks[i].start = data + data_start + (size_t)i * chunk_size;
            ptasks[i].end = data + data_start + (size_t)(i + 1) * chunk_size;
            ptasks[i].quote = ctx->quote_char;
            ptasks[i].parity = 0;
        }
        ctools_persistent_pool *pool = ctools_get_global_pool();
        int pool_ok = 0;
        if (pool != NULL && nseg > 1 &&
            ctools_persistent_pool_submit_batch(pool, cimport_parity_task, ptasks, nseg,
                                                sizeof(CImportParityTask)) == 0) {
            ctools_persistent_pool_wait(pool);
            pool_ok = 1;
        }
        if (!pool_ok) {
            for (int i = 0; i < nseg; i++) cimport_parity_task(&ptasks[i]);
        }

        bool in_quotes = false;
        for (int i = 1; i < num_chunks; i++) {
            in_quotes ^= (ptasks[i - 1].parity != 0);
            size_t target = data_start + (size_t)i * chunk_size;
            /* The previous search found no row start in [its target, its
             * boundary), so if that boundary is past this target it is also
             * this target's next row start. Skipping the rescan keeps a file
             * whose quote never closes from being scanned once per chunk. */
            if (boundaries[i - 1] > target) {
                boundaries[i] = boundaries[i - 1];
                continue;
            }
            const char *row = cimport_find_next_row_strict_from(data + target, data + file_size,
                                                                ctx->quote_char, in_quotes);
            boundaries[i] = (size_t)(row - data);
        }
    } else {
        for (int i = 1; i < num_chunks; i++) {
            boundaries[i] = cimport_loose_boundary(data, file_size, data_start,
                                                   data_start + (size_t)i * chunk_size);
        }
    }

    for (int i = 1; i <= num_chunks; i++) {
        if (boundaries[i] < boundaries[i - 1]) boundaries[i] = boundaries[i - 1];
    }
}

static void *cimport_parse_chunk_parallel(void *arg) {
    CImportChunkParseTask *task = (CImportChunkParseTask *)arg;
    CImportContext *ctx = task->ctx;
    CImportParsedChunk *chunk = task->chunk;

    ctools_arena_init(&chunk->arena, 0);
    chunk->capacity = 16384;
    chunk->rows = ctools_safe_malloc2(chunk->capacity, sizeof(CImportParsedRow *));
    chunk->num_rows = 0;

    chunk->num_col_stats = ctx->num_columns;
    chunk->col_stats = ctools_safe_calloc2(ctx->num_columns, sizeof(CImportColumnParseStats));

    if (!chunk->rows || !chunk->col_stats) {
        atomic_store(&ctx->error_code, 1);
        return NULL;
    }
    cimport_init_col_stats(chunk->col_stats, ctx->num_columns);

    CImportFieldRef field_buf[CTOOLS_MAX_COLUMNS];
    const char *ptr = task->start;
    const char *end = task->end;

    while (ptr < end && atomic_load(&ctx->error_code) == 0) {
        const char *next_row;
        bool had_unmatched;
        int num_fields = cimport_parse_row_fast(ptr, end, ctx->delimiter, ctx->quote_char,
                                                 field_buf, CTOOLS_MAX_COLUMNS, ctx->file_data,
                                                 ctx->bindquotes, &next_row, &had_unmatched);

        /* Handle empty lines (including CRLF blank lines) based on emptylines mode */
        if (ctx->emptylines_mode == CIMPORT_EMPTYLINES_SKIP &&
            cimport_row_is_empty(ctx->file_data, field_buf, num_fields)) {
            ptr = next_row;
            continue;
        }

        bool is_header_row = (ctx->has_header && task->chunk_id == 0 && chunk->num_rows == 0);

        /* Check for unmatched quotes in LOOSE mode (derived from parser state) */
        if (ctx->bindquotes == CIMPORT_BINDQUOTES_LOOSE && !is_header_row) {
            if (had_unmatched) {
                size_t approx_row = chunk->num_rows + 1;
                if (task->chunk_id > 0) {
                    approx_row += task->chunk_id * (ctx->file_size / ctx->num_chunks / 50);
                }
                cimport_record_unmatched_quote(ctx, approx_row);
            }
        }


        if (chunk->num_rows >= chunk->capacity) {
            if (chunk->capacity > SIZE_MAX / 2) {
                atomic_store(&ctx->error_code, 1);
                break;
            }
            size_t new_capacity = chunk->capacity * 2;
            if (new_capacity > SIZE_MAX / sizeof(CImportParsedRow *)) {
                atomic_store(&ctx->error_code, 1);
                break;
            }
            CImportParsedRow **new_rows = realloc(chunk->rows, sizeof(CImportParsedRow *) * new_capacity);
            if (!new_rows) {
                atomic_store(&ctx->error_code, 1);
                break;
            }
            chunk->rows = new_rows;
            chunk->capacity = new_capacity;
        }

        CImportParsedRow *row = ctools_arena_alloc(&chunk->arena, cimport_row_size(num_fields));
        if (!row || !cimport_fill_row(row, field_buf, num_fields)) {
            atomic_store(&ctx->error_code, 1);
            break;
        }
        uint64_t *values = cimport_row_values(row);
        for (int f = 0; f < num_fields; f++) {
            values[f] = !is_header_row && f < chunk->num_col_stats
                ? cimport_accumulate_field(ctx, &chunk->col_stats[f], &field_buf[f], ctx->col_force[f])
                : CIMPORT_NO_VALUE;
        }

        chunk->rows[chunk->num_rows++] = row;

        if (num_fields > chunk->max_fields_in_chunk) {
            chunk->max_fields_in_chunk = num_fields;
        }

        ptr = next_row;
    }

    return NULL;
}

/* ============================================================================
 * Delimiter Auto-Detection (matches Stata's import delimited behavior)
 * ============================================================================ */

/* Native import delimited (Stata 17 and later) counts tab, comma, colon,
 * semicolon and pipe characters outside quotes in the first 50 non-empty
 * lines, header included, and uses the most frequent one. Ties go to comma,
 * then tab, pipe, colon and semicolon; a file with none of them is comma
 * delimited. So `a;b` over rows like `1,5;x` splits on semicolons although
 * each data row has as many commas. Lines end at \n, \r\n or \r, except
 * inside quotes under bindquotes(strict); nobind ignores quotes. first and
 * last (physical lines, 0 when open) limit the lines as native does for
 * varnames(#) and for varnames(nonames) with rowrange(). */
static char cimport_auto_detect_delimiter(const char *data, size_t size,
                                          CImportBindQuotesMode bindquotes,
                                          size_t first, size_t last) {
    static const char candidates[] = {';', ':', '|', '\t', ','};  /* rising tie priority */
    size_t counts[5] = {0};
    const bool quotes = bindquotes != CIMPORT_BINDQUOTES_NOBIND;
    const char *p = data, *end = data + size;
    size_t lines = 0;

    for (int read = 0; read < 50 && !(last && lines >= last); read++) {
        /* Next non-empty line ending on or after the first line */
        const char *start = NULL, *stop = NULL;
        while (p < end && !start) {
            const char *line = p, *line_end = end;
            bool in_quotes = false;
            for (;;) {
                if (p == end) { lines++; break; }
                char c = *p;
                if (c == '\n' || c == '\r') {
                    const char *eol = p;
                    p += (c == '\r' && p + 1 < end && p[1] == '\n') ? 2 : 1;
                    lines++;
                    if (in_quotes && bindquotes == CIMPORT_BINDQUOTES_STRICT) continue;
                    line_end = eol;
                    break;
                }
                if (c == '"' && quotes) in_quotes = !in_quotes;
                p++;
            }
            if (line_end > line && (!first || lines >= first)) {
                start = line;
                stop = line_end;
            }
        }
        if (!start) break;

        bool in_quotes = false;
        for (const char *q = start; q < stop; q++) {
            if (*q == '"' && quotes) in_quotes = !in_quotes;
            if (in_quotes) continue;
            for (int i = 0; i < 5; i++) counts[i] += *q == candidates[i];
        }
    }

    int best = 0;
    for (int i = 1; i < 5; i++) {
        if (counts[i] >= counts[best]) best = i;
    }
    return candidates[best];
}

/* ============================================================================
 * Main Parse Function
 * ============================================================================ */

static CImportContext *cimport_parse_csv(const char *filename, char delimiter, bool has_header, int skip_rows, bool verbose, CImportBindQuotesMode bindquotes,
                                          CImportNumericTypeMode numeric_type_mode, char decimal_sep, char group_sep,
                                          CImportEmptyLinesMode emptylines_mode, int max_quoted_rows,
                                          CImportEncoding requested_encoding,
                                          int *force_numeric_cols, int num_force_numeric,
                                          int *force_string_cols, int num_force_string) {
    CImportContext *ctx = calloc(1, sizeof(CImportContext));
    if (!ctx) return NULL;

    double t_start, t_end;
    char msg[512];

    g_cimport_parse_status = 601;

    char locale[128]="",locale_decimal[8]="",locale_group[8]="";
    SF_macro_use("___cimport_parselocale",locale,sizeof(locale)-1);
    SF_macro_use("___cimport_decimalseparator",locale_decimal,sizeof(locale_decimal)-1);
    SF_macro_use("___cimport_groupseparator",locale_group,sizeof(locale_group)-1);
    int locale_rc=cimport_configure_locale(locale,locale_decimal,locale_group);
    if(locale_rc) {g_cimport_parse_status=locale_rc;free(ctx);free(force_numeric_cols);free(force_string_cols);return NULL;}

    ctx->delimiter = delimiter;
    ctx->quote_char = '"';
    ctx->has_header = has_header;
    ctx->skip_rows = skip_rows;
    ctx->bindquotes = bindquotes;
    char quote_option[16] = "";
    SF_macro_use("___cimport_stripquotes", quote_option, sizeof(quote_option)-1);
    ctx->retain_quotes = strcmp(quote_option, "no") == 0;
    ctx->strip_quotes_all = strcmp(quote_option, "yes") == 0;
    ctx->verbose = verbose;
    atomic_init(&ctx->error_code, 0);

    /* Initialize warning tracking early - must happen before any path that
     * calls cimport_free_context(), which always destroys warning_mutex. */
    ctx->num_unmatched_quote_warnings = 0;
    pthread_mutex_init(&ctx->warning_mutex, NULL);

    /* Take ownership of force lists immediately so cimport_free_context always
     * handles them, preventing double-free if we fail and the caller also frees. */
    ctx->force_numeric_cols = force_numeric_cols;
    ctx->num_force_numeric = num_force_numeric;
    ctx->force_string_cols = force_string_cols;
    ctx->num_force_string = num_force_string;

    ctx->filename = strdup(filename);
    if (ctx->filename == NULL) {
        cimport_free_context(ctx);
        return NULL;
    }

    /* Initialize new options */
    ctx->numeric_type_mode = numeric_type_mode;
    ctx->decimal_separator = decimal_sep;
    ctx->group_separator = group_sep;
    ctx->emptylines_mode = emptylines_mode;
    ctx->max_quoted_rows = max_quoted_rows;

    /* Initialize encoding */
    ctx->requested_encoding = requested_encoding;
    ctx->encoding = CIMPORT_ENC_UTF8;
    ctx->converted_data = NULL;
    ctx->converted_size = 0;
    ctx->time_encoding = 0;

    t_start = ctools_timer_ms();
    if (cimport_mmap_file(ctx, filename) != 0) {
        cimport_display_error(ctx->error_message);
        cimport_free_context(ctx);
        return NULL;
    }
    t_end = ctools_timer_ms();
    ctx->time_mmap = t_end - t_start;

    /* Recognition and conversion execute in C; keep the source charset for r(). */
    t_start = ctools_timer_ms();
    char encoding_name[256] = "";
    SF_macro_use("___cimport_encoding_name", encoding_name,
                 sizeof(encoding_name) - 1);
    bool explicit_encoding = *encoding_name != 0;
    const char *source_encoding =
        explicit_encoding ? encoding_name
        : requested_encoding != CIMPORT_ENC_UNKNOWN
            ? cimport_encoding_name(requested_encoding)
            : cimport_detect_charset(ctx->file_data, ctx->file_size);
    snprintf(ctx->reported_encoding, sizeof(ctx->reported_encoding), "%s",
             !strcmp(source_encoding, "latin1") ? "ISO-8859-1"
                                                : source_encoding);
    ctx->encoding = cimport_parse_encoding_name(source_encoding);
    /* ASCII can retain the original mapping after statistical recognition. */
    bool ascii_fast_path =
        !explicit_encoding &&
        (ctx->encoding == CIMPORT_ENC_UTF8 ||
         ctx->encoding == CIMPORT_ENC_ASCII ||
         ctx->encoding == CIMPORT_ENC_ISO_8859_1 ||
         ctx->encoding == CIMPORT_ENC_ISO_8859_15 ||
         ctx->encoding == CIMPORT_ENC_WINDOWS_1252 ||
         ctx->encoding == CIMPORT_ENC_MACROMAN) &&
        cimport_is_ascii_simd(ctx->file_data, ctx->file_size);
    if (!ascii_fast_path) {
      char *converted = NULL;
      size_t converted_size = 0;
      int rc = cimport_convert_named_to_utf8(ctx->file_data, ctx->file_size,
                                             source_encoding, &converted,
                                             &converted_size);
      if (rc) {
        g_cimport_parse_status = rc;
        cimport_free_context(ctx);
        return NULL;
      }
      ctx->converted_data = converted;
      ctx->converted_size = converted_size;
      ctx->file_data = converted;
      ctx->file_size = converted_size;
    }
    t_end = ctools_timer_ms();
    ctx->time_encoding = t_end - t_start;
    /* Empty input still validates an explicit charset and reports its encoding.
     */
    if (ctx->file_size == 0) {
      if (!ctx->delimiter)
        ctx->delimiter = ',';
      if (verbose)
        cimport_display_msg("(empty file - creating empty dataset)\n");
      return ctx;
    }

    int normalize_rc = cimport_normalize_delimiters(ctx);
    if (!normalize_rc)
      normalize_rc = cimport_check_quoted_rows(ctx);
    if (normalize_rc) {
      g_cimport_parse_status = normalize_rc;
      cimport_display_error("cimport: could not normalize delimiters\n");
      cimport_free_context(ctx);
      return NULL;
    }

    /* Auto-detect delimiter if not specified (delimiter == '\0') */
    if (ctx->delimiter == '\0') {
        /* Native scans from the varnames(#) row, or over the rowrange() lines
         * with varnames(nonames); an automatic header scans from line 1. */
        size_t first = 0, last = 0;
        char autoheader[8] = "";
        SF_macro_use("___cimport_autoheader", autoheader, sizeof(autoheader) - 1);
        if (!ctx->has_header || strcmp(autoheader, "1") != 0) {
            cimport_row_range(&first, &last);
            if (ctx->has_header) first = (size_t)ctx->skip_rows + 1;
        }
        ctx->delimiter = cimport_auto_detect_delimiter(
            ctx->file_data, ctx->file_size, ctx->bindquotes, first, last);
        if (verbose) {
            char delim_name[16];
            if (ctx->delimiter == '\t') snprintf(delim_name, sizeof(delim_name), "tab");
            else if (ctx->delimiter == ' ') snprintf(delim_name, sizeof(delim_name), "space");
            else snprintf(delim_name, sizeof(delim_name), "'%c'", ctx->delimiter);
            snprintf(msg, sizeof(msg), "(delimiter not specified, auto-detected: %s)\n", delim_name);
            cimport_display_msg(msg);
        }
    }

    if (cimport_parse_header(ctx) != 0) {
        cimport_display_error("Failed to parse header");
        cimport_free_context(ctx);
        return NULL;
    }
    /* Only empty lines: no variables, as for an empty file. */
    if (ctx->num_columns == 0) return ctx;
    if (cimport_resolve_col_force(ctx, 0, ctx->num_columns) != 0) {
        cimport_free_context(ctx);
        return NULL;
    }
    /* The single-pass number parser applies only to the default grammar. */
    ctx->plain_numbers = ctx->decimal_separator == '.' && ctx->group_separator == '\0' &&
                         !cimport_locale_active();

#ifdef CIMPORT_WINDOWS
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    int cpu_count = (int)sysinfo.dwNumberOfProcessors;
#else
    int cpu_count = (int)sysconf(_SC_NPROCESSORS_ONLN);
#endif
    if (cpu_count <= 0) cpu_count = 4;
    if (cpu_count > CIMPORT_MAX_THREADS) cpu_count = CIMPORT_MAX_THREADS;
    ctx->num_threads = cpu_count;

    t_start = ctools_timer_ms();

    size_t min_chunk_size = 1 * 1024 * 1024;
    int num_chunks = ctx->num_threads;

    if (ctx->file_size < min_chunk_size * 2) {
        num_chunks = 1;
    }

    ctx->num_chunks = num_chunks;
    ctx->chunks = ctools_safe_calloc2(ctx->num_chunks, sizeof(CImportParsedChunk));
    if (!ctx->chunks) {
        cimport_free_context(ctx);
        return NULL;
    }

    if (num_chunks == 1) {
        CImportParsedChunk *chunk = &ctx->chunks[0];
        ctools_arena_init(&chunk->arena, 0);
        chunk->capacity = 65536;
        chunk->rows = ctools_safe_malloc2(chunk->capacity, sizeof(CImportParsedRow *));
        chunk->num_rows = 0;

        chunk->num_col_stats = ctx->num_columns;
        chunk->col_stats = ctools_safe_calloc2(ctx->num_columns, sizeof(CImportColumnParseStats));

        if (!chunk->rows || !chunk->col_stats) {
            cimport_free_context(ctx);
            return NULL;
        }
        cimport_init_col_stats(chunk->col_stats, ctx->num_columns);

        CImportFieldRef field_buf[CTOOLS_MAX_COLUMNS];
        const char *ptr = ctx->file_data;
        const char *end = ctx->file_data + ctx->file_size;

        /* Skip rows if requested (for varnames(N) where N > 1) */
        if (ctx->skip_rows > 0) {
            ptr = cimport_skip_rows(ptr, end, ctx->skip_rows, ctx->quote_char, ctx->bindquotes);
        }

        while (ptr < end) {
            const char *next_row;
            bool had_unmatched;
            int num_fields = cimport_parse_row_fast(ptr, end, ctx->delimiter, ctx->quote_char,
                                                     field_buf, CTOOLS_MAX_COLUMNS, ctx->file_data,
                                                     ctx->bindquotes, &next_row, &had_unmatched);

            /* Handle empty lines (including CRLF blank lines) based on emptylines mode */
            if (ctx->emptylines_mode == CIMPORT_EMPTYLINES_SKIP &&
                cimport_row_is_empty(ctx->file_data, field_buf, num_fields)) {
                ptr = next_row;
                continue;
            }

            bool is_header_row = (ctx->has_header && chunk->num_rows == 0);

            /* Check for unmatched quotes in LOOSE mode (derived from parser state) */
            if (ctx->bindquotes == CIMPORT_BINDQUOTES_LOOSE && !is_header_row) {
                if (had_unmatched) {
                    size_t row_num = chunk->num_rows + 1;
                    cimport_record_unmatched_quote(ctx, row_num);
                }
            }

            if (chunk->num_rows >= chunk->capacity) {
                if (chunk->capacity > SIZE_MAX / 2) {
                    cimport_free_context(ctx);
                    return NULL;
                }
                size_t new_capacity = chunk->capacity * 2;
                if (new_capacity > SIZE_MAX / sizeof(CImportParsedRow *)) {
                    cimport_free_context(ctx);
                    return NULL;
                }
                CImportParsedRow **new_rows = realloc(chunk->rows, sizeof(CImportParsedRow *) * new_capacity);
                if (!new_rows) {
                    cimport_free_context(ctx);
                    return NULL;
                }
                chunk->rows = new_rows;
                chunk->capacity = new_capacity;
            }

            CImportParsedRow *row = ctools_arena_alloc(&chunk->arena, cimport_row_size(num_fields));
            if (!row || !cimport_fill_row(row, field_buf, num_fields)) {
                cimport_free_context(ctx);
                return NULL;
            }
            uint64_t *values = cimport_row_values(row);
            for (int f = 0; f < num_fields; f++) {
                values[f] = !is_header_row && f < chunk->num_col_stats
                    ? cimport_accumulate_field(ctx, &chunk->col_stats[f], &field_buf[f], ctx->col_force[f])
                    : CIMPORT_NO_VALUE;
            }
            chunk->rows[chunk->num_rows++] = row;

            if (num_fields > chunk->max_fields_in_chunk) {
                chunk->max_fields_in_chunk = num_fields;
            }

            ptr = next_row;
        }
    } else {
        /* Skip rows if requested (for varnames(N) where N > 1) */
        const char *data_start = ctx->file_data;
        size_t data_size = ctx->file_size;
        if (ctx->skip_rows > 0) {
            data_start = cimport_skip_rows(ctx->file_data, ctx->file_data + ctx->file_size,
                                           ctx->skip_rows, ctx->quote_char, ctx->bindquotes);
            data_size = ctx->file_size - (data_start - ctx->file_data);
        }

        size_t chunk_size = data_size / num_chunks;
        size_t start_offset = data_start - ctx->file_data;

        size_t *boundaries = ctools_safe_malloc2((size_t)num_chunks + 1, sizeof(size_t));
        if (!boundaries) {
            cimport_free_context(ctx);
            return NULL;
        }

        cimport_chunk_boundaries(ctx, start_offset, chunk_size, num_chunks, boundaries);

        CImportChunkParseTask tasks[CIMPORT_MAX_THREADS];

        for (int i = 0; i < num_chunks; i++) {
            tasks[i].ctx = ctx;
            tasks[i].chunk_id = i;
            tasks[i].start = ctx->file_data + boundaries[i];
            tasks[i].end = ctx->file_data + boundaries[i + 1];
            tasks[i].chunk = &ctx->chunks[i];
        }

        ctools_persistent_pool *pool = ctools_get_global_pool();
        int pool_ok = 0;
        if (pool != NULL) {
            if (ctools_persistent_pool_submit_batch(pool, cimport_parse_chunk_parallel,
                                                     tasks, num_chunks,
                                                     sizeof(CImportChunkParseTask)) == 0) {
                ctools_persistent_pool_wait(pool);
                pool_ok = 1;
            }
        }
        if (!pool_ok) {
            /* Fallback: raw pthreads */
            pthread_t threads[CIMPORT_MAX_THREADS];
            int threads_created = 0;
            for (int i = 0; i < num_chunks; i++) {
                if (pthread_create(&threads[threads_created], NULL, cimport_parse_chunk_parallel, &tasks[i]) == 0) {
                    threads_created++;
                } else {
                    cimport_parse_chunk_parallel(&tasks[i]);
                }
            }
            for (int i = 0; i < threads_created; i++) {
                pthread_join(threads[i], NULL);
            }
        }

        free(boundaries);

        if (atomic_load(&ctx->error_code) != 0) {
            cimport_free_context(ctx);
            return NULL;
        }
    }

    t_end = ctools_timer_ms();
    ctx->time_parse = t_end - t_start;

    ctx->total_rows = 0;
    for (int c = 0; c < ctx->num_chunks; c++) {
        ctx->total_rows += ctx->chunks[c].num_rows;
    }
    if (ctx->has_header && ctx->total_rows > 0) {
        ctx->total_rows--;
    }

    ctx->max_fields_seen = ctx->num_columns;
    for (int c = 0; c < ctx->num_chunks; c++) {
        if (ctx->chunks[c].max_fields_in_chunk > ctx->max_fields_seen) {
            ctx->max_fields_seen = ctx->chunks[c].max_fields_in_chunk;
        }
    }

    if (ctx->max_fields_seen > ctx->num_columns) {
        int old_num = ctx->num_columns;
        int new_num = ctx->max_fields_seen;

        CImportColumnInfo *new_columns = realloc(ctx->columns, new_num * sizeof(CImportColumnInfo));
        if (!new_columns) {
            cimport_free_context(ctx);
            return NULL;
        }
        ctx->columns = new_columns;

        for (int i = old_num; i < new_num; i++) {
            memset(&ctx->columns[i], 0, sizeof(CImportColumnInfo));
            snprintf(ctx->columns[i].name, CTOOLS_MAX_VARNAME_LEN, "v%d", i + 1);
            ctx->columns[i].type = CIMPORT_COL_UNKNOWN;
            ctx->columns[i].max_strlen = 0;
        }

        ctx->num_columns = new_num;
        if (cimport_resolve_col_force(ctx, old_num, new_num) != 0) {
            cimport_free_context(ctx);
            return NULL;
        }

        /* Expand chunk col_stats arrays and collect type stats for expanded columns.
         * Since chunks were parsed with the old column count, stats for expanded
         * columns were not collected. Re-scan each chunk's rows for the new columns. */
        for (int c = 0; c < ctx->num_chunks; c++) {
            CImportParsedChunk *chunk = &ctx->chunks[c];
            if (chunk->num_col_stats < new_num) {
                CImportColumnParseStats *new_stats = realloc(chunk->col_stats,
                    new_num * sizeof(CImportColumnParseStats));
                if (!new_stats) {
                    cimport_free_context(ctx);
                    return NULL;
                }
                cimport_init_col_stats(new_stats + chunk->num_col_stats, new_num - chunk->num_col_stats);
                chunk->col_stats = new_stats;

                /* Re-scan rows in this chunk for expanded column stats */
                size_t start_row = (c == 0 && ctx->has_header) ? 1 : 0;
                for (size_t r = start_row; r < chunk->num_rows; r++) {
                    CImportParsedRow *row = chunk->rows[r];
                    uint64_t *values = cimport_row_values(row);
                    CImportFieldRef field = {row->offset, 0};
                    for (int f = 0; f < row->num_fields && f < new_num; f++) {
                        field.length = row->lengths[f];
                        if (f >= chunk->num_col_stats)
                            values[f] = cimport_accumulate_field(ctx, &new_stats[f], &field, ctx->col_force[f]);
                        field.offset += CIMPORT_FIELD_LENGTH(field) + 1;
                    }
                }

                chunk->num_col_stats = new_num;
            }
        }
    }

    t_start = ctools_timer_ms();
    cimport_select_rows(ctx);
    cimport_infer_column_types(ctx);
    if (atomic_load(&ctx->error_code)) { cimport_free_context(ctx); return NULL; }
    t_end = ctools_timer_ms();
    ctx->time_type_infer = t_end - t_start;

    ctx->is_loaded = true;
    ctx->cache_ready = false;

    return ctx;
}

/* ============================================================================
 * Column Cache Builder
 * ============================================================================ */

typedef struct {
    CImportContext *ctx;
    int col_start;
    int col_end;
    int thread_id;
} CImportCacheBuildTask;

static void *cimport_build_cache_worker(void *arg) {
    CImportCacheBuildTask *task = (CImportCacheBuildTask *)arg;
    CImportContext *ctx = task->ctx;
    char *field_buf = NULL;
    size_t field_capacity = 0;

    for (int col_idx = task->col_start; col_idx < task->col_end; col_idx++) {
        CImportColumnInfo *col = &ctx->columns[col_idx];
        CImportColumnCache *cache = &ctx->col_cache[col_idx];
        size_t row_idx = 0;

        if (col->type == CIMPORT_COL_STRING) {
            if (cache->string_data == NULL) {
                cache->count = 0;
                continue;
            }

            ctools_arena_init(&cache->string_arena, 0);
            size_t need = (size_t)col->max_strlen + 1;
            if (need > INT_MAX) {atomic_store(&ctx->error_code,109);free(field_buf);return (void*)1;}
            if (need > field_capacity) {
                char *next=realloc(field_buf,need);
                if (!next) {atomic_store(&ctx->error_code,920);free(field_buf);return (void*)1;}
                field_buf=next;field_capacity=need;
            }

            for (int c = 0; c < ctx->num_chunks; c++) {
                CImportParsedChunk *chunk = &ctx->chunks[c];
                size_t start_row = (c == 0 && ctx->has_header) ? 1 : 0;

                /* Per-chunk quote check: skip per-field SIMD scan if
                 * this chunk has no quotes for this column */
                bool chunk_has_quotes = col->column_has_quotes &&
                    chunk->col_stats && col_idx < chunk->num_col_stats &&
                    chunk->col_stats[col_idx].has_quotes;

                for (size_t r = start_row; r < chunk->num_rows; r++) {
                    CImportParsedRow *row = chunk->rows[r];
                    char *str;

                    if (col_idx < row->num_fields) {
                        CImportFieldRef field_ref = cimport_row_field(row, col_idx);
                        CImportFieldRef *field = &field_ref;
                        const char *src = ctx->file_data + field->offset;
                        bool field_has_quote = (field->length & CIMPORT_FIELD_QUOTED_FLAG) != 0;

                        if (ctx->retain_quotes || !chunk_has_quotes || !field_has_quote) {
                            /* Fast path: no quotes in this field.
                             * Copy directly from mmap to arena (single copy). */
                            int src_len = CIMPORT_FIELD_LENGTH(*field);
                            while (src_len > 0 && (src[src_len-1] == '\r' || src[src_len-1] == '\n')) {
                                src_len--;
                            }
                            int copy_len = src_len;
                            str = ctools_arena_alloc(&cache->string_arena, copy_len + 1);
                            if (str) {
                                memcpy(str, src, copy_len);
                                str[copy_len] = '\0';
                            } else {
                                atomic_store(&ctx->error_code,920);free(field_buf);return (void*)1;
                            }
                        } else {
                            /* Field has quotes — extract with quote handling */
                            int len = cimport_extract_field_fast(ctx->file_data, field, field_buf, (int)field_capacity, cimport_extract_quote(ctx));
                            field_buf[len] = '\0';
                            str = ctools_arena_alloc(&cache->string_arena, len + 1);
                            if (str) {
                                memcpy(str, field_buf, len + 1);
                            } else {
                                atomic_store(&ctx->error_code,920);free(field_buf);return (void*)1;
                            }
                        }
                    } else {
                        str = ctools_arena_alloc(&cache->string_arena, 1);
                        if (str) {
                            str[0] = '\0';
                        } else {
                            atomic_store(&ctx->error_code,920);free(field_buf);return (void*)1;
                        }
                    }

                    cache->string_data[row_idx++] = str;
                }
            }
        } else {
            if (cache->numeric_data == NULL) {
                cache->count = 0;
                continue;
            }

            double missing = SV_missval;

            for (int c = 0; c < ctx->num_chunks; c++) {
                CImportParsedChunk *chunk = &ctx->chunks[c];
                size_t start_row = (c == 0 && ctx->has_header) ? 1 : 0;

                for (size_t r = start_row; r < chunk->num_rows; r++) {
                    CImportParsedRow *row = chunk->rows[r];
                    double val;

                    if (col_idx < row->num_fields) {
                        CImportFieldRef field_ref = cimport_row_field(row, col_idx);
                        CImportFieldRef *field = &field_ref;
                        if (!cimport_context_number(ctx,field,&val)) {
                            val = missing;
                        }
                    } else {
                        val = missing;
                    }

                    cache->numeric_data[row_idx++] = val;
                }
            }
        }

        cache->count = row_idx;
    }

    free(field_buf);
    return NULL;
}

static void cimport_build_column_cache(CImportContext *ctx) {
    if (ctx->cache_ready) return;

    /* Handle empty file case - no columns or rows to cache */
    if (ctx->num_columns == 0 || ctx->total_rows == 0) {
        ctx->col_cache = NULL;
        ctx->cache_ready = true;
        ctx->time_cache = 0;
        return;
    }

    double t_start = ctools_timer_ms();

    ctx->col_cache = ctools_safe_calloc2(ctx->num_columns, sizeof(CImportColumnCache));
    if (!ctx->col_cache) return;

    for (int i = 0; i < ctx->num_columns; i++) {
        CImportColumnInfo *col = &ctx->columns[i];
        CImportColumnCache *cache = &ctx->col_cache[i];

        if (col->type == CIMPORT_COL_STRING && col->use_strl) {
            cache->string_data = ctools_safe_calloc2(ctx->total_rows, sizeof(char *));
            if (!cache->string_data) { atomic_store(&ctx->error_code,920);return; }
            cache->numeric_data = NULL;
        } else {
            /* Skip numeric column cache — parsed directly during SPI store */
            cache->numeric_data = NULL;
            cache->string_data = NULL;
        }
        cache->count = 0;
    }

    int num_threads = ctx->num_threads;
    if (num_threads > ctx->num_columns) num_threads = ctx->num_columns;
    if (num_threads < 1) num_threads = 1;

    if (num_threads == 1 || ctx->num_columns < 4) {
        CImportCacheBuildTask task = {ctx, 0, ctx->num_columns, 0};
        cimport_build_cache_worker(&task);
    } else {
        CImportCacheBuildTask tasks[CIMPORT_MAX_THREADS];

        int cols_per_thread = (ctx->num_columns + num_threads - 1) / num_threads;

        for (int t = 0; t < num_threads; t++) {
            tasks[t].ctx = ctx;
            tasks[t].col_start = t * cols_per_thread;
            tasks[t].col_end = (t + 1) * cols_per_thread;
            if (tasks[t].col_end > ctx->num_columns) tasks[t].col_end = ctx->num_columns;
            tasks[t].thread_id = t;
        }

        ctools_persistent_pool *pool = ctools_get_global_pool();
        int pool_ok = 0;
        if (pool != NULL) {
            if (ctools_persistent_pool_submit_batch(pool, cimport_build_cache_worker,
                                                     tasks, num_threads,
                                                     sizeof(CImportCacheBuildTask)) == 0) {
                ctools_persistent_pool_wait(pool);
                pool_ok = 1;
            }
        }
        if (!pool_ok) {
            /* Fallback: raw pthreads */
            pthread_t threads[CIMPORT_MAX_THREADS];
            int threads_created = 0;
            for (int t = 0; t < num_threads; t++) {
                if (pthread_create(&threads[threads_created], NULL, cimport_build_cache_worker, &tasks[t]) == 0) {
                    threads_created++;
                } else {
                    cimport_build_cache_worker(&tasks[t]);
                }
            }
            for (int t = 0; t < threads_created; t++) {
                pthread_join(threads[t], NULL);
            }
        }
    }

    if (atomic_load(&ctx->error_code)) return;
    ctx->cache_ready = true;

    double t_end = ctools_timer_ms();
    ctx->time_cache = t_end - t_start;
}

/* ============================================================================
 * Column List Parsing (for numericcols/stringcols options)
 * ============================================================================ */

/* Parse column list from Stata global macro */
static int *cimport_parse_col_list(const char *macro_name, int *out_count) {
    char buf[4096];
    *out_count = 0;

    if (SF_macro_use((char *)macro_name, buf, sizeof(buf)) != 0 || strlen(buf) == 0) {
        return NULL;
    }

    if (strcmp(buf, "_all") == 0) {
        int *all = calloc(1, sizeof(int));
        if (all) *out_count = 1;
        return all;
    }
    /* Count tokens */
    int count = 0;
    char *p = buf;
    while (*p) {
        while (*p == ' ') p++;
        if (*p == '\0') break;
        count++;
        while (*p && *p != ' ') p++;
    }

    if (count == 0) return NULL;

    int *cols = ctools_safe_malloc2(count, sizeof(int));
    if (!cols) return NULL;

    /* Parse tokens - use strtol directly since ctools_safe_atoi expects null-terminated strings */
    p = buf;
    int idx = 0;
    while (*p && idx < count) {
        while (*p == ' ') p++;
        if (*p == '\0') break;
        char *endptr;
        errno = 0;
        long val = strtol(p, &endptr, 10);
        if (errno == ERANGE || endptr == p || val <= 0 || val > INT_MAX) {
            free(cols);
            return NULL;  /* Invalid column number */
        }
        cols[idx] = (int)val;
        idx++;
        p = endptr;  /* Move to character after the parsed number */
    }

    *out_count = idx;
    return cols;
}

/* ============================================================================
 * SCAN Mode
 * ============================================================================ */

static ST_retcode cimport_do_scan(const char *filename, char delimiter, bool has_header, int skip_rows, bool verbose, CImportBindQuotesMode bindquotes,
                                   CImportNumericTypeMode numeric_type_mode, char decimal_sep, char group_sep,
                                   CImportEmptyLinesMode emptylines_mode, int max_quoted_rows,
                                   CImportEncoding requested_encoding) {
    char macro_val[65536];

    cimport_clear_cached_context();

    /* Parse column type overrides BEFORE parsing CSV so inference runs once with overrides */
    int num_force_num = 0, num_force_str = 0;
    int *force_num = cimport_parse_col_list("CIMPORT_NUMCOLS", &num_force_num);
    int *force_str = cimport_parse_col_list("CIMPORT_STRCOLS", &num_force_str);

    g_cimport_ctx = cimport_parse_csv(filename, delimiter, has_header, skip_rows, verbose, bindquotes,
                                       numeric_type_mode, decimal_sep, group_sep, emptylines_mode, max_quoted_rows,
                                       requested_encoding,
                                       force_num, num_force_num, force_str, num_force_str);

    if (!g_cimport_ctx) {
        /* cimport_parse_csv took ownership of force_num/force_str and freed them
         * on failure via cimport_free_context, so do NOT free them here. */
        return g_cimport_parse_status;
    }

    CImportContext *ctx = g_cimport_ctx;

    ST_retcode report_rc =
        SF_macro_save("_cimport_encoding", ctx->reported_encoding);
    char returned_delimiters[4096] = "";
    if (!report_rc)
      report_rc = SF_macro_use("___cimport_delimiters", returned_delimiters,
                               sizeof(returned_delimiters) - 1);
    if (!*returned_delimiters) {
      returned_delimiters[0] = ctx->delimiter;
      returned_delimiters[1] = 0;
    }
    if (!report_rc)
      report_rc = SF_macro_save("_cimport_delimiters", returned_delimiters);
    if (report_rc) {
      cimport_free_context(g_cimport_ctx);
      g_cimport_ctx = NULL;
      return report_rc;
    }

    snprintf(macro_val, sizeof(macro_val), "%zu", ctx->total_rows);
    SF_macro_save("_cimport_nobs", macro_val);

    snprintf(macro_val, sizeof(macro_val), "%d", ctx->num_columns);
    SF_macro_save("_cimport_nvar", macro_val);

    char *p = macro_val;
    char *end = macro_val + sizeof(macro_val) - 1;
    for (int i = 0; i < ctx->num_columns && p < end; i++) {
        if (i > 0 && p < end) *p++ = ' ';
        size_t name_len = strlen(ctx->columns[i].name);
        size_t space_left = (size_t)(end - p);
        if (name_len > space_left) name_len = space_left;
        memcpy(p, ctx->columns[i].name, name_len);
        p += name_len;
    }
    *p = '\0';
    SF_macro_save("_cimport_varnames", macro_val);

    p = macro_val;
    for (int i = 0; i < ctx->num_columns && p < end; i++) {
        if (i > 0 && p < end) *p++ = ' ';
        if (p < end) *p++ = (ctx->columns[i].type == CIMPORT_COL_STRING) ? '1' : '0';
    }
    *p = '\0';
    SF_macro_save("_cimport_vartypes", macro_val);

    p = macro_val;
    for (int i = 0; i < ctx->num_columns && p < end; i++) {
        if (i > 0 && p < end) *p++ = ' ';
        if (p < end) {
            if (ctx->columns[i].type == CIMPORT_COL_STRING) {
                *p++ = '0';
            } else {
                *p++ = '0' + (int)ctx->columns[i].num_subtype;
            }
        }
    }
    *p = '\0';
    SF_macro_save("_cimport_numtypes", macro_val);

    p = macro_val;
    for (int i = 0; i < ctx->num_columns && p < end; i++) {
        if (i > 0 && p < end) *p++ = ' ';
        int len = ctx->columns[i].type == CIMPORT_COL_STRING ? ctx->columns[i].max_strlen : 0;
        int written = snprintf(p, (size_t)(end - p + 1), "%d", len);
        if (written > 0 && p + written <= end) p += written;
    }
    *p = '\0';
    SF_macro_save("_cimport_strlens", macro_val);
    p=macro_val;
    for(int i=0;i<ctx->num_columns && p<end;i++) {
        if(i>0 && p<end)*p++=' ';
        if(p<end)*p++=ctx->columns[i].use_strl ? '1' : '0';
    }
    *p=0;SF_macro_save("_cimport_longtypes",macro_val);

    /* Transport original headings independently, including quotes and spaces. */
    CImportFieldRef headings[CTOOLS_MAX_COLUMNS];
    int count = 0;
    if (ctx->has_header) {
        count = cimport_parse_header_row(ctx, headings, NULL, NULL);
    }
    for (int i = 0; i < ctx->num_columns; i++) {
        char key[64], value[CTOOLS_MAX_STRING_LEN + 1] = "";
        snprintf(key, sizeof(key), "_cimport_label_%d", i + 1);
        if (i < count) cimport_extract_field_fast(ctx->file_data, &headings[i], value, sizeof(value), cimport_extract_quote(ctx));
        SF_macro_save(key, value);
    }

    /* Display warnings for malformed data (like Stata) */
    cimport_display_warnings(ctx);

    /* Timing data saved to scalars for verbose table in .ado */

    return 0;
}

/* ============================================================================
 * Parallel SPI Store
 * ============================================================================ */

/* Stata keeps observations contiguous, so each task writes every column of a
 * block of rows. Column-parallel stores instead touch one cache line per cell
 * and let neighbouring columns share lines across threads. */
#define CIMPORT_STORE_TILE_ROWS 2048

/* Per-column store kinds */
enum { CIMPORT_STORE_NUMERIC = 0, CIMPORT_STORE_STRING = 1, CIMPORT_STORE_STRL = 2 };

typedef struct {
    CImportContext *ctx;
    const unsigned char *kind;       /* per column: CIMPORT_STORE_* */
    int chunk;
    size_t row_start, row_end;       /* indices into chunk->rows */
    ST_int obs_start;                /* Stata observation of row_start */
    ST_retcode error;
} CImportStoreTile;

static void *cimport_store_tile_worker(void *arg) {
    CImportStoreTile *task = (CImportStoreTile *)arg;
    CImportContext *ctx = task->ctx;
    CImportParsedChunk *chunk = &ctx->chunks[task->chunk];
    const int ncols = ctx->num_columns;
    const double missing = SV_missval;
    const char quote = cimport_extract_quote(ctx);
    ST_IIID store = (_stata_)->safestore;
    char buf[CTOOLS_MAX_STRING_LEN + 1];

    task->error = 0;
    for (size_t r = task->row_start; r < task->row_end; r++) {
        const CImportParsedRow *row = chunk->rows[r];
        const uint64_t *values = cimport_row_values(row);
        ST_int obs = task->obs_start + (ST_int)(r - task->row_start);
        int nfields = row->num_fields;
        CImportFieldRef field = {row->offset, 0};
        for (int col = 0; col < ncols; col++) {
            ST_retcode rc = 0;
            int kind = task->kind[col];
            if (col >= nfields) {
                if (kind == CIMPORT_STORE_NUMERIC) rc = store(col + 1, obs, missing);
                else if (kind == CIMPORT_STORE_STRING) rc = SF_sstore(col + 1, obs, "");
            } else {
                field.length = row->lengths[col];
                if (kind == CIMPORT_STORE_NUMERIC) {
                    double val;
                    if (values[col] != CIMPORT_NO_VALUE) {
                        memcpy(&val, &values[col], sizeof(val));
                    } else if (!cimport_context_number(ctx, &field, &val)) {
                        rc = atomic_load(&ctx->error_code);
                        val = missing;
                    }
                    if (!rc) rc = store(col + 1, obs, val);
                } else if (kind == CIMPORT_STORE_STRING) {
                    /* strL columns are transferred by the ado (cimport blob). */
                    const char *src = ctx->file_data + field.offset;
                    int len = (int)CIMPORT_FIELD_LENGTH(field);
                    if (ctx->retain_quotes || !(field.length & CIMPORT_FIELD_QUOTED_FLAG)) {
                        while (len > 0 && (src[len-1] == '\r' || src[len-1] == '\n')) len--;
                        if (len > CTOOLS_MAX_STRING_LEN) len = CTOOLS_MAX_STRING_LEN;
                        memcpy(buf, src, (size_t)len);
                        buf[len] = '\0';
                    } else {
                        cimport_extract_field_fast(ctx->file_data, &field, buf, sizeof(buf), quote);
                    }
                    rc = SF_sstore(col + 1, obs, buf);
                }
                field.offset += CIMPORT_FIELD_LENGTH(field) + 1;
            }
            if (rc) { task->error = rc; return (void *)1; }
        }
    }
    return NULL;
}

/* ============================================================================
 * LOAD Mode
 * ============================================================================ */

static ST_retcode cimport_do_load(const char *filename, char delimiter, bool has_header, int skip_rows, bool verbose, CImportBindQuotesMode bindquotes,
                                   CImportNumericTypeMode numeric_type_mode, char decimal_sep, char group_sep,
                                   CImportEmptyLinesMode emptylines_mode, int max_quoted_rows,
                                   CImportEncoding requested_encoding) {
    char msg[512];

    CImportContext *ctx = g_cimport_ctx;
    bool used_cache = (ctx != NULL && ctx->is_loaded && ctx->filename && strcmp(ctx->filename, filename) == 0);

    if (!used_cache) {
        cimport_clear_cached_context();

        /* Parse column type overrides BEFORE parsing CSV so inference runs once */
        int num_force_num = 0, num_force_str = 0;
        int *force_num = cimport_parse_col_list("CIMPORT_NUMCOLS", &num_force_num);
        int *force_str = cimport_parse_col_list("CIMPORT_STRCOLS", &num_force_str);

        g_cimport_ctx = cimport_parse_csv(filename, delimiter, has_header, skip_rows, verbose, bindquotes,
                                           numeric_type_mode, decimal_sep, group_sep, emptylines_mode, max_quoted_rows,
                                           requested_encoding,
                                           force_num, num_force_num, force_str, num_force_str);
        ctx = g_cimport_ctx;
        if (!ctx) {
            /* cimport_parse_csv took ownership of force_num/force_str */
            return g_cimport_parse_status;
        }
        cimport_display_warnings(ctx);
    }
    /* else: used_cache — force lists already set during scan, reuse them */

    if (!ctx->cache_ready) {
        /* Switch from MADV_SEQUENTIAL to MADV_NORMAL before multi-threaded
         * cache build + SPI store.  Multiple threads will read different
         * columns from overlapping pages; SEQUENTIAL would discard pages
         * behind one thread that another thread still needs. */
        cimport_madvise_normal(ctx);
        ctx->cache_ready = true;
    }

    double t_load_start = ctools_timer_ms();

    ST_int stata_nvars = SF_nvars();
    ST_int stata_nobs = SF_nobs();

    (void)verbose; /* verbose timing reported via scalars */

    if (stata_nvars < ctx->num_columns) {
        snprintf(msg, sizeof(msg), "Error: Stata has %d variables but CSV has %d columns\n",
                 stata_nvars, ctx->num_columns);
        cimport_display_error(msg);
        cimport_clear_cached_context();
        return 198;
    }

    if ((size_t)stata_nobs < ctx->total_rows) {
        snprintf(msg, sizeof(msg), "Error: Stata has %d observations but CSV has %zu rows\n",
                 stata_nobs, ctx->total_rows);
        cimport_display_error(msg);
        cimport_clear_cached_context();
        return 198;
    }

    /* Parallel SPI store over row tiles; strings are decoded directly from
     * the mapped file, so no intermediate string cache is built. */
    {
        unsigned char *kind = malloc((size_t)ctx->num_columns);
        size_t max_tiles = ctx->num_chunks + ctx->total_rows / CIMPORT_STORE_TILE_ROWS + 1;
        CImportStoreTile *tiles = ctools_safe_malloc2(max_tiles, sizeof(CImportStoreTile));
        if (!kind || !tiles) {
            free(kind); free(tiles);
            cimport_clear_cached_context();
            return 459;
        }

        ST_retcode store_error = 0;
        for (int col = 0; col < ctx->num_columns; col++) {
            ST_int var = col + 1;
            bool is_string = ctx->columns[col].type == CIMPORT_COL_STRING;
            kind[col] = !is_string ? CIMPORT_STORE_NUMERIC
                      : SF_var_is_strl(var) ? CIMPORT_STORE_STRL : CIMPORT_STORE_STRING;
            if ((SF_var_is_string(var) != 0) != is_string) store_error = 198;
        }

        size_t ntiles = 0, obs = 1;
        for (int c = 0; c < ctx->num_chunks && !store_error; c++) {
            CImportParsedChunk *chunk = &ctx->chunks[c];
            size_t r = (c == 0 && ctx->has_header) ? 1 : 0;
            for (; r < chunk->num_rows; r += CIMPORT_STORE_TILE_ROWS) {
                size_t end = r + CIMPORT_STORE_TILE_ROWS;
                if (end > chunk->num_rows) end = chunk->num_rows;
                tiles[ntiles] = (CImportStoreTile){ctx, kind, c, r, end, (ST_int)obs, 0};
                obs += end - r;
                ntiles++;
            }
        }
        if (!store_error && obs - 1 != ctx->total_rows) store_error = 459;

        if (!store_error) {
            ctools_persistent_pool *pool = ntiles > 1 ? ctools_get_global_pool() : NULL;
            if (pool == NULL || ctools_persistent_pool_submit_batch(pool, cimport_store_tile_worker,
                                    tiles, ntiles, sizeof(CImportStoreTile)) != 0) {
                for (size_t t = 0; t < ntiles; t++) cimport_store_tile_worker(&tiles[t]);
            } else {
                ctools_persistent_pool_wait(pool);
            }
            for (size_t t = 0; t < ntiles && !store_error; t++) store_error = tiles[t].error;
        }
        free(tiles);
        free(kind);
        if (store_error) {
            cimport_display_error("cimport: write failed; imported data are incomplete\n");
            cimport_clear_cached_context();
            return store_error;
        }
    }

    double t_load_end = ctools_timer_ms();
    double time_spi = t_load_end - t_load_start;

    /* Timing data saved to scalars for verbose table in .ado */

    /* Save timing and thread diagnostics to Stata scalars (convert ms to seconds) */
    SF_scal_save("_cimport_time_mmap", ctx->time_mmap / 1000.0);
    SF_scal_save("_cimport_time_parse", ctx->time_parse / 1000.0);
    SF_scal_save("_cimport_time_infer", ctx->time_type_infer / 1000.0);
    SF_scal_save("_cimport_time_cache", ctx->time_cache / 1000.0);
    SF_scal_save("_cimport_time_store", time_spi / 1000.0);
    SF_scal_save("_cimport_time_total", (ctx->time_mmap + ctx->time_parse + ctx->time_type_infer + ctx->time_cache + time_spi) / 1000.0);
    CTOOLS_SAVE_THREAD_INFO("_cimport");

    cimport_clear_cached_context();

    return 0;
}

/* ============================================================================
 * Plugin Entry Point
 * ============================================================================ */

/* Structure to hold parsed options for cimport */
typedef struct {
    char *mode;
    char *filename;
    char delimiter;
    bool has_header;
    int header_row;         /* 1-based row number for header (1 = first row, 0 = no header) */
    bool verbose;
    CImportBindQuotesMode bindquotes;
    CImportNumericTypeMode numeric_type_mode;
    char decimal_separator;
    char group_separator;
    CImportEmptyLinesMode emptylines_mode;
    int max_quoted_rows;
    CImportEncoding encoding;
} CImportOptions;

static ST_retcode cimport_blob(int col, int row, size_t offset) {
    CImportContext *ctx=g_cimport_ctx;
    if (!ctx || col<1 || col>ctx->num_columns || row<1 || (size_t)row>ctx->total_rows ||
        ctx->columns[col-1].type!=CIMPORT_COL_STRING) return 198;
    cimport_build_column_cache(ctx);
    if(!ctx->cache_ready)return atomic_load(&ctx->error_code) ? atomic_load(&ctx->error_code) : 920;
    const unsigned char *text=(const unsigned char*)ctx->col_cache[col-1].string_data[row-1];
    size_t length=strlen((const char*)text),take;
    if(offset>length)return 198;
    take=length-offset;if(take>8000)take=8000;
    char hex[16001],next[40];
    const char *digits="0123456789abcdef";
    for(size_t i=0;i<take;i++){hex[2*i]=digits[text[offset+i]>>4];hex[2*i+1]=digits[text[offset+i]&15];}
    hex[2*take]=0;
    snprintf(next,sizeof(next),"%zu",offset+take<length?offset+take:0);
    int rc=SF_macro_save("_cimport_blobhex",hex);
    return rc ? rc : SF_macro_save("_cimport_blobnext",next);
}

ST_retcode cimport_main(const char *args) {
    if (args && !strcmp(args,"clear")) { cimport_clear_cached_context();return 0; }
    if (args && !strncmp(args,"blob ",5)) {
        int col,row;size_t offset;
        if(sscanf(args,"%*s %d %d %zu",&col,&row,&offset)!=3)return 198;
        return cimport_blob(col,row,offset);
    }
    if (args == NULL || strlen(args) == 0) {
        cimport_display_error("cimport: no arguments specified\n");
        return 198;
    }

    /* Parse arguments: mode filename [delimiter] [options...] */
    char args_copy[4096];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    CImportOptions opts = {
        .mode = NULL,
        .filename = NULL,
        .delimiter = '\0',  /* '\0' = auto-detect */
        .has_header = true,
        .header_row = 1,      /* Default: first row is header */
        .verbose = false,
        .bindquotes = CIMPORT_BINDQUOTES_LOOSE,
        .numeric_type_mode = CIMPORT_NUMTYPE_AUTO,
        .decimal_separator = '.',
        .group_separator = '\0',
        .emptylines_mode = CIMPORT_EMPTYLINES_SKIP,
        .max_quoted_rows = 20,
        .encoding = CIMPORT_ENC_UNKNOWN  /* Auto-detect by default */
    };

    char *token = strtok(args_copy, " ");
    int arg_idx = 0;

    while (token != NULL) {
        if (arg_idx == 0) {
            opts.mode = token;
        } else if (arg_idx == 1) {
            opts.filename = token;
        } else {
            if (strcmp(token, "noheader") == 0) {
                opts.has_header = false;
                opts.header_row = 0;
            } else if (strncmp(token, "headerrow=", 10) == 0) {
                if (!ctools_safe_atoi(token + 10, &opts.header_row)) {
                    cimport_display_error("cimport: invalid headerrow value\n");
                    return 198;
                }
                opts.has_header = (opts.header_row > 0);
            } else if (strcmp(token, "verbose") == 0) {
                opts.verbose = true;
            } else if (strcmp(token, "tab") == 0) {
                opts.delimiter = '\t';
            } else if (strcmp(token, "space") == 0) {
                opts.delimiter = ' ';
            } else if (strcmp(token, "auto") == 0) {
                opts.delimiter = '\0';  /* auto-detect */
            } else if (strcmp(token, "bindquotes=strict") == 0) {
                opts.bindquotes = CIMPORT_BINDQUOTES_STRICT;
            } else if (strcmp(token, "bindquotes=loose") == 0) {
                opts.bindquotes = CIMPORT_BINDQUOTES_LOOSE;
            } else if (strcmp(token, "bindquotes=nobind") == 0) {
                opts.bindquotes = CIMPORT_BINDQUOTES_NOBIND;
            } else if (strcmp(token, "asfloat") == 0) {
                opts.numeric_type_mode = CIMPORT_NUMTYPE_FLOAT;
            } else if (strcmp(token, "asdouble") == 0) {
                opts.numeric_type_mode = CIMPORT_NUMTYPE_DOUBLE;
            } else if (strncmp(token, "decimalsep=", 11) == 0) {
                opts.decimal_separator = token[11];
            } else if (strncmp(token, "groupsep=", 9) == 0) {
                /* Handle "space" keyword for space character */
                if (strcmp(token + 9, "space") == 0) {
                    opts.group_separator = ' ';
                } else {
                    opts.group_separator = token[9];
                }
            } else if (strcmp(token, "emptylines=fill") == 0) {
                opts.emptylines_mode = CIMPORT_EMPTYLINES_FILL;
            } else if (strncmp(token, "maxquotedrows=", 14) == 0) {
                if (!ctools_safe_atoi(token + 14, &opts.max_quoted_rows)) {
                    cimport_display_error("cimport: invalid maxquotedrows value\n");
                    return 198;
                }
            } else if (strncmp(token, "encoding=", 9) == 0) {
                opts.encoding = cimport_parse_encoding_name(token + 9);
            } else if (strlen(token) == 1) {
                opts.delimiter = token[0];
            } else if (strlen(token) == 3 && token[0] == '"' && token[2] == '"') {
                opts.delimiter = token[1];
            }
        }
        arg_idx++;
        token = strtok(NULL, " ");
    }

    if (opts.mode == NULL || opts.filename == NULL) {
        cimport_display_error("cimport: mode and filename required\n");
        return 198;
    }

    /* The ado transports the filename separately so spaces/quotes never become
     * command tokens. Keep literal tokens for direct plugin callers. */
    char filename_buf[4096];
    if (strcmp(opts.filename, "@filename") == 0) {
        if (SF_macro_use("___cimport_filename", filename_buf, sizeof(filename_buf) - 1) || !filename_buf[0]) {
            cimport_display_error("cimport: cannot read filename\n");
            return 198;
        }
        opts.filename = filename_buf;
    }

    /* Calculate skip_rows from header_row (header_row=1 means no skip, header_row=3 means skip 2) */
    int skip_rows = (opts.header_row > 1) ? (opts.header_row - 1) : 0;

    if (strcmp(opts.mode, "scan") == 0) {
        return cimport_do_scan(opts.filename, opts.delimiter, opts.has_header, skip_rows, opts.verbose, opts.bindquotes,
                                opts.numeric_type_mode, opts.decimal_separator, opts.group_separator,
                                opts.emptylines_mode, opts.max_quoted_rows, opts.encoding);
    } else if (strcmp(opts.mode, "load") == 0) {
        return cimport_do_load(opts.filename, opts.delimiter, opts.has_header, skip_rows, opts.verbose, opts.bindquotes,
                                opts.numeric_type_mode, opts.decimal_separator, opts.group_separator,
                                opts.emptylines_mode, opts.max_quoted_rows, opts.encoding);
    } else {
        cimport_display_error("cimport: invalid mode. Use 'scan' or 'load'\n");
        return 198;
    }
}

/* ============================================================================
 * Cleanup function for ctools_cleanup system
 * ============================================================================ */

void cimport_cleanup_cache(void)
{
    cio_cleanup();
    cimport_clear_cached_context();
}

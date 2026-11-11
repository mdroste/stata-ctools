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
#include "ctools_file.h"
#include "formats/ctools_dta.h"

/* High-resolution timing */
/* Timer: use ctools_timer_ms() directly from ctools_timer.h */

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/* Local: cimport may use more threads for parsing than general I/O */
#define CIMPORT_MAX_THREADS      32

/* Small inputs keep the legacy row cache. Larger inputs use typed blocks and
 * retain a cached prefix within a hard payload budget; only the suffix is reparsed.
 * Overridable only at build time for differential tests. */
#ifndef CIMPORT_STREAM_THRESHOLD
#define CIMPORT_STREAM_THRESHOLD ((size_t)256 * 1024 * 1024)
#endif

/* Below twice this size an input is parsed in one piece. */
#define CIMPORT_MIN_CHUNK ((size_t)1 << 20)

/* With r(109): a field of CIMPORT_FIELD_MAX_LENGTH bytes or more, often an
 * unclosed quote under bindquotes(strict) that reaches far into the file. */
#define CIMPORT_FIELD_TOO_LONG_MSG \
    "cimport: a field is 2 GB or longer, more than a Stata string holds " \
    "(an unclosed quote under bindquotes(strict)?)\n"

/* Parse chunks follow the thread limit, so threads() and the default cap apply. */
static int cimport_parse_threads(void) {
    int n = ctools_get_max_threads();
    if (n <= 0) n = 1;
    return n > CIMPORT_MAX_THREADS ? CIMPORT_MAX_THREADS : n;
}

/* Global cached context */
static CImportContext *g_cimport_ctx = NULL;
static ST_retcode g_cimport_parse_status = 601;

/* A failed scan after the file opened: report the worker's error (920 for
 * memory, 109 for an oversized field) or the fallback, never 601. */
static void cimport_free_context(CImportContext *ctx);
static void cimport_blobmem_clear(CImportContext *ctx);
static CImportContext *cimport_fail(CImportContext *ctx, ST_retcode fallback) {
    ST_retcode rc = (ST_retcode)atomic_load(&ctx->error_code);
    g_cimport_parse_status = rc ? rc : fallback;
    cimport_free_context(ctx);
    return NULL;
}

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

/* Record a warning for an unmatched quote in the row starting at byte offset
 * (thread-safe). The first CIMPORT_MAX_WARNINGS rows of the file are kept, in
 * file order, whatever order the parallel chunks find them in; the display
 * turns offsets into line numbers. */
static void cimport_record_unmatched_quote(CImportContext *ctx,
                                           size_t offset) {
  pthread_mutex_lock(&ctx->warning_mutex);
  int kept = ctx->num_unmatched_quote_warnings < CIMPORT_MAX_WARNINGS
                 ? ctx->num_unmatched_quote_warnings
                 : CIMPORT_MAX_WARNINGS;
  if (kept < CIMPORT_MAX_WARNINGS || offset < ctx->unmatched_quote_rows[kept - 1]) {
    int i = kept < CIMPORT_MAX_WARNINGS ? kept : kept - 1;
    while (i > 0 && ctx->unmatched_quote_rows[i - 1] > offset) {
      ctx->unmatched_quote_rows[i] = ctx->unmatched_quote_rows[i - 1];
      i--;
    }
    ctx->unmatched_quote_rows[i] = offset;
  }
  ctx->num_unmatched_quote_warnings++;
  pthread_mutex_unlock(&ctx->warning_mutex);
}

static size_t cimport_count_byte(const char *p, const char *end, char byte);

/* The start of the row holding file line n_lines + 1: varnames(#) names the
 * row on line # as native counts lines, so under bindquotes(strict) a field
 * spanning lines above it counts every line, and a row that spans line # is
 * the header. Rows that end before that line are skipped. With loose or
 * nobind quotes every row is one line, so this skips n_lines rows. */
static const char *cimport_skip_rows(const char *data, const char *end, int n_lines, char quote_char, CImportBindQuotesMode bindquotes) {
    const char *ptr = data;
    size_t lines = 0, target = n_lines > 0 ? (size_t)n_lines : 0;
    while (ptr < end && lines < target) {
        const char *next = cimport_find_next_row(ptr, end, quote_char, bindquotes);
        size_t in_row = cimport_count_byte(ptr, next, '\n');
        if (lines + in_row > target) break; /* this row holds line target + 1 */
        lines += in_row;
        ptr = next;
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
/* colrange() as file columns [col_first, cimport_col_stop) */
static int cimport_col_stop(const CImportContext *ctx) {
    return ctx->col_last > 0 && ctx->col_last < ctx->num_columns ? ctx->col_last : ctx->num_columns;
}

static bool cimport_col_selected(const CImportContext *ctx, int col) {
    return col >= ctx->col_first && col < cimport_col_stop(ctx);
}

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

/* Parse one row into *fields, doubling the buffer until the whole row fits:
 * native reads rows of any width (colrange() and maxvar decide what loads).
 * The buffer starts at CTOOLS_MAX_COLUMNS refs. Returns -1 if it cannot grow. */
static int cimport_parse_row_grow(const CImportContext *ctx, const char *row, const char *end,
                                  CImportFieldRef **fields, int *capacity,
                                  const char **next_row, bool *had_unmatched) {
    for (;;) {
        int count = cimport_parse_row_fast(row, end, ctx->delimiter, ctx->quote_char, *fields, *capacity,
                                           ctx->file_data, ctx->bindquotes, next_row, had_unmatched);
        if (count < *capacity) return count;
        if (*capacity > INT_MAX / 2) return -1;
        CImportFieldRef *grown = realloc(*fields, 2 * (size_t)*capacity * sizeof(**fields));
        if (!grown) return -1;
        *fields = grown;
        *capacity *= 2;
    }
}

/* Display warnings for unmatched quotes (Stata-compatible format) */
/* emptylines(include): native counts the empty text after a final line
 * terminator as one more empty line, so the row parsed at that terminator is
 * kept as an empty observation. Only for the part of the file that ends it. */
static bool cimport_trailing_empty_line(const CImportContext *ctx, const char *end) {
    return ctx->emptylines_mode != CIMPORT_EMPTYLINES_SKIP && ctx->file_size > 0 &&
           end == ctx->file_data + ctx->file_size && (end[-1] == '\n' || end[-1] == '\r');
}

static void cimport_display_warnings(CImportContext *ctx) {
    if (ctx->num_unmatched_quote_warnings == 0) return;

    char msg[600];
    int to_display = (ctx->num_unmatched_quote_warnings < CIMPORT_MAX_WARNINGS)
                     ? ctx->num_unmatched_quote_warnings
                     : CIMPORT_MAX_WARNINGS;

    /* Native names the file line the row starts on: blank lines and the
     * rows varnames(#) skips count. Offsets are kept in file order. */
    CImportLineCursor cursor = {ctx->file_data, 0};
    for (int i = 0; i < to_display; i++) {
        size_t line = ctx->file_data
            ? cimport_line_at(&cursor, ctx->file_data + ctx->unmatched_quote_rows[i]) : 0;
        snprintf(msg, sizeof(msg),
            "Note: Unmatched quote while processing row %zu; this can be due to a formatting\n"
            "    problem in the file or because a quoted data element spans multiple\n"
            "    lines. You should carefully inspect your data after importing. Consider\n"
            "    using option bindquote(strict) if quoted data spans multiple lines or\n"
            "    option bindquote(nobind) if quotes are not used for binding data.\n",
            line);
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
    /* As native import does, a column with non-integers is double only when
     * one of its integer values (every value beyond 2^53 is one) leaves the
     * long range; non-integers alone keep float storage. */
    double low = col->is_integer ? col->min_value : col->int_min_value;
    double high = col->is_integer ? col->max_value : col->int_max_value;
    if (low < -2147483647.0 || high > 2147483620.0) {
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
static int cimport_parse_header_row(const CImportContext *ctx, CImportFieldRef **fields, int *capacity,
                                    const char **row_start, const char **next_row) {
    const char *end = ctx->file_data + ctx->file_size;
    const char *row = cimport_skip_rows(ctx->file_data, end, ctx->skip_rows, ctx->quote_char, ctx->bindquotes);
    const char *next = end;
    int count = 0;
    while (row < end) {
        count = cimport_parse_row_grow(ctx, row, end, fields, capacity, &next, NULL);
        if (count < 0) return count == CIMPORT_FIELD_TOO_LONG ? count : -1;
        if (ctx->emptylines_mode != CIMPORT_EMPTYLINES_SKIP ||
            !cimport_row_is_empty(ctx->file_data, *fields, count)) break;
        row = next;
        count = 0;
    }
    if (row_start) *row_start = row;
    if (next_row) *next_row = row < end ? next : end;
    return count;
}

static int cimport_parse_header(CImportContext *ctx) {
    /* The header's field refs (512 KiB, grown for wider rows) live on the
     * heap, off the Stata thread's stack. */
    int capacity = CTOOLS_MAX_COLUMNS;
    CImportFieldRef *fields = ctools_safe_malloc2((size_t)capacity, sizeof(CImportFieldRef));
    const char *end = ctx->file_data + ctx->file_size;
    const char *header_start, *next;
    if (!fields) return -1;

    ctx->num_columns = cimport_parse_header_row(ctx, &fields, &capacity, &header_start, &next);
    if (ctx->num_columns == CIMPORT_FIELD_TOO_LONG) {
        cimport_display_error(CIMPORT_FIELD_TOO_LONG_MSG);
        atomic_store(&ctx->error_code, 109);
    }
    if (ctx->num_columns < 0) { ctx->num_columns = 0; free(fields); return -1; }
    /* Chunked parsing starts at this row, so chunk 0 always holds the header. */
    ctx->first_row_offset = (size_t)(header_start - ctx->file_data);
    /* Only empty rows: native imports nothing, like an empty file. */
    if (ctx->num_columns == 0) { free(fields); return 0; }

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
        /* Only the header's columns are inspected; the parser skips any
         * further fields of a wider row. */
        CImportFieldRef *next_fields = ctools_safe_malloc2((size_t)ctx->num_columns + 1, sizeof(CImportFieldRef));
        bool *candidates = calloc((size_t)ctx->num_columns, sizeof(bool));
        if (!candidates || !next_fields) { free(candidates); free(next_fields); free(fields); return -1; }
        for (int i = 0; i < ctx->num_columns; i++) {
            int len = CIMPORT_FIELD_LENGTH(fields[i]);
            candidates[i] = len && cimport_col_selected(ctx, i) && !cimport_context_numeric(ctx,&fields[i]);
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
                        next_fields, ctx->num_columns + 1, ctx->file_data, ctx->bindquotes, &next, NULL);
            if (count < 0) break; /* a field too long to sample; the parse reports it */
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
        free(candidates);
        free(next_fields);
    }

    ctx->columns = ctools_safe_calloc2(ctx->num_columns, sizeof(CImportColumnInfo));
    if (!ctx->columns) { free(fields); return -1; }

    char name_buf[CTOOLS_MAX_STRING_LEN];
    for (int i = 0; i < ctx->num_columns; i++) {
        if (ctx->has_header) {
            cimport_extract_field_fast(ctx->file_data, &fields[i], name_buf, sizeof(name_buf), cimport_extract_quote(ctx));
            cimport_sanitize_varname(name_buf, strlen(name_buf), ctx->columns[i].name);
            if (!ctx->columns[i].name[0])
                snprintf(ctx->columns[i].name, sizeof(ctx->columns[i].name), "v%d", i + 1 - ctx->col_first);
        } else {
            snprintf(ctx->columns[i].name, CTOOLS_MAX_VARNAME_LEN, "v%d", i + 1 - ctx->col_first);
        }
        ctx->columns[i].type = CIMPORT_COL_UNKNOWN;
        ctx->columns[i].max_strlen = 0;
    }

    free(fields);
    return 0;
}

/* Helper to check if a column is in the force list */
static bool cimport_col_in_list(int col_1based, int *list, int count) {
    for (int i = 0; i < count; i++) {
        if (list[i] == 0 || list[i] == col_1based) return true;
    }
    return false;
}

/* Numeric overrides take precedence over string overrides, as in native. */
static int cimport_column_force(const CImportContext *ctx, int col) {
    if (!cimport_col_selected(ctx, col)) return CIMPORT_FORCE_NONE;
    col -= ctx->col_first;
    if (ctx->force_numeric_cols && cimport_col_in_list(col + 1, ctx->force_numeric_cols, ctx->num_force_numeric))
        return CIMPORT_FORCE_NUMERIC;
    if (ctx->force_string_cols && cimport_col_in_list(col + 1, ctx->force_string_cols, ctx->num_force_string))
        return CIMPORT_FORCE_STRING;
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
        stats[i].int_min = DBL_MAX;
        stats[i].int_max = -DBL_MAX;
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
    else {
        if (val < stats->int_min) stats->int_min = val;
        if (val > stats->int_max) stats->int_max = val;
    }
    if (val < stats->num_min) stats->num_min = val;
    if (val > stats->num_max) stats->num_max = val;
    return cached;
}

static void cimport_infer_column_types(CImportContext *ctx) {
    size_t total_data_rows = 0;
    ctx->strl_bytes = 0;
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
        colinfo->int_min_value = DBL_MAX;
        colinfo->int_max_value = -DBL_MAX;
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
            if (stats->int_min < colinfo->int_min_value) colinfo->int_min_value = stats->int_min;
            if (stats->int_max > colinfo->int_max_value) colinfo->int_max_value = stats->int_max;
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
            if (colinfo->use_strl && cimport_col_selected(ctx, col)) ctx->strl_bytes += textbytes;
            continue;
        }

        colinfo->type = CIMPORT_COL_NUMERIC;
        if (!colinfo->is_integer && ctx->numeric_type_mode == CIMPORT_NUMTYPE_DOUBLE) {
            colinfo->num_subtype = CIMPORT_NUM_DOUBLE;
        } else if (colinfo->min_value == DBL_MAX) {
            /* No number in the column. Text that numericcols() forces to
             * missing takes the type native gives it, by the column's string
             * width: 1 byte, 2 int, 3 long, 4 float, 5 double, wider byte.
             * (Native is not deterministic when text and numbers mix.) */
            static const int by_width[] = {CIMPORT_NUM_BYTE, CIMPORT_NUM_BYTE, CIMPORT_NUM_INT,
                                           CIMPORT_NUM_LONG, CIMPORT_NUM_FLOAT, CIMPORT_NUM_DOUBLE};
            colinfo->num_subtype = force == CIMPORT_FORCE_NUMERIC && is_string && decoded_len <= 5
                ? by_width[decoded_len] : CIMPORT_NUM_BYTE;
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
    cimport_blobmem_clear(ctx);

    if (ctx->chunks) {
        for (int c = 0; c < ctx->num_chunks; c++) {
            CImportParsedChunk *chunk = &ctx->chunks[c];
            free(chunk->rows);
            free(chunk->col_stats);
            ctools_arena_free(&chunk->arena);
            cimport_cache_free(chunk->cache, &ctx->cache_budget);
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

/* "UTF-8" or "UTF8" in any case: names both converters accept for UTF-8. */
static bool cimport_names_utf8(const char *name) {
    char lower[6];
    size_t n = strlen(name);
    if (n > 5) return false;
    for (size_t i = 0; i <= n; i++) lower[i] = (char)tolower((unsigned char)name[i]);
    return !strcmp(lower, "utf-8") || !strcmp(lower, "utf8");
}

/* What reading the input bytes in place needs to know of them: whether any
 * is above 0x7f, whether they are well-formed UTF-8, whether any is zero, and
 * whether a CR ends a line by itself (one not followed by LF). Large inputs
 * are scanned in parallel segments that start at a byte which does not
 * continue a UTF-8 sequence, so each sequence is checked whole by one
 * segment (four continuation bytes in a row are never valid). */
typedef struct {
    const char *start;
    size_t size;
    bool high, valid, nul, lone_cr;
    const char *limit; /* end of the input */
} CImportByteScan;

static void *cimport_byte_scan_task(void *arg) {
    CImportByteScan *scan = arg;
    const char *end = scan->start + scan->size;
    /* A segment starts invalid only after four continuation bytes. */
    if (scan->valid) scan->valid = cimport_scan_utf8(scan->start, scan->size, &scan->high);
    else scan->high = true;
    scan->nul = scan->size && memchr(scan->start, 0, scan->size) != NULL;
    scan->lone_cr = false;
    for (const char *cr = scan->start; cr < end && (cr = memchr(cr, '\r', (size_t)(end - cr))); cr++)
        if (cr + 1 == scan->limit || cr[1] != '\n') { scan->lone_cr = true; break; }
    return NULL;
}

/* Stats of one chunk for the columns a wider row added after it was parsed. */
typedef struct {
    CImportContext *ctx;
    int chunk;
    int new_num;
    int rc;
} CImportExpandTask;

static void *cimport_expand_task(void *arg) {
    CImportExpandTask *task = arg;
    CImportContext *ctx = task->ctx;
    CImportParsedChunk *chunk = &ctx->chunks[task->chunk];
    const int new_num = task->new_num;
    task->rc = 0;
    if (chunk->num_col_stats >= new_num) return NULL;
    CImportColumnParseStats *new_stats = realloc(chunk->col_stats, (size_t)new_num * sizeof(*new_stats));
    if (!new_stats) { task->rc = 920; return NULL; }
    cimport_init_col_stats(new_stats + chunk->num_col_stats, new_num - chunk->num_col_stats);
    chunk->col_stats = new_stats;
    size_t start_row = (task->chunk == 0 && ctx->has_header) ? 1 : 0;
    for (size_t r = start_row; r < chunk->num_rows; r++) {
        CImportParsedRow *row = chunk->rows[r];
        if (row->num_fields <= chunk->num_col_stats) continue;
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
    return NULL;
}

static CImportByteScan cimport_scan_bytes(const char *data, size_t size) {
    enum { SEGMENT_MIN = 8 << 20 };
    CImportByteScan part[CIMPORT_MAX_THREADS];
    int threads = ctools_get_max_threads();
    if (threads > CIMPORT_MAX_THREADS) threads = CIMPORT_MAX_THREADS;
    int n = threads > 1 ? (int)(size / SEGMENT_MIN) : 1;
    if (n > threads) n = threads;
    if (n < 1) n = 1;
    const unsigned char *s = (const unsigned char *)data;
    size_t at = 0;
    for (int i = 0; i < n; i++) {
        size_t end = i == n - 1 ? size : (size / (size_t)n) * (size_t)(i + 1), moved = 0;
        while (end < size && (s[end] & 0xc0) == 0x80 && moved < 4) end++, moved++;
        if (end < at) end = at;
        part[i] = (CImportByteScan){data + at, end - at, false, moved < 4, false, false, data + size};
        at = end;
    }
    ctools_persistent_pool *pool = n > 1 ? ctools_get_global_pool() : NULL;
    if (pool == NULL || ctools_persistent_pool_submit_batch(pool, cimport_byte_scan_task, part, n,
                                                            sizeof(CImportByteScan)) != 0) {
        for (int i = 0; i < n; i++) cimport_byte_scan_task(&part[i]);
    } else {
        ctools_persistent_pool_wait(pool);
    }
    CImportByteScan all = {data, size, false, true, false, false, data + size};
    for (int i = 0; i < n; i++) {
        all.high |= part[i].high;
        all.valid &= part[i].valid;
        all.nul |= part[i].nul;
        all.lone_cr |= part[i].lone_cr;
    }
    return all;
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

/* Limit typed block allocations to a quarter of physical RAM, leaving room for
 * mapped input, Stata's data and its rollback copy. Actual block bytes, not a
 * sampled estimate, consume this budget. Types still use every input value. */
static uint64_t cimport_cache_limit(void)
{
#ifdef CIMPORT_CACHE_BUDGET
    return CIMPORT_CACHE_BUDGET;
#else
    uint64_t memory = ctools_physical_memory_bytes();
    return memory ? memory / 4 : UINT64_C(512) * 1024 * 1024;
#endif
}

static bool cimport_cache_eligible(CImportContext *ctx)
{
    if (ctx->file_size < CIMPORT_STREAM_THRESHOLD) return false;
    size_t first, last;
    cimport_row_range(&first, &last);
    return !first && !last;
}

static bool cimport_should_stream(CImportContext *ctx)
{
    return cimport_cache_eligible(ctx) && !cimport_cache_limit();
}

static bool cimport_flush_cache(CImportContext *ctx, CImportParsedChunk *chunk,
                                CImportCacheBuilder *builder)
{
    for (int col = 0; col < ctx->num_columns; col++)
        builder->numeric[col] = ctx->col_force[col] == CIMPORT_FORCE_NUMERIC ||
            (ctx->col_force[col] != CIMPORT_FORCE_STRING && !chunk->col_stats[col].seen_string);
    CImportCachedBlock *block = cimport_cache_pack(builder, &ctx->cache_budget, SV_missval,
        ctx->numeric_type_mode == CIMPORT_NUMTYPE_FLOAT,
        ctx->numeric_type_mode == CIMPORT_NUMTYPE_DOUBLE,
        ctx->file_data, cimport_extract_quote(ctx));
    if (!block) { cimport_cache_destroy(builder); return false; }
    if (chunk->cache_tail) chunk->cache_tail->next = block;
    else chunk->cache = block;
    chunk->cache_tail = block;
    chunk->cached_rows += block->rows;
    chunk->uncached_start = block->file_end;
    return true;
}

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
    chunk->file_start = (size_t)(task->start - ctx->file_data);
    chunk->file_end = (size_t)(task->end - ctx->file_data);
    chunk->uncached_start = chunk->file_start;
    chunk->capacity = (ctx->stream_rows || ctx->compact_rows) ? 1 : 16384;
    chunk->rows = ctools_safe_malloc2(chunk->capacity, sizeof(CImportParsedRow *));
    chunk->num_rows = 0;

    chunk->num_col_stats = ctx->num_columns;
    chunk->col_stats = ctools_safe_calloc2(ctx->num_columns, sizeof(CImportColumnParseStats));

    if (!chunk->rows || !chunk->col_stats) {
        atomic_store(&ctx->error_code, 920);
        return NULL;
    }
    cimport_init_col_stats(chunk->col_stats, ctx->num_columns);

    /* 512 KiB of field refs (grown for wider rows): on the heap, as a worker
     * stack may be no larger. */
    int field_capacity = CTOOLS_MAX_COLUMNS;
    CImportFieldRef *field_buf = ctools_safe_malloc2((size_t)field_capacity, sizeof(CImportFieldRef));
    if (!field_buf) {
        atomic_store(&ctx->error_code, 920);
        return NULL;
    }
    CImportCacheBuilder builder = {0};
    double *row_values = NULL;
    bool caching = ctx->compact_rows && !cimport_cache_init(&builder, ctx->num_columns);
    if (caching) {
        row_values = ctools_safe_malloc2((size_t)ctx->num_columns, sizeof(double));
        if (!row_values) { cimport_cache_destroy(&builder); caching = false; }
    }
    const char *ptr = task->start;
    const char *end = task->end;
    bool trailing_empty = task->chunk_id == ctx->num_chunks - 1 && cimport_trailing_empty_line(ctx, end);

    while ((ptr < end || trailing_empty) && atomic_load(&ctx->error_code) == 0) {
        if (ptr >= end) {
            ptr = end - 1;
            trailing_empty = false;
        }
        const char *next_row;
        bool had_unmatched;
        int num_fields = cimport_parse_row_grow(ctx, ptr, end, &field_buf, &field_capacity,
                                                &next_row, &had_unmatched);
        if (num_fields < 0) {
            atomic_store(&ctx->error_code, num_fields == CIMPORT_FIELD_TOO_LONG ? 109 : 920);
            break;
        }
        /* Only a row with a line end inside it can hold a quoted field that
         * spans lines. */
        if (ctx->quoted_rows_in_parse && next_row - ptr > 1 &&
            memchr(ptr, '\n', (size_t)(next_row - ptr - 1)) &&
            cimport_quoted_rows_exceed(ptr, next_row, (size_t)ctx->max_quoted_rows)) {
            atomic_store(&ctx->error_code, 5101);
            break;
        }

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
                cimport_record_unmatched_quote(ctx, (size_t)(ptr - ctx->file_data));
            }
        }


        if (!ctx->stream_rows && !ctx->compact_rows && chunk->num_rows >= chunk->capacity) {
            if (chunk->capacity > SIZE_MAX / 2) {
                atomic_store(&ctx->error_code, 920);
                break;
            }
            size_t new_capacity = chunk->capacity * 2;
            if (new_capacity > SIZE_MAX / sizeof(CImportParsedRow *)) {
                atomic_store(&ctx->error_code, 920);
                break;
            }
            CImportParsedRow **new_rows = realloc(chunk->rows, sizeof(CImportParsedRow *) * new_capacity);
            if (!new_rows) {
                atomic_store(&ctx->error_code, 920);
                break;
            }
            chunk->rows = new_rows;
            chunk->capacity = new_capacity;
        }

        if ((ctx->stream_rows || ctx->compact_rows) && !is_header_row) {
            for (int f = 0; f < num_fields && f < chunk->num_col_stats; f++) {
                uint64_t bits = cimport_accumulate_field(ctx, &chunk->col_stats[f], &field_buf[f], ctx->col_force[f]);
                if (caching) {
                    double value = SV_missval;
                    if (bits != CIMPORT_NO_VALUE) memcpy(&value, &bits, sizeof(value));
                    else if ((ctx->col_force[f] == CIMPORT_FORCE_NUMERIC ||
                              !chunk->col_stats[f].seen_string) &&
                             !cimport_context_number(ctx, &field_buf[f], &value)) value = SV_missval;
                    row_values[f] = value;
                }
            }
            if (caching) {
                for (int f = num_fields; f < ctx->num_columns; f++) row_values[f] = SV_missval;
                cimport_cache_append(&builder, field_buf, num_fields, row_values,
                    (uint64_t)(ptr - ctx->file_data), (uint64_t)(next_row - ctx->file_data));
                if (builder.rows == builder.capacity) caching = cimport_flush_cache(ctx, chunk, &builder);
            }
        } else {
            CImportParsedRow *row = ctools_arena_alloc(&chunk->arena, cimport_row_size(num_fields));
            if (!row || !cimport_fill_row(row, field_buf, num_fields)) {
                atomic_store(&ctx->error_code, 920);
                break;
            }
            uint64_t *values = cimport_row_values(row);
            for (int f = 0; f < num_fields; f++) {
                values[f] = !is_header_row && f < chunk->num_col_stats
                    ? cimport_accumulate_field(ctx, &chunk->col_stats[f], &field_buf[f], ctx->col_force[f])
                    : CIMPORT_NO_VALUE;
            }
            chunk->rows[chunk->num_rows] = row;
        }
        chunk->num_rows++;

        if (num_fields > chunk->max_fields_in_chunk) {
            chunk->max_fields_in_chunk = num_fields;
        }

        ptr = next_row;
    }

    if (caching && builder.rows && !atomic_load(&ctx->error_code))
        cimport_flush_cache(ctx, chunk, &builder);
    cimport_cache_destroy(&builder);
    free(row_values);
    free(field_buf);
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
            /* A line past rowrange()'s end is not sampled, also when blank
             * lines before it were skipped. */
            if (last && lines > last) break;
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
    if (!ctx) {
        /* This function owns the force-column lists on every return. */
        free(force_numeric_cols);
        free(force_string_cols);
        g_cimport_parse_status = 920;
        return NULL;
    }

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
    if (ctx->filename == NULL) return cimport_fail(ctx, 920);

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
    /* A Java charset name or alias given in encoding() is reported by its
     * charset's name (latin1 as ISO-8859-1). */
    const char *java_name =
        explicit_encoding ? cimport_java_charset(source_encoding) : NULL;
    snprintf(ctx->reported_encoding, sizeof(ctx->reported_encoding), "%s",
             java_name ? java_name : source_encoding);
    ctx->encoding = cimport_parse_encoding_name(source_encoding);
    CImportByteScan bytes = cimport_scan_bytes(ctx->file_data, ctx->file_size);
    /* These stateless encodings preserve every ASCII byte. The detector can
     * select any of its ISO-8859 models even for entirely numeric ASCII text;
     * keep that reported charset without needlessly transcoding the file.
     * Explicit names still go through converter validation. */
    bool ascii_fast_path =
        !explicit_encoding &&
        (ctx->encoding == CIMPORT_ENC_UTF8 ||
         ctx->encoding == CIMPORT_ENC_ASCII ||
         ctx->encoding == CIMPORT_ENC_ISO_8859_1 ||
         ctx->encoding == CIMPORT_ENC_ISO_8859_15 ||
         ctx->encoding == CIMPORT_ENC_WINDOWS_1252 ||
         ctx->encoding == CIMPORT_ENC_MACROMAN ||
         !strcmp(source_encoding, "ISO-8859-2") ||
         !strcmp(source_encoding, "ISO-8859-5") ||
         !strcmp(source_encoding, "ISO-8859-6") ||
         !strcmp(source_encoding, "ISO-8859-7") ||
         !strcmp(source_encoding, "ISO-8859-8") ||
         !strcmp(source_encoding, "ISO-8859-9") ||
         !strcmp(source_encoding, "KOI8-R")) &&
        !bytes.high;
    /* Well-formed UTF-8 named or detected as UTF-8 converts to itself less a
     * leading BOM, so the mapped file is used in place rather than copied. */
    bool utf8_in_place =
        !ascii_fast_path &&
        cimport_names_utf8(source_encoding) && bytes.valid;
    if (utf8_in_place && ctx->file_size >= 3 &&
        !memcmp(ctx->file_data, "\xef\xbb\xbf", 3)) {
      ctx->file_data += 3;
      ctx->file_size -= 3;
    }
    if (!ascii_fast_path && !utf8_in_place) {
      char *converted = NULL;
      size_t converted_size = 0;
      int rc = cimport_convert_named_to_utf8(ctx->file_data, ctx->file_size,
                                             source_encoding, &converted,
                                             &converted_size);
      if (rc) {
        return cimport_fail(ctx, rc);
      }
      ctx->converted_data = converted;
      ctx->converted_size = converted_size;
      ctx->file_data = converted;
      ctx->file_size = converted_size;
    }
    /* Native import ignores binary zeros, including inside fields. Remove
     * them after decoding (UTF-16/32 contain meaningful zero bytes before
     * decoding) so C string consumers cannot truncate the rest of a field.
     * Keep the mapped input read-only and avoid copying ordinary text. */
    const char *nul = ctx->file_size && (ctx->converted_data || bytes.nul)
        ? memchr(ctx->file_data, 0, ctx->file_size) : NULL;
    if (nul) {
      size_t prefix = (size_t)(nul - ctx->file_data);
      char *clean = ctx->converted_data;
      if (!clean) {
        clean = ctx->file_size < SIZE_MAX ? malloc(ctx->file_size + 1) : NULL;
        if (!clean) {
          return cimport_fail(ctx, 920);
        }
        memcpy(clean, ctx->file_data, prefix);
      }
      size_t used = prefix;
      for (size_t i = prefix; i < ctx->file_size; i++)
        if (ctx->file_data[i]) clean[used++] = ctx->file_data[i];
      clean[used] = 0;
      ctx->file_data = ctx->converted_data = clean;
      ctx->file_size = ctx->converted_size = used;
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

    /* Bytes read in place were scanned for a CR ending a line by itself. */
    int normalize_rc = ctx->converted_data || bytes.lone_cr ? cimport_normalize_line_endings(ctx) : 0;
    if (!normalize_rc)
      normalize_rc = cimport_normalize_delimiters(ctx);
    /* A parse in chunks checks the quoted rows of the rows it reads, in
     * parallel, and the bytes before them; otherwise the whole input is
     * checked here. */
    ctx->quoted_rows_in_parse = cimport_checks_quoted_rows(ctx) && ctx->file_size >= CIMPORT_MIN_CHUNK * 2 &&
                                cimport_parse_threads() > 1;
    if (!normalize_rc && !ctx->quoted_rows_in_parse)
      normalize_rc = cimport_check_quoted_rows(ctx, ctx->file_size);
    if (normalize_rc) {
      cimport_display_error("cimport: could not normalize delimiters\n");
      return cimport_fail(ctx, normalize_rc);
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

    {
        char text[32] = "";
        SF_macro_use("___cimport_startcol", text, sizeof(text) - 1);
        ctx->col_first = atoi(text) > 1 ? atoi(text) - 1 : 0;
        text[0] = 0;
        SF_macro_use("___cimport_endcol", text, sizeof(text) - 1);
        ctx->col_last = atoi(text) > 0 ? atoi(text) : 0;
    }
    if (cimport_parse_header(ctx) != 0) {
        cimport_display_error("Failed to parse header");
        return cimport_fail(ctx, 920);
    }
    /* Only empty lines: no variables, as for an empty file. */
    if (ctx->num_columns == 0) return ctx;
    if (cimport_resolve_col_force(ctx, 0, ctx->num_columns) != 0) {
        return cimport_fail(ctx, 920);
    }
    /* The single-pass number parser applies only to the default grammar. */
    ctx->plain_numbers = ctx->decimal_separator == '.' && ctx->group_separator == '\0' &&
                         !cimport_locale_active();

    /* Chunks follow the thread limit, so threads() and the default cap apply. */
    ctx->num_threads = cimport_parse_threads();

    t_start = ctools_timer_ms();

    int num_chunks = ctx->num_threads;

    if (ctx->file_size < CIMPORT_MIN_CHUNK * 2) {
        num_chunks = 1;
    }

    ctx->num_chunks = num_chunks;
    ctx->stream_rows = cimport_should_stream(ctx);
    ctx->compact_rows = cimport_cache_eligible(ctx) && !ctx->stream_rows;
    ctx->cache_budget.limit = cimport_cache_limit();
    atomic_init(&ctx->cache_budget.used, 0);

parse_rows:
    t_start = ctools_timer_ms();
    ctx->chunks = ctools_safe_calloc2(ctx->num_chunks, sizeof(CImportParsedChunk));
    if (!ctx->chunks) {
        return cimport_fail(ctx, 920);
    }

    if (num_chunks == 1 && !ctx->stream_rows && !ctx->compact_rows) {
        CImportParsedChunk *chunk = &ctx->chunks[0];
        ctools_arena_init(&chunk->arena, 0);
        chunk->capacity = 65536;
        chunk->rows = ctools_safe_malloc2(chunk->capacity, sizeof(CImportParsedRow *));
        chunk->num_rows = 0;

        chunk->num_col_stats = ctx->num_columns;
        chunk->col_stats = ctools_safe_calloc2(ctx->num_columns, sizeof(CImportColumnParseStats));

        if (!chunk->rows || !chunk->col_stats) {
            return cimport_fail(ctx, 920);
        }
        cimport_init_col_stats(chunk->col_stats, ctx->num_columns);

        /* 512 KiB of field refs (grown for wider rows): on the heap, off the
         * Stata thread's stack. */
        int field_capacity = CTOOLS_MAX_COLUMNS;
        CImportFieldRef *field_buf = ctools_safe_malloc2((size_t)field_capacity, sizeof(CImportFieldRef));
        if (!field_buf) {
            return cimport_fail(ctx, 920);
        }
        bool failed = false;
        const char *ptr = ctx->file_data;
        const char *end = ctx->file_data + ctx->file_size;

        /* Skip rows if requested (for varnames(N) where N > 1) */
        if (ctx->skip_rows > 0) {
            ptr = cimport_skip_rows(ptr, end, ctx->skip_rows, ctx->quote_char, ctx->bindquotes);
        }

        bool trailing_empty = cimport_trailing_empty_line(ctx, end);
        while (ptr < end || trailing_empty) {
            if (ptr >= end) {
                ptr = end - 1;
                trailing_empty = false;
            }
            const char *next_row;
            bool had_unmatched;
            int num_fields = cimport_parse_row_grow(ctx, ptr, end, &field_buf, &field_capacity,
                                                    &next_row, &had_unmatched);
            if (num_fields < 0) {
                if (num_fields == CIMPORT_FIELD_TOO_LONG) atomic_store(&ctx->error_code, 109);
                failed = true;
                break;
            }

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
                    cimport_record_unmatched_quote(ctx, (size_t)(ptr - ctx->file_data));
                }
            }

            if (chunk->num_rows >= chunk->capacity) {
                if (chunk->capacity > SIZE_MAX / 2) {
                    failed = true;
                    break;
                }
                size_t new_capacity = chunk->capacity * 2;
                if (new_capacity > SIZE_MAX / sizeof(CImportParsedRow *)) {
                    failed = true;
                    break;
                }
                CImportParsedRow **new_rows = realloc(chunk->rows, sizeof(CImportParsedRow *) * new_capacity);
                if (!new_rows) {
                    failed = true;
                    break;
                }
                chunk->rows = new_rows;
                chunk->capacity = new_capacity;
            }

            CImportParsedRow *row = ctools_arena_alloc(&chunk->arena, cimport_row_size(num_fields));
            if (!row || !cimport_fill_row(row, field_buf, num_fields)) {
                failed = true;
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
        free(field_buf);
        if (failed) {
            if (atomic_load(&ctx->error_code) == 109) cimport_display_error(CIMPORT_FIELD_TOO_LONG_MSG);
            return cimport_fail(ctx, 920);
        }
    } else {
        /* Start at the first row read, after the varnames(#) skip and the
         * empty rows emptylines(skip) drops, so chunk 0 always begins with
         * the header that every consumer skips as its row 0. */
        size_t start_offset = ctx->first_row_offset;
        size_t data_size = ctx->file_size - start_offset;
        size_t chunk_size = data_size / num_chunks;
        if (ctx->quoted_rows_in_parse && cimport_check_quoted_rows(ctx, start_offset))
            return cimport_fail(ctx, 5101);

        size_t *boundaries = ctools_safe_malloc2((size_t)num_chunks + 1, sizeof(size_t));
        if (!boundaries) {
            return cimport_fail(ctx, 920);
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
            pthread_attr_t attr;
            int have_attr = ctools_worker_attr_init(&attr) == 0;
            int threads_created = 0;
            for (int i = 0; i < num_chunks; i++) {
                if (pthread_create(&threads[threads_created], have_attr ? &attr : NULL,
                                   cimport_parse_chunk_parallel, &tasks[i]) == 0) {
                    threads_created++;
                } else {
                    cimport_parse_chunk_parallel(&tasks[i]);
                }
            }
            for (int i = 0; i < threads_created; i++) {
                pthread_join(threads[i], NULL);
            }
            if (have_attr) pthread_attr_destroy(&attr);
        }

        free(boundaries);

        if (atomic_load(&ctx->error_code) == 5101)
            cimport_display_error("Unmatched quote exceeded maxquotedrows()\n");
        if (atomic_load(&ctx->error_code) == 109)
            cimport_display_error(CIMPORT_FIELD_TOO_LONG_MSG);
        if (atomic_load(&ctx->error_code) != 0) {
            return cimport_fail(ctx, 920);
        }
    }

    t_end = ctools_timer_ms();
    ctx->time_parse += t_end - t_start;

    ctx->total_rows = 0;
    for (int c = 0; c < ctx->num_chunks; c++) {
        ctx->total_rows += ctx->chunks[c].num_rows;
    }
    /* Every consumer skips row 0 of chunk 0 as the header, so it must be there. */
    if (ctx->has_header && ctx->total_rows > 0 && ctx->chunks[0].num_rows == 0) {
        return cimport_fail(ctx, 459);
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

    /* Ragged rows adding columns need the general inference pass. */
    if ((ctx->stream_rows || ctx->compact_rows) && ctx->max_fields_seen > ctx->num_columns)
        goto retain_rows;

    if (ctx->max_fields_seen > ctx->num_columns) {
        int old_num = ctx->num_columns;
        int new_num = ctx->max_fields_seen;

        CImportColumnInfo *new_columns = realloc(ctx->columns, new_num * sizeof(CImportColumnInfo));
        if (!new_columns) {
            return cimport_fail(ctx, 920);
        }
        ctx->columns = new_columns;

        for (int i = old_num; i < new_num; i++) {
            memset(&ctx->columns[i], 0, sizeof(CImportColumnInfo));
            snprintf(ctx->columns[i].name, CTOOLS_MAX_VARNAME_LEN, "v%d", i + 1 - ctx->col_first);
            ctx->columns[i].type = CIMPORT_COL_UNKNOWN;
            ctx->columns[i].max_strlen = 0;
        }

        ctx->num_columns = new_num;
        if (cimport_resolve_col_force(ctx, old_num, new_num) != 0) {
            return cimport_fail(ctx, 920);
        }

        /* Chunks were parsed with the old column count, so stats for the
         * expanded columns were not collected: each chunk's rows are scanned
         * again for them, the chunks in parallel. */
        CImportExpandTask expand[CIMPORT_MAX_THREADS];
        int nexpand = ctx->num_chunks < CIMPORT_MAX_THREADS ? ctx->num_chunks : CIMPORT_MAX_THREADS;
        for (int c = 0; c < nexpand; c++) expand[c] = (CImportExpandTask){ctx, c, new_num, 0};
        ctools_persistent_pool *pool = nexpand > 1 ? ctools_get_global_pool() : NULL;
        if (pool == NULL || ctools_persistent_pool_submit_batch(pool, cimport_expand_task, expand, nexpand,
                                                                sizeof(CImportExpandTask)) != 0) {
            for (int c = 0; c < nexpand; c++) cimport_expand_task(&expand[c]);
        } else {
            ctools_persistent_pool_wait(pool);
        }
        for (int c = nexpand; c < ctx->num_chunks; c++) {
            CImportExpandTask task = {ctx, c, new_num, 0};
            cimport_expand_task(&task);
            if (task.rc) return cimport_fail(ctx, task.rc);
        }
        for (int c = 0; c < nexpand; c++)
            if (expand[c].rc) return cimport_fail(ctx, expand[c].rc);
    }

    t_start = ctools_timer_ms();
    cimport_select_rows(ctx);
    cimport_infer_column_types(ctx);
    if (atomic_load(&ctx->error_code)) return cimport_fail(ctx, 920);
    t_end = ctools_timer_ms();
    ctx->time_type_infer += t_end - t_start;
    /* The strL transfers read compact chunks block by block (see
     * CImportStrlBlock), but need rows when only the header was kept. */
    if (ctx->stream_rows) {
        for (int col = 0; col < ctx->num_columns; col++) {
            if (ctx->columns[col].use_strl) goto retain_rows;
        }
    }

    ctx->is_loaded = true;
    ctx->cache_ready = false;

    return ctx;

retain_rows:
    for (int c = 0; c < ctx->num_chunks; c++) {
        free(ctx->chunks[c].rows);
        free(ctx->chunks[c].col_stats);
        ctools_arena_free(&ctx->chunks[c].arena);
        cimport_cache_free(ctx->chunks[c].cache, &ctx->cache_budget);
    }
    free(ctx->chunks);
    ctx->chunks = NULL;
    ctx->stream_rows = false;
    ctx->compact_rows = false;
    ctx->num_unmatched_quote_warnings = 0;
    goto parse_rows;
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

                        if (!chunk_has_quotes || !field_has_quote) {
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

                    /* string_data holds exactly total_rows entries. */
                    if (row_idx >= ctx->total_rows) {atomic_store(&ctx->error_code,459);free(field_buf);return (void*)1;}
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
            pthread_attr_t attr;
            int have_attr = ctools_worker_attr_init(&attr) == 0;
            int threads_created = 0;
            for (int t = 0; t < num_threads; t++) {
                if (pthread_create(&threads[threads_created], have_attr ? &attr : NULL,
                                   cimport_build_cache_worker, &tasks[t]) == 0) {
                    threads_created++;
                } else {
                    cimport_build_cache_worker(&tasks[t]);
                }
            }
            for (int t = 0; t < threads_created; t++) {
                pthread_join(threads[t], NULL);
            }
            if (have_attr) pthread_attr_destroy(&attr);
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

/* Parse a column list from the ado's local <macro_name>. The ado also passes
 * its byte length in <macro_name>_len, so a list that cannot be read whole
 * fails with 198 instead of losing columns. *out_cols stays NULL when no list
 * was given (no or zero length). */
static ST_retcode cimport_parse_col_list(const char *macro_name, int **out_cols, int *out_count) {
    char name[64], length_text[32] = "";
    int length;
    *out_cols = NULL;
    *out_count = 0;

    snprintf(name, sizeof(name), "%s_len", macro_name);
    if (SF_macro_use(name, length_text, sizeof(length_text) - 1) != 0) return 198;
    if (!length_text[0]) return 0;
    if (!ctools_safe_atoi(length_text, &length) || length < 0 || length > INT_MAX - 2) return 198;
    if (length == 0) return 0;

    char *buf = calloc((size_t)length + 2, 1);
    if (!buf) return 920;
    if (SF_macro_use((char *)macro_name, buf, length + 1) != 0 || strlen(buf) != (size_t)length) {
        free(buf);
        return 198;
    }

    if (strcmp(buf, "_all") == 0) {
        free(buf);
        if (!(*out_cols = calloc(1, sizeof(int)))) return 920;
        *out_count = 1;
        return 0;
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

    if (count == 0) { free(buf); return 0; }

    int *cols = ctools_safe_malloc2(count, sizeof(int));
    if (!cols) { free(buf); return 920; }

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
            free(buf);
            return 198;  /* Invalid column number */
        }
        cols[idx] = (int)val;
        idx++;
        p = endptr;  /* Move to character after the parsed number */
    }

    free(buf);
    *out_cols = cols;
    *out_count = idx;
    return 0;
}

/* Read the numericcols() and stringcols() lists; on failure neither is kept. */
static ST_retcode cimport_parse_col_lists(int **force_num, int *num_force_num,
                                          int **force_str, int *num_force_str) {
    ST_retcode rc = cimport_parse_col_list("___cimport_numcols", force_num, num_force_num);
    if (!rc) rc = cimport_parse_col_list("___cimport_strcols", force_str, num_force_str);
    if (rc) {
        free(*force_num);
        *force_num = NULL;
        *num_force_num = 0;
        cimport_display_error("cimport: cannot read numericcols() or stringcols()\n");
    }
    return rc;
}

/* ============================================================================
 * SCAN Mode
 * ============================================================================ */

/* Headers must already be exported before releasing the source. Retain it if
 * any chunk needs a second pass, any text uses references, or final inference
 * requires recovering original spelling/precision from a numeric block. */
static void cimport_release_cached_source(CImportContext *ctx)
{
    if (!ctx->compact_rows || !ctx->file_data) return;
    for (int c = 0; c < ctx->num_chunks; c++) {
        CImportParsedChunk *chunk = &ctx->chunks[c];
        size_t header = c == 0 && ctx->has_header ? 1 : 0;
        if (chunk->cached_rows != chunk->num_rows - header) return;
        for (CImportCachedBlock *block = chunk->cache; block; block = block->next) {
            for (int col = 0; col < ctx->num_columns; col++) {
                int kind = block->columns[col].kind;
                if (ctx->columns[col].type == CIMPORT_COL_STRING) {
                    if (ctx->columns[col].use_strl || kind != CIMPORT_CACHE_INLINE_TEXT) return;
                } else if (kind >= CIMPORT_CACHE_TEXT ||
                           (kind == CIMPORT_CACHE_FLOAT &&
                            ctx->columns[col].num_subtype != CIMPORT_NUM_FLOAT)) return;
            }
        }
    }
    free(ctx->converted_data);
    ctx->converted_data = NULL;
    ctx->converted_size = 0;
    cimport_munmap_file(ctx);
    ctx->file_data = NULL;
}

static ST_retcode cimport_do_scan(const char *filename, char delimiter, bool has_header, int skip_rows, bool verbose, CImportBindQuotesMode bindquotes,
                                   CImportNumericTypeMode numeric_type_mode, char decimal_sep, char group_sep,
                                   CImportEmptyLinesMode emptylines_mode, int max_quoted_rows,
                                   CImportEncoding requested_encoding) {
    char macro_val[65536];

    cimport_clear_cached_context();

    /* Parse column type overrides BEFORE parsing CSV so inference runs once with overrides */
    int num_force_num = 0, num_force_str = 0, *force_num = NULL, *force_str = NULL;
    ST_retcode list_rc = cimport_parse_col_lists(&force_num, &num_force_num, &force_str, &num_force_str);
    if (list_rc) return list_rc;

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

    if (verbose && ctx->stream_rows)
        cimport_display_msg("(row cache exceeds memory budget; using bounded two-pass import)\n");
    if (verbose && ctx->compact_rows) {
        size_t cached = 0;
        for (int c = 0; c < ctx->num_chunks; c++) cached += ctx->chunks[c].cached_rows;
        char cache_msg[160];
        snprintf(cache_msg, sizeof(cache_msg), "(compact row cache: %.2f MiB; %zu of %zu rows cached)\n",
                 (double)atomic_load(&ctx->cache_budget.used) / 1048576.0, cached, ctx->total_rows);
        cimport_display_msg(cache_msg);
    }

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
    SF_macro_save("_cimport_hasheader", ctx->has_header ? "1" : "0");

    /* Only the colrange() columns become variables, numbered from 1. */
    int first = ctx->col_first, stop = cimport_col_stop(ctx);
    if (stop < first) stop = first;
    snprintf(macro_val, sizeof(macro_val), "%d", stop - first);
    SF_macro_save("_cimport_nvar", macro_val);

    /* The lists grow with the column count: names take their bytes plus a
     * space, and a string width at most 11 characters. */
    size_t list_cap = 16 + 12 * (size_t)(stop - first);
    for (int i = first; i < stop; i++) list_cap += strlen(ctx->columns[i].name) + 1;
    char *list = malloc(list_cap);
    if (!list) {
        cimport_free_context(g_cimport_ctx);
        g_cimport_ctx = NULL;
        return 920;
    }
    ST_retcode meta_rc = 0;
    char *p = list;
    char *end = list + list_cap - 1;
    for (int i = first; i < stop && p < end; i++) {
        if (i > first && p < end) *p++ = ' ';
        size_t name_len = strlen(ctx->columns[i].name);
        size_t space_left = (size_t)(end - p);
        if (name_len > space_left) name_len = space_left;
        memcpy(p, ctx->columns[i].name, name_len);
        p += name_len;
    }
    *p = '\0';
    if (!meta_rc) meta_rc = SF_macro_save("_cimport_varnames", list);

    p = list;
    for (int i = first; i < stop && p < end; i++) {
        if (i > first && p < end) *p++ = ' ';
        if (p < end) *p++ = (ctx->columns[i].type == CIMPORT_COL_STRING) ? '1' : '0';
    }
    *p = '\0';
    if (!meta_rc) meta_rc = SF_macro_save("_cimport_vartypes", list);

    p = list;
    for (int i = first; i < stop && p < end; i++) {
        if (i > first && p < end) *p++ = ' ';
        if (p < end) {
            if (ctx->columns[i].type == CIMPORT_COL_STRING) {
                *p++ = '0';
            } else {
                *p++ = '0' + (int)ctx->columns[i].num_subtype;
            }
        }
    }
    *p = '\0';
    if (!meta_rc) meta_rc = SF_macro_save("_cimport_numtypes", list);

    p = list;
    for (int i = first; i < stop && p < end; i++) {
        if (i > first && p < end) *p++ = ' ';
        int len = ctx->columns[i].type == CIMPORT_COL_STRING ? ctx->columns[i].max_strlen : 0;
        int written = snprintf(p, (size_t)(end - p + 1), "%d", len);
        if (written > 0 && p + written <= end) p += written;
    }
    *p = '\0';
    if (!meta_rc) meta_rc = SF_macro_save("_cimport_strlens", list);
    p=list;
    for(int i=first;i<stop && p<end;i++) {
        if(i>first && p<end)*p++=' ';
        if(p<end)*p++=ctx->columns[i].use_strl ? '1' : '0';
    }
    *p=0;
    if (!meta_rc) meta_rc = SF_macro_save("_cimport_longtypes", list);
    snprintf(macro_val, sizeof(macro_val), "%llu", (unsigned long long)ctx->strl_bytes);
    if (!meta_rc) meta_rc = SF_macro_save("_cimport_strlbytes", macro_val);
    free(list);
    if (meta_rc) {
        cimport_free_context(g_cimport_ctx);
        g_cimport_ctx = NULL;
        return meta_rc;
    }

    /* Transport original headings independently, including quotes and spaces.
     * The field refs (512 KiB, grown for wider rows) live on the heap, off the
     * Stata thread's stack. */
    CImportFieldRef *headings = NULL;
    int count = 0;
    if (ctx->has_header) {
        int capacity = CTOOLS_MAX_COLUMNS;
        headings = ctools_safe_malloc2((size_t)capacity, sizeof(CImportFieldRef));
        if (headings) count = cimport_parse_header_row(ctx, &headings, &capacity, NULL, NULL);
        if (!headings || count < 0) {
            free(headings);
            cimport_free_context(g_cimport_ctx);
            g_cimport_ctx = NULL;
            return 920;
        }
    }
    for (int i = first; i < stop; i++) {
        char key[64], value[CTOOLS_MAX_STRING_LEN + 1] = "";
        snprintf(key, sizeof(key), "_cimport_label_%d", i + 1 - first);
        if (i < count) cimport_extract_field_fast(ctx->file_data, &headings[i], value, sizeof(value), cimport_extract_quote(ctx));
        SF_macro_save(key, value);
    }
    free(headings);

    /* Display warnings for malformed data (like Stata) */
    cimport_display_warnings(ctx);

    /* Avoid overlapping the source with Stata's rollback copy and destination. */
    cimport_release_cached_source(ctx);

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
enum { CIMPORT_STORE_NUMERIC = 0, CIMPORT_STORE_STRING = 1, CIMPORT_STORE_STRL = 2,
       CIMPORT_STORE_SKIP = 3 /* outside colrange() */ };

/* Where the store workers put values when they do not store through the
 * SPI: a group of rows for the .dta writer (cimport fulldta). Arrays are per
 * file column, indexed by observation - obs0; a NULL array is not wanted.
 * Text is copied forward into one buffer large enough for the group's bytes,
 * since decoding never lengthens a value. strL values are decoded in full. */
typedef struct {
    double **num;
    ctools_dta_span **str;
    char *text;
    ST_int obs0;
} CImportSink;

typedef struct {
    CImportContext *ctx;
    const unsigned char *kind;       /* per column: CIMPORT_STORE_* */
    int chunk;
    size_t row_start, row_end;       /* indices into chunk->rows */
    ST_int obs_start;                /* Stata observation of row_start */
    ST_retcode error;
    CImportSink *sink;               /* NULL: store through the SPI */
    const CImportCachedBlock *block; /* compact, with a sink: only this block */
} CImportStoreTile;

static inline ST_retcode cimport_put_number(const CImportStoreTile *task, ST_IIID store, int col,
                                            ST_int obs, double value) {
    if (!task->sink) return store(col + 1 - task->ctx->col_first, obs, value);
    if (task->sink->num[col]) task->sink->num[col][obs - task->sink->obs0] = value;
    return 0;
}

/* A str# value, which SF_sstore would cut to the variable's width. */
static inline ST_retcode cimport_put_text(const CImportStoreTile *task, int col, ST_int obs,
                                          const char *text) {
    if (!task->sink) return SF_sstore(col + 1 - task->ctx->col_first, obs, (char *)text);
    CImportSink *sink = task->sink;
    if (!sink->str[col]) return 0;
    size_t n = strlen(text), width = (size_t)task->ctx->columns[col].max_strlen;
    if (width < 1) width = 1;
    if (n > width) n = width;
    memcpy(sink->text, text, n);
    sink->str[col][obs - sink->obs0] = (ctools_dta_span){sink->text, n, 0};
    sink->text += n;
    return 0;
}

/* A strL value, decoded in full as cimport strldta decodes it. */
static inline void cimport_put_strl(const CImportStoreTile *task, int col, ST_int obs,
                                    const CImportFieldRef *field, char quote) {
    CImportSink *sink = task->sink;
    if (!sink || !sink->str[col]) return;
    const char *data = task->ctx->file_data, *src = data + field->offset;
    int n = (int)CIMPORT_FIELD_LENGTH(*field);
    if (!(field->length & CIMPORT_FIELD_QUOTED_FLAG)) {
        while (n > 0 && (src[n - 1] == '\r' || src[n - 1] == '\n')) n--;
        if (n > 0) memcpy(sink->text, src, (size_t)n);
    } else {
        n = cimport_extract_field_fast(data, (CImportFieldRef *)field, sink->text, n + 1, quote);
    }
    if (n <= 0) return;
    sink->str[col][obs - sink->obs0] = (ctools_dta_span){sink->text, (size_t)n, 0};
    sink->text += n;
}

static inline void cimport_store_text(CImportContext *ctx, CImportFieldRef *field, char *buf,
                                      int capacity, char quote)
{
    if (!(field->length & CIMPORT_FIELD_QUOTED_FLAG))
    {
        const char *src = ctx->file_data + field->offset;
        int len = (int)CIMPORT_FIELD_LENGTH(*field);
        while (len > 0 && (src[len - 1] == '\r' || src[len - 1] == '\n'))
            len--;
        if (len >= capacity)
            len = capacity - 1;
        memcpy(buf, src, (size_t)len);
        buf[len] = '\0';
    }
    else
    {
        cimport_extract_field_fast(ctx->file_data, field, buf, capacity, quote);
    }
}

static void *cimport_store_tile_worker(void *arg) {
    CImportStoreTile *task = (CImportStoreTile *)arg;
    CImportContext *ctx = task->ctx;
    CImportParsedChunk *chunk = &ctx->chunks[task->chunk];
    const int ncols = ctx->num_columns;
    const double missing = SV_missval;
    const char quote = cimport_extract_quote(ctx);
    ST_IIID store = task->sink ? NULL : (_stata_)->safestore;
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
                if (kind == CIMPORT_STORE_NUMERIC) rc = cimport_put_number(task, store, col, obs, missing);
                else if (kind == CIMPORT_STORE_STRING) rc = cimport_put_text(task, col, obs, "");
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
                    if (!rc) rc = cimport_put_number(task, store, col, obs, val);
                } else if (kind == CIMPORT_STORE_STRING) {
                    cimport_store_text(ctx, &field, buf, sizeof(buf), quote);
                    rc = cimport_put_text(task, col, obs, buf);
                } else if (kind == CIMPORT_STORE_STRL) {
                    /* Through the SPI, strL columns are transferred by the ado. */
                    cimport_put_strl(task, col, obs, &field, quote);
                }
                field.offset += CIMPORT_FIELD_LENGTH(field) + 1;
            }
            if (rc) { task->error = rc; return (void *)1; }
        }
    }
    return NULL;
}

/* Reparse a disjoint byte range directly into its disjoint Stata observations.
 * Scratch space is O(columns) per worker, independent of the number of rows.
 * Use the same parser, quote decoder, numeric grammar and checked SPI stores as
 * the retained-row path. The first pass still infers from every data value. */
static void *cimport_store_stream_worker(void *arg)
{
    CImportStoreTile *task = (CImportStoreTile *)arg;
    CImportContext *ctx = task->ctx;
    CImportParsedChunk *chunk = &ctx->chunks[task->chunk];
    /* Fields past the columns are skipped: a wider row is an error below. */
    CImportFieldRef *fields = ctools_safe_malloc2((size_t)ctx->num_columns + 1, sizeof(*fields));
    task->error = 0;
    if (!fields)
    {
        task->error = 920;
        return NULL;
    }
    const char *ptr = ctx->file_data + (ctx->compact_rows ? chunk->uncached_start : chunk->file_start);
    const char *end = ctx->file_data + chunk->file_end;
    bool trailing_empty =
        task->chunk == ctx->num_chunks - 1 && cimport_trailing_empty_line(ctx, end);
    size_t r = 0;
    const double missing = SV_missval;
    const char quote = cimport_extract_quote(ctx);
    ST_IIID store = task->sink ? NULL : (_stata_)->safestore;
    char buf[CTOOLS_MAX_STRING_LEN + 1];

    while (ptr < end || trailing_empty)
    {
        if (ptr >= end)
        {
            ptr = end - 1;
            trailing_empty = false;
        }
        const char *next;
        int nfields = cimport_parse_row_fast(ptr, end, ctx->delimiter, ctx->quote_char, fields,
                                             ctx->num_columns + 1, ctx->file_data, ctx->bindquotes,
                                             &next, NULL);
        ptr = next;
        if (ctx->emptylines_mode == CIMPORT_EMPTYLINES_SKIP &&
            cimport_row_is_empty(ctx->file_data, fields, nfields))
            continue;
        if (r++ < task->row_start)
            continue; /* Header in chunk zero only. */
        if (r > task->row_end || nfields > ctx->num_columns)
        {
            task->error = 459;
            break;
        }
        ST_int obs = task->obs_start + (ST_int)(r - task->row_start - 1);
        for (int col = 0; col < ctx->num_columns; col++)
        {
            ST_retcode rc = 0;
            if (task->kind[col] == CIMPORT_STORE_NUMERIC)
            {
                double val = missing;
                if (col < nfields)
                {
                    CImportFieldRef *field = &fields[col];
                    bool plain =
                        ctx->plain_numbers && !(field->length & CIMPORT_FIELD_QUOTED_FLAG) &&
                        cimport_parse_plain_decimal(ctx->file_data + field->offset,
                                                    (int)CIMPORT_FIELD_LENGTH(*field), &val);
                    if (!plain && !cimport_context_number(ctx, field, &val))
                    {
                        rc = atomic_load(&ctx->error_code);
                        val = missing;
                    }
                }
                if (!rc)
                    rc = cimport_put_number(task, store, col, obs, val);
            }
            else if (task->kind[col] == CIMPORT_STORE_STRING)
            {
                buf[0] = '\0';
                if (col < nfields)
                    cimport_store_text(ctx, &fields[col], buf, sizeof(buf), quote);
                rc = cimport_put_text(task, col, obs, buf);
            }
            else if (task->kind[col] == CIMPORT_STORE_STRL && col < nfields)
            {
                cimport_put_strl(task, col, obs, &fields[col], quote);
            }
            if (rc)
            {
                task->error = rc;
                break;
            }
        }
        if (task->error)
            break;
    }
    if (!task->error && r != task->row_end)
        task->error = 459;
    free(fields);
    return NULL;
}

/* Most blocks need no parsing on load. A numeric column can become text, or a
 * provisional float can become double, after its block was cached. Only those
 * columns are recovered from the original bytes; the block boundary suffices
 * to find rows, so no per-row offset or field-length arrays are retained. */
static void *cimport_store_compact_worker(void *arg)
{
    CImportStoreTile *task = arg;
    CImportContext *ctx = task->ctx;
    CImportParsedChunk *chunk = &ctx->chunks[task->chunk];
    unsigned char *repair = malloc((size_t)ctx->num_columns);
    /* Only the columns are stored; the parser skips any further fields. */
    CImportFieldRef *fields = ctools_safe_malloc2((size_t)ctx->num_columns + 1, sizeof(*fields));
    task->error = 0;
    if (!repair || !fields) { free(repair); free(fields); task->error = 920; return NULL; }
    const double missing = SV_missval;
    const char quote = cimport_extract_quote(ctx);
    ST_IIID store = task->sink ? NULL : (_stata_)->safestore;
    char buf[CTOOLS_MAX_STRING_LEN + 1];
    size_t stored = 0;
    /* With a sink, only task->block (a group of the .dta writer). */
    for (const CImportCachedBlock *block = task->block ? task->block : chunk->cache;
         block && !task->error; block = task->block ? NULL : block->next) {
        bool parse = false;
        for (int col = 0; col < ctx->num_columns; col++) {
            int kind = block->columns[col].kind;
            repair[col] = task->kind[col] == CIMPORT_STORE_STRING ||
                          (task->kind[col] == CIMPORT_STORE_STRL && task->sink)
                ? kind < CIMPORT_CACHE_TEXT
                : task->kind[col] == CIMPORT_STORE_NUMERIC &&
                  (kind >= CIMPORT_CACHE_TEXT || (kind == CIMPORT_CACHE_FLOAT &&
                   ctx->columns[col].num_subtype != CIMPORT_NUM_FLOAT));
            parse |= repair[col] != 0;
        }
        const char *ptr = parse ? ctx->file_data + block->file_start : NULL;
        const char *end = parse ? ctx->file_data + block->file_end : NULL;
        for (size_t r = 0; r < block->rows && !task->error; r++, stored++) {
            int nfields = 0;
            if (parse) {
                do {
                    if (ptr >= end) {
                        if (!cimport_trailing_empty_line(ctx, end)) { task->error = 459; break; }
                        ptr = end - 1;
                    }
                    const char *next;
                    nfields = cimport_parse_row_fast(ptr, end, ctx->delimiter, ctx->quote_char,
                        fields, ctx->num_columns + 1, ctx->file_data, ctx->bindquotes, &next, NULL);
                    ptr = next;
                } while (ctx->emptylines_mode == CIMPORT_EMPTYLINES_SKIP &&
                         cimport_row_is_empty(ctx->file_data, fields, nfields));
                if (task->error) break;
            }
            ST_int obs = task->obs_start + (ST_int)stored;
            for (int col = 0; col < ctx->num_columns; col++) {
                ST_retcode rc = 0;
                if (task->kind[col] == CIMPORT_STORE_NUMERIC) {
                    double value = missing;
                    if (!repair[col]) value = cimport_cache_number(block, col, r, missing);
                    else if (col < nfields && !cimport_context_number(ctx, &fields[col], &value)) {
                        rc = atomic_load(&ctx->error_code);
                        value = missing;
                    }
                    if (!rc) rc = cimport_put_number(task, store, col, obs, value);
                } else if (task->kind[col] == CIMPORT_STORE_STRING) {
                    if (!repair[col] && block->columns[col].kind == CIMPORT_CACHE_INLINE_TEXT) {
                        rc = cimport_put_text(task, col, obs, cimport_cache_text(block, col, r));
                        if (rc) { task->error = rc; break; }
                        continue;
                    }
                    buf[0] = '\0';
                    if (!repair[col]) {
                        CImportFieldRef field = cimport_cache_field(block, col, r);
                        cimport_store_text(ctx, &field, buf, sizeof(buf), quote);
                    } else if (col < nfields) cimport_store_text(ctx, &fields[col], buf, sizeof(buf), quote);
                    rc = cimport_put_text(task, col, obs, buf);
                } else if (task->kind[col] == CIMPORT_STORE_STRL && task->sink) {
                    if (!repair[col] && block->columns[col].kind == CIMPORT_CACHE_INLINE_TEXT) {
                        /* Inline text is already decoded, and short. */
                        const char *text = cimport_cache_text(block, col, r);
                        size_t n = strlen(text);
                        if (n && task->sink->str[col]) {
                            memcpy(task->sink->text, text, n);
                            task->sink->str[col][obs - task->sink->obs0] = (ctools_dta_span){task->sink->text, n, 0};
                            task->sink->text += n;
                        }
                    } else if (!repair[col]) {
                        CImportFieldRef field = cimport_cache_field(block, col, r);
                        cimport_put_strl(task, col, obs, &field, quote);
                    } else if (col < nfields) {
                        cimport_put_strl(task, col, obs, &fields[col], quote);
                    }
                }
                if (rc) { task->error = rc; break; }
            }
        }
    }
    free(repair);
    free(fields);
    if (task->block) {
        if (!task->error && stored != task->block->rows) task->error = 459;
        return NULL;
    }
    size_t total = task->row_end - task->row_start;
    if (!task->error && (stored != chunk->cached_rows || stored > total)) task->error = 459;
    if (!task->error && stored < total) {
        CImportStoreTile suffix = *task;
        suffix.obs_start += (ST_int)stored;
        suffix.row_start = stored ? 0 : task->row_start;
        suffix.row_end = total - stored + suffix.row_start;
        cimport_store_stream_worker(&suffix);
        task->error = suffix.error;
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
        int num_force_num = 0, num_force_str = 0, *force_num = NULL, *force_str = NULL;
        ST_retcode list_rc = cimport_parse_col_lists(&force_num, &num_force_num, &force_str, &num_force_str);
        if (list_rc) return list_rc;

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
        cimport_release_cached_source(ctx);
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
    /* A scanned context is healthy (a failed scan frees it), so a code left
     * by an earlier route on this scan, such as an allocation failure in
     * cimport fulldta that the ado fell back from, is stale here. */
    atomic_store(&ctx->error_code, 0);

    double t_load_start = ctools_timer_ms();

    ST_int stata_nvars = SF_nvars();
    ST_int stata_nobs = SF_nobs();

    (void)verbose; /* verbose timing reported via scalars */

    int selected = cimport_col_stop(ctx) - ctx->col_first;
    if (selected < 0) selected = 0;
    if (stata_nvars < selected) {
        snprintf(msg, sizeof(msg), "Error: Stata has %d variables but CSV has %d columns\n",
                 stata_nvars, selected);
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

    /* Parallel SPI store from typed blocks or retained rows. Referenced strings
     * and uncached suffixes are decoded directly from the mapped source. */
    {
        unsigned char *kind = malloc((size_t)ctx->num_columns);
        size_t max_tiles = (ctx->stream_rows || ctx->compact_rows) ? (size_t)ctx->num_chunks
            : ctx->num_chunks + ctx->total_rows / CIMPORT_STORE_TILE_ROWS + 1;
        CImportStoreTile *tiles = ctools_safe_malloc2(max_tiles, sizeof(CImportStoreTile));
        if (!kind || !tiles) {
            free(kind); free(tiles);
            cimport_clear_cached_context();
            return 459;
        }

        ST_retcode store_error = 0;
        for (int col = 0; col < ctx->num_columns; col++) {
            ST_int var = col + 1 - ctx->col_first;
            if (!cimport_col_selected(ctx, col)) {
                kind[col] = CIMPORT_STORE_SKIP;
                continue;
            }
            bool is_string = ctx->columns[col].type == CIMPORT_COL_STRING;
            kind[col] = !is_string ? CIMPORT_STORE_NUMERIC
                      : SF_var_is_strl(var) ? CIMPORT_STORE_STRL : CIMPORT_STORE_STRING;
            if ((SF_var_is_string(var) != 0) != is_string) store_error = 198;
        }

        size_t ntiles = 0, obs = 1;
        for (int c = 0; c < ctx->num_chunks && !store_error; c++) {
            CImportParsedChunk *chunk = &ctx->chunks[c];
            size_t r = (c == 0 && ctx->has_header) ? 1 : 0;
            if (ctx->stream_rows || ctx->compact_rows) {
                tiles[ntiles++] = (CImportStoreTile){ctx, kind, c, r, chunk->num_rows, (ST_int)obs, 0, NULL, NULL};
                obs += chunk->num_rows - r;
                continue;
            }
            for (; r < chunk->num_rows; r += CIMPORT_STORE_TILE_ROWS) {
                size_t end = r + CIMPORT_STORE_TILE_ROWS;
                if (end > chunk->num_rows) end = chunk->num_rows;
                tiles[ntiles] = (CImportStoreTile){ctx, kind, c, r, end, (ST_int)obs, 0, NULL, NULL};
                obs += end - r;
                ntiles++;
            }
        }
        if (!store_error && obs - 1 != ctx->total_rows) store_error = 459;

        if (!store_error) {
            ctools_thread_func worker = ctx->compact_rows ? cimport_store_compact_worker
                : ctx->stream_rows ? cimport_store_stream_worker : cimport_store_tile_worker;
            ctools_persistent_pool *pool = ntiles > 1 ? ctools_get_global_pool() : NULL;
            if (pool == NULL || ctools_persistent_pool_submit_batch(pool, worker,
                                    tiles, ntiles, sizeof(CImportStoreTile)) != 0) {
                for (size_t t = 0; t < ntiles; t++) worker(&tiles[t]);
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
    SF_scal_save("_cimport_time_encoding", ctx->time_encoding / 1000.0);
    SF_scal_save("_cimport_time_parse", ctx->time_parse / 1000.0);
    SF_scal_save("_cimport_time_infer", ctx->time_type_infer / 1000.0);
    SF_scal_save("_cimport_time_cache", ctx->time_cache / 1000.0);
    SF_scal_save("_cimport_time_store", time_spi / 1000.0);
    SF_scal_save("_cimport_time_total", (ctx->time_mmap + ctx->time_encoding + ctx->time_parse + ctx->time_type_infer + ctx->time_cache + time_spi) / 1000.0);
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
    if (ctx) col+=ctx->col_first; /* variable number -> file column (colrange()) */
    if (!ctx || col<1 || !cimport_col_selected(ctx,col-1) || row<1 || (size_t)row>ctx->total_rows ||
        ctx->columns[col-1].type!=CIMPORT_COL_STRING) return 198;
    /* The column cache reads retained rows, which compact chunks lack. */
    if (ctx->compact_rows || ctx->stream_rows) return 459;
    cimport_build_column_cache(ctx);
    if(!ctx->cache_ready)return atomic_load(&ctx->error_code) ? atomic_load(&ctx->error_code) : 920;
    /* Only strL columns keep their text in the column cache. */
    if(!ctx->col_cache || !ctx->col_cache[col-1].string_data || !ctx->col_cache[col-1].string_data[row-1])
        return 198;
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

/* strL transfer. The plugin API cannot write strL (SF_sstore on a strL
 * variable returns 0 and stores nothing), so _cimport_read_blobs in
 * cimport.ado stores them a block of rows at a time, one st_sstore() per
 * column. When the parsed bytes are the file's own (no encoding, delimiter,
 * line-ending or binary-zero rewrite), a value found verbatim in the file
 * travels only as its position: Mata reads the block's byte range of the
 * source once and cuts every such value from it with one substr(). Values
 * that decoding changed (doubled quotes, for example) travel as bytes.
 *
 * File, integers and doubles in host order:
 *   [uint32 version][uint32 columns][uint32 direct][double source bytes]
 *   [uint32 path length][path]
 * then blocks of up to CIMPORT_BLOB_BLOCK_ROWS rows, never spanning chunks:
 *   [uint32 rows b][double region start][double region length]
 *   and per column, in the order of __cimport_longcolumns, [uint32 mode]:
 *     mode 2: every value is empty; nothing follows (a new strL is empty).
 *     mode 0 or 1, plus 4 when at least half the values are empty: then
 *       [uint32 count n][n uint32 rows in block] list the nonempty rows, and
 *       only they follow; otherwise all b rows (n = b) do.
 *     mode 0: [n uint32 lengths][their bytes]
 *     mode 1: [n uint32 positions, 1-based in the region][n uint32 lengths]
 *             [uint32 patches p][p uint32 indices into the n rows]
 *             [p uint32 lengths][their bytes]; a patched row has length 0 in
 *             the first list. */
#define CIMPORT_BLOB_VERSION 3
/* Both limits can be set at build time, as CIMPORT_STREAM_THRESHOLD and
 * CIMPORT_CACHE_BUDGET can, so that a small file exercises the cuts. */
#ifndef CIMPORT_BLOB_BLOCK_ROWS
#define CIMPORT_BLOB_BLOCK_ROWS 65536
#endif
/* A block ends early once its rows span this many source bytes, which also
 * bounds its region and its byte payload (decoding never lengthens text). */
#ifndef CIMPORT_BLOB_BLOCK_BYTES
#define CIMPORT_BLOB_BLOCK_BYTES ((uint64_t)64 << 20)
#endif

enum { CIMPORT_BLOB_EMPTY = 0, CIMPORT_BLOB_DIRECT, CIMPORT_BLOB_BYTES };

/* Source bytes a row spans from its start through its last transferred field. */
static uint64_t cimport_blob_row_end(const CImportParsedRow *row, int max_col) {
    uint64_t end = row->offset;
    int stop = row->num_fields <= max_col ? row->num_fields : max_col + 1;
    for (int f = 0; f < stop; f++) end += (row->lengths[f] & 0x7FFFFFFFu) + 1;
    return end;
}

/* A block of rows of one chunk for the strL transfers; obs0 is the 1-based
 * observation of its first row. Retained rows are chunk->rows[r0, r1). A
 * compact chunk keeps no rows: a block is one of its cached blocks, whose
 * field references serve when every strL column was cached as text, or else
 * its bytes [start, end), parsed again, of which the first skip rows are
 * not data (the header); its rows are then numbered from r0 = 0. */
typedef struct {
    int chunk;
    size_t r0, r1;
    uint64_t obs0;
    uint64_t span; /* source bytes its rows span through the last strL column */
    bool compact;
    const CImportCachedBlock *cached;
    uint64_t start, end;
    size_t skip;
} CImportStrlBlock;

typedef struct {
    int *cols;        /* file columns, 0-based, in transfer order */
    int ncols;
    int max_col;
    bool sorted;
    size_t scratch_cap;
    CImportStrlBlock *blocks;
    size_t nblocks;
} CImportStrlPlan;

static void cimport_strl_plan_free(CImportStrlPlan *plan) {
    free(plan->cols);
    free(plan->blocks);
    memset(plan, 0, sizeof(*plan));
}

static ST_retcode cimport_strl_add(CImportStrlPlan *plan, size_t *cap, CImportStrlBlock block) {
    if (plan->nblocks == *cap) {
        size_t grown_cap = *cap ? 2 * *cap : 64;
        CImportStrlBlock *grown = realloc(plan->blocks, grown_cap * sizeof(*grown));
        if (!grown) return 920;
        plan->blocks = grown;
        *cap = grown_cap;
    }
    plan->blocks[plan->nblocks++] = block;
    return 0;
}

static ST_retcode cimport_strl_split(const CImportContext *ctx, CImportStrlPlan *plan, size_t *cap,
                                     int chunk, uint64_t start, uint64_t end, size_t skip,
                                     size_t nrows, size_t *rows);

/* The strL columns named in local __cimport_longcolumns (variable numbers,
 * byte length in __cimport_longcolumns_len), and the blocks of rows they are
 * transferred in: up to CIMPORT_BLOB_BLOCK_ROWS rows, cut early on the source
 * bytes the rows span, never across chunks. */
static ST_retcode cimport_strl_plan(CImportContext *ctx, CImportStrlPlan *plan) {
    memset(plan, 0, sizeof(*plan));
    ST_retcode rc = cimport_parse_col_list("___cimport_longcolumns", &plan->cols, &plan->ncols);
    if (rc) return rc;
    if (!plan->ncols) rc = 198;
    plan->sorted = true;
    plan->scratch_cap = 1;
    for (int k = 0; k < plan->ncols && !rc; k++) {
        int col = plan->cols[k] + ctx->col_first; /* variable number -> file column (colrange()) */
        if (col < 1 || col > ctx->num_columns || !cimport_col_selected(ctx, col - 1) ||
            ctx->columns[col - 1].type != CIMPORT_COL_STRING) { rc = 198; break; }
        plan->cols[k] = col - 1;
        if (k && plan->cols[k] <= plan->cols[k - 1]) plan->sorted = false;
        if (plan->cols[k] > plan->max_col) plan->max_col = plan->cols[k];
        if ((size_t)ctx->columns[col - 1].max_strlen + 1 > plan->scratch_cap)
            plan->scratch_cap = (size_t)ctx->columns[col - 1].max_strlen + 1;
    }
    if (!rc && plan->scratch_cap > INT_MAX) rc = 109;
    if (!rc && ctx->total_rows && (!ctx->file_data || !ctx->chunks)) rc = 459;

    size_t cap = 0, rows = 0;
#define CIMPORT_STRL_ADD(...) do { \
        if (!rc) rc = cimport_strl_add(plan, &cap, (CImportStrlBlock){__VA_ARGS__}); \
    } while (0)
    /* A compact chunk's cached blocks and uncached suffix become blocks as
     * they are while within the transfer's limits; longer ones are cut. */
    for (int c = 0; c < ctx->num_chunks && !rc && ctx->compact_rows; c++) {
        CImportParsedChunk *chunk = &ctx->chunks[c];
        size_t header = (c == 0 && ctx->has_header) ? 1 : 0;
        for (const CImportCachedBlock *block = chunk->cache; block && !rc; block = block->next) {
            uint64_t span = block->file_end - block->file_start;
            if (block->rows > CIMPORT_BLOB_BLOCK_ROWS || span > CIMPORT_BLOB_BLOCK_BYTES) {
                rc = cimport_strl_split(ctx, plan, &cap, c, block->file_start, block->file_end, 0,
                                        block->rows, &rows);
                continue;
            }
            bool refs = true;
            for (int k = 0; k < plan->ncols && refs; k++) {
                int kind = block->columns[plan->cols[k]].kind;
                refs = kind == CIMPORT_CACHE_TEXT || kind == CIMPORT_CACHE_FIXED_TEXT;
            }
            CIMPORT_STRL_ADD(c, 0, block->rows, (uint64_t)rows + 1, span,
                             true, refs ? block : NULL, block->file_start, block->file_end, 0);
            rows += block->rows;
        }
        size_t rest = chunk->num_rows - header - chunk->cached_rows;
        if (!rc && rest) {
            uint64_t span = chunk->file_end - chunk->uncached_start;
            size_t skip = chunk->cached_rows ? 0 : header;
            if (rest > CIMPORT_BLOB_BLOCK_ROWS || span > CIMPORT_BLOB_BLOCK_BYTES) {
                rc = cimport_strl_split(ctx, plan, &cap, c, chunk->uncached_start, chunk->file_end,
                                        skip, rest, &rows);
            } else {
                CIMPORT_STRL_ADD(c, 0, rest, (uint64_t)rows + 1, span,
                                 true, NULL, chunk->uncached_start, chunk->file_end, skip);
                rows += rest;
            }
        }
    }
    for (int c = 0; c < ctx->num_chunks && !rc && !ctx->compact_rows; c++) {
        CImportParsedChunk *chunk = &ctx->chunks[c];
        size_t r = (c == 0 && ctx->has_header) ? 1 : 0;
        while (r < chunk->num_rows) {
            size_t r0 = r++;
            uint64_t first = chunk->rows[r0]->offset;
            uint64_t span = cimport_blob_row_end(chunk->rows[r0], plan->max_col) - first, next;
            while (r < chunk->num_rows && r - r0 < CIMPORT_BLOB_BLOCK_ROWS &&
                   (next = cimport_blob_row_end(chunk->rows[r], plan->max_col) - first) <= CIMPORT_BLOB_BLOCK_BYTES) {
                span = next;
                r++;
            }
            CIMPORT_STRL_ADD(c, r0, r, (uint64_t)rows + 1, span, false, NULL, 0, 0, 0);
            if (rc) break;
            rows += r - r0;
        }
    }
#undef CIMPORT_STRL_ADD
    if (!rc && rows != ctx->total_rows) rc = 459;
    if (rc) cimport_strl_plan_free(plan);
    return rc;
}

/* The rows of a block, one at a time: the fields of the plan's columns, and
 * how many fields the row has (a column past them is absent). */
typedef struct {
    const CImportContext *ctx;
    const CImportStrlPlan *plan;
    const CImportStrlBlock *block;
    CImportFieldRef *parsed; /* a parsed row's fields */
    const char *ptr, *end;
    size_t i;
} CImportStrlCursor;

static bool cimport_strl_parse_row(CImportStrlCursor *cur, int *nfields) {
    const CImportContext *ctx = cur->ctx;
    do {
        if (cur->ptr >= cur->end) {
            if (!cimport_trailing_empty_line(ctx, cur->end)) return false;
            cur->ptr = cur->end - 1;
        }
        const char *next;
        *nfields = cimport_parse_row_fast(cur->ptr, cur->end, ctx->delimiter, ctx->quote_char,
                                          cur->parsed, ctx->num_columns + 1, ctx->file_data,
                                          ctx->bindquotes, &next, NULL);
        cur->ptr = next;
    } while (ctx->emptylines_mode == CIMPORT_EMPTYLINES_SKIP &&
             cimport_row_is_empty(ctx->file_data, cur->parsed, *nfields));
    return true;
}

static bool cimport_strl_cursor_init(CImportStrlCursor *cur, const CImportContext *ctx,
                                     const CImportStrlPlan *plan, const CImportStrlBlock *block) {
    *cur = (CImportStrlCursor){ctx, plan, block, NULL, NULL, NULL, 0};
    if (!block->compact || block->cached) return true;
    cur->parsed = ctools_safe_malloc2((size_t)ctx->num_columns + 1, sizeof(*cur->parsed));
    if (!cur->parsed) return false;
    cur->ptr = ctx->file_data + block->start;
    cur->end = ctx->file_data + block->end;
    for (size_t s = 0; s < block->skip; s++) {
        int n;
        if (!cimport_strl_parse_row(cur, &n)) return false;
    }
    return true;
}

static bool cimport_strl_next(CImportStrlCursor *cur, CImportFieldRef *fields, int *nfields) {
    const CImportStrlPlan *plan = cur->plan;
    const CImportStrlBlock *block = cur->block;
    const int ncols = plan->ncols;
    size_t i = cur->i++;
    if (!block->compact) {
        const CImportParsedRow *row = cur->ctx->chunks[block->chunk].rows[block->r0 + i];
        *nfields = row->num_fields;
        if (plan->sorted) {
            CImportFieldRef field = {row->offset, 0};
            int k = 0;
            for (int f = 0; f < row->num_fields && f <= plan->max_col && k < ncols; f++) {
                field.length = row->lengths[f];
                if (f == plan->cols[k]) fields[k++] = field;
                field.offset += CIMPORT_FIELD_LENGTH(field) + 1;
            }
        } else {
            for (int k = 0; k < ncols; k++)
                if (plan->cols[k] < row->num_fields) fields[k] = cimport_row_field(row, plan->cols[k]);
        }
        return true;
    }
    if (block->cached) {
        /* A field past the row's end was cached as empty. */
        *nfields = INT_MAX;
        for (int k = 0; k < ncols; k++) fields[k] = cimport_cache_field(block->cached, plan->cols[k], i);
        return true;
    }
    if (!cimport_strl_parse_row(cur, nfields)) return false;
    for (int k = 0; k < ncols; k++)
        if (plan->cols[k] < *nfields) fields[k] = cur->parsed[plan->cols[k]];
    return true;
}

/* Blocks for nrows rows of a compact chunk's bytes [start, end), after its
 * first skip rows (the header), within CIMPORT_BLOB_BLOCK_ROWS rows and about
 * CIMPORT_BLOB_BLOCK_BYTES source bytes each. A cached block or uncached
 * suffix can hold far more (a suffix is unbounded once the cache budget runs
 * out), and the transfers hold a block's region, values and text at once.
 * The rows are parsed here with the strL cursor, exactly as the transfer
 * parses each block again from its own byte range, to find the cuts. */
static ST_retcode cimport_strl_split(const CImportContext *ctx, CImportStrlPlan *plan, size_t *cap,
                                     int chunk, uint64_t start, uint64_t end, size_t skip,
                                     size_t nrows, size_t *rows) {
    CImportStrlBlock whole = {chunk, 0, nrows, 0, end - start, true, NULL, start, end, skip};
    CImportStrlCursor cur;
    ST_retcode rc = cimport_strl_cursor_init(&cur, ctx, plan, &whole) ? 0 : (cur.parsed ? 459 : 920);
    for (size_t done = 0; !rc && done < nrows;) {
        uint64_t first = (uint64_t)(cur.ptr - ctx->file_data);
        size_t b = 0;
        int n;
        while (done + b < nrows && b < CIMPORT_BLOB_BLOCK_ROWS) {
            if (!cimport_strl_parse_row(&cur, &n)) { rc = 459; break; }
            b++;
            if ((uint64_t)(cur.ptr - ctx->file_data) - first >= CIMPORT_BLOB_BLOCK_BYTES) break;
        }
        if (rc) break;
        /* The last block keeps the original end, as the trailing empty line
         * of a file is found from it. */
        uint64_t last = done + b < nrows ? (uint64_t)(cur.ptr - ctx->file_data) : end;
        rc = cimport_strl_add(plan, cap, (CImportStrlBlock){chunk, 0, b, (uint64_t)*rows + 1, last - first,
                                                            true, NULL, first, last, 0});
        *rows += b;
        done += b;
    }
    free(cur.parsed);
    return rc;
}

typedef struct {
    const char *ptr; /* value bytes, or NULL when they are in the patch buffer */
    size_t off;      /* offset in the patch buffer when ptr is NULL */
    uint32_t len;
    unsigned char kind;
    unsigned char t; /* in memory, DIRECT: CIMPORT_BLOB_UNQUOTE/_CRLF applied to ptr[0, len) */
} CImportBlobValue;

/* Transforms Mata applies to source bytes (subinstr(), left to right): ""
 * to ", and CR LF to CR. */
enum { CIMPORT_BLOB_UNQUOTE = 1, CIMPORT_BLOB_CRLF = 2, CIMPORT_BLOB_TEXT = 8 };

/* Whether transform t of src[0, len) gives out[0, n). */
static bool cimport_blob_transform_is(const char *src, size_t len, int t, const char *out, size_t n) {
    size_t i = 0, j = 0;
    while (i < len) {
        char c = src[i];
        if ((t & CIMPORT_BLOB_UNQUOTE) && c == '"' && i + 1 < len && src[i + 1] == '"') i += 2;
        else if ((t & CIMPORT_BLOB_CRLF) && c == '\r' && i + 1 < len && src[i + 1] == '\n') i += 2;
        else i++;
        if (j == n || out[j] != c) return false;
        j++;
    }
    return j == n;
}

/* A block of the in-memory strL transfer (cimport blobmem): the block's
 * values, column by column, as (row in block, position, length, kind): every
 * row of a column more than half full, in order, else only its nonempty
 * ones. A position is 1-based in the block's region of the source, or for
 * kind CIMPORT_BLOB_TEXT in text, the values Mata cannot cut. */
typedef struct {
    uint32_t rows;
    uint64_t start, length; /* the region, in source-file bytes */
    uint32_t *counts;       /* per column */
    int32_t *r, *p, *n;
    unsigned char *t;
    size_t nent;
    char *text;
    size_t ntext;
} CImportBlobMem;

static void cimport_blobmem_free(CImportBlobMem *m) {
    free(m->counts);
    free(m->r);
    free(m->p);
    free(m->n);
    free(m->t);
    free(m->text);
    memset(m, 0, sizeof(*m));
}

typedef struct {
    CImportContext *ctx;
    const int *cols;  /* file columns, 0-based, in transfer order */
    int ncols;
    int max_col;      /* largest of cols */
    bool sorted;      /* cols ascending: one walk over each row's lengths */
    bool direct;      /* file_data[i] is byte base + i of the source file */
    uint64_t base;
    size_t scratch_cap;
    int chunk;
    size_t r0, r1;    /* the block's rows */
    char *out;        /* serialized block */
    size_t out_len;
    ST_retcode rc;
    const CImportStrlPlan *plan;
    const CImportStrlBlock *block;
    CImportBlobMem *mem;  /* non-NULL: build the block in memory, not serialized */
} CImportBlobTask;

static inline void cimport_blob_put(char **p, const void *src, size_t n) {
    memcpy(*p, src, n);
    *p += n;
}

static void *cimport_blob_worker(void *arg) {
    CImportBlobTask *task = arg;
    CImportContext *ctx = task->ctx;
    CImportParsedChunk *chunk = &ctx->chunks[task->chunk];
    const int ncols = task->ncols;
    const size_t b = task->r1 - task->r0;
    const char quote = cimport_extract_quote(ctx);
    const char *data = ctx->file_data;
    CImportBlobValue *values = ctools_safe_malloc2(b * (size_t)ncols, sizeof(*values));
    CImportFieldRef *fields = ctools_safe_malloc2((size_t)ncols, sizeof(*fields));
    bool *chunk_quotes = malloc((size_t)ncols);
    char *scratch = malloc(task->scratch_cap);
    char *patch = NULL;
    size_t patch_len = 0, patch_cap = 0;
    CImportStrlCursor cursor;
    bool cursor_ok = cimport_strl_cursor_init(&cursor, ctx, task->plan, task->block);
    task->out = NULL;
    task->out_len = 0;
    task->rc = 0;
    if (!values || !fields || !chunk_quotes || !scratch || !cursor_ok) { task->rc = 920; goto done; }
    for (int k = 0; k < ncols; k++) {
        int col = task->cols[k];
        chunk_quotes[k] = ctx->columns[col].column_has_quotes && chunk->col_stats &&
                          col < chunk->num_col_stats && chunk->col_stats[col].has_quotes;
    }

    /* Decode each value once, as the column cache did, and note whether the
     * result is a run of the source bytes. */
    uint64_t lo = UINT64_MAX, hi = 0;
    for (size_t i = 0; i < b; i++) {
        int nfields;
        if (!cimport_strl_next(&cursor, fields, &nfields)) { task->rc = 459; goto done; }
        for (int k = 0; k < ncols; k++) {
            CImportBlobValue *v = &values[i * (size_t)ncols + (size_t)k];
            *v = (CImportBlobValue){NULL, 0, 0, CIMPORT_BLOB_EMPTY, 0};
            if (task->cols[k] >= nfields) continue;
            CImportFieldRef *field = &fields[k];
            const char *src = data + field->offset;
            int raw = (int)CIMPORT_FIELD_LENGTH(*field), n;
            if (!chunk_quotes[k] || !(field->length & CIMPORT_FIELD_QUOTED_FLAG)) {
                n = raw;
                while (n > 0 && (src[n - 1] == '\r' || src[n - 1] == '\n')) n--;
                if (!n) continue;
                *v = (CImportBlobValue){src, 0, (uint32_t)n,
                                        task->direct ? CIMPORT_BLOB_DIRECT : CIMPORT_BLOB_BYTES, 0};
            } else {
                n = cimport_extract_field_fast(data, field, scratch, (int)task->scratch_cap, quote);
                if (n <= 0) continue;
                /* A quoted value without doubled quotes is the text inside them. */
                const char *at = NULL;
                if (task->direct) {
                    if (raw >= n + 1 && !memcmp(src + 1, scratch, (size_t)n)) at = src + 1;
                    else if (raw >= n && !memcmp(src, scratch, (size_t)n)) at = src;
                }
                /* In memory, Mata may also cut source bytes and undo doubled
                 * quotes or CR LF in them: the field without its quotes, or
                 * with only the first, or whole, before or after line ends. */
                int t = 0;
                size_t at_len = (size_t)n;
                if (!at && task->direct && task->mem) {
                    size_t ends[2] = {(size_t)raw, (size_t)raw};
                    while (ends[1] > 0 && (src[ends[1] - 1] == '\r' || src[ends[1] - 1] == '\n')) ends[1]--;
                    for (int e = 0; e < 2 && !at; e++) {
                        size_t end = ends[e];
                        const char *from[3] = {src + 1, src + 1, src};
                        size_t len[3] = {end >= 2 ? end - 2 : 0, end >= 1 ? end - 1 : 0, end};
                        bool ok[3] = {end >= 2 && src[0] == '"' && src[end - 1] == '"', end >= 1 && src[0] == '"', true};
                        for (int c = 0; c < 3 && !at; c++)
                            for (int tt = 1; tt <= 3 && ok[c] && !at; tt++)
                                if (cimport_blob_transform_is(from[c], len[c], tt, scratch, (size_t)n)) {
                                    at = from[c];
                                    at_len = len[c];
                                    t = tt;
                                }
                    }
                }
                if (at) {
                    *v = (CImportBlobValue){at, 0, (uint32_t)at_len, CIMPORT_BLOB_DIRECT, (unsigned char)t};
                } else {
                    if (patch_len + (size_t)n > patch_cap) {
                        size_t cap = patch_cap ? patch_cap : 65536;
                        while (cap < patch_len + (size_t)n) cap *= 2;
                        char *next = realloc(patch, cap);
                        if (!next) { task->rc = 920; goto done; }
                        patch = next;
                        patch_cap = cap;
                    }
                    memcpy(patch + patch_len, scratch, (size_t)n);
                    *v = (CImportBlobValue){NULL, patch_len, (uint32_t)n, CIMPORT_BLOB_BYTES, 0};
                    patch_len += (size_t)n;
                }
            }
            if (v->kind == CIMPORT_BLOB_DIRECT) {
                uint64_t start = (uint64_t)(v->ptr - data), end = start + v->len;
                if (start < lo) lo = start;
                if (end > hi) hi = end;
            }
        }
    }
    /* Mata reads the region as one string, which holds at most about 2 GB: a
     * longer region (only a single huge row can make one) sends its values as
     * bytes instead. Values travel verbatim on this route, so that is exact. */
    bool region = lo < hi && hi - lo <= 2147483620u;
    if (task->mem) {
        /* Positions in memory travel as Stata longs. A transformed value
         * has no decoded bytes to fall back on: the ado then uses a file. */
        if (lo < hi && hi - lo > 2147483620u) { task->rc = 459; goto done; }
        if (patch_len > 2147483620u) { task->rc = 459; goto done; }
        CImportBlobMem *m = task->mem;
        /* A column more than half full lists every row, in order, so that
         * Mata stores it as one range; otherwise only its nonempty rows. */
        unsigned char *dense = (unsigned char *)chunk_quotes; /* reuse: one byte per column */
        size_t nent = 0;
        for (int k = 0; k < ncols; k++) {
            size_t full = 0;
            for (size_t i = 0; i < b; i++) full += values[i * (size_t)ncols + (size_t)k].kind != CIMPORT_BLOB_EMPTY;
            dense[k] = 2 * full > b;
            nent += dense[k] ? b : full;
        }
        m->rows = (uint32_t)b;
        m->start = lo < hi ? task->base + lo : 0;
        m->length = lo < hi ? hi - lo : 0;
        m->counts = ctools_safe_calloc2((size_t)ncols, sizeof(*m->counts));
        m->r = ctools_safe_malloc2(nent ? nent : 1, sizeof(*m->r));
        m->p = ctools_safe_malloc2(nent ? nent : 1, sizeof(*m->p));
        m->n = ctools_safe_malloc2(nent ? nent : 1, sizeof(*m->n));
        m->t = malloc(nent ? nent : 1);
        if (!m->counts || !m->r || !m->p || !m->n || !m->t) { task->rc = 920; goto done; }
        size_t e = 0;
        for (int k = 0; k < ncols; k++)
            for (size_t i = 0; i < b; i++) {
                const CImportBlobValue *v = &values[i * (size_t)ncols + (size_t)k];
                if (v->kind == CIMPORT_BLOB_EMPTY && !dense[k]) continue;
                m->r[e] = (int32_t)i + 1;
                m->n[e] = (int32_t)v->len;
                if (v->kind == CIMPORT_BLOB_EMPTY) {
                    m->p[e] = 1;
                    m->n[e] = 0;
                    m->t[e] = 0;
                } else if (v->kind == CIMPORT_BLOB_DIRECT) {
                    m->p[e] = (int32_t)((uint64_t)(v->ptr - data) - lo + 1);
                    m->t[e] = v->t;
                } else if (v->ptr) { /* verbatim bytes outside a region: copied */
                    task->rc = 459;
                    goto done;
                } else {
                    m->p[e] = (int32_t)v->off + 1;
                    m->t[e] = CIMPORT_BLOB_TEXT;
                }
                m->counts[k]++;
                e++;
            }
        m->nent = nent;
        m->text = patch;
        m->ntext = patch_len;
        patch = NULL;
        goto done;
    }
    if (!region) {
        for (size_t i = 0; i < b * (size_t)ncols; i++)
            if (values[i].kind == CIMPORT_BLOB_DIRECT) values[i].kind = CIMPORT_BLOB_BYTES;
        lo = hi = 0;
    }

    /* Size, then serialize. */
    unsigned char *mode = (unsigned char *)chunk_quotes; /* reuse: one byte per column */
    size_t size = sizeof(uint32_t) + 2 * sizeof(double);
    for (int k = 0; k < ncols; k++) {
        size_t direct = 0, patches = 0, bytes = 0;
        for (size_t i = 0; i < b; i++) {
            const CImportBlobValue *v = &values[i * (size_t)ncols + (size_t)k];
            if (v->kind == CIMPORT_BLOB_DIRECT) direct++;
            else if (v->kind == CIMPORT_BLOB_BYTES) { patches++; bytes += v->len; }
        }
        size_t m = direct + patches, n = 2 * m <= b ? m : b;
        mode[k] = !m ? 2 : (direct ? 1 : 0) | (n < b ? 4 : 0);
        size += sizeof(uint32_t);
        if (!m) continue;
        if (n < b) size += sizeof(uint32_t) + n * sizeof(uint32_t);
        size += n * sizeof(uint32_t) + bytes;
        if (mode[k] & 1) size += n * sizeof(uint32_t) + sizeof(uint32_t) + 2 * patches * sizeof(uint32_t);
    }
    char *out = malloc(size ? size : 1), *p = out;
    if (!out) { task->rc = 920; goto done; }
    uint32_t u = (uint32_t)b;
    double start = (double)(task->base + lo), length = (double)(hi - lo);
    cimport_blob_put(&p, &u, sizeof(u));
    cimport_blob_put(&p, &start, sizeof(start));
    cimport_blob_put(&p, &length, sizeof(length));
    for (int k = 0; k < ncols; k++) {
        u = mode[k];
        cimport_blob_put(&p, &u, sizeof(u));
        if (mode[k] == 2) continue;
        const CImportBlobValue *col = values + k;
        const bool sparse = (mode[k] & 4) != 0;
#define CIMPORT_BLOB_AT(i) (&col[(i) * (size_t)ncols])
#define CIMPORT_BLOB_BYTES_OF(v) ((v)->ptr ? (v)->ptr : patch + (v)->off)
#define CIMPORT_BLOB_LISTED(i) (!sparse || CIMPORT_BLOB_AT(i)->kind != CIMPORT_BLOB_EMPTY)
        if (sparse) {
            uint32_t n = 0;
            for (size_t i = 0; i < b; i++) n += CIMPORT_BLOB_AT(i)->kind != CIMPORT_BLOB_EMPTY;
            cimport_blob_put(&p, &n, sizeof(n));
            for (size_t i = 0; i < b; i++)
                if (CIMPORT_BLOB_LISTED(i)) { u = (uint32_t)i + 1; cimport_blob_put(&p, &u, sizeof(u)); }
        }
        if (!(mode[k] & 1)) {
            for (size_t i = 0; i < b; i++)
                if (CIMPORT_BLOB_LISTED(i)) { u = CIMPORT_BLOB_AT(i)->len; cimport_blob_put(&p, &u, sizeof(u)); }
            for (size_t i = 0; i < b; i++) {
                const CImportBlobValue *v = CIMPORT_BLOB_AT(i);
                if (v->len) cimport_blob_put(&p, CIMPORT_BLOB_BYTES_OF(v), v->len);
            }
            continue;
        }
        uint32_t patches = 0;
        for (size_t i = 0; i < b; i++) {
            if (!CIMPORT_BLOB_LISTED(i)) continue;
            const CImportBlobValue *v = CIMPORT_BLOB_AT(i);
            u = v->kind == CIMPORT_BLOB_DIRECT ? (uint32_t)((uint64_t)(v->ptr - data) - lo + 1) : 1;
            cimport_blob_put(&p, &u, sizeof(u));
            patches += v->kind == CIMPORT_BLOB_BYTES;
        }
        for (size_t i = 0; i < b; i++) {
            if (!CIMPORT_BLOB_LISTED(i)) continue;
            const CImportBlobValue *v = CIMPORT_BLOB_AT(i);
            u = v->kind == CIMPORT_BLOB_DIRECT ? v->len : 0;
            cimport_blob_put(&p, &u, sizeof(u));
        }
        cimport_blob_put(&p, &patches, sizeof(patches));
        uint32_t j = 0; /* index among the listed rows */
        for (size_t i = 0; i < b; i++) {
            if (!CIMPORT_BLOB_LISTED(i)) continue;
            j++;
            if (CIMPORT_BLOB_AT(i)->kind == CIMPORT_BLOB_BYTES) cimport_blob_put(&p, &j, sizeof(j));
        }
        for (size_t i = 0; i < b; i++)
            if (CIMPORT_BLOB_AT(i)->kind == CIMPORT_BLOB_BYTES) { u = CIMPORT_BLOB_AT(i)->len; cimport_blob_put(&p, &u, sizeof(u)); }
        for (size_t i = 0; i < b; i++) {
            const CImportBlobValue *v = CIMPORT_BLOB_AT(i);
            if (v->kind == CIMPORT_BLOB_BYTES) cimport_blob_put(&p, CIMPORT_BLOB_BYTES_OF(v), v->len);
        }
#undef CIMPORT_BLOB_AT
#undef CIMPORT_BLOB_BYTES_OF
#undef CIMPORT_BLOB_LISTED
    }
    if ((size_t)(p - out) != size) { free(out); task->rc = 459; goto done; }
    task->out = out;
    task->out_len = size;
done:
    free(values);
    free(fields);
    free(chunk_quotes);
    free(scratch);
    free(patch);
    free(cursor.parsed);
    return NULL;
}


/* A short-lived file. On Windows it is marked temporary, which asks the cache
 * manager to keep it in memory rather than write it out; the ado erases it
 * once read. */
static FILE *cimport_open_temporary(const char *path) {
#ifdef _WIN32
    HANDLE h = ctools_create_file(path, GENERIC_WRITE, FILE_SHARE_READ, CREATE_ALWAYS,
                                  FILE_ATTRIBUTE_TEMPORARY);
    if (h == INVALID_HANDLE_VALUE) return NULL;
    int fd = _open_osfhandle((intptr_t)h, _O_WRONLY | _O_BINARY);
    if (fd < 0) { CloseHandle(h); return NULL; }
    FILE *f = _fdopen(fd, "wb");
    if (!f) _close(fd);
    return f;
#else
    return ctools_fopen(path, "wb");
#endif
}

/* The in-memory strL transfer, for Stata 16 and later (cimport.ado,
 * _cimport_read_blobmem). cimport blobmem builds every block of the strL
 * columns named in local __cimport_longcolumns as CImportBlobMem, held in the
 * context: values found in the source, verbatim or after Mata's transforms,
 * travel as positions, and the rest as text. cimport blobfill g then stores
 * block g's values in the four variables of its call (a frame of the ado's:
 * row, position, length, kind), with its row count, region and text size in
 * global _cimport_blobinfo and the values per column in _cimport_blobcounts;
 * cimport blobtext g offset max puts up to max bytes of the block's text from
 * offset in global _cimport_blobtext and the next offset in
 * _cimport_blobnext. Nothing is written to a file. blobmem reports in global
 * _cimport_blobmemrc and returns 0, so that a refusal (459: the source cannot
 * be cut in place) keeps the scan for the file transfer. */
static void cimport_blobmem_clear(CImportContext *ctx) {
    CImportBlobMem *mem = ctx->blobmem;
    for (size_t g = 0; g < ctx->nblobmem; g++) cimport_blobmem_free(&mem[g]);
    free(mem);
    ctx->blobmem = NULL;
    ctx->nblobmem = 0;
}

static ST_retcode cimport_blobmem(void) {
    CImportContext *ctx = g_cimport_ctx;
    if (!ctx) return 198;
    cimport_blobmem_clear(ctx);
    CImportStrlPlan plan;
    ST_retcode rc = cimport_strl_plan(ctx, &plan);
    if (rc) return rc;
    char flag[8] = "";
    SF_macro_use("___cimport_blobdirect", flag, sizeof(flag) - 1);
    bool direct = strcmp(flag, "0") && !ctx->converted_data && ctx->mmap_data && ctx->filename &&
                  ctx->file_data >= ctx->mmap_data && ctx->file_data <= ctx->mmap_data + ctx->mmap_size;
    uint64_t base = direct ? (uint64_t)(ctx->file_data - ctx->mmap_data) : 0;
    if (!direct) rc = 459;
    size_t ntasks = plan.nblocks;
    CImportBlobTask *tasks = rc ? NULL : ctools_safe_malloc2(ntasks ? ntasks : 1, sizeof(*tasks));
    CImportBlobMem *mem = rc ? NULL : ctools_safe_calloc2(ntasks ? ntasks : 1, sizeof(*mem));
    if (!rc && (!tasks || !mem)) rc = 920;
    for (size_t t = 0; t < ntasks && !rc; t++)
        tasks[t] = (CImportBlobTask){ctx, plan.cols, plan.ncols, plan.max_col, plan.sorted, true, base,
                                     plan.scratch_cap, plan.blocks[t].chunk, plan.blocks[t].r0,
                                     plan.blocks[t].r1, NULL, 0, 0, &plan, &plan.blocks[t], &mem[t]};
    int threads = ctx->num_threads > 0 ? ctx->num_threads : 1;
    ctools_persistent_pool *pool = threads > 1 && ntasks > 1 ? ctools_get_global_pool() : NULL;
    if (!rc && (pool == NULL || ctools_persistent_pool_submit_batch(pool, cimport_blob_worker, tasks, ntasks,
                                                                    sizeof(CImportBlobTask)) != 0)) {
        for (size_t t = 0; t < ntasks; t++) cimport_blob_worker(&tasks[t]);
    } else if (!rc) {
        ctools_persistent_pool_wait(pool);
    }
    size_t most = 0;
    for (size_t t = 0; t < ntasks && !rc; t++) {
        if (!rc) rc = tasks[t].rc;
        if (mem[t].nent > most) most = mem[t].nent;
    }
    if (rc && mem) {
        for (size_t t = 0; t < ntasks; t++) cimport_blobmem_free(&mem[t]);
        free(mem);
        mem = NULL;
    }
    free(tasks);
    ctx->blobmem = mem;
    ctx->nblobmem = rc ? 0 : ntasks;
    ctx->nblobcols = plan.ncols;
    cimport_strl_plan_free(&plan);
    char text[64];
    snprintf(text, sizeof(text), "%zu", ctx->nblobmem);
    ST_retcode mrc = SF_macro_save("_cimport_blobblocks", text);
    snprintf(text, sizeof(text), "%zu", most);
    if (!mrc) mrc = SF_macro_save("_cimport_blobmax", text);
    snprintf(text, sizeof(text), "%.0f", (double)ctx->mmap_size);
    if (!mrc) mrc = SF_macro_save("_cimport_blobsource", text);
    snprintf(text, sizeof(text), "%d", (int)rc);
    if (!mrc) mrc = SF_macro_save("_cimport_blobmemrc", text);
    return mrc;
}

typedef struct {
    const CImportBlobMem *m;
    size_t e0, e1;
    ST_retcode rc;
} CImportBlobFill;

static void *cimport_blobfill_worker(void *arg) {
    CImportBlobFill *f = arg;
    ST_IIID store = (_stata_)->safestore;
    for (size_t e = f->e0; e < f->e1 && !f->rc; e++) {
        ST_int obs = (ST_int)e + 1;
        f->rc = store(1, obs, (double)f->m->r[e]);
        if (!f->rc) f->rc = store(2, obs, (double)f->m->p[e]);
        if (!f->rc) f->rc = store(3, obs, (double)f->m->n[e]);
        if (!f->rc) f->rc = store(4, obs, (double)f->m->t[e]);
    }
    return NULL;
}

static ST_retcode cimport_blobfill(size_t g) {
    CImportContext *ctx = g_cimport_ctx;
    if (!ctx || g >= ctx->nblobmem) return 198;
    const CImportBlobMem *m = &((CImportBlobMem *)ctx->blobmem)[g];
    if (SF_nvars() != 4 || (size_t)SF_nobs() < m->nent) return 198;
    ST_retcode rc = 0;
    enum { PIECES = 32 };
    CImportBlobFill fill[PIECES];
    int threads = ctx->num_threads > 0 ? ctx->num_threads : 1;
    size_t n = m->nent < 65536 || threads < 2 ? 1 : (size_t)(threads < PIECES ? threads : PIECES);
    for (size_t i = 0; i < n; i++)
        fill[i] = (CImportBlobFill){m, m->nent * i / n, m->nent * (i + 1) / n, 0};
    ctools_persistent_pool *pool = n > 1 ? ctools_get_global_pool() : NULL;
    if (pool == NULL || ctools_persistent_pool_submit_batch(pool, cimport_blobfill_worker, fill, n,
                                                            sizeof(CImportBlobFill)) != 0) {
        for (size_t i = 0; i < n; i++) cimport_blobfill_worker(&fill[i]);
    } else {
        ctools_persistent_pool_wait(pool);
    }
    for (size_t i = 0; i < n && !rc; i++) rc = fill[i].rc;
    if (rc) return rc;
    char info[160];
    snprintf(info, sizeof(info), "%u %.0f %.0f %zu %zu", m->rows, (double)m->start, (double)m->length,
             m->ntext, m->nent);
    rc = SF_macro_save("_cimport_blobinfo", info);
    size_t cap = 12 * (size_t)ctx->nblobcols + 2;
    char *list = malloc(cap), *p = list;
    if (!list) return 920;
    *p = 0;
    for (int k = 0; k < ctx->nblobcols; k++) p += snprintf(p, cap - (size_t)(p - list), k ? " %u" : "%u", m->counts[k]);
    if (!rc) rc = SF_macro_save("_cimport_blobcounts", list);
    free(list);
    return rc;
}

static ST_retcode cimport_blobtext(size_t g, size_t offset, size_t max) {
    CImportContext *ctx = g_cimport_ctx;
    if (!ctx || g >= ctx->nblobmem || max < 1) return 198;
    const CImportBlobMem *m = &((CImportBlobMem *)ctx->blobmem)[g];
    if (offset > m->ntext) return 198;
    size_t n = m->ntext - offset < max ? m->ntext - offset : max;
    char *piece = malloc(n + 1);
    if (!piece) return 920;
    memcpy(piece, m->text + offset, n);
    piece[n] = 0; /* the text holds no zero bytes: they were removed */
    ST_retcode rc = SF_macro_save("_cimport_blobtext", piece);
    free(piece);
    char next[32];
    snprintf(next, sizeof(next), "%zu", offset + n);
    if (!rc) rc = SF_macro_save("_cimport_blobnext", next);
    return rc;
}

/* Every strL column named in local __cimport_longcolumns, in that order, to
 * the file named by local __cimport_blobpath, in one plugin call. Local
 * __cimport_blobdirect 0 sends every value as bytes (Mata could not reopen
 * the source). */
static ST_retcode cimport_blobs(void) {
    CImportContext *ctx = g_cimport_ctx;
    if (!ctx) return 198;
    CImportStrlPlan plan;
    ST_retcode rc = cimport_strl_plan(ctx, &plan);
    if (rc) return rc;
    int *cols = plan.cols, ncols = plan.ncols;
    char path[4096] = "", flag[8] = "";
    if (!rc && (SF_macro_use("___cimport_blobpath", path, sizeof(path) - 1) || !*path)) rc = 198;
    SF_macro_use("___cimport_blobdirect", flag, sizeof(flag) - 1);
    bool direct = !rc && strcmp(flag, "0") && !ctx->converted_data && ctx->mmap_data &&
                  ctx->filename && ctx->file_data >= ctx->mmap_data &&
                  ctx->file_data <= ctx->mmap_data + ctx->mmap_size;
    uint64_t base = direct ? (uint64_t)(ctx->file_data - ctx->mmap_data) : 0;

    size_t ntasks = plan.nblocks;
    CImportBlobTask *tasks = rc ? NULL : ctools_safe_malloc2(ntasks ? ntasks : 1, sizeof(*tasks));
    if (!rc && !tasks) rc = 920;
    for (size_t t = 0; t < ntasks && !rc; t++)
        tasks[t] = (CImportBlobTask){ctx, cols, ncols, plan.max_col, plan.sorted, direct, base,
                                     plan.scratch_cap, plan.blocks[t].chunk, plan.blocks[t].r0,
                                     plan.blocks[t].r1, NULL, 0, 0, &plan, &plan.blocks[t], NULL};

    FILE *out = NULL;
    if (!rc && !(out = cimport_open_temporary(path))) rc = 603;
    if (!rc) {
        uint32_t head[3] = {CIMPORT_BLOB_VERSION, (uint32_t)ncols, direct ? 1u : 0u};
        double source = direct ? (double)ctx->mmap_size : 0;
        uint32_t plen = direct ? (uint32_t)strlen(ctx->filename) : 0;
        if (fwrite(head, sizeof(head), 1, out) != 1 || fwrite(&source, sizeof(source), 1, out) != 1 ||
            fwrite(&plen, sizeof(plen), 1, out) != 1 ||
            (plen && fwrite(ctx->filename, 1, plen, out) != plen)) rc = 693;
    }

    /* Blocks are built in parallel a wave at a time and written in order, so
     * at most one wave of serialized blocks is held at once. */
    int threads = ctx->num_threads > 0 ? ctx->num_threads : 1;
    size_t wave = (size_t)threads * 2;
    ctools_persistent_pool *pool = threads > 1 ? ctools_get_global_pool() : NULL;
    for (size_t t0 = 0; t0 < ntasks && !rc; t0 += wave) {
        size_t n = ntasks - t0 < wave ? ntasks - t0 : wave;
        if (n < 2 || pool == NULL ||
            ctools_persistent_pool_submit_batch(pool, cimport_blob_worker, tasks + t0, n,
                                                sizeof(CImportBlobTask)) != 0) {
            for (size_t t = 0; t < n; t++) cimport_blob_worker(&tasks[t0 + t]);
        } else {
            ctools_persistent_pool_wait(pool);
        }
        for (size_t t = 0; t < n; t++) {
            CImportBlobTask *task = &tasks[t0 + t];
            if (!rc) rc = task->rc;
            if (!rc && fwrite(task->out, 1, task->out_len, out) != task->out_len) rc = 693;
            free(task->out);
            task->out = NULL;
        }
    }
    if (out && fclose(out) && !rc) rc = 693;
    free(tasks);
    cimport_strl_plan_free(&plan);
    return rc;
}

/* strL transfer as a Stata dataset. Mata can store a strL only through
 * st_sstore(), several times slower than use reads one, so when the import
 * may replace the dataset outright cimport.ado has the strL columns written
 * as a .dta of N observations, uses it, then adds the other variables and
 * orders them (see _cimport_addvars). The shared writer (ctools_dta_write)
 * pulls the values a block of rows at a time. A block is decoded as in
 * cimport_blob_worker, once for the data section's references and, unless
 * the writer could keep it, once more for the values, which follow every
 * reference in the file. */
typedef struct {
    CImportContext *ctx;
    const CImportStrlPlan *plan;
} CImportDtaSource;

typedef struct {
    char *text;              /* the block's decoded values */
    ctools_dta_span *spans;  /* columns x rows, column by column */
} CImportDtaBlock;

static void cimport_dta_release(void *arg, ctools_dta_group *group) {
    (void)arg;
    CImportDtaBlock *block = group->state;
    if (block) {
        free(block->text);
        free(block->spans);
        free(block);
    }
    group->state = NULL;
}

static int cimport_dta_fill(void *arg, size_t g, int phase, ctools_dta_group *out) {
    (void)phase; /* both passes read only the strL columns */
    CImportDtaSource *source = arg;
    CImportContext *ctx = source->ctx;
    const CImportStrlPlan *plan = source->plan;
    const CImportStrlBlock *block = &plan->blocks[g];
    CImportParsedChunk *chunk = &ctx->chunks[block->chunk];
    const int ncols = plan->ncols;
    const size_t b = block->r1 - block->r0;
    const char quote = cimport_extract_quote(ctx);
    const char *data = ctx->file_data;
    /* Decoding never lengthens a value, so the bytes the rows span bound the
     * values; each also takes a terminating zero. */
    uint64_t bound = block->span + (uint64_t)b * (uint64_t)ncols + 1;
    CImportDtaBlock *state = calloc(1, sizeof(*state));
    CImportFieldRef *fields = ctools_safe_malloc2((size_t)ncols, sizeof(*fields));
    bool *chunk_quotes = malloc((size_t)ncols);
    CImportStrlCursor cursor;
    bool cursor_ok = cimport_strl_cursor_init(&cursor, ctx, plan, block);
    int rc = 0;
    if (state) {
        state->text = bound <= SIZE_MAX ? malloc((size_t)bound) : NULL;
        state->spans = ctools_safe_calloc2(b * (size_t)ncols, sizeof(*state->spans));
    }
    if (!state || !state->text || !state->spans || !fields || !chunk_quotes || !cursor_ok) { rc = 920; goto done; }
    for (int k = 0; k < ncols; k++) {
        int col = plan->cols[k];
        chunk_quotes[k] = ctx->columns[col].column_has_quotes && chunk->col_stats &&
                          col < chunk->num_col_stats && chunk->col_stats[col].has_quotes;
    }
    char *text = state->text;
    for (size_t i = 0; i < b; i++) {
        int nfields;
        if (!cimport_strl_next(&cursor, fields, &nfields)) { rc = 459; goto done; }
        for (int k = 0; k < ncols; k++) {
            if (plan->cols[k] >= nfields) continue;
            CImportFieldRef *field = &fields[k];
            const char *src = data + field->offset;
            int raw = (int)CIMPORT_FIELD_LENGTH(*field), n;
            if (raw >= INT_MAX) { rc = 109; goto done; }
            if (!chunk_quotes[k] || !(field->length & CIMPORT_FIELD_QUOTED_FLAG)) {
                n = raw;
                while (n > 0 && (src[n - 1] == '\r' || src[n - 1] == '\n')) n--;
                if (n) memcpy(text, src, (size_t)n);
            } else {
                n = cimport_extract_field_fast(data, field, text, raw + 1, quote);
            }
            if (n <= 0) continue;
            state->spans[(size_t)k * b + i] = (ctools_dta_span){text, (size_t)n, 0};
            text += (size_t)n + 1;
        }
    }
    for (int k = 0; k < ncols; k++) out->columns[k] = state->spans + (size_t)k * b;
    out->state = state;
    state = NULL;
done:
    if (state) {
        free(state->text);
        free(state->spans);
        free(state);
    }
    free(fields);
    free(chunk_quotes);
    free(cursor.parsed);
    return rc;
}

/* The strL columns named in local __cimport_longcolumns as variables named
 * by local __cimport_strlnames (one name per column, in order; its byte
 * length in __cimport_strlnames_len), written to the .dta named by local
 * __cimport_blobpath. */
static ST_retcode cimport_strldta(void) {
    CImportContext *ctx = g_cimport_ctx;
    if (!ctx) return 198;
    CImportStrlPlan plan;
    ST_retcode rc = cimport_strl_plan(ctx, &plan);
    if (rc) return rc;
    const int ncols = plan.ncols;
    if (ncols > 32767) rc = 198; /* format 118 holds up to 32,767 variables */

    /* Names */
    char length_text[32] = "", *names = NULL, **name = NULL;
    int length = 0;
    if (!rc && (SF_macro_use("___cimport_strlnames_len", length_text, sizeof(length_text) - 1) ||
                !ctools_safe_atoi(length_text, &length) || length < 1 || length > INT_MAX - 2))
        rc = 198;
    if (!rc && (!(names = calloc((size_t)length + 2, 1)) ||
                !(name = ctools_safe_malloc2((size_t)ncols, sizeof(*name))))) rc = 920;
    if (!rc && (SF_macro_use("___cimport_strlnames", names, length + 1) ||
                strlen(names) != (size_t)length)) rc = 198;
    if (!rc) {
        int k = 0;
        for (char *p = names; *p && !rc;) {
            while (*p == ' ') p++;
            if (!*p) break;
            char *e = p;
            while (*e && *e != ' ') e++;
            if (k == ncols || e - p > 128) { rc = 198; break; }
            name[k++] = p;
            if (*e) *e++ = 0;
            p = e;
        }
        if (!rc && k != ncols) rc = 198;
    }
    char path[4096] = "";
    if (!rc && (SF_macro_use("___cimport_blobpath", path, sizeof(path) - 1) || !*path)) rc = 198;


    uint64_t *starts = NULL, *bytes = NULL;
    ctools_dta_column *cols = NULL;
    if (!rc && (!(starts = ctools_safe_malloc2(plan.nblocks + 1, sizeof(*starts))) ||
                !(bytes = ctools_safe_malloc2(plan.nblocks ? plan.nblocks : 1, sizeof(*bytes))) ||
                !(cols = ctools_safe_calloc2((size_t)ncols, sizeof(*cols))))) rc = 920;
    if (!rc) {
        for (size_t t = 0; t < plan.nblocks; t++) {
            const CImportStrlBlock *block = &plan.blocks[t];
            uint64_t cells = (uint64_t)(block->r1 - block->r0) * (uint64_t)ncols;
            starts[t] = block->obs0 - 1;
            bytes[t] = block->span + cells * (1 + sizeof(ctools_dta_span)) + 1;
        }
        starts[plan.nblocks] = ctx->total_rows;
        for (int k = 0; k < ncols; k++)
            cols[k] = (ctools_dta_column){.name = name[k], .type = CTOOLS_DTA_STRL, .format = "%9s",
                                          .source = CTOOLS_DTA_SRC_PRODUCER, .flags = CTOOLS_DTA_COL_TEXT};
        CImportDtaSource source = {ctx, &plan};
        ctools_dta_producer producer = {plan.nblocks, starts, bytes, &source, cimport_dta_fill, cimport_dta_release};
        ctools_dta_spec spec;
        ctools_dta_spec_init(&spec);
        spec.nobs = ctx->total_rows;
        spec.ncols = (size_t)ncols;
        spec.cols = cols;
        spec.producer = &producer;
        spec.threads = ctx->num_threads > 0 ? ctx->num_threads : 1;
        spec.temporary = 1;
        rc = ctools_dta_write(path, &spec, NULL);
    }
    free(starts);
    free(bytes);
    free(cols);
    free(name);
    free(names);
    cimport_strl_plan_free(&plan);
    return rc;
}

/* The whole import as one Stata dataset (option tempdta of cimport.ado):
 * every selected column, named and typed by locals __cimport_fullnames and
 * __cimport_fulltypes (byte lengths in their _len locals), written to the
 * .dta named by local __cimport_blobpath. The ado uses it in place of
 * creating variables, storing through the SPI and transferring strLs. The
 * store workers fill the writer's groups through a CImportSink, so every
 * value is the one the SPI load would store: retained rows a tile at a time,
 * compact chunks a cached block at a time plus the uncached suffix. The
 * outcome goes to global _cimport_fullrc, so that a refusal keeps the scan
 * for the ado's next route: 0, 459 when the rows cannot be grouped that way
 * (only the header was kept, or a long uncached suffix), or the writer's
 * code. The call itself fails only without a scan to write. */
#define CIMPORT_FULL_SUFFIX_ROWS 65536
#define CIMPORT_FULL_SUFFIX_BYTES ((uint64_t)256 << 20)

typedef struct {
    int chunk;
    size_t row_start, row_end;        /* as in CImportStoreTile */
    const CImportCachedBlock *block;  /* compact: a cached block, or NULL for the suffix */
    uint64_t obs0;                    /* 0-based observation of the first row */
    size_t rows;
    uint64_t span;                    /* source bytes the rows span */
} CImportFullGroup;

typedef struct {
    CImportContext *ctx;
    const CImportFullGroup *groups;
    const unsigned char *kind[2];     /* per file column, for each writer phase */
} CImportFullSource;

typedef struct {
    void *arrays;
    char *text;
} CImportFullState;

static void cimport_full_release(void *arg, ctools_dta_group *group) {
    (void)arg;
    CImportFullState *state = group->state;
    if (state) {
        free(state->arrays);
        free(state->text);
        free(state);
    }
    group->state = NULL;
}

static int cimport_full_fill(void *arg, size_t g, int phase, ctools_dta_group *out) {
    CImportFullSource *source = arg;
    CImportContext *ctx = source->ctx;
    const CImportFullGroup *grp = &source->groups[g];
    const unsigned char *kind = source->kind[phase == CTOOLS_DTA_PHASE_STRLS ? 1 : 0];
    const int ncols = ctx->num_columns, first = ctx->col_first;
    size_t bytes = 0;
    for (int col = 0; col < ncols; col++) {
        if (kind[col] == CIMPORT_STORE_NUMERIC) bytes += grp->rows * sizeof(double);
        else if (kind[col] != CIMPORT_STORE_SKIP) bytes += grp->rows * sizeof(ctools_dta_span);
    }
    double **num = ctools_safe_calloc2((size_t)ncols, sizeof(*num));
    ctools_dta_span **str = ctools_safe_calloc2((size_t)ncols, sizeof(*str));
    CImportFullState *state = calloc(1, sizeof(*state));
    int rc = 0;
    if (state) {
        state->arrays = calloc(bytes ? bytes : 1, 1); /* empty spans: "" */
        state->text = grp->span < SIZE_MAX ? malloc((size_t)grp->span + 1) : NULL;
    }
    if (!num || !str || !state || !state->arrays || !state->text) { rc = 920; goto done; }
    char *p = state->arrays;
    for (int col = 0; col < ncols; col++) {
        if (kind[col] == CIMPORT_STORE_NUMERIC) { num[col] = (double *)p; p += grp->rows * sizeof(double); }
        else if (kind[col] != CIMPORT_STORE_SKIP) { str[col] = (ctools_dta_span *)p; p += grp->rows * sizeof(ctools_dta_span); }
    }
    CImportSink sink = {num, str, state->text, (ST_int)(grp->obs0 + 1)};
    CImportStoreTile tile = {ctx, kind, grp->chunk, grp->row_start, grp->row_end, (ST_int)(grp->obs0 + 1),
                             0, &sink, grp->block};
    if (!ctx->compact_rows) cimport_store_tile_worker(&tile);
    else if (grp->block) cimport_store_compact_worker(&tile);
    else cimport_store_stream_worker(&tile);
    rc = tile.error;
    if (!rc && (uint64_t)(sink.text - state->text) > grp->span) rc = 459; /* cannot happen */
    if (rc) goto done;
    for (int col = first; col < cimport_col_stop(ctx); col++)
        out->columns[col - first] = num[col] ? (void *)num[col] : (void *)str[col];
    out->state = state;
    state = NULL;
done:
    if (state) {
        free(state->arrays);
        free(state->text);
        free(state);
    }
    free(num);
    free(str);
    return rc;
}

/* The next space-separated token of a list read from a local; NULL at its end. */
static char *cimport_next_token(char **p) {
    while (**p == ' ') (*p)++;
    if (!**p) return NULL;
    char *t = *p;
    while (**p && **p != ' ') (*p)++;
    if (**p) *(*p)++ = 0;
    return t;
}

static char *cimport_read_local(const char *name) {
    char key[64], length_text[32] = "";
    int length = 0;
    snprintf(key, sizeof(key), "%s_len", name);
    if (SF_macro_use(key, length_text, sizeof(length_text) - 1) ||
        !ctools_safe_atoi(length_text, &length) || length < 1 || length > INT_MAX - 2)
        return NULL;
    char *text = calloc((size_t)length + 2, 1);
    if (text && (SF_macro_use((char *)name, text, length + 1) || strlen(text) != (size_t)length)) {
        free(text);
        text = NULL;
    }
    return text;
}

static ST_retcode cimport_fulldta(void) {
    CImportContext *ctx = g_cimport_ctx;
    if (!ctx || !ctx->total_rows) return 198;
    /* Store workers report allocation failures through error_code. */
    atomic_store(&ctx->error_code, 0);
    double t_start = ctools_timer_ms();
    const int first = ctx->col_first, stop = cimport_col_stop(ctx);
    const int K = stop - first;
    ST_retcode rc = 0;
    if (ctx->stream_rows || !ctx->file_data || !ctx->chunks) rc = 459;
    else if (K < 1 || K > 32767) rc = 198;

    char *names = cimport_read_local("___cimport_fullnames");
    char *types = cimport_read_local("___cimport_fulltypes");
    ctools_dta_column *cols = rc ? NULL : ctools_safe_calloc2((size_t)K, sizeof(*cols));
    unsigned char *kind = rc ? NULL : ctools_safe_malloc2((size_t)ctx->num_columns, 2);
    if (!rc && (!names || !types)) rc = 198;
    else if (!rc && (!cols || !kind)) rc = 920;

    /* Columns: the ado's names and types, which must agree with the scan */
    char *np = names, *tp = types;
    for (int col = 0; col < ctx->num_columns && !rc; col++) {
        kind[col] = kind[ctx->num_columns + col] = CIMPORT_STORE_SKIP;
        if (col < first || col >= stop) continue;
        const CImportColumnInfo *info = &ctx->columns[col];
        char *name = cimport_next_token(&np), *type = cimport_next_token(&tp);
        if (!name || !type) { rc = 198; break; }
        ctools_dta_column *c = &cols[col - first];
        c->name = name;
        c->source = CTOOLS_DTA_SRC_PRODUCER;
        bool string = info->type == CIMPORT_COL_STRING;
        if (!strcmp(type, "strL")) {
            if (!string || !info->use_strl) rc = 198;
            c->type = CTOOLS_DTA_STRL;
            c->flags = CTOOLS_DTA_COL_TEXT; /* binary zeros were removed */
            kind[col] = kind[ctx->num_columns + col] = CIMPORT_STORE_STRL;
        } else if (!strncmp(type, "str", 3)) {
            int width = atoi(type + 3);
            if (!string || info->use_strl || width < 1 || width > 2045) rc = 198;
            c->type = CTOOLS_DTA_STR;
            c->width = width;
            kind[col] = CIMPORT_STORE_STRING;
        } else {
            c->type = !strcmp(type, "byte") ? CTOOLS_DTA_BYTE : !strcmp(type, "int") ? CTOOLS_DTA_INT
                    : !strcmp(type, "long") ? CTOOLS_DTA_LONG : !strcmp(type, "float") ? CTOOLS_DTA_FLOAT
                    : !strcmp(type, "double") ? CTOOLS_DTA_DOUBLE : (ctools_dta_type)0;
            if (string || !c->type) rc = 198;
            kind[col] = CIMPORT_STORE_NUMERIC;
        }
    }
    if (!rc && (cimport_next_token(&np) || cimport_next_token(&tp))) rc = 198;

    /* Groups of rows, each filled by one store worker call */
    CImportFullGroup *groups = NULL;
    size_t ngroups = 0, cap = 0;
    uint64_t obs = 0;
#define CIMPORT_FULL_ADD(...) do { \
        if (ngroups == cap) { \
            cap = cap ? 2 * cap : 256; \
            CImportFullGroup *grown = realloc(groups, cap * sizeof(*grown)); \
            if (!grown) { rc = 920; break; } \
            groups = grown; \
        } \
        groups[ngroups++] = (CImportFullGroup){__VA_ARGS__}; \
    } while (0)
    for (int c = 0; c < ctx->num_chunks && !rc; c++) {
        CImportParsedChunk *chunk = &ctx->chunks[c];
        size_t header = (c == 0 && ctx->has_header) ? 1 : 0;
        if (ctx->compact_rows) {
            for (const CImportCachedBlock *block = chunk->cache; block && !rc; block = block->next) {
                if (!block->rows) continue;
                CIMPORT_FULL_ADD(c, 0, block->rows, block, obs, block->rows, block->file_end - block->file_start);
                obs += block->rows;
            }
            size_t rest = chunk->num_rows - header - chunk->cached_rows;
            if (!rc && rest) {
                uint64_t span = chunk->file_end - chunk->uncached_start;
                if (rest > CIMPORT_FULL_SUFFIX_ROWS || span > CIMPORT_FULL_SUFFIX_BYTES) { rc = 459; break; }
                size_t r0 = chunk->cached_rows ? 0 : header;
                CIMPORT_FULL_ADD(c, r0, r0 + rest, NULL, obs, rest, span);
                obs += rest;
            }
            continue;
        }
        for (size_t r = header; r < chunk->num_rows && !rc; r += CIMPORT_STORE_TILE_ROWS) {
            size_t end = r + CIMPORT_STORE_TILE_ROWS < chunk->num_rows ? r + CIMPORT_STORE_TILE_ROWS : chunk->num_rows;
            uint64_t span = cimport_row_end(chunk->rows[end - 1]) - chunk->rows[r]->offset + 1;
            CIMPORT_FULL_ADD(c, r, end, NULL, obs, end - r, span);
            obs += end - r;
        }
    }
#undef CIMPORT_FULL_ADD
    if (!rc && obs != ctx->total_rows) rc = 459;

    uint64_t *starts = NULL, *bytes = NULL;
    if (!rc && (!(starts = ctools_safe_malloc2(ngroups + 1, sizeof(*starts))) ||
                !(bytes = ctools_safe_malloc2(ngroups ? ngroups : 1, sizeof(*bytes))))) rc = 920;
    char path[4096] = "";
    if (!rc && (SF_macro_use("___cimport_blobpath", path, sizeof(path) - 1) || !*path)) rc = 198;
    if (!rc) {
        for (size_t g = 0; g < ngroups; g++) {
            starts[g] = groups[g].obs0;
            bytes[g] = groups[g].span + 1 + (uint64_t)groups[g].rows * (uint64_t)K * sizeof(ctools_dta_span);
        }
        starts[ngroups] = ctx->total_rows;
        CImportFullSource source = {ctx, groups, {kind, kind + ctx->num_columns}};
        ctools_dta_producer producer = {ngroups, starts, bytes, &source, cimport_full_fill, cimport_full_release};
        ctools_dta_spec spec;
        ctools_dta_spec_init(&spec);
        spec.nobs = ctx->total_rows;
        spec.ncols = (size_t)K;
        spec.cols = cols;
        spec.producer = &producer;
        spec.threads = ctx->num_threads > 0 ? ctx->num_threads : 1;
        spec.temporary = 1;
        ctools_dta_result result;
        memset(&result, 0, sizeof(result));
        rc = ctools_dta_write(path, &spec, &result);
        if (rc && ctx->verbose && *result.message) {
            char msg[320];
            snprintf(msg, sizeof(msg), "(cimport: %s; importing another way)\n", result.message);
            cimport_display_msg(msg);
        }
    }
    if (!rc) {
        double time_write = ctools_timer_ms() - t_start;
        SF_scal_save("_cimport_time_mmap", ctx->time_mmap / 1000.0);
        SF_scal_save("_cimport_time_encoding", ctx->time_encoding / 1000.0);
        SF_scal_save("_cimport_time_parse", ctx->time_parse / 1000.0);
        SF_scal_save("_cimport_time_infer", ctx->time_type_infer / 1000.0);
        SF_scal_save("_cimport_time_cache", ctx->time_cache / 1000.0);
        SF_scal_save("_cimport_time_store", time_write / 1000.0);
        SF_scal_save("_cimport_time_total", (ctx->time_mmap + ctx->time_encoding + ctx->time_parse +
                                             ctx->time_type_infer + ctx->time_cache + time_write) / 1000.0);
        CTOOLS_SAVE_THREAD_INFO("_cimport");
    }
    free(starts);
    free(bytes);
    free(groups);
    free(cols);
    free(kind);
    free(names);
    free(types);
    char status[16];
    snprintf(status, sizeof(status), "%d", (int)rc);
    return SF_macro_save("_cimport_fullrc", status);
}

ST_retcode cimport_main(const char *args) {
    if (args && !strcmp(args,"clear")) { cimport_clear_cached_context();return 0; }
    if (args && !strcmp(args,"blobs")) return cimport_blobs();
    if (args && !strcmp(args,"strldta")) return cimport_strldta();
    if (args && !strcmp(args,"blobmem")) return cimport_blobmem();
    if (args && !strncmp(args,"blobfill ",9)) {
        size_t g;
        if (sscanf(args, "%*s %zu", &g) != 1) return 198;
        return cimport_blobfill(g);
    }
    if (args && !strncmp(args,"blobtext ",9)) {
        size_t g, offset, max;
        if (sscanf(args, "%*s %zu %zu %zu", &g, &offset, &max) != 3) return 198;
        return cimport_blobtext(g, offset, max);
    }
    if (args && !strcmp(args,"fulldta")) return cimport_fulldta();
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
    if (strlen(args) >= sizeof(args_copy)) return 198;
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
    cimport_clear_cached_context();
}

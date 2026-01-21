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
#else
    #include <unistd.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_timer.h"
#include "cimport_impl.h"

/* Helper modules (cimport_context.h includes arena, parse, and mmap headers) */
#include "cimport_context.h"
#include "cimport_mmap.h"

/* High-resolution timing */
#define cimport_get_time_ms() ctools_timer_ms()

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/* Local: cimport may use more threads for parsing than general I/O */
#define CIMPORT_MAX_THREADS      32

/* Global cached context */
static CImportContext *g_cimport_ctx = NULL;

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

static void cimport_display_msg(const char *msg) {
    SF_display((char *)msg);
}

static void cimport_display_error(const char *msg) {
    SF_error((char *)msg);
}

/* Record a warning for unmatched quote (thread-safe) */
static void cimport_record_unmatched_quote(CImportContext *ctx, size_t row_num) {
    pthread_mutex_lock(&ctx->warning_mutex);
    if (ctx->num_unmatched_quote_warnings < CIMPORT_MAX_WARNINGS) {
        ctx->unmatched_quote_rows[ctx->num_unmatched_quote_warnings] = row_num;
    }
    ctx->num_unmatched_quote_warnings++;
    pthread_mutex_unlock(&ctx->warning_mutex);
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

    for (int i = 0; i < len && j < CTOOLS_MAX_VARNAME_LEN - 1; i++) {
        char c = input[i];
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_') {
            output[j++] = c;
        } else if (c >= '0' && c <= '9' && j > 0) {
            output[j++] = c;
        }
    }

    if (j == 0) {
        strcpy(output, "v");
        j = 1;
    }

    output[j] = '\0';
}

/* ============================================================================
 * Type Inference
 * ============================================================================ */

static CImportNumericSubtype cimport_determine_numeric_subtype(CImportColumnInfo *col) {
    if (col->min_value < -2147483647.0 || col->max_value > 2147483620.0) {
        return CIMPORT_NUM_DOUBLE;
    }

    if (!col->is_integer) {
        return CIMPORT_NUM_FLOAT;
    }

    if (col->min_value >= -127 && col->max_value <= 100) {
        return CIMPORT_NUM_BYTE;
    }
    if (col->min_value >= -32767 && col->max_value <= 32740) {
        return CIMPORT_NUM_INT;
    }

    return CIMPORT_NUM_LONG;
}

/* ============================================================================
 * Header and Type Inference
 * ============================================================================ */

static int cimport_parse_header(CImportContext *ctx) {
    CImportFieldRef fields[CTOOLS_MAX_COLUMNS];
    const char *end = cimport_find_next_row(ctx->file_data, ctx->file_data + ctx->file_size, ctx->quote_char, ctx->bindquotes);

    ctx->num_columns = cimport_parse_row_fast(ctx->file_data, end, ctx->delimiter, ctx->quote_char,
                                               fields, CTOOLS_MAX_COLUMNS, ctx->file_data);

    if (ctx->num_columns <= 0) {
        strcpy(ctx->error_message, "No columns found in header");
        return -1;
    }

    ctx->columns = ctools_safe_calloc2(ctx->num_columns, sizeof(CImportColumnInfo));
    if (!ctx->columns) return -1;

    char name_buf[CTOOLS_MAX_STRING_LEN];
    for (int i = 0; i < ctx->num_columns; i++) {
        if (ctx->has_header) {
            cimport_extract_field_fast(ctx->file_data, &fields[i], name_buf, sizeof(name_buf), ctx->quote_char);
            cimport_sanitize_varname(name_buf, strlen(name_buf), ctx->columns[i].name);
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
        if (list[i] == col_1based) return true;
    }
    return false;
}

static void cimport_infer_column_types(CImportContext *ctx) {
    for (int i = 0; i < ctx->num_columns; i++) {
        ctx->columns[i].is_integer = true;
        ctx->columns[i].min_value = DBL_MAX;
        ctx->columns[i].max_value = -DBL_MAX;
        ctx->columns[i].max_decimal_digits = 0;
        ctx->columns[i].num_subtype = CIMPORT_NUM_DOUBLE;
    }

    for (int col = 0; col < ctx->num_columns; col++) {
        CImportColumnInfo *colinfo = &ctx->columns[col];
        int col_1based = col + 1;

        /* Check for forced type overrides first */
        if (ctx->force_string_cols && cimport_col_in_list(col_1based, ctx->force_string_cols, ctx->num_force_string)) {
            /* Force this column to be string */
            colinfo->type = CIMPORT_COL_STRING;
            /* Still need to find max length */
            int max_len = 0;
            for (int c = 0; c < ctx->num_chunks; c++) {
                CImportParsedChunk *chunk = &ctx->chunks[c];
                if (chunk->col_stats && col < chunk->num_col_stats) {
                    if (chunk->col_stats[col].max_field_len > max_len) {
                        max_len = chunk->col_stats[col].max_field_len;
                    }
                }
            }
            colinfo->max_strlen = max_len > 0 ? max_len : 1;
            if (colinfo->max_strlen > CTOOLS_MAX_STRING_LEN) {
                colinfo->max_strlen = CTOOLS_MAX_STRING_LEN;
            }
            continue;
        }

        if (ctx->force_numeric_cols && cimport_col_in_list(col_1based, ctx->force_numeric_cols, ctx->num_force_numeric)) {
            /* Force this column to be numeric */
            colinfo->type = CIMPORT_COL_NUMERIC;
            continue;
        }

        /* Normal type inference */
        bool is_string = false;
        int max_len = 0;

        for (int c = 0; c < ctx->num_chunks; c++) {
            CImportParsedChunk *chunk = &ctx->chunks[c];
            if (chunk->col_stats && col < chunk->num_col_stats) {
                CImportColumnParseStats *stats = &chunk->col_stats[col];
                if (stats->seen_string) {
                    is_string = true;
                }
                if (stats->max_field_len > max_len) {
                    max_len = stats->max_field_len;
                }
            }
        }

        if (is_string) {
            colinfo->type = CIMPORT_COL_STRING;
            colinfo->max_strlen = max_len > 0 ? max_len : 1;
            if (colinfo->max_strlen > CTOOLS_MAX_STRING_LEN) {
                colinfo->max_strlen = CTOOLS_MAX_STRING_LEN;
            }
        } else {
            colinfo->type = CIMPORT_COL_NUMERIC;
        }
    }

    size_t total_data_rows = 0;
    for (int c = 0; c < ctx->num_chunks; c++) {
        size_t start_row = (c == 0 && ctx->has_header) ? 1 : 0;
        total_data_rows += ctx->chunks[c].num_rows - start_row;
    }

    for (int col_idx = 0; col_idx < ctx->num_columns; col_idx++) {
        CImportColumnInfo *col = &ctx->columns[col_idx];
        if (col->type != CIMPORT_COL_NUMERIC) continue;

        size_t sample_count = 0;
        size_t max_samples = 1000;
        size_t sample_step = (total_data_rows > max_samples) ? (total_data_rows / max_samples) : 1;
        size_t row_idx = 0;

        for (int c = 0; c < ctx->num_chunks && sample_count < max_samples; c++) {
            CImportParsedChunk *chunk = &ctx->chunks[c];
            size_t start_row = (c == 0 && ctx->has_header) ? 1 : 0;

            for (size_t r = start_row; r < chunk->num_rows && sample_count < max_samples; r++) {
                if (row_idx % sample_step == 0) {
                    CImportParsedRow *row = chunk->rows[r];
                    if (col_idx < row->num_fields) {
                        CImportFieldRef *field = &row->fields[col_idx];
                        double val;
                        bool is_int;

                        if (cimport_analyze_numeric_fast(ctx->file_data, field, ctx->quote_char, &val, &is_int)) {
                            if (!is_int) {
                                col->is_integer = false;
                            }
                            if (val < col->min_value) {
                                col->min_value = val;
                            }
                            if (val > col->max_value) {
                                col->max_value = val;
                            }
                            sample_count++;
                        }
                    }
                }
                row_idx++;
            }
        }
    }

    for (int i = 0; i < ctx->num_columns; i++) {
        if (ctx->columns[i].type == CIMPORT_COL_UNKNOWN) {
            ctx->columns[i].type = CIMPORT_COL_NUMERIC;
        }
        if (ctx->columns[i].type == CIMPORT_COL_STRING && ctx->columns[i].max_strlen == 0) {
            ctx->columns[i].max_strlen = 1;
        }
        if (ctx->columns[i].max_strlen > CTOOLS_MAX_STRING_LEN) {
            ctx->columns[i].max_strlen = CTOOLS_MAX_STRING_LEN;
        }

        if (ctx->columns[i].type == CIMPORT_COL_NUMERIC) {
            /* Check for forced numeric type mode */
            if (ctx->numeric_type_mode == CIMPORT_NUMTYPE_FLOAT) {
                ctx->columns[i].num_subtype = CIMPORT_NUM_FLOAT;
            } else if (ctx->numeric_type_mode == CIMPORT_NUMTYPE_DOUBLE) {
                ctx->columns[i].num_subtype = CIMPORT_NUM_DOUBLE;
            } else if (ctx->columns[i].min_value == DBL_MAX) {
                ctx->columns[i].num_subtype = CIMPORT_NUM_DOUBLE;
            } else {
                ctx->columns[i].num_subtype = cimport_determine_numeric_subtype(&ctx->columns[i]);
            }
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
            cimport_arena_free(&chunk->arena);
        }
        free(ctx->chunks);
    }

    if (ctx->col_cache) {
        for (int i = 0; i < ctx->num_columns; i++) {
            free(ctx->col_cache[i].numeric_data);
            free(ctx->col_cache[i].string_data);
            cimport_arena_free(&ctx->col_cache[i].string_arena);
        }
        free(ctx->col_cache);
    }
    free(ctx->string_arena);

    free(ctx->columns);
    free(ctx->filename);
    free(ctx->row_offsets);
    free(ctx->force_numeric_cols);
    free(ctx->force_string_cols);
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
static CImportContext *cimport_parse_csv(const char *filename, char delimiter, bool has_header, bool verbose, CImportBindQuotesMode bindquotes,
                                          CImportNumericTypeMode numeric_type_mode, char decimal_sep, char group_sep,
                                          CImportEmptyLinesMode emptylines_mode, int max_quoted_rows);
static void cimport_build_column_cache(CImportContext *ctx);

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

static size_t cimport_find_safe_boundary(const char *data, size_t file_size, size_t target, char quote, CImportBindQuotesMode bindquotes) {
    if (target >= file_size) return file_size;

    const char *ptr = data + target;
    const char *start = data;

    /* Find the nearest newline before target */
    while (ptr > start && ptr[-1] != '\n') {
        ptr--;
    }

    /* In loose mode, just use the nearest newline */
    if (bindquotes == CIMPORT_BINDQUOTES_LOOSE) {
        return ptr - data;
    }

    /* In strict mode, check if we're inside a quoted field */
    size_t scan_limit = 10000;
    const char *scan_start = ptr;
    if ((size_t)(ptr - start) > scan_limit) {
        scan_start = ptr - scan_limit;
        while (scan_start > start && scan_start[-1] != '\n') scan_start--;
    } else {
        scan_start = start;
    }

    int quote_count = 0;
    for (const char *p = scan_start; p < ptr; p++) {
        if (*p == quote) quote_count++;
    }

    if (quote_count % 2 == 1) {
        ptr = cimport_find_next_row_strict(ptr, data + file_size, quote);
    }

    return ptr - data;
}

static void *cimport_parse_chunk_parallel(void *arg) {
    CImportChunkParseTask *task = (CImportChunkParseTask *)arg;
    CImportContext *ctx = task->ctx;
    CImportParsedChunk *chunk = task->chunk;

    cimport_arena_init(&chunk->arena);
    chunk->capacity = 16384;
    chunk->rows = ctools_safe_malloc2(chunk->capacity, sizeof(CImportParsedRow *));
    chunk->num_rows = 0;

    chunk->num_col_stats = ctx->num_columns;
    chunk->col_stats = ctools_safe_calloc2(ctx->num_columns, sizeof(CImportColumnParseStats));

    if (!chunk->rows || !chunk->col_stats) {
        atomic_store(&ctx->error_code, 1);
        return NULL;
    }

    CImportFieldRef field_buf[CTOOLS_MAX_COLUMNS];
    const char *ptr = task->start;
    const char *end = task->end;

    while (ptr < end && atomic_load(&ctx->error_code) == 0) {
        const char *row_end = cimport_find_next_row(ptr, end, ctx->quote_char, ctx->bindquotes);
        if (row_end > end) row_end = end;

        int num_fields = cimport_parse_row_fast(ptr, row_end, ctx->delimiter, ctx->quote_char,
                                                 field_buf, CTOOLS_MAX_COLUMNS, ctx->file_data);

        if (num_fields == 0 || (num_fields == 1 && field_buf[0].length == 0)) {
            ptr = row_end;
            continue;
        }

        bool is_header_row = (ctx->has_header && task->chunk_id == 0 && chunk->num_rows == 0);

        /* Check for unmatched quotes in LOOSE mode and record warning */
        if (ctx->bindquotes == CIMPORT_BINDQUOTES_LOOSE && !is_header_row) {
            if (cimport_row_has_unmatched_quote(ptr, row_end, ctx->quote_char)) {
                size_t approx_row = chunk->num_rows + 1;
                if (task->chunk_id > 0) {
                    approx_row += task->chunk_id * (ctx->file_size / ctx->num_chunks / 50);
                }
                cimport_record_unmatched_quote(ctx, approx_row);
            }
        }

        if (!is_header_row) {
            for (int f = 0; f < num_fields && f < chunk->num_col_stats; f++) {
                CImportColumnParseStats *stats = &chunk->col_stats[f];
                CImportFieldRef *field = &field_buf[f];

                if ((int)field->length > stats->max_field_len) {
                    stats->max_field_len = field->length;
                }

                if (!stats->seen_string && field->length > 0) {
                    stats->seen_non_empty = true;
                    const char *src = ctx->file_data + field->offset;
                    if (!cimport_field_looks_numeric(src, field->length)) {
                        stats->seen_string = true;
                    }
                }
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

        if ((size_t)num_fields > (SIZE_MAX - sizeof(CImportParsedRow)) / sizeof(CImportFieldRef)) {
            atomic_store(&ctx->error_code, 1);
            break;
        }
        size_t row_size = sizeof(CImportParsedRow) + sizeof(CImportFieldRef) * num_fields;
        CImportParsedRow *row = cimport_arena_alloc(&chunk->arena, row_size);
        if (!row) {
            atomic_store(&ctx->error_code, 1);
            break;
        }

        row->num_fields = num_fields;
        memcpy(row->fields, field_buf, sizeof(CImportFieldRef) * num_fields);

        chunk->rows[chunk->num_rows++] = row;

        if (num_fields > chunk->max_fields_in_chunk) {
            chunk->max_fields_in_chunk = num_fields;
        }

        ptr = row_end;
    }

    return NULL;
}

/* ============================================================================
 * Main Parse Function
 * ============================================================================ */

static CImportContext *cimport_parse_csv(const char *filename, char delimiter, bool has_header, bool verbose, CImportBindQuotesMode bindquotes,
                                          CImportNumericTypeMode numeric_type_mode, char decimal_sep, char group_sep,
                                          CImportEmptyLinesMode emptylines_mode, int max_quoted_rows) {
    CImportContext *ctx = calloc(1, sizeof(CImportContext));
    if (!ctx) return NULL;

    double t_start, t_end;

    ctx->delimiter = delimiter;
    ctx->quote_char = '"';
    ctx->has_header = has_header;
    ctx->bindquotes = bindquotes;
    ctx->filename = strdup(filename);
    ctx->verbose = verbose;
    atomic_init(&ctx->error_code, 0);

    /* Initialize new options */
    ctx->numeric_type_mode = numeric_type_mode;
    ctx->decimal_separator = decimal_sep;
    ctx->group_separator = group_sep;
    ctx->emptylines_mode = emptylines_mode;
    ctx->max_quoted_rows = max_quoted_rows;
    ctx->force_numeric_cols = NULL;
    ctx->num_force_numeric = 0;
    ctx->force_string_cols = NULL;
    ctx->num_force_string = 0;

    /* Initialize warning tracking */
    ctx->num_unmatched_quote_warnings = 0;
    pthread_mutex_init(&ctx->warning_mutex, NULL);

    t_start = cimport_get_time_ms();
    if (cimport_mmap_file(ctx, filename) != 0) {
        cimport_display_error(ctx->error_message);
        cimport_free_context(ctx);
        return NULL;
    }
    t_end = cimport_get_time_ms();
    ctx->time_mmap = t_end - t_start;

    if (cimport_parse_header(ctx) != 0) {
        cimport_display_error("Failed to parse header");
        cimport_free_context(ctx);
        return NULL;
    }

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

    t_start = cimport_get_time_ms();

    size_t min_chunk_size = 4 * 1024 * 1024;
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
        cimport_arena_init(&chunk->arena);
        chunk->capacity = 65536;
        chunk->rows = ctools_safe_malloc2(chunk->capacity, sizeof(CImportParsedRow *));
        chunk->num_rows = 0;

        chunk->num_col_stats = ctx->num_columns;
        chunk->col_stats = ctools_safe_calloc2(ctx->num_columns, sizeof(CImportColumnParseStats));

        if (!chunk->rows || !chunk->col_stats) {
            cimport_free_context(ctx);
            return NULL;
        }

        CImportFieldRef field_buf[CTOOLS_MAX_COLUMNS];
        const char *ptr = ctx->file_data;
        const char *end = ctx->file_data + ctx->file_size;

        while (ptr < end) {
            const char *row_end = cimport_find_next_row(ptr, end, ctx->quote_char, ctx->bindquotes);

            int num_fields = cimport_parse_row_fast(ptr, row_end, ctx->delimiter, ctx->quote_char,
                                                     field_buf, CTOOLS_MAX_COLUMNS, ctx->file_data);

            if (num_fields == 0 || (num_fields == 1 && field_buf[0].length == 0)) {
                ptr = row_end;
                continue;
            }

            bool is_header_row = (ctx->has_header && chunk->num_rows == 0);

            if (ctx->bindquotes == CIMPORT_BINDQUOTES_LOOSE && !is_header_row) {
                if (cimport_row_has_unmatched_quote(ptr, row_end, ctx->quote_char)) {
                    size_t row_num = chunk->num_rows + 1;
                    cimport_record_unmatched_quote(ctx, row_num);
                }
            }
            if (!is_header_row) {
                for (int f = 0; f < num_fields && f < chunk->num_col_stats; f++) {
                    CImportColumnParseStats *stats = &chunk->col_stats[f];
                    CImportFieldRef *field = &field_buf[f];

                    if ((int)field->length > stats->max_field_len) {
                        stats->max_field_len = field->length;
                    }

                    if (!stats->seen_string && field->length > 0) {
                        stats->seen_non_empty = true;
                        const char *src = ctx->file_data + field->offset;
                        if (!cimport_field_looks_numeric(src, field->length)) {
                            stats->seen_string = true;
                        }
                    }
                }
            }

            if (chunk->num_rows >= chunk->capacity) {
                chunk->capacity *= 2;
                CImportParsedRow **new_rows = realloc(chunk->rows, sizeof(CImportParsedRow *) * chunk->capacity);
                if (!new_rows) {
                    cimport_free_context(ctx);
                    return NULL;
                }
                chunk->rows = new_rows;
            }

            size_t row_size = sizeof(CImportParsedRow) + sizeof(CImportFieldRef) * num_fields;
            CImportParsedRow *row = cimport_arena_alloc(&chunk->arena, row_size);
            if (!row) {
                cimport_free_context(ctx);
                return NULL;
            }

            row->num_fields = num_fields;
            memcpy(row->fields, field_buf, sizeof(CImportFieldRef) * num_fields);
            chunk->rows[chunk->num_rows++] = row;

            if (num_fields > chunk->max_fields_in_chunk) {
                chunk->max_fields_in_chunk = num_fields;
            }

            ptr = row_end;
        }
    } else {
        size_t chunk_size = ctx->file_size / num_chunks;

        size_t *boundaries = ctools_safe_malloc2((size_t)num_chunks + 1, sizeof(size_t));
        if (!boundaries) {
            cimport_free_context(ctx);
            return NULL;
        }

        boundaries[0] = 0;
        for (int i = 1; i < num_chunks; i++) {
            size_t target = i * chunk_size;
            boundaries[i] = cimport_find_safe_boundary(ctx->file_data, ctx->file_size, target, ctx->quote_char, ctx->bindquotes);
        }
        boundaries[num_chunks] = ctx->file_size;

        pthread_t threads[CIMPORT_MAX_THREADS];
        CImportChunkParseTask tasks[CIMPORT_MAX_THREADS];

        for (int i = 0; i < num_chunks; i++) {
            tasks[i].ctx = ctx;
            tasks[i].chunk_id = i;
            tasks[i].start = ctx->file_data + boundaries[i];
            tasks[i].end = ctx->file_data + boundaries[i + 1];
            tasks[i].chunk = &ctx->chunks[i];

            pthread_create(&threads[i], NULL, cimport_parse_chunk_parallel, &tasks[i]);
        }

        for (int i = 0; i < num_chunks; i++) {
            pthread_join(threads[i], NULL);
        }

        free(boundaries);

        if (atomic_load(&ctx->error_code) != 0) {
            cimport_free_context(ctx);
            return NULL;
        }
    }

    t_end = cimport_get_time_ms();
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
    }

    t_start = cimport_get_time_ms();
    cimport_infer_column_types(ctx);
    t_end = cimport_get_time_ms();
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
    char field_buf[CTOOLS_MAX_STRING_LEN + 1];

    for (int col_idx = task->col_start; col_idx < task->col_end; col_idx++) {
        CImportColumnInfo *col = &ctx->columns[col_idx];
        CImportColumnCache *cache = &ctx->col_cache[col_idx];
        size_t row_idx = 0;

        if (col->type == CIMPORT_COL_STRING) {
            if (cache->string_data == NULL) {
                cache->count = 0;
                continue;
            }

            cimport_arena_init(&cache->string_arena);

            for (int c = 0; c < ctx->num_chunks; c++) {
                CImportParsedChunk *chunk = &ctx->chunks[c];
                size_t start_row = (c == 0 && ctx->has_header) ? 1 : 0;

                for (size_t r = start_row; r < chunk->num_rows; r++) {
                    CImportParsedRow *row = chunk->rows[r];
                    char *str;

                    if (col_idx < row->num_fields) {
                        CImportFieldRef *field = &row->fields[col_idx];
                        const char *src = ctx->file_data + field->offset;

                        int len;
                        if (!cimport_field_contains_quote(src, field->length, ctx->quote_char)) {
                            len = cimport_extract_field_unquoted(ctx->file_data, field, field_buf, CTOOLS_MAX_STRING_LEN);
                        } else {
                            len = cimport_extract_field_fast(ctx->file_data, field, field_buf, CTOOLS_MAX_STRING_LEN, ctx->quote_char);
                        }
                        field_buf[len] = '\0';

                        str = cimport_arena_alloc(&cache->string_arena, len + 1);
                        if (str) {
                            memcpy(str, field_buf, len + 1);
                        } else {
                            str = "";
                        }
                    } else {
                        str = cimport_arena_alloc(&cache->string_arena, 1);
                        if (str) {
                            str[0] = '\0';
                        } else {
                            str = "";
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
                        CImportFieldRef *field = &row->fields[col_idx];
                        const char *src = ctx->file_data + field->offset;

                        if (!ctools_parse_double_fast(src, field->length, &val, missing)) {
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

    return NULL;
}

static void cimport_build_column_cache(CImportContext *ctx) {
    if (ctx->cache_ready) return;

    double t_start = cimport_get_time_ms();

    ctx->col_cache = ctools_safe_calloc2(ctx->num_columns, sizeof(CImportColumnCache));
    if (!ctx->col_cache) return;

    for (int i = 0; i < ctx->num_columns; i++) {
        CImportColumnInfo *col = &ctx->columns[i];
        CImportColumnCache *cache = &ctx->col_cache[i];

        if (col->type == CIMPORT_COL_STRING) {
            cache->string_data = ctools_safe_calloc2(ctx->total_rows, sizeof(char *));
            cache->numeric_data = NULL;
        } else {
            cache->numeric_data = ctools_safe_malloc2(ctx->total_rows, sizeof(double));
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
        pthread_t threads[CIMPORT_MAX_THREADS];
        CImportCacheBuildTask tasks[CIMPORT_MAX_THREADS];

        int cols_per_thread = (ctx->num_columns + num_threads - 1) / num_threads;
        int threads_created = 0;

        for (int t = 0; t < num_threads; t++) {
            tasks[t].ctx = ctx;
            tasks[t].col_start = t * cols_per_thread;
            tasks[t].col_end = (t + 1) * cols_per_thread;
            if (tasks[t].col_end > ctx->num_columns) tasks[t].col_end = ctx->num_columns;
            tasks[t].thread_id = t;

            if (pthread_create(&threads[t], NULL, cimport_build_cache_worker, &tasks[t]) == 0) {
                threads_created++;
            } else {
                cimport_build_cache_worker(&tasks[t]);
            }
        }

        for (int t = 0; t < threads_created; t++) {
            pthread_join(threads[t], NULL);
        }
    }

    ctx->cache_ready = true;

    double t_end = cimport_get_time_ms();
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

    /* Parse tokens */
    p = buf;
    int idx = 0;
    while (*p && idx < count) {
        while (*p == ' ') p++;
        if (*p == '\0') break;
        cols[idx++] = atoi(p);
        while (*p && *p != ' ') p++;
    }

    *out_count = idx;
    return cols;
}

/* ============================================================================
 * SCAN Mode
 * ============================================================================ */

static ST_retcode cimport_do_scan(const char *filename, char delimiter, bool has_header, bool verbose, CImportBindQuotesMode bindquotes,
                                   CImportNumericTypeMode numeric_type_mode, char decimal_sep, char group_sep,
                                   CImportEmptyLinesMode emptylines_mode, int max_quoted_rows) {
    char msg[512];
    char macro_val[65536];

    cimport_clear_cached_context();

    g_cimport_ctx = cimport_parse_csv(filename, delimiter, has_header, verbose, bindquotes,
                                       numeric_type_mode, decimal_sep, group_sep, emptylines_mode, max_quoted_rows);

    /* Parse column type overrides from global macros */
    if (g_cimport_ctx) {
        g_cimport_ctx->force_numeric_cols = cimport_parse_col_list("_CIMPORT_NUMCOLS", &g_cimport_ctx->num_force_numeric);
        g_cimport_ctx->force_string_cols = cimport_parse_col_list("_CIMPORT_STRCOLS", &g_cimport_ctx->num_force_string);
    }
    if (!g_cimport_ctx) {
        return 601;
    }

    CImportContext *ctx = g_cimport_ctx;

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

    /* Display warnings for malformed data (like Stata) */
    cimport_display_warnings(ctx);

    snprintf(msg, sizeof(msg), "Scanned %zu rows, %d columns\n", ctx->total_rows, ctx->num_columns);
    cimport_display_msg(msg);

    if (verbose) {
        snprintf(msg, sizeof(msg), "  [C timing] mmap: %.2f ms, parse: %.2f ms, type_infer: %.2f ms\n",
                 ctx->time_mmap, ctx->time_parse, ctx->time_type_infer);
        cimport_display_msg(msg);
    }

    return 0;
}

/* ============================================================================
 * LOAD Mode
 * ============================================================================ */

static ST_retcode cimport_do_load(const char *filename, char delimiter, bool has_header, bool verbose, CImportBindQuotesMode bindquotes,
                                   CImportNumericTypeMode numeric_type_mode, char decimal_sep, char group_sep,
                                   CImportEmptyLinesMode emptylines_mode, int max_quoted_rows) {
    char msg[512];

    CImportContext *ctx = g_cimport_ctx;
    bool used_cache = (ctx != NULL && ctx->is_loaded && ctx->filename && strcmp(ctx->filename, filename) == 0);

    if (!used_cache) {
        cimport_clear_cached_context();
        g_cimport_ctx = cimport_parse_csv(filename, delimiter, has_header, verbose, bindquotes,
                                           numeric_type_mode, decimal_sep, group_sep, emptylines_mode, max_quoted_rows);
        ctx = g_cimport_ctx;
        if (!ctx) {
            return 601;
        }
        /* Parse column type overrides */
        ctx->force_numeric_cols = cimport_parse_col_list("_CIMPORT_NUMCOLS", &ctx->num_force_numeric);
        ctx->force_string_cols = cimport_parse_col_list("_CIMPORT_STRCOLS", &ctx->num_force_string);
        cimport_display_warnings(ctx);
    }

    if (!ctx->cache_ready) {
        cimport_build_column_cache(ctx);
        if (!ctx->cache_ready) {
            cimport_display_error("Failed to build column cache\n");
            cimport_clear_cached_context();
            return 459;
        }
    }

    double t_load_start = cimport_get_time_ms();
    long total_fields = 0;

    ST_int stata_nvars = SF_nvars();
    ST_int stata_nobs = SF_nobs();

    if (verbose) {
        snprintf(msg, sizeof(msg), "  [debug] SF_nvars()=%d, SF_nobs()=%d, ctx->num_columns=%d, ctx->total_rows=%zu\n",
                 stata_nvars, stata_nobs, ctx->num_columns, ctx->total_rows);
        cimport_display_msg(msg);
    }

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

    for (int col_idx = 0; col_idx < ctx->num_columns; col_idx++) {
        CImportColumnInfo *col = &ctx->columns[col_idx];
        CImportColumnCache *cache = &ctx->col_cache[col_idx];
        ST_int var = col_idx + 1;
        size_t nrows = cache->count;

        if (col->type == CIMPORT_COL_STRING) {
            char **strings = cache->string_data;

            for (size_t p = 0; p < 16 && p < nrows; p++) {
                __builtin_prefetch(strings[p], 0, 0);
            }

            size_t i = 0;
            size_t nrows_aligned = nrows & ~(size_t)3;

            for (; i < nrows_aligned; i += 4) {
                if (i + 16 < nrows) {
                    __builtin_prefetch(&strings[i + 16], 0, 0);
                    __builtin_prefetch(strings[i + 12], 0, 0);
                    __builtin_prefetch(strings[i + 13], 0, 0);
                    __builtin_prefetch(strings[i + 14], 0, 0);
                    __builtin_prefetch(strings[i + 15], 0, 0);
                }

                SF_sstore(var, (ST_int)(i + 1), strings[i]);
                SF_sstore(var, (ST_int)(i + 2), strings[i + 1]);
                SF_sstore(var, (ST_int)(i + 3), strings[i + 2]);
                SF_sstore(var, (ST_int)(i + 4), strings[i + 3]);
            }

            for (; i < nrows; i++) {
                SF_sstore(var, (ST_int)(i + 1), strings[i]);
            }
            total_fields += nrows;
        } else {
            double *values = cache->numeric_data;

            __builtin_prefetch(values, 0, 0);
            __builtin_prefetch(values + 8, 0, 0);
            __builtin_prefetch(values + 16, 0, 0);
            __builtin_prefetch(values + 24, 0, 0);

            size_t i = 0;
            size_t nrows_aligned = nrows & ~(size_t)7;

            for (; i < nrows_aligned; i += 8) {
                __builtin_prefetch(values + i + 32, 0, 0);
                __builtin_prefetch(values + i + 40, 0, 0);

                SF_vstore(var, (ST_int)(i + 1), values[i]);
                SF_vstore(var, (ST_int)(i + 2), values[i + 1]);
                SF_vstore(var, (ST_int)(i + 3), values[i + 2]);
                SF_vstore(var, (ST_int)(i + 4), values[i + 3]);
                SF_vstore(var, (ST_int)(i + 5), values[i + 4]);
                SF_vstore(var, (ST_int)(i + 6), values[i + 5]);
                SF_vstore(var, (ST_int)(i + 7), values[i + 6]);
                SF_vstore(var, (ST_int)(i + 8), values[i + 7]);
            }

            for (; i < nrows; i++) {
                SF_vstore(var, (ST_int)(i + 1), values[i]);
            }
            total_fields += nrows;
        }
    }

    double t_load_end = cimport_get_time_ms();
    double time_spi = t_load_end - t_load_start;

    snprintf(msg, sizeof(msg), "Loaded %zu observations\n", ctx->total_rows);
    cimport_display_msg(msg);

    if (verbose) {
        snprintf(msg, sizeof(msg), "  [C timing] SPI store: %.2f ms (%ld fields, %.0f fields/sec)\n",
                 time_spi, total_fields, total_fields / (time_spi / 1000.0));
        cimport_display_msg(msg);
    }

    /* Save timing and thread diagnostics to Stata scalars */
    SF_scal_save("_cimport_time_mmap", ctx->time_mmap);
    SF_scal_save("_cimport_time_parse", ctx->time_parse);
    SF_scal_save("_cimport_time_infer", ctx->time_type_infer);
    SF_scal_save("_cimport_time_cache", ctx->time_cache);
    SF_scal_save("_cimport_time_store", time_spi);
    SF_scal_save("_cimport_time_total", ctx->time_mmap + ctx->time_parse + ctx->time_type_infer + ctx->time_cache + time_spi);
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
    bool verbose;
    CImportBindQuotesMode bindquotes;
    CImportNumericTypeMode numeric_type_mode;
    char decimal_separator;
    char group_separator;
    CImportEmptyLinesMode emptylines_mode;
    int max_quoted_rows;
} CImportOptions;

ST_retcode cimport_main(const char *args) {
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
        .delimiter = ',',
        .has_header = true,
        .verbose = false,
        .bindquotes = CIMPORT_BINDQUOTES_LOOSE,
        .numeric_type_mode = CIMPORT_NUMTYPE_AUTO,
        .decimal_separator = '.',
        .group_separator = '\0',
        .emptylines_mode = CIMPORT_EMPTYLINES_SKIP,
        .max_quoted_rows = 20
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
            } else if (strcmp(token, "verbose") == 0) {
                opts.verbose = true;
            } else if (strcmp(token, "tab") == 0) {
                opts.delimiter = '\t';
            } else if (strcmp(token, "bindquotes=strict") == 0) {
                opts.bindquotes = CIMPORT_BINDQUOTES_STRICT;
            } else if (strcmp(token, "bindquotes=loose") == 0) {
                opts.bindquotes = CIMPORT_BINDQUOTES_LOOSE;
            } else if (strcmp(token, "asfloat") == 0) {
                opts.numeric_type_mode = CIMPORT_NUMTYPE_FLOAT;
            } else if (strcmp(token, "asdouble") == 0) {
                opts.numeric_type_mode = CIMPORT_NUMTYPE_DOUBLE;
            } else if (strncmp(token, "decimalsep=", 11) == 0) {
                opts.decimal_separator = token[11];
            } else if (strncmp(token, "groupsep=", 9) == 0) {
                opts.group_separator = token[9];
            } else if (strcmp(token, "emptylines=fill") == 0) {
                opts.emptylines_mode = CIMPORT_EMPTYLINES_FILL;
            } else if (strncmp(token, "maxquotedrows=", 14) == 0) {
                opts.max_quoted_rows = atoi(token + 14);
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

    if (strcmp(opts.mode, "scan") == 0) {
        return cimport_do_scan(opts.filename, opts.delimiter, opts.has_header, opts.verbose, opts.bindquotes,
                                opts.numeric_type_mode, opts.decimal_separator, opts.group_separator,
                                opts.emptylines_mode, opts.max_quoted_rows);
    } else if (strcmp(opts.mode, "load") == 0) {
        return cimport_do_load(opts.filename, opts.delimiter, opts.has_header, opts.verbose, opts.bindquotes,
                                opts.numeric_type_mode, opts.decimal_separator, opts.group_separator,
                                opts.emptylines_mode, opts.max_quoted_rows);
    } else {
        cimport_display_error("cimport: invalid mode. Use 'scan' or 'load'\n");
        return 198;
    }
}

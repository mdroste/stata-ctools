/*
 * cimport_context.h
 * Internal context structure for cimport
 *
 * Contains the main CImportContext struct and related types.
 * This header is shared between cimport_impl.c and helper modules.
 */

#ifndef CIMPORT_CONTEXT_H
#define CIMPORT_CONTEXT_H

#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdatomic.h>
#include <pthread.h>

#include "cimport_arena.h"
#include "cimport_parse.h"

#ifdef _WIN32
#include <windows.h>
#define CIMPORT_WINDOWS 1
#endif

/* Warning tracking for malformed data */
#define CIMPORT_MAX_WARNINGS 10

/* Maximum variable name length (from ctools_config.h) */
#ifndef CTOOLS_MAX_VARNAME_LEN
#define CTOOLS_MAX_VARNAME_LEN 128
#endif

/* ============================================================================
 * Type Definitions
 * ============================================================================ */

typedef enum {
    CIMPORT_COL_UNKNOWN = 0,
    CIMPORT_COL_NUMERIC,
    CIMPORT_COL_STRING
} CImportColumnType;

typedef enum {
    CIMPORT_NUM_DOUBLE = 0,
    CIMPORT_NUM_FLOAT,
    CIMPORT_NUM_LONG,
    CIMPORT_NUM_INT,
    CIMPORT_NUM_BYTE
} CImportNumericSubtype;

typedef struct {
    char name[CTOOLS_MAX_VARNAME_LEN];
    CImportColumnType type;
    int max_strlen;
    CImportNumericSubtype num_subtype;
    bool is_integer;
    double min_value;
    double max_value;
    int max_decimal_digits;
} CImportColumnInfo;

typedef struct {
    bool seen_string;
    bool seen_non_empty;
    int max_field_len;
} CImportColumnParseStats;

typedef struct {
    CImportParsedRow **rows;
    size_t num_rows;
    size_t capacity;
    CImportArena arena;
    CImportColumnParseStats *col_stats;
    int num_col_stats;
    int max_fields_in_chunk;
} CImportParsedChunk;

typedef struct {
    double *numeric_data;
    char **string_data;
    CImportArena string_arena;
    size_t count;
} CImportColumnCache;

/* Main context structure */
typedef struct CImportContext {
    char *file_data;
    size_t file_size;
    char *filename;
    CImportColumnInfo *columns;
    int num_columns;
    size_t total_rows;
    char delimiter;
    char quote_char;
    bool has_header;
    CImportBindQuotesMode bindquotes;
    CImportParsedChunk *chunks;
    int num_chunks;
    size_t *row_offsets;
    CImportColumnCache *col_cache;
    char *string_arena;
    size_t string_arena_size;
    size_t string_arena_used;
    int num_threads;
    atomic_int error_code;
    char error_message[256];
    bool is_loaded;
    bool cache_ready;
    bool verbose;
    int max_fields_seen;
    double time_mmap;
    double time_parse;
    double time_type_infer;
    double time_cache;
    /* Warning tracking for unmatched quotes */
    size_t unmatched_quote_rows[CIMPORT_MAX_WARNINGS];
    int num_unmatched_quote_warnings;
    pthread_mutex_t warning_mutex;
#ifdef CIMPORT_WINDOWS
    HANDLE file_handle;
    HANDLE mapping_handle;
#endif
} CImportContext;

#endif /* CIMPORT_CONTEXT_H */

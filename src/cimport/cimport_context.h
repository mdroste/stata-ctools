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

#include "../ctools_arena.h"
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
    ctools_arena arena;
    CImportColumnParseStats *col_stats;
    int num_col_stats;
    int max_fields_in_chunk;
} CImportParsedChunk;

typedef struct {
    double *numeric_data;
    char **string_data;
    ctools_arena string_arena;
    size_t count;
} CImportColumnCache;

/* Numeric type forcing mode */
typedef enum {
    CIMPORT_NUMTYPE_AUTO = 0,   /* Auto-detect optimal type */
    CIMPORT_NUMTYPE_FLOAT,      /* Force all numerics to float */
    CIMPORT_NUMTYPE_DOUBLE      /* Force all numerics to double */
} CImportNumericTypeMode;

/* Empty line handling mode */
typedef enum {
    CIMPORT_EMPTYLINES_SKIP = 0,  /* Skip empty lines (default) */
    CIMPORT_EMPTYLINES_FILL       /* Include empty lines as missing */
} CImportEmptyLinesMode;

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

    /* New options */
    CImportNumericTypeMode numeric_type_mode;  /* asfloat/asdouble */
    char decimal_separator;      /* decimalseparator (default '.') */
    char group_separator;        /* groupseparator (default '\0' = none) */
    CImportEmptyLinesMode emptylines_mode;  /* emptylines handling */
    int max_quoted_rows;         /* maxquotedrows for type inference */

    /* Column type overrides (NULL = no override) */
    int *force_numeric_cols;     /* Array of 1-based column indices to force numeric */
    int num_force_numeric;
    int *force_string_cols;      /* Array of 1-based column indices to force string */
    int num_force_string;

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

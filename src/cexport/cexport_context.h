/*
 * cexport_context.h
 * Internal context structure for cexport
 *
 * Contains the main cexport_context struct and related types.
 * This header is shared between cexport_impl.c and helper modules.
 */

#ifndef CEXPORT_CONTEXT_H
#define CEXPORT_CONTEXT_H

#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>

#include "ctools_types.h"
#include "ctools_arena.h"
#include "cexport_io.h"

/* ========================================================================
   Variable Storage Types (from Stata)
   ======================================================================== */

typedef enum {
    VARTYPE_STRING = 0,
    VARTYPE_BYTE   = 1,
    VARTYPE_INT    = 2,
    VARTYPE_LONG   = 3,
    VARTYPE_FLOAT  = 4,
    VARTYPE_DOUBLE = 5
} vartype_t;

/* ========================================================================
   Export Context and Configuration
   ======================================================================== */

typedef struct {
    /* Output file */
    char *filename;
    FILE *fp;                    /* Legacy: only used for fallback single-threaded path */

    /* Formatting options */
    char delimiter;
    bool write_header;
    bool quote_strings;
    bool quote_if_needed;
    bool verbose;

    /* Line ending configuration (for binary mode with explicit endings) */
    char line_ending[3];         /* "\n" or "\r\n" + null terminator */
    int line_ending_len;         /* 1 for LF, 2 for CRLF */
    bool use_crlf;               /* True for Windows-style CRLF */

    /* I/O backend configuration */
    cexport_io_backend_t io_backend;  /* PWRITE or MMAP */
    cexport_io_flags_t io_flags;      /* I/O flags (DIRECT, NOFSYNC, PREFAULT) */
    bool use_parallel_io;        /* True to use parallel offset writes */

    /* Adaptive sizing (computed at runtime) */
    size_t actual_avg_row_size;  /* Sampled average row size */
    size_t adaptive_chunk_size;  /* Rows per chunk (adaptive) */

    /* Data from Stata */
    ctools_filtered_data filtered;

    /* Observation range for if/in filtering */
    size_t obs1;    /* First observation (1-based Stata index) */
    size_t nobs_loaded;  /* Number of observations loaded */

    /* Variable names (for header) */
    char **varnames;
    size_t nvars;

    /* Variable storage types (for numeric precision) */
    vartype_t *vartypes;

    /* Optimization flags */
    bool all_numeric;    /* True if all variables are numeric (fast path) */

    /* Arena for chunk metadata (format/write args, offsets) */
    ctools_arena chunk_arena;

    /* Timing */
    double time_load;
    double time_format;
    double time_write;
    double time_total;
} cexport_context;

/* ========================================================================
   Thread Argument Structures
   ======================================================================== */

/*
    Thread arguments for parallel row formatting.
*/
typedef struct {
    size_t start_row;        /* First row to format (0-based) */
    size_t end_row;          /* Last row to format (exclusive) */
    char *output_buffer;     /* Output buffer for formatted rows */
    size_t buffer_size;      /* Size of output buffer */
    size_t bytes_written;    /* Actual bytes written to buffer */
    int success;             /* 1 on success, 0 on failure */
} format_chunk_args_t;

/*
    Thread arguments for parallel offset writes.
*/
typedef struct {
    cexport_io_file *file;   /* Output file handle */
    const char *buffer;      /* Data to write */
    size_t len;              /* Number of bytes to write */
    size_t offset;           /* File offset to write at */
    int success;             /* 1 on success, 0 on failure */
} write_chunk_args_t;

/*
    Thread arguments for zero-copy mmap formatting.
    Each thread formats directly into mapped_data + offset.
*/
typedef struct {
    char *dest;              /* Destination pointer in mapped memory */
    size_t dest_size;        /* Available space at destination */
    size_t start_row;        /* First row to format (0-based) */
    size_t end_row;          /* Last row to format (exclusive) */
    size_t bytes_written;    /* Actual bytes written */
    int success;             /* 1 on success, 0 on failure */
} format_mmap_args_t;

#endif /* CEXPORT_CONTEXT_H */

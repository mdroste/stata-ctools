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
    VARTYPE_DOUBLE = 5,
    VARTYPE_FORMATTED = 6
} vartype_t;

/* How the cells of a numeric variable are written */
typedef enum {
    CEXPORT_CELL_PLAIN = 0, /* the number */
    CEXPORT_CELL_TIME  = 1, /* in its %t display format, never quoted */
    CEXPORT_CELL_LABEL = 2  /* its value label's text; the number when unlabeled */
} cexport_cell_t;

/* A compiled %t display format (cexport_format.c) */
typedef struct cexport_timefmt cexport_timefmt;

/* ========================================================================
   Export Context and Configuration
   ======================================================================== */

typedef struct {
    /* Output file */
    char *filename;
    FILE *fp;                    /* Legacy: only used for fallback single-threaded path */

    /* Formatting options */
    char delimiter;
    bool replace;             /* Atomic destination replacement allowed */
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
    cexport_io_flags_t io_flags;      /* I/O flags (DIRECT, PREFAULT) */
    bool use_parallel_io;        /* True to use parallel offset writes */
    bool skip_if;                /* No if condition: select the in-range */
    bool touse;                  /* Last plugin variable is a sample marker */

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

    /* Display handling of numeric variables, [nvars] each; all NULL when
     * every cell is the plain number. timefmt[j] and labelmap[j] are owned
     * and set for the matching cell_kind[j]; cell_reserve[j] bounds the
     * bytes such a cell can take. */
    unsigned char *cell_kind;     /* cexport_cell_t */
    cexport_timefmt **timefmt;
    void **labelmap;              /* cdecode label maps */
    size_t *cell_reserve;
    /* Set with them by their owner (cexport_impl.c): the label of a value
     * in a map (NULL when it has none), and the release of all of the
     * above, which cexport_context_cleanup() calls. */
    const char *(*label_text)(const void *map, double value);
    void (*cells_free)(void *ctx);

    /* Optimization flags */
    bool all_numeric;    /* True if all variables are numeric (fast path) */
    bool stream;         /* Format cells read directly from Stata (no loaded copy) */
    unsigned char *is_string;  /* Per variable, for streaming */

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
    char **buffer_slot;      /* Owning slot, updated when the buffer grows */
    size_t *size_slot;
    size_t bytes_written;    /* Actual bytes written to buffer */
    int success;             /* 1 on success, 0 on failure */
    int read_error;          /* Stata data could not be read */
    int row_too_long;        /* A row exceeded the formatters' INT_MAX limit */
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

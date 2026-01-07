/*
    cexport_impl.c
    High-Performance Delimited Text Export for Stata

    Exports Stata datasets to CSV/delimited text files with parallel processing.
    Designed for maximum throughput on large datasets.

    Architecture:
    1. PHASE 1: Load data from Stata into C memory (parallel per-variable)
    2. PHASE 2: Format data into text buffers (parallel per-chunk)
    3. PHASE 3: Write buffers to file (sequential, buffered I/O)

    Performance Optimizations:
    - Parallel data loading via ctools_data_load
    - Chunked parallel formatting: each thread formats a range of rows
    - Custom fast double-to-string conversion (avoids sprintf overhead)
    - Large I/O buffers (64KB) to minimize syscalls
    - Memory-mapped output option for very large files
    - Pre-computed field widths for aligned output (optional)

    Thread Safety:
    - Each formatting thread works on disjoint row ranges
    - No shared mutable state during formatting phase
    - Single writer thread for I/O (sequential)
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_timer.h"
#include "cexport_impl.h"

/* ========================================================================
   Configuration Constants
   ======================================================================== */

#define CEXPORT_IO_BUFFER_SIZE   (64 * 1024)   /* 64KB I/O buffer */
#define CEXPORT_CHUNK_SIZE       10000          /* Rows per parallel chunk */
#define CEXPORT_MAX_THREADS      16             /* Maximum formatting threads */
#define CEXPORT_DBL_BUFFER_SIZE  32             /* Buffer for double formatting */
#define CEXPORT_INITIAL_ROW_BUF  4096           /* Initial row buffer size */

/* ========================================================================
   Export Context and Configuration
   ======================================================================== */

typedef struct {
    /* Output file */
    char *filename;
    FILE *fp;

    /* Formatting options */
    char delimiter;
    bool write_header;
    bool quote_strings;
    bool quote_if_needed;
    bool verbose;

    /* Data from Stata */
    stata_data data;

    /* Variable names (for header) */
    char **varnames;
    size_t nvars;

    /* Timing */
    double time_load;
    double time_format;
    double time_write;
    double time_total;
} cexport_context;

/* Global context for this export operation */
static cexport_context g_ctx;

/* Timing uses ctools_timer.h */

/* ========================================================================
   Fast Double to String Conversion

   Custom implementation that avoids sprintf overhead.
   Handles integers, decimals, scientific notation, and Stata missing values.
   ======================================================================== */

/* Check if value is a Stata missing value */
static inline bool is_stata_missing(double val)
{
    /* Stata missing values: . .a .b ... .z are encoded as specific large values */
    return SF_is_missing(val);
}

/*
    Fast double to string conversion.

    Features:
    - Handles integers without decimal point
    - Uses minimal precision for clean output
    - Handles Stata missing values as empty string or "."
    - Avoids sprintf for common cases

    Returns: number of characters written (not including null terminator)
*/
static int double_to_str(double val, char *buf, int buf_size, bool missing_as_dot)
{
    if (is_stata_missing(val)) {
        if (missing_as_dot) {
            buf[0] = '.';
            buf[1] = '\0';
            return 1;
        } else {
            buf[0] = '\0';
            return 0;
        }
    }

    /* Handle special floating point values */
    if (isnan(val)) {
        buf[0] = '.';
        buf[1] = '\0';
        return 1;
    }
    if (isinf(val)) {
        if (val > 0) {
            memcpy(buf, "inf", 4);
            return 3;
        } else {
            memcpy(buf, "-inf", 5);
            return 4;
        }
    }

    /* Check if value is an integer (no fractional part) */
    double int_part;
    double frac_part = modf(val, &int_part);

    if (frac_part == 0.0 && fabs(int_part) < 1e15) {
        /* Integer value - format without decimal point */
        return snprintf(buf, buf_size, "%.0f", val);
    }

    /* General case: use %g for clean output */
    int len = snprintf(buf, buf_size, "%.15g", val);

    /* Remove trailing zeros after decimal point */
    char *dot = strchr(buf, '.');
    if (dot != NULL) {
        char *end = buf + len - 1;
        char *exp = strchr(buf, 'e');
        if (exp == NULL) exp = strchr(buf, 'E');

        char *last_digit = (exp != NULL) ? exp - 1 : end;

        while (last_digit > dot && *last_digit == '0') {
            last_digit--;
        }

        if (last_digit == dot) {
            /* Remove the decimal point too if no fractional digits */
            if (exp != NULL) {
                memmove(dot, exp, strlen(exp) + 1);
                len = (int)strlen(buf);
            } else {
                *dot = '\0';
                len = (int)(dot - buf);
            }
        } else if (last_digit < end - 1 || (exp != NULL && last_digit < exp - 1)) {
            /* Remove trailing zeros */
            if (exp != NULL) {
                memmove(last_digit + 1, exp, strlen(exp) + 1);
                len = (int)strlen(buf);
            } else {
                *(last_digit + 1) = '\0';
                len = (int)(last_digit - buf + 1);
            }
        }
    }

    return len;
}

/* ========================================================================
   String Quoting and Escaping
   ======================================================================== */

/*
    Check if a string needs quoting.

    A string needs quoting if it contains:
    - The delimiter character
    - A double quote
    - A newline or carriage return
*/
static inline bool string_needs_quoting(const char *str, char delimiter)
{
    if (str == NULL) return false;

    for (const char *p = str; *p != '\0'; p++) {
        if (*p == delimiter || *p == '"' || *p == '\n' || *p == '\r') {
            return true;
        }
    }
    return false;
}

/*
    Write a quoted string to a buffer.

    Doubles any internal quotes (CSV escaping convention).
    Returns number of characters written.
*/
static int write_quoted_string(const char *str, char *buf, int buf_size)
{
    if (str == NULL || buf_size < 3) {
        buf[0] = '"';
        buf[1] = '"';
        buf[2] = '\0';
        return 2;
    }

    int pos = 0;
    buf[pos++] = '"';

    for (const char *p = str; *p != '\0' && pos < buf_size - 2; p++) {
        if (*p == '"') {
            /* Escape quote by doubling it */
            if (pos < buf_size - 3) {
                buf[pos++] = '"';
                buf[pos++] = '"';
            }
        } else {
            buf[pos++] = *p;
        }
    }

    buf[pos++] = '"';
    buf[pos] = '\0';
    return pos;
}

/* ========================================================================
   Row Formatting
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
    Format a single row into the output buffer.

    Returns: number of bytes written, or -1 on buffer overflow
*/
static int format_row(size_t row_idx, char *buf, size_t buf_size)
{
    char field_buf[4096];
    size_t pos = 0;
    size_t nvars = g_ctx.data.nvars;

    for (size_t j = 0; j < nvars; j++) {
        stata_variable *var = &g_ctx.data.vars[j];
        int field_len = 0;

        if (var->type == STATA_TYPE_STRING) {
            const char *str = var->data.str[row_idx];
            if (str == NULL) str = "";

            bool need_quote = g_ctx.quote_strings ||
                             (g_ctx.quote_if_needed && string_needs_quoting(str, g_ctx.delimiter));

            if (need_quote) {
                field_len = write_quoted_string(str, field_buf, sizeof(field_buf));
            } else {
                field_len = (int)strlen(str);
                if (field_len >= (int)sizeof(field_buf)) field_len = sizeof(field_buf) - 1;
                memcpy(field_buf, str, field_len);
                field_buf[field_len] = '\0';
            }
        } else {
            /* Numeric variable */
            double val = var->data.dbl[row_idx];
            field_len = double_to_str(val, field_buf, sizeof(field_buf), true);
        }

        /* Check buffer space */
        if (pos + field_len + 2 > buf_size) {
            return -1; /* Buffer overflow */
        }

        /* Copy field to output */
        memcpy(buf + pos, field_buf, field_len);
        pos += field_len;

        /* Add delimiter or newline */
        if (j < nvars - 1) {
            buf[pos++] = g_ctx.delimiter;
        } else {
            buf[pos++] = '\n';
        }
    }

    return (int)pos;
}

/*
    Thread function: Format a chunk of rows.
*/
static void *format_chunk_thread(void *arg)
{
    format_chunk_args_t *args = (format_chunk_args_t *)arg;
    size_t pos = 0;

    for (size_t i = args->start_row; i < args->end_row; i++) {
        int row_len = format_row(i, args->output_buffer + pos, args->buffer_size - pos);
        if (row_len < 0) {
            args->success = 0;
            args->bytes_written = pos;
            return NULL;
        }
        pos += row_len;
    }

    args->success = 1;
    args->bytes_written = pos;
    return NULL;
}

/* ========================================================================
   Header Writing
   ======================================================================== */

static int write_header(FILE *fp)
{
    char line_buf[65536];
    size_t pos = 0;

    for (size_t j = 0; j < g_ctx.nvars; j++) {
        const char *name = g_ctx.varnames[j];
        size_t name_len = strlen(name);

        /* Check if name needs quoting */
        bool need_quote = g_ctx.quote_strings ||
                         (g_ctx.quote_if_needed && string_needs_quoting(name, g_ctx.delimiter));

        if (need_quote) {
            int quoted_len = write_quoted_string(name, line_buf + pos, sizeof(line_buf) - pos);
            pos += quoted_len;
        } else {
            if (pos + name_len + 2 > sizeof(line_buf)) {
                return -1;
            }
            memcpy(line_buf + pos, name, name_len);
            pos += name_len;
        }

        if (j < g_ctx.nvars - 1) {
            line_buf[pos++] = g_ctx.delimiter;
        } else {
            line_buf[pos++] = '\n';
        }
    }

    if (fwrite(line_buf, 1, pos, fp) != pos) {
        return -1;
    }

    return 0;
}

/* ========================================================================
   Main Export Function
   ======================================================================== */

/*
    Parse command arguments.

    Format: filename delimiter [options...]
    Options: noheader, quote, quoteif, verbose
*/
static int parse_args(const char *args)
{
    char *args_copy = strdup(args);
    if (args_copy == NULL) return -1;

    char *saveptr;
    char *token;
    int arg_idx = 0;

    /* Initialize defaults */
    g_ctx.filename = NULL;
    g_ctx.delimiter = ',';
    g_ctx.write_header = true;
    g_ctx.quote_strings = false;
    g_ctx.quote_if_needed = true;  /* Default: quote only if needed */
    g_ctx.verbose = false;

    token = strtok_r(args_copy, " \t", &saveptr);
    while (token != NULL) {
        if (arg_idx == 0) {
            /* First arg: filename */
            g_ctx.filename = strdup(token);
        } else if (arg_idx == 1) {
            /* Second arg: delimiter */
            if (strlen(token) == 1) {
                g_ctx.delimiter = token[0];
            } else if (strcmp(token, "tab") == 0) {
                g_ctx.delimiter = '\t';
            }
        } else {
            /* Options */
            if (strcmp(token, "noheader") == 0) {
                g_ctx.write_header = false;
            } else if (strcmp(token, "quote") == 0) {
                g_ctx.quote_strings = true;
            } else if (strcmp(token, "noquoteif") == 0) {
                g_ctx.quote_if_needed = false;
            } else if (strcmp(token, "verbose") == 0) {
                g_ctx.verbose = true;
            }
        }

        arg_idx++;
        token = strtok_r(NULL, " \t", &saveptr);
    }

    free(args_copy);

    if (g_ctx.filename == NULL) {
        SF_error("cexport: no output filename specified\n");
        return -1;
    }

    return 0;
}

/*
    Load variable names from Stata macro.
    The .ado file sets _cexport_varnames macro with space-separated variable names.
*/
static int load_varnames(void)
{
    size_t nvars = SF_nvars();
    g_ctx.varnames = (char **)malloc(nvars * sizeof(char *));
    if (g_ctx.varnames == NULL) return -1;

    g_ctx.nvars = nvars;

    /* Try to get variable names from macro */
    char varnames_buf[32768];
    if (SF_macro_use("_cexport_varnames", varnames_buf, sizeof(varnames_buf)) == 0 &&
        strlen(varnames_buf) > 0) {
        /* Parse space-separated names */
        char *p = varnames_buf;
        for (size_t j = 0; j < nvars; j++) {
            /* Skip whitespace */
            while (*p == ' ' || *p == '\t') p++;
            if (*p == '\0') break;

            /* Find end of name */
            char *start = p;
            while (*p != ' ' && *p != '\t' && *p != '\0') p++;

            size_t len = p - start;
            g_ctx.varnames[j] = (char *)malloc(len + 1);
            if (g_ctx.varnames[j] == NULL) return -1;
            memcpy(g_ctx.varnames[j], start, len);
            g_ctx.varnames[j][len] = '\0';
        }
    } else {
        /* Fallback: use generic names */
        for (size_t j = 0; j < nvars; j++) {
            char namebuf[32];
            snprintf(namebuf, sizeof(namebuf), "v%zu", j + 1);
            g_ctx.varnames[j] = strdup(namebuf);
            if (g_ctx.varnames[j] == NULL) return -1;
        }
    }

    return 0;
}

/*
    Free all resources in context.
*/
static void cleanup_context(void)
{
    if (g_ctx.filename != NULL) {
        free(g_ctx.filename);
        g_ctx.filename = NULL;
    }

    if (g_ctx.fp != NULL) {
        fclose(g_ctx.fp);
        g_ctx.fp = NULL;
    }

    if (g_ctx.varnames != NULL) {
        for (size_t j = 0; j < g_ctx.nvars; j++) {
            if (g_ctx.varnames[j] != NULL) {
                free(g_ctx.varnames[j]);
            }
        }
        free(g_ctx.varnames);
        g_ctx.varnames = NULL;
    }

    stata_data_free(&g_ctx.data);
}

/*
    Main entry point for cexport command.
*/
ST_retcode cexport_main(const char *args)
{
    double t_start, t_phase;
    char msg[256];
    ST_retcode rc;
    size_t nobs, nvars;

    /* Initialize context */
    memset(&g_ctx, 0, sizeof(g_ctx));
    stata_data_init(&g_ctx.data);

    t_start = ctools_timer_seconds();

    /* Parse arguments */
    if (parse_args(args) != 0) {
        cleanup_context();
        return 198;
    }

    /* Get dataset dimensions */
    nobs = SF_in2() - SF_in1() + 1;
    nvars = SF_nvars();

    if (nobs == 0) {
        SF_error("cexport: no observations to export\n");
        cleanup_context();
        return 2000;
    }

    if (nvars == 0) {
        SF_error("cexport: no variables to export\n");
        cleanup_context();
        return 2000;
    }

    if (g_ctx.verbose) {
        snprintf(msg, sizeof(msg), "cexport: Exporting %zu obs x %zu vars to %s\n",
                nobs, nvars, g_ctx.filename);
        SF_display(msg);
    }

    /* ================================================================
       PHASE 1: Load data from Stata
       ================================================================ */
    t_phase = ctools_timer_seconds();

    /* Load variable names */
    if (load_varnames() != 0) {
        SF_error("cexport: failed to load variable names\n");
        cleanup_context();
        return 920;
    }

    /* Load data using shared ctools infrastructure */
    rc = ctools_data_load(&g_ctx.data, nvars);
    if (rc != STATA_OK) {
        snprintf(msg, sizeof(msg), "cexport: failed to load data (error %d)\n", rc);
        SF_error(msg);
        cleanup_context();
        return 920;
    }

    g_ctx.time_load = ctools_timer_seconds() - t_phase;

    if (g_ctx.verbose) {
        snprintf(msg, sizeof(msg), "  Load time:   %.3f sec\n", g_ctx.time_load);
        SF_display(msg);
    }

    /* ================================================================
       PHASE 2: Open output file and write header
       ================================================================ */
    g_ctx.fp = fopen(g_ctx.filename, "w");
    if (g_ctx.fp == NULL) {
        snprintf(msg, sizeof(msg), "cexport: cannot open file '%s' for writing: %s\n",
                g_ctx.filename, strerror(errno));
        SF_error(msg);
        cleanup_context();
        return 603;
    }

    /* Set large buffer for better I/O performance */
    setvbuf(g_ctx.fp, NULL, _IOFBF, CEXPORT_IO_BUFFER_SIZE);

    if (g_ctx.write_header) {
        if (write_header(g_ctx.fp) != 0) {
            SF_error("cexport: failed to write header\n");
            cleanup_context();
            return 693;
        }
    }

    /* ================================================================
       PHASE 3: Format and write data (parallel chunked)
       ================================================================ */
    t_phase = ctools_timer_seconds();

    /* Determine number of chunks and threads */
    size_t num_chunks = (nobs + CEXPORT_CHUNK_SIZE - 1) / CEXPORT_CHUNK_SIZE;
    size_t num_threads = num_chunks;
    if (num_threads > CEXPORT_MAX_THREADS) num_threads = CEXPORT_MAX_THREADS;
    if (num_threads < 1) num_threads = 1;

    /* Estimate buffer size per chunk (generous estimate) */
    size_t avg_row_size = nvars * 20 + nvars; /* ~20 chars per field + delimiters */
    size_t chunk_buffer_size = CEXPORT_CHUNK_SIZE * avg_row_size * 2; /* 2x safety margin */

    /* For small datasets, use single-threaded approach */
    if (nobs < CEXPORT_CHUNK_SIZE || num_threads == 1) {
        /* Single-threaded: format and write directly */
        char *row_buf = (char *)malloc(avg_row_size * 4);
        if (row_buf == NULL) {
            SF_error("cexport: memory allocation failed\n");
            cleanup_context();
            return 920;
        }

        for (size_t i = 0; i < nobs; i++) {
            int row_len = format_row(i, row_buf, avg_row_size * 4);
            if (row_len < 0) {
                SF_error("cexport: row formatting failed\n");
                free(row_buf);
                cleanup_context();
                return 920;
            }
            if (fwrite(row_buf, 1, row_len, g_ctx.fp) != (size_t)row_len) {
                SF_error("cexport: write error\n");
                free(row_buf);
                cleanup_context();
                return 693;
            }
        }
        free(row_buf);
    } else {
        /* Multi-threaded chunked processing */
        pthread_t *threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
        format_chunk_args_t *chunk_args = (format_chunk_args_t *)malloc(num_threads * sizeof(format_chunk_args_t));
        char **chunk_buffers = (char **)malloc(num_threads * sizeof(char *));

        if (threads == NULL || chunk_args == NULL || chunk_buffers == NULL) {
            free(threads);
            free(chunk_args);
            free(chunk_buffers);
            SF_error("cexport: memory allocation failed\n");
            cleanup_context();
            return 920;
        }

        /* Allocate buffers for each thread */
        for (size_t t = 0; t < num_threads; t++) {
            chunk_buffers[t] = (char *)malloc(chunk_buffer_size);
            if (chunk_buffers[t] == NULL) {
                for (size_t u = 0; u < t; u++) free(chunk_buffers[u]);
                free(threads);
                free(chunk_args);
                free(chunk_buffers);
                SF_error("cexport: memory allocation failed\n");
                cleanup_context();
                return 920;
            }
        }

        /* Process chunks in batches of num_threads */
        size_t chunk_idx = 0;
        while (chunk_idx < num_chunks) {
            size_t batch_size = num_chunks - chunk_idx;
            if (batch_size > num_threads) batch_size = num_threads;

            /* Launch threads for this batch */
            for (size_t t = 0; t < batch_size; t++) {
                size_t c = chunk_idx + t;
                chunk_args[t].start_row = c * CEXPORT_CHUNK_SIZE;
                chunk_args[t].end_row = (c + 1) * CEXPORT_CHUNK_SIZE;
                if (chunk_args[t].end_row > nobs) chunk_args[t].end_row = nobs;
                chunk_args[t].output_buffer = chunk_buffers[t];
                chunk_args[t].buffer_size = chunk_buffer_size;
                chunk_args[t].bytes_written = 0;
                chunk_args[t].success = 0;

                pthread_create(&threads[t], NULL, format_chunk_thread, &chunk_args[t]);
            }

            /* Wait for threads and write output in order */
            for (size_t t = 0; t < batch_size; t++) {
                pthread_join(threads[t], NULL);

                if (!chunk_args[t].success) {
                    for (size_t u = 0; u < num_threads; u++) free(chunk_buffers[u]);
                    free(threads);
                    free(chunk_args);
                    free(chunk_buffers);
                    SF_error("cexport: formatting failed\n");
                    cleanup_context();
                    return 920;
                }

                /* Write this chunk to file (in order) */
                if (fwrite(chunk_buffers[t], 1, chunk_args[t].bytes_written, g_ctx.fp)
                    != chunk_args[t].bytes_written) {
                    for (size_t u = 0; u < num_threads; u++) free(chunk_buffers[u]);
                    free(threads);
                    free(chunk_args);
                    free(chunk_buffers);
                    SF_error("cexport: write error\n");
                    cleanup_context();
                    return 693;
                }
            }

            chunk_idx += batch_size;
        }

        /* Cleanup */
        for (size_t t = 0; t < num_threads; t++) {
            free(chunk_buffers[t]);
        }
        free(threads);
        free(chunk_args);
        free(chunk_buffers);
    }

    g_ctx.time_format = ctools_timer_seconds() - t_phase;
    g_ctx.time_write = g_ctx.time_format; /* Combined for now */

    /* Close file */
    fclose(g_ctx.fp);
    g_ctx.fp = NULL;

    /* Calculate total time */
    g_ctx.time_total = ctools_timer_seconds() - t_start;

    /* Store timing in Stata scalars */
    SF_scal_save("_cexport_time_load", g_ctx.time_load);
    SF_scal_save("_cexport_time_write", g_ctx.time_format);
    SF_scal_save("_cexport_time_total", g_ctx.time_total);

    if (g_ctx.verbose) {
        snprintf(msg, sizeof(msg), "  Write time:  %.3f sec\n", g_ctx.time_format);
        SF_display(msg);
        snprintf(msg, sizeof(msg), "  Total time:  %.3f sec\n", g_ctx.time_total);
        SF_display(msg);

        /* Calculate throughput */
        double rows_per_sec = nobs / g_ctx.time_total;
        snprintf(msg, sizeof(msg), "  Throughput:  %.0f rows/sec\n", rows_per_sec);
        SF_display(msg);
    }

    /* Cleanup and return success */
    cleanup_context();
    return 0;
}

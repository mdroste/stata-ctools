/*
    cexport_impl.c
    High-Performance Delimited Text Export for Stata

    Exports Stata datasets to CSV/delimited text files with parallel processing.
    Designed for maximum throughput on large datasets.

    Architecture:
    1. PHASE 1: Load data from Stata into C memory (parallel per-variable)
    2. PHASE 2: Format data into text buffers (parallel per-chunk)
    3. PHASE 3: Write buffers to file (sequential or parallel I/O)

    Performance Optimizations:
    - Parallel data loading via ctools_data_load
    - Chunked parallel formatting: each thread formats a range of rows
    - Custom fast double-to-string conversion (in cexport_format.c)
    - Large I/O buffers to minimize syscalls
    - Memory-mapped output option for very large files
    - Pre-computed field widths for aligned output (optional)

    Thread Safety:
    - Each formatting thread works on disjoint row ranges
    - No shared mutable state during formatting phase
    - Parallel offset writes for I/O (pwrite backend)
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <errno.h>

#include "ctools_threads.h"

#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_timer.h"
#include "cexport_impl.h"
#include "cexport_io.h"
#include "cexport_context.h"
#include "cexport_format.h"
#include "cexport_parse.h"

/* ========================================================================
   Configuration Constants
   ======================================================================== */

#define CEXPORT_SAMPLE_ROWS      100   /* Number of rows to sample for sizing */
#define CEXPORT_ROW_BUFFER_ROWS 1000   /* Rows to buffer before fwrite */

/* ========================================================================
   Global Context
   ======================================================================== */

static cexport_context g_ctx;

/* ========================================================================
   Dynamic Buffer Sizing
   ======================================================================== */

/*
    Sample rows to compute actual average row size.
    Returns average bytes per row, or 0 on error.
*/
static size_t sample_row_sizes(size_t nobs)
{
    char sample_buf[4096];
    size_t sample_count = (nobs < CEXPORT_SAMPLE_ROWS) ? nobs : CEXPORT_SAMPLE_ROWS;
    size_t total_size = 0;
    size_t rows_sampled = 0;
    bool all_numeric = g_ctx.all_numeric;

    /* Sample evenly distributed rows */
    size_t step = (nobs > sample_count) ? (nobs / sample_count) : 1;

    for (size_t i = 0; i < nobs && rows_sampled < sample_count; i += step) {
        /* Check if this observation satisfies the if condition */
        ST_int stata_obs = (ST_int)(g_ctx.obs1 + i);
        if (!SF_ifobs(stata_obs)) {
            continue;
        }

        int row_len;
        if (all_numeric) {
            row_len = cexport_format_row_numeric(&g_ctx, i, sample_buf, sizeof(sample_buf));
        } else {
            row_len = cexport_format_row(&g_ctx, i, sample_buf, sizeof(sample_buf));
        }

        if (row_len > 0) {
            total_size += row_len;
            rows_sampled++;
        }
    }

    if (rows_sampled == 0) {
        return 0;
    }

    /* Return average with 20% safety margin */
    return (total_size / rows_sampled) * 12 / 10 + 64;
}

/*
    Compute adaptive chunk size based on data characteristics.
*/
static size_t compute_adaptive_chunk_size(size_t nobs, bool all_numeric)
{
    size_t base_chunk_size = CTOOLS_EXPORT_CHUNK_SIZE;

    if (all_numeric) {
        /* All-numeric: use 5x larger chunks */
        base_chunk_size = base_chunk_size * 5;
    }

    /* Cap at 10% of total rows or 100K rows */
    size_t max_chunk = nobs / 10;
    if (max_chunk < CTOOLS_EXPORT_CHUNK_SIZE) {
        max_chunk = CTOOLS_EXPORT_CHUNK_SIZE;
    }
    if (max_chunk > 100000) {
        max_chunk = 100000;
    }

    return (base_chunk_size < max_chunk) ? base_chunk_size : max_chunk;
}

/* ========================================================================
   Thread Functions
   ======================================================================== */

/*
    Thread function: Format a chunk of rows.
*/
static void *format_chunk_thread(void *arg)
{
    format_chunk_args_t *args = (format_chunk_args_t *)arg;
    size_t pos = 0;
    bool all_numeric = g_ctx.all_numeric;

    for (size_t i = args->start_row; i < args->end_row; i++) {
        ST_int stata_obs = (ST_int)(g_ctx.obs1 + i);
        if (!SF_ifobs(stata_obs)) {
            continue;
        }

        int row_len;
        if (all_numeric) {
            row_len = cexport_format_row_numeric(&g_ctx, i,
                args->output_buffer + pos, args->buffer_size - pos);
        } else {
            row_len = cexport_format_row(&g_ctx, i,
                args->output_buffer + pos, args->buffer_size - pos);
        }
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

/*
    Thread function: Write a chunk at a specific offset.
*/
static void *write_chunk_thread(void *arg)
{
    write_chunk_args_t *args = (write_chunk_args_t *)arg;

    if (args->len == 0) {
        args->success = 1;
        return NULL;
    }

    ssize_t written = cexport_io_write_at(args->file, args->buffer,
                                           args->len, args->offset);
    args->success = (written == (ssize_t)args->len);

    return args->success ? NULL : (void *)(intptr_t)-1;
}

/*
    Thread function: Format rows directly into mapped memory (zero-copy).
*/
static void *format_mmap_thread(void *arg)
{
    format_mmap_args_t *args = (format_mmap_args_t *)arg;
    size_t pos = 0;
    bool all_numeric = g_ctx.all_numeric;

    for (size_t i = args->start_row; i < args->end_row; i++) {
        ST_int stata_obs = (ST_int)(g_ctx.obs1 + i);
        if (!SF_ifobs(stata_obs)) {
            continue;
        }

        int row_len;
        if (all_numeric) {
            row_len = cexport_format_row_numeric(&g_ctx, i,
                args->dest + pos, args->dest_size - pos);
        } else {
            row_len = cexport_format_row(&g_ctx, i,
                args->dest + pos, args->dest_size - pos);
        }

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
   Single-Threaded Buffered Writer
   ======================================================================== */

static int write_rows_buffered(FILE *fp, size_t nobs, size_t avg_row_size, size_t *rows_written_out)
{
    bool all_numeric = g_ctx.all_numeric;
    size_t rows_written = 0;

    size_t buffer_size = avg_row_size * CEXPORT_ROW_BUFFER_ROWS;
    if (buffer_size < 65536) buffer_size = 65536;
    if (buffer_size > 4 * 1024 * 1024) buffer_size = 4 * 1024 * 1024;

    char *buffer = (char *)malloc(buffer_size);
    if (buffer == NULL) {
        return -1;
    }

    size_t buf_pos = 0;

    for (size_t i = 0; i < nobs; i++) {
        ST_int stata_obs = (ST_int)(g_ctx.obs1 + i);
        if (!SF_ifobs(stata_obs)) {
            continue;
        }

        if (buf_pos > buffer_size - avg_row_size * 2) {
            if (fwrite(buffer, 1, buf_pos, fp) != buf_pos) {
                free(buffer);
                return -1;
            }
            buf_pos = 0;
        }

        int row_len;
        if (all_numeric) {
            row_len = cexport_format_row_numeric(&g_ctx, i, buffer + buf_pos, buffer_size - buf_pos);
        } else {
            row_len = cexport_format_row(&g_ctx, i, buffer + buf_pos, buffer_size - buf_pos);
        }

        if (row_len < 0) {
            if (buf_pos > 0) {
                if (fwrite(buffer, 1, buf_pos, fp) != buf_pos) {
                    free(buffer);
                    return -1;
                }
                buf_pos = 0;
            }

            if (all_numeric) {
                row_len = cexport_format_row_numeric(&g_ctx, i, buffer + buf_pos, buffer_size - buf_pos);
            } else {
                row_len = cexport_format_row(&g_ctx, i, buffer + buf_pos, buffer_size - buf_pos);
            }

            if (row_len < 0) {
                free(buffer);
                return -1;
            }
        }

        buf_pos += row_len;
        rows_written++;
    }

    if (buf_pos > 0) {
        if (fwrite(buffer, 1, buf_pos, fp) != buf_pos) {
            free(buffer);
            return -1;
        }
    }

    free(buffer);
    *rows_written_out = rows_written;
    return 0;
}

/* ========================================================================
   Export Pipeline Implementations
   ======================================================================== */

/*
    Single-threaded export path using buffered FILE* I/O.
*/
static ST_retcode export_single_threaded(size_t nobs, size_t avg_row_size,
                                          char *header_buf, size_t header_len)
{
    char msg[512];

    g_ctx.fp = fopen(g_ctx.filename, "wb");
    if (g_ctx.fp == NULL) {
        snprintf(msg, sizeof(msg), "cexport: cannot open file '%s' for writing: %s\n",
                g_ctx.filename, strerror(errno));
        SF_error(msg);
        free(header_buf);
        return 603;
    }
    setvbuf(g_ctx.fp, NULL, _IOFBF, CTOOLS_IO_BUFFER_SIZE);

    if (header_buf != NULL && header_len > 0) {
        if (fwrite(header_buf, 1, header_len, g_ctx.fp) != header_len) {
            SF_error("cexport: failed to write header\n");
            free(header_buf);
            fclose(g_ctx.fp);
            g_ctx.fp = NULL;
            return 693;
        }
    }
    free(header_buf);

    size_t rows_written = 0;
    if (write_rows_buffered(g_ctx.fp, nobs, avg_row_size, &rows_written) != 0) {
        SF_error("cexport: write error during buffered write\n");
        fclose(g_ctx.fp);
        g_ctx.fp = NULL;
        return 693;
    }

    fclose(g_ctx.fp);
    g_ctx.fp = NULL;
    return 0;
}

/*
    Memory-mapped export path with zero-copy formatting.
*/
static ST_retcode export_mmap(size_t nobs, size_t avg_row_size, size_t chunk_size,
                               size_t num_chunks, char *header_buf, size_t header_len,
                               ctools_persistent_pool *pool)
{
    char msg[512];
    size_t estimated_total = header_len + (nobs * avg_row_size * 15 / 10);

    cexport_io_file outfile;
    cexport_io_init(&outfile);

    if (cexport_io_open(&outfile, g_ctx.filename, CEXPORT_IO_MMAP, g_ctx.io_flags) != 0) {
        snprintf(msg, sizeof(msg), "cexport: cannot open file '%s': %s\n",
                g_ctx.filename, outfile.error_message);
        SF_error(msg);
        free(header_buf);
        return 603;
    }

    if (cexport_io_presize(&outfile, estimated_total) != 0) {
        snprintf(msg, sizeof(msg), "cexport: failed to pre-size file: %s\n",
                outfile.error_message);
        SF_error(msg);
        cexport_io_close(&outfile, 0);
        free(header_buf);
        return 693;
    }

    char *mapped_base = cexport_io_get_mapped_ptr(&outfile, 0);
    if (mapped_base == NULL) {
        SF_error("cexport: failed to get mapped memory pointer\n");
        cexport_io_close(&outfile, 0);
        free(header_buf);
        return 693;
    }

    if (header_buf != NULL && header_len > 0) {
        memcpy(mapped_base, header_buf, header_len);
    }
    free(header_buf);

    size_t *chunk_offsets = (size_t *)ctools_safe_malloc2(num_chunks + 1, sizeof(size_t));
    if (chunk_offsets == NULL) {
        SF_error("cexport: memory allocation failed\n");
        cexport_io_close(&outfile, 0);
        return 920;
    }

    chunk_offsets[0] = header_len;
    for (size_t c = 0; c < num_chunks; c++) {
        size_t rows_in_chunk = chunk_size;
        if ((c + 1) * chunk_size > nobs) {
            rows_in_chunk = nobs - c * chunk_size;
        }
        chunk_offsets[c + 1] = chunk_offsets[c] + (rows_in_chunk * avg_row_size * 15 / 10);
    }

    format_mmap_args_t *mmap_args = (format_mmap_args_t *)ctools_safe_malloc2(num_chunks, sizeof(format_mmap_args_t));
    if (mmap_args == NULL) {
        SF_error("cexport: memory allocation failed\n");
        free(chunk_offsets);
        cexport_io_close(&outfile, 0);
        return 920;
    }

    for (size_t c = 0; c < num_chunks; c++) {
        mmap_args[c].dest = mapped_base + chunk_offsets[c];
        mmap_args[c].dest_size = chunk_offsets[c + 1] - chunk_offsets[c];
        mmap_args[c].start_row = c * chunk_size;
        mmap_args[c].end_row = (c + 1) * chunk_size;
        if (mmap_args[c].end_row > nobs) mmap_args[c].end_row = nobs;
        mmap_args[c].bytes_written = 0;
        mmap_args[c].success = 0;
    }

    for (size_t c = 0; c < num_chunks; c++) {
        ctools_persistent_pool_submit(pool, format_mmap_thread, &mmap_args[c]);
    }
    ctools_persistent_pool_wait(pool);

    int format_failed = 0;
    for (size_t c = 0; c < num_chunks; c++) {
        if (!mmap_args[c].success) {
            format_failed = 1;
            break;
        }
    }

    if (format_failed) {
        SF_error("cexport: formatting to mmap failed\n");
        free(mmap_args);
        free(chunk_offsets);
        cexport_io_close(&outfile, 0);
        return 920;
    }

    /* Compact chunks together */
    size_t write_pos = header_len;
    for (size_t c = 0; c < num_chunks; c++) {
        if (mmap_args[c].bytes_written > 0 && chunk_offsets[c] != write_pos) {
            memmove(mapped_base + write_pos,
                    mapped_base + chunk_offsets[c],
                    mmap_args[c].bytes_written);
        }
        write_pos += mmap_args[c].bytes_written;
    }
    size_t actual_total = write_pos;

    if (cexport_io_close(&outfile, actual_total) != 0) {
        SF_error("cexport: failed to close file\n");
        free(mmap_args);
        free(chunk_offsets);
        return 693;
    }

    free(mmap_args);
    free(chunk_offsets);
    return 0;
}

/*
    Parallel pwrite export path.
*/
static ST_retcode export_pwrite(size_t nobs, size_t avg_row_size, size_t chunk_size,
                                 size_t num_chunks, size_t chunk_buffer_size,
                                 char *header_buf, size_t header_len,
                                 ctools_persistent_pool *pool)
{
    (void)avg_row_size;  /* Used for buffer sizing, computed externally */
    char msg[512];

    format_chunk_args_t *chunk_args = (format_chunk_args_t *)ctools_safe_malloc2(num_chunks, sizeof(format_chunk_args_t));
    char **chunk_buffers = (char **)ctools_safe_malloc2(num_chunks, sizeof(char *));

    if (chunk_args == NULL || chunk_buffers == NULL) {
        free(chunk_args);
        free(chunk_buffers);
        free(header_buf);
        SF_error("cexport: memory allocation failed\n");
        return 920;
    }

    for (size_t c = 0; c < num_chunks; c++) {
        chunk_buffers[c] = NULL;
    }

    for (size_t c = 0; c < num_chunks; c++) {
        chunk_buffers[c] = (char *)malloc(chunk_buffer_size);
        if (chunk_buffers[c] == NULL) {
            for (size_t u = 0; u < num_chunks; u++) {
                if (chunk_buffers[u]) free(chunk_buffers[u]);
            }
            free(chunk_args);
            free(chunk_buffers);
            free(header_buf);
            SF_error("cexport: memory allocation failed for chunk buffers\n");
            return 920;
        }
    }

    for (size_t c = 0; c < num_chunks; c++) {
        chunk_args[c].start_row = c * chunk_size;
        chunk_args[c].end_row = (c + 1) * chunk_size;
        if (chunk_args[c].end_row > nobs) chunk_args[c].end_row = nobs;
        chunk_args[c].output_buffer = chunk_buffers[c];
        chunk_args[c].buffer_size = chunk_buffer_size;
        chunk_args[c].bytes_written = 0;
        chunk_args[c].success = 0;
    }

    for (size_t c = 0; c < num_chunks; c++) {
        ctools_persistent_pool_submit(pool, format_chunk_thread, &chunk_args[c]);
    }
    ctools_persistent_pool_wait(pool);

    for (size_t c = 0; c < num_chunks; c++) {
        if (!chunk_args[c].success) {
            for (size_t u = 0; u < num_chunks; u++) {
                if (chunk_buffers[u]) free(chunk_buffers[u]);
            }
            free(chunk_args);
            free(chunk_buffers);
            free(header_buf);
            SF_error("cexport: formatting failed\n");
            return 920;
        }
    }

    /* Compute file offsets */
    size_t *offsets = (size_t *)ctools_safe_malloc2(num_chunks + 1, sizeof(size_t));
    if (offsets == NULL) {
        for (size_t u = 0; u < num_chunks; u++) {
            if (chunk_buffers[u]) free(chunk_buffers[u]);
        }
        free(chunk_args);
        free(chunk_buffers);
        free(header_buf);
        SF_error("cexport: memory allocation failed for offsets\n");
        return 920;
    }

    offsets[0] = header_len;
    for (size_t c = 0; c < num_chunks; c++) {
        offsets[c + 1] = offsets[c] + chunk_args[c].bytes_written;
    }
    size_t total_file_size = offsets[num_chunks];

    /* Open and pre-size file */
    cexport_io_file outfile;
    cexport_io_init(&outfile);

    if (cexport_io_open(&outfile, g_ctx.filename, g_ctx.io_backend, g_ctx.io_flags) != 0) {
        snprintf(msg, sizeof(msg), "cexport: cannot open file '%s': %s\n",
                g_ctx.filename, outfile.error_message);
        SF_error(msg);
        free(offsets);
        for (size_t u = 0; u < num_chunks; u++) {
            if (chunk_buffers[u]) free(chunk_buffers[u]);
        }
        free(chunk_args);
        free(chunk_buffers);
        free(header_buf);
        return 603;
    }

    if (cexport_io_presize(&outfile, total_file_size) != 0) {
        snprintf(msg, sizeof(msg), "cexport: failed to pre-size file: %s\n",
                outfile.error_message);
        SF_error(msg);
        cexport_io_close(&outfile, 0);
        free(offsets);
        for (size_t u = 0; u < num_chunks; u++) {
            if (chunk_buffers[u]) free(chunk_buffers[u]);
        }
        free(chunk_args);
        free(chunk_buffers);
        free(header_buf);
        return 693;
    }

    /* Write header */
    if (header_buf != NULL && header_len > 0) {
        ssize_t written = cexport_io_write_at(&outfile, header_buf, header_len, 0);
        if (written != (ssize_t)header_len) {
            snprintf(msg, sizeof(msg), "cexport: failed to write header: %s\n",
                    outfile.error_message);
            SF_error(msg);
            cexport_io_close(&outfile, 0);
            free(offsets);
            for (size_t u = 0; u < num_chunks; u++) {
                if (chunk_buffers[u]) free(chunk_buffers[u]);
            }
            free(chunk_args);
            free(chunk_buffers);
            free(header_buf);
            return 693;
        }
    }
    free(header_buf);

    /* Parallel offset writes */
    write_chunk_args_t *write_args = (write_chunk_args_t *)ctools_safe_malloc2(num_chunks, sizeof(write_chunk_args_t));
    if (write_args == NULL) {
        SF_error("cexport: memory allocation failed for write args\n");
        cexport_io_close(&outfile, 0);
        free(offsets);
        for (size_t u = 0; u < num_chunks; u++) {
            if (chunk_buffers[u]) free(chunk_buffers[u]);
        }
        free(chunk_args);
        free(chunk_buffers);
        return 920;
    }

    for (size_t c = 0; c < num_chunks; c++) {
        write_args[c].file = &outfile;
        write_args[c].buffer = chunk_buffers[c];
        write_args[c].len = chunk_args[c].bytes_written;
        write_args[c].offset = offsets[c];
        write_args[c].success = 0;
    }

    for (size_t c = 0; c < num_chunks; c++) {
        ctools_persistent_pool_submit(pool, write_chunk_thread, &write_args[c]);
    }
    ctools_persistent_pool_wait(pool);

    int write_failed = 0;
    for (size_t c = 0; c < num_chunks; c++) {
        if (!write_args[c].success) {
            write_failed = 1;
            break;
        }
    }

    if (cexport_io_close(&outfile, total_file_size) != 0 || write_failed) {
        SF_error("cexport: write error during parallel I/O\n");
        free(write_args);
        free(offsets);
        for (size_t u = 0; u < num_chunks; u++) {
            if (chunk_buffers[u]) free(chunk_buffers[u]);
        }
        free(chunk_args);
        free(chunk_buffers);
        return 693;
    }

    free(write_args);
    free(offsets);
    for (size_t c = 0; c < num_chunks; c++) {
        if (chunk_buffers[c]) free(chunk_buffers[c]);
    }
    free(chunk_args);
    free(chunk_buffers);
    return 0;
}

/* ========================================================================
   Main Export Function
   ======================================================================== */

ST_retcode cexport_main(const char *args)
{
    double t_start, t_phase;
    char msg[512];
    ST_retcode rc;
    size_t nobs, nvars;

    /* Initialize context */
    cexport_context_init(&g_ctx);
    t_start = ctools_timer_seconds();

    /* Parse arguments */
    if (cexport_parse_args(&g_ctx, args) != 0) {
        cexport_context_cleanup(&g_ctx);
        return 198;
    }

    /* Get dataset dimensions */
    g_ctx.obs1 = (size_t)SF_in1();
    nobs = SF_in2() - SF_in1() + 1;
    g_ctx.nobs_loaded = nobs;
    nvars = SF_nvars();

    if (nobs == 0) {
        SF_error("cexport: no observations to export\n");
        cexport_context_cleanup(&g_ctx);
        return 2000;
    }

    if (nvars == 0) {
        SF_error("cexport: no variables to export\n");
        cexport_context_cleanup(&g_ctx);
        return 2000;
    }

    /* ================================================================
       PHASE 1: Load data from Stata
       ================================================================ */
    t_phase = ctools_timer_seconds();

    if (cexport_load_varnames(&g_ctx) != 0) {
        SF_error("cexport: failed to load variable names\n");
        cexport_context_cleanup(&g_ctx);
        return 920;
    }

    if (cexport_load_vartypes(&g_ctx) != 0) {
        SF_error("cexport: failed to load variable types\n");
        cexport_context_cleanup(&g_ctx);
        return 920;
    }

    rc = ctools_data_load(&g_ctx.data, nvars);
    if (rc != STATA_OK) {
        snprintf(msg, sizeof(msg), "cexport: failed to load data (error %d)\n", rc);
        SF_error(msg);
        cexport_context_cleanup(&g_ctx);
        return 920;
    }

    g_ctx.time_load = ctools_timer_seconds() - t_phase;

    /* Detect if all variables are numeric */
    g_ctx.all_numeric = true;
    for (size_t j = 0; j < nvars; j++) {
        if (g_ctx.data.vars[j].type != STATA_TYPE_DOUBLE) {
            g_ctx.all_numeric = false;
            break;
        }
    }

    /* ================================================================
       PHASE 2: Format and write data
       ================================================================ */
    t_phase = ctools_timer_seconds();

    size_t chunk_size = compute_adaptive_chunk_size(nobs, g_ctx.all_numeric);
    g_ctx.adaptive_chunk_size = chunk_size;

    size_t num_chunks = (nobs + chunk_size - 1) / chunk_size;
    size_t num_threads = num_chunks;
    size_t max_threads = (size_t)ctools_get_max_threads();
    if (num_threads > max_threads) num_threads = max_threads;
    if (num_threads < 1) num_threads = 1;

    size_t avg_row_size = sample_row_sizes(nobs);
    if (avg_row_size == 0) {
        if (nvars > SIZE_MAX / 22) {
            SF_error("cexport: variable count too large\n");
            cexport_context_cleanup(&g_ctx);
            return 920;
        }
        avg_row_size = nvars * 22;
    }
    g_ctx.actual_avg_row_size = avg_row_size;

    size_t chunk_buffer_size;
    if (avg_row_size > SIZE_MAX / chunk_size) {
        SF_error("cexport: buffer size overflow\n");
        cexport_context_cleanup(&g_ctx);
        return 920;
    }
    chunk_buffer_size = chunk_size * avg_row_size * 12 / 10 + 4096;

    /* Format header */
    char *header_buf = NULL;
    size_t header_len = 0;
    if (g_ctx.write_header) {
        header_buf = (char *)malloc(65536);
        if (header_buf == NULL) {
            SF_error("cexport: memory allocation failed for header\n");
            cexport_context_cleanup(&g_ctx);
            return 920;
        }
        int hlen = cexport_format_header(&g_ctx, header_buf, 65536);
        if (hlen < 0) {
            SF_error("cexport: failed to format header\n");
            free(header_buf);
            cexport_context_cleanup(&g_ctx);
            return 920;
        }
        header_len = (size_t)hlen;
    }

    /* Choose export path */
    ST_retcode export_rc;

    if (!g_ctx.use_parallel_io || nobs < CTOOLS_EXPORT_CHUNK_SIZE || num_threads == 1) {
        export_rc = export_single_threaded(nobs, avg_row_size, header_buf, header_len);
        g_ctx.time_format = ctools_timer_seconds() - t_phase;
        g_ctx.time_write = g_ctx.time_format;
    } else {
        ctools_persistent_pool *pool = ctools_get_global_pool();
        if (pool == NULL) {
            if (g_ctx.verbose) {
                SF_display("  Thread pool unavailable, using single-threaded fallback\n");
            }
            export_rc = export_single_threaded(nobs, avg_row_size, header_buf, header_len);
            g_ctx.time_format = ctools_timer_seconds() - t_phase;
            g_ctx.time_write = g_ctx.time_format;
        } else if (g_ctx.io_backend == CEXPORT_IO_MMAP) {
            double t_format_start = ctools_timer_seconds();
            export_rc = export_mmap(nobs, avg_row_size, chunk_size, num_chunks,
                                     header_buf, header_len, pool);
            g_ctx.time_format = ctools_timer_seconds() - t_format_start;
            g_ctx.time_write = ctools_timer_seconds() - t_phase - g_ctx.time_format;
        } else {
            double t_format_start = ctools_timer_seconds();
            export_rc = export_pwrite(nobs, avg_row_size, chunk_size, num_chunks,
                                       chunk_buffer_size, header_buf, header_len, pool);
            g_ctx.time_format = ctools_timer_seconds() - t_format_start;
            g_ctx.time_write = ctools_timer_seconds() - t_phase - g_ctx.time_format;
        }
    }

    if (export_rc != 0) {
        cexport_context_cleanup(&g_ctx);
        return export_rc;
    }

    /* Calculate total time */
    g_ctx.time_total = ctools_timer_seconds() - t_start;

    /* Store timing in Stata scalars */
    SF_scal_save("_cexport_time_load", g_ctx.time_load);
    SF_scal_save("_cexport_time_format", g_ctx.time_format);
    SF_scal_save("_cexport_time_write", g_ctx.time_write);
    SF_scal_save("_cexport_time_total", g_ctx.time_total);
    CTOOLS_SAVE_THREAD_INFO("_cexport");

    cexport_context_cleanup(&g_ctx);
    return 0;
}

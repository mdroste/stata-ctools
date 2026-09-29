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
#include "ctools_runtime.h"
#include "cexport_impl.h"
#include "cexport_io.h"
#include "cexport_context.h"
#include "cexport_format.h"
#include "cexport_parse.h"
#include "cexport_xlsx.h"

/* ========================================================================
   Configuration Constants
   ======================================================================== */

#define CEXPORT_SAMPLE_ROWS      500   /* Number of rows to sample for sizing */
#define CEXPORT_ROW_BUFFER_ROWS 1000   /* Rows to buffer before fwrite */

/* export_mmap: chunk windows were too small; redo with export_pwrite */
#define CEXPORT_RETRY_PWRITE (-1)

/* ========================================================================
   Global Context
   ======================================================================== */

static cexport_context g_ctx;

/* ========================================================================
   Row Dispatch
   ======================================================================== */

/* Format the i-th selected observation from Stata or from the loaded copy. */
static inline int format_row_at(size_t i, char *buf, size_t size)
{
    if (g_ctx.stream)
        return cexport_format_row_stream(&g_ctx, (ST_int)g_ctx.filtered.obs_map[i], buf, size);
    if (g_ctx.all_numeric) return cexport_format_row_numeric(&g_ctx, i, buf, size);
    return cexport_format_row(&g_ctx, i, buf, size);
}

/* Room a formatter may reserve for one row: streaming reads each string at
 * full fixed-string width, and quoting can double it. */
static size_t row_reserve_bytes(void)
{
    size_t bytes = 64 + (size_t)g_ctx.line_ending_len;
    for (size_t j = 0; j < g_ctx.nvars; j++) {
        size_t add = 64;
        if (g_ctx.stream) {
            if (g_ctx.is_string[j]) add = 2 * 2046 + 4;
        } else if (j < g_ctx.filtered.data.nvars && g_ctx.filtered.data.vars[j].type == STATA_TYPE_STRING) {
            add = 2 * g_ctx.filtered.data.vars[j].str_maxlen + 4;
        }
        if (add > SIZE_MAX - bytes) return SIZE_MAX;
        bytes += add;
    }
    return bytes;
}

/* ========================================================================
   Dynamic Buffer Sizing
   ======================================================================== */

/*
    Sample rows to compute actual average row size.
    Returns average bytes per row, or 0 on error.
*/
static size_t sample_row_sizes(size_t nobs, bool longtext)
{
    size_t sample_size = row_reserve_bytes();
    char *sample_buf = sample_size == SIZE_MAX ? NULL : malloc(sample_size);
    if (!sample_buf) return 0;
    size_t sample_count = (nobs < CEXPORT_SAMPLE_ROWS) ? nobs : CEXPORT_SAMPLE_ROWS;
    size_t total_size = 0;
    size_t max_row_size = 0;
    size_t rows_sampled = 0;

    /* Sample evenly distributed rows */
    size_t step = (nobs > sample_count) ? (nobs / sample_count) : 1;

    for (size_t i = 0; i < nobs && rows_sampled < sample_count; i += step) {
        int row_len = format_row_at(i, sample_buf, sample_size);

        if (row_len > 0) {
            total_size += row_len;
            if ((size_t)row_len > max_row_size) max_row_size = (size_t)row_len;
            rows_sampled++;
        }
    }

    free(sample_buf);
    if (rows_sampled == 0) {
        return 0;
    }

    /* Use average + 10% margin, but ensure at least max_row_size + padding */
    size_t avg_with_margin = (total_size / rows_sampled) * 11 / 10 + 64;
    size_t min_safe = max_row_size + 64;
    if (longtext && sample_size > min_safe) min_safe = sample_size;
    return (avg_with_margin > min_safe) ? avg_with_margin : min_safe;
}

/* Upper bound on the formatted bytes of one chunk buffer */
#define CEXPORT_MAX_CHUNK_BYTES ((size_t)16 << 20)

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
    args->read_error = 0;

    for (size_t i = args->start_row; i < args->end_row; i++) {
        int row_len;
        while ((row_len = format_row_at(i, args->output_buffer + pos, args->buffer_size - pos)) == -1) {
            /* Rows can exceed the sampled size estimate; grow the chunk
             * buffer rather than failing the export. */
            size_t grown = args->buffer_size * 2;
            char *buf = args->buffer_slot && grown > args->buffer_size && grown <= ((size_t)1 << 36)
                ? (char *)realloc(args->output_buffer, grown) : NULL;
            if (buf == NULL) {
                args->success = 0;
                args->bytes_written = pos;
                return NULL;
            }
            args->output_buffer = *args->buffer_slot = buf;
            args->buffer_size = *args->size_slot = grown;
        }
        if (row_len < 0) {
            args->read_error = 1;
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
    Thread function: Format rows directly into mapped memory (zero-copy).
*/
static void *format_mmap_thread(void *arg)
{
    format_mmap_args_t *args = (format_mmap_args_t *)arg;
    size_t pos = 0;

    for (size_t i = args->start_row; i < args->end_row; i++) {
        int row_len = format_row_at(i, args->dest + pos, args->dest_size - pos);

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
    size_t rows_written = 0;

    /* Keep room for a worst-case row after every flush. */
    size_t reserve = row_reserve_bytes();
    size_t buffer_size = avg_row_size * CEXPORT_ROW_BUFFER_ROWS;
    if (buffer_size < 65536) buffer_size = 65536;
    if (buffer_size > 4 * 1024 * 1024) buffer_size = 4 * 1024 * 1024;
    if (buffer_size < 2 * reserve) buffer_size = 2 * reserve;
    if (buffer_size < avg_row_size * 2 + 64) buffer_size = avg_row_size * 2 + 64;

    char *buffer = (char *)malloc(buffer_size);
    if (buffer == NULL) {
        return -1;
    }

    size_t buf_pos = 0;

    for (size_t i = 0; i < nobs; i++) {
        if (buf_pos > buffer_size - reserve) {
            if (fwrite(buffer, 1, buf_pos, fp) != buf_pos) {
                free(buffer);
                return -1;
            }
            buf_pos = 0;
        }

        int row_len = format_row_at(i, buffer + buf_pos, buffer_size - buf_pos);

        if (row_len == -1) {
            if (buf_pos > 0) {
                if (fwrite(buffer, 1, buf_pos, fp) != buf_pos) {
                    free(buffer);
                    return -1;
                }
                buf_pos = 0;
            }

            row_len = format_row_at(i, buffer + buf_pos, buffer_size - buf_pos);
        }
        if (row_len < 0) {
            free(buffer);
            return -1;
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

    int close_rc = fclose(g_ctx.fp);
    g_ctx.fp = NULL;
    return close_rc ? 693 : 0;
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

    /* Arena-allocated metadata */
    ctools_arena *arena = &g_ctx.chunk_arena;
    ctools_arena_init(arena, 64 * 1024);

    size_t *chunk_offsets = (size_t *)ctools_arena_alloc(arena, (num_chunks + 1) * sizeof(size_t));
    format_mmap_args_t *mmap_args = (format_mmap_args_t *)ctools_arena_alloc(arena, num_chunks * sizeof(format_mmap_args_t));

    if (chunk_offsets == NULL || mmap_args == NULL) {
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
        /* A chunk outgrew its fixed mapped window; the caller rewrites the
         * (unpublished) temporary file through the growable pwrite path. */
        cexport_io_close(&outfile, 0);
        return CEXPORT_RETRY_PWRITE;
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
        return 693;
    }

    /* Arena freed by cexport_context_cleanup() */
    return 0;
}

/*
    Parallel pwrite export path.

    Uses a sliding-window approach: only num_threads buffers are allocated
    (not num_chunks), and chunks are processed in waves. This dramatically
    reduces memory for large datasets.
    Example: 10M rows at 10K/chunk = 1000 chunks * 2.4MB = 2.4GB
             -> 8 threads * 2.4MB = 19MB
*/
static ST_retcode export_pwrite(size_t nobs, size_t avg_row_size, size_t chunk_size,
                                 size_t num_chunks, size_t chunk_buffer_size,
                                 char *header_buf, size_t header_len,
                                 ctools_persistent_pool *pool)
{
    (void)avg_row_size;  /* Used for buffer sizing, computed externally */
    char msg[512];
    ST_retcode ret = 0;

    /* Determine wave size = number of concurrent buffers */
    size_t max_threads = (size_t)ctools_get_max_threads();
    if (max_threads < 1) max_threads = 1;
    size_t wave_size = max_threads;
    if (wave_size > num_chunks) wave_size = num_chunks;

    /* Initialize arena for small metadata allocations */
    ctools_arena *arena = &g_ctx.chunk_arena;
    ctools_arena_init(arena, 64 * 1024);  /* 64KB blocks */

    /* Arena-allocated metadata (sized to total num_chunks for offset tracking) */
    format_chunk_args_t *chunk_args = (format_chunk_args_t *)ctools_arena_alloc(arena, num_chunks * sizeof(format_chunk_args_t));
    cexport_io_batch_entry *batch_entries = (cexport_io_batch_entry *)ctools_arena_alloc(arena, wave_size * sizeof(cexport_io_batch_entry));

    /* Allocate 2x wave_size buffer pointers for double-buffered pipeline */
    size_t total_buffers = wave_size * 2;
    if (total_buffers > num_chunks + wave_size) total_buffers = num_chunks + wave_size;
    char **wave_buffers = (char **)ctools_arena_alloc(arena, total_buffers * sizeof(char *));
    size_t *wave_sizes = (size_t *)ctools_arena_alloc(arena, total_buffers * sizeof(size_t));

    if (chunk_args == NULL || batch_entries == NULL || wave_buffers == NULL || wave_sizes == NULL) {
        free(header_buf);
        SF_error("cexport: memory allocation failed\n");
        return 920;
    }

    for (size_t b = 0; b < total_buffers; b++) {
        wave_buffers[b] = NULL;
    }

    /* Allocate double-buffered wave buffers */
    for (size_t b = 0; b < total_buffers; b++) {
        wave_buffers[b] = (char *)malloc(chunk_buffer_size);
        wave_sizes[b] = chunk_buffer_size;
        if (wave_buffers[b] == NULL) {
            for (size_t u = 0; u < total_buffers; u++) {
                if (wave_buffers[u]) free(wave_buffers[u]);
            }
            free(header_buf);
            SF_error("cexport: memory allocation failed for chunk buffers\n");
            return 920;
        }
    }

    /* Initialize all chunk_args (row ranges only; buffers assigned per wave) */
    for (size_t c = 0; c < num_chunks; c++) {
        chunk_args[c].start_row = c * chunk_size;
        chunk_args[c].end_row = (c + 1) * chunk_size;
        if (chunk_args[c].end_row > nobs) chunk_args[c].end_row = nobs;
        chunk_args[c].buffer_size = chunk_buffer_size;
        chunk_args[c].buffer_slot = NULL;
        chunk_args[c].size_slot = NULL;
        chunk_args[c].bytes_written = 0;
        chunk_args[c].success = 0;
        chunk_args[c].read_error = 0;
    }

    /* Estimate total file size for pre-sizing */
    size_t estimated_total = header_len + (nobs * avg_row_size * 12 / 10);

    /* Open and pre-size file */
    cexport_io_file outfile;
    cexport_io_init(&outfile);

    if (cexport_io_open(&outfile, g_ctx.filename, g_ctx.io_backend, g_ctx.io_flags) != 0) {
        snprintf(msg, sizeof(msg), "cexport: cannot open file '%s': %s\n",
                g_ctx.filename, outfile.error_message);
        SF_error(msg);
        for (size_t u = 0; u < total_buffers; u++) {
            if (wave_buffers[u]) free(wave_buffers[u]);
        }
        free(header_buf);
        return 603;
    }

    if (cexport_io_presize(&outfile, estimated_total) != 0) {
        snprintf(msg, sizeof(msg), "cexport: failed to pre-size file: %s\n",
                outfile.error_message);
        SF_error(msg);
        cexport_io_close(&outfile, 0);
        for (size_t u = 0; u < total_buffers; u++) {
            if (wave_buffers[u]) free(wave_buffers[u]);
        }
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
            for (size_t u = 0; u < total_buffers; u++) {
                if (wave_buffers[u]) free(wave_buffers[u]);
            }
            free(header_buf);
            return 693;
        }
    }
    free(header_buf);

    /* Pipelined format/write: format wave N+1 while writing wave N.
     * Uses alternating buffer sets (A and B) so writes and formats
     * never touch the same buffers simultaneously. */
    size_t file_offset = header_len;
    int cur_buf_set = 0;

    /* Format first wave */
    size_t wave_start = 0;
    size_t wave_end_first = (wave_start + wave_size < num_chunks) ? wave_start + wave_size : num_chunks;
    size_t wave_count = wave_end_first - wave_start;

    for (size_t w = 0; w < wave_count; w++) {
        size_t slot = cur_buf_set * wave_size + w;
        chunk_args[wave_start + w].output_buffer = wave_buffers[slot];
        chunk_args[wave_start + w].buffer_size = wave_sizes[slot];
        chunk_args[wave_start + w].buffer_slot = &wave_buffers[slot];
        chunk_args[wave_start + w].size_slot = &wave_sizes[slot];
        chunk_args[wave_start + w].bytes_written = 0;
        chunk_args[wave_start + w].success = 0;
    }
    for (size_t w = 0; w < wave_count; w++) {
        ctools_persistent_pool_submit(pool, format_chunk_thread, &chunk_args[wave_start + w]);
    }
    ctools_persistent_pool_wait(pool);

    for (size_t w = 0; w < wave_count; w++) {
        if (!chunk_args[wave_start + w].success) {
            for (size_t u = 0; u < total_buffers; u++) {
                if (wave_buffers[u]) free(wave_buffers[u]);
            }
            cexport_io_close(&outfile, 0);
            SF_error(chunk_args[wave_start + w].read_error ? "cexport: failed to read Stata data\n"
                                                           : "cexport: formatting failed\n");
            return 920;
        }
    }

    /* Main pipeline: overlap write of current wave with format of next wave.
     * Writes use batched I/O from the main thread (cexport_io_write_batch)
     * while format tasks run concurrently on the thread pool.
     * On Windows: async WriteFile with pre-allocated events + WaitForMultipleObjects.
     * On POSIX: pwritev for contiguous ranges (single syscall). */
    while (wave_start < num_chunks) {
        /* Build batch entries for current wave's writes */
        size_t wave_bytes = 0;
        for (size_t w = 0; w < wave_count; w++) {
            batch_entries[w].buf = chunk_args[wave_start + w].output_buffer;
            batch_entries[w].len = chunk_args[wave_start + w].bytes_written;
            wave_bytes += chunk_args[wave_start + w].bytes_written;
            batch_entries[w].offset = file_offset;
            file_offset += chunk_args[wave_start + w].bytes_written;
        }

        /* Determine next wave */
        size_t next_start = wave_start + wave_size;
        size_t next_end = (next_start + wave_size < num_chunks) ? next_start + wave_size : num_chunks;
        size_t next_count = (next_start < num_chunks) ? next_end - next_start : 0;
        int next_buf_set = 1 - cur_buf_set;

        /* Assign buffers to next wave from alternate buffer set */
        for (size_t w = 0; w < next_count; w++) {
            size_t slot = next_buf_set * wave_size + w;
            chunk_args[next_start + w].output_buffer = wave_buffers[slot];
            chunk_args[next_start + w].buffer_size = wave_sizes[slot];
            chunk_args[next_start + w].buffer_slot = &wave_buffers[slot];
            chunk_args[next_start + w].size_slot = &wave_sizes[slot];
            chunk_args[next_start + w].bytes_written = 0;
            chunk_args[next_start + w].success = 0;
        }

        /* Submit format tasks for next wave to pool */
        for (size_t w = 0; w < next_count; w++) {
            ctools_persistent_pool_submit(pool, format_chunk_thread, &chunk_args[next_start + w]);
        }

        /* Write current wave from main thread using batched I/O */
        ssize_t batch_written = cexport_io_write_batch(&outfile, batch_entries, wave_count);

        /* Wait for format tasks */
        ctools_persistent_pool_wait(pool);

        /* Check for write failure (a short write would leave a hole) */
        if (batch_written < 0 || (size_t)batch_written != wave_bytes) {
            for (size_t u = 0; u < total_buffers; u++) {
                if (wave_buffers[u]) free(wave_buffers[u]);
            }
            snprintf(msg, sizeof(msg), "cexport: write error during parallel I/O: %s\n",
                     outfile.error_message);
            cexport_io_close(&outfile, 0);
            SF_error(msg);
            return 693;
        }

        /* Check for format failures in next wave */
        for (size_t w = 0; w < next_count; w++) {
            if (!chunk_args[next_start + w].success) {
                for (size_t u = 0; u < total_buffers; u++) {
                    if (wave_buffers[u]) free(wave_buffers[u]);
                }
                cexport_io_close(&outfile, 0);
                SF_error(chunk_args[next_start + w].read_error ? "cexport: failed to read Stata data\n"
                                                               : "cexport: formatting failed\n");
                return 920;
            }
        }

        /* Advance to next wave */
        wave_start = next_start;
        wave_count = next_count;
        cur_buf_set = next_buf_set;
    }

    size_t total_file_size = file_offset;

    if (cexport_io_close(&outfile, total_file_size) != 0) {
        SF_error("cexport: failed to close file\n");
        ret = 693;
    }

    /* Free large data buffers (arena handles metadata) */
    for (size_t u = 0; u < total_buffers; u++) {
        if (wave_buffers[u]) free(wave_buffers[u]);
    }
    return ret;
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

    /* Types are loaded together with the names. */
    if (!g_ctx.vartypes) {
        SF_error("cexport: failed to load variable types\n");
        cexport_context_cleanup(&g_ctx);
        return 920;
    }

    /* Format cells straight from Stata's observation-major storage unless a
     * strL must be read on this thread; the mmap backend keeps fixed-size
     * destinations and so formats from a loaded copy. */
    g_ctx.is_string = (unsigned char *)malloc(nvars);
    if (g_ctx.is_string == NULL) {
        cexport_context_cleanup(&g_ctx);
        return 920;
    }
    bool longtext = false;
    g_ctx.all_numeric = true;
    for (size_t j = 0; j < nvars; j++) {
        ST_int var = (ST_int)(j + 1);
        g_ctx.is_string[j] = SF_var_is_string(var) != 0;
        if (g_ctx.is_string[j]) {
            g_ctx.all_numeric = false;
            if (SF_var_is_strl(var)) longtext = true;
        }
    }
    g_ctx.stream = g_ctx.io_backend != CEXPORT_IO_MMAP && !longtext;

    ctools_filtered_data_init(&g_ctx.filtered);
    rc = g_ctx.stream ? ctools_data_select(&g_ctx.filtered, 0, 0, CTOOLS_LOAD_CHECK_IF)
       : longtext ? cexport_load_text_data(&g_ctx.filtered, SF_nvars(), false)
                  : ctools_data_load(&g_ctx.filtered, NULL, 0, 0, 0, CTOOLS_LOAD_CHECK_IF | CTOOLS_LOAD_NO_SORT_ORDER);
    if (rc != STATA_OK) {
        snprintf(msg, sizeof(msg), "cexport: failed to load data (error %d)\n", rc);
        SF_error(msg);
        cexport_context_cleanup(&g_ctx);
        return longtext ? rc : 920;
    }

    /* Use filtered count: only rows that passed SF_ifobs */
    nobs = g_ctx.filtered.data.nobs;
    g_ctx.nobs_loaded = nobs;

    g_ctx.time_load = ctools_timer_seconds() - t_phase;

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

    size_t avg_row_size = sample_row_sizes(nobs, longtext);
    if (avg_row_size == 0) {
        if (nvars > SIZE_MAX / 22) {
            SF_error("cexport: variable count too large\n");
            cexport_context_cleanup(&g_ctx);
            return 920;
        }
        avg_row_size = nvars * 22;
    }
    g_ctx.actual_avg_row_size = avg_row_size;

    /* Bound each chunk by bytes as well as rows. With wide rows a 50,000-row
     * chunk is hundreds of MB, 2 x threads such buffers are live at once, and
     * a write wave grows past what one write call accepts. */
    if (avg_row_size > 0 && chunk_size > CEXPORT_MAX_CHUNK_BYTES / avg_row_size) {
        chunk_size = CEXPORT_MAX_CHUNK_BYTES / avg_row_size;
        if (chunk_size < 1) chunk_size = 1;
        g_ctx.adaptive_chunk_size = chunk_size;
        num_chunks = (nobs + chunk_size - 1) / chunk_size;
        num_threads = num_chunks;
        if (num_threads > max_threads) num_threads = max_threads;
        if (num_threads < 1) num_threads = 1;
    }

    size_t chunk_buffer_size;
    if (avg_row_size > SIZE_MAX / chunk_size) {
        SF_error("cexport: buffer size overflow\n");
        cexport_context_cleanup(&g_ctx);
        return 920;
    }
    chunk_buffer_size = chunk_size * avg_row_size * 12 / 10 + 4096;
    if (chunk_buffer_size < row_reserve_bytes()) chunk_buffer_size = row_reserve_bytes();

    /* Format header */
    char *header_buf = NULL;
    size_t header_len = 0;
    if (g_ctx.write_header) {
        size_t header_capacity = 3;
        for (size_t j = 0; j < nvars; j++) {
            size_t bytes = 2 * strlen(g_ctx.varnames[j]) + 3;
            if (bytes > SIZE_MAX - header_capacity) {
                cexport_context_cleanup(&g_ctx); return 920;
            }
            header_capacity += bytes;
        }
        header_buf = (char *)malloc(header_capacity);
        if (header_buf == NULL) {
            SF_error("cexport: memory allocation failed for header\n");
            cexport_context_cleanup(&g_ctx);
            return 920;
        }
        int hlen = cexport_format_header(&g_ctx, header_buf, header_capacity);
        if (hlen < 0) {
            SF_error("cexport: failed to format header\n");
            free(header_buf);
            cexport_context_cleanup(&g_ctx);
            return 920;
        }
        header_len = (size_t)hlen;
    }

    cexport_output output;
    rc = cexport_output_prepare(&output, g_ctx.filename, g_ctx.replace);
    if (rc) { free(header_buf); cexport_context_cleanup(&g_ctx); return rc; }
    free(g_ctx.filename);
    g_ctx.filename = strdup(output.temporary);
    if (!g_ctx.filename) {
        free(header_buf); cexport_output_cleanup(&output);
        cexport_context_cleanup(&g_ctx); return 920;
    }

    /* Choose export path */
    ST_retcode export_rc;

    if (!g_ctx.use_parallel_io || nobs < CTOOLS_EXPORT_CHUNK_SIZE || num_threads == 1 || avg_row_size>1024*1024) {
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
            char *header_copy = header_len ? (char *)malloc(header_len) : NULL;
            if (header_len && header_copy == NULL) {
                free(header_buf);
                export_rc = 920;
            } else {
                if (header_copy) memcpy(header_copy, header_buf, header_len);
                export_rc = export_mmap(nobs, avg_row_size, chunk_size, num_chunks,
                                         header_buf, header_len, pool);
                if (export_rc == CEXPORT_RETRY_PWRITE) {
                    ctools_arena_free(&g_ctx.chunk_arena);
                    g_ctx.io_backend = CEXPORT_IO_PWRITE;
                    export_rc = export_pwrite(nobs, avg_row_size, chunk_size, num_chunks,
                                               chunk_buffer_size, header_copy, header_len, pool);
                } else {
                    free(header_copy);
                }
            }
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

    /* Close before publication; a late buffered-write error is fatal. */
    if (g_ctx.fp) {
        if (fclose(g_ctx.fp) && !export_rc) export_rc = 603;
        g_ctx.fp = NULL;
    }
    if (!export_rc) export_rc = cexport_output_commit(&output);
    cexport_output_cleanup(&output);
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

/* ============================================================================
 * Cleanup function for ctools_cleanup system
 * ============================================================================ */

void cexport_cleanup_state(void)
{
    /* cexport uses a static context that is cleaned up at the end of each
     * operation, so this is typically a no-op. Call cleanup for safety
     * in case of interrupted operations. */
    cexport_context_cleanup(&g_ctx);
}

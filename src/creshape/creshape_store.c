/*
 * creshape_store.c - fill the reshaped dataset.
 *
 * Stata keeps one record per observation, so the cost of a store is set by
 * how often each record is fetched. The writers below follow the staged
 * transport (ctools_rowstage.c): each worker owns a contiguous range of
 * result observations and writes every variable of a tile of them before
 * moving on, so a record stays in cache for all of its variables. A tile
 * gathers its source rows first (one fetch per source row, prefetched),
 * then every pass over it runs at regular strides.
 *
 * Long: result observations (r, l), for source row r in order of i and
 * level l, are consecutive, so a tile is a run of source rows times the
 * levels. Wide: result observation g is one group of i; its cells name the
 * source row of each level.
 *
 * Numeric writes use the checked SPI store (see ctools_transport_internal.h:
 * unchecked stores from several threads corrupt columns).
 */
#include "ctools_transport_internal.h"
#include "ctools_rowstage.h"
#include "creshape_internal.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* Bytes of Stata records a tile should cover, and its bounds in result
 * observations (the staged transport's measured values). */
#define STORE_TILE_BYTES ((size_t)64 * 1024)
#define STORE_TILE_MAX_ROWS 256
#define STORE_TILE_MIN_ROWS 16
/* SPI calls a worker should have before another worker is worth waking. */
#define STORE_CELLS_PER_THREAD 65536
/* Bytes of a gathered source row fetched ahead of its use. */
#define STORE_PREFETCH_BYTES 2048

static int store_threads(size_t nobs, size_t nvars)
{
#ifdef _OPENMP
    size_t cells;
    int threads = ctools_get_max_threads();
    if (threads < 1)
        threads = 1;
    if (ctools_safe_mul_size(nobs, nvars, &cells) == 0 &&
        cells / STORE_CELLS_PER_THREAD < (size_t)threads)
        threads = (int)(cells / STORE_CELLS_PER_THREAD);
    return threads < 1 ? 1 : threads;
#else
    (void)nobs;
    (void)nvars;
    return 1;
#endif
}

static void store_share(size_t n, size_t team, size_t tid, size_t *start, size_t *end)
{
    size_t base = n / team, extra = n % team;
    *start = tid * base + (tid < extra ? tid : extra);
    *end = *start + base + (tid < extra);
}

/* Result observations per tile, from the width of a result record. */
static size_t store_tile_rows(const creshape_state *st)
{
    size_t bytes = 0;
    for (size_t v = 0; v < st->out_nvars; v++)
    {
        const creshape_out *o = &st->out[v];
        size_t width = 8;
        if (o->is_string)
        {
            width = 32;
            if (o->kind == CRESHAPE_OUT_COPY)
                width = (size_t)st->staged.width[o->col] + 1;
        }
        bytes += width;
    }
    size_t rows = STORE_TILE_BYTES / (bytes ? bytes : 1);
    if (rows > STORE_TILE_MAX_ROWS)
        rows = STORE_TILE_MAX_ROWS;
    if (rows < STORE_TILE_MIN_ROWS)
        rows = STORE_TILE_MIN_ROWS;
    return rows;
}

static inline void store_prefetch(const char *row, size_t stride)
{
#if defined(__GNUC__) || defined(__clang__)
    size_t bytes = stride < STORE_PREFETCH_BYTES ? stride : STORE_PREFETCH_BYTES;
    for (size_t at = 0; at < bytes; at += CACHE_LINE_SIZE)
        __builtin_prefetch(row + at, 0, 0);
#else
    (void)row;
    (void)stride;
#endif
}

/* ------------------------------------------------------------------------
 * Long
 * ------------------------------------------------------------------------ */

static stata_retcode store_long(const creshape_state *st)
{
    const ctools_rowstage *rs = &st->staged;
    const size_t nobs = st->nobs, nlevels = st->nlevels, nvars = st->out_nvars;
    const size_t stride = rs->stride;
    const perm_idx_t *order = st->order;
    size_t tile_in = store_tile_rows(st) / nlevels;
    if (tile_in < 1)
        tile_in = 1;
    const int threads = store_threads(st->out_nobs, nvars);
    const double missing = SV_missval;
    atomic_int error = 0;
    (void)threads;
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        size_t team = 1, tid = 0, start, end;
#ifdef _OPENMP
        team = (size_t)omp_get_num_threads();
        tid = (size_t)omp_get_thread_num();
#endif
        store_share(nobs, team, tid, &start, &end);
        const char *rows[STORE_TILE_MAX_ROWS];
        ST_IIID vstore = IO_VSTORE_FN;
        ST_IIIS sstore = (_stata_)->sstore;
        char empty[1] = "";
        for (size_t tile = start; tile < end && !atomic_load(&error); tile += tile_in)
        {
            const size_t count = end - tile < tile_in ? end - tile : tile_in;
            for (size_t i = 0; i < count; i++)
            {
                rows[i] = rs->rows + (order ? (size_t)order[tile + i] : tile + i) * stride;
                if (order)
                    store_prefetch(rows[i], stride);
            }
            const ST_int obs0 = (ST_int)(tile * nlevels) + 1;
            for (size_t v = 0; v < nvars; v++)
            {
                const creshape_out *o = &st->out[v];
                const ST_int idx = (ST_int)v + 1;
                ST_int obs = obs0;
                int io_rc = 0;
                if (o->kind == CRESHAPE_OUT_COPY)
                {
                    const size_t off = rs->offset[o->col];
                    if (!o->is_string)
                    {
                        for (size_t i = 0; i < count && !io_rc; i++)
                        {
                            const double value = *(const double *)(const void *)(rows[i] + off);
                            for (size_t l = 0; l < nlevels && !io_rc; l++)
                                io_rc = vstore(idx, obs++, value);
                        }
                    }
                    else
                    {
                        for (size_t i = 0; i < count && !io_rc; i++)
                            for (size_t l = 0; l < nlevels && !io_rc; l++)
                                io_rc = sstore(idx, obs++, (char *)(uintptr_t)(rows[i] + off));
                    }
                }
                else if (o->kind == CRESHAPE_OUT_J)
                {
                    const creshape_level *levels = st->levels;
                    for (size_t i = 0; i < count && !io_rc; i++)
                        for (size_t l = 0; l < nlevels && !io_rc; l++)
                            io_rc = o->is_string ? sstore(idx, obs++, levels[l].text)
                                                 : vstore(idx, obs++, levels[l].num);
                }
                else
                {
                    const int *cols = st->stubs[o->col].cols;
                    if (!o->is_string)
                    {
                        for (size_t i = 0; i < count && !io_rc; i++)
                            for (size_t l = 0; l < nlevels && !io_rc; l++)
                            {
                                const int c = cols[l];
                                io_rc = vstore(
                                    idx, obs++,
                                    c < 0
                                        ? missing
                                        : *(const double *)(const void *)(rows[i] + rs->offset[c]));
                            }
                    }
                    else
                    {
                        for (size_t i = 0; i < count && !io_rc; i++)
                            for (size_t l = 0; l < nlevels && !io_rc; l++)
                            {
                                const int c = cols[l];
                                io_rc = sstore(
                                    idx, obs++,
                                    c < 0 ? empty : (char *)(uintptr_t)(rows[i] + rs->offset[c]));
                            }
                    }
                }
                if (io_rc)
                {
                    ctools_io_record_io_error(&error,
                                              ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE));
                    break;
                }
            }
        }
    }
    ctools_memory_barrier();
    return ctools_io_store_status((stata_retcode)atomic_load(&error));
}

/* ------------------------------------------------------------------------
 * Wide
 * ------------------------------------------------------------------------ */

static stata_retcode store_wide(const creshape_state *st)
{
    const ctools_rowstage *rs = &st->staged;
    const size_t ngroups = st->ngroups, nlevels = st->nlevels, nvars = st->out_nvars;
    const size_t stride = rs->stride;
    const perm_idx_t *order = st->order, *starts = st->group_start, *cell = st->cell;
    const size_t tile_rows = store_tile_rows(st);
    const int threads = store_threads(ngroups, nvars);
    const double missing = SV_missval;
    atomic_int error = 0;
    (void)threads;
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        size_t team = 1, tid = 0, start, end;
#ifdef _OPENMP
        team = (size_t)omp_get_num_threads();
        tid = (size_t)omp_get_thread_num();
#endif
        store_share(ngroups, team, tid, &start, &end);
        const char *first[STORE_TILE_MAX_ROWS];
        ST_IIID vstore = IO_VSTORE_FN;
        ST_IIIS sstore = (_stata_)->sstore;
        char empty[1] = "";
        for (size_t tile = start; tile < end && !atomic_load(&error); tile += tile_rows)
        {
            const size_t count = end - tile < tile_rows ? end - tile : tile_rows;
            /* the first row of each group carries i and the constant
             * variables; fetch the rows of the tile's cells ahead */
            for (size_t i = 0; i < count; i++)
            {
                const size_t p = starts[tile + i];
                first[i] = rs->rows + (order ? (size_t)order[p] : p) * stride;
                if (order)
                {
                    const perm_idx_t *cells = cell + (tile + i) * nlevels;
                    for (size_t l = 0; l < nlevels; l++)
                        if (cells[l] != CTOOLS_ROW_NONE)
                            store_prefetch(rs->rows + (size_t)cells[l] * stride, stride);
                }
            }
            const ST_int obs0 = (ST_int)tile + 1;
            for (size_t v = 0; v < nvars; v++)
            {
                const creshape_out *o = &st->out[v];
                const ST_int idx = (ST_int)v + 1;
                int io_rc = 0;
                if (o->kind == CRESHAPE_OUT_COPY)
                {
                    const size_t off = rs->offset[o->col];
                    if (!o->is_string)
                        for (size_t i = 0; i < count && !io_rc; i++)
                            io_rc = vstore(idx, obs0 + (ST_int)i,
                                           *(const double *)(const void *)(first[i] + off));
                    else
                        for (size_t i = 0; i < count && !io_rc; i++)
                            io_rc =
                                sstore(idx, obs0 + (ST_int)i, (char *)(uintptr_t)(first[i] + off));
                }
                else
                {
                    const size_t off = rs->offset[st->stubs[o->col].long_col];
                    const perm_idx_t *cells = cell + tile * nlevels + (size_t)o->level;
                    if (!o->is_string)
                        for (size_t i = 0; i < count && !io_rc; i++, cells += nlevels)
                            io_rc = vstore(
                                idx, obs0 + (ST_int)i,
                                *cells == CTOOLS_ROW_NONE
                                    ? missing
                                    : *(const double *)(const void *)(rs->rows +
                                                                      (size_t)*cells * stride +
                                                                      off));
                    else
                        for (size_t i = 0; i < count && !io_rc; i++, cells += nlevels)
                            io_rc =
                                sstore(idx, obs0 + (ST_int)i,
                                       *cells == CTOOLS_ROW_NONE
                                           ? empty
                                           : (char *)(uintptr_t)(rs->rows +
                                                                 (size_t)*cells * stride + off));
                }
                if (io_rc)
                {
                    ctools_io_record_io_error(&error,
                                              ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE));
                    break;
                }
            }
        }
    }
    ctools_memory_barrier();
    return ctools_io_store_status((stata_retcode)atomic_load(&error));
}

stata_retcode creshape_store(const creshape_state *st)
{
    /* the dataset must be the one this plan described */
    if ((size_t)SF_nvars() != st->out_nvars || (size_t)SF_nobs() != st->out_nobs)
        return STATA_ERR_INVALID_INPUT;
    if (!st->out_nobs || !st->out_nvars)
        return STATA_OK;
    stata_retcode rc = ctools_io_validate_range(1, st->out_nobs);
    for (size_t v = 0; v < st->out_nvars && !rc; v++)
        rc = ctools_io_validate_variable((int)v + 1, st->out[v].is_string);
    if (rc)
        return rc;
    return st->dir == CRESHAPE_LONG ? store_long(st) : store_wide(st);
}

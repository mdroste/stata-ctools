#include "ctools_transport_internal.h"
#include "ctools_rowstage.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* A tile is the run of observations one worker moves between passes over the
 * variables: small enough that the Stata records and their staged rows stay
 * cache-resident for every pass, large enough to amortize the per-variable
 * loop setup. Measured flat between 64 and 512 rows for 200-byte rows. */
#define ROWSTAGE_TILE_BYTES ((size_t)64 * 1024)
#define ROWSTAGE_TILE_MAX_ROWS 256
#define ROWSTAGE_TILE_MIN_ROWS 16
/* SPI calls a worker should have before another worker is worth waking. */
#ifndef ROWSTAGE_CELLS_PER_THREAD
#define ROWSTAGE_CELLS_PER_THREAD 65536
#endif
/* Reordered stores read staged rows at random: fetch a few rows ahead, and
 * leave the tail of a very wide row to the sequential hardware prefetcher. */
#define ROWSTAGE_PREFETCH_ROWS 8
#define ROWSTAGE_PREFETCH_BYTES 1024
/* Staged strings occupy their declared width. Beyond the loader's bound for
 * contiguous string columns, the column transport packs actual lengths and can
 * need far less memory, so such data are left to it. */
#define ROWSTAGE_MAX_STRING_COLUMN_BYTES ((size_t)2048 * 1024 * 1024)

void ctools_rowstage_init(ctools_rowstage *rs)
{
    if (rs)
        memset(rs, 0, sizeof(*rs));
}

void ctools_rowstage_free(ctools_rowstage *rs)
{
    if (!rs)
        return;
    free(rs->rows);
    free(rs->var_idx);
    free(rs->width);
    free(rs->offset);
    ctools_rowstage_init(rs);
}

static size_t rowstage_tile_rows(size_t stride)
{
    size_t rows = ROWSTAGE_TILE_BYTES / stride;
    if (rows > ROWSTAGE_TILE_MAX_ROWS)
        rows = ROWSTAGE_TILE_MAX_ROWS;
    if (rows < ROWSTAGE_TILE_MIN_ROWS)
        rows = ROWSTAGE_TILE_MIN_ROWS;
    return rows;
}

static int rowstage_threads(size_t nobs, size_t nvars)
{
#ifdef _OPENMP
    size_t cells;
    int threads = ctools_get_max_threads();
    if (threads < 1)
        threads = 1;
    if (ctools_safe_mul_size(nobs, nvars, &cells) == 0 &&
        cells / ROWSTAGE_CELLS_PER_THREAD < (size_t)threads)
        threads = (int)(cells / ROWSTAGE_CELLS_PER_THREAD);
    return threads < 1 ? 1 : threads;
#else
    (void)nobs;
    (void)nvars;
    return 1;
#endif
}

/* Contiguous share of [0, n) for one worker of a team. */
static void rowstage_share(size_t n, size_t team, size_t tid, size_t *start, size_t *end)
{
    size_t base = n / team, extra = n % team;
    *start = tid * base + (tid < extra ? tid : extra);
    *end = *start + base + (tid < extra);
}

/* Lay out one staged row for the listed variables and nobs observations.
 * Allocates the column metadata only; on failure nothing is left allocated. A
 * NULL str_widths reads the caller's _ctools_strw metadata. */
static stata_retcode rowstage_plan(ctools_rowstage *rs, const int *var_indices, size_t nvars,
                                   const int *str_widths, size_t nobs)
{
    ctools_rowstage_init(rs);
    if (!var_indices || nvars == 0)
        return STATA_ERR_INVALID_INPUT;
    int *owned_widths = NULL;
    stata_retcode rc = STATA_OK;
    rs->var_idx = (int *)ctools_safe_malloc2(nvars, sizeof(int));
    rs->width = (int *)ctools_safe_malloc2(nvars, sizeof(int));
    rs->offset = (size_t *)ctools_safe_malloc2(nvars, sizeof(size_t));
    if (!rs->var_idx || !rs->width || !rs->offset)
        rc = STATA_ERR_MEMORY;

    /* Numeric columns come first so every staged double is 8-byte aligned. */
    size_t stride = 0;
    for (int strings = 0; strings <= 1 && !rc; strings++)
    {
        for (size_t j = 0; j < nvars && !rc; j++)
        {
            int idx = var_indices[j];
            if (idx < 1 || idx > SF_nvars())
            {
                rc = STATA_ERR_INVALID_INPUT;
                break;
            }
            int kind = ctools_io_variable_kind(idx);
            if ((kind != IO_NUMERIC) != strings)
                continue;
            size_t bytes = sizeof(double);
            rs->width[j] = 0;
            if (strings)
            {
                if (!str_widths && !owned_widths)
                    str_widths = owned_widths =
                        ctools_io_ctools_read_strw_from_stata((size_t)SF_nvars());
                int width = str_widths ? str_widths[idx - 1] : 0;
                if (kind != IO_FIXED_STRING || width < 1 || width >= STATA_STR_MAXLEN ||
                    nobs > ROWSTAGE_MAX_STRING_COLUMN_BYTES / ((size_t)width + 1))
                {
                    rc = STATA_ERR_UNSUPPORTED_TYPE;
                    break;
                }
                rs->width[j] = width;
                bytes = (size_t)width + 1;
            }
            rs->var_idx[j] = idx;
            rs->offset[j] = stride;
            if (stride > SIZE_MAX - bytes - 8)
                rc = STATA_ERR_MEMORY;
            stride += bytes;
        }
    }
    free(owned_widths);
    if (rc)
    {
        ctools_rowstage_free(rs);
        return rc;
    }
    rs->stride = (stride + 7) & ~(size_t)7;
    rs->nvars = nvars;
    return STATA_OK;
}

int ctools_rowstage_supported(const int *var_indices, size_t nvars, const int *str_widths,
                              size_t nobs)
{
    ctools_rowstage plan;
    stata_retcode rc = rowstage_plan(&plan, var_indices, nvars, str_widths, nobs);
    ctools_rowstage_free(&plan);
    return rc == STATA_OK;
}

stata_retcode ctools_rowstage_load(ctools_rowstage *rs, const int *var_indices, size_t nvars,
                                   const int *str_widths, size_t obs1, size_t nobs)
{
    if (!rs)
        return STATA_ERR_INVALID_INPUT;
    ctools_rowstage_init(rs);
    if (nobs == 0)
        return STATA_ERR_INVALID_INPUT;
    stata_retcode rc = ctools_io_validate_range(obs1, nobs);
    if (!rc)
        rc = rowstage_plan(rs, var_indices, nvars, str_widths, nobs);
    if (rc)
        return rc;
    const size_t stride = rs->stride;

    /* The guard keeps a direct string read at the last rows inside the buffer. */
    size_t bytes;
    if (ctools_safe_mul_size(nobs, stride, &bytes) != 0 || bytes > SIZE_MAX - IO_SLOT_GUARD)
    {
        ctools_rowstage_free(rs);
        return STATA_ERR_MEMORY;
    }
    rs->rows = (char *)malloc(bytes + IO_SLOT_GUARD);
    if (!rs->rows)
    {
        ctools_rowstage_free(rs);
        return STATA_ERR_MEMORY;
    }
    rs->nobs = nobs;

    const size_t tile_rows = rowstage_tile_rows(stride);
    const int threads = rowstage_threads(nobs, nvars);
    atomic_int error = 0;
    (void)threads;
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        size_t team = 1, tid = 0, start, end;
#ifdef _OPENMP
        team = (size_t)omp_get_num_threads();
        tid = (size_t)omp_get_thread_num();
#endif
        rowstage_share(nobs, team, tid, &start, &end);
        /* Each worker owns its rows; only the last one may read into the guard. */
        const char *safe_end = rs->rows + end * stride + (end == nobs ? IO_SLOT_GUARD : 0);
        for (size_t tile = start; tile < end && !atomic_load(&error); tile += tile_rows)
        {
            size_t count = end - tile < tile_rows ? end - tile : tile_rows;
            for (size_t j = 0; j < nvars; j++)
            {
                ST_int idx = (ST_int)rs->var_idx[j];
                ST_int obs = (ST_int)(obs1 + tile);
                char *slot = rs->rows + tile * stride + rs->offset[j];
                int io_rc = 0;
                if (rs->width[j] == 0)
                {
                    ST_IIIDp read_fn = IO_VDATA_FN;
                    for (size_t i = 0; i < count; i++, slot += stride)
                    {
                        int read_rc = read_fn(idx, obs + (ST_int)i, (double *)(void *)slot);
                        if (read_rc)
                        {
                            io_rc = ctools_spi_status(read_rc, STATA_ERR_STATA_READ);
                            break;
                        }
                    }
                }
                else
                {
                    ST_IIIS read_fn = (_stata_)->sdata;
                    size_t width = (size_t)rs->width[j];
                    for (size_t i = 0; i < count && !io_rc; i++, slot += stride)
                        io_rc = ctools_io_read_string_slot_direct(read_fn, idx, obs + (ST_int)i,
                                                                  slot, width, safe_end);
                }
                if (io_rc)
                {
                    ctools_io_record_io_error(&error, io_rc);
                    break;
                }
            }
        }
    }
    ctools_memory_barrier();
    rc = (stata_retcode)atomic_load(&error);
    if (rc)
    {
        ctools_rowstage_free(rs);
        return ctools_io_load_status(rc);
    }
    return STATA_OK;
}

stata_retcode ctools_rowstage_column(const ctools_rowstage *rs, size_t col, stata_variable *out)
{
    if (!rs || !out || !rs->rows || col >= rs->nvars)
        return STATA_ERR_INVALID_INPUT;
    memset(out, 0, sizeof(*out));
    const size_t nobs = rs->nobs, stride = rs->stride;
    const char *first = rs->rows + rs->offset[col];
    const int threads = rowstage_threads(nobs, 4);
    (void)threads;
    out->nobs = nobs;

    if (rs->width[col] == 0)
    {
        double *values = (double *)ctools_safe_cacheline_alloc2(nobs, sizeof(double));
        if (!values)
            return STATA_ERR_MEMORY;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
        for (size_t i = 0; i < nobs; i++)
            values[i] = *(const double *)(const void *)(first + i * stride);
        out->type = STATA_TYPE_DOUBLE;
        out->data.dbl = values;
        return STATA_OK;
    }

    /* Same contiguous slot layout as the loader's fixed-width string columns;
     * calloc keeps every slot zero-padded after its terminator. */
    const size_t width = (size_t)rs->width[col], slot = width + 1;
    char **strings = (char **)ctools_safe_cacheline_alloc2(nobs, sizeof(char *));
    ctools_string_arena *arena = (ctools_string_arena *)calloc(1, sizeof(*arena));
    char *flat = (char *)calloc(nobs + (IO_SLOT_GUARD + slot - 1) / slot, slot);
    if (!strings || !arena || !flat)
    {
        ctools_aligned_free(strings);
        free(arena);
        free(flat);
        return STATA_ERR_MEMORY;
    }
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (size_t i = 0; i < nobs; i++)
    {
        const char *value = first + i * stride;
        strings[i] = flat + i * slot;
        memcpy(strings[i], value, strnlen(value, width));
    }
    arena->base = flat;
    arena->capacity = arena->used = nobs * slot;
    arena->mode = CTOOLS_STRING_ARENA_NO_FALLBACK;
    out->type = STATA_TYPE_STRING;
    out->str_maxlen = width;
    out->data.str = strings;
    out->_arena = arena;
    return STATA_OK;
}

stata_retcode ctools_rowstage_store(const ctools_rowstage *rs, const perm_idx_t *source,
                                    size_t n_out, size_t obs1)
{
    if (!rs || !rs->rows || rs->nvars == 0 || (!source && n_out > rs->nobs))
        return STATA_ERR_INVALID_INPUT;
    if (n_out == 0)
        return STATA_OK;
    stata_retcode rc = ctools_io_validate_range(obs1, n_out);
    for (size_t j = 0; j < rs->nvars && !rc; j++)
        rc = ctools_io_validate_variable(rs->var_idx[j], rs->width[j] > 0);
    if (rc)
        return rc;

    const size_t nobs = rs->nobs, nvars = rs->nvars, stride = rs->stride;
    const int threads = rowstage_threads(n_out, nvars);
    (void)threads;
    if (source)
    {
        int bad = 0;
#pragma omp parallel for schedule(static) reduction(| : bad) num_threads(threads) if (threads > 1)
        for (size_t i = 0; i < n_out; i++)
            bad |= source[i] >= nobs && source[i] != CTOOLS_ROW_NONE;
        if (bad)
            return STATA_ERR_INVALID_INPUT;
    }

    /* What a destination without a source row receives. An in-order store
     * writes straight from the staged rows and allocates nothing, so it can
     * still restore the original observations after a memory failure. */
    char *blank = NULL;
    if (source)
    {
        blank = (char *)calloc(1, stride);
        if (!blank)
            return STATA_ERR_MEMORY;
        for (size_t j = 0; j < nvars; j++)
            if (rs->width[j] == 0)
                *(double *)(void *)(blank + rs->offset[j]) = SV_missval;
    }

    const size_t tile_rows = rowstage_tile_rows(stride);
    const size_t prefetch_bytes =
        stride < ROWSTAGE_PREFETCH_BYTES ? stride : ROWSTAGE_PREFETCH_BYTES;
    atomic_int error = 0;
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        size_t team = 1, tid = 0, start, end;
#ifdef _OPENMP
        team = (size_t)omp_get_num_threads();
        tid = (size_t)omp_get_thread_num();
#endif
        rowstage_share(n_out, team, tid, &start, &end);
        /* Gather a tile of staged rows, then hand each variable to SPI from
         * it: one random row fetch per observation, regular strides after. */
        char *tile_rows_buf = source ? (char *)ctools_safe_malloc2(tile_rows, stride) : NULL;
        if (source && !tile_rows_buf)
            ctools_io_record_io_error(&error, STATA_ERR_MEMORY);
        for (size_t tile = start; (tile_rows_buf || !source) && tile < end && !atomic_load(&error);
             tile += tile_rows)
        {
            size_t count = end - tile < tile_rows ? end - tile : tile_rows;
            char *base = source ? tile_rows_buf : rs->rows + tile * stride;
            for (size_t i = 0; source && i < count; i++)
            {
                size_t row = source[tile + i];
#if defined(__GNUC__) || defined(__clang__)
                if (tile + i + ROWSTAGE_PREFETCH_ROWS < end)
                {
                    perm_idx_t ahead = source[tile + i + ROWSTAGE_PREFETCH_ROWS];
                    if (ahead != CTOOLS_ROW_NONE)
                        for (size_t byte = 0; byte < prefetch_bytes; byte += CACHE_LINE_SIZE)
                            __builtin_prefetch(rs->rows + (size_t)ahead * stride + byte, 0, 0);
                }
#endif
                memcpy(tile_rows_buf + i * stride,
                       row == CTOOLS_ROW_NONE ? blank : rs->rows + row * stride, stride);
            }
            for (size_t j = 0; j < nvars; j++)
            {
                ST_int idx = (ST_int)rs->var_idx[j];
                ST_int obs = (ST_int)(obs1 + tile);
                char *slot = base + rs->offset[j];
                int io_rc = 0;
                if (rs->width[j] == 0)
                {
                    ST_IIID store_fn = IO_VSTORE_FN;
                    for (size_t i = 0; i < count && !io_rc; i++, slot += stride)
                        io_rc = store_fn(idx, obs + (ST_int)i, *(const double *)(const void *)slot);
                }
                else
                {
                    ST_IIIS store_fn = (_stata_)->sstore;
                    for (size_t i = 0; i < count && !io_rc; i++, slot += stride)
                        io_rc = store_fn(idx, obs + (ST_int)i, slot);
                }
                if (io_rc)
                {
                    ctools_io_record_io_error(&error,
                                              ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE));
                    break;
                }
            }
        }
        free(tile_rows_buf);
    }
    free(blank);
    ctools_memory_barrier();
    rc = (stata_retcode)atomic_load(&error);
    /* An allocation failure before a worker's first write is still reported as
     * a failed store: other workers may already have written. */
    return ctools_io_store_status(rc);
}

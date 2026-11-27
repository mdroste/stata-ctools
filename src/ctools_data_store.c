#include "ctools_transport_internal.h"
#include <math.h>

stata_retcode ctools_io_store_status(int rc)
{
    if (rc)
        SF_error("ctools: write failed; data may be partially modified\n");
    return (stata_retcode)rc;
}

int ctools_io_store_single_variable(stata_variable *var, int var_idx, size_t obs1, size_t nobs,
                                    const perm_idx_t *order)
{
    size_t i;
    ST_int stata_var_idx = (ST_int)var_idx;

    /* Fuse gather and SPI store. Loaded input columns remain immutable, so
     * sequential destination writes need no full-column permutation buffer. */
    if (order)
    {
        if (var->type == STATA_TYPE_DOUBLE)
        {
            const double *values = var->data.dbl;
            ST_IIID store_fn = IO_VSTORE_FN;
            /* For source columns up to 2 MiB, directly consume the gather
             * instead of copying it through a second buffer. Use the full
             * source size, not the size of this tile/chunk. Large random
             * gathers retain their bounded staging buffer. */
            if (var->nobs <= IO_DIRECT_GATHER_MAX_ROWS)
            {
                for (i = 0; i < nobs; i++)
                {
#if defined(__GNUC__) || defined(__clang__)
                    if (i + 16 < nobs)
                        __builtin_prefetch(values + order[i + 16], 0, 0);
#endif
                    int io_rc = store_fn(stata_var_idx, (ST_int)(i + obs1), values[order[i]]);
                    if (io_rc)
                        return ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
                }
                return STATA_OK;
            }
            double gathered[1024];
            for (size_t start = 0; start < nobs; start += 1024)
            {
                size_t count = nobs - start;
                if (count > 1024)
                    count = 1024;
                for (i = 0; i < count; i++)
                {
#if defined(__GNUC__) || defined(__clang__)
                    if (start + i + 16 < nobs)
                        __builtin_prefetch(values + order[start + i + 16], 0, 0);
#endif
                    gathered[i] = values[order[start + i]];
                }
                for (i = 0; i < count; i++)
                {
                    {
                        int io_rc =
                            store_fn(stata_var_idx, (ST_int)(start + i + obs1), gathered[i]);
                        if (io_rc)
                            return ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
                    }
                }
            }
        }
        else
        {
            char *const *values = var->data.str;
            ST_IIIS store_fn = (_stata_)->sstore;
            /* Keep the gathered strings cache-resident while SPI consumes
             * them. Scratch space is bounded independently of row count. */
            if (var->str_maxlen > 0 && var->str_maxlen < STATA_STR_MAXLEN)
            {
                char gathered[65536];
                size_t stride = var->str_maxlen + 1;
                size_t batch = sizeof(gathered) / stride;
                for (size_t start = 0; start < nobs; start += batch)
                {
                    size_t count = nobs - start;
                    if (count > batch)
                        count = batch;
                    for (i = 0; i < count; i++)
                    {
#if defined(__GNUC__) || defined(__clang__)
                        if (start + i + 16 < nobs)
                        {
                            const char *next = values[order[start + i + 16]];
                            if (next)
                                __builtin_prefetch(next, 0, 0);
                        }
#endif
                        const char *value = values[order[start + i]];
                        size_t len = value ? strlen(value) : 0;
                        if (len >= stride)
                            return STATA_ERR_STATA_WRITE;
                        if (len)
                            memcpy(gathered + i * stride, value, len);
                        gathered[i * stride + len] = '\0';
                    }
                    for (i = 0; i < count; i++)
                    {
                        {
                            int io_rc = store_fn(stata_var_idx, (ST_int)(start + i + obs1),
                                                 gathered + i * stride);
                            if (io_rc)
                                return ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
                        }
                    }
                }
                return STATA_OK;
            }
            for (i = 0; i < nobs; i++)
            {
#if defined(__GNUC__) || defined(__clang__)
                if (i + 16 < nobs)
                {
                    const char *next = values[order[i + 16]];
                    if (next)
                        __builtin_prefetch(next, 0, 0);
                }
#endif
                char *value = values[order[i]];
                {
                    int io_rc = store_fn(stata_var_idx, (ST_int)(i + obs1), value ? value : "");
                    if (io_rc)
                        return ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
                }
            }
        }
        return STATA_OK;
    }

    if (var->type == STATA_TYPE_DOUBLE)
    {
        const double *restrict dbl_data = var->data.dbl;

        /* Cache SPI function pointer — avoid reloading _stata_ per iteration */
        ST_IIID store_fn = IO_VSTORE_FN;
        for (i = 0; i + 4 <= nobs; i += 4)
        {
            {
                int io_rc = store_fn(stata_var_idx, (ST_int)(i + 0 + obs1), dbl_data[i + 0]);
                if (io_rc)
                    return ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
            }
            {
                int io_rc = store_fn(stata_var_idx, (ST_int)(i + 1 + obs1), dbl_data[i + 1]);
                if (io_rc)
                    return ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
            }
            {
                int io_rc = store_fn(stata_var_idx, (ST_int)(i + 2 + obs1), dbl_data[i + 2]);
                if (io_rc)
                    return ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
            }
            {
                int io_rc = store_fn(stata_var_idx, (ST_int)(i + 3 + obs1), dbl_data[i + 3]);
                if (io_rc)
                    return ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
            }
        }
        for (; i < nobs; i++)
        {
            {
                int io_rc = store_fn(stata_var_idx, (ST_int)(i + 0 + obs1), dbl_data[i + 0]);
                if (io_rc)
                    return ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
            }
        }
    }
    else
    {
        /* String variable */
        char *const *restrict str_data = var->data.str;
        size_t str_width = var->str_maxlen;

        /* Cache SPI function pointer — avoid reloading _stata_ per iteration */
        ST_IIIS sstore_fn = (_stata_)->sstore;

        /* Replace only the old repack path. Keep the existing direct fallback
         * for unknown widths and columns above its 256 MiB limit; neither
         * needs a new per-cell width branch or length scan. */
        if (str_width > 0 && str_width < STATA_STR_MAXLEN &&
            nobs <= (256ULL * 1024 * 1024) / (str_width + 1))
        {
            for (i = 0; i < nobs; i++)
            {
#if defined(__GNUC__) || defined(__clang__)
                if (i + 16 < nobs && str_data[i + 16])
                    __builtin_prefetch(str_data[i + 16], 0, 0);
#endif
                const char *value = str_data[i] ? str_data[i] : "";
                if (strlen(value) > str_width)
                    return STATA_ERR_STATA_WRITE;
                {
                    int io_rc = sstore_fn(stata_var_idx, (ST_int)(i + obs1), (char *)value);
                    if (io_rc)
                        return ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
                }
            }
            return STATA_OK;
        }

        for (i = 0; i < nobs; i++)
        {
#if defined(__GNUC__) || defined(__clang__)
            if (i + 16 < nobs && str_data[i + 16])
                __builtin_prefetch(str_data[i + 16], 0, 0);
#endif
            {
                int io_rc =
                    sstore_fn(stata_var_idx, (ST_int)(i + obs1), str_data[i] ? str_data[i] : "");
                if (io_rc)
                    return ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
            }
        }
    }
    return STATA_OK;
}

void *ctools_io_store_variable_thread(void *arg)
{
    ctools_var_io_args *args = (ctools_var_io_args *)arg;
    args->error = ctools_io_store_single_variable(args->var, args->var_idx, args->obs1, args->nobs,
                                                  args->store_order);
    return args->error ? (void *)1 : NULL;
}

void *ctools_io_store_chunk_thread(void *arg)
{
    ctools_store_chunk *a = (ctools_store_chunk *)arg;
    if (atomic_load(a->error))
        return NULL;
    stata_variable view = *a->var;
    if (!a->order)
    {
        if (view.type == STATA_TYPE_DOUBLE)
        {
            view.data.dbl += a->dest;
        }
        else
        {
            view.data.str += a->dest;
            size_t width = view.str_maxlen;
            if (width > 0 && width < STATA_STR_MAXLEN &&
                a->total_nobs > (256ULL * 1024 * 1024) / (width + 1))
                view.str_maxlen = 0;
        }
    }
    int rc = ctools_io_store_single_variable(&view, a->var_idx, a->obs1, a->count,
                                             a->order ? a->order + a->dest : NULL);
    if (rc)
    {
        ctools_io_record_io_error(a->error, rc);
        return (void *)1;
    }
    return NULL;
}

stata_retcode ctools_io_data_store_impl(stata_data *data, int *var_indices, size_t nvars,
                                        size_t obs1, const perm_idx_t *order)
{
    if (data == NULL)
    {
        return STATA_ERR_INVALID_INPUT;
    }

    size_t nobs = data->nobs;

    /* If no var_indices, store all variables sequentially */
    if (var_indices == NULL)
    {
        nvars = data->nvars;
    }
    else if (nvars == 0)
    {
        return STATA_ERR_INVALID_INPUT;
    }

    if (!data->vars || nvars > data->nvars)
        return STATA_ERR_INVALID_INPUT;
    stata_retcode rc = ctools_io_validate_range(obs1, nobs);
    if (rc)
        return rc;
    for (size_t j = 0; j < nvars; j++)
    {
        stata_variable *var = &data->vars[j];
        if (var->type != STATA_TYPE_DOUBLE && var->type != STATA_TYPE_STRING)
            return STATA_ERR_UNSUPPORTED_TYPE;
        if (var->nobs < nobs || (nobs && !var->data.dbl))
            return STATA_ERR_INVALID_INPUT;
        rc = ctools_io_validate_variable(var_indices ? var_indices[j] : (int)j + 1,
                                         var->type == STATA_TYPE_STRING);
        if (rc)
            return rc;
    }
    if (order)
    {
        for (size_t i = 0; i < nobs; i++)
            if (order[i] >= nobs)
                return STATA_ERR_INVALID_INPUT;
    }

    /* Wide rows: keep a bounded row range hot while writing each column.
     * Each worker owns a disjoint destination range, including sorted writes. */
    size_t string_bytes = 0;
    int store_tile_ok = 0;
#ifdef _OPENMP
    store_tile_ok = nvars >= IO_MIN_TILED_VARS && nobs >= MIN_OBS_PER_THREAD * 2;
#endif
    for (size_t j = 0; j < nvars && store_tile_ok; j++)
    {
        if (data->vars[j].type == STATA_TYPE_STRING)
        {
            size_t width = data->vars[j].str_maxlen;
            if (!width || width >= STATA_STR_MAXLEN)
                store_tile_ok = 0;
            else
                string_bytes += width;
        }
    }
    int few_columns = 0;
#ifdef _OPENMP
    few_columns = nvars >= 2 && nvars <= 4 && nobs >= MIN_OBS_PER_THREAD * 2 &&
                  ctools_get_max_threads() >= (int)(2 * nvars);
#endif
    /* Row tiles whenever the columns can be tiled: each observation's
     * memory is then written once for all its columns, not once per column.
     * Measured 2026-11-27 (1M-10M rows, 5-24 columns): numeric columns took
     * 0.27-0.41 of the column-chunk store's time, with string columns
     * among them 0.43-0.56 (string_bytes is no longer a condition, unless
     * IO_AUTO_TILES is 0: the rule before, for A/B measurements). */
    if (few_columns || (store_tile_ok && (IO_AUTO_TILES || string_bytes >= IO_MIN_TILED_STRING_BYTES)))
    {
        atomic_int row_error = 0;
#pragma omp parallel
        {
            size_t team = 1;
            size_t tid = 0;
#ifdef _OPENMP
            team = (size_t)omp_get_num_threads();
            tid = (size_t)omp_get_thread_num();
#endif
            size_t base = nobs / team, extra = nobs % team;
            size_t start = tid * base + (tid < extra ? tid : extra);
            size_t count = base + (tid < extra);
            for (size_t tile = start; tile < start + count; tile += IO_STORE_TILE_ROWS)
            {
                size_t tile_count = start + count - tile;
                if (tile_count > IO_STORE_TILE_ROWS)
                    tile_count = IO_STORE_TILE_ROWS;
                for (size_t j = 0; j < nvars; j++)
                {
                    stata_variable view = data->vars[j];
                    if (!order)
                    {
                        if (view.type == STATA_TYPE_DOUBLE)
                            view.data.dbl += tile;
                        else
                        {
                            view.data.str += tile;
                            size_t width = view.str_maxlen;
                            if (width > 0 && width < STATA_STR_MAXLEN &&
                                nobs > (256ULL * 1024 * 1024) / (width + 1))
                                view.str_maxlen = 0;
                        }
                    }
                    int rc = ctools_io_store_single_variable(
                        &view, var_indices ? var_indices[j] : (int)j + 1, obs1 + tile, tile_count,
                        order ? order + tile : NULL);
                    if (rc)
                    {
                        ctools_io_record_io_error(&row_error, rc);
                        break;
                    }
                }
            }
        }
        ctools_memory_barrier();
        return ctools_io_store_status(atomic_load(&row_error));
    }

    /*
        nvars == 1: row-parallel store (OpenMP over observations).
        For numeric variables, this is a simple parallel SF_vstore loop.
        For string variables, use the existing ctools_io_store_single_variable which
        already handles checked pointer-array writes.
    */
    if (nvars == 1)
    {
        stata_variable *var = &data->vars[0];
        int var_idx = var_indices ? var_indices[0] : 1;
        ST_int stata_var = (ST_int)var_idx;

        atomic_int error = 0;
        if (var->type == STATA_TYPE_DOUBLE)
        {
            const double *restrict dbl_data = var->data.dbl;
            /* Cache SPI function pointer — avoid reloading _stata_ per iteration */
            ST_IIID store_fn = IO_VSTORE_FN;
#pragma omp parallel for schedule(static) if (nobs >= MIN_OBS_PER_THREAD * 2)
            for (size_t i = 0; i < nobs; i++)
            {
                size_t row = order ? order[i] : i;
                {
                    int io_rc = store_fn(stata_var, (ST_int)(i + obs1), dbl_data[row]);
                    if (io_rc)
                        ctools_io_record_io_error(&error,
                                                  ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE));
                }
            }
        }
        else
        {
            /* String store: row-parallel from pointer array */
            char *const *restrict str_data = var->data.str;
            ST_IIIS sstore_fn = (_stata_)->sstore;
#pragma omp parallel for schedule(static) if (nobs >= MIN_OBS_PER_THREAD * 2)
            for (size_t i = 0; i < nobs; i++)
            {
                size_t row = order ? order[i] : i;
                {
                    int io_rc = sstore_fn(stata_var, (ST_int)(i + obs1),
                                          str_data[row] ? str_data[row] : "");
                    if (io_rc)
                        ctools_io_record_io_error(&error,
                                                  ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE));
                }
            }
        }

        ctools_memory_barrier();
        return ctools_io_store_status(atomic_load(&error));
    }

    /* nvars >= 2: chunked column store when enough rows are available for
     * the queue to balance unequal or indivisible column counts. */
    size_t per_col = ctools_io_chunks_per_column(nvars, nobs, MIN_OBS_PER_THREAD);
    ctools_persistent_pool *chunk_pool = per_col >= 2 ? ctools_get_global_pool() : NULL;
    size_t store_tasks;
    if (chunk_pool && ctools_safe_mul_size(nvars, per_col, &store_tasks) == 0)
    {
        ctools_store_chunk *chunks =
            (ctools_store_chunk *)ctools_safe_malloc2(store_tasks, sizeof(ctools_store_chunk));
        if (chunks)
        {
            atomic_int error = 0;
            size_t built = 0;
            for (size_t j = 0; j < nvars; j++)
            {
                size_t share = nobs / per_col, extra = nobs % per_col;
                size_t dest = 0;
                for (size_t c = 0; c < per_col; c++)
                {
                    size_t count = share + (c < extra);
                    ctools_store_chunk *ck = &chunks[built++];
                    ck->var = &data->vars[j];
                    ck->var_idx = var_indices ? var_indices[j] : (int)j + 1;
                    ck->obs1 = obs1 + dest;
                    ck->dest = dest;
                    ck->count = count;
                    ck->total_nobs = nobs;
                    ck->order = order;
                    ck->error = &error;
                    dest += count;
                }
            }
            if (ctools_persistent_pool_submit_batch(chunk_pool, ctools_io_store_chunk_thread,
                                                    chunks, built, sizeof(ctools_store_chunk)) != 0)
                ctools_io_record_io_error(&error, STATA_ERR_MEMORY);
            else if (ctools_persistent_pool_wait(chunk_pool))
                ctools_io_record_io_error(&error, STATA_ERR_MEMORY);
            free(chunks);
            ctools_memory_barrier();
            return ctools_io_store_status(atomic_load(&error));
        }
    }

    /* Column-parallel store via thread pool */
    ctools_var_io_args *thread_args =
        (ctools_var_io_args *)ctools_safe_malloc2(nvars, sizeof(ctools_var_io_args));
    if (thread_args == NULL)
    {
        return STATA_ERR_MEMORY;
    }

    ctools_io_init_io_thread_args(thread_args, data, var_indices, nvars, obs1, nobs, IO_MODE_STORE,
                                  NULL);
    for (size_t j = 0; j < nvars; j++)
        thread_args[j].store_order = order;
    rc = ctools_io_execute_io_parallel(thread_args, nvars, ctools_io_store_variable_thread);
    free(thread_args);

    ctools_memory_barrier();
    return ctools_io_store_status(rc);
}

stata_retcode ctools_data_store_ex(stata_data *data, int *var_indices, size_t nvars, size_t obs1)
{
    return ctools_io_data_store_impl(data, var_indices, nvars, obs1, NULL);
}

stata_retcode ctools_data_store_sorted(stata_data *data, size_t obs1)
{
    if (!data || (data->nobs && !data->sort_order))
        return STATA_ERR_INVALID_INPUT;
    return ctools_io_data_store_impl(data, NULL, 0, obs1, data->sort_order);
}

stata_retcode ctools_data_store(stata_data *data, size_t obs1)
{
    return ctools_data_store_ex(data, NULL, 0, obs1);
}

stata_retcode ctools_data_store_selective(stata_data *data, int *var_indices, size_t nvars,
                                          size_t obs1)
{
    if (var_indices == NULL || nvars == 0)
    {
        return STATA_ERR_INVALID_INPUT;
    }
    return ctools_data_store_ex(data, var_indices, nvars, obs1);
}

stata_retcode ctools_stream_var_permuted(int var_idx, int64_t *source_rows, size_t output_nobs,
                                         size_t obs1)
{
    stata_retcode rc = ctools_io_validate_variable(var_idx, -1);
    if (rc)
        return rc;
    rc = ctools_io_validate_range(obs1, output_nobs);
    if (rc || (output_nobs && !source_rows))
        return STATA_ERR_INVALID_INPUT;
    if (!output_nobs)
        return STATA_OK;
    size_t available = (size_t)SF_nobs();
    for (size_t i = 0; i < output_nobs; i++)
    {
        if (source_rows[i] < -1 ||
            (source_rows[i] >= 0 && (uint64_t)source_rows[i] > available - obs1))
            return STATA_ERR_INVALID_INPUT;
    }

    if (SF_var_is_string(var_idx))
    {
        rc = ctools_io_validate_variable(var_idx, 1);
        if (rc)
            return rc;
        char **buf = (char **)ctools_safe_cacheline_alloc2(output_nobs, sizeof(char *));
        if (!buf)
            return STATA_ERR_MEMORY;
        memset(buf, 0, output_nobs * sizeof(char *));
        ctools_string_arena *arena =
            output_nobs <= SIZE_MAX / 64
                ? ctools_string_arena_create(output_nobs * 64, CTOOLS_STRING_ARENA_STRDUP_FALLBACK)
                : NULL;
        for (size_t i = 0; i < output_nobs; i++)
        {
            char strbuf[STATA_STR_MAXLEN + 1] = "";
            int io_rc = source_rows[i] >= 0
                            ? SF_sdata(var_idx, (ST_int)(source_rows[i] + obs1), strbuf)
                            : 0;
            if (io_rc)
            {
                rc = ctools_spi_status(io_rc, STATA_ERR_STATA_READ);
                break;
            }
            buf[i] = ctools_string_arena_strdup(arena, strbuf);
            if (!buf[i])
            {
                rc = STATA_ERR_MEMORY;
                break;
            }
        }
        /* Complete the gather before any writes, including in-place permutations. */
        if (!rc)
        {
            for (size_t i = 0; i < output_nobs; i++)
            {
                int io_rc = SF_sstore(var_idx, (ST_int)(obs1 + i), buf[i]);
                if (io_rc)
                {
                    rc = ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
                    break;
                }
            }
        }
        for (size_t i = 0; i < output_nobs; i++)
            if (buf[i] && !ctools_string_arena_owns(arena, buf[i]))
                free(buf[i]);
        ctools_string_arena_free(arena);
        ctools_aligned_free(buf);
    }
    else
    {
        double *buf = (double *)ctools_safe_cacheline_alloc2(output_nobs, sizeof(double));
        if (!buf)
            return STATA_ERR_MEMORY;
        for (size_t i = 0; i < output_nobs; i++)
        {
            buf[i] = SV_missval;
            int io_rc = source_rows[i] >= 0
                            ? IO_VDATA_FN(var_idx, (ST_int)(source_rows[i] + obs1), &buf[i])
                            : 0;
            if (io_rc)
            {
                rc = ctools_spi_status(io_rc, STATA_ERR_STATA_READ);
                break;
            }
        }
        if (!rc)
        {
            for (size_t i = 0; i < output_nobs; i++)
            {
                int io_rc = IO_VSTORE_FN(var_idx, (ST_int)(obs1 + i), buf[i]);
                if (io_rc)
                {
                    rc = ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
                    break;
                }
            }
        }
        ctools_aligned_free(buf);
    }
    return rc == STATA_ERR_STATA_WRITE ? ctools_io_store_status(rc) : rc;
}

stata_retcode ctools_store_filtered(double *values, size_t n_filtered, int var_idx,
                                    perm_idx_t *obs_map)
{
    if (!n_filtered)
        return STATA_OK;
    if (!values)
        return STATA_ERR_INVALID_INPUT;
    stata_retcode rc = ctools_io_validate_variable(var_idx, 0);
    if (rc)
        return rc;
    rc = ctools_io_validate_obs_map(obs_map, n_filtered);
    if (rc)
        return rc;
    for (size_t i = 0; i < n_filtered; i++)
    {
        {
            int io_rc = IO_VSTORE_FN(var_idx, (ST_int)obs_map[i], values[i]);
            if (io_rc)
                return ctools_io_store_status(ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE));
        }
    }
    return STATA_OK;
}

stata_retcode ctools_store_filtered_rowpar(double *values, size_t n_filtered, int var_idx,
                                           perm_idx_t *obs_map)
{
    if (!n_filtered)
        return STATA_OK;
    if (!values)
        return STATA_ERR_INVALID_INPUT;
    stata_retcode rc = ctools_io_validate_variable(var_idx, 0);
    if (rc)
        return rc;
    int parallel = 0, distinct = 0;
#ifdef _OPENMP
    parallel = n_filtered >= MIN_OBS_PER_THREAD * 2 && ctools_get_max_threads() > 1;
#endif
    if (!parallel)
    {
        /* Serial writes already preserve repeated destinations; only bounds
         * need checking for small transfers and builds without OpenMP. */
        rc = ctools_io_validate_obs_map(obs_map, n_filtered);
        if (rc)
            return rc;
    }
    else
    {
        if (!obs_map)
            return STATA_ERR_INVALID_INPUT;
        const perm_idx_t available = (perm_idx_t)SF_nobs();
        /* Independent OR reductions let the compiler vectorize both proofs.
         * Unsigned subtraction rejects zero as well as values above nobs.
         * Examine every destination before permitting any worker writes. */
        unsigned invalid = (obs_map[0] - 1u >= available);
        unsigned nonincreasing = 0;
        for (size_t i = 1; i < n_filtered; i++)
        {
            const perm_idx_t obs = obs_map[i];
            invalid |= (obs - 1u >= available);
            nonincreasing |= (obs <= obs_map[i - 1]);
        }
        if (invalid)
            return STATA_ERR_INVALID_INPUT;
        distinct = !nonincreasing;
        if (!distinct)
        {
            size_t min_obs = (size_t)obs_map[0], max_obs = min_obs;
            int decreasing = 1;
            for (size_t i = 1; i < n_filtered; i++)
            {
                size_t obs = (size_t)obs_map[i];
                if (obs_map[i] >= obs_map[i - 1])
                    decreasing = 0;
                if (obs < min_obs)
                    min_obs = obs;
                if (obs > max_obs)
                    max_obs = obs;
            }
            distinct = decreasing;
            /* By-group commands supply unordered unique permutations.
             * Prove uniqueness with scratch no larger than the existing
             * map; sparse ranges or OOM safely retain sequential writes. */
            size_t bytes = (max_obs - min_obs) / 8 + 1;
            size_t map_cells = bytes / sizeof(perm_idx_t) + (bytes % sizeof(perm_idx_t) != 0);
            if (!distinct && map_cells <= n_filtered)
            {
                unsigned char *seen = (unsigned char *)calloc(bytes, 1);
                if (seen)
                {
                    distinct = 1;
                    for (size_t i = 0; i < n_filtered; i++)
                    {
                        size_t bit = (size_t)obs_map[i] - min_obs;
                        unsigned char mask = (unsigned char)(1u << (bit & 7));
                        if (seen[bit >> 3] & mask)
                        {
                            distinct = 0;
                            break;
                        }
                        seen[bit >> 3] |= mask;
                    }
                    free(seen);
                }
            }
        }
    }
    ST_IIID store_fn = IO_VSTORE_FN;
    if (!parallel || !distinct)
    {
        /* Repeated destinations require input-order last-write-wins. */
        for (size_t i = 0; i < n_filtered; i++)
        {
            int io_rc = store_fn(var_idx, (ST_int)obs_map[i], values[i]);
            if (io_rc)
                return ctools_io_store_status(ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE));
        }
        return STATA_OK;
    }
    /* The existing OpenMP callback contract is unchanged: checked numeric
     * SPI writes to disjoint destinations. The map proof establishes that
     * disjointness before any worker may modify Stata's data. */
    atomic_int error = 0;
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_filtered; i++)
    {
        int io_rc = store_fn(var_idx, (ST_int)obs_map[i], values[i]);
        if (io_rc)
            ctools_io_record_io_error(&error, ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE));
    }
    return ctools_io_store_status(atomic_load(&error));
}

stata_retcode ctools_store_filtered_cols(const double *const *columns, const double *scale,
                                         const double *fill, const int *var_idx, size_t ncols,
                                         size_t n_filtered, const perm_idx_t *obs_map,
                                         size_t fill_rows, int *unstorable)
{
    if (unstorable)
        *unstorable = 0;
    if (!ncols)
        return STATA_OK;
    if (!columns || !var_idx || (n_filtered && !obs_map))
        return STATA_ERR_INVALID_INPUT;
    for (size_t j = 0; j < ncols; j++)
    {
        stata_retcode rc = ctools_io_validate_variable(var_idx[j], 0);
        if (rc)
            return rc;
    }
    const size_t available = (size_t)SF_nobs();
    if (fill_rows > available)
        return STATA_ERR_INVALID_INPUT;
    /* Strictly increasing destinations from 1 to the last observation (or to
     * fill_rows) are distinct and in range: checked before any write, with
     * an OR reduction the compiler can vectorize. */
    if (n_filtered)
    {
        const size_t last = fill_rows ? fill_rows : available;
        unsigned bad = (obs_map[0] < 1) | ((size_t)obs_map[n_filtered - 1] > last);
        for (size_t i = 1; i < n_filtered; i++)
            bad |= (obs_map[i] <= obs_map[i - 1]);
        if (bad)
            return STATA_ERR_INVALID_INPUT;
    }
    const size_t rows = fill_rows ? fill_rows : n_filtered;
    if (!rows)
        return STATA_OK;
    const double missval = SV_missval;
    ST_IIID store_fn = IO_VSTORE_FN;
    int threads = 1;
#ifdef _OPENMP
    if (rows >= MIN_OBS_PER_THREAD * 2)
        threads = ctools_get_max_threads();
#endif
    if (threads < 1)
        threads = 1;
    atomic_int error = 0;
    int bad_value = 0;
    /* Each worker owns a contiguous range of rows (filtered indices, or with
     * fill_rows the observations themselves) and writes every column of a
     * tile of rows before the next tile: each record is fetched once while
     * all its outputs are written, not once per column. */
#pragma omp parallel num_threads(threads) if (threads > 1) reduction(| : bad_value)
    {
        size_t team = 1, tid = 0;
#ifdef _OPENMP
        team = (size_t)omp_get_num_threads();
        tid = (size_t)omp_get_thread_num();
#endif
        const size_t base = rows / team, extra = rows % team;
        const size_t start = tid * base + (tid < extra ? tid : extra);
        const size_t stop = start + base + (tid < extra);
        size_t q = start;
        if (fill_rows)
        {
            /* the first filtered row at or after observation start + 1 */
            size_t lo = 0, hi = n_filtered;
            while (lo < hi)
            {
                size_t mid = lo + (hi - lo) / 2;
                if ((size_t)obs_map[mid] < start + 1)
                    lo = mid + 1;
                else
                    hi = mid;
            }
            q = lo;
        }
        for (size_t tile = start; tile < stop && !atomic_load(&error); tile += IO_STORE_TILE_ROWS)
        {
            const size_t tile_end = stop - tile > IO_STORE_TILE_ROWS ? tile + IO_STORE_TILE_ROWS : stop;
            size_t next = q;
            for (size_t j = 0; j < ncols; j++)
            {
                const double *col = columns[j];
                const double s = scale ? scale[j] : 1.0;
                const double other = fill ? fill[j] : missval;
                const int var = var_idx[j];
                size_t p = q;
                int io_rc = 0;
                if (fill_rows)
                {
                    for (size_t o = tile + 1; o <= tile_end && !io_rc; o++)
                    {
                        double v = other;
                        if (p < n_filtered && (size_t)obs_map[p] == o)
                        {
                            v = col ? col[p] * s : s;
                            if (!(v < missval && v > -INFINITY))
                            {
                                bad_value = 1;
                                v = missval;
                            }
                            p++;
                        }
                        io_rc = store_fn(var, (ST_int)o, v);
                    }
                }
                else
                {
                    for (; p < tile_end && !io_rc; p++)
                    {
                        double v = col ? col[p] * s : s;
                        if (!(v < missval && v > -INFINITY))
                        {
                            bad_value = 1;
                            v = missval;
                        }
                        io_rc = store_fn(var, (ST_int)obs_map[p], v);
                    }
                }
                if (io_rc)
                {
                    ctools_io_record_io_error(&error, ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE));
                    break;
                }
                next = p;
            }
            q = next;
        }
    }
    ctools_memory_barrier();
    if (unstorable)
        *unstorable = bad_value;
    return ctools_io_store_status(atomic_load(&error));
}

stata_retcode ctools_store_filtered_str(char **strings, size_t n_filtered, int var_idx,
                                        perm_idx_t *obs_map)
{
    if (!n_filtered)
        return STATA_OK;
    if (!strings)
        return STATA_ERR_INVALID_INPUT;
    stata_retcode rc = ctools_io_validate_variable(var_idx, 1);
    if (rc)
        return rc;
    /* Duplicate destinations require sequential last-write-wins behavior.
     * Small/no-OpenMP transfers retain the existing validation and store loop. */
    int ordered = 0;
#ifdef _OPENMP
    if (n_filtered >= MIN_OBS_PER_THREAD * 2 && ctools_get_max_threads() > 1)
    {
        if (!obs_map)
            return STATA_ERR_INVALID_INPUT;
        size_t available = (size_t)SF_nobs();
        ordered = 1;
        for (size_t i = 0; i < n_filtered; i++)
        {
            if (obs_map[i] < 1 || obs_map[i] > available)
                return STATA_ERR_INVALID_INPUT;
            if (i && obs_map[i] <= obs_map[i - 1])
                ordered = 0;
        }
    }
    else
#endif
    {
        rc = ctools_io_validate_obs_map(obs_map, n_filtered);
        if (rc)
            return rc;
    }
    ST_IIIS store_fn = (_stata_)->sstore;
    if (ordered && n_filtered >= MIN_OBS_PER_THREAD * 2)
    {
        atomic_int error = 0;
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < n_filtered; i++)
        {
            int io_rc = store_fn(var_idx, (ST_int)obs_map[i], strings[i] ? strings[i] : "");
            if (io_rc)
                ctools_io_record_io_error(&error, ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE));
        }
        return ctools_io_store_status(atomic_load(&error));
    }
    for (size_t i = 0; i < n_filtered; i++)
    {
        int io_rc = store_fn(var_idx, (ST_int)obs_map[i], strings[i] ? strings[i] : "");
        if (io_rc)
            return ctools_io_store_status(ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE));
    }
    return STATA_OK;
}

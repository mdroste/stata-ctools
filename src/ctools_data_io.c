#include "ctools_transport_internal.h"

int ctools_io_load_single_variable(stata_variable *var, int var_idx, size_t obs1,
                                 size_t nobs, int is_string, int str_width)
{
    size_t i;

    var->nobs = nobs;

    /* Edge case: if no observations, allocate minimal buffer to avoid NULL pointer issues.
     * This prevents crashes when downstream code assumes data.dbl/data.str is non-NULL
     * after successful load. Allocating 1 element is minimal overhead for safety. */
    if (nobs == 0) {
        if (is_string) {
            var->type = STATA_TYPE_STRING;
            var->str_maxlen = 0;
            var->_arena = NULL;
            /* Preserve the aligned allocation contract even for an empty column. */
            var->data.str = (char **)ctools_cacheline_alloc(sizeof(char *));
            if (var->data.str) var->data.str[0] = NULL;
            if (!var->data.str) return STATA_ERR_MEMORY;
        } else {
            var->type = STATA_TYPE_DOUBLE;
            var->_arena = NULL;
            /* Allocate minimal buffer to avoid NULL */
            var->data.dbl = (double *)ctools_cacheline_alloc(sizeof(double));
            if (!var->data.dbl) return STATA_ERR_MEMORY;
            var->data.dbl[0] = SV_missval;  /* Initialize to missing */
        }
        return 0;  /* Success with empty data */
    }

    if (is_string) {
        /* Storage type is constant for a column throughout the transfer. */
        ST_IIIS read_fn = is_string == IO_LONG_STRING
            ? ctools_io_read_stata_strl : (_stata_)->sdata;
        var->type = STATA_TYPE_STRING;

        /* ---- Flat buffer path: known string width ----
           When the actual Stata variable width is known (e.g. 17 for str17),
           allocate a single contiguous nobs×stride buffer. Read through a
           bounded temporary buffer: ado width hints can be stale. Benefits:
           - No per-string arena reservation
           - Fixed strides for strings
           - Better cache locality for subsequent sort comparisons
           - O(1) cleanup via ctools_string_arena_free wrapper */
        if (str_width > 0 && str_width < STATA_STR_MAXLEN) {
            size_t stride = (size_t)str_width + 1;
            size_t flat_size = nobs * stride;

            /* Overflow and size guard (2 GB cap) */
            if (flat_size / stride == nobs && flat_size <= (2048ULL * 1024 * 1024)) {
                /* Whole guard slots keep the (count, stride) allocation shape. */
                size_t guard_slots = (IO_SLOT_GUARD + stride - 1) / stride;
                char *flat_buf = (char *)calloc(nobs + guard_slots, stride);
                if (flat_buf) {
                    /* Allocate pointer array */
                    size_t str_array_size;
                    if (ctools_safe_mul_size(nobs, sizeof(char *), &str_array_size) != 0) {
                        free(flat_buf);
                        return STATA_ERR_MEMORY;
                    }
                    char **str_ptrs = (char **)ctools_cacheline_alloc(str_array_size);
                    if (!str_ptrs) {
                        free(flat_buf);
                        return STATA_ERR_MEMORY;
                    }

                    /* One worker owns the whole column, so every slot may be
                     * read directly; the guard tail bounds the final rows. */
                    const char *flat_safe_end = flat_buf + (nobs + guard_slots) * stride;
                    for (i = 0; i < nobs; i++) {
                        str_ptrs[i] = flat_buf + i * stride;
                        int read_rc = ctools_io_read_string_slot_direct(read_fn, var_idx, (ST_int)(i + obs1), str_ptrs[i], (size_t)str_width, flat_safe_end);
                        if (read_rc) {
                            ctools_aligned_free(str_ptrs);
                            free(flat_buf);
                            return read_rc;
                        }
                    }

                    /* Wrap flat buffer as a ctools_string_arena for compatible cleanup.
                       stata_data_free checks has_fallback=0 → takes O(1) bulk free path:
                       ctools_string_arena_free does free(base) + free(arena). */
                    ctools_string_arena *arena = (ctools_string_arena *)malloc(sizeof(ctools_string_arena));
                    if (!arena) {
                        ctools_aligned_free(str_ptrs);
                        free(flat_buf);
                        return STATA_ERR_MEMORY;
                    }
                    arena->base = flat_buf;
                    arena->capacity = flat_size;
                    arena->used = flat_size;
                    arena->mode = CTOOLS_STRING_ARENA_NO_FALLBACK;
                    arena->has_fallback = 0;

                    var->str_maxlen = (size_t)str_width;
                    var->data.str = str_ptrs;
                    var->_arena = arena;
                    return 0;
                }
                /* calloc failed — fall through to arena path */
            }
        }

        /* ---- Arena path: width unknown or flat buffer not feasible ---- */
        var->str_maxlen = STATA_STR_MAXLEN;
        var->_arena = NULL;

        /* Allocate pointer array (cache-line aligned, zero-initialized for safe cleanup) */
        size_t str_array_size;
        if (ctools_safe_mul_size(nobs, sizeof(char *), &str_array_size) != 0) {
            return STATA_ERR_MEMORY;  /* Overflow */
        }
        var->data.str = (char **)ctools_cacheline_alloc(str_array_size);
        if (var->data.str == NULL) {
            return STATA_ERR_MEMORY;
        }
        memset(var->data.str, 0, str_array_size);

        char **str_ptr = var->data.str;
        char strbuf[STATA_STR_MAXLEN + 1];

        /* Create arena for all strings - estimate avg 64 bytes per string */
        /* Overflow check: nobs * 64 - skip arena if overflow, rely on strdup fallback */
        ctools_string_arena *arena = NULL;
        if (nobs <= SIZE_MAX / 64) {
            arena = ctools_io_create_private_string_arena(nobs, str_width);
        }
        if (arena != NULL) {
            var->_arena = arena;
        }

        /* Load strings */
        ST_IIIS sdata_fn2 = read_fn;
        for (i = 0; i < nobs; i++) {
            int read_rc = sdata_fn2((ST_int)var_idx, (ST_int)(i + obs1), strbuf);
            if (read_rc) return read_fn == ctools_io_read_stata_strl ? read_rc : ctools_spi_status(read_rc, STATA_ERR_STATA_READ);
            str_ptr[i] = ctools_io_copy_private_string(arena, strbuf);
            if (str_ptr[i] == NULL) {
                /* Cleanup on allocation failure:
                 * - Free fallback strings (not owned by arena)
                 * - Free the string pointer array
                 * - Free the arena itself */
                if (arena != NULL) {
                    for (size_t j = 0; j < i; j++) {
                        if (str_ptr[j] != NULL && !ctools_string_arena_owns(arena, str_ptr[j])) {
                            free(str_ptr[j]);
                        }
                    }
                    ctools_string_arena_free(arena);
                    var->_arena = NULL;
                } else {
                    /* No arena - all strings were strdup'd */
                    for (size_t j = 0; j < i; j++) {
                        free(str_ptr[j]);
                    }
                }
                ctools_aligned_free(var->data.str);
                var->data.str = NULL;
                return STATA_ERR_MEMORY;
            }
        }
    } else {
        /* Numeric variable - use cache-line aligned allocation */
        var->type = STATA_TYPE_DOUBLE;
        var->_arena = NULL;
        size_t dbl_array_size;
        if (ctools_safe_mul_size(nobs, sizeof(double), &dbl_array_size) != 0) {
            return STATA_ERR_MEMORY;  /* Overflow */
        }
        var->data.dbl = (double *)ctools_cacheline_alloc(dbl_array_size);

        if (var->data.dbl == NULL) {
            return STATA_ERR_MEMORY;
        }

        double * restrict dbl_ptr = var->data.dbl;

        /* Cache SPI function pointer — avoid reloading _stata_ per iteration */
        ST_IIIDp vdata_fn = IO_VDATA_FN;
        for (i = 0; i + 4 <= nobs; i += 4) {
            { int io_rc = vdata_fn((ST_int)var_idx, (ST_int)(i + 0 + obs1), &dbl_ptr[i + 0]); if (io_rc) return ctools_spi_status(io_rc, STATA_ERR_STATA_READ); }
            { int io_rc = vdata_fn((ST_int)var_idx, (ST_int)(i + 1 + obs1), &dbl_ptr[i + 1]); if (io_rc) return ctools_spi_status(io_rc, STATA_ERR_STATA_READ); }
            { int io_rc = vdata_fn((ST_int)var_idx, (ST_int)(i + 2 + obs1), &dbl_ptr[i + 2]); if (io_rc) return ctools_spi_status(io_rc, STATA_ERR_STATA_READ); }
            { int io_rc = vdata_fn((ST_int)var_idx, (ST_int)(i + 3 + obs1), &dbl_ptr[i + 3]); if (io_rc) return ctools_spi_status(io_rc, STATA_ERR_STATA_READ); }
        }
        for (; i < nobs; i++) {
            { int io_rc = vdata_fn((ST_int)var_idx, (ST_int)(i + 0 + obs1), &dbl_ptr[i + 0]); if (io_rc) return ctools_spi_status(io_rc, STATA_ERR_STATA_READ); }
        }
    }

    return 0;
}

void *ctools_io_load_variable_thread(void *arg)
{
    ctools_var_io_args *args = (ctools_var_io_args *)arg;
    args->error = ctools_io_load_single_variable(args->var, args->var_idx, args->obs1,
                                       args->nobs, args->is_string, args->str_width);
    return args->error ? (void *)1 : NULL;
}

void *ctools_io_load_chunk_thread(void *arg)
{
    ctools_io_chunk *a = (ctools_io_chunk *)arg;
    if (atomic_load(a->error)) return NULL;
    double *dst = a->var->data.dbl + a->dest;
    ST_IIIDp fn = IO_VDATA_FN;
    for (size_t i = 0; i < a->count; i++) {
        ST_int obs = a->map ? (ST_int)a->map[a->dest + i] : (ST_int)(a->obs1 + i);
        int io_rc = fn(a->var_idx, obs, dst + i);
        if (io_rc) {
            ctools_io_record_io_error(a->error, ctools_spi_status(io_rc, STATA_ERR_STATA_READ));
            return (void *)1;
        }
    }
    return NULL;
}

int ctools_io_try_chunked_load(stata_data *data, int *var_indices, size_t nvars,
                            const perm_idx_t *obs_map, size_t nobs,
                            int was_filtered, stata_retcode *rc_out)
{
    *rc_out = STATA_OK;
    for (size_t j = 0; j < nvars; j++)
        if (SF_var_is_string(var_indices[j])) return 0;
    size_t per_col = ctools_io_chunks_per_column(nvars, nobs, IO_CHUNK_MIN_NUMERIC_ROWS);
    if (nvars < 2 || per_col < 2) return 0;
    ctools_persistent_pool *pool = ctools_get_global_pool();
    if (!pool) return 0;

    size_t total_tasks;
    if (ctools_safe_mul_size(nvars, per_col, &total_tasks) != 0) return 0;
    ctools_io_chunk *chunks = (ctools_io_chunk *)
        ctools_safe_malloc2(total_tasks, sizeof(ctools_io_chunk));
    if (!chunks) { *rc_out = STATA_ERR_MEMORY; return 1; }

    atomic_int error = 0;
    size_t built = 0;
    for (size_t j = 0; j < nvars; j++) {
        stata_variable *var = &data->vars[j];
        var->nobs = nobs;
        var->type = STATA_TYPE_DOUBLE;
        var->_arena = NULL;
        size_t bytes;
        if (ctools_safe_mul_size(nobs, sizeof(double), &bytes) != 0 ||
            !(var->data.dbl = (double *)ctools_cacheline_alloc(bytes))) {
            free(chunks); *rc_out = STATA_ERR_MEMORY; return 1;
        }
        size_t share = nobs / per_col, extra = nobs % per_col;
        size_t dest = 0;
        for (size_t c = 0; c < per_col; c++) {
            size_t count = share + (c < extra);
            ctools_io_chunk *ck = &chunks[built++];
            ck->var = var;
            ck->var_idx = var_indices[j];
            ck->map = was_filtered ? obs_map : NULL;
            ck->obs1 = was_filtered ? 0 : (size_t)obs_map[0] + dest;
            ck->dest = dest;
            ck->count = count;
            ck->error = &error;
            dest += count;
        }
    }

    /* A failed batch submit may already have queued chunks; wait() inside
     * submit_batch drains them, so treat it as the allocation failure it is. */
    if (ctools_persistent_pool_submit_batch(pool, ctools_io_load_chunk_thread,
                chunks, built, sizeof(ctools_io_chunk)) != 0)
        ctools_io_record_io_error(&error, STATA_ERR_MEMORY);
    else if (ctools_persistent_pool_wait(pool))
        ctools_io_record_io_error(&error, STATA_ERR_MEMORY);
    free(chunks);
    ctools_memory_barrier();
    *rc_out = (stata_retcode)atomic_load(&error);
    return 1;
}

int ctools_io_load_row_parallel(stata_variable *var, int var_idx,
                              perm_idx_t *obs_map, size_t n_filtered,
                              int str_width)
{
    int is_string = ctools_io_variable_kind(var_idx);

    if (is_string == IO_LONG_STRING)
        return ctools_io_load_filtered_variable(var, var_idx, obs_map, n_filtered,
                                      is_string, str_width);

    var->nobs = n_filtered;

    /* Edge case: empty filtered set */
    if (n_filtered == 0) {
        if (is_string) {
            var->type = STATA_TYPE_STRING;
            var->str_maxlen = 0;
            var->_arena = NULL;
            var->data.str = (char **)ctools_cacheline_alloc(sizeof(char *));
            if (var->data.str) var->data.str[0] = NULL;
            return var->data.str ? 0 : -1;
        } else {
            var->type = STATA_TYPE_DOUBLE;
            var->_arena = NULL;
            var->data.dbl = (double *)ctools_cacheline_alloc(sizeof(double));
            if (!var->data.dbl) return STATA_ERR_MEMORY;
            var->data.dbl[0] = SV_missval;
            return 0;
        }
    }

    if (is_string) {
        /* Storage type is constant for a column throughout the transfer. */
        ST_IIIS read_fn = (_stata_)->sdata;
        var->type = STATA_TYPE_STRING;
        var->str_maxlen = STATA_STR_MAXLEN;
        var->_arena = NULL;

        size_t str_array_size;
        if (ctools_safe_mul_size(n_filtered, sizeof(char *), &str_array_size) != 0) {
            return STATA_ERR_MEMORY;
        }
        var->data.str = (char **)ctools_cacheline_alloc(str_array_size);
        if (var->data.str == NULL) return STATA_ERR_MEMORY;
        memset(var->data.str, 0, str_array_size);

        /* Fast path: flat buffer when string width is known */
        int used_flat_buffer = 0;
        if (str_width > 0 && str_width < STATA_STR_MAXLEN) {
            size_t stride = (size_t)str_width + 1;
            size_t flat_size;
            if (ctools_safe_mul_size(n_filtered, stride, &flat_size) == 0) {
                size_t guard_slots = (IO_SLOT_GUARD + stride - 1) / stride;
                char *flat_buf = n_filtered <= SIZE_MAX / stride - guard_slots
                    ? (char *)calloc(n_filtered + guard_slots, stride) : NULL;
                if (flat_buf != NULL) {
                    atomic_int load_error = 0;
                    #pragma omp parallel if(n_filtered >= MIN_OBS_PER_THREAD * 2)
                    {
                        size_t team = 1, tid = 0;
                        #ifdef _OPENMP
                        team = (size_t)omp_get_num_threads();
                        tid = (size_t)omp_get_thread_num();
                        #endif
                        size_t share = n_filtered / team, extra = n_filtered % team;
                        size_t start = tid * share + (tid < extra ? tid : extra);
                        size_t count = share + (tid < extra);
                        /* Workers own disjoint row ranges; only the final
                         * range extends into the guard tail. */
                        const char *region_end = flat_buf + (start + count) * stride
                            + (start + count == n_filtered ? guard_slots * stride : 0);
                        for (size_t i = start; i < start + count; i++) {
                            if (atomic_load(&load_error)) break;
                            char *slot = flat_buf + i * stride;
                            int read_rc = ctools_io_read_string_slot_direct(read_fn, var_idx, (ST_int)obs_map[i], slot, (size_t)str_width, region_end);
                            ctools_io_record_io_error(&load_error, read_rc);
                            var->data.str[i] = slot;
                        }
                    }

                    if (atomic_load(&load_error)) {
                        free(flat_buf);
                        ctools_aligned_free(var->data.str);
                        var->data.str = NULL;
                        return atomic_load(&load_error);
                    }

                    /* Wrap flat buffer as arena for compatible cleanup */
                    ctools_string_arena *wrapper = (ctools_string_arena *)calloc(1, sizeof(ctools_string_arena));
                    if (wrapper != NULL) {
                        wrapper->base = flat_buf;
                        wrapper->capacity = flat_size;
                        wrapper->used = flat_size;
                        wrapper->has_fallback = 0;
                        var->_arena = wrapper;
                        used_flat_buffer = 1;
                    } else {
                        free(flat_buf);
                        ctools_aligned_free(var->data.str);
                        var->data.str = NULL;
                        return STATA_ERR_MEMORY;
                    }
                }
            }
        }

        if (!used_flat_buffer) {
            /* Fallback: arena with per-thread slabs */
            ctools_string_arena *arena = NULL;
            if (n_filtered <= SIZE_MAX / 64) {
                arena = ctools_io_create_private_string_arena(n_filtered, str_width);
            }
            if (arena != NULL) {
                var->_arena = arena;
            }

            atomic_int load_error = 0;

            #pragma omp parallel if(n_filtered >= MIN_OBS_PER_THREAD * 2)
            {
                /* Partition arena into per-thread slabs */
                char *slab_ptr = NULL;
                char *slab_end = NULL;
                if (arena != NULL) {
                    int nthreads = 1;
                    int tid = 0;
                    #ifdef _OPENMP
                    nthreads = omp_get_num_threads();
                    tid = omp_get_thread_num();
                    #endif
                    size_t slab_size = arena->capacity / (size_t)nthreads;
                    slab_ptr = arena->base + (size_t)tid * slab_size;
                    slab_end = (tid == nthreads - 1)
                        ? arena->base + arena->capacity
                        : slab_ptr + slab_size;
                }

                ST_IIIS sdata_fn2 = read_fn;
                #pragma omp for schedule(static)
                for (size_t i = 0; i < n_filtered; i++) {
                    if (atomic_load(&load_error)) continue;
                    char strbuf[STATA_STR_MAXLEN + 1];
                    int read_rc = sdata_fn2((ST_int)var_idx, (ST_int)obs_map[i], strbuf);
                    if (read_rc) {
                        ctools_io_record_io_error(&load_error, ctools_spi_status(read_rc, STATA_ERR_STATA_READ));
                        continue;
                    }
                    size_t len = strlen(strbuf) + 1;

                    char *s;
                    if (slab_ptr != NULL && slab_ptr + len <= slab_end) {
                        memcpy(slab_ptr, strbuf, len);
                        s = slab_ptr;
                        slab_ptr += len;
                    } else {
                        s = strdup(strbuf);
                        if (arena != NULL) ARENA_ATOMIC_STORE_INT(&arena->has_fallback, 1);
                        if (s == NULL) {
                            ctools_io_record_io_error(&load_error, STATA_ERR_MEMORY);
                            continue;
                        }
                    }
                    var->data.str[i] = s;
                }

                #pragma omp single
                {
                    if (arena != NULL) arena->used = arena->capacity;
                }
            }

            if (atomic_load(&load_error)) {
                /* Cleanup: free fallback strings, arena, pointer array */
                if (arena != NULL) {
                    for (size_t j = 0; j < n_filtered; j++) {
                        if (var->data.str[j] != NULL &&
                            !ctools_string_arena_owns(arena, var->data.str[j])) {
                            free(var->data.str[j]);
                        }
                    }
                    ctools_string_arena_free(arena);
                    var->_arena = NULL;
                } else {
                    for (size_t j = 0; j < n_filtered; j++) {
                        free(var->data.str[j]);
                    }
                }
                ctools_aligned_free(var->data.str);
                var->data.str = NULL;
                return atomic_load(&load_error);
            }
        }
    } else {
        /* Numeric variable */
        var->type = STATA_TYPE_DOUBLE;
        var->_arena = NULL;

        size_t dbl_array_size;
        if (ctools_safe_mul_size(n_filtered, sizeof(double), &dbl_array_size) != 0) {
            return STATA_ERR_MEMORY;
        }
        var->data.dbl = (double *)ctools_cacheline_alloc(dbl_array_size);
        if (var->data.dbl == NULL) return STATA_ERR_MEMORY;

        double * restrict dbl_ptr = var->data.dbl;

        /* Cache SPI function pointer — avoid reloading _stata_ per iteration */
        ST_IIIDp vdata_fn = IO_VDATA_FN;
        atomic_int load_error = 0;
        #pragma omp parallel for schedule(static) if(n_filtered >= MIN_OBS_PER_THREAD * 2)
        for (size_t i = 0; i < n_filtered; i++) {
            { int io_rc = vdata_fn((ST_int)var_idx, (ST_int)obs_map[i], &dbl_ptr[i]); if (io_rc) ctools_io_record_io_error(&load_error, ctools_spi_status(io_rc, STATA_ERR_STATA_READ)); }
        }
        return atomic_load(&load_error);
    }

    return 0;
}

void ctools_filtered_data_init(ctools_filtered_data *fd)
{
    if (fd == NULL) return;
    stata_data_init(&fd->data);
    fd->obs_map = NULL;
    fd->n_range = 0;
    fd->was_filtered = 0;
}

void ctools_filtered_data_free(ctools_filtered_data *fd)
{
    if (fd == NULL) return;
    stata_data_free(&fd->data);
    if (fd->obs_map != NULL) {
        ctools_aligned_free(fd->obs_map);
        fd->obs_map = NULL;
    }
    fd->n_range = 0;
    fd->was_filtered = 0;
}

int ctools_io_load_filtered_variable(stata_variable *var, int var_idx,
                                   perm_idx_t *obs_map, size_t n_filtered,
                                   int is_string, int str_width)
{
    size_t i;

    var->nobs = n_filtered;

    /* Edge case: empty filtered set */
    if (n_filtered == 0) {
        if (is_string) {
            var->type = STATA_TYPE_STRING;
            var->str_maxlen = 0;
            var->_arena = NULL;
            var->data.str = (char **)ctools_cacheline_alloc(sizeof(char *));
            if (var->data.str) var->data.str[0] = NULL;
            if (!var->data.str) return STATA_ERR_MEMORY;
        } else {
            var->type = STATA_TYPE_DOUBLE;
            var->_arena = NULL;
            var->data.dbl = (double *)ctools_cacheline_alloc(sizeof(double));
            if (!var->data.dbl) return STATA_ERR_MEMORY;
            var->data.dbl[0] = SV_missval;
        }
        return 0;
    }

    if (is_string) {
        /* Storage type is constant for a column throughout the transfer. */
        ST_IIIS read_fn = is_string == IO_LONG_STRING
            ? ctools_io_read_stata_strl : (_stata_)->sdata;
        var->type = STATA_TYPE_STRING;

        /* ---- Flat buffer path: known string width ---- */
        if (str_width > 0 && str_width < STATA_STR_MAXLEN) {
            size_t stride = (size_t)str_width + 1;
            size_t flat_size = n_filtered * stride;

            if (flat_size / stride == n_filtered && flat_size <= (2048ULL * 1024 * 1024)) {
                size_t guard_slots = (IO_SLOT_GUARD + stride - 1) / stride;
                char *flat_buf = (char *)calloc(n_filtered + guard_slots, stride);
                if (flat_buf) {
                    size_t str_array_size;
                    if (ctools_safe_mul_size(n_filtered, sizeof(char *), &str_array_size) != 0) {
                        free(flat_buf);
                        return STATA_ERR_MEMORY;
                    }
                    char **str_ptrs = (char **)ctools_cacheline_alloc(str_array_size);
                    if (!str_ptrs) {
                        free(flat_buf);
                        return STATA_ERR_MEMORY;
                    }

                    const char *flat_safe_end = flat_buf + (n_filtered + guard_slots) * stride;
                    for (i = 0; i < n_filtered; i++) {
                        str_ptrs[i] = flat_buf + i * stride;
                        int read_rc = ctools_io_read_string_slot_direct(read_fn, var_idx, (ST_int)obs_map[i], str_ptrs[i], (size_t)str_width, flat_safe_end);
                        if (read_rc) {
                            ctools_aligned_free(str_ptrs);
                            free(flat_buf);
                            return read_rc;
                        }
                    }

                    ctools_string_arena *arena = (ctools_string_arena *)malloc(sizeof(ctools_string_arena));
                    if (!arena) {
                        ctools_aligned_free(str_ptrs);
                        free(flat_buf);
                        return STATA_ERR_MEMORY;
                    }
                    arena->base = flat_buf;
                    arena->capacity = flat_size;
                    arena->used = flat_size;
                    arena->mode = CTOOLS_STRING_ARENA_NO_FALLBACK;
                    arena->has_fallback = 0;

                    var->str_maxlen = (size_t)str_width;
                    var->data.str = str_ptrs;
                    var->_arena = arena;
                    return 0;
                }
            }
        }

        /* ---- Arena path: width unknown or flat buffer not feasible ---- */
        var->str_maxlen = STATA_STR_MAXLEN;
        var->_arena = NULL;

        /* Allocate pointer array */
        size_t str_array_size;
        if (ctools_safe_mul_size(n_filtered, sizeof(char *), &str_array_size) != 0) {
            return STATA_ERR_MEMORY;
        }
        var->data.str = (char **)ctools_cacheline_alloc(str_array_size);
        if (var->data.str == NULL) {
            return STATA_ERR_MEMORY;
        }
        memset(var->data.str, 0, str_array_size);

        char **str_ptr = var->data.str;
        char strbuf[STATA_STR_MAXLEN + 1];

        /* Create arena for string storage */
        ctools_string_arena *arena = NULL;
        if (n_filtered <= SIZE_MAX / 64) {
            arena = ctools_io_create_private_string_arena(n_filtered, str_width);
        }
        if (arena != NULL) {
            var->_arena = arena;
        }

        /* Load strings using obs_map */
        ST_IIIS sdata_fn3 = read_fn;
        for (i = 0; i < n_filtered; i++) {
            int read_rc = sdata_fn3((ST_int)var_idx, (ST_int)obs_map[i], strbuf);
            if (read_rc) return read_fn == ctools_io_read_stata_strl ? read_rc : ctools_spi_status(read_rc, STATA_ERR_STATA_READ);
            str_ptr[i] = ctools_io_copy_private_string(arena, strbuf);
            if (str_ptr[i] == NULL) {
                /* Cleanup on failure */
                if (arena != NULL) {
                    for (size_t j = 0; j < i; j++) {
                        if (str_ptr[j] != NULL && !ctools_string_arena_owns(arena, str_ptr[j])) {
                            free(str_ptr[j]);
                        }
                    }
                    ctools_string_arena_free(arena);
                    var->_arena = NULL;
                } else {
                    for (size_t j = 0; j < i; j++) {
                        free(str_ptr[j]);
                    }
                }
                ctools_aligned_free(var->data.str);
                var->data.str = NULL;
                return STATA_ERR_MEMORY;
            }
        }
    } else {
        /* Numeric variable */
        var->type = STATA_TYPE_DOUBLE;
        var->_arena = NULL;

        size_t dbl_array_size;
        if (ctools_safe_mul_size(n_filtered, sizeof(double), &dbl_array_size) != 0) {
            return STATA_ERR_MEMORY;
        }
        var->data.dbl = (double *)ctools_cacheline_alloc(dbl_array_size);
        if (var->data.dbl == NULL) {
            return STATA_ERR_MEMORY;
        }

        double * restrict dbl_ptr = var->data.dbl;

        /* Cache SPI function pointer — avoid reloading _stata_ per iteration */
        ST_IIIDp vdata_fn = IO_VDATA_FN;
        /* Load using obs_map - obs_map contains 1-based Stata obs numbers */
        for (i = 0; i + 4 <= n_filtered; i += 4) {
            { int io_rc = vdata_fn((ST_int)var_idx, (ST_int)obs_map[i + 0], &dbl_ptr[i + 0]); if (io_rc) return ctools_spi_status(io_rc, STATA_ERR_STATA_READ); }
            { int io_rc = vdata_fn((ST_int)var_idx, (ST_int)obs_map[i + 1], &dbl_ptr[i + 1]); if (io_rc) return ctools_spi_status(io_rc, STATA_ERR_STATA_READ); }
            { int io_rc = vdata_fn((ST_int)var_idx, (ST_int)obs_map[i + 2], &dbl_ptr[i + 2]); if (io_rc) return ctools_spi_status(io_rc, STATA_ERR_STATA_READ); }
            { int io_rc = vdata_fn((ST_int)var_idx, (ST_int)obs_map[i + 3], &dbl_ptr[i + 3]); if (io_rc) return ctools_spi_status(io_rc, STATA_ERR_STATA_READ); }
        }
        for (; i < n_filtered; i++) {
            { int io_rc = vdata_fn((ST_int)var_idx, (ST_int)obs_map[i + 0], &dbl_ptr[i + 0]); if (io_rc) return ctools_spi_status(io_rc, STATA_ERR_STATA_READ); }
        }
    }

    return 0;
}

void *ctools_io_load_filtered_variable_thread(void *arg)
{
    filtered_var_io_args *args = (filtered_var_io_args *)arg;
    args->error = ctools_io_load_filtered_variable(args->var, args->var_idx, args->obs_map,
                                         args->n_filtered, args->is_string, args->str_width);
    return args->error ? (void *)1 : NULL;
}

int ctools_io_load_tiled_variables(stata_data *data, int *indices,
                                 const perm_idx_t *map, const int *widths)
{
    const size_t n = data->nobs, k = data->nvars;
    for (size_t j = 0; j < k; j++) {
        stata_variable *v = data->vars + j;
        int idx = indices[j];
        v->nobs = n;
        if (SF_var_is_string(idx)) {
            size_t stride = (size_t)widths[idx - 1] + 1;
            v->type = STATA_TYPE_STRING;
            v->str_maxlen = stride - 1;
            v->data.str = (char **)ctools_safe_cacheline_alloc2(n, sizeof(char *));
            if (!v->data.str) return STATA_ERR_MEMORY;
            memset(v->data.str, 0, n * sizeof(char *));
            ctools_string_arena *arena = calloc(1, sizeof(*arena));
            if (!arena) return STATA_ERR_MEMORY;
            v->_arena = arena;
            arena->base = calloc(n + (IO_SLOT_GUARD + stride - 1) / stride, stride);
            if (!arena->base) return STATA_ERR_MEMORY;
            arena->capacity = arena->used = n * stride;
            arena->mode = CTOOLS_STRING_ARENA_NO_FALLBACK;
        } else {
            v->type = STATA_TYPE_DOUBLE;
            v->data.dbl = ctools_safe_cacheline_alloc2(n, sizeof(double));
            if (!v->data.dbl) return STATA_ERR_MEMORY;
        }
    }
    atomic_int error = 0;
    #pragma omp parallel for schedule(static)
    for (size_t tile = 0; tile < n; tile += IO_LOAD_TILE_ROWS) {
        size_t end = tile + IO_LOAD_TILE_ROWS;
        if (end > n) end = n;
        for (size_t j = 0; j < k; j++) {
            stata_variable *v = data->vars + j;
            ST_int idx = indices[j];
            if (v->type == STATA_TYPE_DOUBLE) {
                ST_IIIDp read_fn = IO_VDATA_FN;
                for (size_t i = tile; i < end; i++)
                    { int read_rc = read_fn(idx, (ST_int)map[i], v->data.dbl + i);
                    if (read_rc) {
                        ctools_io_record_io_error(&error, ctools_spi_status(read_rc, STATA_ERR_STATA_READ));
                        break;
                    }
                    }
            } else {
                ctools_string_arena *arena = v->_arena;
                size_t stride = v->str_maxlen + 1;
                ST_IIIS read_fn = (_stata_)->sdata;
                /* Tiles may belong to different workers, so direct reads are
                 * bounded by this tile; the guard tail covers the last one. */
                const char *tile_safe_end = arena->base + end * stride
                    + (end == n ? ((IO_SLOT_GUARD + stride - 1) / stride) * stride : 0);
                for (size_t i = tile; i < end; i++) {
                    char *slot = arena->base + i * stride;
                    int rc = ctools_io_read_string_slot_direct(read_fn, idx, (ST_int)map[i], slot, stride - 1, tile_safe_end);
                    v->data.str[i] = slot;
                    if (rc) {
                        ctools_io_record_io_error(&error, rc);
                        break;
                    }
                }
            }
        }
    }
    return atomic_load(&error);
}

stata_retcode ctools_data_load_ex(ctools_filtered_data *result,
                                   int *var_indices, size_t nvars,
                                   size_t obs_start, size_t obs_end,
                                   int flags, const int *str_widths)
{
    if (!result) return STATA_ERR_INVALID_INPUT;
    ctools_filtered_data_init(result);
    ctools_load_plan plan;
    stata_retcode rc = ctools_io_plan_load(&plan, var_indices, nvars,
                                           obs_start, obs_end, flags, str_widths);
    if (rc != STATA_OK || plan.nvars == 0) return rc;
    var_indices = plan.indices;
    nvars = plan.nvars;
    str_widths = plan.widths;
    int *auto_indices = plan.owned_indices, *auto_str_widths = plan.owned_widths;
    perm_idx_t *obs_map = plan.obs_map;
    size_t n_range = plan.n_range, n_filtered = plan.nobs;
    int was_filtered = plan.was_filtered;

    result->n_range = n_range;
    result->obs_map = obs_map;
    result->was_filtered = was_filtered;

    /* Handle empty result */
    if (n_filtered == 0) {
        if (n_range > 0) {
            /* Some obs in range but all filtered out */
            rc = ctools_io_init_data_structure(&result->data, nvars, 0, flags);
        }
        free(auto_str_widths);
        free(auto_indices);
        return ctools_io_load_status(rc);
    }

    /* Allocate data structure */
    rc = ctools_io_init_data_structure(&result->data, nvars, n_filtered, flags);
    if (rc != STATA_OK) {
        ctools_filtered_data_free(result);
        free(auto_str_widths);
        free(auto_indices);
        return ctools_io_load_status(rc);
    }

    int tile_ok = 0;
    #ifdef _OPENMP
    tile_ok = nvars >= IO_MIN_TILED_VARS &&
        (n_filtered >= MIN_OBS_PER_THREAD * 2 ||
         (nvars >= IO_NUMERIC_TILE_MIN_VARS && n_filtered >= IO_NUMERIC_TILE_MIN_ROWS));
    #endif
    size_t string_bytes = 0;
    for (size_t j = 0; j < nvars && tile_ok; j++) {
        int idx = var_indices[j];
        if (SF_var_is_string(idx)) {
            int width = str_widths ? str_widths[idx - 1] : 0;
            if (width > 0) string_bytes += (size_t)width;
            if (SF_var_is_strl(idx) || width < 1 || width >= STATA_STR_MAXLEN ||
                n_filtered > (2048ULL * 1024 * 1024) / ((size_t)width + 1)) tile_ok = 0;
        }
    }
    if (tile_ok &&
        ((n_filtered >= MIN_OBS_PER_THREAD * 2 && string_bytes >= IO_MIN_TILED_STRING_BYTES) ||
         (string_bytes == 0 && nvars >= IO_NUMERIC_TILE_MIN_VARS && n_filtered >= IO_NUMERIC_TILE_MIN_ROWS))) {
        rc = ctools_io_load_tiled_variables(&result->data, var_indices, obs_map, str_widths);
        if (rc != STATA_ERR_MEMORY) {
            if (rc) ctools_filtered_data_free(result);
            free(auto_str_widths); free(auto_indices);
            return ctools_io_load_status(rc);
        }
        /* Speculative flat allocations must not remove the existing arena
         * fallback when address space or available memory is constrained. */
        stata_data_free(&result->data);
        rc = ctools_io_init_data_structure(&result->data, nvars, n_filtered, flags);
        if (rc) {
            ctools_filtered_data_free(result);
            free(auto_str_widths); free(auto_indices);
            return ctools_io_load_status(rc);
        }
    }

    /*
        Auto-dispatch: row-parallel for single variable, column-parallel for multiple.
        For one column, or two/three columns with ample idle workers, split
        observations across threads. Otherwise assign columns to workers.
    */
    int few_columns = 0;
    #ifdef _OPENMP
    few_columns = nvars >= 2 && nvars <= 3 && n_filtered >= MIN_OBS_PER_THREAD * 2 &&
                  ctools_get_max_threads() >= (int)(2 * nvars);
    #endif
    int chunk_handled = 0;
    if (!(nvars == 1 || few_columns)) {
        chunk_handled = ctools_io_try_chunked_load(&result->data, var_indices, nvars,
                                         obs_map, n_filtered, was_filtered, &rc);
        if (chunk_handled && rc) {
            ctools_filtered_data_free(result);
            free(auto_str_widths);
            free(auto_indices);
            return ctools_io_load_status(rc);
        }
    }
    if (chunk_handled) {
        /* loaded via chunked scheduler */
    } else if (nvars == 1 || few_columns) {
        for (size_t j = 0; j < nvars && !rc; j++) {
            int width = str_widths ? str_widths[var_indices[j] - 1] : 0;
            rc = ctools_io_load_row_parallel(&result->data.vars[j], var_indices[j], obs_map, n_filtered, width);
        }
        if (rc) {
            ctools_filtered_data_free(result);
            free(auto_str_widths);
            free(auto_indices);
            return ctools_io_load_status(rc);
        }
    } else if (!was_filtered) {
        /*
            Identity case (nvars >= 2): all observations pass, use contiguous load.
            obs_map[0] gives the 1-based start observation.
        */
        ctools_var_io_args *thread_args = (ctools_var_io_args *)
            ctools_safe_malloc2(nvars, sizeof(ctools_var_io_args));
        if (thread_args == NULL) {
            ctools_filtered_data_free(result);
            free(auto_str_widths);
            free(auto_indices);
            return STATA_ERR_MEMORY;
        }

        ST_int obs1 = (ST_int)obs_map[0];
        ctools_io_init_io_thread_args(thread_args, &result->data, var_indices, nvars,
                            obs1, n_filtered, IO_MODE_LOAD, str_widths);

        rc = ctools_io_execute_io_parallel(thread_args, nvars, ctools_io_load_variable_thread);
        if (rc != 0) {
            free(thread_args);
            ctools_filtered_data_free(result);
            free(auto_str_widths);
            free(auto_indices);
            return ctools_io_load_status(rc);
        }
        free(thread_args);
    } else {
        /*
            Filtered case (nvars >= 2): use obs_map for gather.
        */
        filtered_var_io_args *thread_args = (filtered_var_io_args *)
            ctools_safe_malloc2(nvars, sizeof(filtered_var_io_args));
        if (thread_args == NULL) {
            ctools_filtered_data_free(result);
            free(auto_str_widths);
            free(auto_indices);
            return STATA_ERR_MEMORY;
        }

        for (size_t j = 0; j < nvars; j++) {
            thread_args[j].var = &result->data.vars[j];
            thread_args[j].var_idx = var_indices[j];
            thread_args[j].obs_map = obs_map;
            thread_args[j].n_filtered = n_filtered;
            thread_args[j].is_string = ctools_io_variable_kind(var_indices[j]);
            thread_args[j].str_width = str_widths ? str_widths[var_indices[j] - 1] : 0;
            thread_args[j].error = STATA_OK;
        }

        size_t parallel_count = 0;
        for (size_t j = 0; j < nvars; j++) {
            if (thread_args[j].is_string == IO_LONG_STRING) {
                ctools_io_load_filtered_variable_thread(&thread_args[j]);
                if (thread_args[j].error) {
                    rc = thread_args[j].error;
                    free(thread_args);
                    ctools_filtered_data_free(result);
                    free(auto_str_widths);
                    free(auto_indices);
                    return ctools_io_load_status(rc);
                }
            } else {
                thread_args[parallel_count++] = thread_args[j];
            }
        }

        size_t row_work = 0;
        for (size_t j = 0; j < parallel_count; j++) {
            int width = thread_args[j].str_width;
            row_work = ctools_io_add_work(row_work, ctools_io_column_work(thread_args[j].is_string,
                                  width > 0 ? (size_t)width : 0));
        }
        int serial = ctools_io_serial_work(n_filtered, parallel_count, row_work);
        ctools_persistent_pool *pool = parallel_count >= 2 && !serial ? ctools_get_global_pool() : NULL;
        if (pool != NULL) {
            if (ctools_persistent_pool_submit_batch(pool, ctools_io_load_filtered_variable_thread,
                                                     thread_args, parallel_count,
                                                     sizeof(filtered_var_io_args)) != 0) {
                free(thread_args);
                ctools_filtered_data_free(result);
                free(auto_str_widths);
                free(auto_indices);
                return STATA_ERR_MEMORY;
            }
            if (ctools_persistent_pool_wait(pool)) rc = STATA_ERR_MEMORY;
        } else {
            for (size_t j = 0; j < parallel_count; j++) ctools_io_load_filtered_variable_thread(&thread_args[j]);
        }
        for (size_t j = 0; j < parallel_count; j++) {
            if (thread_args[j].error) { rc = thread_args[j].error; break; }
        }
        if (rc) {
            free(thread_args);
            ctools_filtered_data_free(result);
            free(auto_str_widths);
            free(auto_indices);
            return ctools_io_load_status(rc);
        }

        free(thread_args);
    }

    ctools_memory_barrier();

    free(auto_str_widths);
    free(auto_indices);
    return STATA_OK;
}

stata_retcode ctools_data_load(ctools_filtered_data *result,
                                         int *var_indices, size_t nvars,
                                         size_t obs_start, size_t obs_end,
                                         int flags)
{
    return ctools_data_load_ex(result, var_indices, nvars, obs_start, obs_end,
                                flags, NULL);
}

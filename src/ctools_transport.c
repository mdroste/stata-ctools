#include "ctools_transport_internal.h"

stata_retcode ctools_io_validate_variable(int var_idx, int expected_string)
{
    if (var_idx < 1 || var_idx > SF_nvars())
        return STATA_ERR_INVALID_INPUT;
    int is_string = SF_var_is_string(var_idx);
    if (expected_string >= 0 && is_string != expected_string)
        return STATA_ERR_UNSUPPORTED_TYPE;
    if (expected_string >= 0 && is_string && SF_var_is_strl(var_idx))
    {
        SF_error("ctools: writing strL is unsupported; recast to str2045 when lossless\n");
        return STATA_ERR_UNSUPPORTED_TYPE;
    }
    return STATA_OK;
}

stata_retcode ctools_io_validate_range(size_t obs1, size_t nobs)
{
    size_t available = (size_t)SF_nobs();
    if (obs1 < 1 || obs1 > available + 1 || nobs > available - (obs1 - 1))
        return STATA_ERR_INVALID_INPUT;
    return STATA_OK;
}

stata_retcode ctools_io_validate_obs_map(const perm_idx_t *obs_map, size_t nobs)
{
    if (nobs && !obs_map)
        return STATA_ERR_INVALID_INPUT;
    size_t available = (size_t)SF_nobs();
    for (size_t i = 0; i < nobs; i++)
        if (obs_map[i] < 1 || obs_map[i] > available)
            return STATA_ERR_INVALID_INPUT;
    return STATA_OK;
}

int ctools_io_read_stata_strl(ST_int var_idx, ST_int obs, char *value)
{
    if (SF_var_is_binary(var_idx, obs))
        return STATA_ERR_UNSUPPORTED_TYPE;
    ST_int length = SF_sdatalen(var_idx, obs);
    if (length < 0)
        return STATA_ERR_STATA_READ;
    if (length > STATA_STR_MAXLEN)
        return STATA_ERR_UNSUPPORTED_TYPE;
    /* SF_strldata returns bytes copied (not a zero-on-success status). */
    if (SF_strldata(var_idx, obs, value, length + 1) != length)
        return STATA_ERR_STATA_READ;
    value[length] = '\0';
    return STATA_OK;
}

int ctools_io_variable_kind(int var_idx)
{
    if (!SF_var_is_string(var_idx))
        return IO_NUMERIC;
    return SF_var_is_strl(var_idx) ? IO_LONG_STRING : IO_FIXED_STRING;
}

int ctools_io_read_string_slot(ST_IIIS read_fn, int var_idx, ST_int obs, char *slot, size_t width)
{
    char value[STATA_STR_MAXLEN + 1];
    int rc = read_fn(var_idx, obs, value);
    if (rc)
        return read_fn == ctools_io_read_stata_strl ? rc
                                                    : ctools_spi_status(rc, STATA_ERR_STATA_READ);
    size_t len = strlen(value);
    if (len > width)
        return STATA_ERR_STATA_READ;
    memcpy(slot, value, len + 1);
    return STATA_OK;
}

int ctools_io_read_string_slot_direct(ST_IIIS read_fn, int var_idx, ST_int obs, char *slot,
                                      size_t width, const char *safe_end)
{
    if (width >= IO_DIRECT_READ_MIN_WIDTH && (size_t)(safe_end - slot) >= IO_SLOT_GUARD)
    {
        int rc = read_fn(var_idx, obs, slot);
        if (rc)
            return read_fn == ctools_io_read_stata_strl
                       ? rc
                       : ctools_spi_status(rc, STATA_ERR_STATA_READ);
        if (strlen(slot) > width)
            return STATA_ERR_STATA_READ;
        return STATA_OK;
    }
    return ctools_io_read_string_slot(read_fn, var_idx, obs, slot, width);
}

stata_retcode ctools_io_load_status(int rc)
{
    if (rc == STATA_ERR_UNSUPPORTED_TYPE)
        SF_error("ctools: binary strL or text exceeding 2045 bytes is unsupported\n");
    else if (rc == STATA_ERR_STATA_READ)
        SF_error("ctools: failed to read Stata data\n");
    return (stata_retcode)rc;
}

void ctools_io_record_io_error(atomic_int *error, int rc)
{
    if (rc)
    {
        int expected = 0;
        atomic_compare_exchange_strong(error, &expected, rc);
    }
}

/* Longest __ctools_strw read: above Stata's macro length limit (c(macrolen)). */
#define IO_STRW_MAX_BYTES ((size_t)8 * 1024 * 1024)

int *ctools_io_ctools_read_strw_from_stata(size_t nvars)
{
    /* SF_macro_use convention: prefix local name with underscore.
       Local "__ctools_strw" -> SPI name "___ctools_strw".
       SF_macro_use cuts a longer value to the buffer without saying so, and a
       cut entry reads as a narrower width ("2045" as "20"). A value that
       reaches the end of the buffer may have been cut, so it is read again
       into a larger one; metadata that never fits is no hint. */
    char stack_buf[16384]; /* ~3200 variables at max width "2045," */
    char *buf = stack_buf, *heap = NULL;
    size_t cap = sizeof(stack_buf);
    for (;;)
    {
        buf[0] = buf[cap - 1] = '\0';
        ST_int rc = SF_macro_use("___ctools_strw", buf, (ST_int)(cap - 1));
        if (rc != 0 || buf[0] == '\0')
        {
            free(heap);
            return NULL; /* Macro not set — .ado didn't call _ctools_strw */
        }
        if (strlen(buf) + 2 < cap)
            break; /* the whole value */
        free(heap);
        heap = NULL;
        if (cap > IO_STRW_MAX_BYTES / 2)
            return NULL;
        cap *= 2;
        heap = (char *)malloc(cap);
        if (!heap)
            return NULL;
        buf = heap;
    }

    int *widths = (int *)calloc(nvars, sizeof(int));
    if (!widths)
    {
        free(heap);
        return NULL;
    }

    const char *p = buf;
    for (size_t j = 0; j < nvars && *p != '\0'; j++)
    {
        char *endptr;
        long val = strtol(p, &endptr, 10);
        widths[j] = (val > 0 && val <= STATA_STR_MAXLEN) ? (int)val : 0;
        p = endptr;
        if (*p == ',')
            p++;
    }

    free(heap);
    return widths;
}

char *ctools_io_copy_private_string(ctools_string_arena *arena, const char *value)
{
    size_t len = strlen(value) + 1;
    if (arena && len <= arena->capacity - arena->used)
    {
        char *result = arena->base + arena->used;
        arena->used += len;
        memcpy(result, value, len);
        return result;
    }
    char *result = strdup(value);
    if (result && arena)
        arena->has_fallback = 1;
    return result;
}

ctools_string_arena *ctools_io_create_private_string_arena(size_t nobs, int str_width)
{
    const size_t original_capacity = nobs * 64;
    size_t capacity = original_capacity;
    if (str_width == STATA_STR_MAXLEN)
    {
        const size_t cap = (size_t)2048 * 1024 * 1024;
        const size_t stride = STATA_STR_MAXLEN + 1;
        capacity = nobs <= cap / stride ? nobs * stride : cap;
    }
    ctools_string_arena *arena =
        ctools_string_arena_create(capacity, CTOOLS_STRING_ARENA_STRDUP_FALLBACK);
    if (!arena && capacity != original_capacity)
        arena = ctools_string_arena_create(original_capacity, CTOOLS_STRING_ARENA_STRDUP_FALLBACK);
    return arena;
}

stata_retcode ctools_io_init_data_structure(stata_data *data, size_t nvars, size_t nobs, int flags)
{
    /* Initialize the data structure */
    stata_data_init(data);
    data->nobs = nobs;
    data->nvars = nvars;

    /* Allocate array of variables (cache-line aligned, overflow-safe) */
    data->vars = (stata_variable *)ctools_safe_cacheline_alloc2(nvars, sizeof(stata_variable));
    if (data->vars == NULL)
    {
        return STATA_ERR_MEMORY;
    }
    memset(data->vars, 0, nvars * sizeof(stata_variable));

    /* Read-only/column-transform callers can omit the unused row order. */
    if (flags & CTOOLS_LOAD_NO_SORT_ORDER)
        return STATA_OK;

    /* Allocate sort order array (cache-line aligned, overflow-safe) */
    data->sort_order = (perm_idx_t *)ctools_safe_cacheline_alloc2(nobs, sizeof(perm_idx_t));
    if (data->sort_order == NULL)
    {
        stata_data_free(data);
        return STATA_ERR_MEMORY;
    }

/* Initialize sort order to identity permutation (parallel for large datasets) */
#pragma omp parallel for schedule(static) if (nobs >= MIN_OBS_PER_THREAD * 2)
    for (size_t i = 0; i < nobs; i++)
    {
        data->sort_order[i] = (perm_idx_t)i;
    }

    return STATA_OK;
}

void ctools_io_init_io_thread_args(ctools_var_io_args *args, stata_data *data, int *var_indices,
                                   size_t nvars, size_t obs1, size_t nobs, io_mode_t mode,
                                   const int *str_widths)
{
    for (size_t j = 0; j < nvars; j++)
    {
        args[j].var = &data->vars[j];
        args[j].var_idx = var_indices ? var_indices[j] : (int)(j + 1);
        args[j].obs1 = obs1;
        args[j].nobs = nobs;
        args[j].is_string = (mode == IO_MODE_LOAD) ? ctools_io_variable_kind(args[j].var_idx)
                                                   : 0; /* Store uses var->str_maxlen instead */
        args[j].str_width =
            (str_widths && mode == IO_MODE_LOAD) ? str_widths[args[j].var_idx - 1] : 0;
        args[j].error = STATA_OK;
        args[j].store_order = NULL;
    }
}

size_t ctools_io_column_work(int is_string, size_t width)
{
    if (!is_string)
        return 1;
    if (width < 1 || width > STATA_STR_MAXLEN)
        width = STATA_STR_MAXLEN;
    return 2 + (width + 7) / 8;
}

size_t ctools_io_add_work(size_t total, size_t column)
{
    return total > SIZE_MAX - column ? SIZE_MAX : total + column;
}

int ctools_io_serial_work(size_t nobs, size_t nvars, size_t row_work)
{
    if (nvars < 2 || row_work == 0)
        return 0;
    size_t budget = nvars > (SIZE_MAX - IO_SERIAL_BASE_WORK) / IO_SERIAL_COLUMN_WORK
                        ? SIZE_MAX
                        : IO_SERIAL_BASE_WORK + IO_SERIAL_COLUMN_WORK * nvars;
    return nobs <= budget / row_work;
}

int ctools_io_execute_io_parallel(ctools_var_io_args *args, size_t nvars, io_thread_func func)
{
    if (func == ctools_io_load_variable_thread)
    {
        size_t parallel_count = 0;
        for (size_t j = 0; j < nvars; j++)
        {
            if (args[j].is_string == IO_LONG_STRING)
            {
                func(&args[j]);
                if (args[j].error)
                    return args[j].error;
            }
            else
            {
                args[parallel_count++] = args[j];
            }
        }
        nvars = parallel_count;
    }
    size_t row_work = 0;
    for (size_t j = 0; j < nvars; j++)
    {
        int is_string = func == ctools_io_load_variable_thread
                            ? args[j].is_string
                            : args[j].var->type == STATA_TYPE_STRING;
        size_t width = func == ctools_io_load_variable_thread
                           ? (args[j].str_width > 0 ? (size_t)args[j].str_width : 0)
                           : args[j].var->str_maxlen;
        row_work = ctools_io_add_work(row_work, ctools_io_column_work(is_string, width));
    }
    int serial = nvars >= 2 && ctools_io_serial_work(args[0].nobs, nvars, row_work);
    ctools_persistent_pool *pool = nvars >= 2 && !serial ? ctools_get_global_pool() : NULL;
    if (pool)
    {
        if (ctools_persistent_pool_submit_batch(pool, func, args, nvars,
                                                sizeof(ctools_var_io_args)) != 0)
            return STATA_ERR_MEMORY;
        int pool_result = ctools_persistent_pool_wait(pool);
        for (size_t j = 0; j < nvars; j++)
            if (args[j].error)
                return args[j].error;
        if (pool_result)
            return STATA_ERR_MEMORY;
    }
    else
    {
        for (size_t j = 0; j < nvars; j++)
        {
            func(&args[j]);
            if (args[j].error)
                return args[j].error;
        }
    }
    return STATA_OK;
}

size_t ctools_io_chunks_per_column(size_t nvars, size_t nobs, size_t min_chunk_rows)
{
    size_t nthreads = (size_t)ctools_get_max_threads();
    size_t per_col = (4 * nthreads + nvars - 1) / nvars;
    size_t max_per_col = nobs / min_chunk_rows;
    if (per_col > max_per_col)
        per_col = max_per_col;
    return per_col;
}

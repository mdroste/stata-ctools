/*
    ctools_data_io.c
    Stata <-> C Data Transfer Module

    Handles bidirectional data transfer between Stata's data space and C memory
    using the Stata Plugin Interface (SPI). This module is designed to be
    reusable for any operation requiring bulk data transfer, not just sorting.

    Architecture:
    - Column-major storage: Each variable is stored in a contiguous array
    - Parallel I/O: One thread per variable when nvars >= 2
    - Uses SF_vdata()/SF_sdata() for reading, SF_vstore()/SF_sstore() for writing
    - SPI convention: variable index first, then observation index (both 1-based)

    Thread Safety:
    - Each thread operates on a different Stata variable (column)
    - No shared state during I/O phase
    - Assumes Stata's SPI is thread-safe for column-level parallelism
    - Note: Custom streaming code (e.g., in cmerge) should use sequential I/O
      to avoid potential race conditions with complex access patterns

    Performance Optimizations:
    - Parallel variable I/O overlaps operations across columns
    - Checked SPI access remains enabled even in SD_FASTMODE builds
    - Cache-line aligned allocations for optimal memory access
    - String arena allocator for reduced malloc overhead
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <stdatomic.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_threads.h"
#include "ctools_config.h"
#include "ctools_arena.h"

/* Maximum string buffer size for Stata string variables */
#define STATA_STR_MAXLEN 2045

/* Validate plugin-visible indices before type queries or unchecked string SPI.
 * Writes require fixed strings; bounded textual strL reads use the strL API. */
static stata_retcode validate_variable(int var_idx, int expected_string)
{
    if (var_idx < 1 || var_idx > SF_nvars()) return STATA_ERR_INVALID_INPUT;
    int is_string = SF_var_is_string(var_idx);
    if (expected_string >= 0 && is_string != expected_string)
        return STATA_ERR_UNSUPPORTED_TYPE;
    if (expected_string >= 0 && is_string && SF_var_is_strl(var_idx)) {
        SF_error("ctools: writing strL is unsupported; recast to str2045 when lossless\n");
        return STATA_ERR_UNSUPPORTED_TYPE;
    }
    return STATA_OK;
}

static stata_retcode validate_range(size_t obs1, size_t nobs)
{
    size_t available = (size_t)SF_nobs();
    if (obs1 < 1 || obs1 > available + 1 || nobs > available - (obs1 - 1))
        return STATA_ERR_INVALID_INPUT;
    return STATA_OK;
}

static stata_retcode validate_obs_map(const perm_idx_t *obs_map, size_t nobs)
{
    if (nobs && !obs_map) return STATA_ERR_INVALID_INPUT;
    size_t available = (size_t)SF_nobs();
    for (size_t i = 0; i < nobs; i++)
        if (obs_map[i] < 1 || obs_map[i] > available) return STATA_ERR_INVALID_INPUT;
    return STATA_OK;
}

/* Read bounded text strL through SPI's length-aware interface. Binary and
 * oversized values are rejected before touching the fixed-width buffer.
 * See https://www.stata.com/plugins/ (Routines for handling data). */
static int read_stata_string(ST_int var_idx, ST_int obs, char *value)
{
    if (!SF_var_is_strl(var_idx))
        return SF_sdata(var_idx, obs, value) ? STATA_ERR_STATA_READ : STATA_OK;
    if (SF_var_is_binary(var_idx, obs)) return STATA_ERR_UNSUPPORTED_TYPE;
    ST_int length = SF_sdatalen(var_idx, obs);
    if (length < 0) return STATA_ERR_STATA_READ;
    if (length > STATA_STR_MAXLEN) return STATA_ERR_UNSUPPORTED_TYPE;
    /* SF_strldata returns bytes copied (not a zero-on-success status). */
    if (SF_strldata(var_idx, obs, value, length + 1) != length) return STATA_ERR_STATA_READ;
    value[length] = '\0';
    return STATA_OK;
}

/* Width hints are not trusted as buffer bounds: a stale hint must never allow
 * SPI to overwrite the next slot. */
static int read_string_slot(int var_idx, ST_int obs, char *slot, size_t width)
{
    char value[STATA_STR_MAXLEN + 1];
    int rc = read_stata_string(var_idx, obs, value);
    if (rc) return rc;
    size_t len = strlen(value);
    if (len > width) return STATA_ERR_STATA_READ;
    memcpy(slot, value, len + 1);
    return STATA_OK;
}

static stata_retcode load_status(int rc)
{
    if (rc == STATA_ERR_UNSUPPORTED_TYPE)
        SF_error("ctools: binary strL or text exceeding 2045 bytes is unsupported\n");
    else if (rc == STATA_ERR_STATA_READ)
        SF_error("ctools: failed to read Stata data\n");
    return (stata_retcode)rc;
}

/* Only workers write the error slot; compare-exchange preserves the first
 * failure without volatile data races. The caller reads it after joining. */
static void record_io_error(atomic_int *error, int rc)
{
    if (rc) {
        int expected = 0;
        atomic_compare_exchange_strong(error, &expected, rc);
    }
}

static stata_retcode store_status(int rc)
{
    if (rc) SF_error("ctools: write failed; data may be partially modified\n");
    return (stata_retcode)rc;
}

/* ===========================================================================
   String Width Auto-Detection via Stata Local Macro

   When the calling .ado runs `_ctools_strw varlist` before `plugin call`,
   it sets a Stata local `__ctools_strw` with comma-separated string widths
   (0 = numeric, >0 = actual width like 17 for str17).

   ctools_data_load reads this local automatically via SF_macro_use to
   enable flat buffer string I/O without per-command changes.
   =========================================================================== */

/*
    Read and parse the __ctools_strw Stata local set by _ctools_strw.ado.
    Returns malloc'd array of nvars ints (0 = numeric), or NULL if not available.
    Caller must free.
*/
static int *ctools_read_strw_from_stata(size_t nvars)
{
    /* SF_macro_use convention: prefix local name with underscore.
       Local "__ctools_strw" -> SPI name "___ctools_strw" */
    char buf[16384];  /* ~3200 variables at max width "2045," */
    ST_int rc = SF_macro_use("___ctools_strw", buf, sizeof(buf) - 1);
    if (rc != 0 || buf[0] == '\0') {
        return NULL;  /* Macro not set — .ado didn't call _ctools_strw */
    }

    int *widths = (int *)calloc(nvars, sizeof(int));
    if (!widths) return NULL;

    const char *p = buf;
    for (size_t j = 0; j < nvars && *p != '\0'; j++) {
        char *endptr;
        long val = strtol(p, &endptr, 10);
        widths[j] = (val > 0 && val < STATA_STR_MAXLEN) ? (int)val : 0;
        p = endptr;
        if (*p == ',') p++;
    }

    return widths;
}

/* ===========================================================================
   Data Loading (Stata -> C)
   =========================================================================== */

/*
    Load a single variable from Stata into C memory (internal helper).
    Used by both single-var and multi-var loading.

    When str_width > 0 (actual Stata variable width, e.g. 17 for str17),
    uses a flat contiguous buffer for strings instead of the arena allocator.
    This eliminates per-string strlen, memcpy, and atomic CAS overhead.
*/
static int load_single_variable(stata_variable *var, int var_idx, size_t obs1,
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
            /* Allocate minimal buffer to avoid NULL - calloc ensures it's zeroed */
            var->data.str = (char **)calloc(1, sizeof(char *));
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
        var->type = STATA_TYPE_STRING;

        /* ---- Flat buffer path: known string width ----
           When the actual Stata variable width is known (e.g. 17 for str17),
           allocate a single contiguous nobs×stride buffer and have SF_sdata
           read directly into it. Benefits:
           - Single copy (SPI → flat buffer) instead of two (SPI → stack → arena)
           - No strlen per string (stride is known)
           - No atomic CAS overhead (no arena allocation)
           - Better cache locality for subsequent sort comparisons
           - O(1) cleanup via ctools_string_arena_free wrapper */
        if (str_width > 0 && str_width < STATA_STR_MAXLEN) {
            size_t stride = (size_t)str_width + 1;
            size_t flat_size = nobs * stride;

            /* Overflow and size guard (2 GB cap) */
            if (flat_size / stride == nobs && flat_size <= (2048ULL * 1024 * 1024)) {
                char *flat_buf = (char *)calloc(nobs, stride);
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

                    /* Read directly from Stata into flat buffer — single copy */
                    for (i = 0; i < nobs; i++) {
                        str_ptrs[i] = flat_buf + i * stride;
                        int read_rc = read_string_slot(var_idx, (ST_int)(i + obs1), str_ptrs[i], (size_t)str_width);
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
            size_t arena_capacity = nobs * 64;
            arena = ctools_string_arena_create(arena_capacity, CTOOLS_STRING_ARENA_STRDUP_FALLBACK);
        }
        if (arena != NULL) {
            var->_arena = arena;
        }

        /* Load strings */
        ST_IIIS sdata_fn2 = read_stata_string;
        for (i = 0; i < nobs; i++) {
            int read_rc = sdata_fn2((ST_int)var_idx, (ST_int)(i + obs1), strbuf);
            if (read_rc) return read_rc;
            str_ptr[i] = ctools_string_arena_strdup(arena, strbuf);
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
        ST_IIIDp vdata_fn = (_stata_)->safevdata;
        for (i = 0; i < nobs; i++) {
            if (vdata_fn((ST_int)var_idx, (ST_int)(i + obs1), &dbl_ptr[i])) return STATA_ERR_STATA_READ;
        }
    }

    return 0;
}

/*
    Thread function: Load a single variable from Stata into C memory.

    Allocates memory and reads all observations for one variable using
    SF_vdata (numeric) or SF_sdata (string).

    @param arg  Pointer to ctools_var_io_args with input/output parameters
    @return     NULL on success, non-NULL on failure (for ctools_threads)
*/
static void *load_variable_thread(void *arg)
{
    ctools_var_io_args *args = (ctools_var_io_args *)arg;
    args->error = load_single_variable(args->var, args->var_idx, args->obs1,
                                       args->nobs, args->is_string, args->str_width);
    return args->error ? (void *)1 : NULL;
}

/* ===========================================================================
   Internal Helper Functions
   =========================================================================== */

/*
    Data structure initialization helper.
    Allocates vars array and sort_order with proper alignment.
    Returns STATA_OK on success, error code on failure.
*/
static stata_retcode init_data_structure(stata_data *data, size_t nvars, size_t nobs)
{
    /* Initialize the data structure */
    stata_data_init(data);
    data->nobs = nobs;
    data->nvars = nvars;

    /* Allocate array of variables (cache-line aligned, overflow-safe) */
    data->vars = (stata_variable *)ctools_safe_cacheline_alloc2(nvars, sizeof(stata_variable));
    if (data->vars == NULL) {
        return STATA_ERR_MEMORY;
    }
    memset(data->vars, 0, nvars * sizeof(stata_variable));

    /* Allocate sort order array (cache-line aligned, overflow-safe) */
    data->sort_order = (perm_idx_t *)ctools_safe_cacheline_alloc2(nobs, sizeof(perm_idx_t));
    if (data->sort_order == NULL) {
        stata_data_free(data);
        return STATA_ERR_MEMORY;
    }

    /* Initialize sort order to identity permutation (parallel for large datasets) */
    #pragma omp parallel for schedule(static) if(nobs >= MIN_OBS_PER_THREAD * 2)
    for (size_t i = 0; i < nobs; i++) {
        data->sort_order[i] = (perm_idx_t)i;
    }

    return STATA_OK;
}

/*
    I/O mode enum for thread args initialization.
*/
typedef enum { IO_MODE_LOAD, IO_MODE_STORE } io_mode_t;

/*
    Thread args batch initialization helper.
    Initializes an array of ctools_var_io_args for load or store operations.

    @param args         [out] Pre-allocated array of nvars args
    @param data         [in]  stata_data structure
    @param var_indices  [in]  Array of 1-based variable indices, or NULL for sequential
    @param nvars        [in]  Number of variables
    @param obs1         [in]  First observation (1-based)
    @param nobs         [in]  Number of observations
    @param mode         [in]  IO_MODE_LOAD or IO_MODE_STORE
    @param str_widths   [in]  Optional array of string widths per variable position
                              (0 = numeric/unknown). NULL = no width hints.
*/
static void init_io_thread_args(ctools_var_io_args *args, stata_data *data,
                                int *var_indices, size_t nvars,
                                size_t obs1, size_t nobs, io_mode_t mode,
                                const int *str_widths)
{
    for (size_t j = 0; j < nvars; j++) {
        args[j].var = &data->vars[j];
        args[j].var_idx = var_indices ? var_indices[j] : (int)(j + 1);
        args[j].obs1 = obs1;
        args[j].nobs = nobs;
        args[j].is_string = (mode == IO_MODE_LOAD)
            ? SF_var_is_string((ST_int)args[j].var_idx)
            : 0;  /* Store uses var->str_maxlen instead */
        args[j].str_width = (str_widths && mode == IO_MODE_LOAD)
            ? str_widths[args[j].var_idx - 1] : 0;
        args[j].error = STATA_OK;
    }
}

/*
    Thread function type for I/O operations.
*/
typedef void *(*io_thread_func)(void *);

/*
    Thread pool orchestration helper.
    Executes I/O operations either in parallel (using thread pool) or sequentially.

    @param args         [in]  Array of thread arguments
    @param nvars        [in]  Number of variables (and args)
    @param func         [in]  Thread function to execute

    @return 0 on success, non-zero on failure
*/
static int execute_io_parallel(ctools_var_io_args *args, size_t nvars,
                               io_thread_func func)
{
    ctools_persistent_pool *pool = nvars >= 2 ? ctools_get_global_pool() : NULL;
    if (pool) {
        if (ctools_persistent_pool_submit_batch(pool, func, args, nvars,
                                               sizeof(ctools_var_io_args)) != 0)
            return STATA_ERR_MEMORY;
        int pool_result = ctools_persistent_pool_wait(pool);
        for (size_t j = 0; j < nvars; j++)
            if (args[j].error) return args[j].error;
        if (pool_result) return STATA_ERR_MEMORY;
    } else {
        for (size_t j = 0; j < nvars; j++) {
            func(&args[j]);
            if (args[j].error) return args[j].error;
        }
    }
    return STATA_OK;
}

/* ===========================================================================
   Data Storing (C -> Stata)
   =========================================================================== */

/*
    Store a single variable from C memory to Stata (internal helper).
    Used by both single-var and multi-var storing.

    For string variables, var->str_maxlen holds the actual Stata variable
    width (e.g. 17 for str17) when set by the caller from .ado metadata.
    When available, we pack scattered arena pointers into a contiguous flat
    buffer so SF_sstore reads are sequential and cache-friendly.
*/
static int store_single_variable(stata_variable *var, int var_idx,
                                  size_t obs1, size_t nobs)
{
    size_t i;
    ST_int stata_var_idx = (ST_int)var_idx;

    if (var->type == STATA_TYPE_DOUBLE) {
        const double * restrict dbl_data = var->data.dbl;

        /* Cache SPI function pointer — avoid reloading _stata_ per iteration */
        ST_IIID store_fn = (_stata_)->safestore;
        for (i = 0; i < nobs; i++) {
            if (store_fn(stata_var_idx, (ST_int)(i + obs1), dbl_data[i])) return STATA_ERR_STATA_WRITE;
        }

    } else {
        /* String variable */
        char * const * restrict str_data = var->data.str;
        size_t str_width = var->str_maxlen;

        /* Cache SPI function pointer — avoid reloading _stata_ per iteration */
        ST_IIIS sstore_fn = (_stata_)->sstore;

        /* Use flat buffer when we have actual variable width (not default 2045)
           and the buffer fits in 256 MB */
        if (str_width > 0 && str_width < STATA_STR_MAXLEN) {
            size_t stride = str_width + 1;

            /* Repack path: strings may be scattered after permutation,
               copy into new flat buffer for sequential SF_sstore */
            size_t flat_size = nobs * stride;
            char *flat_buf = NULL;
            if (flat_size / stride == nobs &&  /* overflow check */
                flat_size <= (256ULL * 1024 * 1024)) {
                flat_buf = (char *)calloc(nobs, stride);
            }

            if (flat_buf) {
                for (i = 0; i < nobs; i++) {
                    const char *value = str_data[i];
                    if (value) {
                        size_t len = strlen(value);
                        if (len > str_width) {
                            free(flat_buf);
                            return STATA_ERR_STATA_WRITE;
                        }
                        memcpy(flat_buf + i * stride, value, len + 1);
                    }
                }

                /* Store: sequential scan through flat buffer */
                for (i = 0; i < nobs; i++) {
                    if (sstore_fn(stata_var_idx, (ST_int)(i + obs1), flat_buf + i * stride)) {
                        free(flat_buf);
                        return STATA_ERR_STATA_WRITE;
                    }
                }

                free(flat_buf);
                return STATA_OK;
            }
        }

        /* Fallback: width unknown, buffer too large, or alloc failed */
        for (i = 0; i < nobs; i++) {
            if (sstore_fn(stata_var_idx, (ST_int)(i + obs1), str_data[i] ? str_data[i] : "")) return STATA_ERR_STATA_WRITE;
        }
    }
    return STATA_OK;
}

/*
    Thread function: Store a single variable from C memory to Stata.

    Writes all observations for one variable using SF_vstore (numeric) or
    SF_sstore (string).

    @param arg  Pointer to ctools_var_io_args with input/output parameters
    @return     NULL on success (for ctools_threads compatibility)
*/
static void *store_variable_thread(void *arg)
{
    ctools_var_io_args *args = (ctools_var_io_args *)arg;
    args->error = store_single_variable(args->var, args->var_idx, args->obs1, args->nobs);
    return args->error ? (void *)1 : NULL;
}

/*
    Unified data store with auto-dispatch.

    Writes variables from C memory to Stata. Auto-dispatches between:
    - nvars == 1: row-parallel store (OpenMP over observations)
    - nvars >= 2: column-parallel store (thread pool, 1 thread per variable)

    @param data        [in]  stata_data structure containing data to write
    @param var_indices [in]  Array of 1-based Stata variable indices, or NULL
                             for sequential 1..data->nvars
    @param nvars       [in]  Number of variables to store (ignored if var_indices is NULL)
    @param obs1        [in]  First observation in Stata (1-based)

    @return STATA_OK on success, error code otherwise
*/
stata_retcode ctools_data_store_ex(stata_data *data, int *var_indices,
                                    size_t nvars, size_t obs1)
{
    if (data == NULL) {
        return STATA_ERR_INVALID_INPUT;
    }

    size_t nobs = data->nobs;

    /* If no var_indices, store all variables sequentially */
    if (var_indices == NULL) {
        nvars = data->nvars;
    } else if (nvars == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    if (!data->vars || nvars > data->nvars) return STATA_ERR_INVALID_INPUT;
    stata_retcode rc = validate_range(obs1, nobs);
    if (rc) return rc;
    for (size_t j = 0; j < nvars; j++) {
        stata_variable *var = &data->vars[j];
        if (var->type != STATA_TYPE_DOUBLE && var->type != STATA_TYPE_STRING) return STATA_ERR_UNSUPPORTED_TYPE;
        if (var->nobs < nobs || (nobs && !var->data.dbl)) return STATA_ERR_INVALID_INPUT;
        rc = validate_variable(var_indices ? var_indices[j] : (int)j + 1,
                               var->type == STATA_TYPE_STRING);
        if (rc) return rc;
    }

    /*
        nvars == 1: row-parallel store (OpenMP over observations).
        For numeric variables, this is a simple parallel SF_vstore loop.
        For string variables, use the existing store_single_variable which
        already handles flat buffer packing.
    */
    if (nvars == 1) {
        stata_variable *var = &data->vars[0];
        int var_idx = var_indices ? var_indices[0] : 1;
        ST_int stata_var = (ST_int)var_idx;

        atomic_int error = 0;
        if (var->type == STATA_TYPE_DOUBLE) {
            const double * restrict dbl_data = var->data.dbl;
            /* Cache SPI function pointer — avoid reloading _stata_ per iteration */
            ST_IIID store_fn = (_stata_)->safestore;
            #pragma omp parallel for schedule(static) if(nobs >= MIN_OBS_PER_THREAD * 2)
            for (size_t i = 0; i < nobs; i++) {
                if (store_fn(stata_var, (ST_int)(i + obs1), dbl_data[i])) record_io_error(&error, STATA_ERR_STATA_WRITE);
            }
        } else {
            /* String store: row-parallel from pointer array */
            char * const * restrict str_data = var->data.str;
            ST_IIIS sstore_fn = (_stata_)->sstore;
            #pragma omp parallel for schedule(static) if(nobs >= MIN_OBS_PER_THREAD * 2)
            for (size_t i = 0; i < nobs; i++) {
                if (sstore_fn(stata_var, (ST_int)(i + obs1), str_data[i] ? str_data[i] : "")) record_io_error(&error, STATA_ERR_STATA_WRITE);
            }
        }

        ctools_memory_barrier();
        return store_status(atomic_load(&error));
    }

    /* nvars >= 2: column-parallel store via thread pool */
    ctools_var_io_args *thread_args = (ctools_var_io_args *)
        ctools_safe_malloc2(nvars, sizeof(ctools_var_io_args));
    if (thread_args == NULL) {
        return STATA_ERR_MEMORY;
    }

    init_io_thread_args(thread_args, data, var_indices, nvars, obs1, nobs,
                        IO_MODE_STORE, NULL);
    rc = execute_io_parallel(thread_args, nvars, store_variable_thread);
    free(thread_args);

    ctools_memory_barrier();
    return store_status(rc);
}

/*
    Write all variables from C memory back to Stata's data space.
    Thin wrapper around ctools_data_store_ex with var_indices=NULL.
*/
stata_retcode ctools_data_store(stata_data *data, size_t obs1)
{
    return ctools_data_store_ex(data, NULL, 0, obs1);
}

/*
    Store variables to specified Stata variable indices.
    Thin wrapper around ctools_data_store_ex.
*/
stata_retcode ctools_data_store_selective(stata_data *data, int *var_indices,
                                           size_t nvars, size_t obs1)
{
    if (var_indices == NULL || nvars == 0) {
        return STATA_ERR_INVALID_INPUT;
    }
    return ctools_data_store_ex(data, var_indices, nvars, obs1);
}

/* ===========================================================================
   Streaming Permuted Variable Write (C -> Stata with reordering)
   =========================================================================== */

/*
    Write a single variable to Stata with row permutation.

    For each output row i, reads from Stata row source_rows[i] and writes
    to Stata row i. This enables in-place reordering without loading the
    entire variable into C memory first.

    Optimizations:
    - Cache-line aligned buffer allocation
    - String arena for reduced malloc overhead

    @param var_idx      [in] 1-based Stata variable index
    @param source_rows  [in] Array mapping output row -> source row (0-based)
                             -1 indicates missing value should be written
    @param output_nobs  [in] Number of output rows
    @param obs1         [in] First output observation in Stata (1-based)

    @return STATA_OK on success, error code otherwise
*/
stata_retcode ctools_stream_var_permuted(int var_idx, int64_t *source_rows,
                                          size_t output_nobs, size_t obs1)
{
    stata_retcode rc = validate_variable(var_idx, -1);
    if (rc) return rc;
    rc = validate_range(obs1, output_nobs);
    if (rc || (output_nobs && !source_rows)) return STATA_ERR_INVALID_INPUT;
    if (!output_nobs) return STATA_OK;
    size_t available = (size_t)SF_nobs();
    for (size_t i = 0; i < output_nobs; i++) {
        if (source_rows[i] < -1 ||
            (source_rows[i] >= 0 && (uint64_t)source_rows[i] > available - obs1))
            return STATA_ERR_INVALID_INPUT;
    }

    if (SF_var_is_string(var_idx)) {
        rc = validate_variable(var_idx, 1);
        if (rc) return rc;
        char **buf = (char **)ctools_safe_cacheline_alloc2(output_nobs, sizeof(char *));
        if (!buf) return STATA_ERR_MEMORY;
        memset(buf, 0, output_nobs * sizeof(char *));
        ctools_string_arena *arena = output_nobs <= SIZE_MAX / 64
            ? ctools_string_arena_create(output_nobs * 64, CTOOLS_STRING_ARENA_STRDUP_FALLBACK)
            : NULL;
        for (size_t i = 0; i < output_nobs; i++) {
            char strbuf[STATA_STR_MAXLEN + 1] = "";
            if (source_rows[i] >= 0 && SF_sdata(var_idx, (ST_int)(source_rows[i] + obs1), strbuf)) {
                rc = STATA_ERR_STATA_READ;
                break;
            }
            buf[i] = ctools_string_arena_strdup(arena, strbuf);
            if (!buf[i]) { rc = STATA_ERR_MEMORY; break; }
        }
        /* Complete the gather before any writes, including in-place permutations. */
        if (!rc) {
            for (size_t i = 0; i < output_nobs; i++) {
                if (SF_sstore(var_idx, (ST_int)(obs1 + i), buf[i])) {
                    rc = STATA_ERR_STATA_WRITE;
                    break;
                }
            }
        }
        for (size_t i = 0; i < output_nobs; i++)
            if (buf[i] && !ctools_string_arena_owns(arena, buf[i])) free(buf[i]);
        ctools_string_arena_free(arena);
        ctools_aligned_free(buf);
    } else {
        double *buf = (double *)ctools_safe_cacheline_alloc2(output_nobs, sizeof(double));
        if (!buf) return STATA_ERR_MEMORY;
        for (size_t i = 0; i < output_nobs; i++) {
            buf[i] = SV_missval;
            if (source_rows[i] >= 0 && (_stata_)->safevdata(var_idx, (ST_int)(source_rows[i] + obs1), &buf[i])) {
                rc = STATA_ERR_STATA_READ;
                break;
            }
        }
        if (!rc) {
            for (size_t i = 0; i < output_nobs; i++) {
                if ((_stata_)->safestore(var_idx, (ST_int)(obs1 + i), buf[i])) {
                    rc = STATA_ERR_STATA_WRITE;
                    break;
                }
            }
        }
        ctools_aligned_free(buf);
    }
    return rc == STATA_ERR_STATA_WRITE ? store_status(rc) : rc;
}

/* ===========================================================================
   Filtered Data Loading (Stata -> C with if/in filtering)
   =========================================================================== */

/* Forward declaration */
static int probe_for_identity(ST_int obs1, ST_int obs2, size_t n_range);

/*
    Build observation map from if/in filtering.

    Resolves the observation range, applies SF_ifobs filtering (unless
    CTOOLS_LOAD_SKIP_IF), and builds obs_map. This logic is shared by
    the column-parallel and row-parallel load paths.

    @param obs_start      [in]  First observation (1-based), 0 = use SF_in1()
    @param obs_end        [in]  Last observation (1-based), 0 = use SF_in2()
    @param flags          [in]  CTOOLS_LOAD_CHECK_IF or CTOOLS_LOAD_SKIP_IF
    @param obs_map_out    [out] Allocated obs_map (1-based Stata obs indices)
                                NULL if n_filtered == 0
    @param n_filtered_out [out] Number of observations passing filter
    @param n_range_out    [out] Total observations in range (before filtering)
    @param was_filtered_out [out] 1 if any obs excluded, 0 if identity

    @return STATA_OK on success, STATA_ERR_MEMORY on allocation failure
*/
static stata_retcode build_obs_map(
    size_t obs_start, size_t obs_end, int flags,
    perm_idx_t **obs_map_out, size_t *n_filtered_out,
    size_t *n_range_out, int *was_filtered_out)
{
    ST_int obs1, obs2;
    size_t n_range, n_filtered;
    perm_idx_t *obs_map = NULL;
    int is_identity = 0;
    int skip_if_check = (flags & CTOOLS_LOAD_SKIP_IF) != 0;

    if (obs_start > INT_MAX || obs_end > INT_MAX ||
        (obs_start && obs_end && obs_end < obs_start)) return STATA_ERR_INVALID_INPUT;

    /* Resolve observation range */
    if (obs_start == 0 || obs_end == 0) {
        obs1 = SF_in1();
        obs2 = SF_in2();
        if (obs1 < 1 || obs2 < obs1) {
            *obs_map_out = NULL;
            *n_filtered_out = 0;
            *n_range_out = 0;
            *was_filtered_out = 0;
            return STATA_OK;  /* Empty but valid */
        }
    } else {
        obs1 = (ST_int)obs_start;
        obs2 = (ST_int)obs_end;
    }

    if (validate_range((size_t)obs1, (size_t)(obs2 - obs1 + 1))) return STATA_ERR_INVALID_INPUT;
    n_range = (size_t)(obs2 - obs1 + 1);

    /* Fast path: skip SF_ifobs checks entirely */
    if (skip_if_check) {
        is_identity = 1;
        n_filtered = n_range;
        obs_map = (perm_idx_t *)ctools_safe_cacheline_alloc2(n_filtered, sizeof(perm_idx_t));
        if (obs_map == NULL) return STATA_ERR_MEMORY;
        for (size_t i = 0; i < n_filtered; i++) {
            obs_map[i] = (perm_idx_t)(obs1 + (ST_int)i);
        }
        *obs_map_out = obs_map;
        *n_filtered_out = n_filtered;
        *n_range_out = n_range;
        *was_filtered_out = 0;
        return STATA_OK;
    }

    /* Quick identity detection for small ranges */
    if (n_range <= 32) {
        is_identity = probe_for_identity(obs1, obs2, n_range);
        if (is_identity) {
            n_filtered = n_range;
            obs_map = (perm_idx_t *)ctools_safe_cacheline_alloc2(n_filtered, sizeof(perm_idx_t));
            if (obs_map == NULL) return STATA_ERR_MEMORY;
            for (size_t i = 0; i < n_filtered; i++) {
                obs_map[i] = (perm_idx_t)(obs1 + (ST_int)i);
            }
            *obs_map_out = obs_map;
            *n_filtered_out = n_filtered;
            *n_range_out = n_range;
            *was_filtered_out = 0;
            return STATA_OK;
        }
    }

    /* Standard path: two-pass SF_ifobs */
    /* Pass 1: count */
    n_filtered = 0;
    for (ST_int obs = obs1; obs <= obs2; obs++) {
        if (SF_ifobs(obs)) n_filtered++;
    }

    if (n_filtered == 0) {
        *obs_map_out = NULL;
        *n_filtered_out = 0;
        *n_range_out = n_range;
        *was_filtered_out = 1;
        return STATA_OK;
    }

    is_identity = (n_filtered == n_range);

    /* Pass 2: build obs_map */
    obs_map = (perm_idx_t *)ctools_safe_cacheline_alloc2(n_filtered, sizeof(perm_idx_t));
    if (obs_map == NULL) return STATA_ERR_MEMORY;

    if (is_identity) {
        for (size_t i = 0; i < n_filtered; i++) {
            obs_map[i] = (perm_idx_t)(obs1 + (ST_int)i);
        }
    } else {
        size_t idx = 0;
        for (ST_int obs = obs1; obs <= obs2; obs++) {
            if (SF_ifobs(obs)) obs_map[idx++] = (perm_idx_t)obs;
        }
    }

    *obs_map_out = obs_map;
    *n_filtered_out = n_filtered;
    *n_range_out = n_range;
    *was_filtered_out = !is_identity;
    return STATA_OK;
}

/*
    Row-parallel single-variable loader.

    Loads a single variable using OpenMP row-parallelism (splitting observations
    across threads). For single-variable loads, the column-parallel path is
    entirely sequential; this function provides parallelism instead.

    @param var        [out] Pre-allocated stata_variable to populate
    @param var_idx    [in]  1-based Stata variable index
    @param obs_map    [in]  Observation map (1-based Stata obs indices)
    @param n_filtered [in]  Number of observations
    @param str_width  [in]  Known string width (0 = numeric or unknown width)

    @return STATA_OK on success, a shared I/O error code on failure
*/
static int load_row_parallel(stata_variable *var, int var_idx,
                              perm_idx_t *obs_map, size_t n_filtered,
                              int str_width)
{
    int is_string = SF_var_is_string((ST_int)var_idx);

    var->nobs = n_filtered;

    /* Edge case: empty filtered set */
    if (n_filtered == 0) {
        if (is_string) {
            var->type = STATA_TYPE_STRING;
            var->str_maxlen = 0;
            var->_arena = NULL;
            var->data.str = (char **)calloc(1, sizeof(char *));
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
                char *flat_buf = (char *)calloc(n_filtered, stride);
                if (flat_buf != NULL) {
                    atomic_int load_error = 0;
                    #pragma omp parallel for schedule(static) if(n_filtered >= MIN_OBS_PER_THREAD * 2)
                    for (size_t i = 0; i < n_filtered; i++) {
                        if (atomic_load(&load_error)) continue;
                        char *slot = flat_buf + i * stride;
                        int read_rc = read_string_slot(var_idx, (ST_int)obs_map[i], slot, (size_t)str_width);
                        record_io_error(&load_error, read_rc);
                        var->data.str[i] = slot;
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
                size_t arena_capacity = n_filtered * 64;
                arena = ctools_string_arena_create(arena_capacity,
                                                   CTOOLS_STRING_ARENA_STRDUP_FALLBACK);
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

                ST_IIIS sdata_fn2 = read_stata_string;
                #pragma omp for schedule(static)
                for (size_t i = 0; i < n_filtered; i++) {
                    if (atomic_load(&load_error)) continue;
                    char strbuf[STATA_STR_MAXLEN + 1];
                    int read_rc = sdata_fn2((ST_int)var_idx, (ST_int)obs_map[i], strbuf);
                    if (read_rc) {
                        record_io_error(&load_error, read_rc);
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
                            record_io_error(&load_error, STATA_ERR_MEMORY);
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
        ST_IIIDp vdata_fn = (_stata_)->safevdata;
        atomic_int load_error = 0;
        #pragma omp parallel for schedule(static) if(n_filtered >= MIN_OBS_PER_THREAD * 2)
        for (size_t i = 0; i < n_filtered; i++) {
            if (vdata_fn((ST_int)var_idx, (ST_int)obs_map[i], &dbl_ptr[i])) record_io_error(&load_error, STATA_ERR_STATA_READ);
        }
        return atomic_load(&load_error);
    }

    return 0;
}

/*
    Initialize filtered data structure to safe empty state.
*/
void ctools_filtered_data_init(ctools_filtered_data *fd)
{
    if (fd == NULL) return;
    stata_data_init(&fd->data);
    fd->obs_map = NULL;
    fd->n_range = 0;
    fd->was_filtered = 0;
}

/*
    Free all memory associated with filtered data.
*/
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

/*
    Quick identity detection: probe first and last observations.
    Returns 1 if likely no filtering needed, 0 if filtering detected.
*/
static int probe_for_identity(ST_int obs1, ST_int obs2, size_t n_range)
{
    /* Probe size: 16 at each end, or full range if small */
    size_t probe_size = 16;
    if (n_range <= 32) {
        /* Small range: check all observations */
        for (ST_int obs = obs1; obs <= obs2; obs++) {
            if (!SF_ifobs(obs)) {
                return 0;  /* Filtering detected */
            }
        }
        return 1;  /* Identity confirmed */
    }

    /* Probe first 16 observations */
    for (size_t i = 0; i < probe_size; i++) {
        if (!SF_ifobs(obs1 + (ST_int)i)) {
            return 0;  /* Filtering detected */
        }
    }

    /* Probe last 16 observations */
    for (size_t i = 0; i < probe_size; i++) {
        if (!SF_ifobs(obs2 - (ST_int)i)) {
            return 0;  /* Filtering detected */
        }
    }

    return 1;  /* Likely identity (will verify during full pass if needed) */
}

/*
    Load a single variable from Stata for filtered observations only.
    Uses obs_map to read only the observations that passed filtering.
    When str_width > 0, uses flat buffer for string variables.
*/
static int load_filtered_variable(stata_variable *var, int var_idx,
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
            var->data.str = (char **)calloc(1, sizeof(char *));
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
        var->type = STATA_TYPE_STRING;

        /* ---- Flat buffer path: known string width ---- */
        if (str_width > 0 && str_width < STATA_STR_MAXLEN) {
            size_t stride = (size_t)str_width + 1;
            size_t flat_size = n_filtered * stride;

            if (flat_size / stride == n_filtered && flat_size <= (2048ULL * 1024 * 1024)) {
                char *flat_buf = (char *)calloc(n_filtered, stride);
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

                    for (i = 0; i < n_filtered; i++) {
                        str_ptrs[i] = flat_buf + i * stride;
                        int read_rc = read_string_slot(var_idx, (ST_int)obs_map[i], str_ptrs[i], (size_t)str_width);
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
            size_t arena_capacity = n_filtered * 64;
            arena = ctools_string_arena_create(arena_capacity, CTOOLS_STRING_ARENA_STRDUP_FALLBACK);
        }
        if (arena != NULL) {
            var->_arena = arena;
        }

        /* Load strings using obs_map */
        ST_IIIS sdata_fn3 = read_stata_string;
        for (i = 0; i < n_filtered; i++) {
            int read_rc = sdata_fn3((ST_int)var_idx, (ST_int)obs_map[i], strbuf);
            if (read_rc) return read_rc;
            str_ptr[i] = ctools_string_arena_strdup(arena, strbuf);
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
        ST_IIIDp vdata_fn = (_stata_)->safevdata;
        /* Load using obs_map - obs_map contains 1-based Stata obs numbers */
        for (i = 0; i < n_filtered; i++) {
            if (vdata_fn((ST_int)var_idx, (ST_int)obs_map[i], &dbl_ptr[i])) return STATA_ERR_STATA_READ;
        }
    }

    return 0;
}

/*
    Thread function for filtered variable loading.
*/
typedef struct {
    stata_variable *var;
    int var_idx;
    perm_idx_t *obs_map;
    size_t n_filtered;
    int is_string;
    int str_width;
    int error;
} filtered_var_io_args;

static void *load_filtered_variable_thread(void *arg)
{
    filtered_var_io_args *args = (filtered_var_io_args *)arg;
    args->error = load_filtered_variable(args->var, args->var_idx, args->obs_map,
                                         args->n_filtered, args->is_string, args->str_width);
    return args->error ? (void *)1 : NULL;
}

/*
    Extended data load with optional string width hints for flat buffer optimization.

    Auto-dispatches between:
    - nvars == 1: row-parallel path (OpenMP over observations)
    - nvars >= 2: column-parallel path (thread pool, 1 thread per variable)
*/
stata_retcode ctools_data_load_ex(ctools_filtered_data *result,
                                   int *var_indices, size_t nvars,
                                   size_t obs_start, size_t obs_end,
                                   int flags, const int *str_widths)
{
    size_t n_range, n_filtered;
    perm_idx_t *obs_map = NULL;
    int was_filtered = 0;
    int *auto_indices = NULL;

    if (result == NULL) {
        return STATA_ERR_INVALID_INPUT;
    }

    ctools_filtered_data_init(result);

    /* "Load all variables" mode */
    if (var_indices == NULL) {
        nvars = (size_t)SF_nvars();
        if (nvars == 0) {
            return STATA_OK;
        }
        auto_indices = (int *)malloc(nvars * sizeof(int));
        if (auto_indices == NULL) {
            return STATA_ERR_MEMORY;
        }
        for (size_t i = 0; i < nvars; i++) {
            auto_indices[i] = (int)(i + 1);
        }
        var_indices = auto_indices;
    } else if (nvars == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    for (size_t j = 0; j < nvars; j++) {
        stata_retcode valid = validate_variable(var_indices[j], -1);
        if (valid) { free(auto_indices); return valid; }
    }

    /* Auto-detect string widths from Stata local */
    int *auto_str_widths = NULL;
    if (str_widths == NULL) {
        auto_str_widths = ctools_read_strw_from_stata((size_t)SF_nvars());
        str_widths = auto_str_widths;
    }

    /* Build observation map (shared filtering logic) */
    stata_retcode rc = build_obs_map(obs_start, obs_end, flags,
                                      &obs_map, &n_filtered, &n_range,
                                      &was_filtered);
    if (rc != STATA_OK) {
        free(auto_str_widths);
        free(auto_indices);
        return load_status(rc);
    }

    result->n_range = n_range;
    result->obs_map = obs_map;
    result->was_filtered = was_filtered;

    /* Handle empty result */
    if (n_filtered == 0) {
        if (n_range > 0) {
            /* Some obs in range but all filtered out */
            rc = init_data_structure(&result->data, nvars, 0);
        }
        free(auto_str_widths);
        free(auto_indices);
        return load_status(rc);
    }

    /* Allocate data structure */
    rc = init_data_structure(&result->data, nvars, n_filtered);
    if (rc != STATA_OK) {
        ctools_filtered_data_free(result);
        free(auto_str_widths);
        free(auto_indices);
        return load_status(rc);
    }

    /*
        Auto-dispatch: row-parallel for single variable, column-parallel for multiple.
        For nvars==1, OpenMP splits observations across threads.
        For nvars>=2, thread pool assigns one thread per variable.
    */
    if (nvars == 1) {
        /* Row-parallel path */
        int str_width = str_widths ? str_widths[var_indices[0] - 1] : 0;
        rc = load_row_parallel(&result->data.vars[0], var_indices[0], obs_map, n_filtered, str_width);
        if (rc != 0) {
            ctools_filtered_data_free(result);
            free(auto_str_widths);
            free(auto_indices);
            return load_status(rc);
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
        init_io_thread_args(thread_args, &result->data, var_indices, nvars,
                            obs1, n_filtered, IO_MODE_LOAD, str_widths);

        rc = execute_io_parallel(thread_args, nvars, load_variable_thread);
        if (rc != 0) {
            free(thread_args);
            ctools_filtered_data_free(result);
            free(auto_str_widths);
            free(auto_indices);
            return load_status(rc);
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
            thread_args[j].is_string = SF_var_is_string((ST_int)var_indices[j]);
            thread_args[j].str_width = str_widths ? str_widths[var_indices[j] - 1] : 0;
            thread_args[j].error = STATA_OK;
        }

        ctools_persistent_pool *pool = ctools_get_global_pool();
        if (pool != NULL) {
            if (ctools_persistent_pool_submit_batch(pool, load_filtered_variable_thread,
                                                     thread_args, nvars,
                                                     sizeof(filtered_var_io_args)) != 0) {
                free(thread_args);
                ctools_filtered_data_free(result);
                free(auto_str_widths);
                free(auto_indices);
                return STATA_ERR_MEMORY;
            }
            if (ctools_persistent_pool_wait(pool)) rc = STATA_ERR_MEMORY;
        } else {
            for (size_t j = 0; j < nvars; j++) load_filtered_variable_thread(&thread_args[j]);
        }
        for (size_t j = 0; j < nvars; j++) {
            if (thread_args[j].error) { rc = thread_args[j].error; break; }
        }
        if (rc) {
            free(thread_args);
            ctools_filtered_data_free(result);
            free(auto_str_widths);
            free(auto_indices);
            return load_status(rc);
        }

        free(thread_args);
    }

    ctools_memory_barrier();

    free(auto_str_widths);
    free(auto_indices);
    return STATA_OK;
}

/*
    Backward-compatible wrapper: calls ctools_data_load_ex with no string width hints.
*/
stata_retcode ctools_data_load(ctools_filtered_data *result,
                                         int *var_indices, size_t nvars,
                                         size_t obs_start, size_t obs_end,
                                         int flags)
{
    return ctools_data_load_ex(result, var_indices, nvars, obs_start, obs_end,
                                flags, NULL);
}

/* ===========================================================================
   Filtered Data Storing (C -> Stata with obs_map)
   =========================================================================== */

/*
    Write filtered numeric values back to Stata using obs_map.
*/
stata_retcode ctools_store_filtered(double *values, size_t n_filtered,
                                     int var_idx, perm_idx_t *obs_map)
{
    if (!n_filtered) return STATA_OK;
    if (!values) return STATA_ERR_INVALID_INPUT;
    stata_retcode rc = validate_variable(var_idx, 0);
    if (rc) return rc;
    rc = validate_obs_map(obs_map, n_filtered);
    if (rc) return rc;
    for (size_t i = 0; i < n_filtered; i++) {
        if ((_stata_)->safestore(var_idx, (ST_int)obs_map[i], values[i])) return store_status(STATA_ERR_STATA_WRITE);
    }
    return STATA_OK;
}

/*
    Row-parallel store for a single numeric variable.
    Counterpart to ctools_data_load_single_var_rowpar().
    Validates the entire observation map before starting worker writes.
*/
stata_retcode ctools_store_filtered_rowpar(double *values, size_t n_filtered,
                                            int var_idx, perm_idx_t *obs_map)
{
    if (!n_filtered) return STATA_OK;
    if (!values) return STATA_ERR_INVALID_INPUT;
    stata_retcode rc = validate_variable(var_idx, 0);
    if (rc) return rc;
    rc = validate_obs_map(obs_map, n_filtered);
    if (rc) return rc;
    atomic_int error = 0;
    #pragma omp parallel for schedule(static) if(n_filtered >= MIN_OBS_PER_THREAD * 2)
    for (size_t i = 0; i < n_filtered; i++) {
        if ((_stata_)->safestore(var_idx, (ST_int)obs_map[i], values[i])) record_io_error(&error, STATA_ERR_STATA_WRITE);
    }
    return store_status(atomic_load(&error));
}

/*
    Write filtered string values back to Stata using obs_map.
*/
stata_retcode ctools_store_filtered_str(char **strings, size_t n_filtered,
                                         int var_idx, perm_idx_t *obs_map)
{
    if (!n_filtered) return STATA_OK;
    if (!strings) return STATA_ERR_INVALID_INPUT;
    stata_retcode rc = validate_variable(var_idx, 1);
    if (rc) return rc;
    rc = validate_obs_map(obs_map, n_filtered);
    if (rc) return rc;
    for (size_t i = 0; i < n_filtered; i++) {
        if (SF_sstore(var_idx, (ST_int)obs_map[i], strings[i] ? strings[i] : "")) return store_status(STATA_ERR_STATA_WRITE);
    }
    return STATA_OK;
}

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
    - SD_FASTMODE (compile flag) disables SPI bounds checking
    - Cache-line aligned allocations for optimal memory access
    - String arena allocator for reduced malloc overhead
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_threads.h"
#include "ctools_config.h"
#include "ctools_arena.h"

/* Maximum string buffer size for Stata string variables */
#define STATA_STR_MAXLEN 2045

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
    ST_int rc = SF_macro_use("___ctools_strw", buf, sizeof(buf));
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
            if (!var->data.str) return -1;
        } else {
            var->type = STATA_TYPE_DOUBLE;
            var->_arena = NULL;
            /* Allocate minimal buffer to avoid NULL */
            var->data.dbl = (double *)ctools_cacheline_alloc(sizeof(double));
            if (!var->data.dbl) return -1;
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
                        return -1;
                    }
                    char **str_ptrs = (char **)ctools_cacheline_alloc(str_array_size);
                    if (!str_ptrs) {
                        free(flat_buf);
                        return -1;
                    }

                    /* Read directly from Stata into flat buffer — single copy */
                    for (i = 0; i < nobs; i++) {
                        str_ptrs[i] = flat_buf + i * stride;
                        SF_sdata((ST_int)var_idx, (ST_int)(i + obs1), str_ptrs[i]);
                    }

                    /* Wrap flat buffer as a ctools_string_arena for compatible cleanup.
                       stata_data_free checks has_fallback=0 → takes O(1) bulk free path:
                       ctools_string_arena_free does free(base) + free(arena). */
                    ctools_string_arena *arena = (ctools_string_arena *)malloc(sizeof(ctools_string_arena));
                    if (!arena) {
                        ctools_aligned_free(str_ptrs);
                        free(flat_buf);
                        return -1;
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
            return -1;  /* Overflow */
        }
        var->data.str = (char **)ctools_cacheline_alloc(str_array_size);
        if (var->data.str == NULL) {
            return -1;
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
        for (i = 0; i < nobs; i++) {
            SF_sdata((ST_int)var_idx, (ST_int)(i + obs1), strbuf);
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
                return -1;
            }
        }
    } else {
        /* Numeric variable - use cache-line aligned allocation */
        var->type = STATA_TYPE_DOUBLE;
        var->_arena = NULL;
        size_t dbl_array_size;
        if (ctools_safe_mul_size(nobs, sizeof(double), &dbl_array_size) != 0) {
            return -1;  /* Overflow */
        }
        var->data.dbl = (double *)ctools_cacheline_alloc(dbl_array_size);

        if (var->data.dbl == NULL) {
            return -1;
        }

        double * restrict dbl_ptr = var->data.dbl;

        for (i = 0; i < nobs; i++) {
            SF_vdata((ST_int)var_idx, (ST_int)(i + obs1), &dbl_ptr[i]);
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
    args->success = 0;

    if (load_single_variable(args->var, args->var_idx, args->obs1,
                              args->nobs, args->is_string,
                              args->str_width) != 0) {
        return (void *)1;
    }

    args->success = 1;
    return NULL;
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
        args[j].success = 0;
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
    @param check_error  [in]  Non-zero to check for errors in load mode

    @return 0 on success, non-zero on failure
*/
static int execute_io_parallel(ctools_var_io_args *args, size_t nvars,
                               io_thread_func func, int check_error)
{
    int use_parallel = (nvars >= 2);

    if (use_parallel) {
        ctools_persistent_pool *pool = ctools_get_global_pool();

        if (pool != NULL) {
            /* Submit batch to thread pool */
            if (ctools_persistent_pool_submit_batch(pool, func, args, nvars,
                                                     sizeof(ctools_var_io_args)) != 0) {
                return -1;
            }

            int pool_result = ctools_persistent_pool_wait(pool);
            if (check_error && pool_result != 0) {
                return -1;
            }
        } else {
            /* Fallback to sequential */
            for (size_t j = 0; j < nvars; j++) {
                void *result = func(&args[j]);
                if (check_error && result != NULL) {
                    return -1;
                }
            }
        }
    } else {
        /* Sequential execution */
        for (size_t j = 0; j < nvars; j++) {
            void *result = func(&args[j]);
            if (check_error && result != NULL) {
                return -1;
            }
        }
    }

    return 0;
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
static void store_single_variable(stata_variable *var, int var_idx,
                                  size_t obs1, size_t nobs)
{
    size_t i;
    ST_int stata_var_idx = (ST_int)var_idx;

    if (var->type == STATA_TYPE_DOUBLE) {
        const double * restrict dbl_data = var->data.dbl;

        for (i = 0; i < nobs; i++) {
            SF_vstore(stata_var_idx, (ST_int)(i + obs1), dbl_data[i]);
        }

    } else {
        /* String variable */
        char * const * restrict str_data = var->data.str;
        size_t str_width = var->str_maxlen;

        /* Use flat buffer when we have actual variable width (not default 2045)
           and the buffer fits in 256 MB */
        if (str_width > 0 && str_width < STATA_STR_MAXLEN) {
            size_t stride = str_width + 1;
            ctools_string_arena *arena = (ctools_string_arena *)var->_arena;

            /* Fast path: sequential flat buffer (after flat-buffer permutation).
               Detect by checking if first pointer matches arena base — if so,
               the data is already contiguous and we can SF_sstore directly
               without any repacking. */
            if (arena != NULL && !arena->has_fallback &&
                nobs > 0 && str_data[0] == arena->base) {
                char *base = arena->base;
                for (i = 0; i < nobs; i++) {
                    SF_sstore(stata_var_idx, (ST_int)(i + obs1),
                              base + i * stride);
                }
                return;
            }

            /* Repack path: strings are scattered, copy into new flat buffer */
            size_t flat_size = nobs * stride;
            char *flat_buf = NULL;
            if (flat_size / stride == nobs &&  /* overflow check */
                flat_size <= (256ULL * 1024 * 1024)) {
                flat_buf = (char *)calloc(nobs, stride);
            }

            if (flat_buf) {
                /* Use stride memcpy when source is flat buffer (no strlen needed) */
                if (arena != NULL && !arena->has_fallback) {
                    for (i = 0; i < nobs; i++) {
                        const char *s = str_data[i];
                        if (s) memcpy(flat_buf + i * stride, s, stride);
                    }
                } else {
                    /* Non-flat source: strlen + bounded memcpy */
                    for (i = 0; i < nobs; i++) {
                        const char *s = str_data[i];
                        if (s) {
                            size_t len = strlen(s);
                            if (len >= stride) len = stride - 1;
                            memcpy(flat_buf + i * stride, s, len);
                        }
                    }
                }

                /* Store: sequential scan through flat buffer */
                for (i = 0; i < nobs; i++) {
                    SF_sstore(stata_var_idx, (ST_int)(i + obs1),
                              flat_buf + i * stride);
                }

                free(flat_buf);
                return;
            }
        }

        /* Fallback: strL, width unknown, buffer too large, or alloc failed */
        for (i = 0; i < nobs; i++) {
            SF_sstore(stata_var_idx, (ST_int)(i + obs1), str_data[i]);
        }
    }
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
    store_single_variable(args->var, args->var_idx, args->obs1, args->nobs);
    args->success = 1;
    return NULL;
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

        if (var->type == STATA_TYPE_DOUBLE) {
            const double * restrict dbl_data = var->data.dbl;
            #pragma omp parallel for schedule(static) if(nobs >= MIN_OBS_PER_THREAD * 2)
            for (size_t i = 0; i < nobs; i++) {
                SF_vstore(stata_var, (ST_int)(i + obs1), dbl_data[i]);
            }
        } else {
            /* String store: row-parallel when possible */
            char * const * restrict str_data = var->data.str;
            size_t str_width = var->str_maxlen;
            ctools_string_arena *arena = (ctools_string_arena *)var->_arena;

            /* Fast path: sequential flat buffer — row-parallel SF_sstore
               without any repacking (data is already contiguous) */
            if (arena != NULL && !arena->has_fallback &&
                str_width > 0 && str_width < STATA_STR_MAXLEN &&
                nobs > 0 && str_data[0] == arena->base) {
                size_t stride = str_width + 1;
                char *base = arena->base;
                #pragma omp parallel for schedule(static) if(nobs >= MIN_OBS_PER_THREAD * 2)
                for (size_t i = 0; i < nobs; i++) {
                    SF_sstore(stata_var, (ST_int)(i + obs1),
                              base + i * stride);
                }
            } else {
                /* Non-sequential or non-flat: row-parallel from pointer array */
                #pragma omp parallel for schedule(static) if(nobs >= MIN_OBS_PER_THREAD * 2)
                for (size_t i = 0; i < nobs; i++) {
                    SF_sstore(stata_var, (ST_int)(i + obs1), str_data[i]);
                }
            }
        }

        ctools_memory_barrier();
        return STATA_OK;
    }

    /* nvars >= 2: column-parallel store via thread pool */
    ctools_var_io_args *thread_args = (ctools_var_io_args *)
        ctools_safe_malloc2(nvars, sizeof(ctools_var_io_args));
    if (thread_args == NULL) {
        return STATA_ERR_MEMORY;
    }

    init_io_thread_args(thread_args, data, var_indices, nvars, obs1, nobs,
                        IO_MODE_STORE, NULL);
    execute_io_parallel(thread_args, nvars, store_variable_thread, 0);
    free(thread_args);

    ctools_memory_barrier();
    return STATA_OK;
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
    size_t i;
    ST_int stata_var = (ST_int)var_idx;
    int is_string = SF_var_is_string(stata_var);
    ST_int nobs_max = SF_nobs();

    if (is_string) {
        /* String variable: allocate buffer, gather, scatter (overflow-safe) */
        char **buf = (char **)ctools_safe_cacheline_alloc2(output_nobs, sizeof(char *));
        if (!buf) return STATA_ERR_MEMORY;
        memset(buf, 0, output_nobs * sizeof(char *));  /* Zero-init for safe cleanup */

        char strbuf[2048];

        /* Create arena for string storage (estimate 64 bytes avg per string) */
        /* Overflow check: output_nobs * 64 - skip arena if overflow, rely on strdup fallback */
        ctools_string_arena *arena = NULL;
        if (output_nobs <= SIZE_MAX / 64) {
            size_t arena_capacity = output_nobs * 64;
            arena = ctools_string_arena_create(arena_capacity, CTOOLS_STRING_ARENA_STRDUP_FALLBACK);
        }
        /* Arena failure is ok - we fall back to strdup */

        /* GATHER: Read from source positions */
        for (i = 0; i < output_nobs; i++) {
            if (source_rows[i] >= 0) {
                /* Overflow check: ensure addition doesn't overflow int64_t */
                if (source_rows[i] > INT64_MAX - (int64_t)obs1) {
                    char errbuf[128];
                    snprintf(errbuf, sizeof(errbuf),
                             "ctools: source_rows[%zu]=%lld + obs1=%zu would overflow\n",
                             i, (long long)source_rows[i], obs1);
                    SF_error(errbuf);
                    /* Cleanup */
                    for (size_t j = 0; j < i; j++) {
                        if (buf[j] != NULL && !ctools_string_arena_owns(arena, buf[j])) {
                            free(buf[j]);
                        }
                    }
                    ctools_string_arena_free(arena);
                    ctools_aligned_free(buf);
                    return STATA_ERR_INVALID_INPUT;
                }
                /* Bounds check: source observation must be within valid Stata range */
                int64_t src_obs = source_rows[i] + (int64_t)obs1;
                if (src_obs < 1 || src_obs > nobs_max) {
                    char errbuf[128];
                    snprintf(errbuf, sizeof(errbuf),
                             "ctools: source_rows[%zu]=%lld -> obs %lld out of bounds [1,%d]\n",
                             i, (long long)source_rows[i], (long long)src_obs, (int)nobs_max);
                    SF_error(errbuf);
                    /* Cleanup */
                    for (size_t j = 0; j < i; j++) {
                        if (buf[j] != NULL && !ctools_string_arena_owns(arena, buf[j])) {
                            free(buf[j]);
                        }
                    }
                    ctools_string_arena_free(arena);
                    ctools_aligned_free(buf);
                    return STATA_ERR_INVALID_INPUT;
                }
                SF_sdata(stata_var, (ST_int)src_obs, strbuf);
                buf[i] = ctools_string_arena_strdup(arena, strbuf);
            } else {
                buf[i] = ctools_string_arena_strdup(arena, "");  /* Missing -> empty string */
            }
            if (!buf[i]) {
                /* Cleanup on allocation failure - free any fallback strings first */
                for (size_t j = 0; j < i; j++) {
                    if (buf[j] != NULL && !ctools_string_arena_owns(arena, buf[j])) {
                        free(buf[j]);
                    }
                }
                ctools_string_arena_free(arena);
                ctools_aligned_free(buf);
                return STATA_ERR_MEMORY;
            }
        }

        /* SCATTER: Write to sequential output positions */
        for (i = 0; i < output_nobs; i++) {
            SF_sstore(stata_var, (ST_int)(i + obs1), buf[i]);
        }

        /* Free fallback strings (those not owned by arena), then free arena */
        for (i = 0; i < output_nobs; i++) {
            if (buf[i] != NULL && !ctools_string_arena_owns(arena, buf[i])) {
                free(buf[i]);
            }
        }
        ctools_string_arena_free(arena);
        ctools_aligned_free(buf);

    } else {
        /* Numeric variable: allocate aligned buffer, gather, scatter (overflow-safe) */
        double *buf = (double *)ctools_safe_cacheline_alloc2(output_nobs, sizeof(double));
        if (!buf) return STATA_ERR_MEMORY;

        /* GATHER: Read from source positions */
        for (i = 0; i < output_nobs; i++) {
            if (source_rows[i] >= 0) {
                /* Overflow check: ensure addition doesn't overflow int64_t */
                if (source_rows[i] > INT64_MAX - (int64_t)obs1) {
                    char errbuf[128];
                    snprintf(errbuf, sizeof(errbuf),
                             "ctools: source_rows[%zu]=%lld + obs1=%zu would overflow\n",
                             i, (long long)source_rows[i], obs1);
                    SF_error(errbuf);
                    ctools_aligned_free(buf);
                    return STATA_ERR_INVALID_INPUT;
                }
                /* Bounds check: source observation must be within valid Stata range */
                int64_t src_obs = source_rows[i] + (int64_t)obs1;
                if (src_obs < 1 || src_obs > nobs_max) {
                    char errbuf[128];
                    snprintf(errbuf, sizeof(errbuf),
                             "ctools: source_rows[%zu]=%lld -> obs %lld out of bounds [1,%d]\n",
                             i, (long long)source_rows[i], (long long)src_obs, (int)nobs_max);
                    SF_error(errbuf);
                    ctools_aligned_free(buf);
                    return STATA_ERR_INVALID_INPUT;
                }
                SF_vdata(stata_var, (ST_int)src_obs, &buf[i]);
            } else {
                buf[i] = SV_missval;  /* Missing value */
            }
        }

        /* SCATTER: Write to sequential output positions */
        for (i = 0; i < output_nobs; i++) {
            SF_vstore(stata_var, (ST_int)(i + obs1), buf[i]);
        }

        ctools_aligned_free(buf);
    }

    return STATA_OK;
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

    @return 0 on success, -1 on failure
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
            if (!var->data.dbl) return -1;
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
            return -1;
        }
        var->data.str = (char **)ctools_cacheline_alloc(str_array_size);
        if (var->data.str == NULL) return -1;
        memset(var->data.str, 0, str_array_size);

        /* Fast path: flat buffer when string width is known */
        int used_flat_buffer = 0;
        if (str_width > 0 && str_width < STATA_STR_MAXLEN) {
            size_t stride = (size_t)str_width + 1;
            size_t flat_size;
            if (ctools_safe_mul_size(n_filtered, stride, &flat_size) == 0) {
                char *flat_buf = (char *)calloc(n_filtered, stride);
                if (flat_buf != NULL) {
                    volatile int load_error = 0;
                    #pragma omp parallel for schedule(static) if(n_filtered >= MIN_OBS_PER_THREAD * 2)
                    for (size_t i = 0; i < n_filtered; i++) {
                        if (load_error) continue;
                        char *slot = flat_buf + i * stride;
                        if (SF_sdata((ST_int)var_idx, (ST_int)obs_map[i], slot) != 0) {
                            load_error = 1;
                        }
                        var->data.str[i] = slot;
                    }

                    if (load_error) {
                        free(flat_buf);
                        ctools_aligned_free(var->data.str);
                        var->data.str = NULL;
                        return -1;
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
                        return -1;
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

            volatile int load_error = 0;

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

                #pragma omp for schedule(static)
                for (size_t i = 0; i < n_filtered; i++) {
                    if (load_error) continue;
                    char strbuf[STATA_STR_MAXLEN + 1];
                    SF_sdata((ST_int)var_idx, (ST_int)obs_map[i], strbuf);
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
                            load_error = 1;
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

            if (load_error) {
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
                return -1;
            }
        }
    } else {
        /* Numeric variable */
        var->type = STATA_TYPE_DOUBLE;
        var->_arena = NULL;

        size_t dbl_array_size;
        if (ctools_safe_mul_size(n_filtered, sizeof(double), &dbl_array_size) != 0) {
            return -1;
        }
        var->data.dbl = (double *)ctools_cacheline_alloc(dbl_array_size);
        if (var->data.dbl == NULL) return -1;

        double * restrict dbl_ptr = var->data.dbl;

        #pragma omp parallel for schedule(static) if(n_filtered >= MIN_OBS_PER_THREAD * 2)
        for (size_t i = 0; i < n_filtered; i++) {
            SF_vdata((ST_int)var_idx, (ST_int)obs_map[i], &dbl_ptr[i]);
        }
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
            if (!var->data.str) return -1;
        } else {
            var->type = STATA_TYPE_DOUBLE;
            var->_arena = NULL;
            var->data.dbl = (double *)ctools_cacheline_alloc(sizeof(double));
            if (!var->data.dbl) return -1;
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
                        return -1;
                    }
                    char **str_ptrs = (char **)ctools_cacheline_alloc(str_array_size);
                    if (!str_ptrs) {
                        free(flat_buf);
                        return -1;
                    }

                    for (i = 0; i < n_filtered; i++) {
                        str_ptrs[i] = flat_buf + i * stride;
                        SF_sdata((ST_int)var_idx, (ST_int)obs_map[i], str_ptrs[i]);
                    }

                    ctools_string_arena *arena = (ctools_string_arena *)malloc(sizeof(ctools_string_arena));
                    if (!arena) {
                        ctools_aligned_free(str_ptrs);
                        free(flat_buf);
                        return -1;
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
            return -1;
        }
        var->data.str = (char **)ctools_cacheline_alloc(str_array_size);
        if (var->data.str == NULL) {
            return -1;
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
        for (i = 0; i < n_filtered; i++) {
            SF_sdata((ST_int)var_idx, (ST_int)obs_map[i], strbuf);
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
                return -1;
            }
        }
    } else {
        /* Numeric variable */
        var->type = STATA_TYPE_DOUBLE;
        var->_arena = NULL;

        size_t dbl_array_size;
        if (ctools_safe_mul_size(n_filtered, sizeof(double), &dbl_array_size) != 0) {
            return -1;
        }
        var->data.dbl = (double *)ctools_cacheline_alloc(dbl_array_size);
        if (var->data.dbl == NULL) {
            return -1;
        }

        double * restrict dbl_ptr = var->data.dbl;

        /* Load using obs_map - obs_map contains 1-based Stata obs numbers */
        for (i = 0; i < n_filtered; i++) {
            SF_vdata((ST_int)var_idx, (ST_int)obs_map[i], &dbl_ptr[i]);
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
    int success;
} filtered_var_io_args;

static void *load_filtered_variable_thread(void *arg)
{
    filtered_var_io_args *args = (filtered_var_io_args *)arg;
    args->success = 0;

    if (load_filtered_variable(args->var, args->var_idx, args->obs_map,
                                args->n_filtered, args->is_string,
                                args->str_width) != 0) {
        return (void *)1;
    }

    args->success = 1;
    return NULL;
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
        return rc;
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
        return rc;
    }

    /* Allocate data structure */
    rc = init_data_structure(&result->data, nvars, n_filtered);
    if (rc != STATA_OK) {
        ctools_filtered_data_free(result);
        free(auto_str_widths);
        free(auto_indices);
        return rc;
    }

    /*
        Auto-dispatch: row-parallel for single variable, column-parallel for multiple.
        For nvars==1, OpenMP splits observations across threads.
        For nvars>=2, thread pool assigns one thread per variable.
    */
    if (nvars == 1) {
        /* Row-parallel path */
        int str_width = str_widths ? str_widths[var_indices[0] - 1] : 0;
        if (load_row_parallel(&result->data.vars[0], var_indices[0],
                               obs_map, n_filtered, str_width) != 0) {
            ctools_filtered_data_free(result);
            free(auto_str_widths);
            free(auto_indices);
            return STATA_ERR_MEMORY;
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

        if (execute_io_parallel(thread_args, nvars, load_variable_thread, 1) != 0) {
            free(thread_args);
            ctools_filtered_data_free(result);
            free(auto_str_widths);
            free(auto_indices);
            return STATA_ERR_MEMORY;
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
            thread_args[j].success = 0;
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
            int pool_result = ctools_persistent_pool_wait(pool);
            if (pool_result != 0) {
                free(thread_args);
                ctools_filtered_data_free(result);
                free(auto_str_widths);
                free(auto_indices);
                return STATA_ERR_MEMORY;
            }
        } else {
            for (size_t j = 0; j < nvars; j++) {
                void *res = load_filtered_variable_thread(&thread_args[j]);
                if (res != NULL) {
                    free(thread_args);
                    ctools_filtered_data_free(result);
                    free(auto_str_widths);
                    free(auto_indices);
                    return STATA_ERR_MEMORY;
                }
            }
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
    if (values == NULL || obs_map == NULL || n_filtered == 0) {
        return (values == NULL && n_filtered == 0) ? STATA_OK : STATA_ERR_INVALID_INPUT;
    }

    ST_int stata_var = (ST_int)var_idx;
    ST_int nobs_max = SF_nobs();

    /* Write values using obs_map for indexing */
    for (size_t i = 0; i < n_filtered; i++) {
        /* Bounds check: obs_map values must be valid 1-based Stata observation indices */
        perm_idx_t obs = obs_map[i];
        if (obs < 1 || obs > (perm_idx_t)nobs_max) {
            char buf[128];
            snprintf(buf, sizeof(buf),
                     "ctools: obs_map[%zu]=%u out of bounds [1,%d]\n",
                     i, (unsigned)obs, (int)nobs_max);
            SF_error(buf);
            return STATA_ERR_INVALID_INPUT;
        }
        SF_vstore(stata_var, (ST_int)obs, values[i]);
    }

    return STATA_OK;
}

/*
    Row-parallel store for a single numeric variable.
    Counterpart to ctools_data_load_single_var_rowpar().
    Skips per-element bounds checks (obs_map is trusted from our own loader).
*/
stata_retcode ctools_store_filtered_rowpar(double *values, size_t n_filtered,
                                            int var_idx, perm_idx_t *obs_map)
{
    if (values == NULL || obs_map == NULL || n_filtered == 0) {
        return (values == NULL && n_filtered == 0) ? STATA_OK : STATA_ERR_INVALID_INPUT;
    }

    ST_int stata_var = (ST_int)var_idx;

    #pragma omp parallel for schedule(static) if(n_filtered >= MIN_OBS_PER_THREAD * 2)
    for (size_t i = 0; i < n_filtered; i++) {
        SF_vstore(stata_var, (ST_int)obs_map[i], values[i]);
    }

    return STATA_OK;
}

/*
    Write filtered string values back to Stata using obs_map.
*/
stata_retcode ctools_store_filtered_str(char **strings, size_t n_filtered,
                                         int var_idx, perm_idx_t *obs_map)
{
    if (strings == NULL || obs_map == NULL || n_filtered == 0) {
        return (strings == NULL && n_filtered == 0) ? STATA_OK : STATA_ERR_INVALID_INPUT;
    }

    ST_int stata_var = (ST_int)var_idx;
    ST_int nobs_max = SF_nobs();

    /* Write strings using obs_map for indexing */
    for (size_t i = 0; i < n_filtered; i++) {
        /* Bounds check: obs_map values must be valid 1-based Stata observation indices */
        perm_idx_t obs = obs_map[i];
        if (obs < 1 || obs > (perm_idx_t)nobs_max) {
            char buf[128];
            snprintf(buf, sizeof(buf),
                     "ctools: obs_map[%zu]=%u out of bounds [1,%d]\n",
                     i, (unsigned)obs, (int)nobs_max);
            SF_error(buf);
            return STATA_ERR_INVALID_INPUT;
        }
        SF_sstore(stata_var, (ST_int)obs, strings[i] ? strings[i] : "");
    }

    return STATA_OK;
}

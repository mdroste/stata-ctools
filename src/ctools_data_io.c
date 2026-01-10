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

    Performance Optimizations:
    - 8x loop unrolling reduces loop overhead and improves instruction pipelining
    - Parallel variable I/O overlaps operations across columns
    - SD_FASTMODE (compile flag) disables SPI bounds checking

    Thread Safety:
    - Each thread operates on a different Stata variable (column)
    - No shared state during I/O phase
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_threads.h"

/* Maximum string buffer size for Stata string variables */
#define STATA_STR_MAXLEN 2045

/* ===========================================================================
   Data Loading (Stata -> C)
   =========================================================================== */

/*
    Thread function: Load a single variable from Stata into C memory.

    Allocates memory and reads all observations for one variable using
    SF_vdata (numeric) or SF_sdata (string). Uses 8x unrolled loop for
    numeric data to improve throughput.

    @param arg  Pointer to ctools_var_io_args with input/output parameters
    @return     NULL on success, non-NULL on failure (for ctools_threads)
*/
static void *load_variable_thread(void *arg)
{
    ctools_var_io_args *args = (ctools_var_io_args *)arg;
    size_t i;
    size_t nobs = args->nobs;
    size_t obs1 = args->obs1;
    ST_int var_idx = (ST_int)args->var_idx;
    ST_double val0, val1, val2, val3, val4, val5, val6, val7;
    char strbuf[STATA_STR_MAXLEN + 1];

    args->var->nobs = nobs;
    args->success = 0;  /* Assume failure until complete */

    if (args->is_string) {
        /* String variable */
        args->var->type = STATA_TYPE_STRING;
        args->var->str_maxlen = STATA_STR_MAXLEN;
        args->var->data.str = (char **)calloc(nobs, sizeof(char *));

        if (args->var->data.str == NULL) {
            return (void *)1;  /* Signal failure */
        }

        char **str_ptr = args->var->data.str;
        for (i = 0; i < nobs; i++) {
            SF_sdata(var_idx, (ST_int)(i + obs1), strbuf);
            str_ptr[i] = strdup(strbuf);
            if (str_ptr[i] == NULL) {
                return (void *)1;  /* Signal failure */
            }
        }
    } else {
        /* Numeric variable */
        args->var->type = STATA_TYPE_DOUBLE;
        args->var->data.dbl = (double *)malloc(nobs * sizeof(double));

        if (args->var->data.dbl == NULL) {
            return (void *)1;  /* Signal failure */
        }

        double *dbl_ptr = args->var->data.dbl;
        size_t i_end = nobs - (nobs % 8);

        /* Unrolled loop - process 8 at a time */
        for (i = 0; i < i_end; i += 8) {
            SF_vdata(var_idx, (ST_int)(i + obs1), &val0);
            SF_vdata(var_idx, (ST_int)(i + 1 + obs1), &val1);
            SF_vdata(var_idx, (ST_int)(i + 2 + obs1), &val2);
            SF_vdata(var_idx, (ST_int)(i + 3 + obs1), &val3);
            SF_vdata(var_idx, (ST_int)(i + 4 + obs1), &val4);
            SF_vdata(var_idx, (ST_int)(i + 5 + obs1), &val5);
            SF_vdata(var_idx, (ST_int)(i + 6 + obs1), &val6);
            SF_vdata(var_idx, (ST_int)(i + 7 + obs1), &val7);

            dbl_ptr[i] = val0;
            dbl_ptr[i + 1] = val1;
            dbl_ptr[i + 2] = val2;
            dbl_ptr[i + 3] = val3;
            dbl_ptr[i + 4] = val4;
            dbl_ptr[i + 5] = val5;
            dbl_ptr[i + 6] = val6;
            dbl_ptr[i + 7] = val7;
        }
        /* Handle remaining elements */
        for (; i < nobs; i++) {
            SF_vdata(var_idx, (ST_int)(i + obs1), &val0);
            dbl_ptr[i] = val0;
        }
    }

    args->success = 1;
    return NULL;  /* Success */
}

/*
    Load all variables from Stata's data space into C memory.

    This is the core data transfer function - it performs pure data loading
    without any sort-specific setup. Use this for general-purpose data transfer
    between Stata and C.

    @param data   [out] stata_data structure to populate (caller frees)
    @param nvars  [in]  Total number of variables in dataset

    @return STATA_OK on success, or:
            STATA_ERR_INVALID_INPUT if data is NULL or nvars is 0
            STATA_ERR_MEMORY on allocation failure

    Memory Layout After Return:
    - data->vars[j].data.dbl: contiguous array of nobs doubles (numeric)
    - data->vars[j].data.str: array of nobs char* pointers (string)
    - data->sort_order: identity permutation [0, 1, 2, ..., nobs-1]
*/
stata_retcode ctools_data_load(stata_data *data, size_t nvars)
{
    size_t i;
    size_t nobs;
    size_t obs1;
    ctools_var_io_args *thread_args = NULL;
    size_t j;
    int use_parallel;

    if (data == NULL || nvars == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Cache observation bounds */
    obs1 = (size_t)SF_in1();
    nobs = (size_t)(SF_in2() - SF_in1() + 1);

    /* Initialize the data structure */
    stata_data_init(data);
    data->nobs = (size_t)nobs;
    data->nvars = nvars;

    /* Allocate array of variables */
    data->vars = (stata_variable *)calloc(nvars, sizeof(stata_variable));
    if (data->vars == NULL) {
        return STATA_ERR_MEMORY;
    }

    /* Allocate sort order array */
    data->sort_order = (size_t *)malloc(nobs * sizeof(size_t));
    if (data->sort_order == NULL) {
        stata_data_free(data);
        return STATA_ERR_MEMORY;
    }

    /* Initialize sort order to identity permutation - unrolled */
    size_t idx_end = (size_t)(nobs - (nobs % 8));
    for (i = 0; (size_t)i < idx_end; i += 8) {
        data->sort_order[i] = (size_t)i;
        data->sort_order[i + 1] = (size_t)(i + 1);
        data->sort_order[i + 2] = (size_t)(i + 2);
        data->sort_order[i + 3] = (size_t)(i + 3);
        data->sort_order[i + 4] = (size_t)(i + 4);
        data->sort_order[i + 5] = (size_t)(i + 5);
        data->sort_order[i + 6] = (size_t)(i + 6);
        data->sort_order[i + 7] = (size_t)(i + 7);
    }
    for (; i < nobs; i++) {
        data->sort_order[i] = (size_t)i;
    }

    /* Allocate thread arguments (used for both parallel and sequential) */
    thread_args = (ctools_var_io_args *)malloc(nvars * sizeof(ctools_var_io_args));
    if (thread_args == NULL) {
        stata_data_free(data);
        return STATA_ERR_MEMORY;
    }

    /* Initialize thread arguments */
    for (j = 0; j < nvars; j++) {
        thread_args[j].var = &data->vars[j];
        thread_args[j].var_idx = (int)(j + 1);
        thread_args[j].obs1 = obs1;
        thread_args[j].nobs = nobs;
        thread_args[j].is_string = SF_var_is_string((ST_int)(j + 1));
        thread_args[j].success = 0;
    }

    /* Decide whether to use parallel loading */
    use_parallel = (nvars >= 2);

    if (use_parallel) {
        /* Parallel loading using thread pool */
        ctools_thread_pool pool;
        if (ctools_pool_init(&pool, nvars, thread_args, sizeof(ctools_var_io_args)) != 0) {
            free(thread_args);
            stata_data_free(data);
            return STATA_ERR_MEMORY;
        }

        int pool_result = ctools_pool_run(&pool, load_variable_thread);
        ctools_pool_free(&pool);
        free(thread_args);

        if (pool_result != 0) {
            stata_data_free(data);
            return STATA_ERR_MEMORY;
        }
    } else {
        /* Sequential loading - reuse thread function directly */
        for (j = 0; j < nvars; j++) {
            if (load_variable_thread(&thread_args[j]) != NULL) {
                free(thread_args);
                stata_data_free(data);
                return STATA_ERR_MEMORY;
            }
        }
        free(thread_args);
    }

    return STATA_OK;
}

/* ===========================================================================
   Data Storing (C -> Stata)
   =========================================================================== */

/*
    Thread function: Store a single variable from C memory to Stata.

    Writes all observations for one variable using SF_vstore (numeric) or
    SF_sstore (string). Uses 8x unrolled loop for numeric data. Data is
    read sequentially from C memory (assumes permutation already applied).

    @param arg  Pointer to ctools_var_io_args with input/output parameters
    @return     NULL on success (for ctools_threads compatibility)
*/
static void *store_variable_thread(void *arg)
{
    ctools_var_io_args *args = (ctools_var_io_args *)arg;
    size_t i;
    size_t nobs = args->nobs;
    size_t obs1 = args->obs1;
    ST_int var_idx = (ST_int)args->var_idx;

    if (args->var->type == STATA_TYPE_DOUBLE) {
        /* Numeric variable - sequential store */
        double *dbl_data = args->var->data.dbl;
        size_t i_end = nobs - (nobs % 8);

        /* Main loop with 8x unrolling */
        for (i = 0; i < i_end; i += 8) {
            SF_vstore(var_idx, (ST_int)(i + obs1), dbl_data[i]);
            SF_vstore(var_idx, (ST_int)(i + 1 + obs1), dbl_data[i + 1]);
            SF_vstore(var_idx, (ST_int)(i + 2 + obs1), dbl_data[i + 2]);
            SF_vstore(var_idx, (ST_int)(i + 3 + obs1), dbl_data[i + 3]);
            SF_vstore(var_idx, (ST_int)(i + 4 + obs1), dbl_data[i + 4]);
            SF_vstore(var_idx, (ST_int)(i + 5 + obs1), dbl_data[i + 5]);
            SF_vstore(var_idx, (ST_int)(i + 6 + obs1), dbl_data[i + 6]);
            SF_vstore(var_idx, (ST_int)(i + 7 + obs1), dbl_data[i + 7]);
        }
        /* Handle remaining elements */
        for (; i < nobs; i++) {
            SF_vstore(var_idx, (ST_int)(i + obs1), dbl_data[i]);
        }
    } else {
        /* String variable - sequential store */
        char **str_data = args->var->data.str;
        for (i = 0; i < nobs; i++) {
            SF_sstore(var_idx, (ST_int)(i + obs1), str_data[i]);
        }
    }

    args->success = 1;
    return NULL;  /* Success */
}

/*
    Write all variables from C memory back to Stata's data space.

    This is the main entry point for the C-to-Stata data transfer module.
    It writes all variables back to Stata in their current order. The data
    must already be in final sorted order (permutation applied by the sort
    module).

    @param data  [in] stata_data structure with data to write; data arrays
                      must be in final order (sequential read, no permutation)
    @param obs1  [in] First observation in Stata (1-based, from SF_in1)

    @return STATA_OK on success, or:
            STATA_ERR_INVALID_INPUT if data is NULL
            STATA_ERR_MEMORY on allocation failure (parallel mode only)

    Parallelization:
    - Uses parallel writes when nvars >= 2 (one thread per variable)
    - Falls back to sequential writes for single-variable datasets
*/
stata_retcode ctools_data_store(stata_data *data, size_t obs1)
{
    size_t j;
    size_t nobs, nvars;
    ctools_var_io_args *thread_args = NULL;
    int use_parallel;

    if (data == NULL) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Cache frequently accessed values */
    nobs = data->nobs;
    nvars = data->nvars;

    /* Allocate thread arguments (used for both parallel and sequential) */
    thread_args = (ctools_var_io_args *)malloc(nvars * sizeof(ctools_var_io_args));
    if (thread_args == NULL) {
        return STATA_ERR_MEMORY;
    }

    /* Initialize thread arguments */
    for (j = 0; j < nvars; j++) {
        thread_args[j].var = &data->vars[j];
        thread_args[j].var_idx = (int)(j + 1);
        thread_args[j].obs1 = obs1;
        thread_args[j].nobs = nobs;
        thread_args[j].is_string = 0;  /* Not used by store */
        thread_args[j].success = 0;
    }

    /* Decide whether to use parallel storing */
    use_parallel = (nvars >= 2);

    if (use_parallel) {
        /* Parallel storing using thread pool */
        ctools_thread_pool pool;
        if (ctools_pool_init(&pool, nvars, thread_args, sizeof(ctools_var_io_args)) != 0) {
            free(thread_args);
            return STATA_ERR_MEMORY;
        }

        ctools_pool_run(&pool, store_variable_thread);
        ctools_pool_free(&pool);
    } else {
        /* Sequential storing - reuse thread function directly */
        for (j = 0; j < nvars; j++) {
            store_variable_thread(&thread_args[j]);
        }
    }

    free(thread_args);
    return STATA_OK;
}

/* ===========================================================================
   Selective Data Loading (Stata -> C) - Load only specified variables
   =========================================================================== */

/*
    Load only specified variables from Stata into C memory.

    This function loads a subset of variables, identified by their 1-based
    Stata indices. This is more memory-efficient when only specific columns
    are needed (e.g., key variables for merge operations).

    @param data        [out] stata_data structure to populate
    @param var_indices [in]  Array of 1-based Stata variable indices
    @param nvars       [in]  Number of variables to load
    @param obs_start   [in]  First observation (1-based), 0 = use SF_in1()
    @param obs_end     [in]  Last observation (1-based), 0 = use SF_in2()

    @return STATA_OK on success, error code otherwise
*/
stata_retcode ctools_data_load_selective(stata_data *data, int *var_indices,
                                          size_t nvars, size_t obs_start, size_t obs_end)
{
    size_t i;
    size_t nobs;
    size_t obs1;
    ctools_var_io_args *thread_args = NULL;
    size_t j;
    int use_parallel;

    if (data == NULL || var_indices == NULL || nvars == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Determine observation range */
    if (obs_start == 0 || obs_end == 0) {
        obs1 = (size_t)SF_in1();
        nobs = (size_t)(SF_in2() - SF_in1() + 1);
    } else {
        obs1 = obs_start;
        nobs = obs_end - obs_start + 1;
    }

    /* Initialize the data structure */
    stata_data_init(data);
    data->nobs = nobs;
    data->nvars = nvars;

    /* Allocate array of variables */
    data->vars = (stata_variable *)calloc(nvars, sizeof(stata_variable));
    if (data->vars == NULL) {
        return STATA_ERR_MEMORY;
    }

    /* Allocate sort order array */
    data->sort_order = (size_t *)malloc(nobs * sizeof(size_t));
    if (data->sort_order == NULL) {
        stata_data_free(data);
        return STATA_ERR_MEMORY;
    }

    /* Initialize sort order to identity permutation - unrolled */
    size_t idx_end = nobs - (nobs % 8);
    for (i = 0; i < idx_end; i += 8) {
        data->sort_order[i] = i;
        data->sort_order[i + 1] = i + 1;
        data->sort_order[i + 2] = i + 2;
        data->sort_order[i + 3] = i + 3;
        data->sort_order[i + 4] = i + 4;
        data->sort_order[i + 5] = i + 5;
        data->sort_order[i + 6] = i + 6;
        data->sort_order[i + 7] = i + 7;
    }
    for (; i < nobs; i++) {
        data->sort_order[i] = i;
    }

    /* Allocate thread arguments */
    thread_args = (ctools_var_io_args *)malloc(nvars * sizeof(ctools_var_io_args));
    if (thread_args == NULL) {
        stata_data_free(data);
        return STATA_ERR_MEMORY;
    }

    /* Initialize thread arguments with SPECIFIED variable indices */
    for (j = 0; j < nvars; j++) {
        thread_args[j].var = &data->vars[j];
        thread_args[j].var_idx = var_indices[j];  /* Use specified index (1-based) */
        thread_args[j].obs1 = obs1;
        thread_args[j].nobs = nobs;
        thread_args[j].is_string = SF_var_is_string((ST_int)var_indices[j]);
        thread_args[j].success = 0;
    }

    /* Decide whether to use parallel loading */
    use_parallel = (nvars >= 2);

    if (use_parallel) {
        /* Parallel loading using thread pool */
        ctools_thread_pool pool;
        if (ctools_pool_init(&pool, nvars, thread_args, sizeof(ctools_var_io_args)) != 0) {
            free(thread_args);
            stata_data_free(data);
            return STATA_ERR_MEMORY;
        }

        int pool_result = ctools_pool_run(&pool, load_variable_thread);
        ctools_pool_free(&pool);
        free(thread_args);

        if (pool_result != 0) {
            stata_data_free(data);
            return STATA_ERR_MEMORY;
        }
    } else {
        /* Sequential loading */
        for (j = 0; j < nvars; j++) {
            if (load_variable_thread(&thread_args[j]) != NULL) {
                free(thread_args);
                stata_data_free(data);
                return STATA_ERR_MEMORY;
            }
        }
        free(thread_args);
    }

    return STATA_OK;
}

/* ===========================================================================
   Streaming Permuted Variable Write (C -> Stata with reordering)
   =========================================================================== */

/*
    Write a single variable to Stata with row permutation.

    For each output row i, reads from Stata row source_rows[i] and writes
    to Stata row i. This enables in-place reordering without loading the
    entire variable into C memory first.

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

    if (is_string) {
        /* String variable: allocate buffer, gather, scatter */
        char **buf = (char **)malloc(output_nobs * sizeof(char *));
        if (!buf) return STATA_ERR_MEMORY;

        char strbuf[2048];

        /* GATHER: Read from source positions */
        for (i = 0; i < output_nobs; i++) {
            if (source_rows[i] >= 0) {
                SF_sdata(stata_var, (ST_int)(source_rows[i] + obs1), strbuf);
                buf[i] = strdup(strbuf);
            } else {
                buf[i] = strdup("");  /* Missing -> empty string */
            }
            if (!buf[i]) {
                /* Cleanup on allocation failure */
                for (size_t k = 0; k < i; k++) free(buf[k]);
                free(buf);
                return STATA_ERR_MEMORY;
            }
        }

        /* SCATTER: Write to sequential output positions */
        for (i = 0; i < output_nobs; i++) {
            SF_sstore(stata_var, (ST_int)(i + obs1), buf[i]);
            free(buf[i]);
        }
        free(buf);

    } else {
        /* Numeric variable: allocate buffer, gather, scatter */
        double *buf = (double *)malloc(output_nobs * sizeof(double));
        if (!buf) return STATA_ERR_MEMORY;

        /* GATHER: Read from source positions */
        for (i = 0; i < output_nobs; i++) {
            if (source_rows[i] >= 0) {
                SF_vdata(stata_var, (ST_int)(source_rows[i] + obs1), &buf[i]);
            } else {
                buf[i] = SV_missval;  /* Missing value */
            }
        }

        /* SCATTER: Write to sequential output positions */
        size_t i_end = output_nobs - (output_nobs % 8);
        for (i = 0; i < i_end; i += 8) {
            SF_vstore(stata_var, (ST_int)(i + obs1), buf[i]);
            SF_vstore(stata_var, (ST_int)(i + 1 + obs1), buf[i + 1]);
            SF_vstore(stata_var, (ST_int)(i + 2 + obs1), buf[i + 2]);
            SF_vstore(stata_var, (ST_int)(i + 3 + obs1), buf[i + 3]);
            SF_vstore(stata_var, (ST_int)(i + 4 + obs1), buf[i + 4]);
            SF_vstore(stata_var, (ST_int)(i + 5 + obs1), buf[i + 5]);
            SF_vstore(stata_var, (ST_int)(i + 6 + obs1), buf[i + 6]);
            SF_vstore(stata_var, (ST_int)(i + 7 + obs1), buf[i + 7]);
        }
        for (; i < output_nobs; i++) {
            SF_vstore(stata_var, (ST_int)(i + obs1), buf[i]);
        }

        free(buf);
    }

    return STATA_OK;
}

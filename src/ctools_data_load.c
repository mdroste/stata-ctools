/*
    ctools_data_load.c
    Stata to C Data Transfer Module

    Reads data from Stata's data space into C memory using the Stata Plugin
    Interface (SPI). This module is designed to be reusable for any operation
    requiring bulk data transfer from Stata, not just sorting.

    Architecture:
    - Column-major storage: Each variable is stored in a contiguous array
    - Parallel loading: One thread per variable when nvars >= 2
    - Uses SF_vdata() for numeric and SF_sdata() for string variables
    - SPI convention: variable index first, then observation index (both 1-based)

    Performance Optimizations:
    - 8x loop unrolling reduces loop overhead and improves instruction pipelining
    - Parallel variable loading overlaps I/O across columns
    - SD_FASTMODE (compile flag) disables SPI bounds checking

    Thread Safety:
    - Each thread reads from a different Stata variable (column)
    - No shared state during reading phase
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_threads.h"

/* Maximum string buffer size for Stata string variables */
#define STATA_STR_MAXLEN 2045

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

/*
    ctools_data_store.c
    C to Stata Data Transfer Module

    Writes data from C memory back to Stata's data space using the Stata
    Plugin Interface (SPI). This module is designed to be reusable for any
    operation requiring bulk data transfer to Stata, not just sorting.

    Precondition:
    - Data in C memory must already be in final order (i.e., permutation
      applied during sort phase). This module performs purely sequential writes.

    Architecture:
    - Parallel writing: One thread per variable when nvars >= 2
    - Uses SF_vstore() for numeric and SF_sstore() for string variables
    - SPI convention: variable index first, then observation index (both 1-based)

    Performance Optimizations:
    - 8x loop unrolling reduces loop overhead and improves instruction pipelining
    - Parallel variable writing overlaps I/O across columns
    - Sequential memory access pattern (no permutation indirection)
    - SD_FASTMODE (compile flag) disables SPI bounds checking

    Thread Safety:
    - Each thread writes to a different Stata variable (column)
    - No shared state during writing phase
*/

#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_threads.h"

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

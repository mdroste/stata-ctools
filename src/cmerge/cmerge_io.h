/*
 * cmerge_io.h
 * Data streaming and I/O functions for cmerge
 *
 * Provides optimized variable-level I/O with prefetching and batching.
 */

#ifndef CMERGE_IO_H
#define CMERGE_IO_H

#include <stdlib.h>
#include <stdint.h>
#include "stplugin.h"
#include "cmerge_join.h"
#include "../ctools_types.h"

/* Forward declaration for using cache */
typedef struct {
    ctools_filtered_data keys;
    ctools_filtered_data keepusing;
    size_t nobs;
    int nkeys;
    int n_keepusing;
    int loaded;
    int merge_by_n;
} cmerge_using_cache_t;

/* External reference to global cache (defined in cmerge_impl.c) */
extern cmerge_using_cache_t g_using_cache;

/* Thread arguments for keepusing variable write */
typedef struct {
    int keepusing_idx;              /* Index in keepusing array */
    ST_int dest_idx;                /* Destination Stata variable index */
    int is_shared;                  /* Is this a shared variable? */
    int update_mode;                /* Update missing master values */
    int replace_mode;               /* Replace all master values */
    cmerge_output_spec_t *specs;    /* Output specifications */
    size_t output_nobs;             /* Number of output observations */
    int success;                    /* Thread result */
} cmerge_keepusing_write_args_t;

/*
 * Thread function to write keepusing variables from cache to Stata.
 * Handles both numeric and string variables with update/replace modes.
 *
 * @param arg  Pointer to cmerge_keepusing_write_args_t
 * @return     NULL on success, non-NULL on error
 */
void *cmerge_write_keepusing_var_thread(void *arg);

#endif /* CMERGE_IO_H */

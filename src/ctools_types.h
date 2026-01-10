/*
    ctools_types.h
    common type definitions for ctools modules

    Defines shared data structures used across Stata plugin modules (load, sort,
    store). These types enable modular, reusable components for Stata plugins.

    Memory Model:
    - Column-major: each variable stored as contiguous array
    - Numeric: double[] (8 bytes per observation)
    - String: char*[] (pointer array, each pointing to heap-allocated string)
*/

#ifndef CTOOLS_TYPES_H
#define CTOOLS_TYPES_H

#include <stdint.h>
#include <stddef.h>

/*
    Return codes for all Stata-C operations.
    STATA_OK (0) indicates success; non-zero indicates error.
*/
typedef enum {
    STATA_OK = 0,               /* Success */
    STATA_ERR_MEMORY = 1,       /* Memory allocation failed */
    STATA_ERR_INVALID_INPUT = 2,/* NULL pointer or invalid parameter */
    STATA_ERR_STATA_READ = 3,   /* Stata plugin read error */
    STATA_ERR_STATA_WRITE = 4,  /* Stata plugin write error */
    STATA_ERR_UNSUPPORTED_TYPE = 5 /* Unsupported variable type */
} stata_retcode;

// Variable type: numeric (double) or string
typedef enum {
    STATA_TYPE_DOUBLE = 0,      /* Numeric variable (8-byte double) */
    STATA_TYPE_STRING = 1       /* String variable (char* array) */
} stata_vartype;

//  Single variable's data in C memory
// Storage is contiguous for cache-friendly access patterns
typedef struct {
    stata_vartype type;         /* STATA_TYPE_DOUBLE or STATA_TYPE_STRING */
    size_t nobs;                /* Number of observations */
    union {
        double *dbl;            /* Numeric: contiguous double[nobs] */
        char **str;             /* String: char*[nobs], each heap-allocated */
    } data;
    size_t str_maxlen;          /* Max string length (string vars only) */
} stata_variable;

/*
    Complete dataset in C memory.

    Layout:
    - vars[0..nvars-1]: array of stata_variable, one per Stata variable
    - sort_order[0..nobs-1]: permutation array (0-based), used internally by sort

    Lifecycle operation example with sorting (see ctools_data_io.c):
    1. ctools_data_load: populates vars[], sort_order = identity
    2. ctools_sort_radix_lsd: sorts data, applies permutation to all vars[]
    3. ctools_data_store: writes vars[] to Stata (sequential)
    4. stata_data_free: releases all allocated memory

*/
typedef struct {
    size_t nobs;                /* Number of observations */
    size_t nvars;               /* Number of variables (columns) */
    stata_variable *vars;       /* Array of all variables [nvars] */
    size_t *sort_order;         /* Permutation array (0-based) [nobs] */
} stata_data;

// Performance timing
typedef struct {
    double load_time;           /* Stata → C transfer time (seconds) */
    double sort_time;           /* Radix sort time (seconds) */
    double store_time;          /* C → Stata transfer time (seconds) */
    double total_time;          /* Wall-clock total (seconds) */
} stata_timer;

// Initialize stata_data to safe empty state
// Sets all pointers to NULL, all counts to 0
void stata_data_init(stata_data *data);

// Free all memory associated with stata_data.
// Safe to call multiple times; resets structure to empty state
void stata_data_free(stata_data *data);

/* ---------------------------------------------------------------------------
   Data I/O Operations
   --------------------------------------------------------------------------- */

// Load all variables from Stata into C memory
stata_retcode ctools_data_load(stata_data *data, size_t nvars);

// Load only specified variables from Stata into C memory
// var_indices: array of 1-based Stata variable indices to load
// nvars: number of variables to load
// obs_start, obs_end: observation range (1-based, inclusive), 0 means use SF_in1/SF_in2
stata_retcode ctools_data_load_selective(stata_data *data, int *var_indices,
                                          size_t nvars, size_t obs_start, size_t obs_end);

// Write all variables from C memory back to Stata
stata_retcode ctools_data_store(stata_data *data, size_t obs1);

// Write specific output variable values to Stata using a source row mapping
// For each output row i, reads from source_rows[i] (0-based, -1 = missing)
// var_idx: 1-based Stata variable index
// Handles both numeric and string variables
stata_retcode ctools_stream_var_permuted(int var_idx, int64_t *source_rows,
                                          size_t output_nobs, size_t obs1);

// Sort data using parallel LSD radix sort
stata_retcode ctools_sort_radix_lsd(stata_data *data, int *sort_vars, size_t nsort);

// Sort data and return permutation mapping (sorted_idx -> original_idx)
// If perm_out is non-NULL, fills it with the permutation before applying
// perm_out must be pre-allocated with at least data->nobs elements
stata_retcode ctools_sort_radix_lsd_with_perm(stata_data *data, int *sort_vars,
                                               size_t nsort, size_t *perm_out);

// Sort data using parallel MSD radix sort
// MSD (Most Significant Digit) radix sort is optimized for variable-length strings
// where it can short-circuit on unique prefixes. Best for:
// - Variable-length strings with common prefixes
// - Data with low entropy in high-order bits/characters
stata_retcode ctools_sort_radix_msd(stata_data *data, int *sort_vars, size_t nsort);

// Sort data using Timsort
// Timsort is an adaptive, stable, hybrid sort that combines merge sort and
// insertion sort. Best for:
// - Partially sorted data (panel data, time series)
// - Data with natural runs
// - When O(N) best case for nearly-sorted data is desired
stata_retcode ctools_sort_timsort(stata_data *data, int *sort_vars, size_t nsort);

/* ---------------------------------------------------------------------------
   Sort Algorithm Selection
   --------------------------------------------------------------------------- */

// Sort algorithm enumeration
typedef enum {
    SORT_ALG_LSD = 0,    // LSD radix sort (default) - best for fixed-width keys
    SORT_ALG_MSD = 1,    // MSD radix sort - best for variable-length strings
    SORT_ALG_TIMSORT = 2 // Timsort - best for partially sorted data
} sort_algorithm_t;

// Sort data using specified algorithm
// Wrapper function that dispatches to the appropriate sort implementation
stata_retcode ctools_sort(stata_data *data, int *sort_vars, size_t nsort,
                          sort_algorithm_t algorithm);

#endif /* CTOOLS_TYPES_H */

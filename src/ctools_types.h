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

    Lifecycle operation example with sorting:
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

// Write all variables from C memory back to Stata
stata_retcode ctools_data_store(stata_data *data, size_t obs1);

// Sort data using parallel LSD radix sort
stata_retcode ctools_sort_radix_lsd(stata_data *data, int *sort_vars, size_t nsort);

#endif /* CTOOLS_TYPES_H */

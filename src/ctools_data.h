#ifndef CTOOLS_DATA_H
#define CTOOLS_DATA_H
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>

/*
    Return codes for all Stata-C operations.
    STATA_OK (0) indicates success; non-zero indicates error.
*/
typedef enum
{
    STATA_OK = 0,                   /* Success */
    STATA_ERR_MEMORY = 1,           /* Memory allocation failed */
    STATA_ERR_INVALID_INPUT = 2,    /* NULL pointer or invalid parameter */
    STATA_ERR_STATA_READ = 3,       /* Stata plugin read error */
    STATA_ERR_STATA_WRITE = 4,      /* Stata plugin write error */
    STATA_ERR_UNSUPPORTED_TYPE = 5, /* Unsupported variable type */
    STATA_ERR_CANCELLED = 6,
    STATA_ERR_SINGULAR = 7,
    STATA_ERR_NONCONVERGENCE = 8
} stata_retcode;

/* Internal codes must never escape as raw Stata r() values. Call only at
 * command boundaries, never on an already-native SPI return code. */
static inline int ctools_stata_rc(stata_retcode status)
{
    switch (status)
    {
    case STATA_OK:
        return 0;
    case STATA_ERR_MEMORY:
        return 920;
    case STATA_ERR_INVALID_INPUT:
        return 198;
    case STATA_ERR_STATA_READ:
    case STATA_ERR_STATA_WRITE:
        return 459;
    case STATA_ERR_UNSUPPORTED_TYPE:
        return 109;
    case STATA_ERR_CANCELLED:
        return 1;
    case STATA_ERR_SINGULAR:
        return 506;
    case STATA_ERR_NONCONVERGENCE:
        return 430;
    default:
        return 459; /* Unknown internal failure, never accidental success. */
    }
}

/* SPI r(1) means user cancellation, not allocation failure. Keep this
 * conversion at the read/write callback, while the original rc is available. */
static inline stata_retcode ctools_spi_status(int rc, stata_retcode failure)
{
    return rc == 0 ? STATA_OK : rc == 1 ? STATA_ERR_CANCELLED : failure;
}

/*
    Permutation index type for sort operations.

    Stata's SPI limits observations to < 2^31 (SF_in2 returns int).
    Using uint32_t instead of size_t saves 50% memory for permutation arrays
    and improves cache utilization during sorting.

    Memory savings for 10M rows:
    - size_t (64-bit): 80MB per permutation array
    - uint32_t (32-bit): 40MB per permutation array

    The PERM_IDX_MAX constant can be used for validation.
*/
typedef uint32_t perm_idx_t;
#define PERM_IDX_MAX UINT32_MAX

// Variable type: numeric (double) or string
typedef enum
{
    STATA_TYPE_DOUBLE = 0, /* Numeric variable (8-byte double) */
    STATA_TYPE_STRING = 1  /* String variable (char* array) */
} stata_vartype;

/* Zero initialization preserves the default owning/aligned storage contract.
 * Borrowed contents must outlive every store that reads the column. */
typedef enum
{
    CTOOLS_BUFFER_ALIGNED = 0,
    CTOOLS_BUFFER_HEAP,
    CTOOLS_BUFFER_BORROWED
} ctools_buffer_storage;
typedef enum
{
    CTOOLS_STRINGS_OWNED = 0,
    CTOOLS_STRINGS_BORROWED
} ctools_string_storage;

//  Single variable's data in C memory
// Storage is contiguous for cache-friendly access patterns
typedef struct
{
    stata_vartype type; /* STATA_TYPE_DOUBLE or STATA_TYPE_STRING */
    size_t nobs;        /* Number of observations */
    union
    {
        double *dbl; /* Numeric: contiguous double[nobs] */
        char **str;  /* String: char*[nobs], each heap-allocated */
    } data;
    size_t str_maxlen; /* Max string length (string vars only) */
    void *_arena;      /* Internal: string arena for fast bulk free (NULL if not used) */
    ctools_buffer_storage buffer_storage;
    ctools_string_storage string_storage;
} stata_variable;

/*
    Complete dataset in C memory.

    Layout:
    - vars[0..nvars-1]: array of stata_variable, one per Stata variable
    - sort_order[0..nobs-1]: permutation array (0-based), used internally by sort

    Lifecycle operation example with sorting (see ctools_data_io.c):
    1. ctools_data_load: populates vars[], sort_order = identity
    2. ctools_sort_dispatch: computes sort_order, leaves vars[] unchanged
    3. ctools_data_store_sorted: writes vars[] to Stata in sort_order
    4. stata_data_free: releases all allocated memory

*/
typedef struct
{
    size_t nobs;          /* Number of observations */
    size_t nvars;         /* Number of variables (columns) */
    stata_variable *vars; /* Array of all variables [nvars] */
    perm_idx_t
        *sort_order; /* Permutation array (0-based) [nobs] - uses uint32_t for 50% memory savings */
} stata_data;

// Performance timing
typedef struct
{
    double load_time;  /* Stata → C transfer time (seconds) */
    double sort_time;  /* Radix sort time (seconds) */
    double store_time; /* C → Stata transfer time (seconds) */
    double total_time; /* Wall-clock total (seconds) */
} stata_timer;

// Initialize stata_data to safe empty state
// Sets all pointers to NULL, all counts to 0
void stata_data_init(stata_data *data);

// Free all memory associated with stata_data.
// Safe to call multiple times; resets structure to empty state
void stata_data_free(stata_data *data);
void stata_variable_free(stata_variable *var);

#endif /* CTOOLS_DATA_H */

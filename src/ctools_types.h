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
#include <stdbool.h>
#include <string.h>

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
    void *_arena;               /* Internal: string arena for fast bulk free (NULL if not used) */
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

// Sort data using Parallel Sample Sort
// Sample sort is the gold standard for parallel sorting. It achieves near-linear
// speedup by ensuring each processor handles ~n/p elements. Best for:
// - Large datasets (millions of observations)
// - High core counts (8+ cores)
// - Random or semi-random data distribution
stata_retcode ctools_sort_sample(stata_data *data, int *sort_vars, size_t nsort);

// Sort data using Parallel Counting Sort
// Counting sort is optimal for integer data with a small range of values. Best for:
// - Integer data with range < 1,000,000
// - Categorical variables (year, state codes, etc.)
// Returns STATA_ERR_UNSUPPORTED_TYPE if data is not suitable
stata_retcode ctools_sort_counting(stata_data *data, int *sort_vars, size_t nsort);

// Check if counting sort is suitable for a variable
// Returns 1 if suitable (integer data with small range), 0 otherwise
int ctools_counting_sort_suitable(stata_data *data, int var_idx);

// Sort data using Parallel Merge Sort
// A stable, predictable O(n log n) parallel algorithm. Best for:
// - When stable sort is required
// - When predictable performance is needed
// - General-purpose parallel sorting
stata_retcode ctools_sort_merge(stata_data *data, int *sort_vars, size_t nsort);

// Sort data using IPS4o (In-place Parallel Super Scalar Samplesort)
// A highly efficient parallel sorting algorithm combining:
// - In-place sorting with minimal auxiliary memory
// - Branchless bucket classification for superscalar performance
// - Block-based data movement for cache efficiency
// Best for:
// - Large datasets where memory efficiency matters
// - High core counts with good cache utilization
// - General-purpose high-performance sorting
stata_retcode ctools_sort_ips4o(stata_data *data, int *sort_vars, size_t nsort);

// Sort data using IPS4o and return permutation mapping (sorted_idx -> original_idx)
// If perm_out is non-NULL, fills it with the permutation before applying
// perm_out must be pre-allocated with at least data->nobs elements
stata_retcode ctools_sort_ips4o_with_perm(stata_data *data, int *sort_vars,
                                           size_t nsort, size_t *perm_out);

/* ---------------------------------------------------------------------------
   Sort Algorithm Selection
   --------------------------------------------------------------------------- */

// Sort algorithm enumeration
typedef enum {
    SORT_ALG_LSD = 0,      // LSD radix sort (default) - best for fixed-width keys
    SORT_ALG_MSD = 1,      // MSD radix sort - best for variable-length strings
    SORT_ALG_TIMSORT = 2,  // Timsort - best for partially sorted data
    SORT_ALG_SAMPLE = 3,   // Sample sort - best for large datasets with many cores
    SORT_ALG_COUNTING = 4, // Counting sort - best for integer data with small range
    SORT_ALG_MERGE = 5,    // Parallel merge sort - stable, predictable O(n log n)
    SORT_ALG_IPS4O = 6     // IPS4o - in-place parallel super scalar samplesort
} sort_algorithm_t;

// Sort data using specified algorithm
// Wrapper function that dispatches to the appropriate sort implementation
stata_retcode ctools_sort(stata_data *data, int *sort_vars, size_t nsort,
                          sort_algorithm_t algorithm);

/* ---------------------------------------------------------------------------
   Type Conversion Utilities

   Shared functions for fast numeric/string conversion used across ctools.
   These avoid sprintf/sscanf overhead for performance-critical code paths.
   --------------------------------------------------------------------------- */

/*
    Convert double to sortable uint64 for radix sort.

    IEEE 754 doubles don't sort correctly as raw bit patterns because:
    - Negative numbers have sign bit set but are "less than" positive
    - Negative numbers sort in reverse order when compared as unsigned

    This function transforms the bit pattern so that:
    - Positive numbers: flip sign bit (0x8000... becomes 0x0000...)
    - Negative numbers: flip all bits (preserves ordering)
    - Missing values: map to UINT64_MAX (sort to end)

    @param d  Input double value
    @param is_missing  If true, treat as Stata missing value
    @return   Sortable uint64 representation
*/
static inline uint64_t ctools_double_to_sortable(double d, bool is_missing)
{
    uint64_t bits;
    memcpy(&bits, &d, sizeof(bits));

    /* Sort missing values to the end (largest possible value) */
    if (is_missing) {
        return UINT64_MAX;
    }

    /* If negative (sign bit set), flip all bits */
    /* If positive, flip only the sign bit */
    if (bits & ((uint64_t)1 << 63)) {
        bits = ~bits;
    } else {
        bits ^= ((uint64_t)1 << 63);
    }

    return bits;
}

/*
    Fast unsigned integer to string conversion.
    Writes digits in reverse order then reverses.

    @param val  Value to convert
    @param buf  Output buffer (must hold at least 21 chars)
    @return     Number of characters written (not including null terminator)
*/
int ctools_uint64_to_str(uint64_t val, char *buf);

/*
    Fast signed integer to string conversion.

    @param val  Value to convert
    @param buf  Output buffer (must hold at least 22 chars for sign + digits)
    @return     Number of characters written (not including null terminator)
*/
int ctools_int64_to_str(int64_t val, char *buf);

/*
    Fast unsigned integer to string using digit pair lookup table.
    Processes two digits at a time for ~2x speedup on large numbers.

    @param val  Value to convert
    @param buf  Output buffer (must hold at least 21 chars)
    @return     Number of characters written (not including null terminator)
*/
int ctools_uint64_to_str_fast(uint64_t val, char *buf);

/*
    Fast string to double parser.
    Handles integers, decimals, scientific notation, and common missing values.
    Much faster than strtod/atof for typical numeric data.

    Recognized missing value representations:
    - Empty string or whitespace only
    - "." (Stata missing)
    - "NA", "na", "NaN", "nan"

    @param str     Input string (not necessarily null-terminated)
    @param len     Length of input string
    @param result  Output: parsed double value (SV_missval for missing)
    @param missval The value to use for missing (typically SV_missval from stplugin.h)
    @return        true if successfully parsed, false if invalid format
*/
bool ctools_parse_double_fast(const char *str, int len, double *result, double missval);

/*
    Power of 10 lookup table for fast float parsing/formatting.
    Covers 10^0 through 10^22 (full double precision range without overflow).
*/
extern const double ctools_pow10_table[23];

#endif /* CTOOLS_TYPES_H */

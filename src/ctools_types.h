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
    perm_idx_t *sort_order;     /* Permutation array (0-based) [nobs] - uses uint32_t for 50% memory savings */
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

// Write variables from C memory to specified Stata variable indices
// var_indices[j] is the 1-based Stata variable index for data->vars[j]
stata_retcode ctools_data_store_selective(stata_data *data, int *var_indices,
                                           size_t nvars, size_t obs1);

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

// Sort data using IPS4o but only compute sort_order (don't apply permutation)
// After this call, data->sort_order contains the permutation but data is unchanged.
// Call ctools_apply_permutation() separately to apply the permutation.
stata_retcode ctools_sort_ips4o_order_only(stata_data *data, int *sort_vars, size_t nsort);

// Sort data using LSD radix sort but only compute sort_order (don't apply permutation)
stata_retcode ctools_sort_radix_lsd_order_only(stata_data *data, int *sort_vars, size_t nsort);

// Sort data using MSD radix sort but only compute sort_order (don't apply permutation)
stata_retcode ctools_sort_radix_msd_order_only(stata_data *data, int *sort_vars, size_t nsort);

// Sort data using Timsort but only compute sort_order (don't apply permutation)
stata_retcode ctools_sort_timsort_order_only(stata_data *data, int *sort_vars, size_t nsort);

// Sort data using Sample sort but only compute sort_order (don't apply permutation)
stata_retcode ctools_sort_sample_order_only(stata_data *data, int *sort_vars, size_t nsort);

// Sort data using Counting sort but only compute sort_order (don't apply permutation)
stata_retcode ctools_sort_counting_order_only(stata_data *data, int *sort_vars, size_t nsort);

// Sort data using Merge sort but only compute sort_order (don't apply permutation)
stata_retcode ctools_sort_merge_order_only(stata_data *data, int *sort_vars, size_t nsort);

/* ---------------------------------------------------------------------------
   Permutation Application
   --------------------------------------------------------------------------- */

// Apply permutation (sort_order) to all variables in the dataset.
// After calling, data is physically reordered and sort_order is reset to identity.
// This is used when sorting was done with _order_only variants.
stata_retcode ctools_apply_permutation(stata_data *data);

/* ---------------------------------------------------------------------------
   Sort Algorithm Selection
   --------------------------------------------------------------------------- */

// Sort algorithm enumeration
typedef enum {
    SORT_ALG_LSD = 0,      // LSD radix sort - best for fixed-width keys
    SORT_ALG_MSD = 1,      // MSD radix sort - best for variable-length strings
    SORT_ALG_TIMSORT = 2,  // Timsort - best for partially sorted data
    SORT_ALG_SAMPLE = 3,   // Sample sort - best for large datasets with many cores
    SORT_ALG_COUNTING = 4, // Counting sort - best for integer data with small range
    SORT_ALG_MERGE = 5,    // Parallel merge sort - stable, predictable O(n log n)
    SORT_ALG_IPS4O = 6     // IPS4o (default) - in-place parallel super scalar samplesort
} sort_algorithm_t;

/*
    Unified sort dispatcher for order-only sort operations.

    Dispatches to the appropriate *_order_only sort function based on the
    algorithm enum. Computes sort_order but does NOT apply the permutation.
    Call ctools_apply_permutation() separately to reorder the data.

    For SORT_ALG_COUNTING, if the data is unsuitable (non-integer or range too
    large), automatically falls back to SORT_ALG_LSD.

    @param data       [in/out] stata_data with allocated sort_order
    @param sort_vars  [in] Array of 1-based variable indices specifying sort keys
    @param nsort      [in] Number of sort key variables
    @param algorithm  [in] Sort algorithm to use
    @return           STATA_OK on success, or error code
*/
stata_retcode ctools_sort_dispatch(stata_data *data, int *sort_vars, size_t nsort,
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
    Fast double parser with custom decimal and group separators.
    Handles European formats like "1.234,56" (decimal=',', group='.').

    @param str       Input string (need not be null-terminated)
    @param len       Length of input
    @param result    [out] Parsed double value
    @param missval   The value to use for missing
    @param dec_sep   Decimal separator character (default '.')
    @param grp_sep   Group/thousands separator character ('\0' = none)
    @return          true if successfully parsed, false if invalid format
*/
bool ctools_parse_double_with_separators(const char *str, int len, double *result, double missval,
                                          char dec_sep, char grp_sep);

/*
    Power of 10 lookup table for fast float parsing/formatting.
    Covers 10^0 through 10^22 (full double precision range without overflow).
*/
extern const double ctools_pow10_table[23];

/*
    Two-digit lookup table for fast digit pair output.
    "00", "01", "02", ... "99" stored as pairs of characters.
    Usage: buf[0] = CTOOLS_DIGIT_PAIRS[val*2]; buf[1] = CTOOLS_DIGIT_PAIRS[val*2+1];
*/
extern const char CTOOLS_DIGIT_PAIRS[200];

/* ---------------------------------------------------------------------------
   Safe String-to-Number Parsing

   These functions replace unsafe atoi()/atof() which return 0 on parse failure,
   making it impossible to distinguish "0" from invalid input like "abc".
   --------------------------------------------------------------------------- */

#include <stdlib.h>
#include <errno.h>
#include <limits.h>

/*
    Safe string to int conversion.
    Uses strtol internally with full error checking.

    @param str     Input string to parse
    @param result  Output: parsed integer value (unchanged on failure)
    @return        true if successfully parsed, false if:
                   - str is NULL or empty
                   - contains non-numeric characters
                   - value overflows int range
*/
static inline bool ctools_safe_atoi(const char *str, int *result)
{
    if (str == NULL || *str == '\0') {
        return false;
    }

    char *endptr;
    errno = 0;
    long val = strtol(str, &endptr, 10);

    /* Check for conversion errors */
    if (errno == ERANGE || val > INT_MAX || val < INT_MIN) {
        return false;  /* Overflow */
    }
    if (endptr == str) {
        return false;  /* No digits found */
    }
    /* Allow trailing whitespace but not other characters */
    while (*endptr == ' ' || *endptr == '\t') {
        endptr++;
    }
    if (*endptr != '\0') {
        return false;  /* Trailing non-whitespace characters */
    }

    *result = (int)val;
    return true;
}

/*
    Safe string to long conversion.

    @param str     Input string to parse
    @param result  Output: parsed long value (unchanged on failure)
    @return        true if successfully parsed, false on error
*/
static inline bool ctools_safe_atol(const char *str, long *result)
{
    if (str == NULL || *str == '\0') {
        return false;
    }

    char *endptr;
    errno = 0;
    long val = strtol(str, &endptr, 10);

    if (errno == ERANGE) {
        return false;
    }
    if (endptr == str) {
        return false;
    }
    while (*endptr == ' ' || *endptr == '\t') {
        endptr++;
    }
    if (*endptr != '\0') {
        return false;
    }

    *result = val;
    return true;
}

/*
    Safe string to size_t conversion.

    @param str     Input string to parse
    @param result  Output: parsed size_t value (unchanged on failure)
    @return        true if successfully parsed, false on error or negative value
*/
static inline bool ctools_safe_atozu(const char *str, size_t *result)
{
    if (str == NULL || *str == '\0') {
        return false;
    }

    /* Reject negative numbers for size_t */
    const char *p = str;
    while (*p == ' ' || *p == '\t') p++;
    if (*p == '-') {
        return false;
    }

    char *endptr;
    errno = 0;
    unsigned long long val = strtoull(str, &endptr, 10);

    if (errno == ERANGE || val > SIZE_MAX) {
        return false;
    }
    if (endptr == str) {
        return false;
    }
    while (*endptr == ' ' || *endptr == '\t') {
        endptr++;
    }
    if (*endptr != '\0') {
        return false;
    }

    *result = (size_t)val;
    return true;
}

/*
    Safe string to double conversion.

    @param str     Input string to parse
    @param result  Output: parsed double value (unchanged on failure)
    @return        true if successfully parsed, false on error
*/
static inline bool ctools_safe_atof(const char *str, double *result)
{
    if (str == NULL || *str == '\0') {
        return false;
    }

    char *endptr;
    errno = 0;
    double val = strtod(str, &endptr);

    if (errno == ERANGE) {
        return false;
    }
    if (endptr == str) {
        return false;
    }
    while (*endptr == ' ' || *endptr == '\t') {
        endptr++;
    }
    if (*endptr != '\0') {
        return false;
    }

    *result = val;
    return true;
}

#endif /* CTOOLS_TYPES_H */

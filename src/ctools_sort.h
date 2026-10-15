#ifndef CTOOLS_SORT_H
#define CTOOLS_SORT_H
#include "ctools_data.h"

/* ---------------------------------------------------------------------------
   Sort Algorithms

   Each *_order_only function computes data->sort_order without moving the
   loaded columns. Call ctools_apply_permutation() or ctools_data_store_sorted()
   afterwards. Most callers should use ctools_sort_dispatch() below.
   --------------------------------------------------------------------------- */

// LSD radix sort: best for fixed-width keys. Numeric keys and string keys of
// up to CTOOLS_SORT_RECORD_STRING_MAX bytes sort together as 64-bit words
// (a string as one or two); longer strings get one pass per byte position.
stata_retcode ctools_sort_radix_lsd_order_only(stata_data *data, int *sort_vars, size_t nsort);
#define CTOOLS_SORT_RECORD_STRING_MAX 16

// MSD radix sort: optimized for variable-length strings, where it can
// short-circuit on unique prefixes
stata_retcode ctools_sort_radix_msd_order_only(stata_data *data, int *sort_vars, size_t nsort);

// Timsort: adaptive, stable merge/insertion hybrid; O(N) on data with long
// natural runs (panel data, time series)
stata_retcode ctools_sort_timsort_order_only(stata_data *data, int *sort_vars, size_t nsort);

// Parallel sample sort: splitter-based partitioning for large datasets and
// high core counts
stata_retcode ctools_sort_sample_order_only(stata_data *data, int *sort_vars, size_t nsort);

// Parallel counting sort: integer keys with a small range (year, state codes).
// Returns STATA_ERR_UNSUPPORTED_TYPE if the data is not suitable.
stata_retcode ctools_sort_counting_order_only(stata_data *data, int *sort_vars, size_t nsort);

// Parallel merge sort: stable, predictable O(n log n)
stata_retcode ctools_sort_merge_order_only(stata_data *data, int *sort_vars, size_t nsort);

// IPS4o (In-place Parallel Super Scalar Samplesort)
stata_retcode ctools_sort_ips4o_order_only(stata_data *data, int *sort_vars, size_t nsort);

// IPS4o that also physically reorders data->vars and returns the permutation
// (perm_out[sorted_idx] = original_idx; pre-allocated with data->nobs elements).
// Used by cmerge.
stata_retcode ctools_sort_ips4o_with_perm(stata_data *data, int *sort_vars, size_t nsort,
                                          size_t *perm_out);

/* ---------------------------------------------------------------------------
   Permutation Application
   --------------------------------------------------------------------------- */

// Apply permutation (sort_order) to all variables in the dataset.
// Requires owned aligned column buffers; views/heap output vectors use
// order-only sorting plus ctools_data_store_sorted instead.
// After calling, data is physically reordered and sort_order is reset to identity.
// This is used when sorting was done with _order_only variants.
stata_retcode ctools_apply_permutation(stata_data *data);

/* ---------------------------------------------------------------------------
   Sort Algorithm Selection
   --------------------------------------------------------------------------- */

// Sort algorithm enumeration
typedef enum
{
    SORT_ALG_LSD = 0,      // LSD radix sort - best for fixed-width keys
    SORT_ALG_MSD = 1,      // MSD radix sort - best for variable-length strings
    SORT_ALG_TIMSORT = 2,  // Timsort - best for partially sorted data
    SORT_ALG_SAMPLE = 3,   // Sample sort - best for large datasets with many cores
    SORT_ALG_COUNTING = 4, // Counting sort - best for integer data with small range
    SORT_ALG_MERGE = 5,    // Parallel merge sort - stable, predictable O(n log n)
    SORT_ALG_IPS4O = 6,    // IPS4o - in-place parallel super scalar samplesort
    SORT_ALG_AUTO = 7      // Auto-select best algorithm based on data characteristics
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

#endif /* CTOOLS_SORT_H */

/*
 * creghdfe_utils.h
 *
 * Utility functions: timing, hash tables, union-find
 * Part of the ctools Stata plugin suite
 */

#ifndef CREGHDFE_UTILS_H
#define CREGHDFE_UTILS_H

#include "creghdfe_types.h"

/* ========================================================================
 * High-resolution timing
 * ======================================================================== */

double get_time_sec(void);

/* ========================================================================
 * Hash table functions
 * ======================================================================== */

IntHashTable *hash_create(ST_int expected_size);
void hash_destroy(IntHashTable *ht);
ST_int hash_insert_or_get(IntHashTable *ht, ST_int key);

/* Count connected components in a bipartite graph */
ST_int count_connected_components(
    const ST_int *fe1_levels,
    const ST_int *fe2_levels,
    ST_int N,
    ST_int num_levels1,
    ST_int num_levels2
);

/* ========================================================================
 * Sort-based value remapping (replaces hash table approach)
 * ======================================================================== */

/*
 * Remap values to contiguous levels using sorting.
 * Much faster than hash tables for large N with many unique values.
 *
 * @param values     Input: array of N double values (e.g., FE variable values)
 * @param N          Number of observations
 * @param levels_out Output: array of N level assignments (1-indexed), must be pre-allocated
 * @param num_levels Output: number of unique levels found
 * @return 0 on success, non-zero on error
 */
int remap_values_sorted(const double *values, ST_int N, ST_int *levels_out, ST_int *num_levels);

/* ========================================================================
 * Optimized FE remap with fused counting
 * ======================================================================== */

/*
 * Remap FE values to contiguous 1-based levels with fused counting.
 * Uses O(N) counting-based approach for integer values (common case),
 * falls back to O(N log N) sort for non-integer or very sparse values.
 *
 * @param values      Input: N FE values (already filtered, no missing)
 * @param N           Number of observations
 * @param levels_out  Output: 1-based level per obs (pre-allocated, N elements)
 * @param num_levels  Output: number of unique levels
 * @param counts_out  Output: per-level counts, allocated here (caller frees)
 * @param weights     Optional: per-obs weights (NULL if no weights)
 * @param wcounts_out Output: weighted counts, allocated if weights != NULL (caller frees)
 * @return 0 on success, -1 on error
 */
int remap_and_count(const double *values, ST_int N, ST_int *levels_out,
                    ST_int *num_levels, ST_int **counts_out,
                    const double *weights, ST_double **wcounts_out);

#endif /* CREGHDFE_UTILS_H */

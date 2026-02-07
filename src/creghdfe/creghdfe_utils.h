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

/*
 * Sort observations by cluster_id and compute cluster boundaries.
 * Used for streaming VCE computation.
 *
 * @param cluster_ids  Input: array of N cluster IDs (0-indexed)
 * @param N            Number of observations
 * @param num_clusters Number of unique clusters
 * @param sort_perm    Output: permutation array (sort_perm[i] = original index), must be pre-allocated
 * @param boundaries   Output: cluster boundaries (boundaries[c] = start of cluster c), size num_clusters+1
 * @return 0 on success, non-zero on error
 */
int sort_by_cluster(const ST_int *cluster_ids, ST_int N, ST_int num_clusters,
                    ST_int *sort_perm, ST_int *boundaries);

#endif /* CREGHDFE_UTILS_H */

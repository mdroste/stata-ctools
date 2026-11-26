/*
 * hdfe_remap.h
 *
 * Utility functions: timing, connected components, FE value remapping
 * Part of the ctools Stata plugin suite
 */

#ifndef CTOOLS_HDFE_REMAP_H
#define CTOOLS_HDFE_REMAP_H

#include "hdfe_types.h"

/* ========================================================================
 * High-resolution timing
 * ======================================================================== */

double get_time_sec(void);

/* ========================================================================
 * Connected components
 * ======================================================================== */

/* Count connected components in a bipartite graph */
ST_int count_connected_components(const ST_int *fe1_levels, const ST_int *fe2_levels, ST_int N,
                                  ST_int num_levels1, ST_int num_levels2);

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
 * Remap FE values to contiguous 1-based levels, numbered in order of first
 * appearance, with optional counts.
 * Uses O(N) counting-based approach for integer values (common case) and a
 * hash for non-integer or very sparse values; a table that cannot be
 * allocated is an error. The direct-table path runs in parallel over row
 * blocks (same levels) outside parallel regions.
 *
 * @param values      Input: N FE values (already filtered, no missing)
 * @param N           Number of observations
 * @param levels_out  Output: 1-based level per obs (pre-allocated, N elements)
 * @param num_levels  Output: number of unique levels
 * @param counts_out  Output: per-level counts, allocated here (caller frees);
 *                    NULL skips the counts
 * @param weights     Optional: per-obs weights (NULL if no weights)
 * @param wcounts_out Output: weighted counts, allocated if weights != NULL (caller frees)
 * @return 0 on success, -1 on error
 */
/* On failure, num_levels is zero and count outputs are NULL; levels_out
 * may be partially written and must not be read until the call succeeds. */
int remap_and_count(const double *values, ST_int N, ST_int *levels_out, ST_int *num_levels,
                    ST_int **counts_out, const double *weights, ST_double **wcounts_out);

/*
 * The levels of remap_and_count() without counts, for callers that
 * renumber them (hdfe_build_factors) or need only the partition they define
 * (singleton removal): sparse and non-integer codes are then also hashed in
 * parallel, and numbered in an order that depends on the values and N but
 * is not first appearance. Returns 0, or -1 with num_levels zero.
 */
int hdfe_remap_levels(const double *values, ST_int N, ST_int *levels_out, ST_int *num_levels);

#endif /* CTOOLS_HDFE_REMAP_H */

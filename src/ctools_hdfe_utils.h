/*
    ctools_hdfe_utils.h
    Shared utilities for HDFE regression commands (creghdfe, civreghdfe)

    Provides:
    - Singleton detection and removal
    - FE level remapping to contiguous indices
    - Cluster ID remapping
*/

#ifndef CTOOLS_HDFE_UTILS_H
#define CTOOLS_HDFE_UTILS_H

#include "stplugin.h"

/*
    Detect singletons in a single FE factor.

    A singleton is an observation where the FE level appears only once.

    Parameters:
    - levels: FE level array (N x 1), 1-based indices
    - counts: Level count array (must be pre-allocated to max_level+1)
    - mask: Keep mask (1=keep, 0=drop), modified in place
    - N: Number of observations
    - max_level: Maximum level value in levels array

    Returns: Number of singletons found in this factor
*/
ST_int ctools_find_singletons_factor(
    const ST_int *levels,
    ST_int *counts,
    ST_int *mask,
    ST_int N,
    ST_int max_level
);

/*
    Iteratively detect and remove singletons across multiple FE factors.

    Singletons can cascade: removing one singleton may create new singletons
    in other factors. This function iterates until no more singletons are found.

    Parameters:
    - fe_levels: Array of G FE level arrays, each of length N
    - G: Number of FE factors
    - N: Number of observations (updated on return)
    - mask: Output keep mask (1=keep, 0=drop), must be pre-allocated to N
    - max_iter: Maximum iterations (safety limit)
    - verbose: Print progress messages

    Returns: Total number of singletons removed
*/
ST_int ctools_remove_singletons(
    ST_int **fe_levels,
    ST_int G,
    ST_int N,
    ST_int *mask,
    ST_int max_iter,
    ST_int verbose
);

/*
    Remap FE levels to contiguous 1-based indices.

    After singleton removal, FE levels may have gaps. This function
    remaps them to 1, 2, 3, ..., num_levels for efficient array indexing.

    Parameters:
    - levels: FE level array (N x 1), modified in place
    - N: Number of observations
    - num_levels: Output - number of unique levels after remapping

    Returns: 0 on success, -1 on memory allocation failure
*/
ST_int ctools_remap_fe_levels(
    ST_int *levels,
    ST_int N,
    ST_int *num_levels
);

/*
    Remap cluster IDs to contiguous 1-based indices.

    Parameters:
    - cluster_ids: Cluster ID array (N x 1), modified in place
    - N: Number of observations
    - num_clusters: Output - number of unique clusters

    Returns: 0 on success, -1 on memory allocation failure
*/
ST_int ctools_remap_cluster_ids(
    ST_int *cluster_ids,
    ST_int N,
    ST_int *num_clusters
);

/*
    Compact arrays by removing flagged observations.

    Copies data from source to dest for observations where mask[i] == 1.

    Parameters:
    - src: Source array (N_src elements)
    - dest: Destination array (must be pre-allocated to N_dest)
    - mask: Keep mask (1=keep, 0=drop)
    - N_src: Source array length
    - N_dest: Expected destination length (for validation)

    Returns: Actual number of elements copied
*/
ST_int ctools_compact_array_double(
    const ST_double *src,
    ST_double *dest,
    const ST_int *mask,
    ST_int N_src,
    ST_int N_dest
);

ST_int ctools_compact_array_int(
    const ST_int *src,
    ST_int *dest,
    const ST_int *mask,
    ST_int N_src,
    ST_int N_dest
);

/*
    Compact a column-major matrix by removing flagged observations.

    Parameters:
    - src: Source matrix (N_src x K), column-major
    - dest: Destination matrix (N_dest x K), column-major
    - mask: Keep mask (1=keep, 0=drop)
    - N_src: Source row count
    - N_dest: Destination row count
    - K: Number of columns

    Returns: Actual number of rows copied
*/
ST_int ctools_compact_matrix_double(
    const ST_double *src,
    ST_double *dest,
    const ST_int *mask,
    ST_int N_src,
    ST_int N_dest,
    ST_int K
);

/*
    Union-Find (Disjoint Set Union) for connected components.
    Used to efficiently count mobility groups in bipartite graphs.
*/
typedef struct {
    ST_int *parent;
    ST_int *rank;
    ST_int size;
} ctools_UnionFind;

ctools_UnionFind *ctools_uf_create(ST_int size);
void ctools_uf_destroy(ctools_UnionFind *uf);
ST_int ctools_uf_find(ctools_UnionFind *uf, ST_int x);
void ctools_uf_union(ctools_UnionFind *uf, ST_int x, ST_int y);

/*
    Count connected components in a bipartite graph.

    Given two FE level vectors of length N, where:
      - fe1_levels[i] is in [1, num_levels1]
      - fe2_levels[i] is in [1, num_levels2]

    Each observation creates an edge between fe1_levels[i] and fe2_levels[i].
    Returns the number of connected components (mobility groups).
    Returns -1 on error.
*/
ST_int ctools_count_connected_components(
    const ST_int *fe1_levels,
    const ST_int *fe2_levels,
    ST_int N,
    ST_int num_levels1,
    ST_int num_levels2
);

#endif /* CTOOLS_HDFE_UTILS_H */

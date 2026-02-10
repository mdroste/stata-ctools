/*
    ctools_hdfe_utils.h
    Shared utilities for HDFE regression commands (creghdfe, civreghdfe)

    Provides:
    - FE_Factor and HDFE_State type definitions
    - Singleton detection and removal
    - FE level remapping to contiguous indices
    - Cluster ID remapping
    - FE nesting detection
    - HDFE state buffer allocation and cleanup
*/

#ifndef CTOOLS_HDFE_UTILS_H
#define CTOOLS_HDFE_UTILS_H

#include "stplugin.h"

/* ========================================================================
 * Fixed Effect Factor structure
 * ======================================================================== */

typedef struct {
    ST_int num_levels;           /* Number of unique levels (for df calculation) */
    ST_int max_level;            /* Maximum level value (for array indexing) */
    ST_int has_intercept;
    ST_int num_slopes;
    ST_int *levels;              /* Level assignment per obs (1-indexed) */
    ST_double *counts;           /* Unweighted counts per level */
    ST_double *inv_counts;       /* Precomputed 1/counts for fast division */
    ST_double *weighted_counts;  /* Sum of weights per level (NULL if no weights) */
    ST_double *inv_weighted_counts; /* Precomputed 1/weighted_counts */
    ST_double *means;
    /* Sorted observation indices for cache-friendly projection */
    ST_int *sorted_indices;      /* indices[N]: observation indices sorted by level */
    ST_int *sorted_levels;       /* levels[N]: level values in sorted order */
    ST_int *level_offsets;       /* offsets[num_levels+1]: start position per level in sorted arrays */
    ST_int sorted_initialized;   /* Flag: 1 if sorted format is built */
} FE_Factor;

/* ========================================================================
 * Global HDFE State structure
 * Holds all state needed for the CG solver and factor operations
 * ======================================================================== */

typedef struct {
    ST_int G;
    ST_int N;
    ST_int K;
    ST_int in1;
    ST_int in2;
    ST_int has_weights;
    ST_int weight_type;          /* 0=none, 1=aweight, 2=fweight, 3=pweight */
    ST_double *weights;
    ST_double sum_weights;       /* Sum of weights (for fweight df calculation) */
    FE_Factor *factors;
    ST_int maxiter;
    ST_double tolerance;
    ST_int verbose;
    /* Per-thread working buffers for parallel column processing */
    ST_double **thread_cg_r;
    ST_double **thread_cg_u;
    ST_double **thread_cg_v;
    ST_double **thread_proj;      /* Not used in creghdfe but needed by cqreg/civreghdfe */
    ST_double **thread_fe_means;  /* Per-thread means buffers for each FE */
    ST_int num_threads;
    /* Cached data from HDFE init for reuse in partial_out */
    ST_int factors_initialized;  /* Flag indicating factors are ready */
    ST_int df_a;                 /* Degrees of freedom absorbed */
    ST_int mobility_groups;      /* Number of mobility groups (connected components) */
} HDFE_State;

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
    Convert string array to integer cluster IDs using sort-based grouping.

    Equivalent of Stata's `egen group()` but in C. Assigns consecutive
    non-negative integer IDs to each unique string. NULL or empty strings
    are assigned -1 (missing sentinel).

    Parameters:
    - strings: Array of N string pointers (from ctools_data_load string var)
    - N: Number of observations
    - cluster_ids: Output array (pre-allocated to N), will contain group IDs
                   Non-missing: 0, 1, 2, ... (consecutive)
                   Missing (NULL or ""): -1
    - num_groups: Output - number of unique non-missing groups

    Returns: 0 on success, -1 on memory allocation failure
*/
ST_int ctools_strings_to_cluster_ids(
    char **strings,
    ST_int N,
    ST_int *cluster_ids,
    ST_int *num_groups
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

/*
    Check if a fixed effect is nested within a cluster variable.

    An FE is nested if every FE level maps to exactly one cluster.
    Uses a first-seen mapping: for each FE level, records the first
    cluster ID observed and checks subsequent observations match.

    Parameters:
    - fe_levels: FE level array (N x 1), 1-based indices
    - num_fe_levels: Number of unique FE levels
    - cluster_ids: Cluster ID array (N x 1), 0-based indices
    - N: Number of observations

    Returns: 1 if nested, 0 if not, -1 on memory allocation error
*/
ST_int ctools_fe_nested_in_cluster(
    const ST_int *fe_levels,
    ST_int num_fe_levels,
    const ST_int *cluster_ids,
    ST_int N
);

/*
    Compute degrees of freedom absorbed by fixed effects and mobility groups.

    Algorithm:
    1. Sum all FE levels → df_a
    2. For G >= 2: count connected components between first two FEs
       and subtract from df_a
    3. For G > 2: subtract (G-2) additional, add to mobility_groups

    Parameters:
    - factors: Array of G FE_Factor structs (levels must be 1-based, contiguous)
    - G: Number of FE factors
    - N: Number of observations
    - df_a: Output - degrees of freedom absorbed
    - mobility_groups: Output - number of mobility groups

    Returns: 0 on success, -1 on error
*/
ST_int ctools_compute_hdfe_dof(
    const FE_Factor *factors,
    ST_int G,
    ST_int N,
    ST_int *df_a,
    ST_int *mobility_groups
);

/*
    Build a sorted permutation for cache-friendly FE projection (scatter-gather).

    Performs a counting sort on f->levels[0..N-1] (1-based contiguous values)
    to produce sorted_indices and level_offsets arrays. This enables sequential
    memory access during the Kaczmarz scatter-gather, replacing random access
    to the means[] array with sequential iteration by level.

    Only builds the permutation when num_levels > 64 (small factors have no
    cache pressure and don't benefit).

    Memory cost: ~8N bytes for sorted_indices + 8*num_levels for level_offsets.

    Parameters:
    - f: FE_Factor with levels already set (1-based, contiguous after remap)
    - N: Number of observations

    Returns: 0 on success, -1 on memory allocation failure
*/
ST_int ctools_build_sorted_permutation(FE_Factor *f, ST_int N);

/*
    Allocate per-thread CG solver buffers and compute inv_counts/inv_weighted_counts.

    Call this after factors[g].levels, counts, weighted_counts, and num_levels
    are already set up. This function:
    1. Computes inv_counts and inv_weighted_counts for each factor
    2. Allocates thread_cg_r/u/v (and optionally thread_proj)
    3. Allocates thread_fe_means

    Only min(num_threads, max_columns) buffer sets are allocated, since the CG
    solver parallelizes over columns and never uses more threads than columns.

    Parameters:
    - state: HDFE_State with factors already initialized (levels, counts, num_levels set)
    - alloc_proj: 1 to allocate thread_proj buffers (needed by civreghdfe/cqreg), 0 to skip
    - max_columns: maximum number of columns that will be partialled out (K)

    Returns: 0 on success, -1 on memory allocation failure
*/
ST_int ctools_hdfe_alloc_buffers(HDFE_State *state, ST_int alloc_proj, ST_int max_columns);

/*
    Free all dynamically allocated memory in an HDFE_State.

    Frees factors (levels, counts, inv_counts, weighted_counts, inv_weighted_counts,
    means, sorted_indices, sorted_levels), weights, thread buffers, and the
    factors array itself.

    Does NOT free the HDFE_State struct itself (caller may have stack-allocated it).
    Does NOT free state->weights (caller manages weight ownership).
*/
void ctools_hdfe_state_cleanup(HDFE_State *state);

#endif /* CTOOLS_HDFE_UTILS_H */

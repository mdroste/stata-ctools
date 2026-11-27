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

#include <stdint.h>
#include <string.h>
#include "stplugin.h"
#ifdef _OPENMP
#include <omp.h>
#endif

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
    ST_double *slope;            /* Optional single heterogeneous slope (owned) */
    ST_double *slope_center;     /* Weighted within-group means when intercept is absorbed */
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
    ST_double **thread_proj;      /* Per-thread N-buffers, only when alloc_proj=1 (no current caller) */
    ST_double **thread_fe_means;  /* Per-thread means buffers for each FE */
    ST_int num_threads;
    /* Cached data from HDFE init for reuse in partial_out */
    ST_int factors_initialized;  /* Flag indicating factors are ready */
    ST_int df_a;                 /* Degrees of freedom absorbed */
    ST_int mobility_groups;      /* Number of mobility groups (connected components) */
    /* 1: the iteration on the rows takes a residual it has brought to the
     * rounding floor as converged whatever the range of the positive
     * weights; 0 (the default): a range wider than that floor can certify
     * makes it a breakdown (hdfe_solver.c, cg_breakdown_result). Set by
     * cppmlhdfe, whose IRLS weights fall without bound for zero outcomes. */
    ST_int floor_any_weight_span;
} HDFE_State;

/*
    Iteratively detect and remove singletons across multiple FE factors.

    A singleton is an observation alone in one of its FE levels. Singletons
    can cascade: removing one may create new singletons in other factors.
    Removal always runs to the fixed point (no iteration cap), which is the
    sample reghdfe keeps.

    Parameters:
    - fe_levels: Array of G FE level arrays, each of length N (1-based levels)
    - G: Number of FE factors
    - N: Number of observations
    - mask: Output keep mask (1=keep, 0=drop), must be pre-allocated to N
    - fweights: Frequency weights (N x 1, positive integers) or NULL. With
      fweights a level is a singleton only when its total fweight is 1, as in
      reghdfe; pass NULL for unweighted, aweight and pweight estimation.
    - verbose: Print progress messages

    Returns: Total number of singletons removed (observations, not weights),
    or -1 on allocation failure or a negative level (mask then unspecified)
*/
ST_int ctools_remove_singletons(
    ST_int **fe_levels,
    ST_int G,
    ST_int N,
    ST_int *mask,
    const ST_double *fweights,
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
/* Dense zero-based IDs from full numeric values; missing values receive -1. */
ST_int ctools_numeric_to_cluster_ids(const ST_double *values, ST_int N,
                                    ST_int *ids, ST_int *num_groups);

ST_int ctools_strings_to_cluster_ids(
    char **strings,
    ST_int N,
    ST_int *cluster_ids,
    ST_int *num_groups
);

/*
    Group IDs of the pairs (a[i], b[i]): the intersection of two clusterings.

    Memory is O(N) whatever the number of levels of a and b: a dense a-by-b
    table is used only while it is small relative to N, and a sort otherwise;
    both give the same IDs. IDs are dense, 0-based and numbered in order of
    first appearance. a and b must be nonnegative.

    Parameters:
    - a, b: ID arrays (N x 1)
    - N: Number of observations
    - ids: Output array (pre-allocated to N); may not overlap a or b
    - num_groups: Output - number of distinct pairs

    Returns: 0 on success, -1 on memory allocation failure or a negative ID
*/
ST_int ctools_intersect_cluster_ids(const ST_int *a, const ST_int *b, ST_int N,
                                    ST_int *ids, ST_int *num_groups);

/* The same groups numbered in ascending (a, b) order instead. */
ST_int ctools_intersect_cluster_ids_sorted(const ST_int *a, const ST_int *b, ST_int N,
                                           ST_int *ids, ST_int *num_groups);

/*
    Dense 0-based group IDs of arbitrary integer codes, numbered in order of
    first appearance. O(N) memory whatever the range or number of codes: a
    table over the code range while that is small relative to N, a sort
    otherwise. ids may not overlap codes.

    Returns: 0 on success, -1 on memory allocation failure
*/
ST_int ctools_dense_cluster_ids(const ST_int *codes, ST_int N, ST_int *ids, ST_int *num_groups);

/* The same groups numbered in ascending code order instead (what egen
   group() gives, minus one). */
ST_int ctools_dense_cluster_ids_sorted(const ST_int *codes, ST_int N, ST_int *ids, ST_int *num_groups);

ST_int ctools_compact_array_int(
    const ST_int *src,
    ST_int *dest,
    const ST_int *mask,
    ST_int N_src,
    ST_int N_dest
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
    Compute degrees of freedom absorbed by fixed effects and mobility groups,
    from the first pair of effects.

    dof_adjust_type selects reghdfe's dofadjustments():
    - 0 (all) and 3 (pairwise): sum of levels, minus the connected components
      of the first two effects, minus one more for each later effect; the
      mobility groups are those components plus one per later effect
    - 2 (firstpair): as above without the later effects
    - 1 (none): one coefficient is redundant when G >= 2; mobility groups = 1
    With a single effect (types other than none) nothing is redundant and
    the mobility groups are 0.

    Parameters:
    - factors: Array of G FE_Factor structs (levels must be 1-based, contiguous)
    - G: Number of FE factors
    - N: Number of observations
    - dof_adjust_type: 0=all, 1=none, 2=firstpair, 3=pairwise
    - df_a: Output - degrees of freedom absorbed
    - mobility_groups: Output - number of mobility groups

    Returns: 0 (an allocation failure in the component count falls back to
    one mobility group)
*/
ST_int ctools_compute_hdfe_dof(
    const FE_Factor *factors,
    ST_int G,
    ST_int N,
    ST_int dof_adjust_type,
    ST_int *df_a,
    ST_int *mobility_groups
);

/* Mobility-group adjustments of ctools_hdfe_redundant_intercepts(). */
#define CTOOLS_DOF_PAIRWISE  1  /* every pair of effects */
#define CTOOLS_DOF_FIRSTPAIR 2  /* the first pair only */

/*
    Redundant intercept coefficients of each absorbed effect, counted as
    reghdfe's estimate_dof() counts them.

    Effects flagged in nested (those nested within a cluster variable) are
    wholly redundant and take no further part. Of the remaining effects, each
    one after the first has at least one redundant coefficient, and so does
    the first when any effect is nested. With CTOOLS_DOF_PAIRWISE an effect's
    count is the largest number of connected components it forms with any
    earlier remaining effect; with CTOOLS_DOF_FIRSTPAIR only the first such
    pair is examined; with neither no pair is.

    The absorbed degrees of freedom are the levels minus these counts.

    Parameters:
    - factors: Array of G FE_Factor structs (levels must be 1-based, contiguous)
    - G: Number of FE factors
    - N: Number of observations
    - nested: G flags, or NULL when no effect is nested
    - flags: CTOOLS_DOF_PAIRWISE, CTOOLS_DOF_FIRSTPAIR or 0
    - redundant: Output - G counts

    Returns: 0 on success, -1 on memory allocation error
*/
/*
    Component of every level of the first effect in the bipartite graph of
    two effects' levels (rows link them): root[l] (num_levels1 entries) is
    the smallest level of the first effect in the component of level l,
    0-based. Levels are 1-based in fe1_levels/fe2_levels, every level of the
    second effect on some row. Returns 0, or -1 on memory allocation error.
*/
ST_int ctools_connected_component_roots(
    const ST_int *fe1_levels,
    const ST_int *fe2_levels,
    ST_int N,
    ST_int num_levels1,
    ST_int num_levels2,
    ST_int *root
);

ST_int ctools_hdfe_redundant_intercepts(
    const FE_Factor *factors,
    ST_int G,
    ST_int N,
    const ST_int *nested,
    ST_int flags,
    ST_int *redundant
);

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
    - alloc_proj: 1 to also allocate thread_proj buffers, 0 to skip (all current callers pass 0)
    - max_columns: maximum number of columns that will be partialled out (K)

    Returns: 0 on success, -1 on memory allocation failure
*/
ST_int ctools_hdfe_alloc_buffers(HDFE_State *state, ST_int alloc_proj, ST_int max_columns);

/*
    Free all dynamically allocated memory in an HDFE_State.

    Frees factors (levels, counts, inv_counts, weighted_counts, inv_weighted_counts,
    means, slope, slope_center), weights, thread buffers, and the
    factors array itself.

    Does NOT free the HDFE_State struct itself (caller may have stack-allocated it).
    Does NOT free state->weights (caller manages weight ownership).
*/
void ctools_hdfe_state_cleanup(HDFE_State *state);

/*
    Row blocks for the parallel sample-construction passes (FE remap, factor
    renumbering, connected components, nesting): one per OpenMP thread, with
    at least CTOOLS_HDFE_BLOCK_ROWS rows each. Returns 1, meaning serial, inside
    any parallel region, active or not (bootstrap replicates run there), or
    for few rows.
    Callers share the blocks out with work-sharing loops, so a runtime that
    supplies fewer threads still covers every block.
*/
#ifndef CTOOLS_HDFE_BLOCK_ROWS /* tests set a small one */
#define CTOOLS_HDFE_BLOCK_ROWS 65536
#endif
static inline int ctools_hdfe_row_blocks(ST_int N)
{
#ifdef _OPENMP
    /* Inside any parallel region, active or not (a bootstrap batch of one
       replicate): regions nested there would each need a new thread team */
    if (omp_get_level() > 0) return 1;
    int threads = omp_get_max_threads();
    long long most = (long long)N / CTOOLS_HDFE_BLOCK_ROWS;
    if (threads < 2 || most < 2) return 1;
    return most < threads ? (int)most : threads;
#else
    (void)N;
    return 1;
#endif
}

/* Rows [*i0, *i1) of block b of B over N rows */
static inline void ctools_hdfe_block_range(ST_int N, int B, int b, ST_int *i0, ST_int *i1)
{
    *i0 = (ST_int)(((long long)N * b) / B);
    *i1 = (ST_int)(((long long)N * (b + 1)) / B);
}

/*
    Per-level work over many levels in rows of no order, in T ranges of levels
    with a thread each: every thread goes down the rows and takes those whose
    level is in its range, so that no two threads write near each other, none
    waits alone for the whole table, and each level meets its rows in row
    order. ctools_hdfe_level_range() gives the first level of range t of T
    over levels 1..L and the number of levels in it;
    CTOOLS_HDFE_ROWS_IN_RANGE runs ROW, with `row` set, for the rows from
    `start` on whose level is one of the `width` levels from `lo`. The rows
    are screened 64 at a time (a byte per row, read as eight words: byte j of
    a word is bits 8j and up, as on every platform the plugin is built for).
*/
static inline void ctools_hdfe_level_range(ST_int L, int T, int t, uint32_t *lo, uint32_t *width)
{
    *lo = 1 + (uint32_t)(((uint64_t)L * (uint64_t)t) / (uint64_t)T);
    *width = 1 + (uint32_t)(((uint64_t)L * ((uint64_t)t + 1)) / (uint64_t)T) - *lo;
}

#define CTOOLS_HDFE_ROWS_IN_RANGE(lev, N, lo, width, start, ROW) \
    do { \
        unsigned char in_[64]; \
        ST_int i_ = (start); \
        for (; i_ <= (N) - 64; i_ += 64) { \
            for (int j_ = 0; j_ < 64; j_++) \
                in_[j_] = (uint32_t)(lev)[i_ + j_] - (lo) < (width); \
            for (int q_ = 0; q_ < 8; q_++) { \
                uint64_t word_; \
                memcpy(&word_, in_ + 8 * q_, 8); \
                while (word_) { \
                    const ST_int row = i_ + 8 * q_ + (ST_int)(__builtin_ctzll(word_) >> 3); \
                    word_ &= word_ - 1; \
                    ROW; \
                } \
            } \
        } \
        for (; i_ < (N); i_++) \
            if ((uint32_t)(lev)[i_] - (lo) < (width)) { \
                const ST_int row = i_; \
                ROW; \
            } \
    } while (0)

/* Relaxed atomic access to the int tables these passes share (GCC and
   Clang builtins; plain loads and stores on x86-64 and ARM64) */
#define CTOOLS_HDFE_LOAD(p) __atomic_load_n((p), __ATOMIC_RELAXED)
#define CTOOLS_HDFE_STORE(p, v) __atomic_store_n((p), (v), __ATOMIC_RELAXED)
#define CTOOLS_HDFE_CAS(p, expected, v) \
    __atomic_compare_exchange_n((p), (expected), (v), 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)

#endif /* CTOOLS_HDFE_UTILS_H */

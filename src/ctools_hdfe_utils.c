/*
    ctools_hdfe_utils.c
    Shared utilities for HDFE regression commands (creghdfe, civreghdfe)

    Implements singleton detection, FE level remapping, and cluster ID remapping.
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "ctools_hdfe_utils.h"

/*
    Detect singletons in a single FE factor.
*/
ST_int ctools_find_singletons_factor(
    const ST_int *levels,
    ST_int *counts,
    ST_int *mask,
    ST_int N,
    ST_int max_level
)
{
    ST_int i, num_singletons = 0;

    /* First pass: count levels (only for non-masked observations) */
    memset(counts, 0, (max_level + 1) * sizeof(ST_int));
    for (i = 0; i < N; i++) {
        if (mask[i]) {
            ST_int level = levels[i];
            /* Bounds check to prevent out-of-bounds access */
            if (level >= 0 && level <= max_level) {
                counts[level]++;
            }
        }
    }

    /* Second pass: mark singletons */
    for (i = 0; i < N; i++) {
        if (mask[i] == 0) continue;  /* Already dropped */
        ST_int level = levels[i];
        /* Bounds check to prevent out-of-bounds access */
        if (level >= 0 && level <= max_level && counts[level] == 1) {
            mask[i] = 0;  /* Mark for dropping */
            num_singletons++;
        }
    }

    return num_singletons;
}

/*
    Iteratively detect and remove singletons across multiple FE factors.
*/
ST_int ctools_remove_singletons(
    ST_int **fe_levels,
    ST_int G,
    ST_int N,
    ST_int *mask,
    ST_int max_iter,
    ST_int verbose
)
{
    ST_int g, i;
    ST_int num_singletons_total = 0;
    ST_int num_singletons_iter;
    ST_int iter = 0;

    if (G == 0) {
        /* No FE factors means no singletons, but still initialize mask */
        for (i = 0; i < N; i++) mask[i] = 1;
        return 0;
    }

    /* Find max level for each factor */
    ST_int *max_levels = (ST_int *)malloc(G * sizeof(ST_int));
    ST_int **counts = (ST_int **)malloc(G * sizeof(ST_int *));

    if (!max_levels || !counts) {
        free(max_levels);
        free(counts);
        return -1;
    }

    for (g = 0; g < G; g++) {
        max_levels[g] = 0;
        for (i = 0; i < N; i++) {
            if (fe_levels[g][i] > max_levels[g]) {
                max_levels[g] = fe_levels[g][i];
            }
        }
        counts[g] = (ST_int *)calloc(max_levels[g] + 1, sizeof(ST_int));
        if (!counts[g]) {
            for (ST_int gg = 0; gg < g; gg++) free(counts[gg]);
            free(counts);
            free(max_levels);
            return -1;
        }
    }

    /* Initialize mask to all 1s (keep all) */
    for (i = 0; i < N; i++) mask[i] = 1;

    /* Iterate until no more singletons */
    do {
        num_singletons_iter = 0;

        for (g = 0; g < G; g++) {
            ST_int found = ctools_find_singletons_factor(
                fe_levels[g], counts[g], mask, N, max_levels[g]
            );
            num_singletons_iter += found;
        }

        if (num_singletons_iter > 0) {
            num_singletons_total += num_singletons_iter;
        }

        iter++;
    } while (num_singletons_iter > 0 && iter < max_iter);

    if (verbose && num_singletons_total > 0) {
        char buf[256];
        snprintf(buf, sizeof(buf), "ctools: Dropped %d singletons in %d iterations\n",
                 (int)num_singletons_total, iter);
        SF_display(buf);
    }

    /* Cleanup */
    for (g = 0; g < G; g++) free(counts[g]);
    free(counts);
    free(max_levels);

    return num_singletons_total;
}

/*
    Remap cluster IDs to contiguous 1-based indices.
*/
ST_int ctools_remap_cluster_ids(
    ST_int *cluster_ids,
    ST_int N,
    ST_int *num_clusters
)
{
    ST_int i;

    /* Find max cluster ID and check for invalid (negative) values */
    ST_int max_id = 0;
    for (i = 0; i < N; i++) {
        /* Defensive check: negative IDs indicate missing values or corruption */
        if (cluster_ids[i] < 0) {
            return -1;  /* Error: invalid cluster ID */
        }
        if (cluster_ids[i] > max_id) max_id = cluster_ids[i];
    }

    /* Create remap array */
    ST_int *remap = (ST_int *)calloc(max_id + 1, sizeof(ST_int));
    if (!remap) return -1;

    /* First pass: assign contiguous IDs starting from 1 */
    ST_int next_id = 1;
    for (i = 0; i < N; i++) {
        ST_int old_id = cluster_ids[i];
        if (remap[old_id] == 0) {
            remap[old_id] = next_id++;
        }
    }

    *num_clusters = next_id - 1;

    /* Second pass: remap in place */
    for (i = 0; i < N; i++) {
        cluster_ids[i] = remap[cluster_ids[i]];
    }

    free(remap);
    return 0;
}

/*
    Compact a double array by removing flagged observations.
*/
ST_int ctools_compact_array_double(
    const ST_double *src,
    ST_double *dest,
    const ST_int *mask,
    ST_int N_src,
    ST_int N_dest
)
{
    ST_int i, idx = 0;
    (void)N_dest;  /* Used for validation in debug builds */

    for (i = 0; i < N_src; i++) {
        if (mask[i]) {
            dest[idx++] = src[i];
        }
    }

    return idx;
}

/*
    Compact an int array by removing flagged observations.
*/
ST_int ctools_compact_array_int(
    const ST_int *src,
    ST_int *dest,
    const ST_int *mask,
    ST_int N_src,
    ST_int N_dest
)
{
    ST_int i, idx = 0;
    (void)N_dest;

    for (i = 0; i < N_src; i++) {
        if (mask[i]) {
            dest[idx++] = src[i];
        }
    }

    return idx;
}

/*
    Compact a column-major matrix by removing flagged observations.
*/
ST_int ctools_compact_matrix_double(
    const ST_double *src,
    ST_double *dest,
    const ST_int *mask,
    ST_int N_src,
    ST_int N_dest,
    ST_int K
)
{
    ST_int i, k, idx = 0;

    for (i = 0; i < N_src; i++) {
        if (mask[i]) {
            for (k = 0; k < K; k++) {
                dest[k * N_dest + idx] = src[k * N_src + i];
            }
            idx++;
        }
    }

    return idx;
}

/* ========================================================================
 * Build sorted permutation for cache-friendly FE projection
 * ======================================================================== */

ST_int ctools_build_sorted_permutation(FE_Factor *f, ST_int N)
{
    /* The sorted path trades sequential ans[] access (N-sized) for sequential
     * means[] access (num_levels-sized). This is only worthwhile when means[]
     * exceeds L2 cache (~256KB), i.e. num_levels > ~32K. For typical FE factors,
     * means[] fits in L1/L2 and the original sequential-scan of ans[] is far
     * faster because the hardware prefetcher handles it perfectly. Disable
     * by using a threshold that effectively never triggers. */
    (void)N;
    f->sorted_initialized = 0;
    f->sorted_indices = NULL;
    f->sorted_levels = NULL;
    f->level_offsets = NULL;
    return 0;

#if 0  /* Sorted permutation disabled — see comment above */
    if (num_levels <= 64) {
        f->sorted_initialized = 0;
        return 0;
    }

    /* Allocate sorted_indices and level_offsets */
    f->sorted_indices = (ST_int *)malloc((size_t)N * sizeof(ST_int));
    f->level_offsets = (ST_int *)malloc(((size_t)num_levels + 1) * sizeof(ST_int));

    if (!f->sorted_indices || !f->level_offsets) {
        if (f->sorted_indices) { free(f->sorted_indices); f->sorted_indices = NULL; }
        if (f->level_offsets) { free(f->level_offsets); f->level_offsets = NULL; }
        f->sorted_initialized = 0;
        return -1;
    }

    /* Counting sort: levels are 1-based contiguous [1..num_levels] */

    /* Step 1: Count occurrences per level (reuse level_offsets as temp) */
    memset(f->level_offsets, 0, ((size_t)num_levels + 1) * sizeof(ST_int));
    for (i = 0; i < N; i++) {
        f->level_offsets[f->levels[i]]++;  /* levels[i] is 1-based */
    }

    /* Step 2: Compute prefix sums (exclusive) to get starting offsets.
     * level_offsets[lev] = start position for level lev (1-based). */
    {
        ST_int cumsum = 0;
        for (lev = 1; lev <= num_levels; lev++) {
            ST_int count = f->level_offsets[lev];
            f->level_offsets[lev] = cumsum;
            cumsum += count;
        }
        f->level_offsets[0] = 0;  /* Sentinel / unused for 1-based levels */
    }

    /* Step 3: Place indices into sorted order.
     * After this, sorted_indices[level_offsets[lev]..level_offsets[lev]+count-1]
     * contains the observation indices for level lev. */
    {
        /* Use a temporary copy of offsets as write cursors */
        ST_int *cursors = (ST_int *)malloc(((size_t)num_levels + 1) * sizeof(ST_int));
        if (!cursors) {
            free(f->sorted_indices); f->sorted_indices = NULL;
            free(f->level_offsets); f->level_offsets = NULL;
            f->sorted_initialized = 0;
            return -1;
        }
        memcpy(cursors, f->level_offsets, ((size_t)num_levels + 1) * sizeof(ST_int));

        for (i = 0; i < N; i++) {
            ST_int l = f->levels[i];  /* 1-based */
            f->sorted_indices[cursors[l]++] = i;
        }
        free(cursors);
    }

    /* Shift level_offsets to be 0-based indexed: offsets[0..num_levels]
     * where offsets[lev] is the start for 0-based level index lev,
     * and offsets[num_levels] = N. */
    for (lev = 0; lev < num_levels; lev++) {
        f->level_offsets[lev] = f->level_offsets[lev + 1];
    }
    f->level_offsets[num_levels] = N;

    /* Don't need sorted_levels — the level for sorted_indices[j] where
     * level_offsets[lev] <= j < level_offsets[lev+1] is simply lev. */
    f->sorted_levels = NULL;

    f->sorted_initialized = 1;
    return 0;
#endif  /* sorted permutation disabled */
}

/* ========================================================================
 * Union-Find implementation for connected components
 * ======================================================================== */

ctools_UnionFind *ctools_uf_create(ST_int size)
{
    ctools_UnionFind *uf = (ctools_UnionFind *)malloc(sizeof(ctools_UnionFind));
    if (!uf) return NULL;

    uf->parent = (ST_int *)malloc(size * sizeof(ST_int));
    uf->rank = (ST_int *)calloc(size, sizeof(ST_int));
    uf->size = size;

    if (!uf->parent || !uf->rank) {
        if (uf->parent) free(uf->parent);
        if (uf->rank) free(uf->rank);
        free(uf);
        return NULL;
    }

    for (ST_int i = 0; i < size; i++) {
        uf->parent[i] = i;
    }

    return uf;
}

void ctools_uf_destroy(ctools_UnionFind *uf)
{
    if (uf) {
        if (uf->parent) free(uf->parent);
        if (uf->rank) free(uf->rank);
        free(uf);
    }
}

ST_int ctools_uf_find(ctools_UnionFind *uf, ST_int x)
{
    if (uf->parent[x] != x) {
        uf->parent[x] = ctools_uf_find(uf, uf->parent[x]);
    }
    return uf->parent[x];
}

void ctools_uf_union(ctools_UnionFind *uf, ST_int x, ST_int y)
{
    ST_int root_x = ctools_uf_find(uf, x);
    ST_int root_y = ctools_uf_find(uf, y);

    if (root_x == root_y) return;

    if (uf->rank[root_x] < uf->rank[root_y]) {
        uf->parent[root_x] = root_y;
    } else if (uf->rank[root_x] > uf->rank[root_y]) {
        uf->parent[root_y] = root_x;
    } else {
        uf->parent[root_y] = root_x;
        uf->rank[root_x]++;
    }
}

ST_int ctools_count_connected_components(
    const ST_int *fe1_levels,
    const ST_int *fe2_levels,
    ST_int N,
    ST_int num_levels1,
    ST_int num_levels2
)
{
    ST_int total_nodes = num_levels1 + num_levels2;
    ctools_UnionFind *uf = ctools_uf_create(total_nodes);
    if (!uf) return -1;

    ST_int i;

    for (i = 0; i < N; i++) {
        ST_int node1 = fe1_levels[i] - 1;
        ST_int node2 = num_levels1 + fe2_levels[i] - 1;
        ctools_uf_union(uf, node1, node2);
    }

    ST_int num_components = 0;
    ST_int *seen_roots = (ST_int *)calloc(total_nodes, sizeof(ST_int));
    if (!seen_roots) {
        ctools_uf_destroy(uf);
        return -1;
    }

    for (i = 0; i < num_levels1; i++) {
        ST_int root = ctools_uf_find(uf, i);
        if (!seen_roots[root]) {
            seen_roots[root] = 1;
            num_components++;
        }
    }

    free(seen_roots);
    ctools_uf_destroy(uf);

    return num_components;
}

/*
    Check if a fixed effect is nested within a cluster variable.
*/
ST_int ctools_fe_nested_in_cluster(
    const ST_int *fe_levels,
    ST_int num_fe_levels,
    const ST_int *cluster_ids,
    ST_int N)
{
    ST_int *fe_to_cluster = (ST_int *)malloc((size_t)num_fe_levels * sizeof(ST_int));
    if (!fe_to_cluster) return -1;

    for (ST_int i = 0; i < num_fe_levels; i++)
        fe_to_cluster[i] = -1;

    ST_int is_nested = 1;
    for (ST_int i = 0; i < N && is_nested; i++) {
        ST_int fe_level = fe_levels[i] - 1;  /* Convert 1-based to 0-based */
        ST_int clust_id = cluster_ids[i];
        if (fe_to_cluster[fe_level] == -1) {
            fe_to_cluster[fe_level] = clust_id;
        } else if (fe_to_cluster[fe_level] != clust_id) {
            is_nested = 0;
        }
    }

    free(fe_to_cluster);
    return is_nested;
}

/*
    Compute degrees of freedom absorbed by fixed effects and mobility groups.
*/
ST_int ctools_compute_hdfe_dof(
    const FE_Factor *factors,
    ST_int G,
    ST_int N,
    ST_int *df_a,
    ST_int *mobility_groups)
{
    ST_int dfa = 0;
    ST_int mg = 0;

    for (ST_int g = 0; g < G; g++)
        dfa += factors[g].num_levels;

    if (G >= 2) {
        mg = ctools_count_connected_components(
            factors[0].levels, factors[1].levels,
            N, factors[0].num_levels, factors[1].num_levels);

        if (mg < 0) mg = 1;  /* Fallback on error */
        dfa -= mg;

        if (G > 2) {
            ST_int extra = G - 2;
            dfa -= extra;
            mg += extra;
        }
    }

    *df_a = dfa;
    *mobility_groups = mg;
    return 0;
}

/*
    Allocate per-thread CG solver buffers and compute inv_counts/inv_weighted_counts.
    Caller must set state->num_threads, state->N, state->G, state->has_weights,
    and state->factors[g].{num_levels, counts, weighted_counts} before calling.

    Only min(num_threads, max_columns) N-sized buffer sets are allocated, since the
    CG solver parallelizes over columns and never uses more threads than columns.
    This also caps state->num_threads so the OMP pragmas in partial_out_columns
    don't launch more threads than we have buffers for.
*/
ST_int ctools_hdfe_alloc_buffers(HDFE_State *state, ST_int alloc_proj, ST_int max_columns)
{
    ST_int G = state->G;
    ST_int N = state->N;
    ST_int num_threads = state->num_threads;

    /* Cap threads to number of columns — the CG solver parallelizes over columns,
       so we never need more buffer sets than columns to partial out */
    if (max_columns > 0 && num_threads > max_columns) {
        num_threads = max_columns;
        state->num_threads = num_threads;
    }

    /* Compute inv_counts and inv_weighted_counts for each factor */
    for (ST_int g = 0; g < G; g++) {
        ST_int num_lev = state->factors[g].num_levels;

        state->factors[g].inv_counts = (ST_double *)malloc((size_t)num_lev * sizeof(ST_double));
        if (!state->factors[g].inv_counts) {
            return -1;  /* Allocation failed — caller must clean up */
        }
        for (ST_int lev = 0; lev < num_lev; lev++) {
            state->factors[g].inv_counts[lev] =
                (state->factors[g].counts[lev] > 0) ? 1.0 / state->factors[g].counts[lev] : 0.0;
        }

        if (state->has_weights && state->factors[g].weighted_counts) {
            state->factors[g].inv_weighted_counts = (ST_double *)malloc((size_t)num_lev * sizeof(ST_double));
            if (!state->factors[g].inv_weighted_counts) {
                return -1;  /* Allocation failed — caller must clean up */
            }
            for (ST_int lev = 0; lev < num_lev; lev++) {
                state->factors[g].inv_weighted_counts[lev] =
                    (state->factors[g].weighted_counts[lev] > 0) ? 1.0 / state->factors[g].weighted_counts[lev] : 0.0;
            }
        } else {
            state->factors[g].inv_weighted_counts = NULL;
        }
    }

    /* Allocate thread buffer arrays */
    state->thread_cg_r = (ST_double **)calloc((size_t)num_threads, sizeof(ST_double *));
    state->thread_cg_u = (ST_double **)calloc((size_t)num_threads, sizeof(ST_double *));
    state->thread_cg_v = (ST_double **)calloc((size_t)num_threads, sizeof(ST_double *));
    state->thread_proj = alloc_proj ? (ST_double **)calloc((size_t)num_threads, sizeof(ST_double *)) : NULL;
    state->thread_fe_means = (ST_double **)calloc((size_t)num_threads * G, sizeof(ST_double *));

    if (!state->thread_cg_r || !state->thread_cg_u || !state->thread_cg_v ||
        !state->thread_fe_means || (alloc_proj && !state->thread_proj)) {
        return -1;
    }

    /* Allocate per-thread buffers */
    for (ST_int t = 0; t < num_threads; t++) {
        state->thread_cg_r[t] = (ST_double *)malloc((size_t)N * sizeof(ST_double));
        state->thread_cg_u[t] = (ST_double *)malloc((size_t)N * sizeof(ST_double));
        state->thread_cg_v[t] = (ST_double *)malloc((size_t)N * sizeof(ST_double));
        if (!state->thread_cg_r[t] || !state->thread_cg_u[t] || !state->thread_cg_v[t])
            return -1;

        if (alloc_proj) {
            state->thread_proj[t] = (ST_double *)malloc((size_t)N * sizeof(ST_double));
            if (!state->thread_proj[t]) return -1;
        }

        for (ST_int g = 0; g < G; g++) {
            state->thread_fe_means[t * G + g] = (ST_double *)malloc(
                (size_t)state->factors[g].num_levels * sizeof(ST_double));
            if (!state->thread_fe_means[t * G + g]) return -1;
        }
    }

    return 0;
}

/*
    Free all dynamically allocated memory inside an HDFE_State.
    Does NOT free the HDFE_State struct itself.
    Does NOT free state->weights (caller manages weight ownership).
*/
void ctools_hdfe_state_cleanup(HDFE_State *state)
{
    if (!state) return;

    if (state->factors) {
        for (ST_int g = 0; g < state->G; g++) {
            if (state->factors[g].levels) free(state->factors[g].levels);
            if (state->factors[g].counts) free(state->factors[g].counts);
            if (state->factors[g].inv_counts) free(state->factors[g].inv_counts);
            if (state->factors[g].weighted_counts) free(state->factors[g].weighted_counts);
            if (state->factors[g].inv_weighted_counts) free(state->factors[g].inv_weighted_counts);
            if (state->factors[g].means) free(state->factors[g].means);
            if (state->factors[g].sorted_indices) free(state->factors[g].sorted_indices);
            if (state->factors[g].sorted_levels) free(state->factors[g].sorted_levels);
            if (state->factors[g].level_offsets) free(state->factors[g].level_offsets);
        }
        free(state->factors);
        state->factors = NULL;
    }

    /* Free per-thread CG buffers */
    if (state->thread_cg_r) {
        for (ST_int t = 0; t < state->num_threads; t++)
            if (state->thread_cg_r[t]) free(state->thread_cg_r[t]);
        free(state->thread_cg_r);
        state->thread_cg_r = NULL;
    }
    if (state->thread_cg_u) {
        for (ST_int t = 0; t < state->num_threads; t++)
            if (state->thread_cg_u[t]) free(state->thread_cg_u[t]);
        free(state->thread_cg_u);
        state->thread_cg_u = NULL;
    }
    if (state->thread_cg_v) {
        for (ST_int t = 0; t < state->num_threads; t++)
            if (state->thread_cg_v[t]) free(state->thread_cg_v[t]);
        free(state->thread_cg_v);
        state->thread_cg_v = NULL;
    }
    if (state->thread_proj) {
        for (ST_int t = 0; t < state->num_threads; t++)
            if (state->thread_proj[t]) free(state->thread_proj[t]);
        free(state->thread_proj);
        state->thread_proj = NULL;
    }
    if (state->thread_fe_means) {
        for (ST_int t = 0; t < state->num_threads * state->G; t++)
            if (state->thread_fe_means[t]) free(state->thread_fe_means[t]);
        free(state->thread_fe_means);
        state->thread_fe_means = NULL;
    }
}

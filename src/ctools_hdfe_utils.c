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
            counts[levels[i]]++;
        }
    }

    /* Second pass: mark singletons */
    for (i = 0; i < N; i++) {
        if (mask[i] == 0) continue;  /* Already dropped */
        if (counts[levels[i]] == 1) {
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

    if (G == 0) return 0;

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
    Remap FE levels to contiguous 1-based indices.
*/
ST_int ctools_remap_fe_levels(
    ST_int *levels,
    ST_int N,
    ST_int *num_levels
)
{
    ST_int i;

    /* Find max level */
    ST_int max_level = 0;
    for (i = 0; i < N; i++) {
        if (levels[i] > max_level) max_level = levels[i];
    }

    /* Create remap array */
    ST_int *remap = (ST_int *)calloc(max_level + 1, sizeof(ST_int));
    if (!remap) return -1;

    /* First pass: assign contiguous IDs starting from 1 */
    ST_int next_level = 1;
    for (i = 0; i < N; i++) {
        ST_int old_level = levels[i];
        if (remap[old_level] == 0) {
            remap[old_level] = next_level++;
        }
    }

    *num_levels = next_level - 1;

    /* Second pass: remap in place */
    for (i = 0; i < N; i++) {
        levels[i] = remap[levels[i]];
    }

    free(remap);
    return 0;
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

    /* Find max cluster ID */
    ST_int max_id = 0;
    for (i = 0; i < N; i++) {
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

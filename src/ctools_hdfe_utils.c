/*
    ctools_hdfe_utils.c
    Shared utilities for HDFE regression commands (creghdfe, civreghdfe)

    Implements singleton detection, FE level remapping, and cluster ID remapping.
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "ctools_hdfe_utils.h"

/*
    `row` is the only row still counted in one of its levels. Drop it unless
    it is already dropped or, with fweights, its weight is not 1: reghdfe
    (FE_Factor::drop_singletons) drops a level only when its total fweight is
    1, i.e. one row of weight 1 since fweights are positive integers. Dropped
    rows wait on a stack until their level counts are decremented; the stack
    grows by copying, and each row is pushed at most once.
    Returns -1 on allocation failure.
*/
static int singleton_drop(ST_int row, ST_int *mask, const ST_double *fweights,
                          ST_int **stack, size_t *len, size_t *cap)
{
    if (!mask[row] || (fweights && fweights[row] != 1.0)) return 0;
    if (*len == *cap) {
        size_t grown_cap = *cap ? 2 * *cap : 1024;
        ST_int *grown = (ST_int *)malloc(grown_cap * sizeof(ST_int));
        if (!grown) return -1;
        if (*len) memcpy(grown, *stack, *len * sizeof(ST_int));
        free(*stack);
        *stack = grown;
        *cap = grown_cap;
    }
    mask[row] = 0;
    (*stack)[(*len)++] = row;
    return 0;
}

/*
    Iteratively detect and remove singletons across multiple FE factors.

    Peels to the fixed point, which is the sample reghdfe's cyclic loop
    reaches whatever the order of removals. Each level keeps the count of its
    undropped rows and the XOR of their row indices, so a level left with one
    row names that row directly, and dropping a row touches only its G
    levels: total work is O(N*G) however long the singleton chains are.
*/
ST_int ctools_remove_singletons(
    ST_int **fe_levels,
    ST_int G,
    ST_int N,
    ST_int *mask,
    const ST_double *fweights,
    ST_int verbose
)
{
    ST_int g, i;
    ST_int num_singletons = 0;
    ST_int result = -1;
    ST_int **counts = NULL, **row_xor = NULL;
    ST_int *max_levels = NULL, *stack = NULL;
    size_t stack_len = 0, stack_cap = 0;

    /* Initialize mask to all 1s (keep all) */
    for (i = 0; i < N; i++) mask[i] = 1;
    if (G <= 0 || N <= 0) return 0;

    max_levels = (ST_int *)calloc((size_t)G, sizeof(ST_int));
    counts = (ST_int **)calloc((size_t)G, sizeof(ST_int *));
    row_xor = (ST_int **)calloc((size_t)G, sizeof(ST_int *));
    if (!max_levels || !counts || !row_xor) goto cleanup;

    for (g = 0; g < G; g++) {
        for (i = 0; i < N; i++) {
            if (fe_levels[g][i] < 0) goto cleanup;  /* levels are 1-based */
            if (fe_levels[g][i] > max_levels[g]) max_levels[g] = fe_levels[g][i];
        }
        counts[g] = (ST_int *)calloc((size_t)max_levels[g] + 1, sizeof(ST_int));
        row_xor[g] = (ST_int *)calloc((size_t)max_levels[g] + 1, sizeof(ST_int));
        if (!counts[g] || !row_xor[g]) goto cleanup;
        for (i = 0; i < N; i++) {
            counts[g][fe_levels[g][i]]++;
            row_xor[g][fe_levels[g][i]] ^= i;
        }
    }

    /* Seed with the rows that are alone in a level */
    for (g = 0; g < G; g++) {
        for (ST_int lev = 0; lev <= max_levels[g]; lev++) {
            if (counts[g][lev] == 1 &&
                singleton_drop(row_xor[g][lev], mask, fweights, &stack, &stack_len, &stack_cap))
                goto cleanup;
        }
    }

    /* Removing a row can leave a single row in its other levels. When that
     * row is itself dropped but not yet popped, the level empties later. */
    while (stack_len > 0) {
        ST_int row = stack[--stack_len];
        num_singletons++;
        for (g = 0; g < G; g++) {
            ST_int lev = fe_levels[g][row];
            counts[g][lev]--;
            row_xor[g][lev] ^= row;
            if (counts[g][lev] == 1 &&
                singleton_drop(row_xor[g][lev], mask, fweights, &stack, &stack_len, &stack_cap))
                goto cleanup;
        }
    }
    result = num_singletons;

    if (verbose && num_singletons > 0) {
        char buf[256];
        snprintf(buf, sizeof(buf), "ctools: Dropped %d singletons\n", (int)num_singletons);
        SF_display(buf);
    }

cleanup:
    for (g = 0; g < G; g++) {
        if (counts) free(counts[g]);
        if (row_xor) free(row_xor[g]);
    }
    free(counts);
    free(row_xor);
    free(max_levels);
    free(stack);
    return result;
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
    Compare function for string-index pairs (used by qsort in ctools_strings_to_cluster_ids).
*/
typedef struct {
    const char *str;
    ST_int index;
} StringIndexPair;

static int compare_string_index(const void *a, const void *b)
{
    const StringIndexPair *pa = (const StringIndexPair *)a;
    const StringIndexPair *pb = (const StringIndexPair *)b;
    return strcmp(pa->str, pb->str);
}

/*
    Convert string array to integer cluster IDs using sort-based grouping.
*/
ST_int ctools_strings_to_cluster_ids(
    char **strings,
    ST_int N,
    ST_int *cluster_ids,
    ST_int *num_groups
)
{
    if (N <= 0 || !strings || !cluster_ids || !num_groups) return -1;

    /* First pass: count non-missing strings and mark missing */
    ST_int n_valid = 0;
    for (ST_int i = 0; i < N; i++) {
        if (strings[i] == NULL || strings[i][0] == '\0') {
            cluster_ids[i] = -1;  /* Missing sentinel */
        } else {
            cluster_ids[i] = 0;  /* Placeholder, will be filled */
            n_valid++;
        }
    }

    if (n_valid == 0) {
        *num_groups = 0;
        return 0;
    }

    /* Build (string, index) pairs for non-missing strings */
    StringIndexPair *pairs = (StringIndexPair *)malloc(n_valid * sizeof(StringIndexPair));
    if (!pairs) return -1;

    ST_int j = 0;
    for (ST_int i = 0; i < N; i++) {
        if (cluster_ids[i] != -1) {
            pairs[j].str = strings[i];
            pairs[j].index = i;
            j++;
        }
    }

    /* Sort by string value */
    qsort(pairs, n_valid, sizeof(StringIndexPair), compare_string_index);

    /* Assign consecutive group IDs (0-based) */
    ST_int current_group = 0;
    cluster_ids[pairs[0].index] = 0;

    for (ST_int i = 1; i < n_valid; i++) {
        if (strcmp(pairs[i].str, pairs[i - 1].str) != 0) {
            current_group++;
        }
        cluster_ids[pairs[i].index] = current_group;
    }

    *num_groups = current_group + 1;
    free(pairs);
    return 0;
}

/* (a, b) packed into one key, with the row it came from */
typedef struct { unsigned long long key; ST_int index; } ctools_pair_key;

static int ctools_compare_pair_key(const void *x, const void *y)
{
    const ctools_pair_key *p = (const ctools_pair_key *)x;
    const ctools_pair_key *q = (const ctools_pair_key *)y;
    if (p->key != q->key) return (p->key > q->key) - (p->key < q->key);
    return (p->index > q->index) - (p->index < q->index);
}

/* Dense 0-based IDs, numbered in order of first appearance, for N rows whose
   keys are in pairs (index = row). Sorts pairs. Returns 0, or -1 on
   allocation failure. */
static ST_int ctools_ids_from_pair_keys(ctools_pair_key *pairs, ST_int N,
                                        ST_int *ids, ST_int *num_groups)
{
    qsort(pairs, (size_t)N, sizeof(*pairs), ctools_compare_pair_key);

    /* Sorted-order group of every row */
    ST_int groups = 0;
    for (ST_int i = 0; i < N; i++) {
        if (i == 0 || pairs[i].key != pairs[i - 1].key) groups++;
        ids[pairs[i].index] = groups - 1;
    }

    /* Renumber in order of first appearance */
    ST_int *first_seen = (ST_int *)malloc((size_t)groups * sizeof(ST_int));
    if (!first_seen) return -1;
    for (ST_int grp = 0; grp < groups; grp++) first_seen[grp] = -1;
    ST_int next = 0;
    for (ST_int i = 0; i < N; i++) {
        if (first_seen[ids[i]] < 0) first_seen[ids[i]] = next++;
        ids[i] = first_seen[ids[i]];
    }
    free(first_seen);
    *num_groups = groups;
    return 0;
}

/* Pair IDs through a dense (a, b) table, when that table is small relative
   to N: sorted != 0 numbers the pairs in ascending (a, b) order, otherwise in
   order of first appearance, exactly as the sort-based route does. Returns 1
   if done, 0 if the table would be too large or cannot be allocated. */
static int ctools_pair_ids_by_table(const ST_int *a, const ST_int *b, ST_int N,
                                    ST_int max_a, ST_int max_b, int sorted,
                                    ST_int *ids, ST_int *num_groups)
{
    const size_t width = (size_t)max_b + 1;
    const size_t cells = ((size_t)max_a + 1) * width;
    if (cells / width != (size_t)max_a + 1 || cells > 8 * (size_t)N + 65536) return 0;
    ST_int *table = (ST_int *)calloc(cells, sizeof(ST_int));
    if (!table) return 0;

    ST_int groups = 0;
    if (sorted) {
        for (ST_int i = 0; i < N; i++) table[(size_t)a[i] * width + (size_t)b[i]] = 1;
        for (size_t c = 0; c < cells; c++) {
            if (table[c]) table[c] = ++groups;
        }
        for (ST_int i = 0; i < N; i++) ids[i] = table[(size_t)a[i] * width + (size_t)b[i]] - 1;
    } else {
        for (ST_int i = 0; i < N; i++) {
            ST_int *cell = &table[(size_t)a[i] * width + (size_t)b[i]];
            if (*cell == 0) *cell = ++groups;
            ids[i] = *cell - 1;
        }
    }
    free(table);
    *num_groups = groups;
    return 1;
}

static ST_int ctools_pair_ids(const ST_int *a, const ST_int *b, ST_int N, int sorted,
                              ST_int *ids, ST_int *num_groups)
{
    if (N < 0 || !a || !b || !ids || !num_groups) return -1;
    *num_groups = 0;
    if (N == 0) return 0;

    ST_int max_a = 0, max_b = 0;
    for (ST_int i = 0; i < N; i++) {
        if (a[i] < 0 || b[i] < 0) return -1;
        if (a[i] > max_a) max_a = a[i];
        if (b[i] > max_b) max_b = b[i];
    }
    if (ctools_pair_ids_by_table(a, b, N, max_a, max_b, sorted, ids, num_groups)) return 0;

    ctools_pair_key *pairs = (ctools_pair_key *)malloc((size_t)N * sizeof(*pairs));
    if (!pairs) return -1;
    for (ST_int i = 0; i < N; i++) {
        pairs[i].key = ((unsigned long long)(unsigned int)a[i] << 32) | (unsigned int)b[i];
        pairs[i].index = i;
    }
    ST_int result = 0;
    if (sorted) {
        qsort(pairs, (size_t)N, sizeof(*pairs), ctools_compare_pair_key);
        ST_int groups = 0;
        for (ST_int i = 0; i < N; i++) {
            if (i == 0 || pairs[i].key != pairs[i - 1].key) groups++;
            ids[pairs[i].index] = groups - 1;
        }
        *num_groups = groups;
    } else {
        result = ctools_ids_from_pair_keys(pairs, N, ids, num_groups);
    }
    free(pairs);
    return result;
}

/*
    Group IDs of the pairs (a[i], b[i]), numbered in order of first appearance.
*/
ST_int ctools_intersect_cluster_ids(const ST_int *a, const ST_int *b, ST_int N,
                                    ST_int *ids, ST_int *num_groups)
{
    return ctools_pair_ids(a, b, N, 0, ids, num_groups);
}

/*
    Group IDs of the pairs (a[i], b[i]), numbered in ascending (a, b) order.
*/
ST_int ctools_intersect_cluster_ids_sorted(const ST_int *a, const ST_int *b, ST_int N,
                                           ST_int *ids, ST_int *num_groups)
{
    return ctools_pair_ids(a, b, N, 1, ids, num_groups);
}

/* Dense IDs of integer codes; sorted != 0 numbers them in ascending code
   order, otherwise in order of first appearance. A presence table over
   [min, max] is used while it is small relative to N, a sort otherwise;
   both give the same IDs. */
static ST_int ctools_code_ids(const ST_int *codes, ST_int N, int sorted,
                              ST_int *ids, ST_int *num_groups)
{
    if (N < 0 || !codes || !ids || !num_groups) return -1;
    *num_groups = 0;
    if (N == 0) return 0;

    ST_int lo = codes[0], hi = codes[0];
    for (ST_int i = 1; i < N; i++) {
        if (codes[i] < lo) lo = codes[i];
        if (codes[i] > hi) hi = codes[i];
    }
    const unsigned long long span = (unsigned long long)((long long)hi - (long long)lo) + 1ULL;
    if (span <= 8ULL * (unsigned long long)N + 65536ULL) {
        ST_int *table = (ST_int *)calloc((size_t)span, sizeof(ST_int));
        if (table) {
            ST_int groups = 0;
            if (sorted) {
                for (ST_int i = 0; i < N; i++) table[(size_t)((long long)codes[i] - lo)] = 1;
                for (size_t c = 0; c < (size_t)span; c++) {
                    if (table[c]) table[c] = ++groups;
                }
                for (ST_int i = 0; i < N; i++) ids[i] = table[(size_t)((long long)codes[i] - lo)] - 1;
            } else {
                for (ST_int i = 0; i < N; i++) {
                    ST_int *cell = &table[(size_t)((long long)codes[i] - lo)];
                    if (*cell == 0) *cell = ++groups;
                    ids[i] = *cell - 1;
                }
            }
            free(table);
            *num_groups = groups;
            return 0;
        }
    }

    ctools_pair_key *pairs = (ctools_pair_key *)malloc((size_t)N * sizeof(*pairs));
    if (!pairs) return -1;
    for (ST_int i = 0; i < N; i++) {
        /* offset by the minimum so that keys order like the signed codes */
        pairs[i].key = (unsigned long long)((long long)codes[i] - (long long)lo);
        pairs[i].index = i;
    }
    ST_int result = 0;
    if (sorted) {
        qsort(pairs, (size_t)N, sizeof(*pairs), ctools_compare_pair_key);
        ST_int groups = 0;
        for (ST_int i = 0; i < N; i++) {
            if (i == 0 || pairs[i].key != pairs[i - 1].key) groups++;
            ids[pairs[i].index] = groups - 1;
        }
        *num_groups = groups;
    } else {
        result = ctools_ids_from_pair_keys(pairs, N, ids, num_groups);
    }
    free(pairs);
    return result;
}

/*
    Group IDs of integer codes, numbered in order of first appearance.
*/
ST_int ctools_dense_cluster_ids(const ST_int *codes, ST_int N, ST_int *ids, ST_int *num_groups)
{
    return ctools_code_ids(codes, N, 0, ids, num_groups);
}

/*
    Group IDs of integer codes, numbered in ascending code order.
*/
ST_int ctools_dense_cluster_ids_sorted(const ST_int *codes, ST_int N, ST_int *ids, ST_int *num_groups)
{
    return ctools_code_ids(codes, N, 1, ids, num_groups);
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

/* ========================================================================
 * Build sorted permutation for cache-friendly FE projection
 * ======================================================================== */

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
    ST_int dof_adjust_type,
    ST_int *df_a,
    ST_int *mobility_groups)
{
    ST_int dfa = 0;
    ST_int mg = 0;

    for (ST_int g = 0; g < G; g++)
        dfa += factors[g].num_levels;

    if (dof_adjust_type == 1) {
        /* none: no mobility calculation; one coefficient is still redundant
           once two effects are absorbed */
        mg = 1;
        if (G >= 2) dfa -= 1;
    } else if (G >= 2) {
        mg = ctools_count_connected_components(
            factors[0].levels, factors[1].levels,
            N, factors[0].num_levels, factors[1].num_levels);

        if (mg < 0) mg = 1;  /* Fallback on error */
        dfa -= mg;

        if (G > 2 && (dof_adjust_type == 0 || dof_adjust_type == 3)) {
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
    Redundant intercept coefficients of each absorbed effect (reghdfe's rule).
*/
ST_int ctools_hdfe_redundant_intercepts(
    const FE_Factor *factors,
    ST_int G,
    ST_int N,
    const ST_int *nested,
    ST_int flags,
    ST_int *redundant)
{
    ST_int any_nested = 0, first = -1, pairs = 0;

    for (ST_int g = 0; g < G; g++)
        if (nested && nested[g]) any_nested = 1;

    for (ST_int g = 0; g < G; g++) {
        if (nested && nested[g]) {
            redundant[g] = factors[g].num_levels;
            continue;
        }
        redundant[g] = (first >= 0 || any_nested) ? 1 : 0;
        if (first < 0) first = g;
        if (!(flags & (CTOOLS_DOF_PAIRWISE | CTOOLS_DOF_FIRSTPAIR))) continue;

        for (ST_int h = 0; h < g; h++) {
            if (nested && nested[h]) continue;
            if ((flags & CTOOLS_DOF_FIRSTPAIR) && pairs) break;
            ST_int groups = ctools_count_connected_components(
                factors[h].levels, factors[g].levels,
                N, factors[h].num_levels, factors[g].num_levels);
            if (groups < 0) return -1;
            pairs++;
            if (groups > redundant[g]) redundant[g] = groups;
        }
    }
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
            free(state->factors[g].slope);
            free(state->factors[g].slope_center);
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

/* Retain the original double partition before producing integer IDs. */
typedef struct { ST_double value; ST_int index; } ctools_numeric_cluster_pair;
static int ctools_compare_numeric_cluster(const void *a, const void *b)
{
    ST_double x = ((const ctools_numeric_cluster_pair *)a)->value;
    ST_double y = ((const ctools_numeric_cluster_pair *)b)->value;
    return (x > y) - (x < y);
}

ST_int ctools_numeric_to_cluster_ids(const ST_double *values, ST_int N,
                                    ST_int *ids, ST_int *num_groups)
{
    if (N < 0 || !values || !ids || !num_groups) return -1;
    *num_groups = 0;
    if (N == 0) return 0;

    /* Integer codes in a compact range (the usual cluster variable) are
     * ranked by counting: same dense IDs in ascending value order as the
     * sort below, in O(N + range). */
    {
        ST_double lo = 0.0, hi = 0.0;
        int compact = 1, any = 0;
        for (ST_int i = 0; i < N && compact; i++) {
            ST_double v = values[i];
            if (SF_is_missing(v) || !isfinite(v)) continue;
            if (v != floor(v) || fabs(v) > 4.0e15) compact = 0;
            else if (!any) { lo = hi = v; any = 1; }
            else { if (v < lo) lo = v; if (v > hi) hi = v; }
        }
        if (compact && any && hi - lo < 4.0 * (ST_double)N + 1024.0) {
            size_t range = (size_t)(hi - lo) + 1;
            ST_int *rank = (ST_int *)calloc(range, sizeof(ST_int));
            if (rank) {
                for (ST_int i = 0; i < N; i++) {
                    ST_double v = values[i];
                    if (!SF_is_missing(v) && isfinite(v)) rank[(size_t)(v - lo)] = 1;
                }
                ST_int groups = 0;
                for (size_t r = 0; r < range; r++) {
                    if (rank[r]) rank[r] = ++groups;
                }
                for (ST_int i = 0; i < N; i++) {
                    ST_double v = values[i];
                    ids[i] = (!SF_is_missing(v) && isfinite(v)) ? rank[(size_t)(v - lo)] - 1 : -1;
                }
                free(rank);
                *num_groups = groups;
                return 0;
            }
        }
        if (compact && !any) {
            for (ST_int i = 0; i < N; i++) ids[i] = -1;
            return 0;
        }
    }

    ctools_numeric_cluster_pair *pairs = malloc((size_t)N * sizeof(*pairs));
    if (!pairs) return -1;
    ST_int count = 0;
    for (ST_int i = 0; i < N; i++) {
        ids[i] = -1;
        if (!SF_is_missing(values[i]) && isfinite(values[i])) {
            pairs[count].value = values[i];
            pairs[count++].index = i;
        }
    }
    qsort(pairs, (size_t)count, sizeof(*pairs), ctools_compare_numeric_cluster);
    ST_int group = -1;
    for (ST_int i = 0; i < count; i++) {
        if (i == 0 || pairs[i].value != pairs[i - 1].value) group++;
        ids[pairs[i].index] = group;
    }
    *num_groups = group + 1;
    free(pairs);
    return 0;
}

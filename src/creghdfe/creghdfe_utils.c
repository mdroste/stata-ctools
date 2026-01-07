/*
 * creghdfe_utils.c
 *
 * Utility functions: timing, hash tables, union-find
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_utils.h"
#include "ctools_timer.h"

/* ========================================================================
 * High-resolution timing for profiling
 * Uses shared ctools_timer module for cross-platform compatibility.
 * ======================================================================== */

double get_time_sec(void) {
    return ctools_timer_seconds();
}

/* ========================================================================
 * Hash table implementation
 * ======================================================================== */

static ST_int hash_int(ST_int key, ST_int capacity)
{
    /* Simple multiplicative hash */
    unsigned int h = (unsigned int)key * 2654435761u;
    return (ST_int)(h % (unsigned int)capacity);
}

IntHashTable *hash_create(ST_int expected_size)
{
    IntHashTable *ht = (IntHashTable *)malloc(sizeof(IntHashTable));
    if (!ht) return NULL;

    /* Size to maintain load factor */
    ht->capacity = (ST_int)(expected_size / HASH_LOAD_FACTOR) + 1;
    if (ht->capacity < 16) ht->capacity = 16;

    ht->keys = (ST_int *)malloc(ht->capacity * sizeof(ST_int));
    ht->values = (ST_int *)malloc(ht->capacity * sizeof(ST_int));
    ht->size = 0;

    if (!ht->keys || !ht->values) {
        if (ht->keys) free(ht->keys);
        if (ht->values) free(ht->values);
        free(ht);
        return NULL;
    }

    /* Initialize all slots as empty */
    for (ST_int i = 0; i < ht->capacity; i++) {
        ht->keys[i] = HASH_EMPTY;
    }

    return ht;
}

void hash_destroy(IntHashTable *ht)
{
    if (ht) {
        if (ht->keys) free(ht->keys);
        if (ht->values) free(ht->values);
        free(ht);
    }
}

/* Insert or get existing level for a key. Returns the level (1-indexed) */
ST_int hash_insert_or_get(IntHashTable *ht, ST_int key)
{
    ST_int idx = hash_int(key, ht->capacity);

    /* Linear probing */
    while (ht->keys[idx] != HASH_EMPTY) {
        if (ht->keys[idx] == key) {
            return ht->values[idx];  /* Already exists */
        }
        idx = (idx + 1) % ht->capacity;
    }

    /* Insert new key */
    ht->keys[idx] = key;
    ht->values[idx] = ++ht->size;  /* 1-indexed level */
    return ht->size;
}

/* ========================================================================
 * Union-Find implementation
 * ======================================================================== */

UnionFind *uf_create(ST_int size)
{
    UnionFind *uf = (UnionFind *)malloc(sizeof(UnionFind));
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

    /* Initialize: each element is its own parent */
    for (ST_int i = 0; i < size; i++) {
        uf->parent[i] = i;
    }

    return uf;
}

void uf_destroy(UnionFind *uf)
{
    if (uf) {
        if (uf->parent) free(uf->parent);
        if (uf->rank) free(uf->rank);
        free(uf);
    }
}

/* Find with path compression */
ST_int uf_find(UnionFind *uf, ST_int x)
{
    if (uf->parent[x] != x) {
        uf->parent[x] = uf_find(uf, uf->parent[x]);
    }
    return uf->parent[x];
}

/* Union by rank */
void uf_union(UnionFind *uf, ST_int x, ST_int y)
{
    ST_int root_x = uf_find(uf, x);
    ST_int root_y = uf_find(uf, y);

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

/* Count connected components in a bipartite graph
 *
 * Given two FE level vectors of length N, where:
 *   - fe1_levels[i] is in [1, num_levels1]
 *   - fe2_levels[i] is in [1, num_levels2]
 *
 * Each observation (i) creates an edge between fe1_levels[i] and fe2_levels[i]
 * in a bipartite graph. We count the number of connected components.
 *
 * We use Union-Find on a combined node space:
 *   - Nodes 0 to num_levels1-1 represent FE1 levels
 *   - Nodes num_levels1 to num_levels1+num_levels2-1 represent FE2 levels
 *
 * Returns: number of connected components (mobility groups)
 */
ST_int count_connected_components(
    const ST_int *fe1_levels,
    const ST_int *fe2_levels,
    ST_int N,
    ST_int num_levels1,
    ST_int num_levels2
)
{
    ST_int total_nodes = num_levels1 + num_levels2;
    UnionFind *uf = uf_create(total_nodes);
    if (!uf) return -1;

    ST_int i;

    /* For each observation, union the FE1 and FE2 nodes */
    for (i = 0; i < N; i++) {
        ST_int node1 = fe1_levels[i] - 1;                   /* Convert to 0-indexed */
        ST_int node2 = num_levels1 + fe2_levels[i] - 1;     /* Offset for FE2 */
        uf_union(uf, node1, node2);
    }

    /* Count unique roots among FE1 nodes (this gives us the number of components) */
    /* We only need to count among one side of the bipartite graph */
    ST_int num_components = 0;
    ST_int *seen_roots = (ST_int *)calloc(total_nodes, sizeof(ST_int));
    if (!seen_roots) {
        uf_destroy(uf);
        return -1;
    }

    for (i = 0; i < num_levels1; i++) {
        ST_int root = uf_find(uf, i);
        if (!seen_roots[root]) {
            seen_roots[root] = 1;
            num_components++;
        }
    }

    free(seen_roots);
    uf_destroy(uf);

    return num_components;
}

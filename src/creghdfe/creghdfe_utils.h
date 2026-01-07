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

/* ========================================================================
 * Union-Find functions
 * ======================================================================== */

UnionFind *uf_create(ST_int size);
void uf_destroy(UnionFind *uf);
ST_int uf_find(UnionFind *uf, ST_int x);
void uf_union(UnionFind *uf, ST_int x, ST_int y);

/* Count connected components in a bipartite graph */
ST_int count_connected_components(
    const ST_int *fe1_levels,
    const ST_int *fe2_levels,
    ST_int N,
    ST_int num_levels1,
    ST_int num_levels2
);

#endif /* CREGHDFE_UTILS_H */

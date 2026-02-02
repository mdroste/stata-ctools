/*
 * ctools_sort_pairs.h
 * Order pair sorting for ctools
 *
 * Provides parallel LSD radix sort for (order_key, orig_idx) pairs.
 * Used by cmerge and other commands that need stable ordering.
 *
 * The sort is O(n) for bounded integer keys with parallel execution.
 * Falls back to qsort for small arrays where parallel overhead dominates.
 */

#ifndef CTOOLS_SORT_PAIRS_H
#define CTOOLS_SORT_PAIRS_H

#include <stdlib.h>
#include <stddef.h>

/*
 * Order pair for sorting with original index preservation.
 * Enables stable sort by comparing orig_idx when order_keys are equal.
 */
typedef struct {
    size_t order_key;    /* Primary sort key */
    size_t orig_idx;    /* Original array index (for stability) */
} ctools_order_pair_t;

/*
 * Sort an array of order pairs by order_key using parallel LSD radix sort.
 * Falls back to qsort for small arrays (n < 10000) where parallel overhead
 * would dominate.
 *
 * Sort is stable: pairs with equal order_key values are ordered by orig_idx.
 *
 * @param pairs  Array of order pairs to sort (modified in place)
 * @param n      Number of pairs
 */
void ctools_sort_order_pairs(ctools_order_pair_t *pairs, size_t n);

/*
 * Comparison function for qsort of order pairs.
 * Compares by order_key first, then orig_idx for stability.
 *
 * @param a  Pointer to first ctools_order_pair_t
 * @param b  Pointer to second ctools_order_pair_t
 * @return   -1 if a < b, 1 if a > b, 0 if equal
 */
int ctools_compare_order_pairs(const void *a, const void *b);

/* ============================================================================
 * Compatibility layer for cmerge
 *
 * These typedefs and macros allow cmerge to continue using its existing
 * naming conventions without modification.
 * ============================================================================ */

/* Type compatibility */
typedef ctools_order_pair_t cmerge_order_pair_t;

/* Function compatibility */
#define cmerge_radix_sort_order_pairs ctools_sort_order_pairs
#define cmerge_compare_order_pairs ctools_compare_order_pairs

#endif /* CTOOLS_SORT_PAIRS_H */

/*
 * cmerge_radix_sort.h
 * Parallel LSD Radix Sort for order pairs
 *
 * O(n) stable sort optimized for bounded integer keys with parallel execution.
 * Uses 4 passes of 8-bit radix sort (sufficient for Stata's 2^31 obs limit).
 */

#ifndef CMERGE_RADIX_SORT_H
#define CMERGE_RADIX_SORT_H

#include <stdlib.h>
#include <stdint.h>

/* Order pair for preserve_order sorting */
typedef struct {
    size_t order_key;
    size_t orig_idx;
} cmerge_order_pair_t;

/*
 * Sort an array of order pairs by order_key using parallel LSD radix sort.
 * Falls back to qsort for small arrays where parallel overhead dominates.
 *
 * @param pairs  Array of order pairs to sort (modified in place)
 * @param n      Number of pairs
 */
void cmerge_radix_sort_order_pairs(cmerge_order_pair_t *pairs, size_t n);

/* Comparison function for qsort fallback */
int cmerge_compare_order_pairs(const void *a, const void *b);

#endif /* CMERGE_RADIX_SORT_H */

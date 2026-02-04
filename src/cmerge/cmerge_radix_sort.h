/*
 * cmerge_radix_sort.h
 * Parallel LSD Radix Sort for order pairs
 *
 * This header now wraps ctools_sort_pairs.h for backwards compatibility.
 * All functionality has been consolidated into the unified sort pairs module.
 */

#ifndef CMERGE_RADIX_SORT_H
#define CMERGE_RADIX_SORT_H

/* Include the unified sort pairs header which provides:
 * - cmerge_order_pair_t (typedef for ctools_order_pair_t)
 * - cmerge_radix_sort_order_pairs (macro for ctools_sort_order_pairs)
 * - cmerge_compare_order_pairs (macro for ctools_compare_order_pairs)
 */
#include "../ctools_sort_pairs.h"

#endif /* CMERGE_RADIX_SORT_H */

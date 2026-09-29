/*
 * cmerge_group_search.h
 * Group boundary detection for sorted merge operations
 *
 * Provides hybrid linear/binary search for finding the end of groups
 * of matching keys efficiently.
 */

#ifndef CMERGE_GROUP_SEARCH_H
#define CMERGE_GROUP_SEARCH_H

#include <stdlib.h>
#include "../ctools_types.h"

/*
 * Find the end of a group of matching keys starting at 'start'.
 * Uses a hybrid approach:
 * - Linear scan when groups are likely small
 * - Binary search when groups appear large
 *
 * @param data        Dataset containing key variables
 * @param start       Starting index of the group
 * @param nobs        Total number of observations
 * @param nkeys       Number of key variables
 * @param all_numeric True if all keys are numeric (enables fast path)
 * @return            First index after the group (index where keys differ)
 */
size_t cmerge_find_group_end_hybrid(stata_data *data, size_t start,
                                     size_t nobs, int nkeys, int all_numeric);

#endif /* CMERGE_GROUP_SEARCH_H */

#ifndef CTOOLS_GROUP_H
#define CTOOLS_GROUP_H
#include "ctools_data.h"
/* Runs refer to positions in order (zero based), or physical rows if NULL.
 * keys are zero-based columns, or NULL for the first nkeys columns.
 * Caller provides nobs ranges. No allocation failure escapes: the parallel
 * boundary bitmap falls back to serial construction. Exact numeric equality,
 * including distinct Stata missings; NULL string and empty string are equal. */
typedef struct
{
    size_t start, count;
} ctools_group_range;
void ctools_group_ranges(const stata_data *data, const perm_idx_t *order, const int *keys,
                         size_t nkeys, ctools_group_range *groups, size_t *ngroups);
#endif

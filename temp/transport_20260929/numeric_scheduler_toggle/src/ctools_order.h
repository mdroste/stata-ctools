#ifndef CTOOLS_ORDER_H
#define CTOOLS_ORDER_H
#include "ctools_types.h"
/* Exact lexicographic comparison, including distinct extended missing values.
 * Column indices are zero based. Strings compare by their stored bytes. */
int ctools_compare_rows(const stata_data *a, size_t i, const int *ak,
                       const stata_data *b, size_t j, const int *bk, size_t nk);
/* Stable index sort. Does not move data or coalesce extended missing keys. */
stata_retcode ctools_order_stable(stata_data *data, const int *keys, size_t nkeys);
/* Split sorted rows into runs of equal keys, like Stata's by-groups: rows are
 * equal when ctools_compare_rows returns 0 (exact strings, distinct missing
 * codes, -0 == +0). order lists the rows in sorted order, or is NULL when the
 * rows of d are already sorted. Writes the first position of each run to
 * starts[0..n) and starts[n] = d->nobs; starts must hold d->nobs + 1 entries.
 * Returns the number of runs n (0 when d has no rows). */
size_t ctools_group_starts(const stata_data *d, const perm_idx_t *order,
                           const int *keys, size_t nkeys, size_t *starts);
#endif

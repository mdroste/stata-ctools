/*
 * cmerge_group_search.c
 * Group boundary detection for sorted merge operations
 */

#include "cmerge_group_search.h"
#include "cmerge_keys.h"

/* Linear scan threshold: below this, always use linear scan */
#define LINEAR_SCAN_THRESHOLD 16

size_t cmerge_find_group_end_hybrid(stata_data *data, size_t start,
                                     size_t nobs, int nkeys, int all_numeric)
{
    /* For small remaining range, just use linear scan */
    if (nobs - start <= LINEAR_SCAN_THRESHOLD) {
        size_t idx = start + 1;
        while (idx < nobs) {
            int cmp;
            if (all_numeric) {
                cmp = cmerge_compare_keys_numeric_same(data, start, idx, nkeys);
            } else {
                cmp = cmerge_compare_keys_same(data, start, idx, nkeys);
            }
            if (cmp != 0) break;
            idx++;
        }
        return idx;
    }

    /* Try linear scan first (optimistic: small groups are common) */
    size_t probe_end = start + LINEAR_SCAN_THRESHOLD;
    if (probe_end > nobs) probe_end = nobs;

    size_t idx = start + 1;
    while (idx < probe_end) {
        int cmp;
        if (all_numeric) {
            cmp = cmerge_compare_keys_numeric_same(data, start, idx, nkeys);
        } else {
            cmp = cmerge_compare_keys_same(data, start, idx, nkeys);
        }
        if (cmp != 0) return idx;
        idx++;
    }

    /* Group extends past linear probe, switch to binary search */
    /* lo is the last known matching index, hi is the first unknown */
    size_t lo = probe_end - 1;  /* Last index we know matches */
    size_t hi = nobs;           /* First index beyond matches */

    while (lo + 1 < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int cmp;
        if (all_numeric) {
            cmp = cmerge_compare_keys_numeric_same(data, start, mid, nkeys);
        } else {
            cmp = cmerge_compare_keys_same(data, start, mid, nkeys);
        }

        if (cmp == 0) {
            lo = mid;  /* Key matches, extend known-match range */
        } else {
            hi = mid;  /* Key differs, narrow search range */
        }
    }

    /* lo is the last known matching index, hi is the first known differing index.
     * Return the first index after the matching group. */
    return lo + 1;
}

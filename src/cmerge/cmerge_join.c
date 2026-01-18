/*
 * cmerge_join.c
 * Sorted merge join algorithm for cmerge
 */

#include <stdlib.h>
#include <string.h>
#include "cmerge_join.h"
#include "cmerge_keys.h"
#include "cmerge_group_search.h"

int64_t cmerge_sorted_join(
    stata_data *master_keys, stata_data *using_keys,
    int nkeys, cmerge_type_t merge_type,
    cmerge_output_spec_t **output_specs_out)
{
    size_t m_nobs = master_keys->nobs;
    size_t u_nobs = using_keys->nobs;

    /* Determine if we can use the numeric fast path */
    int master_numeric = cmerge_all_keys_numeric(master_keys, nkeys);
    int using_numeric = cmerge_all_keys_numeric(using_keys, nkeys);
    int all_numeric = master_numeric && using_numeric;

    /* OPTIMIZATION 1: Better initial capacity estimation
       Based on merge type and dataset sizes to minimize reallocations */
    size_t capacity;
    switch (merge_type) {
        case MERGE_1_1:
            /* Max possible: all unique = m + u, typical: max(m,u) */
            capacity = m_nobs + u_nobs;
            break;
        case MERGE_M_1:
            /* Output = master rows (each gets one using match or none) */
            capacity = m_nobs + u_nobs;
            break;
        case MERGE_1_M:
            /* Output can expand: each master can match multiple using */
            capacity = m_nobs + u_nobs + u_nobs / 2;
            break;
        case MERGE_M_M:
        default:
            /* Sequential pairing with expansion: output count = max of group counts */
            capacity = m_nobs + u_nobs;
            break;
    }
    /* Round up to power of 2 for better realloc behavior */
    size_t pow2 = 1024;
    while (pow2 < capacity) pow2 *= 2;
    capacity = pow2;

    cmerge_output_spec_t *specs = malloc(capacity * sizeof(cmerge_output_spec_t));
    if (!specs) return -1;

    size_t count = 0;
    size_t m_idx = 0;
    size_t u_idx = 0;

    /* Macro to ensure capacity - reduces code duplication */
    #define ENSURE_CAPACITY(needed) \
        do { \
            if (count + (needed) > capacity) { \
                while (capacity < count + (needed)) capacity *= 2; \
                cmerge_output_spec_t *new_specs = realloc(specs, capacity * sizeof(cmerge_output_spec_t)); \
                if (!new_specs) { free(specs); return -1; } \
                specs = new_specs; \
            } \
        } while (0)

    while (m_idx < m_nobs || u_idx < u_nobs) {
        if (m_idx >= m_nobs) {
            /* Master exhausted - bulk add remaining using rows */
            size_t remaining = u_nobs - u_idx;
            ENSURE_CAPACITY(remaining);
            for (size_t i = 0; i < remaining; i++) {
                specs[count].master_sorted_row = -1;
                specs[count].using_sorted_row = (int64_t)(u_idx + i);
                specs[count].merge_result = MERGE_RESULT_USING_ONLY;
                count++;
            }
            u_idx = u_nobs;
        }
        else if (u_idx >= u_nobs) {
            /* Using exhausted - bulk add remaining master rows */
            size_t remaining = m_nobs - m_idx;
            ENSURE_CAPACITY(remaining);
            for (size_t i = 0; i < remaining; i++) {
                specs[count].master_sorted_row = (int64_t)(m_idx + i);
                specs[count].using_sorted_row = -1;
                specs[count].merge_result = MERGE_RESULT_MASTER_ONLY;
                count++;
            }
            m_idx = m_nobs;
        }
        else {
            /* OPTIMIZATION 2: Use specialized numeric comparator when possible */
            int cmp;
            if (all_numeric) {
                cmp = cmerge_compare_keys_numeric(master_keys, m_idx, using_keys, u_idx, nkeys);
            } else {
                cmp = cmerge_compare_keys(master_keys, m_idx, using_keys, u_idx, nkeys);
            }

            if (cmp < 0) {
                /* Master only */
                ENSURE_CAPACITY(1);
                specs[count].master_sorted_row = (int64_t)m_idx;
                specs[count].using_sorted_row = -1;
                specs[count].merge_result = MERGE_RESULT_MASTER_ONLY;
                count++;
                m_idx++;
            }
            else if (cmp > 0) {
                /* Using only */
                ENSURE_CAPACITY(1);
                specs[count].master_sorted_row = -1;
                specs[count].using_sorted_row = (int64_t)u_idx;
                specs[count].merge_result = MERGE_RESULT_USING_ONLY;
                count++;
                u_idx++;
            }
            else {
                /* Keys match - find group boundaries */
                size_t m_start = m_idx;
                size_t u_start = u_idx;

                /* OPTIMIZATION 3: Use binary search for group boundaries
                   instead of linear scan when groups might be large */
                size_t m_end = cmerge_find_group_end_hybrid(master_keys, m_start,
                                                            m_nobs, nkeys, master_numeric);
                size_t u_end = cmerge_find_group_end_hybrid(using_keys, u_start,
                                                            u_nobs, nkeys, using_numeric);

                size_t m_count = m_end - m_start;
                size_t u_count = u_end - u_start;

                /* Pre-calculate output size for this group */
                size_t group_output;
                switch (merge_type) {
                    case MERGE_1_1:
                        group_output = 1;
                        break;
                    case MERGE_M_1:
                        group_output = m_count;
                        break;
                    case MERGE_1_M:
                        group_output = u_count;
                        break;
                    case MERGE_M_M:
                    default:
                        group_output = (m_count > u_count) ? m_count : u_count;  /* Max for sequential pairing */
                        break;
                }
                ENSURE_CAPACITY(group_output);

                /* Generate output based on merge type */
                switch (merge_type) {
                    case MERGE_1_1:
                        specs[count].master_sorted_row = (int64_t)m_start;
                        specs[count].using_sorted_row = (int64_t)u_start;
                        specs[count].merge_result = MERGE_RESULT_BOTH;
                        count++;
                        break;

                    case MERGE_M_1:
                        /* Multiple master rows map to one using row */
                        for (size_t i = 0; i < m_count; i++) {
                            specs[count].master_sorted_row = (int64_t)(m_start + i);
                            specs[count].using_sorted_row = (int64_t)u_start;
                            specs[count].merge_result = MERGE_RESULT_BOTH;
                            count++;
                        }
                        break;

                    case MERGE_1_M:
                        /* One master row maps to multiple using rows */
                        for (size_t i = 0; i < u_count; i++) {
                            specs[count].master_sorted_row = (int64_t)m_start;
                            specs[count].using_sorted_row = (int64_t)(u_start + i);
                            specs[count].merge_result = MERGE_RESULT_BOTH;
                            count++;
                        }
                        break;

                    case MERGE_M_M:
                        /* Sequential pairing: pair row 1 with row 1, row 2 with row 2, etc.
                         * If one side has more rows, excess rows get the last row of the other side.
                         * This matches Stata's m:m merge behavior. */
                        {
                            size_t max_count = (m_count > u_count) ? m_count : u_count;
                            for (size_t i = 0; i < max_count; i++) {
                                /* For master: use row i if available, else use last row */
                                size_t mi = (i < m_count) ? i : (m_count - 1);
                                /* For using: use row i if available, else use last row */
                                size_t ui = (i < u_count) ? i : (u_count - 1);

                                specs[count].master_sorted_row = (int64_t)(m_start + mi);
                                specs[count].using_sorted_row = (int64_t)(u_start + ui);
                                specs[count].merge_result = MERGE_RESULT_BOTH;
                                count++;
                            }
                        }
                        break;
                }

                m_idx = m_end;
                u_idx = u_end;
            }
        }
    }

    #undef ENSURE_CAPACITY

    *output_specs_out = specs;
    return (int64_t)count;
}

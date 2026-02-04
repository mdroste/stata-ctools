/*
 * cmerge_join.c
 * Sorted merge join algorithm for cmerge
 */

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "cmerge_join.h"
#include "cmerge_keys.h"
#include "cmerge_group_search.h"
#include "../ctools_config.h"  /* For ctools_safe_mul_size */

/* Prefetch hint macro for read-ahead optimization */
#if defined(__GNUC__) || defined(__clang__)
    /* Prefetch for read (0), low temporal locality (0) - data won't be reused soon */
    #define CMERGE_PREFETCH(addr) __builtin_prefetch((addr), 0, 0)
    /* Prefetch distance: how many rows ahead to prefetch */
    #define CMERGE_PREFETCH_DISTANCE 32
#else
    #define CMERGE_PREFETCH(addr) ((void)0)
    #define CMERGE_PREFETCH_DISTANCE 0
#endif

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

    /* Overflow-safe capacity calculation */
    /* Maximum safe capacity to avoid overflow in size calculations */
    const size_t MAX_SAFE_CAPACITY = SIZE_MAX / sizeof(cmerge_output_spec_t) / 2;

    switch (merge_type) {
        case MERGE_1_1:
        case MERGE_M_1:
        case MERGE_M_M:
        default:
            /* Max possible: m + u - check for overflow */
            if (m_nobs > SIZE_MAX - u_nobs) {
                /* Overflow would occur - cap at safe maximum */
                capacity = MAX_SAFE_CAPACITY;
            } else {
                capacity = m_nobs + u_nobs;
            }
            break;
        case MERGE_1_M:
            /* Output can expand: m + u + u/2 - check for overflow */
            if (m_nobs > SIZE_MAX - u_nobs) {
                capacity = MAX_SAFE_CAPACITY;
            } else {
                size_t base = m_nobs + u_nobs;
                size_t extra = u_nobs / 2;
                if (base > SIZE_MAX - extra) {
                    capacity = MAX_SAFE_CAPACITY;
                } else {
                    capacity = base + extra;
                }
            }
            break;
    }

    /* Cap capacity at safe maximum */
    if (capacity > MAX_SAFE_CAPACITY) {
        capacity = MAX_SAFE_CAPACITY;
    }

    /* Round up to power of 2 for better realloc behavior - with overflow protection */
    size_t pow2 = 1024;
    while (pow2 < capacity && pow2 <= SIZE_MAX / 2) {
        pow2 *= 2;
    }
    /* If pow2 would overflow, use capacity directly */
    if (pow2 >= capacity) {
        capacity = pow2;
    }

    /* Safe allocation size calculation */
    size_t alloc_size;
    if (ctools_safe_mul_size(capacity, sizeof(cmerge_output_spec_t), &alloc_size) != 0) {
        return -1;  /* Overflow in allocation size */
    }
    cmerge_output_spec_t *specs = malloc(alloc_size);
    if (!specs) return -1;

    size_t count = 0;
    size_t m_idx = 0;
    size_t u_idx = 0;

    /* Macro to ensure capacity - with overflow protection */
    #define ENSURE_CAPACITY(needed) \
        do { \
            size_t _needed = (needed); \
            /* Check for overflow in count + needed */ \
            if (_needed > SIZE_MAX - count) { \
                free(specs); return -1; /* Overflow */ \
            } \
            if (count + _needed > capacity) { \
                size_t new_capacity = capacity; \
                size_t target = count + _needed; \
                /* Double capacity until sufficient, with overflow check */ \
                while (new_capacity < target) { \
                    if (new_capacity > SIZE_MAX / 2) { \
                        /* Can't double - try exact fit */ \
                        new_capacity = target; \
                        break; \
                    } \
                    new_capacity *= 2; \
                } \
                /* Check allocation size for overflow */ \
                if (new_capacity > SIZE_MAX / sizeof(cmerge_output_spec_t)) { \
                    free(specs); return -1; /* Allocation would overflow */ \
                } \
                cmerge_output_spec_t *new_specs = realloc(specs, new_capacity * sizeof(cmerge_output_spec_t)); \
                if (!new_specs) { free(specs); return -1; } \
                specs = new_specs; \
                capacity = new_capacity; \
            } \
        } while (0)

    while (m_idx < m_nobs || u_idx < u_nobs) {
        /* OPTIMIZATION: Prefetch ahead for better cache utilization */
        #if CMERGE_PREFETCH_DISTANCE > 0
        if (all_numeric && nkeys > 0) {
            /* Prefetch master key data ahead */
            if (m_idx + CMERGE_PREFETCH_DISTANCE < m_nobs) {
                CMERGE_PREFETCH(&master_keys->vars[0].data.dbl[m_idx + CMERGE_PREFETCH_DISTANCE]);
            }
            /* Prefetch using key data ahead */
            if (u_idx + CMERGE_PREFETCH_DISTANCE < u_nobs) {
                CMERGE_PREFETCH(&using_keys->vars[0].data.dbl[u_idx + CMERGE_PREFETCH_DISTANCE]);
            }
        }
        #endif

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

/*
 * cmerge_join.h
 * Sorted merge join algorithm for cmerge
 *
 * Implements the core merge join logic for 1:1, m:1, 1:m, and m:m merge types.
 */

#ifndef CMERGE_JOIN_H
#define CMERGE_JOIN_H

#include <stdlib.h>
#include <stdint.h>
#include "../ctools_types.h"

/* Merge types */
typedef enum {
    MERGE_1_1 = 0,
    MERGE_M_1 = 1,
    MERGE_1_M = 2,
    MERGE_M_M = 3
} cmerge_type_t;

/* Merge result codes */
typedef enum {
    MERGE_RESULT_MASTER_ONLY = 1,
    MERGE_RESULT_USING_ONLY = 2,
    MERGE_RESULT_BOTH = 3
} cmerge_result_t;

/* Output specification for each row in merged output.
 * Uses int32_t for row indices since Stata max obs is INT32_MAX.
 * This halves the struct size (24 -> 12 bytes) for better cache density. */
typedef struct {
    int32_t master_sorted_row;  /* Row index in sorted master (-1 if using-only) */
    int32_t using_sorted_row;   /* Row index in sorted using (-1 if master-only) */
    int8_t merge_result;        /* 1=master only, 2=using only, 3=matched */
} cmerge_output_spec_t;

/*
 * Perform a sorted merge join between master and using datasets.
 *
 * @param master_keys      Sorted master key data
 * @param using_keys       Sorted using key data
 * @param nkeys            Number of key variables
 * @param merge_type       Type of merge (1:1, m:1, 1:m, m:m)
 * @param output_specs_out Pointer to receive allocated output specs array
 * @return                 Number of output rows, or -1 on error
 */
int64_t cmerge_sorted_join(
    stata_data *master_keys, stata_data *using_keys,
    int nkeys, cmerge_type_t merge_type,
    cmerge_output_spec_t **output_specs_out);

#endif /* CMERGE_JOIN_H */

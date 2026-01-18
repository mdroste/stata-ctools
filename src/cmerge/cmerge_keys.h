/*
 * cmerge_keys.h
 * Key comparison functions for cmerge
 *
 * Provides optimized key comparison with specialized numeric fast path
 * and general path for mixed string/numeric keys.
 */

#ifndef CMERGE_KEYS_H
#define CMERGE_KEYS_H

#include "../ctools_types.h"

/* Check if all keys are numeric (for fast path selection) */
int cmerge_all_keys_numeric(stata_data *data, int nkeys);

/* Fast path: compare all-numeric keys without type checking in loop */
int cmerge_compare_keys_numeric(stata_data *data_a, size_t row_a,
                                 stata_data *data_b, size_t row_b,
                                 int nkeys);

/* Fast path: compare same-dataset numeric keys */
int cmerge_compare_keys_numeric_same(stata_data *data, size_t row_a,
                                      size_t row_b, int nkeys);

/* General path: handles mixed string/numeric keys */
int cmerge_compare_keys(stata_data *data_a, size_t row_a,
                         stata_data *data_b, size_t row_b,
                         int nkeys);

/* General path: same-dataset comparison */
int cmerge_compare_keys_same(stata_data *data, size_t row_a, size_t row_b, int nkeys);

#endif /* CMERGE_KEYS_H */

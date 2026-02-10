/*
 * creghdfe_utils.c
 *
 * Utility functions: timing, hash tables, union-find
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_utils.h"
#include "ctools_runtime.h"
#include "../ctools_hdfe_utils.h"

/* ========================================================================
 * High-resolution timing for profiling
 * Uses shared ctools_timer module for cross-platform compatibility.
 * ======================================================================== */

double get_time_sec(void) {
    return ctools_timer_seconds();
}

/* ========================================================================
 * Hash table implementation
 * ======================================================================== */

static ST_int hash_int(ST_int key, ST_int capacity)
{
    /* Simple multiplicative hash */
    unsigned int h = (unsigned int)key * 2654435761u;
    return (ST_int)(h % (unsigned int)capacity);
}

IntHashTable *hash_create(ST_int expected_size)
{
    IntHashTable *ht = (IntHashTable *)malloc(sizeof(IntHashTable));
    if (!ht) return NULL;

    /* Size to maintain load factor */
    ht->capacity = (ST_int)(expected_size / HASH_LOAD_FACTOR) + 1;
    if (ht->capacity < 16) ht->capacity = 16;

    ht->keys = (ST_int *)malloc(ht->capacity * sizeof(ST_int));
    ht->values = (ST_int *)malloc(ht->capacity * sizeof(ST_int));
    ht->size = 0;

    if (!ht->keys || !ht->values) {
        if (ht->keys) free(ht->keys);
        if (ht->values) free(ht->values);
        free(ht);
        return NULL;
    }

    /* Initialize all slots as empty */
    for (ST_int i = 0; i < ht->capacity; i++) {
        ht->keys[i] = HASH_EMPTY;
    }

    return ht;
}

void hash_destroy(IntHashTable *ht)
{
    if (ht) {
        if (ht->keys) free(ht->keys);
        if (ht->values) free(ht->values);
        free(ht);
    }
}

/* Insert or get existing level for a key. Returns the level (1-indexed),
 * or -1 if the table is full (prevents infinite loop on overload). */
ST_int hash_insert_or_get(IntHashTable *ht, ST_int key)
{
    ST_int idx = hash_int(key, ht->capacity);

    /* Prevent infinite loop: probe at most capacity slots */
    ST_int probes = 0;

    /* Linear probing */
    while (ht->keys[idx] != HASH_EMPTY) {
        if (ht->keys[idx] == key) {
            return ht->values[idx];  /* Already exists */
        }
        idx = (idx + 1) % ht->capacity;
        if (++probes >= ht->capacity) {
            return -1;  /* Table full — should not happen with correct sizing */
        }
    }

    /* Check load factor before inserting */
    if (ht->size >= ht->capacity - 1) {
        return -1;  /* Table full */
    }

    /* Insert new key */
    ht->keys[idx] = key;
    ht->values[idx] = ++ht->size;  /* 1-indexed level */
    return ht->size;
}

/* Delegate to shared implementation */
ST_int count_connected_components(
    const ST_int *fe1_levels,
    const ST_int *fe2_levels,
    ST_int N,
    ST_int num_levels1,
    ST_int num_levels2
)
{
    return ctools_count_connected_components(fe1_levels, fe2_levels, N, num_levels1, num_levels2);
}

/* ========================================================================
 * Sort-based value remapping
 * ======================================================================== */

/* Structure for sorting (value, index) pairs */
typedef struct {
    double value;   /* FE value (may be non-integer) */
    ST_int index;   /* Original observation index */
} ValueIndexPair;

/* Comparison function for qsort */
static int compare_value_index_inline(const ValueIndexPair *a, const ValueIndexPair *b)
{
    if (a->value < b->value) return -1;
    if (a->value > b->value) return 1;
    return 0;
}

/*
 * Remap values to contiguous levels using sort.
 * Algorithm:
 *   1. Create (value, index) pairs
 *   2. Sort by value (qsort, handles doubles including non-integer FE values)
 *   3. Linear scan: assign level 1 to first value, increment when value changes
 *   4. Write levels to output array using original indices
 */
int remap_values_sorted(const double *values, ST_int N, ST_int *levels_out, ST_int *num_levels)
{
    if (N <= 0 || !values || !levels_out || !num_levels) {
        return -1;
    }

    /* Allocate (value, index) pairs */
    ValueIndexPair *pairs = (ValueIndexPair *)malloc(N * sizeof(ValueIndexPair));
    if (!pairs) {
        return -1;
    }

    /* Fill pairs - check for missing values to avoid undefined cast behavior */
    for (ST_int i = 0; i < N; i++) {
        /* Missing values must be filtered out before calling this function */
        if (SF_is_missing(values[i])) {
            free(pairs);
            return -1;  /* Error: unexpected missing value */
        }
        pairs[i].value = values[i];
        pairs[i].index = i;
    }

    /* Sort by value - use qsort for double values */
    qsort(pairs, N, sizeof(ValueIndexPair),
          (int (*)(const void*, const void*))compare_value_index_inline);

    /* Linear scan to assign levels */
    ST_int current_level = 1;
    double prev_value = pairs[0].value;
    levels_out[pairs[0].index] = current_level;

    for (ST_int i = 1; i < N; i++) {
        if (pairs[i].value != prev_value) {
            current_level++;
            prev_value = pairs[i].value;
        }
        levels_out[pairs[i].index] = current_level;
    }

    *num_levels = (ST_int)current_level;
    free(pairs);
    return 0;
}

/*
 * Sort observations by cluster_id using counting sort.
 * Since cluster_ids are already 0-indexed and contiguous, we can use
 * counting sort which is O(N + num_clusters) instead of O(N log N).
 *
 * Algorithm:
 *   1. Count observations per cluster
 *   2. Compute prefix sums to get boundaries
 *   3. Fill sort_perm by placing each observation in its cluster's slot
 */
int sort_by_cluster(const ST_int *cluster_ids, ST_int N, ST_int num_clusters,
                    ST_int *sort_perm, ST_int *boundaries)
{
    if (N <= 0 || !cluster_ids || !sort_perm || !boundaries || num_clusters <= 0) {
        return -1;
    }

    /* Count observations per cluster */
    ST_int *counts = (ST_int *)calloc(num_clusters, sizeof(ST_int));
    if (!counts) {
        return -1;
    }

    for (ST_int i = 0; i < N; i++) {
        ST_int c = cluster_ids[i];
        if (c < 0 || c >= num_clusters) { free(counts); return -1; }
        counts[c]++;
    }

    /* Compute prefix sums to get starting positions */
    boundaries[0] = 0;
    for (ST_int c = 0; c < num_clusters; c++) {
        boundaries[c + 1] = boundaries[c] + counts[c];
    }

    /* Reset counts to use as insertion pointers */
    for (ST_int c = 0; c < num_clusters; c++) {
        counts[c] = 0;
    }

    /* Fill sort_perm */
    for (ST_int i = 0; i < N; i++) {
        ST_int c = cluster_ids[i];
        ST_int pos = boundaries[c] + counts[c];
        sort_perm[pos] = i;
        counts[c]++;
    }

    free(counts);
    return 0;
}

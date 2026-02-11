/*
 * creghdfe_utils.c
 *
 * Utility functions: timing, hash tables, union-find
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_utils.h"
#include "ctools_runtime.h"
#include "../ctools_hdfe_utils.h"
#include <stdint.h>

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

/* ========================================================================
 * Optimized FE remap with fused counting
 * ======================================================================== */

/*
 * Counting-based remap for integer FE values with manageable range: O(N + range).
 * Two-pass: (1) assign levels via direct-mapped array, (2) count with exact allocation.
 */
static int remap_counting_impl(const double *values, ST_int N, ST_int *levels_out,
                                ST_int *num_levels, ST_int **counts_out,
                                const double *weights, ST_double **wcounts_out,
                                int64_t vmin, int64_t range)
{
    /* Remap table: maps (value - vmin) -> 1-based level, 0 = unseen */
    ST_int *remap = (ST_int *)calloc((size_t)range, sizeof(ST_int));
    if (!remap) return -1;

    /* Pass 1: Assign contiguous 1-based levels */
    ST_int level = 0;
    for (ST_int i = 0; i < N; i++) {
        ST_int offset = (ST_int)((int64_t)values[i] - vmin);
        if (remap[offset] == 0) {
            remap[offset] = ++level;
        }
        levels_out[i] = remap[offset];
    }
    free(remap);
    *num_levels = level;

    /* Pass 2: Count per level with exact-size allocation */
    ST_int *counts = (ST_int *)calloc(level, sizeof(ST_int));
    if (!counts) return -1;

    if (weights && wcounts_out) {
        ST_double *wcounts = (ST_double *)calloc(level, sizeof(ST_double));
        if (!wcounts) { free(counts); return -1; }
        for (ST_int i = 0; i < N; i++) {
            ST_int lev = levels_out[i] - 1;
            counts[lev]++;
            wcounts[lev] += weights[i];
        }
        *wcounts_out = wcounts;
    } else {
        for (ST_int i = 0; i < N; i++) {
            counts[levels_out[i] - 1]++;
        }
    }

    *counts_out = counts;
    return 0;
}

/*
 * Sort-based remap fallback for non-integer or very sparse values: O(N log N).
 * Same algorithm as remap_values_sorted but without SF_is_missing check
 * (data is already filtered) and with fused counting.
 */
static int remap_sort_impl(const double *values, ST_int N, ST_int *levels_out,
                            ST_int *num_levels, ST_int **counts_out,
                            const double *weights, ST_double **wcounts_out)
{
    ValueIndexPair *pairs = (ValueIndexPair *)malloc(N * sizeof(ValueIndexPair));
    if (!pairs) return -1;

    for (ST_int i = 0; i < N; i++) {
        pairs[i].value = values[i];
        pairs[i].index = i;
    }

    qsort(pairs, N, sizeof(ValueIndexPair),
          (int (*)(const void*, const void*))compare_value_index_inline);

    ST_int level = 1;
    double prev_value = pairs[0].value;
    levels_out[pairs[0].index] = level;

    for (ST_int i = 1; i < N; i++) {
        if (pairs[i].value != prev_value) {
            level++;
            prev_value = pairs[i].value;
        }
        levels_out[pairs[i].index] = level;
    }
    free(pairs);
    *num_levels = level;

    /* Count per level */
    ST_int *counts = (ST_int *)calloc(level, sizeof(ST_int));
    if (!counts) return -1;

    if (weights && wcounts_out) {
        ST_double *wcounts = (ST_double *)calloc(level, sizeof(ST_double));
        if (!wcounts) { free(counts); return -1; }
        for (ST_int i = 0; i < N; i++) {
            ST_int lev = levels_out[i] - 1;
            counts[lev]++;
            wcounts[lev] += weights[i];
        }
        *wcounts_out = wcounts;
    } else {
        for (ST_int i = 0; i < N; i++) {
            counts[levels_out[i] - 1]++;
        }
    }

    *counts_out = counts;
    return 0;
}

int remap_and_count(const double *values, ST_int N, ST_int *levels_out,
                    ST_int *num_levels, ST_int **counts_out,
                    const double *weights, ST_double **wcounts_out)
{
    if (N <= 0 || !values || !levels_out || !num_levels || !counts_out)
        return -1;

    *counts_out = NULL;
    if (wcounts_out) *wcounts_out = NULL;

    /* Phase 1: Check if all values are integers and find min/max */
    int all_integer = 1;
    int64_t vmin = INT64_MAX, vmax = INT64_MIN;

    for (ST_int i = 0; i < N; i++) {
        double v = values[i];
        int64_t iv = (int64_t)v;
        if ((double)iv != v) {
            all_integer = 0;
            break;
        }
        if (iv < vmin) vmin = iv;
        if (iv > vmax) vmax = iv;
    }

    if (all_integer) {
        int64_t range = vmax - vmin + 1;
        /* Use counting if range is manageable: max(4*N, 16M) entries */
        int64_t threshold = 4 * (int64_t)N;
        if (threshold < 16 * 1024 * 1024) threshold = 16 * 1024 * 1024;
        if (range <= threshold) {
            return remap_counting_impl(values, N, levels_out, num_levels,
                                       counts_out, weights, wcounts_out, vmin, range);
        }
    }

    /* Fallback: sort-based approach */
    return remap_sort_impl(values, N, levels_out, num_levels,
                           counts_out, weights, wcounts_out);
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
static int sort_by_cluster(const ST_int *cluster_ids, ST_int N, ST_int num_clusters,
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

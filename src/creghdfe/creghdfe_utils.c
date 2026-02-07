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

/* Insert or get existing level for a key. Returns the level (1-indexed) */
ST_int hash_insert_or_get(IntHashTable *ht, ST_int key)
{
    ST_int idx = hash_int(key, ht->capacity);

    /* Linear probing */
    while (ht->keys[idx] != HASH_EMPTY) {
        if (ht->keys[idx] == key) {
            return ht->values[idx];  /* Already exists */
        }
        idx = (idx + 1) % ht->capacity;
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
    ST_int value;   /* Integer value (cast from double) */
    ST_int index;   /* Original observation index */
} ValueIndexPair;

/* Forward declaration for fallback */
static int compare_value_index_inline(const ValueIndexPair *a, const ValueIndexPair *b);

/*
 * LSD Radix sort for ValueIndexPair array, sorting by value field.
 * O(N × 4) for 32-bit integers, much faster than qsort's O(N log N).
 * Uses 8-bit digits (256 buckets) for 4 passes.
 */
static void radix_sort_pairs(ValueIndexPair *pairs, ST_int N)
{
    /* Allocate temp buffer for ping-pong */
    ValueIndexPair *temp = (ValueIndexPair *)malloc(N * sizeof(ValueIndexPair));
    if (!temp) {
        /* Fallback to qsort if malloc fails */
        qsort(pairs, N, sizeof(ValueIndexPair),
              (int (*)(const void*, const void*))compare_value_index_inline);
        return;
    }

    ValueIndexPair *src = pairs;
    ValueIndexPair *dst = temp;
    ST_int counts[256];

    /* 4 passes for 32-bit integers, 8 bits per pass */
    for (int pass = 0; pass < 4; pass++) {
        int shift = pass * 8;

        /* Count occurrences of each digit */
        memset(counts, 0, sizeof(counts));
        for (ST_int i = 0; i < N; i++) {
            /* Handle signed integers: flip sign bit on last pass */
            ST_int val = src[i].value;
            if (pass == 3) val ^= 0x80000000;  /* Flip sign bit for signed sort */
            int digit = (val >> shift) & 0xFF;
            counts[digit]++;
        }

        /* Convert counts to starting positions (prefix sum) */
        ST_int total = 0;
        for (int d = 0; d < 256; d++) {
            ST_int count = counts[d];
            counts[d] = total;
            total += count;
        }

        /* Distribute elements to destination */
        for (ST_int i = 0; i < N; i++) {
            ST_int val = src[i].value;
            if (pass == 3) val ^= 0x80000000;
            int digit = (val >> shift) & 0xFF;
            dst[counts[digit]++] = src[i];
        }

        /* Swap src and dst for next pass */
        ValueIndexPair *swap = src;
        src = dst;
        dst = swap;
    }

    /* If result is in temp, copy back to pairs */
    if (src != pairs) {
        memcpy(pairs, src, N * sizeof(ValueIndexPair));
    }

    free(temp);
}

/* Inline comparison for fallback qsort */
static int compare_value_index_inline(const ValueIndexPair *a, const ValueIndexPair *b)
{
    if (a->value < b->value) return -1;
    if (a->value > b->value) return 1;
    return 0;
}

/*
 * Remap values to contiguous levels using radix sort.
 * Algorithm:
 *   1. Create (value, index) pairs
 *   2. Radix sort by value - O(4N) for 32-bit integers
 *   3. Linear scan: assign level 1 to first value, increment when value changes
 *   4. Write levels to output array using original indices
 *
 * Complexity: O(N) time (4 passes), O(N) space
 * Much faster than qsort O(N log N) for large N.
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
        pairs[i].value = (ST_int)values[i];
        pairs[i].index = i;
    }

    /* Radix sort by value - O(4N) instead of O(N log N) */
    radix_sort_pairs(pairs, N);

    /* Linear scan to assign levels */
    ST_int current_level = 1;
    ST_int prev_value = pairs[0].value;
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
        counts[cluster_ids[i]]++;
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

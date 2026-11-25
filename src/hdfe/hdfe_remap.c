/*
 * hdfe_remap.c
 *
 * Utility functions: timing, hash tables, union-find
 * Part of the ctools Stata plugin suite
 */

#include "hdfe_remap.h"
#include "ctools_runtime.h"
#include "../ctools_hdfe_utils.h"
#include <stdint.h>
#include <string.h>

/* ========================================================================
 * High-resolution timing for profiling
 * Uses shared ctools_timer module for cross-platform compatibility.
 * ======================================================================== */

double get_time_sec(void)
{
    return ctools_timer_seconds();
}

/* ========================================================================
 * Connected components (delegates to the shared implementation)
 * ======================================================================== */

ST_int count_connected_components(const ST_int *fe1_levels, const ST_int *fe2_levels, ST_int N,
                                  ST_int num_levels1, ST_int num_levels2)
{
    return ctools_count_connected_components(fe1_levels, fe2_levels, N, num_levels1, num_levels2);
}

/* ========================================================================
 * Sort-based value remapping
 * ======================================================================== */

/* Structure for sorting (value, index) pairs */
typedef struct
{
    double value; /* FE value (may be non-integer) */
    ST_int index; /* Original observation index */
} ValueIndexPair;

/* Comparison function for qsort */
static int compare_value_index_inline(const ValueIndexPair *a, const ValueIndexPair *b)
{
    if (a->value < b->value)
        return -1;
    if (a->value > b->value)
        return 1;
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
    if (N <= 0 || !values || !levels_out || !num_levels)
    {
        return -1;
    }

    /* Allocate (value, index) pairs */
    ValueIndexPair *pairs = (ValueIndexPair *)malloc(N * sizeof(ValueIndexPair));
    if (!pairs)
    {
        return -1;
    }

    /* Fill pairs - check for missing values to avoid undefined cast behavior */
    for (ST_int i = 0; i < N; i++)
    {
        /* Missing values must be filtered out before calling this function */
        if (SF_is_missing(values[i]))
        {
            free(pairs);
            return -1; /* Error: unexpected missing value */
        }
        pairs[i].value = values[i];
        pairs[i].index = i;
    }

    /* Sort by value - use qsort for double values */
    qsort(pairs, N, sizeof(ValueIndexPair),
          (int (*)(const void *, const void *))compare_value_index_inline);

    /* Linear scan to assign levels */
    ST_int current_level = 1;
    double prev_value = pairs[0].value;
    levels_out[pairs[0].index] = current_level;

    for (ST_int i = 1; i < N; i++)
    {
        if (pairs[i].value != prev_value)
        {
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
                               ST_int *num_levels, ST_int **counts_out, const double *weights,
                               ST_double **wcounts_out, int64_t vmin, int64_t range)
{
    /* Remap table: maps (value - vmin) -> 1-based level, 0 = unseen */
    ST_int *remap = (ST_int *)calloc((size_t)range, sizeof(ST_int));
    if (!remap)
        return -1;

    /* Pass 1: Assign contiguous 1-based levels */
    ST_int level = 0;
    for (ST_int i = 0; i < N; i++)
    {
        size_t offset = (size_t)((int64_t)values[i] - vmin); /* range may exceed 2^31 */
        if (remap[offset] == 0)
        {
            remap[offset] = ++level;
        }
        levels_out[i] = remap[offset];
    }
    free(remap);

    /* Pass 2: Count per level with exact-size allocation */
    ST_int *counts = (ST_int *)calloc(level, sizeof(ST_int));
    if (!counts)
        return -1;

    if (weights && wcounts_out)
    {
        ST_double *wcounts = (ST_double *)calloc(level, sizeof(ST_double));
        if (!wcounts)
        {
            free(counts);
            return -1;
        }
        for (ST_int i = 0; i < N; i++)
        {
            ST_int lev = levels_out[i] - 1;
            counts[lev]++;
            wcounts[lev] += weights[i];
        }
        *wcounts_out = wcounts;
    }
    else
    {
        for (ST_int i = 0; i < N; i++)
        {
            counts[levels_out[i] - 1]++;
        }
    }

    *counts_out = counts;
    *num_levels = level;
    return 0;
}

/*
 * Hash remap for non-integer or very sparse values: O(N) expected.
 * Levels are numbered in order of first appearance, as in remap_counting_impl.
 * Equal values share a level, including signed zeros. The caller filters
 * out missing codes before remapping.
 */
static inline uint64_t remap_mix64(uint64_t x)
{
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

static int remap_hash_impl(const double *values, ST_int N, ST_int *levels_out, ST_int *num_levels,
                           ST_int **counts_out, const double *weights, ST_double **wcounts_out)
{
    size_t cap = 1u << 12;
    uint64_t *keys = (uint64_t *)malloc(cap * sizeof(uint64_t));
    ST_int *slots = (ST_int *)calloc(cap, sizeof(ST_int)); /* 0 = empty, else level */
    if (!keys || !slots)
    {
        free(keys);
        free(slots);
        return -1;
    }
    ST_int level = 0;
    for (ST_int i = 0; i < N; i++)
    {
        double v = values[i];
        if (v == 0.0)
            v = 0.0;
        uint64_t key;
        memcpy(&key, &v, sizeof key);
        size_t mask = cap - 1, h = (size_t)remap_mix64(key) & mask;
        while (slots[h] && keys[h] != key)
            h = (h + 1) & mask;
        if (slots[h])
        {
            levels_out[i] = slots[h];
            continue;
        }
        keys[h] = key;
        slots[h] = ++level;
        levels_out[i] = level;
        if ((size_t)level * 2 > cap)
        {
            /* Keep the load factor at most 1/2 */
            size_t ncap = cap * 2;
            uint64_t *nkeys = (uint64_t *)malloc(ncap * sizeof(uint64_t));
            ST_int *nslots = (ST_int *)calloc(ncap, sizeof(ST_int));
            if (!nkeys || !nslots)
            {
                free(nkeys);
                free(nslots);
                free(keys);
                free(slots);
                return -1;
            }
            for (size_t s = 0; s < cap; s++)
            {
                if (!slots[s])
                    continue;
                size_t nh = (size_t)remap_mix64(keys[s]) & (ncap - 1);
                while (nslots[nh])
                    nh = (nh + 1) & (ncap - 1);
                nkeys[nh] = keys[s];
                nslots[nh] = slots[s];
            }
            free(keys);
            free(slots);
            keys = nkeys;
            slots = nslots;
            cap = ncap;
        }
    }
    free(keys);
    free(slots);

    ST_int *counts = (ST_int *)calloc(level, sizeof(ST_int));
    if (!counts)
        return -1;
    if (weights && wcounts_out)
    {
        ST_double *wcounts = (ST_double *)calloc(level, sizeof(ST_double));
        if (!wcounts)
        {
            free(counts);
            return -1;
        }
        for (ST_int i = 0; i < N; i++)
        {
            ST_int lev = levels_out[i] - 1;
            counts[lev]++;
            wcounts[lev] += weights[i];
        }
        *wcounts_out = wcounts;
    }
    else
    {
        for (ST_int i = 0; i < N; i++)
            counts[levels_out[i] - 1]++;
    }
    *counts_out = counts;
    *num_levels = level;
    return 0;
}

int remap_and_count(const double *values, ST_int N, ST_int *levels_out, ST_int *num_levels,
                    ST_int **counts_out, const double *weights, ST_double **wcounts_out)
{
    if (num_levels)
        *num_levels = 0;
    if (counts_out)
        *counts_out = NULL;
    if (wcounts_out)
        *wcounts_out = NULL;
    if (N <= 0 || !values || !levels_out || !num_levels || !counts_out)
        return -1;

    /* Phase 1: Check if all values are integers and find min/max */
    int all_integer = 1;
    int64_t vmin = INT64_MAX, vmax = INT64_MIN;

    for (ST_int i = 0; i < N; i++)
    {
        double v = values[i];
        /* Keep counting codes within exact double integers. This also makes
         * the later signed range subtraction safe. Other labels use hashing. */
        if (!isfinite(v) || v < -0x1p53 || v > 0x1p53 || trunc(v) != v)
        {
            all_integer = 0;
            break;
        }
        int64_t iv = (int64_t)v;
        if ((double)iv != v)
        {
            all_integer = 0;
            break;
        }
        if (iv < vmin)
            vmin = iv;
        if (iv > vmax)
            vmax = iv;
    }

    if (all_integer)
    {
        int64_t range = vmax - vmin + 1;
        /* Use counting if range is manageable: max(4*N, 16M) entries */
        int64_t threshold = 4 * (int64_t)N;
        if (threshold < 16 * 1024 * 1024)
            threshold = 16 * 1024 * 1024;
        if (range <= threshold)
        {
            return remap_counting_impl(values, N, levels_out, num_levels, counts_out, weights,
                                       wcounts_out, vmin, range);
        }
    }

    return remap_hash_impl(values, N, levels_out, num_levels, counts_out, weights, wcounts_out);
}

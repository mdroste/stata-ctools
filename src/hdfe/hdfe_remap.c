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
 * FE remap: levels, then optional counts
 * ======================================================================== */

/*
 * Counting-based remap for integer FE values with manageable range: O(N + range).
 * One pass assigns levels in order of first appearance via a direct-mapped array.
 */
static int remap_counting_impl(const double *values, ST_int N, ST_int *levels_out,
                               ST_int *num_levels, int64_t vmin, int64_t range)
{
    /* Remap table: maps (value - vmin) -> 1-based level, 0 = unseen */
    ST_int *remap = (ST_int *)calloc((size_t)range, sizeof(ST_int));
    if (!remap)
        return -1;

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
    *num_levels = level;
    return 0;
}

/*
 * Integer codes in a manageable range, in B row blocks, numbered in
 * ascending code order. The values are read once: the first pass checks
 * them and keeps the low 32 bits of each code in levels_out (the offset from
 * the smallest code is then exact modulo 2^32, since the range is below
 * 2^32); the presence table is marked, ranked, and read through that copy.
 * Returns 1 (levels_out unspecified) when the codes are not such integers.
 */
static int remap_dense_parallel(const double *values, ST_int N, ST_int *levels_out,
                                ST_int *num_levels, int B)
{
    int all_integer = 1, b;
    int64_t vmin = INT64_MAX, vmax = INT64_MIN;
    uint32_t *code = (uint32_t *)levels_out;
#pragma omp parallel for schedule(static) reduction(&& : all_integer) reduction(min : vmin) \
    reduction(max : vmax)
    for (b = 0; b < B; b++)
    {
        ST_int i0, i1;
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1 && all_integer; i++)
        {
            double v = values[i];
            /* exact double integers only (the comparisons fail for NaN) */
            if (!(v >= -0x1p53 && v <= 0x1p53))
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
            code[i] = (uint32_t)iv;
        }
    }
    if (!all_integer)
        return 1;
    const int64_t range = vmax - vmin + 1;
    int64_t threshold = 4 * (int64_t)N;
    if (threshold < 16 * 1024 * 1024)
        threshold = 16 * 1024 * 1024;
    if (range > threshold || range >= ((int64_t)1 << 32))
        return 1;
    const uint32_t lo = (uint32_t)vmin;

    ST_int *table = (ST_int *)calloc((size_t)range, sizeof(ST_int));
    ST_int *base = (ST_int *)calloc((size_t)B + 1, sizeof(ST_int));
    if (!table || !base)
    {
        free(table);
        free(base);
        return -1;
    }
    /* Mark the codes present (racing stores all write 1) */
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        ST_int i0, i1;
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
        {
            ST_int *slot = &table[(uint32_t)(code[i] - lo)];
            if (!CTOOLS_HDFE_LOAD(slot))
                CTOOLS_HDFE_STORE(slot, 1);
        }
    }
    /* Rank them, the range cut into B pieces */
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        const size_t r0 = (size_t)((range * b) / B), r1 = (size_t)((range * (b + 1)) / B);
        ST_int n = 0;
        for (size_t r = r0; r < r1; r++)
            n += table[r];
        base[b + 1] = n;
    }
    for (b = 0; b < B; b++)
        base[b + 1] += base[b];
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        const size_t r0 = (size_t)((range * b) / B), r1 = (size_t)((range * (b + 1)) / B);
        ST_int next = base[b];
        for (size_t r = r0; r < r1; r++)
            if (table[r])
                table[r] = ++next;
    }
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        ST_int i0, i1;
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
            levels_out[i] = table[(uint32_t)(code[i] - lo)];
    }
    *num_levels = base[B];
    free(table);
    free(base);
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

static inline uint64_t remap_key(double v)
{
    if (v == 0.0)
        v = 0.0;
    uint64_t key;
    memcpy(&key, &v, sizeof key);
    return key;
}

static int remap_hash_impl(const double *values, ST_int N, ST_int *levels_out, ST_int *num_levels)
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
        uint64_t key = remap_key(values[i]);
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
    *num_levels = level;
    return 0;
}

/*
 * Hash remap in B row blocks through 2^bits hash partitions (by the top bits
 * of the hash): the codes are copied out partition by partition in row
 * order, each partition's distinct codes go into a table of its own (small,
 * built by one thread), numbered in order of first appearance within the
 * partition, the partitions one after another; then every row looks its
 * code up, in row order. The numbering depends on the values and N only,
 * not on the threads; it is not row order, which callers of
 * hdfe_remap_levels accept. Returns 1 when an allocation fails, so the
 * caller can fall back.
 */
typedef struct
{
    uint64_t *keys;
    ST_int *ids; /* 0 = empty, else the level within the partition */
    size_t cap;  /* power of two, at least twice the codes held */
} remap_part_table;

static int remap_part_insert(remap_part_table *t, uint64_t key, ST_int *distinct)
{
    size_t mask = t->cap - 1, h = (size_t)remap_mix64(key) & mask;
    while (t->ids[h] && t->keys[h] != key)
        h = (h + 1) & mask;
    if (t->ids[h])
        return 0;
    t->keys[h] = key;
    t->ids[h] = ++*distinct;
    if ((size_t)*distinct * 2 <= t->cap)
        return 0;
    size_t ncap = t->cap * 2;
    uint64_t *nkeys = (uint64_t *)calloc(ncap, sizeof(uint64_t));
    ST_int *nids = (ST_int *)calloc(ncap, sizeof(ST_int));
    if (!nkeys || !nids)
    {
        free(nkeys);
        free(nids);
        return -1;
    }
    for (size_t s = 0; s < t->cap; s++)
    {
        if (!t->ids[s])
            continue;
        size_t nh = (size_t)remap_mix64(t->keys[s]) & (ncap - 1);
        while (nids[nh])
            nh = (nh + 1) & (ncap - 1);
        nkeys[nh] = t->keys[s];
        nids[nh] = t->ids[s];
    }
    free(t->keys);
    free(t->ids);
    t->keys = nkeys;
    t->ids = nids;
    t->cap = ncap;
    return 0;
}

static int remap_sparse_parallel(const double *values, ST_int N, ST_int *levels_out,
                                 ST_int *num_levels, int B)
{
    int bits = 4;
    while (bits < 12 && ((size_t)16384 << bits) < (size_t)N)
        bits++;
    const size_t P = (size_t)1 << bits;
    const int shift = 64 - bits;
    size_t *pos = (size_t *)calloc((size_t)B * P, sizeof(size_t));
    size_t *start = (size_t *)malloc((P + 1) * sizeof(size_t));
    uint64_t *keys = (uint64_t *)malloc((size_t)N * sizeof(uint64_t));
    ST_int *base = (ST_int *)calloc(P + 1, sizeof(ST_int));
    remap_part_table *tab = (remap_part_table *)calloc(P, sizeof(remap_part_table));
    int failed = 0, b;
    long long p;
    const long long np = (long long)P;
    if (!pos || !start || !keys || !base || !tab)
    {
        failed = 1;
        goto done;
    }

    /* Rows of each partition per block */
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        size_t *h = pos + (size_t)b * P;
        ST_int i0, i1;
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
            h[remap_mix64(remap_key(values[i])) >> shift]++;
    }
    /* Partition-major offsets: within a partition, block b follows b - 1 */
    size_t off = 0;
    for (size_t q = 0; q < P; q++)
    {
        start[q] = off;
        for (b = 0; b < B; b++)
        {
            size_t n = pos[(size_t)b * P + q];
            pos[(size_t)b * P + q] = off;
            off += n;
        }
    }
    start[P] = off;
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        size_t *h = pos + (size_t)b * P;
        ST_int i0, i1;
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
        {
            uint64_t key = remap_key(values[i]);
            keys[h[remap_mix64(key) >> shift]++] = key;
        }
    }

    /* Each partition's distinct codes, in order of first appearance in it */
#pragma omp parallel for schedule(dynamic, 4) reduction(| : failed)
    for (p = 0; p < np; p++)
    {
        remap_part_table *t = &tab[p];
        t->cap = 64;
        t->keys = (uint64_t *)calloc(t->cap, sizeof(uint64_t));
        t->ids = (ST_int *)calloc(t->cap, sizeof(ST_int));
        if (!t->keys || !t->ids)
        {
            failed = 1;
            continue;
        }
        ST_int distinct = 0;
        for (size_t r = start[p]; r < start[p + 1]; r++)
        {
            if (remap_part_insert(t, keys[r], &distinct))
            {
                failed = 1;
                break;
            }
        }
        base[p + 1] = distinct;
    }
    if (failed)
        goto done;
    for (size_t q = 0; q < P; q++)
        base[q + 1] += base[q];

    /* Every row's level, in row order */
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        ST_int i0, i1;
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
        {
            uint64_t key = remap_key(values[i]), hash = remap_mix64(key);
            const size_t q = (size_t)(hash >> shift);
            const remap_part_table *t = &tab[q];
            size_t mask = t->cap - 1, h = (size_t)hash & mask;
            /* every code is in its table, and no empty slot precedes it */
            while (t->keys[h] != key || !t->ids[h])
                h = (h + 1) & mask;
            levels_out[i] = base[q] + t->ids[h];
        }
    }
    *num_levels = base[P];

done:
    if (tab)
        for (size_t q = 0; q < P; q++)
        {
            free(tab[q].keys);
            free(tab[q].ids);
        }
    free(tab);
    free(pos);
    free(start);
    free(keys);
    free(base);
    return failed ? 1 : 0;
}

/*
 * Levels 1..L of N values without missing codes: integer codes in a
 * manageable range by a direct table, others hashed
 * when the hash table cannot be allocated. Serially, both the table and the
 * hash number the levels in order of first appearance. With any_order the
 * call runs in parallel over row blocks (unless it is already inside a
 * parallel region or N is small), numbering integer codes in ascending
 * order and hashed codes by hash partition.
 */
static int remap_levels_impl(const double *values, ST_int N, ST_int *levels_out,
                             ST_int *num_levels, int any_order)
{
    const int B = any_order ? ctools_hdfe_row_blocks(N) : 1;
    int rc;

    if (B > 1)
    {
        rc = remap_dense_parallel(values, N, levels_out, num_levels, B);
        if (rc != 1)
            return rc;
        rc = remap_sparse_parallel(values, N, levels_out, num_levels, B);
        if (rc != 1)
            return rc;
    }
    else
    {
        /* Integer codes and their range. Keep counting codes within exact
         * double integers: this also makes the later signed range subtraction
         * safe (the comparisons fail for NaN and infinities). */
        int all_integer = 1;
        int64_t vmin = INT64_MAX, vmax = INT64_MIN;
        for (ST_int i = 0; i < N; i++)
        {
            double v = values[i];
            if (!(v >= -0x1p53 && v <= 0x1p53))
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
                return remap_counting_impl(values, N, levels_out, num_levels, vmin, range);
        }
    }

    /* Sparse or non-integer codes: hash them. A table that cannot be
     * allocated is an error; there is no other path */
    return remap_hash_impl(values, N, levels_out, num_levels);
}

int hdfe_remap_levels(const double *values, ST_int N, ST_int *levels_out, ST_int *num_levels)
{
    if (num_levels)
        *num_levels = 0;
    if (N <= 0 || !values || !levels_out || !num_levels)
        return -1;
    if (remap_levels_impl(values, N, levels_out, num_levels, 1) != 0)
    {
        *num_levels = 0;
        return -1;
    }
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
    if (N <= 0 || !values || !levels_out || !num_levels)
        return -1;

    ST_int level = 0;
    if (remap_levels_impl(values, N, levels_out, &level, 0) != 0)
        return -1;
    if (!counts_out)
    {
        *num_levels = level;
        return 0;
    }

    /* Count per level with exact-size allocation; weights in row order */
    ST_int *counts = (ST_int *)calloc((size_t)(level > 0 ? level : 1), sizeof(ST_int));
    if (!counts)
        return -1;
    if (weights && wcounts_out)
    {
        ST_double *wcounts = (ST_double *)calloc((size_t)(level > 0 ? level : 1), sizeof(ST_double));
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

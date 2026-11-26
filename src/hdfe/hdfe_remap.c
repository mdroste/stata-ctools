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
static int remap_dense_table(ST_int *levels_out, ST_int N, int64_t vmin, int64_t range,
                             ST_int *num_levels, int B);

/* The first pass of remap_dense_parallel: the low 32 bits of every code in
 * code[], the smallest code and the size of the range. Returns 1 when the
 * values are not integers in a manageable range. */
static int remap_dense_codes(const double *values, ST_int N, uint32_t *code, int B,
                             int64_t *vmin_out, int64_t *range_out)
{
    int all_integer = 1, b;
    int64_t vmin = INT64_MAX, vmax = INT64_MIN;
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
    *vmin_out = vmin;
    *range_out = range;
    return 0;
}

static int remap_dense_parallel(const double *values, ST_int N, ST_int *levels_out,
                                ST_int *num_levels, int B)
{
    int64_t vmin = 0, range = 0;
    if (remap_dense_codes(values, N, (uint32_t *)levels_out, B, &vmin, &range))
        return 1;
    return remap_dense_table(levels_out, N, vmin, range, num_levels, B);
}

/* The rest of remap_dense_parallel: levels_out holds the codes' low 32 bits */
static int remap_dense_table(ST_int *levels_out, ST_int N, int64_t vmin, int64_t range,
                             ST_int *num_levels, int B)
{
    int b;
    uint32_t *code = (uint32_t *)levels_out;
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

/* ========================================================================
 * Levels in order of first appearance, in parallel
 * ======================================================================== */

/* Per-block bitsets of remap_first_order may take this much, or two bytes a
 * row when that is more */
#define REMAP_ORDER_MIN_BUDGET ((size_t)64 << 20)
/* Codes from which unordered rows hand their numbers over in ranges */
#ifndef REMAP_SCATTER_CODES /* tests set a small one */
#define REMAP_SCATTER_CODES 65536
#endif

/* A code and the level it takes (sign bit: the code is in one row) */
typedef struct
{
    uint32_t code;
    ST_int id;
} remap_code_id;

/* a when take is 1, b when it is 0, without a branch */
static inline void *remap_pick(void *a, void *b, uintptr_t take)
{
    return (void *)((uintptr_t)b ^ (((uintptr_t)a ^ (uintptr_t)b) & ((uintptr_t)0 - take)));
}

static int remap_order_fits(size_t range, ST_int N, int B)
{
    size_t budget = 2 * (size_t)N;
    if (budget < REMAP_ORDER_MIN_BUDGET)
        budget = REMAP_ORDER_MIN_BUDGET;
    return (range + 63) / 64 <= budget / 16 / (size_t)B;
}

/*
 * levels_out holds a code in lo..lo+range-1 for every row (modulo 2^32).
 * Renumber the codes present 1..L in order of first appearance, in B row
 * blocks and without atomic operations:
 *   1. each block marks its codes in a bitset of its own, marks in a second
 *      one those it meets again, and flags the rows where it meets a code
 *      for the first time
 *   2. word by word, a code stays marked only in the first block that has
 *      it; a code is in one row only when one block has it and that block
 *      met it once. Counting each block's marks gives the number its first
 *      new level takes
 *   3. each block numbers its marked codes down its flagged rows, which is
 *      the order of their first rows
 *   4. every row reads its level
 * Also returned: the first row of every level, and the rows alone in their
 * level (ascending). Returns 0, or -1 when an allocation fails; levels_out
 * then still holds the codes.
 */
static int remap_first_order(ST_int *levels_out, ST_int N, uint32_t lo, size_t range, int B,
                             ST_int *num_levels, HDFE_RemapOrder *order)
{
    const uint32_t *code = (const uint32_t *)levels_out;
    const size_t words = (range + 63) / 64;
    /* A block's marks: marks[b * stride + 2 * word] = its codes, the next
     * word those it met more than once (after step 2, block 0's second
     * words hold the codes in one row only). At least two cache lines a
     * block, so that the blocks of a small range never write the same one */
    const size_t stride = words * 2 < 16 ? 16 : words * 2;
    /* The code words in B ranges: range p holds the words w with
     * (w * part_mult) >> 32 == p */
    const uint64_t part_mult = (((uint64_t)B) << 32) / (uint64_t)words;
    const size_t place_stride = (size_t)B + 32;
    int rc = -1, b;
    long long w;
    ST_int *first_row = NULL, *single_rows = NULL, *bucket = NULL, *places = NULL;
    remap_code_id *pairs = NULL;
    uint64_t *marks = (uint64_t *)calloc((size_t)B * stride, sizeof(uint64_t));
    /* fresh: a bit per row, set where its block meets a code for the first
     * time; block b's rows start at word (first row / 64) + b */
    uint64_t *fresh = (uint64_t *)malloc(((size_t)N / 64 + (size_t)B + 1) * sizeof(uint64_t));
    ST_int *id = (ST_int *)malloc(range * sizeof(ST_int));
    ST_int *base = (ST_int *)calloc((size_t)B + 1, sizeof(ST_int));
    ST_int *sbase = (ST_int *)calloc((size_t)B + 1, sizeof(ST_int));
    if (!marks || !fresh || !id || !base || !sbase)
        goto done;

    /* 1. Codes of each block, and the rows where it first meets one */
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        uint64_t *m = marks + (size_t)b * stride;
        ST_int i0, i1;
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        uint64_t *fr = fresh + ((size_t)i0 >> 6) + (size_t)b;
        const ST_int rows = i1 - i0;
        uint64_t acc = 0;
        if (words == 1)
        {
            /* one word: kept in registers */
            uint64_t mine = 0, again = 0;
            for (ST_int r = 0; r < rows; r++)
            {
                const uint64_t bit = (uint64_t)1 << ((code[i0 + r] - lo) & 63), had = mine & bit;
                again |= had;
                mine |= bit;
                acc |= (uint64_t)!had << (r & 63);
                if ((r & 63) == 63)
                {
                    fr[r >> 6] = acc;
                    acc = 0;
                }
            }
            m[0] = mine;
            m[1] = again;
        }
        else
        {
            for (ST_int r = 0; r < rows; r++)
            {
                const uint32_t c = code[i0 + r] - lo;
                uint64_t *slot = &m[2 * (size_t)(c >> 6)];
                const uint64_t bit = (uint64_t)1 << (c & 63), had = slot[0] & bit;
                slot[1] |= had;
                slot[0] |= bit;
                acc |= (uint64_t)!had << (r & 63);
                if ((r & 63) == 63)
                {
                    fr[r >> 6] = acc;
                    acc = 0;
                }
            }
        }
        if (rows & 63)
            fr[rows >> 6] = acc;
    }

    /* 2. First block of every code, and the codes in one row */
    long long shared_codes = 0; /* codes that more than one block has */
#pragma omp parallel for schedule(static) reduction(+ : shared_codes)
    for (w = 0; w < (long long)words; w++)
    {
        uint64_t seen = 0, several = 0, shared = 0;
        for (int t = 0; t < B; t++)
        {
            uint64_t *slot = &marks[(size_t)t * stride + (size_t)w * 2];
            shared |= seen & slot[0];
            several |= (seen & slot[0]) | slot[1];
            const uint64_t mine = slot[0] & ~seen;
            seen |= slot[0];
            slot[0] = mine;
        }
        shared_codes += __builtin_popcountll(shared);
        /* each block's second word: its new codes that are in one row */
        const uint64_t once = seen & ~several;
        for (int t = 0; t < B; t++)
        {
            uint64_t *slot = &marks[(size_t)t * stride + (size_t)w * 2];
            slot[1] = slot[0] & once;
        }
    }
    /* New codes of each block in each range of code words, and its codes in
     * one row */
    bucket = (ST_int *)calloc((size_t)B * (size_t)B + 1, sizeof(ST_int));
    if (!bucket)
        goto done;
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        const uint64_t *m = marks + (size_t)b * stride;
        ST_int *count = bucket + (size_t)b * (size_t)B;
        ST_int s = 0;
        for (size_t k = 0; k < words; k++)
        {
            count[(k * part_mult) >> 32] += __builtin_popcountll(m[2 * k]);
            s += __builtin_popcountll(m[2 * k + 1]);
        }
        sbase[b + 1] = s;
    }
    /* Codes that are not new to the block with the most new codes of their
     * range: when they are many (rows in no order of the codes), threads
     * numbering their blocks would keep writing the same parts of id[] */
    ST_int mixed = 0;
    for (int p = 0; p < B; p++)
    {
        ST_int all = 0, most = 0;
        for (b = 0; b < B; b++)
        {
            const ST_int n = bucket[(size_t)b * (size_t)B + (size_t)p];
            all += n;
            if (n > most)
                most = n;
        }
        mixed += all - most;
    }
    {
        /* bucket[b * B + p] becomes the place of the first such code */
        ST_int total = 0;
        for (b = 0; b < B; b++)
        {
            base[b] = total;
            for (int p = 0; p < B; p++)
            {
                const ST_int n = bucket[(size_t)b * (size_t)B + (size_t)p];
                bucket[(size_t)b * (size_t)B + (size_t)p] = total;
                total += n;
            }
            sbase[b + 1] += sbase[b];
        }
        base[B] = total;
        bucket[(size_t)B * (size_t)B] = total;
    }
    const ST_int L = base[B], singles = sbase[B];
    /* Then a block hands the numbers of its codes to the thread of their
     * range, and steps down its rows without branches (most of the rows
     * where it first meets a code are not that code's first) */
    const int scatter = L >= REMAP_SCATTER_CODES && mixed > L / 8;
    first_row = (ST_int *)malloc((size_t)(L > 0 ? L : 1) * sizeof(ST_int));
    single_rows = (ST_int *)malloc((size_t)(singles > 0 ? singles : 1) * sizeof(ST_int));
    if (scatter)
    {
        pairs = (remap_code_id *)malloc((size_t)L * sizeof(remap_code_id));
        /* each block's places apart from the next block's: they are
         * written at every row */
        places = (ST_int *)malloc((size_t)B * place_stride * sizeof(ST_int));
        for (b = 0; places && b < B; b++)
            memcpy(places + (size_t)b * place_stride, bucket + (size_t)b * (size_t)B,
                   (size_t)B * sizeof(ST_int));
    }
    if (!first_row || !single_rows || (scatter && (!pairs || !places)))
        goto done;

    /* 3. Number each block's new codes in the order of their first rows;
     * the sign bit flags a code in one row */
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        const uint64_t *m = marks + (size_t)b * stride;
        ST_int i0, i1, next = base[b];
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        const uint64_t *fr = fresh + ((size_t)i0 >> 6) + (size_t)b;
        const size_t fresh_words = ((size_t)(i1 - i0) + 63) / 64;
        if (scatter)
        {
            /* a row whose code another block met first writes to the spare
             * entries instead */
            ST_int *place = places + (size_t)b * place_stride, spare_row = 0;
            remap_code_id spare_pair;
            for (size_t k = 0; k < fresh_words; k++)
            {
                uint64_t rows = fr[k];
                while (rows)
                {
                    const ST_int i = i0 + (ST_int)(k * 64) + __builtin_ctzll(rows);
                    rows &= rows - 1;
                    const uint32_t c = code[i] - lo;
                    const uint64_t *slot = &m[2 * (size_t)(c >> 6)];
                    const uintptr_t first = (uintptr_t)((slot[0] >> (c & 63)) & 1);
                    const uint32_t once = (uint32_t)((slot[1] >> (c & 63)) & 1);
                    ST_int *at = &place[((uint64_t)(c >> 6) * part_mult) >> 32];
                    ST_int *row_to = (ST_int *)remap_pick(&first_row[next], &spare_row, first);
                    remap_code_id *pair = (remap_code_id *)remap_pick(&pairs[*at], &spare_pair, first);
                    *row_to = i;
                    next += (ST_int)first;
                    *at += (ST_int)first;
                    pair->code = c;
                    pair->id = (ST_int)((uint32_t)next | (once << 31));
                }
            }
            (void)spare_row;
            (void)spare_pair;
            continue;
        }
        for (size_t k = 0; k < fresh_words; k++)
        {
            uint64_t rows = fr[k];
            while (rows)
            {
                const ST_int i = i0 + (ST_int)(k * 64) + __builtin_ctzll(rows);
                rows &= rows - 1;
                const uint32_t c = code[i] - lo;
                const uint64_t *slot = &m[2 * (size_t)(c >> 6)];
                const uint64_t bit = (uint64_t)1 << (c & 63);
                if (!(slot[0] & bit))
                    continue;
                first_row[next] = i;
                ++next;
                id[c] = (slot[1] & bit) ? (ST_int)((uint32_t)next | 0x80000000u) : next;
            }
        }
    }
    if (scatter)
    {
        int p;
#pragma omp parallel for schedule(static)
        for (p = 0; p < B; p++)
        {
            for (int t = 0; t < B; t++)
            {
                const size_t at = (size_t)t * (size_t)B + (size_t)p;
                /* the next bucket starts where this one ends */
                const ST_int stop = bucket[at + 1];
                for (ST_int k = bucket[at]; k < stop; k++)
                    id[pairs[k].code] = pairs[k].id;
            }
        }
    }

    /* 4. Levels of the rows; a row alone in its level is in the block that
     * numbered the level */
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        ST_int i0, i1, s = sbase[b];
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
        {
            const ST_int v = id[(uint32_t)levels_out[i] - lo];
            levels_out[i] = v & 0x7fffffff;
            if (v < 0)
                single_rows[s++] = i;
        }
    }
    *num_levels = L;
    order->first_row = first_row;
    order->single_rows = single_rows;
    order->num_single = singles;
    order->local = shared_codes <= L / 16;
    first_row = single_rows = NULL;
    rc = 0;

done:
    free(first_row);
    free(single_rows);
    free(bucket);
    free(places);
    free(pairs);
    free(marks);
    free(fresh);
    free(id);
    free(base);
    free(sbase);
    return rc;
}

int hdfe_remap_levels_ordered(const double *values, ST_int N, ST_int *levels_out,
                              ST_int *num_levels, HDFE_RemapOrder *order)
{
    if (num_levels)
        *num_levels = 0;
    if (order)
    {
        order->first_row = NULL;
        order->single_rows = NULL;
        order->num_single = -1;
        order->local = 0;
    }
    if (N <= 0 || !values || !levels_out || !num_levels || !order)
        return -1;

    const int B = ctools_hdfe_row_blocks(N);
    ST_int L = 0;
    if (B == 1)
    {
        /* Serially the table and the hash both number in row order */
        if (remap_levels_impl(values, N, levels_out, &L, 0) != 0)
            return -1;
        *num_levels = L;
        return 0;
    }

    int64_t vmin = 0, range = 0;
    if (remap_dense_codes(values, N, (uint32_t *)levels_out, B, &vmin, &range) == 0)
    {
        if (remap_order_fits((size_t)range, N, B) &&
            remap_first_order(levels_out, N, (uint32_t)vmin, (size_t)range, B, &L, order) == 0)
        {
            *num_levels = L;
            return 0;
        }
        /* A range too wide for the bitsets: ascending levels first */
        if (remap_dense_table(levels_out, N, vmin, range, &L, B) != 0)
            return -1;
    }
    else
    {
        int rc = remap_sparse_parallel(values, N, levels_out, &L, B);
        if (rc != 0)
        {
            if (remap_hash_impl(values, N, levels_out, &L) != 0)
                return -1;
            *num_levels = L;
            return 0;
        }
    }
    /* Levels 1..L in another order: renumber them */
    *num_levels = L;
    if (L > 0 && remap_order_fits((size_t)L, N, B) &&
        remap_first_order(levels_out, N, 1, (size_t)L, B, &L, order) == 0)
        return 0;
    return 1;
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

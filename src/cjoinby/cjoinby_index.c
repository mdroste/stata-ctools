/*
    cjoinby_index.c
    Using-side key index for cjoinby

    The using observations are put in stable key order, so the using rows of
    a key group keep their order in the file; the runs of equal keys become
    groups, and the distinct key values are copied into compact per-key
    arrays. Master rows find their group through one of two lookups:

    - direct: when every key is an integer (or a missing value) and the
      product of the key ranges is small, a key tuple is a mixed-radix
      position in a table of groups, one read per master row. The index is
      then built from those positions alone: one radix sort of 8-byte
      (position, row) words gives the key order, the groups and the table;
    - hashed: otherwise the rows are ordered by csort's stable engines
      (counting or radix passes; data already in key order are only
      checked), and the groups go into an open-addressing table, built with
      compare-and-swap inserts in parallel (the keys are distinct, so an
      insert only needs a free slot) and probed in parallel.

    Keys compare as joinby's merge compares them: numbers as doubles (each
    extended missing value is its own value, and -0 equals 0), strings by
    their bytes.
*/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ctools_config.h"
#include "ctools_order.h"
#include "ctools_sort.h"
#include "ctools_release.h"
#include "cjoinby_internal.h"

#define CJ_HASH_EMPTY UINT64_MAX
/* Rows a worker should have before another worker is worth waking */
#define CJ_ROWS_PER_THREAD 16384
/* A direct table may hold this many entries however few rows there are
 * (16 MB), and beyond that two entries per using row or one per master row */
#define CJ_DIRECT_MIN_ENTRIES ((uint64_t)1 << 22)
#define CJ_DIRECT_MAX_ENTRIES ((uint64_t)1 << 30)
/* Stata's missing values as doubles: . and .a-.z are 27 consecutive codes */
#define CJ_MISSING_BITS UINT64_C(0x7FE0000000000000)
#define CJ_MISSING_CODES 27
/* Integers a double holds exactly */
#define CJ_EXACT_INT 9007199254740992.0

static int cj_threads_for(size_t n)
{
    int threads = ctools_get_max_threads();
    if (threads < 1)
        threads = 1;
    if (n / CJ_ROWS_PER_THREAD < (size_t)threads)
        threads = (int)(n / CJ_ROWS_PER_THREAD);
    return threads < 1 ? 1 : threads;
}

/* 0-26 for ., .a-.z; -1 for anything else at or above the missing range */
static inline int cj_missing_code(double x)
{
    uint64_t bits;
    memcpy(&bits, &x, sizeof(bits));
    uint64_t delta = bits - CJ_MISSING_BITS;
    uint64_t code = delta >> 40;
    if (code >= CJ_MISSING_CODES || (delta & ((UINT64_C(1) << 40) - 1)))
        return -1;
    return (int)code;
}

static inline uint64_t cj_mix(uint64_t h)
{
    h ^= h >> 33;
    h *= UINT64_C(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h *= UINT64_C(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    return h;
}

static inline uint64_t cj_hash_double(uint64_t h, double x)
{
    uint64_t bits;
    if (x == 0.0)
        x = 0.0; /* -0 and 0 are one key */
    memcpy(&bits, &x, sizeof(bits));
    return cj_mix(h ^ bits) + UINT64_C(0x9e3779b97f4a7c15);
}

static inline uint64_t cj_hash_string(uint64_t h, const char *s)
{
    uint64_t f = UINT64_C(0xcbf29ce484222325);
    for (const unsigned char *p = (const unsigned char *)s; *p; p++)
        f = (f ^ *p) * UINT64_C(0x100000001b3);
    return cj_mix(h ^ f) + UINT64_C(0x9e3779b97f4a7c15);
}

void cj_index_init(cj_index *ix)
{
    if (!ix)
        return;
    memset(ix, 0, sizeof(*ix));
    stata_data_init(&ix->keys);
}

static void cj_index_drop(cj_index *ix, int later)
{
    if (!ix)
        return;
    if (ix->keytab)
    {
        for (int k = 0; k < ix->nkeys; k++)
        {
            if (later)
            {
                ctools_release_later(ix->keytab[k].num);
                ctools_release_later(ix->keytab[k].str);
            }
            else
            {
                free(ix->keytab[k].num);
                free(ix->keytab[k].str);
            }
        }
        free(ix->keytab);
    }
    if (later)
    {
        ctools_release_later(ix->direct.table);
        ctools_release_later(ix->hash.slots);
        ctools_release_later(ix->gstart);
    }
    else
    {
        free(ix->direct.table);
        free(ix->hash.slots);
        free(ix->gstart);
    }
    free(ix->direct.base);
    free(ix->direct.top);
    free(ix->direct.finite);
    free(ix->direct.mult);
    free(ix->direct.has_missing);
    free(ix->keycol);
    stata_data_free(&ix->keys);
    cj_index_init(ix);
}

void cj_index_free(cj_index *ix)
{
    cj_index_drop(ix, 0);
}

void cj_index_release(cj_index *ix)
{
    cj_index_drop(ix, 1);
}

/* The longest value of a string column (at most STATA's str# maximum) */
static size_t cj_longest(const stata_variable *var)
{
    const size_t n = var->nobs;
    char *const *s = var->data.str;
    size_t longest = 0;
    int threads = cj_threads_for(n);
    (void)threads;
#pragma omp parallel for schedule(static) reduction(max : longest)                                 \
    num_threads(threads) if (threads > 1)
    for (size_t i = 0; i < n; i++)
    {
        size_t len = s[i] ? strnlen(s[i], 2045) : 0;
        if (len > longest)
            longest = len;
    }
    return longest;
}

/* Copy each group's key values into compact arrays, in group order */
static stata_retcode cj_index_keytab(cj_index *ix)
{
    const size_t ng = ix->ngroups;
    const perm_idx_t *order = ix->keys.sort_order;
    const uint32_t *gstart = ix->gstart;
    int threads = cj_threads_for(ng);
    (void)threads;

    ix->keytab = (cj_keycol *)calloc((size_t)ix->nkeys, sizeof(cj_keycol));
    if (!ix->keytab)
        return STATA_ERR_MEMORY;
    for (int k = 0; k < ix->nkeys; k++)
    {
        const stata_variable *var = &ix->keys.vars[k];
        cj_keycol *kc = &ix->keytab[k];
        if (var->type == STATA_TYPE_DOUBLE)
        {
            const double *x = var->data.dbl;
            double *num = (double *)ctools_safe_malloc2(ng ? ng : 1, sizeof(double));
            if (!num)
                return STATA_ERR_MEMORY;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
            for (size_t g = 0; g < ng; g++)
                num[g] = x[order[gstart[g]]];
            kc->num = num;
        }
        else
        {
            char *const *s = var->data.str;
            const size_t slot = cj_longest(var) + 1;
            char *str = (char *)ctools_safe_calloc2(ng ? ng : 1, slot);
            if (!str)
                return STATA_ERR_MEMORY;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
            for (size_t g = 0; g < ng; g++)
            {
                const char *value = s[order[gstart[g]]];
                if (value)
                    memcpy(str + g * slot, value, strnlen(value, slot - 1));
            }
            kc->is_str = 1;
            kc->str = str;
            kc->slot = slot;
        }
    }
    return STATA_OK;
}

/* Decide whether the keys qualify for the direct lookup, from every using
 * row: all numeric, every value an integer or a missing value, and a product
 * of the key ranges within the table limit. With *eligible set, d->base, top,
 * finite, mult, has_missing and size describe the table. */
static stata_retcode cj_direct_scan(cj_index *ix, size_t master_hint, int *eligible)
{
    const size_t n = ix->nobs;
    const int nk = ix->nkeys;
    cj_direct *d = &ix->direct;

    *eligible = 0;
    for (int k = 0; k < nk; k++)
        if (ix->keys.vars[k].type != STATA_TYPE_DOUBLE)
            return STATA_OK;

    d->base = (double *)calloc((size_t)nk, sizeof(double));
    d->top = (double *)calloc((size_t)nk, sizeof(double));
    d->finite = (uint64_t *)calloc((size_t)nk, sizeof(uint64_t));
    d->mult = (uint64_t *)calloc((size_t)nk, sizeof(uint64_t));
    d->has_missing = (int *)calloc((size_t)nk, sizeof(int));
    if (!d->base || !d->top || !d->finite || !d->mult || !d->has_missing)
        return STATA_ERR_MEMORY;

    uint64_t limit = (uint64_t)n * 2;
    if (limit < (uint64_t)master_hint)
        limit = (uint64_t)master_hint;
    if (limit < CJ_DIRECT_MIN_ENTRIES)
        limit = CJ_DIRECT_MIN_ENTRIES;
    if (limit > CJ_DIRECT_MAX_ENTRIES)
        limit = CJ_DIRECT_MAX_ENTRIES;

    int threads = cj_threads_for(n);
    (void)threads;
    uint64_t size = 1;
    for (int k = nk - 1; k >= 0; k--)
    {
        const double *x = ix->keys.vars[k].data.dbl;
        double lo = 0, hi = 0;
        int any_finite = 0, any_missing = 0, bad = 0;
#pragma omp parallel num_threads(threads) if (threads > 1)
        {
            double tlo = 0, thi = 0;
            int tfinite = 0, tmissing = 0, tbad = 0;
#pragma omp for schedule(static) nowait
            for (size_t i = 0; i < n; i++)
            {
                double v = x[i];
                if (v >= SV_missval)
                {
                    if (cj_missing_code(v) < 0)
                        tbad = 1;
                    tmissing = 1;
                }
                else if (!(v >= -CJ_EXACT_INT && v <= CJ_EXACT_INT) || v != (double)(int64_t)v)
                {
                    tbad = 1;
                }
                else if (!tfinite)
                {
                    tlo = thi = v;
                    tfinite = 1;
                }
                else
                {
                    if (v < tlo)
                        tlo = v;
                    if (v > thi)
                        thi = v;
                }
            }
#pragma omp critical(cjoinby_direct_range)
            {
                bad |= tbad;
                any_missing |= tmissing;
                if (tfinite)
                {
                    if (!any_finite || tlo < lo)
                        lo = tlo;
                    if (!any_finite || thi > hi)
                        hi = thi;
                    any_finite = 1;
                }
            }
        }
        if (bad)
            return STATA_OK;
        double span = any_finite ? hi - lo + 1.0 : 0.0;
        if (span > (double)CJ_DIRECT_MAX_ENTRIES)
            return STATA_OK;
        uint64_t codes = (uint64_t)span + (any_missing ? CJ_MISSING_CODES : 0);
        if (codes == 0)
            codes = 1;
        d->base[k] = lo;
        d->top[k] = any_finite ? hi : lo - 1.0;
        d->finite[k] = (uint64_t)span;
        d->has_missing[k] = any_missing;
        d->mult[k] = size;
        if (codes > limit / size)
            return STATA_OK;
        size *= codes;
    }
    d->size = (size_t)size;
    *eligible = 1;
    return STATA_OK;
}

/* Stable sort of packed (table position << 32 | row) words by position: LSD
 * radix passes over the position bytes that can differ. The rows are packed
 * in ascending order, so they stay ascending within a position. */
static stata_retcode cj_sort_packed(uint64_t *packed, size_t n, uint64_t size)
{
    int threads = cj_threads_for(n);
    uint64_t *tmp = (uint64_t *)ctools_safe_malloc2(n, sizeof(uint64_t));
    size_t *hist = (size_t *)ctools_safe_malloc2((size_t)threads * 256, sizeof(size_t));
    if (!tmp || !hist)
    {
        free(tmp);
        free(hist);
        return STATA_ERR_MEMORY;
    }
    uint64_t *src = packed, *dst = tmp;
    for (unsigned shift = 32; shift < 64 && ((size - 1) >> (shift - 32)) != 0; shift += 8)
    {
#pragma omp parallel num_threads(threads) if (threads > 1)
        {
            size_t team = 1, tid = 0;
#ifdef _OPENMP
            team = (size_t)omp_get_num_threads();
            tid = (size_t)omp_get_thread_num();
#endif
            const size_t begin = n * tid / team, end = n * (tid + 1) / team;
            size_t *h = hist + tid * 256;
            memset(h, 0, 256 * sizeof(size_t));
            for (size_t i = begin; i < end; i++)
                h[(src[i] >> shift) & 255]++;
#pragma omp barrier
#pragma omp single
            {
                size_t offset = 0;
                for (size_t b = 0; b < 256; b++)
                {
                    for (size_t t = 0; t < team; t++)
                    {
                        size_t *slot = hist + t * 256 + b, count = *slot;
                        *slot = offset;
                        offset += count;
                    }
                }
            }
            for (size_t i = begin; i < end; i++)
                dst[h[(src[i] >> shift) & 255]++] = src[i];
        }
        uint64_t *swap = src;
        src = dst;
        dst = swap;
    }
    if (src != packed)
        memcpy(packed, src, n * sizeof(uint64_t));
    free(tmp);
    free(hist);
    return STATA_OK;
}

/* Index the using rows by table position: the key order, the groups and the
 * table come from one sort of 8-byte words, with no key comparison. */
static stata_retcode cj_index_direct(cj_index *ix)
{
    const size_t n = ix->nobs;
    const int nk = ix->nkeys;
    cj_direct *d = &ix->direct;
    const stata_variable *keys = ix->keys.vars;
    int threads = cj_threads_for(n);
    (void)threads;

    uint64_t *packed = (uint64_t *)ctools_safe_malloc2(n, sizeof(uint64_t));
    perm_idx_t *order = (perm_idx_t *)ctools_safe_cacheline_alloc2(n, sizeof(perm_idx_t));
    if (!packed || !order)
    {
        free(packed);
        ctools_aligned_free(order);
        return STATA_ERR_MEMORY;
    }
    ix->keys.sort_order = order;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (size_t i = 0; i < n; i++)
    {
        uint64_t pos = 0;
        for (int k = 0; k < nk; k++)
        {
            double v = keys[k].data.dbl[i];
            uint64_t code = v >= SV_missval ? d->finite[k] + (uint64_t)cj_missing_code(v)
                                            : (uint64_t)(v - d->base[k]);
            pos += code * d->mult[k];
        }
        packed[i] = (pos << 32) | (uint64_t)i;
    }
    int sorted = 1;
#pragma omp parallel for schedule(static) reduction(&& : sorted)                                   \
    num_threads(threads) if (threads > 1)
    for (size_t i = 1; i < n; i++)
        sorted = sorted && (packed[i - 1] >> 32) <= (packed[i] >> 32);
    if (!sorted)
    {
        stata_retcode rc = cj_sort_packed(packed, n, (uint64_t)d->size);
        if (rc)
        {
            free(packed);
            return rc;
        }
    }

    size_t ng = 1;
#pragma omp parallel for schedule(static) reduction(+ : ng) num_threads(threads) if (threads > 1)
    for (size_t p = 1; p < n; p++)
        ng += (packed[p - 1] >> 32) != (packed[p] >> 32);
    uint32_t *gstart = (uint32_t *)ctools_safe_malloc2(ng + 1, sizeof(uint32_t));
    uint32_t *table = (uint32_t *)ctools_safe_malloc2(d->size, sizeof(uint32_t));
    if (!gstart || !table)
    {
        free(packed);
        free(gstart);
        free(table);
        return STATA_ERR_MEMORY;
    }
    memset(table, 0xFF, d->size * sizeof(uint32_t));
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (size_t p = 0; p < n; p++)
        order[p] = (perm_idx_t)(packed[p] & UINT64_C(0xFFFFFFFF));
    size_t g = 0;
    uint64_t previous = 0;
    for (size_t p = 0; p < n; p++)
    {
        const uint64_t pos = packed[p] >> 32;
        if (p > 0 && pos == previous)
            continue;
        gstart[g] = (uint32_t)p;
        table[pos] = (uint32_t)g++;
        previous = pos;
    }
    gstart[ng] = (uint32_t)n;
    free(packed);
    free(ix->gstart);
    ix->gstart = gstart;
    ix->ngroups = ng;
    d->table = table;
    ix->use_direct = 1;
    return STATA_OK;
}

static inline uint64_t cj_hash_group(const cj_index *ix, size_t g)
{
    uint64_t h = 0;
    for (int k = 0; k < ix->nkeys; k++)
    {
        const cj_keycol *kc = &ix->keytab[k];
        h = kc->is_str ? cj_hash_string(h, kc->str + g * kc->slot) : cj_hash_double(h, kc->num[g]);
    }
    return h;
}

static stata_retcode cj_index_hash(cj_index *ix)
{
    const size_t ng = ix->ngroups;
    size_t capacity = 16;
    while (capacity < ng + ng / 2 + 1)
    {
        if (capacity > SIZE_MAX / 2)
            return STATA_ERR_MEMORY;
        capacity *= 2;
    }
    uint64_t *slots = (uint64_t *)ctools_safe_malloc2(capacity, sizeof(uint64_t));
    if (!slots)
        return STATA_ERR_MEMORY;
    memset(slots, 0xFF, capacity * sizeof(uint64_t));
    const size_t mask = capacity - 1;
    int threads = cj_threads_for(ng);
    (void)threads;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (size_t g = 0; g < ng; g++)
    {
        uint64_t h = cj_hash_group(ix, g);
        uint64_t entry = (h & UINT64_C(0xFFFFFFFF00000000)) | (uint64_t)g;
        size_t i = (size_t)h & mask;
        for (;;)
        {
            uint64_t expected = CJ_HASH_EMPTY;
            if (__atomic_compare_exchange_n(&slots[i], &expected, entry, 0, __ATOMIC_RELAXED,
                                            __ATOMIC_RELAXED))
                break;
            i = (i + 1) & mask;
        }
    }
    ix->hash.slots = slots;
    ix->hash.mask = mask;
    return STATA_OK;
}

/* Put ix->keys.sort_order in stable key order. Sorted using data (panels
 * usually are) cost one pass of comparisons. Otherwise the dispatcher picks
 * a stable engine by the keys: counting or LSD radix passes for numbers and
 * short strings, MSD radix for long strings. It chooses by str_maxlen, so
 * string keys loaded without a declared width get their longest value. */
static stata_retcode cj_index_order(cj_index *ix)
{
    const size_t n = ix->nobs;
    const int nk = ix->nkeys;
    int threads = cj_threads_for(n);
    (void)threads;

    int sorted = 1;
#pragma omp parallel for schedule(static) reduction(&& : sorted)                                   \
    num_threads(threads) if (threads > 1)
    for (size_t i = 1; i < n; i++)
        sorted = sorted && ctools_compare_rows(&ix->keys, i - 1, ix->keycol, &ix->keys, i,
                                               ix->keycol, (size_t)nk) <= 0;
    if (sorted)
        return STATA_OK;

    int *sortvars = (int *)malloc((size_t)nk * sizeof(int));
    if (!sortvars)
        return STATA_ERR_MEMORY;
    for (int k = 0; k < nk; k++)
    {
        stata_variable *var = &ix->keys.vars[k];
        sortvars[k] = k + 1;
        if (var->type == STATA_TYPE_STRING && (var->str_maxlen == 0 || var->str_maxlen >= 2045))
            var->str_maxlen = cj_longest(var);
    }
    stata_retcode rc = ctools_sort_dispatch(&ix->keys, sortvars, (size_t)nk, SORT_ALG_AUTO);
    free(sortvars);
    return rc;
}

stata_retcode cj_index_build(cj_index *ix, const cj_table *utab, int nkeys, size_t master_hint)
{
    if (!ix || !utab || nkeys < 1 || (size_t)nkeys > utab->nvars)
        return STATA_ERR_INVALID_INPUT;
    cj_index_init(ix);
    const size_t n = utab->nobs;
    stata_retcode rc = STATA_OK;

    ix->nkeys = nkeys;
    ix->nobs = n;
    ix->keycol = (int *)malloc((size_t)nkeys * sizeof(int));
    ix->keys.vars =
        (stata_variable *)ctools_cacheline_alloc((size_t)nkeys * sizeof(stata_variable));
    if (!ix->keycol || !ix->keys.vars)
        return STATA_ERR_MEMORY;
    memset(ix->keys.vars, 0, (size_t)nkeys * sizeof(stata_variable));
    ix->keys.nvars = (size_t)nkeys;
    ix->keys.nobs = n;
    for (int k = 0; k < nkeys && !rc; k++)
    {
        ix->keycol[k] = k;
        rc = cj_table_column(utab, (size_t)k, &ix->keys.vars[k]);
    }
    if (rc)
        return rc;

    ix->gstart = (uint32_t *)malloc(sizeof(uint32_t));
    if (!ix->gstart)
        return STATA_ERR_MEMORY;
    ix->gstart[0] = 0;
    if (n == 0)
        return cj_index_keytab(ix);

    /* Integer keys of a small combined range: sorted by table position */
    int direct = 0;
    rc = cj_direct_scan(ix, master_hint, &direct);
    if (!rc && direct)
        rc = cj_index_direct(ix);
    if (rc)
        return rc;
    if (ix->use_direct)
        return cj_index_keytab(ix);

    /* Any other keys: stable key order, then the runs of equal keys */
    perm_idx_t *order = (perm_idx_t *)ctools_safe_cacheline_alloc2(n, sizeof(perm_idx_t));
    if (!order)
        return STATA_ERR_MEMORY;
    int threads = cj_threads_for(n);
    (void)threads;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (size_t i = 0; i < n; i++)
        order[i] = (perm_idx_t)i;
    ix->keys.sort_order = order;
    rc = cj_index_order(ix);
    if (rc)
        return rc;

    size_t *starts = (size_t *)ctools_safe_malloc2(n + 1, sizeof(size_t));
    if (!starts)
        return STATA_ERR_MEMORY;
    const size_t ng = ctools_group_starts_parallel(&ix->keys, ix->keys.sort_order, ix->keycol,
                                                   (size_t)nkeys, starts);
    uint32_t *gstart = (uint32_t *)ctools_safe_malloc2(ng + 1, sizeof(uint32_t));
    if (!gstart)
    {
        free(starts);
        return STATA_ERR_MEMORY;
    }
    for (size_t g = 0; g <= ng; g++)
        gstart[g] = (uint32_t)starts[g];
    free(starts);
    free(ix->gstart);
    ix->gstart = gstart;
    ix->ngroups = ng;

    rc = cj_index_keytab(ix);
    if (!rc)
        rc = cj_index_hash(ix);
    return rc;
}

/* 1 when master row i and its predecessor hold the same key values */
static inline int cj_same_as_previous(const stata_variable *mk, int nk, size_t i)
{
    for (int k = 0; k < nk; k++)
    {
        if (mk[k].type == STATA_TYPE_DOUBLE)
        {
            if (mk[k].data.dbl[i] != mk[k].data.dbl[i - 1])
                return 0;
        }
        else if (strcmp(mk[k].data.str[i], mk[k].data.str[i - 1]) != 0)
        {
            return 0;
        }
    }
    return 1;
}

static inline uint32_t cj_lookup_direct(const cj_index *ix, const stata_variable *mk, size_t i)
{
    const cj_direct *d = &ix->direct;
    uint64_t pos = 0;
    for (int k = 0; k < ix->nkeys; k++)
    {
        double v = mk[k].data.dbl[i];
        uint64_t code;
        if (v >= SV_missval)
        {
            int m = d->has_missing[k] ? cj_missing_code(v) : -1;
            if (m < 0)
                return CJ_NONE;
            code = d->finite[k] + (uint64_t)m;
        }
        else
        {
            if (!(v >= d->base[k] && v <= d->top[k]))
                return CJ_NONE;
            /* Test v itself: v - base can round a value one step from an
             * integer onto that integer (1e-15 - (-127) == 127 exactly).
             * base and top are integers within CJ_EXACT_INT, so the cast is
             * defined and, for an integer v, the difference is exact. */
            if (v != (double)(int64_t)v)
                return CJ_NONE;
            code = (uint64_t)(v - d->base[k]);
        }
        pos += code * d->mult[k];
    }
    return d->table[pos];
}

static inline uint32_t cj_lookup_hash(const cj_index *ix, const stata_variable *mk, size_t i)
{
    const int nk = ix->nkeys;
    uint64_t h = 0;
    for (int k = 0; k < nk; k++)
        h = mk[k].type == STATA_TYPE_STRING ? cj_hash_string(h, mk[k].data.str[i])
                                            : cj_hash_double(h, mk[k].data.dbl[i]);
    const uint64_t tag = h & UINT64_C(0xFFFFFFFF00000000);
    const uint64_t *slots = ix->hash.slots;
    const size_t mask = ix->hash.mask;
    for (size_t s = (size_t)h & mask;; s = (s + 1) & mask)
    {
        const uint64_t entry = slots[s];
        if (entry == CJ_HASH_EMPTY)
            return CJ_NONE;
        if ((entry & UINT64_C(0xFFFFFFFF00000000)) != tag)
            continue;
        const size_t g = (size_t)(entry & UINT64_C(0xFFFFFFFF));
        int equal = 1;
        for (int k = 0; k < nk && equal; k++)
        {
            const cj_keycol *kc = &ix->keytab[k];
            equal = kc->is_str ? strcmp(kc->str + g * kc->slot, mk[k].data.str[i]) == 0
                               : kc->num[g] == mk[k].data.dbl[i];
        }
        if (equal)
            return (uint32_t)g;
    }
}

stata_retcode cj_index_probe(const cj_index *ix, const stata_data *mkeys, uint32_t *gid,
                             unsigned char *matched)
{
    if (!ix || !mkeys || (mkeys->nobs && !gid) || mkeys->nvars < (size_t)ix->nkeys)
        return STATA_ERR_INVALID_INPUT;
    const size_t n = mkeys->nobs;
    const int nk = ix->nkeys;
    const stata_variable *mk = mkeys->vars;
    for (int k = 0; k < nk; k++)
        if (mk[k].type != ix->keys.vars[k].type)
            return STATA_ERR_INVALID_INPUT;
    if (n == 0)
        return STATA_OK;
    if (ix->ngroups == 0)
    {
        memset(gid, 0xFF, n * sizeof(uint32_t));
        return STATA_OK;
    }

    const int direct = ix->use_direct;
    int threads = cj_threads_for(n);
    (void)threads;
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        size_t team = 1, tid = 0;
#ifdef _OPENMP
        team = (size_t)omp_get_num_threads();
        tid = (size_t)omp_get_thread_num();
#endif
        const size_t base = n / team, extra = n % team;
        const size_t start = tid * base + (tid < extra ? tid : extra);
        const size_t end = start + base + (tid < extra);
        for (size_t i = start; i < end; i++)
        {
            uint32_t g;
            /* Runs of one key (sorted or grouped masters) look up once */
            if (i > start && cj_same_as_previous(mk, nk, i))
            {
                gid[i] = gid[i - 1];
                continue;
            }
            g = direct ? cj_lookup_direct(ix, mk, i) : cj_lookup_hash(ix, mk, i);
            gid[i] = g;
            if (g != CJ_NONE && matched && !__atomic_load_n(&matched[g], __ATOMIC_RELAXED))
                __atomic_store_n(&matched[g], 1, __ATOMIC_RELAXED);
        }
    }
    return STATA_OK;
}

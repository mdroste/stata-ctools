/*
 * ccollapse_group.c
 *
 * Observations -> groups, numbered in the sort order of the by() values.
 *
 * Four methods, tried in the order of their cost:
 *
 *   ordered  The rows are already in key order (panel data, the result of an
 *            earlier sort). One parallel pass of adjacent comparisons finds
 *            the group boundaries; nothing is sorted or hashed.
 *   dense    Every key is a whole number (missing values allowed) and the
 *            keys span a small box: each row's cell number in the box is its
 *            mixed-radix code, and the occupied cells, taken in order, are
 *            the groups in sort order.
 *   hash     Workers number the keys of contiguous blocks of rows in private
 *            tables, the block tables are merged, and only the distinct keys
 *            are sorted. Given up when a block meets too many distinct keys.
 *   sort     A stable order of all rows by key (ctools_order_stable), for
 *            data with about as many groups as rows.
 *
 * Every method gives the same numbering: groups in ascending key order, and
 * first[g] the lowest row of group g. Rows keep their order within a group,
 * which the order-dependent statistics (sums, first, last) rely on.
 */

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ccollapse_internal.h"
#include "ctools_alloc.h"
#include "ctools_config.h"
#include "ctools_order.h"

/* Rows per block of a parallel pass, and the least rows worth a team */
#define CC_GROUP_BLOCK ((size_t)65536)
#define CC_GROUP_PARALLEL_MIN ((size_t)131072)

/* Stata's missing values are 2^1023 * (1 + code/4096) for code 0 (.) to 26 (.z) */
#define CC_MISSING_BITS UINT64_C(0x7FE0000000000000)
#define CC_MISSING_CODES 27

void cc_groups_free(cc_groups *groups)
{
    if (!groups)
        return;
    free(groups->gid);
    free(groups->first);
    free(groups->starts);
    free(groups->perm);
    memset(groups, 0, sizeof(*groups));
}

static inline int cc_team(int threads, size_t nobs)
{
    if (threads < 1 || nobs < CC_GROUP_PARALLEL_MIN)
        return 1;
    size_t blocks = nobs / CC_GROUP_BLOCK;
    return (size_t)threads > blocks ? (int)blocks : threads;
}

static inline const char *cc_str(const stata_variable *v, size_t row)
{
    const char *s = v->data.str[row];
    return s ? s : "";
}

/* Keys of rows a and b: <0, 0, >0 as ctools_compare_rows, NULL strings as "" */
static inline int cc_compare_keys(const stata_data *d, size_t nby, size_t a, size_t b)
{
    for (size_t k = 0; k < nby; k++)
    {
        const stata_variable *v = &d->vars[k];
        int c;
        if (v->type == STATA_TYPE_STRING)
            c = strcmp(cc_str(v, a), cc_str(v, b));
        else
        {
            double x = v->data.dbl[a], y = v->data.dbl[b];
            c = (x > y) - (x < y);
        }
        if (c)
            return c;
    }
    return 0;
}

static inline int cc_equal_keys(const stata_data *d, size_t nby, size_t a, size_t b)
{
    for (size_t k = 0; k < nby; k++)
    {
        const stata_variable *v = &d->vars[k];
        if (v->type == STATA_TYPE_STRING)
        {
            if (strcmp(cc_str(v, a), cc_str(v, b)) != 0)
                return 0;
        }
        else if (v->data.dbl[a] != v->data.dbl[b])
            return 0;
    }
    return 1;
}

/* ===========================================================================
   One group
   =========================================================================== */

static stata_retcode cc_group_all(size_t nobs, cc_groups *out)
{
    out->gid = calloc(nobs ? nobs : 1, sizeof(*out->gid));
    out->first = malloc(sizeof(*out->first));
    out->starts = malloc(2 * sizeof(*out->starts));
    if (!out->gid || !out->first || !out->starts)
        return STATA_ERR_MEMORY;
    out->first[0] = 0;
    out->starts[0] = 0;
    out->starts[1] = nobs;
    out->ngroups = 1;
    out->contiguous = 1;
    out->method = "one group";
    out->method_code = 1;
    return STATA_OK;
}

/* ===========================================================================
   Ordered rows
   =========================================================================== */

/*
    Compare each row of [lo, hi) with the row before it. Returns 0 at the
    first row that sorts before its predecessor; otherwise sets start[i] to 1
    where row i begins a group (i == 0 included) and returns 1.
*/
static int cc_mark_starts(const stata_data *d, size_t nby, size_t lo, size_t hi,
                          unsigned char *start)
{
    if (nby == 1 && d->vars[0].type == STATA_TYPE_DOUBLE)
    {
        const double *x = d->vars[0].data.dbl;
        for (size_t i = lo; i < hi; i++)
        {
            if (i == 0)
            {
                start[0] = 1;
                continue;
            }
            if (x[i - 1] > x[i])
                return 0;
            start[i] = x[i - 1] != x[i];
        }
        return 1;
    }
    for (size_t i = lo; i < hi; i++)
    {
        if (i == 0)
        {
            start[0] = 1;
            continue;
        }
        int c = cc_compare_keys(d, nby, i - 1, i);
        if (c > 0)
            return 0;
        start[i] = c != 0;
    }
    return 1;
}

/*
    Groups of rows that are already in key order. Returns STATA_OK with
    out->ngroups == 0 when they are not.
*/
static stata_retcode cc_group_ordered(const stata_data *d, size_t nby, int threads,
                                      cc_groups *out)
{
    size_t n = d->nobs;
    unsigned char *start = malloc(n);
    if (!start)
        return STATA_ERR_MEMORY;

    /* Unordered data usually show it within the first rows */
    size_t probe = n < 4096 ? n : 4096;
    int ordered = cc_mark_starts(d, nby, 0, probe, start);
    if (ordered && probe < n)
    {
        int team = cc_team(threads, n);
        size_t nblocks = (n - probe + CC_GROUP_BLOCK - 1) / CC_GROUP_BLOCK;
        int unordered = 0;
        #pragma omp parallel for schedule(dynamic, 4) num_threads(team) if (team > 1)
        for (size_t b = 0; b < nblocks; b++)
        {
            int seen;
            #pragma omp atomic read
            seen = unordered;
            if (seen)
                continue;
            size_t lo = probe + b * CC_GROUP_BLOCK;
            size_t hi = lo + CC_GROUP_BLOCK < n ? lo + CC_GROUP_BLOCK : n;
            if (!cc_mark_starts(d, nby, lo, hi, start))
            {
                #pragma omp atomic write
                unordered = 1;
            }
        }
        ordered = !unordered;
    }
    if (!ordered)
    {
        free(start);
        return STATA_OK;
    }

    size_t ngroups = 0;
    for (size_t i = 0; i < n; i++)
        ngroups += start[i];
    out->gid = malloc(n * sizeof(*out->gid));
    out->first = malloc(ngroups * sizeof(*out->first));
    out->starts = malloc((ngroups + 1) * sizeof(*out->starts));
    if (!out->gid || !out->first || !out->starts)
    {
        free(start);
        return STATA_ERR_MEMORY;
    }
    size_t g = 0;
    for (size_t i = 0; i < n; i++)
    {
        if (start[i])
        {
            out->first[g] = (perm_idx_t)i;
            out->starts[g] = i;
            g++;
        }
        out->gid[i] = (uint32_t)(g - 1);
    }
    out->starts[ngroups] = n;
    out->ngroups = ngroups;
    out->contiguous = 1;
    out->method = "ordered";
    out->method_code = 2;
    free(start);
    return STATA_OK;
}

/* ===========================================================================
   Dense integer keys
   =========================================================================== */

typedef struct
{
    double min, max; /* over the nonmissing values */
    int any;         /* a nonmissing value was seen */
    int missing;     /* a missing value was seen */
    int usable;      /* whole numbers and canonical missing values only */
} cc_key_range;

static void cc_scan_key(const double *x, size_t lo, size_t hi, double missval, cc_key_range *r)
{
    double min = r->min, max = r->max;
    int any = r->any, missing = r->missing, usable = r->usable;
    for (size_t i = lo; i < hi && usable; i++)
    {
        double v = x[i];
        if (v < missval)
        {
            /* whole numbers of magnitude below 2^52 convert exactly */
            if (!(v > -4503599627370496.0 && v < 4503599627370496.0) ||
                (double)(int64_t)v != v)
                usable = 0;
            if (!any || v < min)
                min = v;
            if (!any || v > max)
                max = v;
            any = 1;
        }
        else
        {
            uint64_t bits;
            memcpy(&bits, &v, sizeof(bits));
            uint64_t offset = bits - CC_MISSING_BITS;
            if ((offset & UINT64_C(0xFFFFFFFFFF)) || (offset >> 40) >= CC_MISSING_CODES)
                usable = 0;
            missing = 1;
        }
    }
    r->min = min;
    r->max = max;
    r->any = any;
    r->missing = missing;
    r->usable = usable;
}

static cc_key_range cc_key_extent(const double *x, size_t n, double missval, int threads)
{
    cc_key_range total = {0.0, 0.0, 0, 0, 1};
    int team = cc_team(threads, n);
    if (team <= 1)
    {
        cc_scan_key(x, 0, n, missval, &total);
        return total;
    }
    size_t nblocks = (n + CC_GROUP_BLOCK - 1) / CC_GROUP_BLOCK;
    #pragma omp parallel num_threads(team)
    {
        cc_key_range local = {0.0, 0.0, 0, 0, 1};
        #pragma omp for schedule(static) nowait
        for (size_t b = 0; b < nblocks; b++)
        {
            size_t lo = b * CC_GROUP_BLOCK;
            size_t hi = lo + CC_GROUP_BLOCK < n ? lo + CC_GROUP_BLOCK : n;
            cc_scan_key(x, lo, hi, missval, &local);
        }
        #pragma omp critical(ccollapse_key_extent)
        {
            if (local.any)
            {
                if (!total.any || local.min < total.min)
                    total.min = local.min;
                if (!total.any || local.max > total.max)
                    total.max = local.max;
                total.any = 1;
            }
            total.missing |= local.missing;
            total.usable &= local.usable;
        }
    }
    return total;
}

/*
    lowest[c] = the lowest row whose code is c, CC_NONE when no row has it.
    Workers take contiguous shares of the rows and private tables while those
    stay small; the first share that saw a code holds its lowest row.
*/
static stata_retcode cc_lowest_rows(const uint32_t *code, size_t n, size_t ncodes,
                                    perm_idx_t *lowest, int threads)
{
    int team = cc_team(threads, n);
    if (team > 1 && ncodes > ((size_t)4 << 20) / (size_t)team)
        team = 1;
    perm_idx_t *shared = NULL;
    if (team > 1)
    {
        shared = malloc((size_t)team * ncodes * sizeof(*shared));
        if (!shared)
            team = 1;
    }
    if (team <= 1)
    {
        memset(lowest, 0xFF, ncodes * sizeof(*lowest));
        for (size_t i = 0; i < n; i++)
        {
            if (lowest[code[i]] == CC_NONE)
                lowest[code[i]] = (perm_idx_t)i;
        }
        return STATA_OK;
    }
    memset(shared, 0xFF, (size_t)team * ncodes * sizeof(*shared));
    /* A static schedule of one iteration per share keeps the shares in row
     * order whatever team the runtime grants */
    #pragma omp parallel for schedule(static, 1) num_threads(team)
    for (int t = 0; t < team; t++)
    {
        perm_idx_t *mine = shared + (size_t)t * ncodes;
        size_t lo = n * (size_t)t / (size_t)team, hi = n * (size_t)(t + 1) / (size_t)team;
        for (size_t i = lo; i < hi; i++)
        {
            if (mine[code[i]] == CC_NONE)
                mine[code[i]] = (perm_idx_t)i;
        }
    }
    #pragma omp parallel for schedule(static) num_threads(team)
    for (size_t c = 0; c < ncodes; c++)
    {
        perm_idx_t row = CC_NONE;
        for (int t = 0; t < team && row == CC_NONE; t++)
            row = shared[(size_t)t * ncodes + c];
        lowest[c] = row;
    }
    free(shared);
    return STATA_OK;
}

/*
    Groups of whole-number keys in a small box. Returns STATA_OK with
    out->ngroups == 0 when the keys do not qualify.
*/
static stata_retcode cc_group_dense(const stata_data *d, size_t nby, int threads, cc_groups *out)
{
    enum { CC_DENSE_MAX_KEYS = 8 };
    size_t n = d->nobs;
    const double missval = SV_missval;
    double base[CC_DENSE_MAX_KEYS];
    size_t span[CC_DENSE_MAX_KEYS], width[CC_DENSE_MAX_KEYS];

    if (nby > CC_DENSE_MAX_KEYS)
        return STATA_OK;
    for (size_t k = 0; k < nby; k++)
        if (d->vars[k].type != STATA_TYPE_DOUBLE)
            return STATA_OK;

    /* The cell table may be as large as the row count: beyond that a hash
     * table of the distinct keys is the smaller structure */
    size_t limit = n > 65536 ? n : 65536;
    if (limit > UINT32_MAX - 1)
        limit = UINT32_MAX - 1;
    size_t cells = 1;
    for (size_t k = 0; k < nby; k++)
    {
        cc_key_range r = cc_key_extent(d->vars[k].data.dbl, n, missval, threads);
        if (!r.usable)
            return STATA_OK;
        double extent = r.any ? r.max - r.min + 1.0 : 0.0;
        if (extent > (double)limit)
            return STATA_OK;
        base[k] = r.any ? r.min : 0.0;
        span[k] = (size_t)extent;
        width[k] = span[k] + (r.missing ? CC_MISSING_CODES : 0);
        if (width[k] == 0 || cells > limit / width[k])
            return STATA_OK;
        cells *= width[k];
    }

    stata_retcode rc = STATA_OK;
    uint32_t *code = malloc(n * sizeof(*code));
    perm_idx_t *lowest = malloc(cells * sizeof(*lowest));
    if (!code || !lowest)
    {
        free(code);
        free(lowest);
        return STATA_ERR_MEMORY;
    }

    int team = cc_team(threads, n);
    /* The first key is the most significant digit, as in a sort */
    #pragma omp parallel for schedule(static) num_threads(team) if (team > 1)
    for (size_t i = 0; i < n; i++)
    {
        size_t c = 0;
        for (size_t k = 0; k < nby; k++)
        {
            double v = d->vars[k].data.dbl[i];
            size_t digit;
            if (v < missval)
                digit = (size_t)(v - base[k]);
            else
            {
                uint64_t bits;
                memcpy(&bits, &v, sizeof(bits));
                digit = span[k] + (size_t)((bits - CC_MISSING_BITS) >> 40);
            }
            c = c * width[k] + digit;
        }
        code[i] = (uint32_t)c;
    }

    rc = cc_lowest_rows(code, n, cells, lowest, threads);
    if (rc != STATA_OK)
    {
        free(code);
        free(lowest);
        return rc;
    }

    /* Occupied cells in order are the groups in sort order; lowest[] turns
     * from a row into the cell's group number */
    size_t ngroups = 0;
    for (size_t c = 0; c < cells; c++)
        ngroups += lowest[c] != CC_NONE;
    out->first = malloc((ngroups ? ngroups : 1) * sizeof(*out->first));
    if (!out->first)
    {
        free(code);
        free(lowest);
        return STATA_ERR_MEMORY;
    }
    size_t g = 0;
    for (size_t c = 0; c < cells; c++)
    {
        if (lowest[c] == CC_NONE)
            continue;
        out->first[g] = lowest[c];
        lowest[c] = (perm_idx_t)g++;
    }
    #pragma omp parallel for schedule(static) num_threads(team) if (team > 1)
    for (size_t i = 0; i < n; i++)
        code[i] = lowest[code[i]];

    free(lowest);
    out->gid = code;
    out->ngroups = ngroups;
    out->contiguous = 0;
    out->method = "dense";
    out->method_code = 3;
    return STATA_OK;
}

/* ===========================================================================
   Hashed keys
   =========================================================================== */

static inline uint64_t cc_mix(uint64_t x)
{
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

static inline uint64_t cc_hash_keys(const stata_data *d, size_t nby, size_t row)
{
    uint64_t h = UINT64_C(0x9e3779b97f4a7c15);
    for (size_t k = 0; k < nby; k++)
    {
        const stata_variable *v = &d->vars[k];
        uint64_t bits = 0;
        if (v->type == STATA_TYPE_DOUBLE)
        {
            double x = v->data.dbl[row];
            /* negative and positive zero are one key */
            if (x != 0)
                memcpy(&bits, &x, sizeof(bits));
        }
        else
        {
            bits = UINT64_C(14695981039346656037);
            for (const unsigned char *s = (const unsigned char *)cc_str(v, row); *s; s++)
                bits = (bits ^ *s) * UINT64_C(1099511628211);
        }
        h = cc_mix(h ^ bits) + UINT64_C(0x9e3779b97f4a7c15);
    }
    return h;
}

/* Open addressing, linear probing, at most half full. A slot holds the
 * number of its key plus one (0: empty) and the key's hash. */
typedef struct
{
    uint64_t hash;
    uint32_t id1;
} cc_slot;

typedef struct
{
    cc_slot *slots;
    size_t capacity;     /* power of two */
    size_t count;        /* keys numbered so far */
    perm_idx_t *rows;    /* [count] lowest row of each key */
    uint64_t *hashes;    /* [count] */
    size_t reserved;     /* entries rows and hashes can hold */
} cc_table;

static void cc_table_free(cc_table *t)
{
    free(t->slots);
    free(t->rows);
    free(t->hashes);
    memset(t, 0, sizeof(*t));
}

static int cc_table_init(cc_table *t, size_t capacity)
{
    memset(t, 0, sizeof(*t));
    t->slots = calloc(capacity, sizeof(*t->slots));
    t->reserved = capacity / 2;
    t->rows = malloc(t->reserved * sizeof(*t->rows));
    t->hashes = malloc(t->reserved * sizeof(*t->hashes));
    if (!t->slots || !t->rows || !t->hashes)
    {
        cc_table_free(t);
        return -1;
    }
    t->capacity = capacity;
    return 0;
}

static int cc_table_grow(cc_table *t)
{
    if (t->capacity > SIZE_MAX / 2 / sizeof(*t->slots))
        return -1;
    size_t capacity = t->capacity * 2;
    cc_slot *slots = calloc(capacity, sizeof(*slots));
    perm_idx_t *rows = realloc(t->rows, (capacity / 2) * sizeof(*rows));
    if (rows)
        t->rows = rows;
    uint64_t *hashes = realloc(t->hashes, (capacity / 2) * sizeof(*hashes));
    if (hashes)
        t->hashes = hashes;
    if (!slots || !rows || !hashes)
    {
        free(slots);
        return -1;
    }
    for (size_t k = 0; k < t->count; k++)
    {
        size_t s = t->hashes[k] & (capacity - 1);
        while (slots[s].id1)
            s = (s + 1) & (capacity - 1);
        slots[s].hash = t->hashes[k];
        slots[s].id1 = (uint32_t)(k + 1);
    }
    free(t->slots);
    t->slots = slots;
    t->capacity = capacity;
    t->reserved = capacity / 2;
    return 0;
}

/* The number of row's key in t, adding the key when it is new. -1: no memory */
static inline int64_t cc_table_number(cc_table *t, const stata_data *d, size_t nby, size_t row,
                                      uint64_t hash)
{
    size_t mask = t->capacity - 1, s = hash & mask;
    while (t->slots[s].id1)
    {
        if (t->slots[s].hash == hash &&
            cc_equal_keys(d, nby, t->rows[t->slots[s].id1 - 1], row))
            return (int64_t)t->slots[s].id1 - 1;
        s = (s + 1) & mask;
    }
    if (t->count >= t->reserved)
    {
        if (cc_table_grow(t))
            return -1;
        mask = t->capacity - 1;
        s = hash & mask;
        while (t->slots[s].id1)
            s = (s + 1) & mask;
    }
    t->slots[s].hash = hash;
    t->slots[s].id1 = (uint32_t)(t->count + 1);
    t->rows[t->count] = (perm_idx_t)row;
    t->hashes[t->count] = hash;
    return (int64_t)t->count++;
}

/*
    Sort order of the distinct keys: order[r] is the key (an index into rows)
    of rank r. The keys are gathered into columns of their own and ordered by
    the shared stable sorter.
*/
static stata_retcode cc_order_keys(const stata_data *d, size_t nby, const perm_idx_t *rows,
                                   size_t nkeys, perm_idx_t *order)
{
    stata_data keys;
    stata_data_init(&keys);
    stata_retcode rc = STATA_OK;
    int *columns = malloc(nby * sizeof(*columns));
    keys.vars = calloc(nby, sizeof(*keys.vars));
    if (!columns || !keys.vars)
        rc = STATA_ERR_MEMORY;
    for (size_t k = 0; k < nby && rc == STATA_OK; k++)
    {
        const stata_variable *v = &d->vars[k];
        stata_variable *out = &keys.vars[k];
        columns[k] = (int)k;
        out->type = v->type;
        out->nobs = nkeys;
        if (v->type == STATA_TYPE_DOUBLE)
        {
            out->data.dbl = malloc((nkeys ? nkeys : 1) * sizeof(double));
            if (!out->data.dbl)
                rc = STATA_ERR_MEMORY;
            else
                for (size_t j = 0; j < nkeys; j++)
                    out->data.dbl[j] = v->data.dbl[rows[j]];
        }
        else
        {
            out->data.str = malloc((nkeys ? nkeys : 1) * sizeof(char *));
            if (!out->data.str)
                rc = STATA_ERR_MEMORY;
            else
                for (size_t j = 0; j < nkeys; j++)
                    out->data.str[j] = (char *)cc_str(v, rows[j]);
        }
    }
    if (rc == STATA_OK)
    {
        keys.nobs = nkeys;
        keys.nvars = nby;
        keys.sort_order = order;
        for (size_t j = 0; j < nkeys; j++)
            order[j] = (perm_idx_t)j;
        rc = ctools_order_stable(&keys, columns, nby);
    }
    /* The strings are borrowed and order belongs to the caller */
    if (keys.vars)
    {
        for (size_t k = 0; k < nby; k++)
            free(d->vars[k].type == STATA_TYPE_DOUBLE ? (void *)keys.vars[k].data.dbl
                                                      : (void *)keys.vars[k].data.str);
    }
    free(keys.vars);
    free(columns);
    return rc;
}

/*
    Groups through hash tables. Returns STATA_OK with out->ngroups == 0 when
    a share of the rows holds more distinct keys than is worth hashing.
*/
static stata_retcode cc_group_hash(const stata_data *d, size_t nby, int threads, cc_groups *out)
{
    size_t n = d->nobs;
    int team = cc_team(threads, n);
    stata_retcode rc = STATA_OK;
    int failed = 0, crowded = 0;

    /* Numbering more keys than this per share costs more than sorting rows */
    size_t share = n / (size_t)team;
    size_t cap = share / 4 < n / 64 ? share / 4 : n / 64;
    if (cap < 4096)
        cap = 4096;

    uint32_t *gid = malloc(n * sizeof(*gid));
    cc_table *tables = calloc((size_t)team, sizeof(*tables));
    uint32_t **number = calloc((size_t)team, sizeof(*number));
    cc_table merged;
    memset(&merged, 0, sizeof(merged));
    perm_idx_t *order = NULL;
    uint32_t *rank = NULL;
    if (!gid || !tables || !number)
    {
        rc = STATA_ERR_MEMORY;
        goto done;
    }

    #pragma omp parallel for schedule(static, 1) num_threads(team) if (team > 1)
    for (int t = 0; t < team; t++)
    {
        cc_table *table = &tables[t];
        size_t lo = n * (size_t)t / (size_t)team, hi = n * (size_t)(t + 1) / (size_t)team;
        if (cc_table_init(table, 1024))
        {
            #pragma omp atomic write
            failed = 1;
            continue;
        }
        for (size_t i = lo; i < hi; i++)
        {
            int64_t id = cc_table_number(table, d, nby, i, cc_hash_keys(d, nby, i));
            if (id < 0)
            {
                #pragma omp atomic write
                failed = 1;
                break;
            }
            gid[i] = (uint32_t)id;
            if (table->count > cap)
            {
                #pragma omp atomic write
                crowded = 1;
                break;
            }
            if ((i & 4095) == 0)
            {
                int stop;
                #pragma omp atomic read
                stop = crowded;
                if (stop)
                    break;
            }
        }
    }
    if (failed)
    {
        rc = STATA_ERR_MEMORY;
        goto done;
    }
    if (crowded)
        goto done;

    /* Merge the shares in row order, so that a key's row stays its lowest */
    size_t total = 0;
    for (int t = 0; t < team; t++)
        total += tables[t].count;
    size_t capacity = 1024;
    while (capacity / 2 < tables[0].count)
        capacity *= 2;
    if (cc_table_init(&merged, capacity))
    {
        rc = STATA_ERR_MEMORY;
        goto done;
    }
    for (int t = 0; t < team; t++)
    {
        number[t] = malloc((tables[t].count ? tables[t].count : 1) * sizeof(**number));
        if (!number[t])
        {
            rc = STATA_ERR_MEMORY;
            goto done;
        }
        for (size_t k = 0; k < tables[t].count; k++)
        {
            int64_t id = cc_table_number(&merged, d, nby, tables[t].rows[k], tables[t].hashes[k]);
            if (id < 0)
            {
                rc = STATA_ERR_MEMORY;
                goto done;
            }
            number[t][k] = (uint32_t)id;
        }
    }
    (void)total;

    size_t ngroups = merged.count;
    order = malloc((ngroups ? ngroups : 1) * sizeof(*order));
    rank = malloc((ngroups ? ngroups : 1) * sizeof(*rank));
    out->first = malloc((ngroups ? ngroups : 1) * sizeof(*out->first));
    if (!order || !rank || !out->first)
    {
        rc = STATA_ERR_MEMORY;
        goto done;
    }
    rc = cc_order_keys(d, nby, merged.rows, ngroups, order);
    if (rc != STATA_OK)
        goto done;
    for (size_t r = 0; r < ngroups; r++)
    {
        rank[order[r]] = (uint32_t)r;
        out->first[r] = merged.rows[order[r]];
    }
    /* share number -> merged number -> rank, composed once per share */
    for (int t = 0; t < team; t++)
        for (size_t k = 0; k < tables[t].count; k++)
            number[t][k] = rank[number[t][k]];
    #pragma omp parallel for schedule(static, 1) num_threads(team) if (team > 1)
    for (int t = 0; t < team; t++)
    {
        size_t lo = n * (size_t)t / (size_t)team, hi = n * (size_t)(t + 1) / (size_t)team;
        const uint32_t *map = number[t];
        for (size_t i = lo; i < hi; i++)
            gid[i] = map[gid[i]];
    }
    out->gid = gid;
    gid = NULL;
    out->ngroups = ngroups;
    out->contiguous = 0;
    out->method = "hash";
    out->method_code = 4;

done:
    if (tables)
        for (int t = 0; t < team; t++)
            cc_table_free(&tables[t]);
    if (number)
        for (int t = 0; t < team; t++)
            free(number[t]);
    cc_table_free(&merged);
    free(tables);
    free(number);
    free(order);
    free(rank);
    free(gid);
    if (rc != STATA_OK || !out->gid)
    {
        free(out->first);
        out->first = NULL;
    }
    return rc;
}

/* ===========================================================================
   Sorted rows
   =========================================================================== */

static stata_retcode cc_group_sort(const stata_data *d, size_t nby, int threads, cc_groups *out)
{
    size_t n = d->nobs;
    stata_retcode rc = STATA_OK;
    int *columns = malloc(nby * sizeof(*columns));
    perm_idx_t *perm = malloc(n * sizeof(*perm));
    unsigned char *start = malloc(n);
    if (!columns || !perm || !start)
    {
        rc = STATA_ERR_MEMORY;
        goto done;
    }
    for (size_t k = 0; k < nby; k++)
        columns[k] = (int)k;
    for (size_t i = 0; i < n; i++)
        perm[i] = (perm_idx_t)i;

    /* The by() columns are the first columns of d; the view borrows them */
    stata_data view = *d;
    view.sort_order = perm;
    rc = ctools_order_stable(&view, columns, nby);
    if (rc != STATA_OK)
        goto done;

    int team = cc_team(threads, n);
    start[0] = 1;
    #pragma omp parallel for schedule(static) num_threads(team) if (team > 1)
    for (size_t p = 1; p < n; p++)
        start[p] = !cc_equal_keys(d, nby, perm[p - 1], perm[p]);

    size_t ngroups = 0;
    for (size_t p = 0; p < n; p++)
        ngroups += start[p];
    out->gid = malloc(n * sizeof(*out->gid));
    out->first = malloc(ngroups * sizeof(*out->first));
    out->starts = malloc((ngroups + 1) * sizeof(*out->starts));
    if (!out->gid || !out->first || !out->starts)
    {
        rc = STATA_ERR_MEMORY;
        goto done;
    }
    size_t g = 0;
    for (size_t p = 0; p < n; p++)
    {
        if (start[p])
        {
            /* a stable order puts a group's lowest row first */
            out->first[g] = perm[p];
            out->starts[g] = p;
            g++;
        }
        out->gid[perm[p]] = (uint32_t)(g - 1);
    }
    out->starts[ngroups] = n;
    out->ngroups = ngroups;
    out->perm = perm;
    perm = NULL;
    out->contiguous = 0;
    out->method = "sort";
    out->method_code = 5;

done:
    free(columns);
    free(perm);
    free(start);
    return rc;
}

/* ===========================================================================
   Entry points
   =========================================================================== */

stata_retcode cc_group_rows(const stata_data *data, size_t nby, int threads, cc_groups *out)
{
    memset(out, 0, sizeof(*out));
    out->nobs = data->nobs;
    if (data->nobs == 0)
        return STATA_ERR_INVALID_INPUT;
    if (data->nobs >= UINT32_MAX)
        return STATA_ERR_UNSUPPORTED_TYPE;

    stata_retcode rc;
    if (nby == 0)
        rc = cc_group_all(data->nobs, out);
    else
    {
        rc = cc_group_ordered(data, nby, threads, out);
        if (rc == STATA_OK && !out->ngroups)
            rc = cc_group_dense(data, nby, threads, out);
        if (rc == STATA_OK && !out->ngroups)
            rc = cc_group_hash(data, nby, threads, out);
        if (rc == STATA_OK && !out->ngroups)
            rc = cc_group_sort(data, nby, threads, out);
    }
    if (rc != STATA_OK)
        cc_groups_free(out);
    return rc;
}

stata_retcode cc_groups_order(cc_groups *groups, int threads)
{
    size_t n = groups->nobs, ngroups = groups->ngroups;
    if (groups->starts && (groups->contiguous || groups->perm))
        return STATA_OK;

    const uint32_t *gid = groups->gid;
    size_t *starts = groups->starts;
    if (!starts)
    {
        starts = calloc(ngroups + 1, sizeof(*starts));
        if (!starts)
            return STATA_ERR_MEMORY;
        for (size_t i = 0; i < n; i++)
            starts[gid[i] + 1]++;
        for (size_t g = 0; g < ngroups; g++)
            starts[g + 1] += starts[g];
        groups->starts = starts;
    }
    if (groups->contiguous)
        return STATA_OK;

    perm_idx_t *perm = malloc(n * sizeof(*perm));
    if (!perm)
        return STATA_ERR_MEMORY;

    /* A stable counting sort of the rows by group. Workers take contiguous
     * shares of the rows; share t writes each group's rows after those of
     * the shares before it, so that rows stay ascending within a group. */
    int team = cc_team(threads, n);
    if (team > 1 && ngroups > ((size_t)4 << 20) / (size_t)team)
        team = 1;
    size_t *cursor = team > 1 ? calloc((size_t)team * ngroups, sizeof(*cursor)) : NULL;
    if (team > 1 && cursor)
    {
        #pragma omp parallel for schedule(static, 1) num_threads(team)
        for (int t = 0; t < team; t++)
        {
            size_t *mine = cursor + (size_t)t * ngroups;
            size_t lo = n * (size_t)t / (size_t)team, hi = n * (size_t)(t + 1) / (size_t)team;
            for (size_t i = lo; i < hi; i++)
                mine[gid[i]]++;
        }
        #pragma omp parallel for schedule(static) num_threads(team)
        for (size_t g = 0; g < ngroups; g++)
        {
            size_t at = starts[g];
            for (int t = 0; t < team; t++)
            {
                size_t count = cursor[(size_t)t * ngroups + g];
                cursor[(size_t)t * ngroups + g] = at;
                at += count;
            }
        }
        #pragma omp parallel for schedule(static, 1) num_threads(team)
        for (int t = 0; t < team; t++)
        {
            size_t *mine = cursor + (size_t)t * ngroups;
            size_t lo = n * (size_t)t / (size_t)team, hi = n * (size_t)(t + 1) / (size_t)team;
            for (size_t i = lo; i < hi; i++)
                perm[mine[gid[i]]++] = (perm_idx_t)i;
        }
        free(cursor);
    }
    else
    {
        free(cursor);
        cursor = malloc((ngroups ? ngroups : 1) * sizeof(*cursor));
        if (!cursor)
        {
            free(perm);
            return STATA_ERR_MEMORY;
        }
        memcpy(cursor, starts, ngroups * sizeof(*cursor));
        for (size_t i = 0; i < n; i++)
            perm[cursor[gid[i]]++] = (perm_idx_t)i;
        free(cursor);
    }
    groups->perm = perm;
    return STATA_OK;
}

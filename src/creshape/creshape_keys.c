/*
 * creshape_keys.c - order, groups and checks on the staged observations.
 *
 * Everything here reads the staged rows (ctools_rowstage) and never calls
 * the SPI for data, so the loops run on every thread. Rows compare as Stata
 * sorts them: numbers as doubles (the missing values are the largest
 * doubles, in order), strings by their bytes.
 *
 * The result of each parallel loop is the same for any number of threads:
 * flags are combined with |, sets are merged and then sorted.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "stplugin.h"
#include "creshape_internal.h"
#include "ctools_config.h"
#include "ctools_sort.h"

/* Rows below this run on the calling thread. */
#define KEYS_PARALLEL_ROWS 50000
/* A direct table from an integer j value to its level, for ranges this wide. */
#define KEYS_DENSE_RANGE 65536

static int keys_threads(size_t nobs)
{
    int nt = ctools_get_max_threads();
    return nobs < KEYS_PARALLEL_ROWS || nt < 1 ? 1 : nt;
}

static inline const char *keys_row(const ctools_rowstage *rs, size_t row)
{
    return rs->rows + row * rs->stride;
}

static inline double keys_num(const ctools_rowstage *rs, size_t row, int col)
{
    return *(const double *)(const void *)(keys_row(rs, row) + rs->offset[col]);
}

static inline const char *keys_str(const ctools_rowstage *rs, size_t row, int col)
{
    return keys_row(rs, row) + rs->offset[col];
}

static inline int keys_cmp(const ctools_rowstage *rs, size_t a, size_t b, const int *cols,
                           size_t ncols)
{
    const char *ra = keys_row(rs, a), *rb = keys_row(rs, b);
    for (size_t k = 0; k < ncols; k++)
    {
        const size_t off = rs->offset[cols[k]];
        if (rs->width[cols[k]] == 0)
        {
            double x = *(const double *)(const void *)(ra + off);
            double y = *(const double *)(const void *)(rb + off);
            if (x < y)
                return -1;
            if (x > y)
                return 1;
        }
        else
        {
            int c = strcmp(ra + off, rb + off);
            if (c)
                return c;
        }
    }
    return 0;
}

/* Compare rows on the i variables and then on j; *i_differs says whether
 * the i variables already decide (the later row starts a group of i). */
static inline int keys_cmp_ij(const creshape_state *st, size_t a, size_t b, int *i_differs)
{
    int c = keys_cmp(&st->staged, a, b, st->icols, st->ni);
    *i_differs = c != 0;
    return c ? c : keys_cmp(&st->staged, a, b, &st->jcol, 1);
}

/* ------------------------------------------------------------------------
 * Order
 * ------------------------------------------------------------------------ */

int creshape_order_rows(creshape_state *st, int *unique)
{
    const ctools_rowstage *rs = &st->staged;
    const size_t nobs = st->nobs;
    *unique = 1;
    st->order = NULL;
    if (nobs < 2)
        return 0;

    /* the keys: i, and for long data j */
    size_t nkeys = st->ni + (st->dir == CRESHAPE_WIDE ? 1 : 0);
    int *keys = (int *)malloc(nkeys * sizeof(int));
    if (!keys)
        return 920;
    memcpy(keys, st->icols, st->ni * sizeof(int));
    if (st->dir == CRESHAPE_WIDE)
        keys[st->ni] = st->jcol;

    /* Data that are already in order (the usual case for a panel) need no
     * sort, and the same pass shows whether a key repeats. For long data
     * creshape_wide_levels made that pass together with its own. */
    const int nt = keys_threads(nobs);
    int unsorted = st->scan_unsorted, equal = st->scan_equal;
    (void)nt;
    if (!st->first)
    {
#pragma omp parallel for schedule(static) reduction(| : unsorted, equal) num_threads(nt) if (nt > 1)
        for (size_t r = 1; r < nobs; r++)
        {
            if (unsorted)
                continue;
            int c = keys_cmp(rs, r - 1, r, keys, nkeys);
            unsorted |= c > 0;
            equal |= c == 0;
        }
    }
    if (!unsorted)
    {
        *unique = !equal;
        free(keys);
        return 0;
    }

    /* Sort copies of the key columns with the shared sorter (stable, the
     * order csort gives). */
    stata_data d;
    stata_data_init(&d);
    int rc = 0;
    int *sort_vars = (int *)malloc(nkeys * sizeof(int));
    d.vars = (stata_variable *)ctools_safe_cacheline_alloc2(nkeys, sizeof(stata_variable));
    d.sort_order = (perm_idx_t *)ctools_safe_cacheline_alloc2(nobs, sizeof(perm_idx_t));
    if (!sort_vars || !d.vars || !d.sort_order)
        rc = 920;
    if (!rc)
    {
        memset(d.vars, 0, nkeys * sizeof(stata_variable));
        d.nobs = nobs;
        d.nvars = nkeys;
#pragma omp parallel for schedule(static) num_threads(nt) if (nt > 1)
        for (size_t r = 0; r < nobs; r++)
            d.sort_order[r] = (perm_idx_t)r;
        for (size_t k = 0; k < nkeys && !rc; k++)
        {
            sort_vars[k] = (int)k + 1;
            rc = ctools_stata_rc(ctools_rowstage_column(rs, (size_t)keys[k], &d.vars[k]));
        }
    }
    if (!rc)
        rc = ctools_stata_rc(ctools_sort_dispatch(&d, sort_vars, nkeys, SORT_ALG_AUTO));
    if (!rc)
    {
        st->order = d.sort_order;
        d.sort_order = NULL;
    }
    free(sort_vars);
    stata_data_free(&d);
    if (rc)
    {
        free(keys);
        return rc;
    }

    const perm_idx_t *order = st->order;
    unsigned char *first = st->first;
    int repeated = 0;
    if (first)
    {
        /* long data: where i changes in the new order, for the groups */
#pragma omp parallel for schedule(static) reduction(| : repeated) num_threads(nt) if (nt > 1)
        for (size_t p = 1; p < nobs; p++)
        {
            int i_differs;
            repeated |= keys_cmp_ij(st, order[p - 1], order[p], &i_differs) == 0;
            first[p] = (unsigned char)i_differs;
        }
    }
    else
    {
#pragma omp parallel for schedule(static) reduction(| : repeated) num_threads(nt) if (nt > 1)
        for (size_t p = 1; p < nobs; p++)
            repeated |= keys_cmp(rs, order[p - 1], order[p], keys, nkeys) == 0;
    }
    *unique = !repeated;
    free(keys);
    return 0;
}

/* ------------------------------------------------------------------------
 * Wide: the j values in the data
 * ------------------------------------------------------------------------ */

/* A set of 8-byte keys: the bits of a double, or a pointer to a staged
 * string. Open addressing; never holds KEYS_EMPTY. */
#define KEYS_EMPTY UINT64_C(0x7FF8000000C70015)

typedef struct
{
    uint64_t *slots;
    size_t cap, n;
    int strings; /* keys are const char *, compared by content */
} keys_set;

static inline uint64_t keys_mix(uint64_t x)
{
    x ^= x >> 33;
    x *= UINT64_C(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x *= UINT64_C(0xc4ceb9fe1a85ec53);
    return x ^ (x >> 33);
}

static inline uint64_t keys_hash(const keys_set *s, uint64_t key)
{
    if (!s->strings)
        return keys_mix(key);
    uint64_t h = UINT64_C(1469598103934665603);
    for (const unsigned char *p = (const unsigned char *)(uintptr_t)key; *p; p++)
        h = (h ^ *p) * UINT64_C(1099511628211);
    return keys_mix(h);
}

static inline int keys_same(const keys_set *s, uint64_t a, uint64_t b)
{
    return a == b ||
           (s->strings && !strcmp((const char *)(uintptr_t)a, (const char *)(uintptr_t)b));
}

static int keys_set_init(keys_set *s, int strings)
{
    s->cap = 64;
    s->n = 0;
    s->strings = strings;
    s->slots = (uint64_t *)malloc(s->cap * sizeof(uint64_t));
    if (!s->slots)
        return 920;
    for (size_t i = 0; i < s->cap; i++)
        s->slots[i] = KEYS_EMPTY;
    return 0;
}

static int keys_set_add(keys_set *s, uint64_t key);

static int keys_set_grow(keys_set *s)
{
    keys_set bigger = {NULL, s->cap * 2, 0, s->strings};
    if (s->cap > SIZE_MAX / 16)
        return 920;
    bigger.slots = (uint64_t *)malloc(bigger.cap * sizeof(uint64_t));
    if (!bigger.slots)
        return 920;
    for (size_t i = 0; i < bigger.cap; i++)
        bigger.slots[i] = KEYS_EMPTY;
    for (size_t i = 0; i < s->cap; i++)
        if (s->slots[i] != KEYS_EMPTY)
            (void)keys_set_add(&bigger, s->slots[i]);
    free(s->slots);
    *s = bigger;
    return 0;
}

static int keys_set_add(keys_set *s, uint64_t key)
{
    if (2 * (s->n + 1) > s->cap && keys_set_grow(s))
        return 920;
    size_t at = (size_t)keys_hash(s, key) & (s->cap - 1);
    while (s->slots[at] != KEYS_EMPTY)
    {
        if (keys_same(s, s->slots[at], key))
            return 0;
        at = (at + 1) & (s->cap - 1);
    }
    s->slots[at] = key;
    s->n++;
    return 0;
}

static inline uint64_t keys_bits(double v)
{
    uint64_t b;
    if (v == 0)
        v = 0.0; /* -0 and +0 are one value */
    memcpy(&b, &v, sizeof(b));
    return b;
}

static int keys_cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static int keys_cmp_text(const void *a, const void *b)
{
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

/* One pass over the long data for everything that reads each row once:
 * the distinct values of j, in *set; what reshape refuses in j (*flags: bit
 * 0 a missing value, bit 1 leading or trailing blanks); whether the rows
 * are in order of (i, j) and whether a pair repeats; and first[r], set
 * where row r starts a group of i. */
static int keys_scan_j(creshape_state *st, keys_set *set, int *flags)
{
    const ctools_rowstage *rs = &st->staged;
    const size_t nobs = st->nobs;
    const int col = st->jcol, strings = st->string_j;
    const double missing = SV_missval;
    const int nt = keys_threads(nobs);
    unsigned char *first = st->first;
    int bad = 0, failed = 0, unsorted = 0, equal = 0;
    (void)nt;
#pragma omp parallel num_threads(nt) if (nt > 1) reduction(| : bad, failed, unsorted, equal)
    {
        keys_set mine;
        int oom = keys_set_init(&mine, strings);
        uint64_t last = KEYS_EMPTY;
#pragma omp for schedule(static) nowait
        for (size_t r = 0; r < nobs; r++)
        {
            if (r)
            {
                int i_differs;
                int c = keys_cmp_ij(st, r - 1, r, &i_differs);
                first[r] = (unsigned char)i_differs;
                unsorted |= c > 0;
                equal |= c == 0;
            }
            else
                first[0] = 1;
            if (oom)
                continue;
            uint64_t key;
            if (strings)
            {
                const char *s = keys_str(rs, r, col);
                size_t len = strlen(s);
                size_t blank = strspn(s, " ");
                if (blank == len)
                {
                    bad |= 1;
                    continue;
                }
                if (blank || s[len - 1] == ' ')
                    bad |= 2;
                /* runs of one value (sorted data) cost one comparison */
                if (last != KEYS_EMPTY && !strcmp((const char *)(uintptr_t)last, s))
                    continue;
                key = (uint64_t)(uintptr_t)s;
            }
            else
            {
                double v = keys_num(rs, r, col);
                if (v >= missing)
                {
                    bad |= 1;
                    continue;
                }
                key = keys_bits(v);
                if (key == last)
                    continue;
            }
            last = key;
            oom = keys_set_add(&mine, key);
        }
#pragma omp critical(creshape_levels)
        {
            for (size_t i = 0; !oom && mine.slots && i < mine.cap; i++)
                if (mine.slots[i] != KEYS_EMPTY)
                    oom = keys_set_add(set, mine.slots[i]);
        }
        failed |= oom;
        free(mine.slots);
    }
    *flags = bad;
    st->scan_unsorted = unsorted;
    st->scan_equal = equal;
    return failed ? 920 : 0;
}

int creshape_wide_levels(creshape_state *st)
{
    keys_set set;
    int flags = 0;
    if (keys_set_init(&set, st->string_j))
        return 920;
    if (st->nobs && !(st->first = (unsigned char *)malloc(st->nobs)))
    {
        free(set.slots);
        return 920;
    }
    int rc = st->nobs ? keys_scan_j(st, &set, &flags) : 0;
    if (rc)
    {
        free(set.slots);
        return rc;
    }
    if (flags & 1)
    {
        creshape_error("variable {bf:%s} contains missing values", st->jname);
        free(set.slots);
        return 498;
    }
    if (flags & 2)
    {
        creshape_error("variable {bf:%s} has leading or trailing blanks", st->jname);
        free(set.slots);
        return 498;
    }

    /* the values, sorted */
    const size_t n = set.n;
    double *nums = NULL;
    const char **texts = NULL;
    if (st->string_j)
        texts = (const char **)malloc((n ? n : 1) * sizeof(*texts));
    else
        nums = (double *)malloc((n ? n : 1) * sizeof(*nums));
    if (!nums && !texts)
    {
        free(set.slots);
        return 920;
    }
    size_t k = 0;
    for (size_t i = 0; i < set.cap; i++)
    {
        if (set.slots[i] == KEYS_EMPTY)
            continue;
        if (st->string_j)
            texts[k++] = (const char *)(uintptr_t)set.slots[i];
        else
            memcpy(&nums[k++], &set.slots[i], sizeof(double));
    }
    free(set.slots);
    if (st->string_j)
        qsort(texts, n, sizeof(*texts), keys_cmp_text);
    else
        qsort(nums, n, sizeof(*nums), keys_cmp_double);

    if (!st->user_levels)
    {
        st->levels = (creshape_level *)calloc(n ? n : 1, sizeof(creshape_level));
        if (!st->levels)
            rc = 920;
        for (size_t l = 0; l < n && !rc; l++)
        {
            creshape_level *lv = &st->levels[l];
            if (st->string_j)
            {
                if (strchr(texts[l], ' '))
                {
                    char line[512];
                    snprintf(line, sizeof(line),
                             "{p 0 4 2}string j variable {bf:%s} not allowed to contain "
                             "spaces;{break}{bf:encode} string variable {bf:%s} into a numeric "
                             "variable and use the new numeric variable as the j "
                             "variable{p_end}",
                             st->jname, st->jname);
                    creshape_error_text(line);
                    rc = 111;
                    break;
                }
                lv->text = creshape_strdup(texts[l]);
            }
            else
            {
                char number[40];
                creshape_format_number(nums[l], number, sizeof(number));
                lv->num = nums[l];
                lv->text = creshape_strdup(number);
            }
            if (!lv->text)
                rc = 920;
            lv->found = 1;
            st->nlevels = l + 1;
        }
        if (!rc)
        {
            /* (j = 1 2 3) */
            size_t need = 16;
            for (size_t l = 0; l < st->nlevels; l++)
                need += strlen(st->levels[l].text) + 1;
            char *line = (char *)malloc(need);
            if (line)
            {
                size_t at = (size_t)snprintf(line, need, "{txt}(j =");
                for (size_t l = 0; l < st->nlevels; l++)
                    at += (size_t)snprintf(line + at, need - at, " %s", st->levels[l].text);
                snprintf(line + at, need - at, ")\n");
                SF_display(line);
                free(line);
            }
        }
    }
    else
    {
        /* j() listed the values: those without observations are reported
         * and make no variables; the order of the list is kept */
        size_t kept = 0;
        for (size_t l = 0; l < st->nlevels; l++)
        {
            creshape_level *lv = &st->levels[l];
            const void *hit = st->string_j
                                  ? bsearch(&lv->text, texts, n, sizeof(*texts), keys_cmp_text)
                                  : bsearch(&lv->num, nums, n, sizeof(*nums), keys_cmp_double);
            if (!hit)
            {
                creshape_note("no data for %s == %s", st->jname, lv->text);
                free(lv->text);
                continue;
            }
            lv->found = 1;
            st->levels[kept++] = *lv;
        }
        st->nlevels = kept;
    }
    free(nums);
    free(texts);
    return rc ? rc : creshape_levels_report(st);
}

/* ------------------------------------------------------------------------
 * Wide: groups of i, the row of each cell, and constant xi variables
 * ------------------------------------------------------------------------ */

/* From a j value to its position in levels[]. */
typedef struct
{
    const creshape_state *st;
    int32_t *dense; /* owned: numeric j, levels[] position by value - lo, -1 for none */
    double lo;
    size_t range;
} keys_lookup;

static int keys_lookup_init(keys_lookup *lk, const creshape_state *st)
{
    lk->st = st;
    lk->dense = NULL;
    lk->lo = 0;
    lk->range = 0;
    if (st->string_j || !st->nreport)
        return 0;
    const double lo = st->levels[st->report_levels[0]].num;
    const double hi = st->levels[st->report_levels[st->nreport - 1]].num;
    if (!(hi - lo < KEYS_DENSE_RANGE))
        return 0;
    for (size_t l = 0; l < st->nlevels; l++)
        if (st->levels[l].num != (double)(long long)st->levels[l].num)
            return 0;
    lk->range = (size_t)(hi - lo) + 1;
    lk->dense = (int32_t *)malloc(lk->range * sizeof(int32_t));
    if (!lk->dense)
        return 920;
    memset(lk->dense, 0xff, lk->range * sizeof(int32_t));
    for (size_t l = 0; l < st->nlevels; l++)
        lk->dense[(size_t)(st->levels[l].num - lo)] = (int32_t)l;
    lk->lo = lo;
    return 0;
}

static inline long keys_lookup_row(const keys_lookup *lk, const ctools_rowstage *rs, size_t row)
{
    const creshape_state *st = lk->st;
    if (lk->dense)
    {
        double d = keys_num(rs, row, st->jcol) - lk->lo;
        if (!(d >= 0) || d >= (double)lk->range || d != (double)(size_t)d)
            return -1;
        return lk->dense[(size_t)d];
    }
    /* binary search over the levels in order of value */
    size_t a = 0, b = st->nreport;
    if (st->string_j)
    {
        const char *s = keys_str(rs, row, st->jcol);
        while (a < b)
        {
            size_t mid = a + (b - a) / 2;
            int c = strcmp(st->levels[st->report_levels[mid]].text, s);
            if (!c)
                return (long)st->report_levels[mid];
            if (c < 0)
                a = mid + 1;
            else
                b = mid;
        }
        return -1;
    }
    double v = keys_num(rs, row, st->jcol);
    while (a < b)
    {
        size_t mid = a + (b - a) / 2;
        double x = st->levels[st->report_levels[mid]].num;
        if (x == v)
            return (long)st->report_levels[mid];
        if (x < v)
            a = mid + 1;
        else
            b = mid;
    }
    return -1;
}

int creshape_wide_groups(creshape_state *st, unsigned char **varying_out)
{
    const ctools_rowstage *rs = &st->staged;
    const size_t nobs = st->nobs, nlevels = st->nlevels, nxi = st->nxi;
    const perm_idx_t *order = st->order;
    const int nt = keys_threads(nobs);
    *varying_out = NULL;
    st->ngroups = 0;
    (void)nt;

    /* group starts: the positions flagged while the rows were scanned (or
     * sorted), gathered by each thread for its share and laid end to end */
    const unsigned char *first = st->first;
    size_t ngroups = 0;
    size_t *counts = (size_t *)calloc((size_t)nt + 1, sizeof(size_t));
    if (!counts)
        return 920;
    int no_memory = 0;
#pragma omp parallel num_threads(nt) if (nt > 1)
    {
        size_t team = 1, tid = 0;
#ifdef _OPENMP
        team = (size_t)omp_get_num_threads();
        tid = (size_t)omp_get_thread_num();
#endif
        const size_t base = nobs / team, extra = nobs % team;
        const size_t p0 = tid * base + (tid < extra ? tid : extra);
        const size_t p1 = p0 + base + (tid < extra);
        size_t mine = 0;
        for (size_t p = p0; p < p1; p++)
            mine += first[p];
        counts[tid + 1] = mine;
#pragma omp barrier
#pragma omp single
        {
            for (size_t t = 0; t < team; t++)
                counts[t + 1] += counts[t];
            ngroups = counts[team];
            st->group_start = (perm_idx_t *)malloc((ngroups + 1) * sizeof(perm_idx_t));
            if (!st->group_start)
                no_memory = 1;
        }
        if (!no_memory)
        {
            perm_idx_t *out = st->group_start + counts[tid];
            for (size_t p = p0; p < p1; p++)
                if (first[p])
                    *out++ = (perm_idx_t)p;
        }
    }
    free(counts);
    if (no_memory)
        return 920;
    st->group_start[ngroups] = (perm_idx_t)nobs;
    st->ngroups = ngroups;
    st->out_nobs = ngroups;
    free(st->first);
    st->first = NULL;

    size_t ncells;
    if (ctools_safe_mul_size(ngroups, nlevels, &ncells) != 0 ||
        ncells > SIZE_MAX / sizeof(perm_idx_t))
        return 920;
    if (ncells)
    {
        st->cell = (perm_idx_t *)malloc(ncells * sizeof(perm_idx_t));
        if (!st->cell)
            return 920;
    }
    unsigned char *varying = (unsigned char *)calloc(nxi ? nxi : 1, 1);
    keys_lookup lk;
    if (!varying || keys_lookup_init(&lk, st))
    {
        free(varying);
        return 920;
    }

    const perm_idx_t *starts = st->group_start;
    perm_idx_t *cell = st->cell;
    const int *xicols = st->xicols;
    int failed = 0;
#pragma omp parallel num_threads(nt) if (nt > 1) reduction(| : failed)
    {
        unsigned char *mine = (unsigned char *)calloc(nxi ? nxi : 1, 1);
        if (!mine)
            failed = 1;
#pragma omp for schedule(dynamic, 1024) nowait
        for (size_t grp = 0; grp < ngroups; grp++)
        {
            if (!mine)
                continue;
            perm_idx_t *cells = cell ? cell + grp * nlevels : NULL;
            for (size_t l = 0; l < nlevels; l++)
                cells[l] = CTOOLS_ROW_NONE;
            const size_t p0 = starts[grp], p1 = starts[grp + 1];
            const size_t base = order ? order[p0] : p0;
            const char *base_row = keys_row(rs, base);
            for (size_t p = p0; p < p1; p++)
            {
                const size_t row = order ? order[p] : p;
                if (cells)
                {
                    long l = keys_lookup_row(&lk, rs, row);
                    if (l >= 0)
                        cells[l] = (perm_idx_t)row;
                }
                if (p == p0)
                    continue;
                /* the other variables must not vary within i */
                const char *this_row = keys_row(rs, row);
                for (size_t x = 0; x < nxi; x++)
                {
                    const size_t off = rs->offset[xicols[x]];
                    if (rs->width[xicols[x]] == 0)
                    {
                        if (*(const double *)(const void *)(this_row + off) !=
                            *(const double *)(const void *)(base_row + off))
                            mine[x] = 1;
                    }
                    else if (strcmp(this_row + off, base_row + off))
                        mine[x] = 1;
                }
            }
        }
#pragma omp critical(creshape_varying)
        {
            for (size_t x = 0; mine && x < nxi; x++)
                varying[x] |= mine[x];
        }
        free(mine);
    }
    free(lk.dense);
    if (failed)
    {
        free(varying);
        return 920;
    }
    *varying_out = varying;
    return 0;
}

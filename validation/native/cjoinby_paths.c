/* cjoinby's key index and output plan against a brute-force join.
 *
 * Covers both lookups (direct and hashed), numeric, string and mixed keys,
 * missing values and -0, every unmatched() mask, and the placement of
 * using-only groups in a sorted master. The sort engines have suites of
 * their own: here the dispatcher is the stable merge order of ctools_order.c.
 */
#undef NDEBUG
#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_order.h"
#include "cjoinby/cjoinby_internal.h"

static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static int threads = 1;
int ctools_get_max_threads(void)
{
    return threads;
}
/* ctools_release_queue comes from src/ctools_release.c, which the harness
 * links into every native test; outside a session it frees nothing. */
void stata_data_init(stata_data *data)
{
    memset(data, 0, sizeof(*data));
}
void stata_data_free(stata_data *data)
{
    /* every key column of these tests is borrowed from the fixture */
    ctools_aligned_free(data->vars);
    ctools_aligned_free(data->sort_order);
    stata_data_init(data);
}
stata_retcode ctools_sort_dispatch(stata_data *data, int *sort_vars, size_t nsort,
                                   sort_algorithm_t algorithm)
{
    int keys[8];
    (void)algorithm;
    assert(nsort <= 8);
    for (size_t k = 0; k < nsort; k++)
        keys[k] = sort_vars[k] - 1;
    return ctools_order_stable(data, keys, nsort);
}

#include "ctools_order.c"
#include "cjoinby/cjoinby_index.c"
#include "cjoinby/cjoinby_plan.c"

/* ---- fixtures: columns of doubles or strings ---------------------------- */

enum
{
    MAXCOLS = 6
};
typedef struct
{
    size_t n;
    int ncols;
    int is_str[MAXCOLS];
    double *num[MAXCOLS];
    char **str[MAXCOLS];
    stata_variable vars[MAXCOLS];
    int width[MAXCOLS];
    int var_idx[MAXCOLS];
    size_t offset[MAXCOLS];
} fixture;

static uint64_t rng_state = 88172645463325252ULL;
static uint32_t rng(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (uint32_t)(rng_state >> 16);
}

static double missing_value(int code) /* 0 for ., 1-26 for .a-.z */
{
    uint64_t bits = UINT64_C(0x7FE0000000000000) + ((uint64_t)code << 40);
    double x;
    memcpy(&x, &bits, sizeof(x));
    return x;
}

static void fixture_init(fixture *f, size_t n, int ncols, const int *is_str)
{
    memset(f, 0, sizeof(*f));
    f->n = n;
    f->ncols = ncols;
    for (int c = 0; c < ncols; c++)
    {
        f->is_str[c] = is_str[c];
        f->var_idx[c] = c + 1;
        f->width[c] = is_str[c] ? 12 : 0;
        if (is_str[c])
        {
            f->str[c] = calloc(n ? n : 1, sizeof(char *));
            assert(f->str[c]);
            for (size_t i = 0; i < n; i++)
            {
                f->str[c][i] = calloc(1, 13);
                assert(f->str[c][i]);
            }
        }
        else
        {
            f->num[c] = calloc(n ? n : 1, sizeof(double));
            assert(f->num[c]);
        }
    }
}

static void fixture_bind(fixture *f)
{
    for (int c = 0; c < f->ncols; c++)
    {
        stata_variable *v = &f->vars[c];
        memset(v, 0, sizeof(*v));
        v->nobs = f->n;
        v->type = f->is_str[c] ? STATA_TYPE_STRING : STATA_TYPE_DOUBLE;
        v->buffer_storage = CTOOLS_BUFFER_BORROWED;
        v->string_storage = CTOOLS_STRINGS_BORROWED;
        if (f->is_str[c])
        {
            v->data.str = f->str[c];
            v->str_maxlen = 12;
        }
        else
        {
            v->data.dbl = f->num[c];
        }
    }
}

static void fixture_free(fixture *f)
{
    for (int c = 0; c < f->ncols; c++)
    {
        if (f->is_str[c])
        {
            for (size_t i = 0; i < f->n; i++)
                free(f->str[c][i]);
            free(f->str[c]);
        }
        else
        {
            free(f->num[c]);
        }
    }
    memset(f, 0, sizeof(*f));
}

/* The using keys as cjoinby's table, held as columns */
static void fixture_table(fixture *f, cj_table *t)
{
    fixture_bind(f);
    memset(t, 0, sizeof(*t));
    t->rowmajor = 0;
    t->nobs = f->n;
    t->nvars = (size_t)f->ncols;
    t->stride = 8;
    t->var_idx = f->var_idx;
    t->width = f->width;
    t->offset = f->offset;
    t->cols.data.nobs = f->n;
    t->cols.data.nvars = (size_t)f->ncols;
    t->cols.data.vars = f->vars;
}

stata_retcode cj_table_column(const cj_table *t, size_t col, stata_variable *out)
{
    *out = t->cols.data.vars[col];
    return STATA_OK;
}

static int cell_cmp(const fixture *a, int ca, size_t i, const fixture *b, int cb, size_t j)
{
    if (a->is_str[ca])
        return strcmp(a->str[ca][i], b->str[cb][j]);
    double u = a->num[ca][i], v = b->num[cb][j];
    return (u > v) - (u < v);
}

static int keys_cmp(const fixture *a, size_t i, const fixture *b, size_t j, int nk)
{
    for (int k = 0; k < nk; k++)
    {
        int c = cell_cmp(a, k, i, b, k, j);
        if (c)
            return c;
    }
    return 0;
}

/* Key values of scenario `kind` in the first nk columns */
enum
{
    K_INT,
    K_TWO,
    K_FRACTION,
    K_SPARSE,
    K_STRING,
    K_MIXED,
    K_KINDS
};

static int kind_keys(int kind)
{
    return kind == K_TWO || kind == K_MIXED ? 2 : 1;
}

static void fill_keys(fixture *f, int kind, uint32_t spread)
{
    for (size_t i = 0; i < f->n; i++)
    {
        uint32_t r = rng() % spread, q = rng() % 100;
        switch (kind)
        {
        case K_INT:
            f->num[0][i] = q < 6 ? missing_value(0) : q < 10 ? missing_value(1) : (double)r - 3;
            break;
        case K_TWO:
            f->num[0][i] = q < 5 ? missing_value(26) : (double)(r % 6) - 2;
            f->num[1][i] = (double)(rng() % 7);
            break;
        case K_FRACTION:
            f->num[0][i] = q < 5 ? -0.0 : q < 10 ? 0.0 : q < 14 ? missing_value(2) : r / 4.0 + .1;
            break;
        case K_SPARSE:
            f->num[0][i] = q < 5 ? -5.0 : 1e9 + 7919.0 * r;
            break;
        case K_STRING:
            if (q < 6)
                f->str[0][i][0] = '\0';
            else
                snprintf(f->str[0][i], 13, q < 12 ? "k%u " : "k%u", (unsigned)r);
            break;
        case K_MIXED:
            snprintf(f->str[0][i], 13, q < 6 ? "" : "g%u", (unsigned)(r % 4));
            f->num[1][i] = (double)(rng() % 5) + .5;
            break;
        default:
            abort();
        }
    }
}

static void kind_types(int kind, int *is_str, int extra)
{
    const int nk = kind_keys(kind);
    for (int c = 0; c < MAXCOLS; c++)
        is_str[c] = 0;
    if (kind == K_STRING || kind == K_MIXED)
        is_str[0] = 1;
    if (extra)
        is_str[nk + 1] = 1; /* nk: a number, nk + 1: a string */
}

/* ---- the index ----------------------------------------------------------- */

static void check_index(const cj_index *ix, const fixture *u, int nk, int expect_direct)
{
    const size_t n = u->n;
    const perm_idx_t *order = ix->keys.sort_order;
    unsigned char *seen = calloc(n ? n : 1, 1);
    assert(seen);
    size_t g = 0;
    assert(ix->nobs == n && ix->gstart[0] == 0);
    for (size_t p = 0; p < n; p++)
    {
        assert(order[p] < n && !seen[order[p]]);
        seen[order[p]] = 1;
        if (p == 0)
            continue;
        int c = keys_cmp(u, order[p - 1], u, order[p], nk);
        assert(c <= 0);
        if (c == 0)
            assert(order[p - 1] < order[p]); /* stable */
        else
            assert(ix->gstart[++g] == p);
    }
    assert(ix->ngroups == (n ? g + 1 : 0));
    assert(ix->gstart[ix->ngroups] == n);
    if (expect_direct >= 0 && n)
        assert(ix->use_direct == expect_direct);
    free(seen);
}

/* ---- the plan, without placement ----------------------------------------- */

static void check_plan(const cj_index *ix, const fixture *u, const fixture *m, int nk,
                       int unmatched)
{
    stata_data mk = {m->n, (size_t)m->ncols, (stata_variable *)m->vars, NULL};
    cj_plan p;
    assert(cj_plan_build(&p, ix, m->n ? &mk : NULL, m->n, unmatched, NULL) == 0);
    const int keep_master = unmatched & CJ_KEEP_MASTER, keep_using = unmatched & CJ_KEEP_USING;
    const perm_idx_t *order = ix->keys.sort_order;
    unsigned char *hit = calloc(ix->ngroups ? ix->ngroups : 1, 1);
    assert(hit && p.nm == m->n && p.off[0] == 0 && !p.row && !p.placed);
    size_t both = 0, lone = 0, stage_from = m->n;
    for (size_t i = 0; i < m->n; i++)
    {
        size_t count = 0, first = 0;
        for (size_t j = 0; j < u->n; j++)
        {
            if (keys_cmp(m, i, u, j, nk))
                continue;
            if (!count++)
                first = j;
        }
        const uint32_t g = p.gid[i];
        if (!count)
        {
            assert(g == CJ_NONE);
            lone++;
            count = keep_master ? 1 : 0;
        }
        else
        {
            assert(g != CJ_NONE && cj_group_count(ix, g) == count);
            assert(order[ix->gstart[g]] == first);
            for (size_t k = 0; k < count; k++)
                assert(!keys_cmp(m, i, u, order[ix->gstart[g] + k], nk));
            hit[g] = 1;
            both += count;
        }
        assert(p.off[i + 1] - p.off[i] == count);
        if (stage_from == m->n && count != 1)
            stage_from = i;
        assert(cj_plan_master_row(&p, i) == i);
    }
    size_t d = m->n, tail = 0;
    for (size_t g = 0; g < ix->ngroups && keep_using; g++)
    {
        if (hit[g])
            continue;
        assert(d < p.ndrv && p.gid[d] == g && cj_plan_master_row(&p, d) == CJ_NONE);
        assert(p.off[d + 1] - p.off[d] == cj_group_count(ix, (uint32_t)g));
        tail += cj_group_count(ix, (uint32_t)g);
        d++;
    }
    assert(d == p.ndrv && p.nout == p.off[p.ndrv]);
    assert(p.n_both == both && p.n_using == tail && p.n_master == (keep_master ? lone : 0));
    assert(p.nout == both + tail + (keep_master ? lone : 0));
    assert(p.stage_from == stage_from);
    free(hit);
    cj_plan_free(&p);
}

/* ---- placement in a sorted master ---------------------------------------- */

static const fixture *sort_fixture;
static const int *sort_cols;
static int sort_ncols;

static int row_cmp(const void *a, const void *b)
{
    const size_t i = *(const size_t *)a, j = *(const size_t *)b;
    for (int s = 0; s < sort_ncols; s++)
    {
        int c = cell_cmp(sort_fixture, sort_cols[s], i, sort_fixture, sort_cols[s], j);
        if (c)
            return c;
    }
    return (i > j) - (i < j);
}

/* Reorder the rows of f by the columns cols */
static void fixture_sort(fixture *f, const int *cols, int ncols)
{
    size_t *rows = malloc((f->n ? f->n : 1) * sizeof(size_t));
    assert(rows);
    for (size_t i = 0; i < f->n; i++)
        rows[i] = i;
    sort_fixture = f;
    sort_cols = cols;
    sort_ncols = ncols;
    qsort(rows, f->n, sizeof(size_t), row_cmp);
    for (int c = 0; c < f->ncols; c++)
    {
        if (f->is_str[c])
        {
            char **moved = malloc((f->n ? f->n : 1) * sizeof(char *));
            assert(moved);
            for (size_t i = 0; i < f->n; i++)
                moved[i] = f->str[c][rows[i]];
            free(f->str[c]);
            f->str[c] = moved;
        }
        else
        {
            double *moved = malloc((f->n ? f->n : 1) * sizeof(double));
            assert(moved);
            for (size_t i = 0; i < f->n; i++)
                moved[i] = f->num[c][rows[i]];
            free(f->num[c]);
            f->num[c] = moved;
        }
    }
    free(rows);
    fixture_bind(f);
}

/* Sort variable s of driver d: the master's value, or for a using-only group
 * its key value, and a missing value where the variable is not a key */
static int driver_cmp(const cj_plan *p, const cj_index *ix, const fixture *u, const fixture *m,
                      int nk, const int *cols, int ncols, size_t a, size_t b)
{
    const perm_idx_t *order = ix->keys.sort_order;
    for (int s = 0; s < ncols; s++)
    {
        const int col = cols[s];
        const uint32_t ra = cj_plan_master_row(p, a), rb = cj_plan_master_row(p, b);
        int c;
        if (m->is_str[col])
        {
            const char *x = ra != CJ_NONE ? m->str[col][ra]
                            : col < nk    ? u->str[col][order[ix->gstart[p->gid[a]]]]
                                          : "";
            const char *y = rb != CJ_NONE ? m->str[col][rb]
                            : col < nk    ? u->str[col][order[ix->gstart[p->gid[b]]]]
                                          : "";
            c = strcmp(x, y);
        }
        else
        {
            double x = ra != CJ_NONE ? m->num[col][ra]
                       : col < nk    ? u->num[col][order[ix->gstart[p->gid[a]]]]
                                     : mock.missval;
            double y = rb != CJ_NONE ? m->num[col][rb]
                       : col < nk    ? u->num[col][order[ix->gstart[p->gid[b]]]]
                                     : mock.missval;
            c = (x > y) - (x < y);
        }
        if (c)
            return c;
    }
    return 0;
}

static void check_placement(const cj_index *ix, const fixture *u, fixture *m, int nk,
                            const int *cols, int ncols)
{
    fixture_sort(m, cols, ncols);
    stata_data mk = {m->n, (size_t)m->ncols, m->vars, NULL};
    cj_sortspec spec = {ncols, cols};
    cj_plan p;
    assert(cj_plan_build(&p, ix, m->n ? &mk : NULL, m->n, CJ_KEEP_MASTER | CJ_KEEP_USING, &spec) ==
           0);
    unsigned char *seen = calloc(ix->ngroups ? ix->ngroups : 1, 1);
    assert(seen);
    size_t next_master = 0, tails = 0, tail_rows = 0;
    for (size_t d = 0; d < p.ndrv; d++)
    {
        const uint32_t row = cj_plan_master_row(&p, d);
        if (row != CJ_NONE)
        {
            assert(row == next_master++);
        }
        else
        {
            assert(p.gid[d] < ix->ngroups && !seen[p.gid[d]]);
            seen[p.gid[d]] = 1;
            tails++;
            tail_rows += cj_group_count(ix, p.gid[d]);
        }
        if (d == 0)
            continue;
        const int c = driver_cmp(&p, ix, u, m, nk, cols, ncols, d - 1, d);
        const uint32_t before = cj_plan_master_row(&p, d - 1);
        assert(c <= 0);
        if (c == 0)
        {
            /* ties: master observations first, then groups in key order */
            assert(!(before == CJ_NONE && row != CJ_NONE));
            if (before == CJ_NONE && row == CJ_NONE)
                assert(p.gid[d - 1] < p.gid[d]);
        }
    }
    assert(next_master == m->n && p.ndrv == m->n + tails && p.n_using == tail_rows);
    assert(!!p.placed == (tails > 0) && (p.row != NULL) == (tails > 0));
    /* the same join as the unplaced layout */
    cj_plan plain;
    assert(cj_plan_build(&plain, ix, m->n ? &mk : NULL, m->n, CJ_KEEP_MASTER | CJ_KEEP_USING,
                         NULL) == 0);
    assert(plain.nout == p.nout && plain.n_both == p.n_both && plain.n_master == p.n_master);
    for (size_t d = 0; d < p.ndrv; d++)
    {
        const uint32_t row = cj_plan_master_row(&p, d);
        if (row != CJ_NONE)
            assert(p.gid[d] == plain.gid[row]);
        else
            assert(!keys_cmp(u, ix->keys.sort_order[ix->gstart[p.gid[d]]], u,
                             ix->keys.sort_order[ix->gstart[p.gid[d]]], nk));
        if (d >= p.stage_from)
            continue;
        assert(row == d && p.off[d] == d && p.off[d + 1] == d + 1);
    }
    cj_plan_free(&plain);
    cj_plan_free(&p);
    free(seen);
}

/* ---- scenarios ------------------------------------------------------------ */

static void small_joins(void)
{
    for (int kind = 0; kind < K_KINDS; kind++)
    {
        const int nk = kind_keys(kind);
        for (int round = 0; round < 6; round++)
        {
            const size_t nu = round == 0 ? 0 : 40 + rng() % 220;
            const size_t nm = round == 1 ? 0 : 30 + rng() % 260;
            int types[MAXCOLS];
            fixture u, m;
            cj_table table;
            cj_index ix;
            kind_types(kind, types, 0);
            fixture_init(&u, nu, nk, types);
            fill_keys(&u, kind, 14);
            kind_types(kind, types, 1);
            fixture_init(&m, nm, nk + 2, types);
            fill_keys(&m, kind, 18);
            for (size_t i = 0; i < nm; i++)
            {
                m.num[nk][i] = rng() % 10 < 2 ? missing_value(rng() % 3) : (double)(rng() % 5);
                snprintf(m.str[nk + 1][i], 13, rng() % 5 ? "s%u" : "", (unsigned)(rng() % 4));
            }
            fixture_bind(&m);
            fixture_table(&u, &table);
            assert(cj_index_build(&ix, &table, nk, nm) == STATA_OK);
            check_index(&ix, &u, nk, kind == K_INT || kind == K_TWO ? 1 : 0);
            for (int unmatched = 0; unmatched < 4; unmatched++)
                check_plan(&ix, &u, &m, nk, unmatched);

            /* the master's sort orders: keys, part of them, other variables */
            const int by_keys[2] = {0, 1}, by_other[1] = {nk}, by_text[1] = {nk + 1};
            const int by_first_other[2] = {0, nk}, by_text_first[2] = {nk + 1, 0};
            const int by_last[2] = {nk - 1, nk + 1};
            check_placement(&ix, &u, &m, nk, by_keys, nk);
            check_placement(&ix, &u, &m, nk, by_keys, 1);
            check_placement(&ix, &u, &m, nk, by_other, 1);
            check_placement(&ix, &u, &m, nk, by_text, 1);
            check_placement(&ix, &u, &m, nk, by_first_other, 2);
            check_placement(&ix, &u, &m, nk, by_text_first, 2);
            check_placement(&ix, &u, &m, nk, by_last, 2);

            /* a key that is a string in one dataset and a number in the other */
            if (nm && nu)
            {
                cj_plan p;
                stata_data mk = {m.n, (size_t)m.ncols, m.vars, NULL};
                stata_variable held = m.vars[0];
                m.vars[0].type =
                    held.type == STATA_TYPE_STRING ? STATA_TYPE_DOUBLE : STATA_TYPE_STRING;
                assert(cj_plan_build(&p, &ix, &mk, m.n, 0, NULL) == 106);
                m.vars[0] = held;
            }
            cj_index_free(&ix);
            fixture_free(&u);
            fixture_free(&m);
        }
    }
}

static int double_cmp(const void *a, const void *b)
{
    const double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* Enough rows for every parallel kernel; the oracle is a sorted copy */
static void large_joins(void)
{
    for (int hashed = 0; hashed <= 1; hashed++)
    {
        const size_t nu = 70000, nm = 90000;
        const int types[MAXCOLS] = {0};
        fixture u, m;
        cj_table table;
        cj_index ix;
        cj_plan p;
        fixture_init(&u, nu, 1, types);
        fixture_init(&m, nm, 1, types);
        for (size_t i = 0; i < nu; i++)
            u.num[0][i] = hashed ? (rng() % 40000) * 1000003.0 + .5 : (double)(rng() % 40000);
        for (size_t i = 0; i < nm; i++)
            m.num[0][i] = hashed ? (rng() % 50000) * 1000003.0 + .5 : (double)(rng() % 50000);
        /* runs of one key, as in a sorted master */
        for (size_t i = 1; i < nm; i++)
            if (i % 7)
                m.num[0][i] = m.num[0][i - 1];
        fixture_bind(&m);
        fixture_table(&u, &table);
        assert(cj_index_build(&ix, &table, 1, nm) == STATA_OK);
        check_index(&ix, &u, 1, !hashed);
        stata_data mk = {m.n, 1, m.vars, NULL};
        assert(cj_plan_build(&p, &ix, &mk, nm, CJ_KEEP_MASTER | CJ_KEEP_USING, NULL) == 0);

        double *sorted = malloc(nu * sizeof(double));
        assert(sorted);
        memcpy(sorted, u.num[0], nu * sizeof(double));
        qsort(sorted, nu, sizeof(double), double_cmp);
        size_t both = 0, lone = 0;
        for (size_t i = 0; i < nm; i++)
        {
            const double *found = bsearch(&m.num[0][i], sorted, nu, sizeof(double), double_cmp);
            const uint32_t g = p.gid[i];
            assert((found != NULL) == (g != CJ_NONE));
            if (g == CJ_NONE)
            {
                lone++;
                assert(p.off[i + 1] - p.off[i] == 1);
                continue;
            }
            assert(u.num[0][ix.keys.sort_order[ix.gstart[g]]] == m.num[0][i]);
            assert(p.off[i + 1] - p.off[i] == cj_group_count(&ix, g));
            both += cj_group_count(&ix, g);
        }
        assert(p.n_both == both && p.n_master == lone && p.nout == both + lone + p.n_using);
        free(sorted);
        cj_plan_free(&p);
        cj_index_free(&ix);
        fixture_free(&u);
        fixture_free(&m);
    }
}

/* Master keys one rounding step from an integer, against using keys whose
 * direct table starts far below zero: v - base rounds them onto an integer
 * (1e-15 + 127 == 127), so the lookup must test v itself. */
static void near_integer_joins(void)
{
    static const double bases[] = {-1, -127, -32767, -1e6};
    const double eps = 2.220446049250313e-16; /* epsdouble() */
    const double near[] = {1e-15, -1e-15, 2 - eps, 2 - 2 * eps, 2 + 2 * eps, 3 - eps,
                           3 + 4 * eps, 5 + 1e-12, 5 - 1e-12, 5 + 1e-11, 0.5, -0.0,
                           0, 2, 5, 6, -1, -127, -32767, -1e6, -1e6 + 1e-9, 1e-300};
    const size_t nnear = sizeof(near) / sizeof(near[0]);
    for (size_t b = 0; b < sizeof(bases) / sizeof(bases[0]); b++)
    {
        const double ukeys[] = {bases[b], 0, 2, 5, missing_value(0)};
        const size_t nu = sizeof(ukeys) / sizeof(ukeys[0]);
        const int types[MAXCOLS] = {0};
        for (size_t nm = 1; nm <= nnear; nm += nnear - 1)
        {
            fixture u, m;
            cj_table table;
            cj_index ix;
            fixture_init(&u, nu, 1, types);
            fixture_init(&m, nm, 1, types);
            memcpy(u.num[0], ukeys, sizeof(ukeys));
            /* one row (each value alone) and all of them together */
            for (size_t start = 0; start + nm <= nnear; start += nm)
            {
                memcpy(m.num[0], near + start, nm * sizeof(double));
                fixture_bind(&m);
                fixture_table(&u, &table);
                assert(cj_index_build(&ix, &table, 1, nm) == STATA_OK);
                check_index(&ix, &u, 1, 1);
                for (int unmatched = 0; unmatched < 4; unmatched++)
                    check_plan(&ix, &u, &m, 1, unmatched);
                cj_index_free(&ix);
            }
            fixture_free(&u);
            fixture_free(&m);
        }
    }
}

int main(void)
{
    mock.missval = missing_value(0);
    for (threads = 1; threads <= 4; threads += 3)
    {
        rng_state = 88172645463325252ULL;
        small_joins();
        large_joins();
        near_integer_joins();
#ifndef _OPENMP
        break;
#endif
    }
    puts("cjoinby index and plan agree with a brute-force join");
    return 0;
}

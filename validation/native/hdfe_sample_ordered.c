/* Sample construction that shares its steps' work (hdfe_load_factors_ordered,
 * hdfe_drop_singletons, hdfe_build_factors_ordered, hdfe_count_factors)
 * against an independent reference and against the separate steps
 * (hdfe_load_factors, ctools_remove_singletons, hdfe_build_factors): the kept
 * rows, the levels in order of first appearance among them, the rows per
 * level and the weights per level summed in row order must match bit for
 * bit. The thresholds are small, so that a few hundred rows reach the row
 * blocks, the ranges of codes and levels, the counting pass of the singleton
 * step, its fallback, and both parallel counts. */
#undef NDEBUG
#ifndef CTOOLS_HDFE_BLOCK_ROWS
#define CTOOLS_HDFE_BLOCK_ROWS 16
#endif
#ifndef HDFE_SINGLETON_PROBE_ROWS
#define HDFE_SINGLETON_PROBE_ROWS 8
#endif
#ifndef HDFE_LISTED_MIN_DROPPED
#define HDFE_LISTED_MIN_DROPPED 8
#endif
#ifndef HDFE_ORDER_RANGE_LEVELS
#define HDFE_ORDER_RANGE_LEVELS 4
#endif
#ifndef REMAP_SCATTER_CODES
#define REMAP_SCATTER_CODES 4
#endif
#ifndef HDFE_COUNT_FAR_SLOTS
#define HDFE_COUNT_FAR_SLOTS 8
#endif
#ifndef HDFE_COUNT_BLOCK_ENTRIES
#define HDFE_COUNT_BLOCK_ENTRIES 64
#endif
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "stplugin.h"

static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static ST_boolean mock_missing(ST_double x) { return x >= 8.988465674311579e307; }
static ST_retcode mock_out(char *s)
{
    (void)s;
    return 0;
}
double ctools_timer_seconds(void) { return 0; }

#include "ctools_hdfe_utils.c"
#include "hdfe/hdfe_remap.c"
#include "hdfe/hdfe_sample.c"
#include "cppmlhdfe/cppmlhdfe_separation.c"

static uint64_t rs = 20261126;
static uint64_t rnd(void)
{
    uint64_t z = (rs += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static uint64_t below(uint64_t n) { return n ? rnd() % n : 0; }

static void threads(int n)
{
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(n);
#else
    (void)n;
#endif
}

#define MAXG 3

typedef struct {
    ST_int N;                 /* kept rows */
    ST_int *keep;             /* flags of the loaded rows */
    ST_int L[MAXG];
    ST_int *lev[MAXG];        /* levels of the kept rows */
    double *cnt[MAXG], *wsum[MAXG];
} sample;

static void sample_free(sample *s, int G)
{
    free(s->keep);
    for (int g = 0; g < G; g++) {
        free(s->lev[g]);
        free(s->cnt[g]);
        free(s->wsum[g]);
    }
    memset(s, 0, sizeof(*s));
}

/* The reference: equal values by comparison, singletons by rounds of
 * recounting, levels numbered at their first kept row, sums in row order */
static void reference(sample *s, int G, const double *const *v, ST_int N, const double *fw,
                      const double *w)
{
    ST_int *group[MAXG];
    memset(s, 0, sizeof(*s));
    assert((s->keep = malloc((size_t)(N ? N : 1) * sizeof(ST_int))));
    for (int g = 0; g < G; g++) {
        /* group[g][i] = first row with row i's value */
        assert((group[g] = malloc((size_t)(N ? N : 1) * sizeof(ST_int))));
        for (ST_int i = 0; i < N; i++) {
            ST_int j = 0;
            while (v[g][j] != v[g][i]) j++;
            group[g][i] = j;
        }
    }
    for (ST_int i = 0; i < N; i++) s->keep[i] = 1;
    ST_int *cnt = calloc((size_t)(N ? N : 1), sizeof(ST_int));
    assert(cnt);
    for (;;) {
        ST_int round = 0;
        for (int g = 0; g < G; g++) {
            memset(cnt, 0, (size_t)N * sizeof(ST_int));
            for (ST_int i = 0; i < N; i++) cnt[group[g][i]] += s->keep[i];
            for (ST_int i = 0; i < N; i++)
                if (s->keep[i] && cnt[group[g][i]] == 1 && (!fw || fw[i] == 1.0)) {
                    s->keep[i] = 0;
                    round++;
                }
        }
        if (!round) break;
    }
    for (ST_int i = 0; i < N; i++) s->N += s->keep[i];
    for (int g = 0; g < G; g++) {
        ST_int next = 0, idx = 0;
        assert((s->lev[g] = malloc((size_t)(s->N ? s->N : 1) * sizeof(ST_int))));
        memset(cnt, 0, (size_t)N * sizeof(ST_int)); /* new number of each group */
        for (ST_int i = 0; i < N; i++) {
            if (!s->keep[i]) continue;
            if (!cnt[group[g][i]]) cnt[group[g][i]] = ++next;
            s->lev[g][idx++] = cnt[group[g][i]];
        }
        s->L[g] = next;
        assert((s->cnt[g] = calloc((size_t)(next ? next : 1), sizeof(double))));
        assert((s->wsum[g] = calloc((size_t)(next ? next : 1), sizeof(double))));
        idx = 0;
        for (ST_int i = 0; i < N; i++) {
            if (!s->keep[i]) continue;
            s->cnt[g][s->lev[g][idx] - 1] += 1.0;
            if (w) s->wsum[g][s->lev[g][idx] - 1] += w[i];
            idx++;
        }
        free(group[g]);
    }
    free(cnt);
}

static void from_state(sample *s, int G, const HDFE_State *S, const ST_int *keep, ST_int N,
                       ST_int dropped, int weighted)
{
    memset(s, 0, sizeof(*s));
    s->N = N - dropped;
    assert(S->N == s->N && S->G == G);
    assert((s->keep = malloc((size_t)(N ? N : 1) * sizeof(ST_int))));
    for (ST_int i = 0; i < N; i++) s->keep[i] = keep ? keep[i] : 1;
    for (int g = 0; g < G; g++) {
        const FE_Factor *f = &S->factors[g];
        const size_t L = (size_t)f->num_levels;
        s->L[g] = f->num_levels;
        assert(f->max_level == f->num_levels && f->has_intercept == 1);
        assert((s->lev[g] = malloc((size_t)(s->N ? s->N : 1) * sizeof(ST_int))));
        memcpy(s->lev[g], f->levels, (size_t)s->N * sizeof(ST_int));
        assert((s->cnt[g] = calloc(L ? L : 1, sizeof(double))));
        assert((s->wsum[g] = calloc(L ? L : 1, sizeof(double))));
        memcpy(s->cnt[g], f->counts, L * sizeof(double));
        assert(!weighted == !f->weighted_counts);
        if (weighted) memcpy(s->wsum[g], f->weighted_counts, L * sizeof(double));
    }
}

/* 0 on success, 1 when no row is left (nothing is built then) */
static int separate_steps(sample *s, int G, const double *const *v, ST_int N, const double *fw,
                          const double *w)
{
    HDFE_LoadedFactor fac[MAXG];
    HDFE_State S;
    ST_int *lv[MAXG], lmax[MAXG], *mask;
    memset(fac, 0, sizeof(fac));
    memset(&S, 0, sizeof(S));
    assert((mask = malloc((size_t)N * sizeof(ST_int))));
    assert(hdfe_load_factors(fac, G, v, N) == 0);
    for (int g = 0; g < G; g++) lv[g] = fac[g].levels, lmax[g] = fac[g].num_levels;
    ST_int dropped = ctools_remove_singletons(lv, G, N, mask, fw, 0);
    assert(dropped >= 0);
    int none = dropped == N;
    if (!none) {
        assert(hdfe_build_factors(&S, G, lv, lmax, N, mask, N - dropped, w, 0) == 0);
        from_state(s, G, &S, mask, N, dropped, w != NULL);
        ctools_hdfe_state_cleanup(&S);
    }
    for (int g = 0; g < G; g++) free(fac[g].levels);
    free(mask);
    return none;
}

static int shared_listed, shared_dropped;

/* the counts either inside the build or left to the caller, job by job */
static int shared_steps(sample *s, int G, const double *const *v, ST_int N, const double *fw,
                        const double *w, int jobs)
{
    HDFE_LoadedFactor fac[MAXG];
    HDFE_FactorOrder ord[MAXG];
    HDFE_Singletons sing;
    HDFE_State S;
    ST_int *lv[MAXG], lmax[MAXG], pending = 0;
    memset(fac, 0, sizeof(fac));
    memset(ord, 0, sizeof(ord));
    memset(&sing, 0, sizeof(sing));
    memset(&S, 0, sizeof(S));
    assert(hdfe_load_factors_ordered(fac, ord, G, v, N) == 0);
    for (int g = 0; g < G; g++) {
        lv[g] = fac[g].levels, lmax[g] = fac[g].num_levels;
        /* small samples always come out in order of first appearance */
        assert(ord[g].ordered);
        if (ord[g].rows.num_single >= 0) {
            /* the rows alone in their level, ascending */
            ST_int *c = calloc((size_t)lmax[g] + 1, sizeof(ST_int)), k = 0;
            assert(c);
            for (ST_int i = 0; i < N; i++) c[lv[g][i]]++;
            for (ST_int i = 0; i < N; i++)
                if (c[lv[g][i]] == 1) assert(k < ord[g].rows.num_single && ord[g].rows.single_rows[k++] == i);
            assert(k == ord[g].rows.num_single);
            /* and every level's first row */
            memset(c, 0, ((size_t)lmax[g] + 1) * sizeof(ST_int));
            for (ST_int i = 0; i < N; i++)
                if (!c[lv[g][i]]++) assert(ord[g].rows.first_row[lv[g][i] - 1] == i);
            free(c);
        }
    }
    assert(hdfe_drop_singletons(&sing, lv, lmax, ord, G, N, fw, 0) == 0);
    assert(!sing.keep == !sing.num_dropped);
    shared_listed = sing.listed;
    shared_dropped = sing.num_dropped;
    int none = sing.num_dropped == N;
    if (!none) {
        assert(hdfe_build_factors_ordered(&S, G, fac, ord, &sing, N, w, 0, jobs ? &pending : NULL) == 0);
        if (pending) {
            if (jobs == 1) {
                hdfe_count_factors(&S, ord, N, sing.keep, w);
            } else {
                /* as a caller with helper threads would: serial passes
                 * first, then those that take the row blocks */
                for (int g = G - 1; g >= 0; g--) {
                    hdfe_count_weights(&S, g, N, sing.keep, w);
                    if (!hdfe_count_rows_parallel(&S, g, ord[g].rows.local)) hdfe_count_rows_serial(&S, g);
                }
                for (int g = 0; g < G; g++)
                    if (hdfe_count_rows_parallel(&S, g, ord[g].rows.local)) hdfe_count_rows(&S, g, ord[g].rows.local);
            }
        }
        from_state(s, G, &S, sing.keep, N, sing.num_dropped, w != NULL);
        ctools_hdfe_state_cleanup(&S);
    }
    for (int g = 0; g < G; g++) free(fac[g].levels);
    hdfe_factor_order_free(ord, G);
    hdfe_singletons_free(&sing);
    return none;
}

static void same(const sample *a, const sample *b, int G, ST_int N, int weighted);

static int keep_listed;

/* Keep flags made elsewhere (any rows, not only singletons): the factors
 * hdfe_build_factors() makes of them against those made from the list of
 * dropped rows (hdfe_singletons_from_keep) */
static void keep_steps(int G, const double *const *v, ST_int N, const double *w, ST_int *keep,
                       ST_int dropped)
{
    HDFE_LoadedFactor fac[MAXG];
    HDFE_FactorOrder ord[MAXG];
    HDFE_Singletons sing;
    HDFE_State S, T;
    ST_int *lv[MAXG], lmax[MAXG];
    sample a, b;
    memset(fac, 0, sizeof(fac));
    memset(ord, 0, sizeof(ord));
    memset(&S, 0, sizeof(S));
    memset(&T, 0, sizeof(T));
    assert(hdfe_load_factors_ordered(fac, ord, G, v, N) == 0);
    for (int g = 0; g < G; g++) lv[g] = fac[g].levels, lmax[g] = fac[g].num_levels;
    assert(hdfe_build_factors(&S, G, lv, lmax, N, keep, N - dropped, w, 0) == 0);
    from_state(&a, G, &S, keep, N, dropped, w != NULL);
    assert(hdfe_singletons_from_keep(&sing, lv, lmax, G, N, keep, dropped) == 0);
    assert(sing.keep == keep && sing.num_dropped == dropped);
    keep_listed = sing.listed;
    if (sing.listed && dropped) {
        ST_int k = 0;
        for (ST_int i = 0; i < N; i++)
            if (!keep[i]) assert(sing.dropped[k++] == i);
        assert(k == dropped);
        for (int g = 0; g < G; g++)
            for (ST_int a2 = 0; a2 < sing.affected[g].n; a2++) {
                ST_int left = 0;
                for (ST_int i = 0; i < N; i++) left += keep[i] && lv[g][i] == sing.affected[g].level[a2];
                assert(!left == !sing.affected[g].left[a2]);
            }
    }
    assert(hdfe_build_factors_ordered(&T, G, fac, ord, &sing, N, w, 0, NULL) == 0);
    from_state(&b, G, &T, keep, N, dropped, w != NULL);
    same(&a, &b, G, N, w != NULL);
    sample_free(&a, G);
    sample_free(&b, G);
    ctools_hdfe_state_cleanup(&S);
    ctools_hdfe_state_cleanup(&T);
    for (int g = 0; g < G; g++) free(fac[g].levels);
    hdfe_factor_order_free(ord, G);
    sing.keep = NULL; /* the caller's */
    hdfe_singletons_free(&sing);
}

/* The PPML selection keeps every row exactly when the one-pass test says so */
static void whole_steps(int G, const double *const *v, ST_int N)
{
    HDFE_LoadedFactor fac[MAXG];
    HDFE_FactorOrder ord[MAXG];
    ST_int *lv[MAXG], lmax[MAXG], single[MAXG], slope_only[MAXG], known = 1, ns = 0, nz = 0;
    memset(fac, 0, sizeof(fac));
    memset(ord, 0, sizeof(ord));
    assert(hdfe_load_factors_ordered(fac, ord, G, v, N) == 0);
    for (int g = 0; g < G; g++) {
        lv[g] = fac[g].levels, lmax[g] = fac[g].num_levels;
        single[g] = ord[g].rows.num_single;
        known &= single[g] >= 0;
        slope_only[g] = below(4) == 0;
    }
    double *y = malloc((size_t)N * sizeof(double));
    ST_int *mask = malloc((size_t)N * sizeof(ST_int));
    assert(y && mask);
    const int zeros = (int)below(4); /* 0: none */
    for (ST_int i = 0; i < N; i++) {
        y[i] = (zeros && below(zeros == 1 ? 2 : 6) == 0) ? 0.0 : 1.0 + (double)below(5);
        if (zeros == 3 && v[0][i] == v[0][0]) y[i] = 0.0; /* a level with no positive outcome */
        mask[i] = 1;
    }
    const ST_int keep_singletons = below(4) == 0, check = below(4) != 0;
    const ST_int *so = below(2) ? slope_only : NULL;
    const int whole = ppml_fe_sample_is_whole(y, lv, lmax, single, G, N, so, keep_singletons, check);
    assert(ppml_select_fe_sample_ex(y, lv, lmax, G, N, mask, &ns, &nz, so, keep_singletons, check) == 0);
    if (whole) assert(ns + nz == 0);
    if (known) assert(whole == (ns + nz == 0));
    free(y);
    free(mask);
    for (int g = 0; g < G; g++) free(fac[g].levels);
    hdfe_factor_order_free(ord, G);
}

static void same(const sample *a, const sample *b, int G, ST_int N, int weighted)
{
    assert(a->N == b->N);
    assert(!memcmp(a->keep, b->keep, (size_t)N * sizeof(ST_int)));
    for (int g = 0; g < G; g++) {
        assert(a->L[g] == b->L[g]);
        assert(!memcmp(a->lev[g], b->lev[g], (size_t)a->N * sizeof(ST_int)));
        assert(!memcmp(a->cnt[g], b->cnt[g], (size_t)a->L[g] * sizeof(double)));
        if (weighted) assert(!memcmp(a->wsum[g], b->wsum[g], (size_t)a->L[g] * sizeof(double)));
    }
}

/* One factor's values: kind picks the shape */
static void fill(double *x, ST_int N, int kind)
{
    const uint64_t L = 1 + below(kind == 0 ? 3 : kind == 1 ? (uint64_t)N : kind == 2 ? (uint64_t)N * 4
                                                                                      : (uint64_t)N / 2 + 1);
    const double base = below(3) == 0 ? -50 : below(2) ? 1 : 1e9;
    const ST_int run = 1 + (ST_int)below(6);
    for (ST_int i = 0; i < N; i++) {
        if (kind == 3) x[i] = (double)(i / run);                          /* sorted runs */
        else if (kind == 4) x[i] = (double)below(L) * 0.5;                /* not integers */
        else if (kind == 5) x[i] = (double)(below(L) * 1000003ull) * 1e3; /* wide range */
        else if (kind == 6) x[i] = (double)i;                             /* all alone */
        else if (kind == 7) x[i] = (double)((i + 1) / 2);                 /* a chain with kind 8 */
        else if (kind == 8) x[i] = (double)(i / 2);
        else if (kind == 9) x[i] = below(5) == 0 ? -7.0 : (double)(i / run); /* sorted, one level everywhere */
        else if (kind == 10) x[i] = below(2) ? 0.0 : -0.0;                /* signed zeros: one level */
        else x[i] = base + (double)below(L);
    }
}

int main(void)
{
    mock.ismissing = mock_missing;
    mock.spoutsml = mock_out;
    int listed = 0, with_drops = 0, listed_drops = 0, fallback = 0, kept_listed = 0, kept_renumbered = 0;
    for (int c = 0; c < 2500; c++) {
        const ST_int N = 1 + (ST_int)below(below(4) == 0 ? 40 : 1200);
        const int G = 1 + (int)below(MAXG);
        double *v[MAXG], *fw = NULL, *w = NULL;
        for (int g = 0; g < MAXG; g++) assert((v[g] = malloc((size_t)N * sizeof(double))));
        if (below(8) == 0 && G >= 2) {
            fill(v[0], N, 7);
            fill(v[1], N, 8);
            for (int g = 2; g < G; g++) fill(v[g], N, 0);
        } else {
            for (int g = 0; g < G; g++) fill(v[g], N, (int)below(11));
        }
        if (below(4) == 0) {
            assert((fw = malloc((size_t)N * sizeof(double))));
            for (ST_int i = 0; i < N; i++) fw[i] = 1 + (double)below(3);
        }
        if (below(3) == 0) {
            assert((w = malloc((size_t)N * sizeof(double))));
            for (ST_int i = 0; i < N; i++) w[i] = 0.25 + (double)below(1000) / 997.0;
        }
        const double *const vals[MAXG] = {v[0], v[1], v[2]};
        sample ref, old, new1;
        reference(&ref, G, vals, N, fw, w);
        for (int t = 0; t < 3; t++) {
            static const int nthreads[3] = {1, 3, 8};
            threads(nthreads[t]);
            if (separate_steps(&old, G, vals, N, fw, w)) {
                assert(ref.N == 0);
            } else {
                same(&ref, &old, G, N, w != NULL);
                sample_free(&old, G);
            }
            for (int jobs = 0; jobs < 3; jobs++) {
                if (shared_steps(&new1, G, vals, N, fw, w, jobs)) {
                    assert(ref.N == 0);
                } else {
                    same(&ref, &new1, G, N, w != NULL);
                    sample_free(&new1, G);
                }
                assert(shared_dropped == N - ref.N);
            }
            if (t == 2) {
                listed += shared_listed;
                if (shared_dropped) with_drops++, listed_drops += shared_listed;
                if (shared_dropped && !shared_listed) fallback++;
            }
            /* any rows dropped: none, one, a few, many */
            {
                ST_int *keep = malloc((size_t)N * sizeof(ST_int)), dropped = 0;
                const uint64_t rate = below(4) == 0 ? 2 : below(2) ? 40 : 400;
                const int mode = (int)below(5);
                assert(keep);
                for (ST_int i = 0; i < N; i++) {
                    keep[i] = mode == 0 ? 1 : mode == 1 ? i != N / 2 : below(rate) != 0;
                    dropped += !keep[i];
                }
                if (dropped < N) {
                    keep_steps(G, vals, N, w, keep, dropped);
                    if (dropped && t == 2) keep_listed ? kept_listed++ : kept_renumbered++;
                }
                free(keep);
            }
            whole_steps(G, vals, N);
        }
        sample_free(&ref, G);
        for (int g = 0; g < MAXG; g++) free(v[g]);
        free(fw);
        free(w);
    }
#ifdef _OPENMP
    /* both ways of removing singletons were used, and both ways of building
     * the factors from keep flags */
    assert(listed_drops > 100 && fallback > 100);
    assert(kept_listed > 100 && kept_renumbered > 100);
#endif
    printf("hdfe_sample_ordered: listed %d, with drops %d (listed %d, every level counted %d); "
           "keep flags listed %d, renumbered %d\n",
           listed, with_drops, listed_drops, fallback, kept_listed, kept_renumbered);
    return 0;
}

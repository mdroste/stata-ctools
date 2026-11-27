/* Parallel sample construction (FE remap, singleton removal, PPML sample
 * selection, factor renumbering, connected components, nesting) against the
 * serial paths and independent references. */
#undef NDEBUG
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

/* Singletons by rounds until nothing changes (reghdfe's loop): the same
 * fixed point as any order of removals */
static ST_int ref_singletons(ST_int *const *lev, const ST_int *lmax, int G, ST_int N,
                             const double *fw, ST_int *mask)
{
    ST_int *cnt[3] = {0}, dropped = 0;
    for (int g = 0; g < G; g++) assert((cnt[g] = calloc((size_t)lmax[g] + 1, sizeof(ST_int))));
    for (ST_int i = 0; i < N; i++) mask[i] = 1;
    for (;;) {
        for (int g = 0; g < G; g++) memset(cnt[g], 0, ((size_t)lmax[g] + 1) * sizeof(ST_int));
        for (ST_int i = 0; i < N; i++)
            for (int g = 0; g < G; g++) cnt[g][lev[g][i]] += mask[i];
        ST_int round = 0;
        for (ST_int i = 0; i < N; i++) {
            if (!mask[i] || (fw && fw[i] != 1.0)) continue;
            for (int g = 0; g < G; g++)
                if (cnt[g][lev[g][i]] == 1) {
                    mask[i] = 0;
                    round++;
                    break;
                }
        }
        if (!round) break;
        dropped += round;
    }
    for (int g = 0; g < G; g++) free(cnt[g]);
    return dropped;
}

/* PPML sample selection as before its counts were updated incrementally:
 * every round recounts the kept rows and positive outcomes */
static ST_int ref_ppml_select(const double *y, ST_int *const *lev, const ST_int *nl, int G, ST_int N,
                              ST_int *mask, ST_int *ns, ST_int *nz, const ST_int *slope_only,
                              ST_int keep_singletons, ST_int check_separation)
{
    ST_int *cnt[3] = {0};
    unsigned char *pos[3] = {0};
    *ns = *nz = 0;
    for (int g = 0; g < G; g++) {
        assert((cnt[g] = calloc((size_t)nl[g], sizeof(ST_int))));
        assert((pos[g] = calloc((size_t)nl[g], 1)));
    }
    for (;;) {
        for (int g = 0; g < G; g++) {
            memset(cnt[g], 0, (size_t)nl[g] * sizeof(ST_int));
            memset(pos[g], 0, (size_t)nl[g]);
            for (ST_int i = 0; i < N; i++) {
                if (!mask[i]) continue;
                cnt[g][lev[g][i] - 1]++;
                if (y[i] > 0) pos[g][lev[g][i] - 1] = 1;
            }
        }
        ST_int removed = 0;
        for (ST_int i = 0; i < N; i++) {
            if (!mask[i]) continue;
            int singleton = 0, separated = 0;
            for (int g = 0; g < G; g++) {
                singleton |= !keep_singletons && cnt[g][lev[g][i] - 1] == 1;
                separated |= check_separation && !(slope_only && slope_only[g]) && !pos[g][lev[g][i] - 1];
            }
            if (singleton || separated) {
                mask[i] = 0;
                removed++;
                if (singleton) (*ns)++;
                else (*nz)++;
            }
        }
        if (!removed) break;
    }
    for (int g = 0; g < G; g++) {
        free(cnt[g]);
        free(pos[g]);
    }
    return 0;
}

static uint64_t rs = 20261125;
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

/* The same partition of the rows: a bijection between the two numberings */
static int same_partition(const ST_int *a, const ST_int *b, ST_int N, ST_int La, ST_int Lb)
{
    if (La != Lb) return 0;
    ST_int *ab = calloc((size_t)La + 1, sizeof(ST_int)), *ba = calloc((size_t)Lb + 1, sizeof(ST_int));
    assert(ab && ba);
    int ok = 1;
    for (ST_int i = 0; i < N && ok; i++) {
        if (a[i] < 1 || a[i] > La || b[i] < 1 || b[i] > Lb) ok = 0;
        else if (!ab[a[i]] && !ba[b[i]]) ab[a[i]] = b[i], ba[b[i]] = a[i];
        else if (ab[a[i]] != b[i] || ba[b[i]] != a[i]) ok = 0;
    }
    free(ab);
    free(ba);
    return ok;
}

/* Components holding an FE1 level, with the rank-based union-find */
static ST_int ref_components(const ST_int *x, const ST_int *y, ST_int N, ST_int L1, ST_int L2)
{
    ctools_UnionFind *uf = ctools_uf_create(L1 + L2);
    char *seen = calloc((size_t)L1 + L2 + 1, 1);
    assert(uf && seen);
    for (ST_int i = 0; i < N; i++) ctools_uf_union(uf, x[i] - 1, L1 + y[i] - 1);
    ST_int n = 0;
    for (ST_int i = 0; i < L1; i++) {
        ST_int r = ctools_uf_find(uf, i);
        if (!seen[r]) seen[r] = 1, n++;
    }
    free(seen);
    ctools_uf_destroy(uf);
    return n;
}

/* Codes of one absorbed variable; kind 11 places them in the components
 * that comp[] assigns, so that pairs of such factors are disconnected */
static void gen_factor(double *v, ST_int N, int kind, const ST_int *comp)
{
    ST_int K = 1 + (ST_int)below((uint64_t)N / 4 + 2);
    for (ST_int i = 0; i < N; i++) {
        switch (kind) {
        case 0: v[i] = 1 + (double)below(K); break;                           /* dense */
        case 1: v[i] = 3e9 + (double)below(K); break;                         /* above 2^31 */
        case 2: v[i] = -(double)(K / 2) + (double)below(K); break;            /* negative */
        case 3: v[i] = below(3) == 0 ? -0.0 : (below(2) ? 0.0 : 1 + (double)below(5)); break;
        case 4: v[i] = 0.5 + (double)below(K); break;                         /* non-integer */
        case 5: v[i] = (double)((below(K) * 0x9E3779B97F4A7C15ULL) >> 24); break; /* sparse */
        case 6: v[i] = below(2) ? 0 : 0x1p52; break;                          /* huge range */
        case 7: v[i] = 7; break;                                              /* one level */
        case 8: v[i] = (double)i * 3 - 5; break;                              /* unique */
        case 9: v[i] = -2147483648.0 + (double)below(K); break;               /* int32 minimum */
        case 10: v[i] = 1990 + (double)below(20); break;                      /* few */
        default: v[i] = (double)comp[i] * 1000 + (double)below(7); break;
        }
    }
}

static void run_case(int cs, int T)
{
    static const ST_int sizes[] = {131072, 131073, 140000, 262145};
    const ST_int N = sizes[below(4)];
    const int G = 1 + (int)below(3);
    double *v[3] = {0};
    ST_int *comp = malloc((size_t)N * sizeof(ST_int));
    assert(comp);
    const ST_int ncomp = 1 + (ST_int)below(40);
    for (ST_int i = 0; i < N; i++) comp[i] = (ST_int)below((uint64_t)ncomp);
    for (int g = 0; g < G; g++) {
        v[g] = malloc((size_t)N * sizeof(double));
        assert(v[g]);
        gen_factor(v[g], N, below(3) == 0 ? 11 : (int)below(11), comp);
    }
    const ST_double *vals[3] = {v[0], v[1], v[2]};

    /* Load: the same partition serially and in parallel; remap_and_count
     * keeps first-appearance numbering, which the serial remap gives */
    HDFE_LoadedFactor f1[3], fT[3];
    memset(f1, 0, sizeof f1);
    memset(fT, 0, sizeof fT);
    threads(1);
    assert(hdfe_load_factors(f1, G, vals, N) == 0);
    threads(T);
    assert(hdfe_load_factors(fT, G, vals, N) == 0);
    for (int g = 0; g < G; g++) {
        assert(fT[g].counts == NULL);
        assert(same_partition(f1[g].levels, fT[g].levels, N, f1[g].num_levels, fT[g].num_levels));
    }
    {
        ST_int *lv = malloc((size_t)N * sizeof(ST_int)), L = 0, *counts = NULL;
        assert(lv);
        assert(remap_and_count(v[0], N, lv, &L, &counts, NULL, NULL) == 0);
        assert(L == f1[0].num_levels && !memcmp(lv, f1[0].levels, (size_t)N * sizeof(ST_int)));
        ST_int total = 0;
        for (ST_int l = 0; l < L; l++) total += counts[l];
        assert(total == N);
        free(counts);
        assert(remap_and_count(v[0], N, lv, &L, NULL, NULL, NULL) == 0 && L == f1[0].num_levels);
        free(lv);
    }

    /* Renumbering on all rows, after singleton removal, or after random drops */
    ST_int *keep = malloc((size_t)N * sizeof(ST_int)), Nk = N;
    double *w = malloc((size_t)N * sizeof(double));
    assert(keep && w);
    ST_int *lv[3], lmax[3];
    for (int g = 0; g < G; g++) {
        lv[g] = fT[g].levels;
        lmax[g] = fT[g].num_levels + (ST_int)below(3);
    }

    /* Singleton removal, without and with fweights of 1 and 2: serial,
     * parallel and by rounds */
    {
        ST_int *m1 = malloc((size_t)N * sizeof(ST_int)), *mT = malloc((size_t)N * sizeof(ST_int));
        ST_int *mr = malloc((size_t)N * sizeof(ST_int));
        double *fw = malloc((size_t)N * sizeof(double));
        assert(m1 && mT && mr && fw);
        for (ST_int i = 0; i < N; i++) fw[i] = below(4) ? 1.0 : 2.0;
        for (int weighted = 0; weighted < 2; weighted++) {
            const double *f = weighted ? fw : NULL;
            threads(1);
            const ST_int s1 = ctools_remove_singletons(lv, G, N, m1, f, 0);
            threads(T);
            const ST_int sT = ctools_remove_singletons(lv, G, N, mT, f, 0);
            const ST_int sr = ref_singletons(lv, lmax, G, N, f, mr);
            assert(s1 == sr && sT == sr);
            assert(!memcmp(m1, mr, (size_t)N * sizeof(ST_int)) && !memcmp(mT, mr, (size_t)N * sizeof(ST_int)));
        }
        free(m1); free(mT); free(mr); free(fw);
    }

    /* PPML selection (singletons and all-zero groups), against recounting
     * every round: outcomes with zeros, a fifth of comp groups all zero */
    {
        double *y = malloc((size_t)N * sizeof(double));
        ST_int *m1 = malloc((size_t)N * sizeof(ST_int)), *mT = malloc((size_t)N * sizeof(ST_int));
        ST_int *mr = malloc((size_t)N * sizeof(ST_int));
        assert(y && m1 && mT && mr);
        for (ST_int i = 0; i < N; i++) y[i] = (comp[i] % 5 == 0 || below(3) == 0) ? 0 : 1 + (double)below(4);
        const ST_int slope_only[3] = {0, (ST_int)below(2), 0};
        for (int variant = 0; variant < 3; variant++) {
            const ST_int ks = variant == 1, cs = variant != 2;
            ST_int n1, z1, nT, zT, nr, zr;
            for (ST_int i = 0; i < N; i++) m1[i] = mT[i] = mr[i] = below(50) != 0;
            const ST_int *nl = lmax;
            ref_ppml_select(y, lv, nl, G, N, mr, &nr, &zr, slope_only, ks, cs);
            threads(1);
            assert(ppml_select_fe_sample_ex(y, lv, nl, G, N, m1, &n1, &z1, slope_only, ks, cs) == 0);
            threads(T);
            assert(ppml_select_fe_sample_ex(y, lv, nl, G, N, mT, &nT, &zT, slope_only, ks, cs) == 0);
            assert(n1 == nr && z1 == zr && nT == nr && zT == zr);
            assert(!memcmp(m1, mr, (size_t)N * sizeof(ST_int)) && !memcmp(mT, mr, (size_t)N * sizeof(ST_int)));
        }
        free(y); free(m1); free(mT); free(mr);
    }
    const int keepmode = (int)below(3);
    if (keepmode == 1) {
        Nk = N - ctools_remove_singletons(lv, G, N, keep, NULL, 0);
    } else {
        Nk = 0;
        for (ST_int i = 0; i < N; i++) Nk += keep[i] = keepmode == 0 || below(10) != 0;
    }
    for (ST_int i = 0; i < N; i++) w[i] = 0.25 + (double)below(1u << 20) / 7.0;
    for (int variant = 0; variant < 3 && Nk > 0; variant++) {
        const ST_int *kp = keepmode == 0 && below(2) ? NULL : keep;
        HDFE_State S1, ST;
        memset(&S1, 0, sizeof S1);
        memset(&ST, 0, sizeof ST);
        threads(1);
        assert(hdfe_build_factors(&S1, G, lv, lmax, N, kp, Nk, variant == 1 ? w : NULL, variant == 2) == 0);
        threads(T);
        assert(hdfe_build_factors(&ST, G, lv, lmax, N, kp, Nk, variant == 1 ? w : NULL, variant == 2) == 0);
        for (int g = 0; g < G; g++) {
            const FE_Factor *a = &S1.factors[g], *b = &ST.factors[g];
            const ST_int L = a->num_levels;
            assert(L == b->num_levels && a->max_level == b->max_level);
            assert(!memcmp(a->levels, b->levels, (size_t)Nk * sizeof(ST_int)));
            assert(!memcmp(a->counts, b->counts, (size_t)L * sizeof(double)));
            assert(!a->weighted_counts == !b->weighted_counts);
            if (a->weighted_counts)
                assert(!memcmp(a->weighted_counts, b->weighted_counts, (size_t)L * sizeof(double)));
        }
        if (variant == 0) {
            for (int g = 0; g < G; g++)
                for (int h = 0; h < g; h++) {
                    const FE_Factor *x = &ST.factors[h], *y = &ST.factors[g];
                    const ST_int ref = ref_components(x->levels, y->levels, Nk, x->num_levels, y->num_levels);
                    threads(1);
                    assert(ctools_count_connected_components(x->levels, y->levels, Nk, x->num_levels,
                                                             y->num_levels) == ref);
                    threads(T);
                    assert(ctools_count_connected_components(x->levels, y->levels, Nk, x->num_levels,
                                                             y->num_levels) == ref);
                }
            /* nested, nested but for one row, unrelated */
            for (int mode = 0; mode < 3; mode++) {
                const FE_Factor *x = &ST.factors[0];
                ST_int *cl = malloc((size_t)Nk * sizeof(ST_int));
                assert(cl);
                for (ST_int i = 0; i < Nk; i++) cl[i] = mode == 2 ? (ST_int)below(5) : (x->levels[i] * 7) % 13;
                if (mode == 1) cl[below((uint64_t)Nk)] += 1;
                threads(1);
                const ST_int n1 = ctools_fe_nested_in_cluster(x->levels, x->num_levels, cl, Nk);
                threads(T);
                const ST_int nT = ctools_fe_nested_in_cluster(x->levels, x->num_levels, cl, Nk);
                assert(n1 == nT && (mode != 0 || n1 == 1));
                free(cl);
            }
        }
        ctools_hdfe_state_cleanup(&S1);
        ctools_hdfe_state_cleanup(&ST);
    }
    for (int g = 0; g < G; g++) {
        free(f1[g].levels);
        free(fT[g].levels);
        free(v[g]);
    }
    free(keep);
    free(w);
    free(comp);
    (void)cs;
}

int main(void)
{
    mock.ismissing = mock_missing;
    mock.spoutsml = mock_out;
    /* Inside a parallel region every path runs serially */
#ifdef _OPENMP
#pragma omp parallel num_threads(2)
    {
        assert(ctools_hdfe_row_blocks(1 << 20) == 1);
    }
    threads(4);
    assert(ctools_hdfe_row_blocks(1 << 20) == 4 && ctools_hdfe_row_blocks(100000) == 1);
#else
    assert(ctools_hdfe_row_blocks(1 << 20) == 1);
#endif
    static const int team[] = {4, 3, 7};
    for (int cs = 0; cs < 18; cs++) run_case(cs, team[cs % 3]);
    /* Components of two effects of few levels: every block joins them all
     * within its first rows, and the count stops there. Then the same joined
     * by the last row alone (nothing can stop early), and not joined. */
    {
        const ST_int N = 400000;
        ST_int *a = malloc((size_t)N * sizeof(ST_int)), *b = malloc((size_t)N * sizeof(ST_int));
        assert(a && b);
        for (int T = 1; T <= 4; T += 3) {
            threads(T);
            for (ST_int i = 0; i < N; i++) a[i] = 1 + i % 7, b[i] = 1 + (i / 7) % 5;
            assert(ref_components(a, b, N, 7, 5) == 1);
            assert(ctools_count_connected_components(a, b, N, 7, 5) == 1);
            for (ST_int i = 0; i < N; i++) a[i] = b[i] = 1 + (i >= N / 2);
            assert(ctools_count_connected_components(a, b, N, 2, 2) == 2);
            b[N - 1] = 1;
            assert(ref_components(a, b, N, 2, 2) == 1);
            assert(ctools_count_connected_components(a, b, N, 2, 2) == 1);
            /* a level without rows is a component of its own in neither
             * count, and keeps the links below the number that stops */
            assert(ctools_count_connected_components(a, b, N, 2, 3) == ref_components(a, b, N, 2, 3));
        }
        free(a);
        free(b);
    }
    puts("hdfe sample parallel: ok");
    return 0;
}

/* Bootstrap engine tests (src/hdfe/hdfe_boot.c, the binomial sampler in
 * src/ctools_rng.h). Select one CTOOLS_CASE_* per compilation. */

#if defined(CTOOLS_CASE_boot_engine)
#undef NDEBUG
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ctools_rng.h"
#include "hdfe/hdfe_boot.h"

#define CHECK(cond)                                                                            \
    do                                                                                         \
    {                                                                                          \
        if (!(cond))                                                                           \
        {                                                                                      \
            fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);                    \
            exit(1);                                                                           \
        }                                                                                      \
    } while (0)

/* Generated with StataNow 19: set seed; the 2U uniforms bsample takes (draws,
 * then tie-breaks) from Mata's runiform(); then, after the same seed,
 * bsample, weight() and each unit's weight in bsample's unit order. Designs:
 * 13 observations; 13 observations in strata of 7 and 6; 6 clusters;
 * 7 clusters in strata of 4 and 3 (seeds 101, 202, 303, 404). */
static const double obs_u[] = {.061053244277654239,.032920169748520522,.63396126686271803,.52764955969457394,.018586553417547891,.022901864477512413,.84096607935979784,.32294726756864212,.037794883659786671,.016016725797258213,.17792536653141589,.17379473030012116,.051987988885491165,.099095640166035781,.35795708085949352,.18151687900760949,.93275291209140032,.05825301767940172,.57034812327556994,.95133266780128156,.58093136765799069,.72909883608844039,.81399382351418892,.98744512558341879,.70432907848358528,.24459987474433154};
static const int obs_m[] = {0,0,1,1,0,0,1,1,0,7,2,0,0};
static const double strobs_u[] = {.85359254995840528,.62178324566980236,.73564495324344559,.97932994129175388,.78594145940650118,.8819989837873955,.65581769473577134,.25045660453630358,.80187077600096701,.66076623440643523,.061447069276230004,.37344214738153225,.18190492558333193,.76539872241096318,.47731374212440902,.22795852144324213,.65353415858105912,.43276056793602191,.78256467463261614,.40542058211817855,.80896933617745936,.29294603383915863,.40710992671213952,.66060940408087354,.91084610663446031,.096653298098169027};
static const int strobs_m[] = {3,2,0,0,0,2,0,2,1,1,1,1,0};
static const double clust_u[] = {.1710413512579626,.41931434979881088,.52574401336470367,.9610704594344962,.79451348600617455,.74974710184200644,.23280161144629175,.27494222987288308,.20467883116433228,.016575012247240228,.91911667819092047,.43423829994557128};
static const int clust_m[] = {1,1,1,1,2,0};
static const double strclust_u[] = {.48334573796315949,.30360157243852626,.58371127395577516,.64408935271340073,.61918193833083157,.03129329392997493,.16602012270601108,.64215038710828576,.25290021009782337,.9750615025814372,.93406378909587173,.94077862279841706,.66075916290133552,.29079357603664724};
static const int strclust_m[] = {2,0,2,0,1,2,0};

static void check_bsample(const double *u, const int *expected, ST_int U, ST_int S,
                          const ST_int *starts, int cluster_mode)
{
    hdfe_boot_design d = {U, S, starts, NULL, 0};
    ST_int m[16];
    CHECK(hdfe_boot_design_check(&d) == 0);
    CHECK(hdfe_boot_map_bsample(&d, u, u + U, cluster_mode, m) == 0);
    for (ST_int i = 0; i < U; i++)
        CHECK(m[i] == expected[i]);
}

static void test_bsample_mapping(void)
{
    static const ST_int strobs_starts[] = {0, 7, 13};
    static const ST_int strclust_starts[] = {0, 4, 7};
    check_bsample(obs_u, obs_m, 13, 1, NULL, 0);
    check_bsample(strobs_u, strobs_m, 13, 2, strobs_starts, 0);
    check_bsample(clust_u, clust_m, 6, 1, NULL, 1);
    check_bsample(strclust_u, strclust_m, 7, 2, strclust_starts, 1);
    /* uniforms must lie in [0, 1) */
    hdfe_boot_design d = {2, 1, NULL, NULL, 0};
    double bad[4] = {0.5, 1.0, 0.1, 0.2};
    ST_int m[2];
    CHECK(hdfe_boot_map_bsample(&d, bad, bad + 2, 0, m) == 198);
}

static void test_design_check(void)
{
    static const ST_int empty[] = {0, 3, 3};
    static const ST_int wrong_end[] = {0, 2, 4};
    hdfe_boot_design d = {3, 2, empty, NULL, 0};
    CHECK(hdfe_boot_design_check(&d) == 2001);
    d.stratum_start = wrong_end;
    CHECK(hdfe_boot_design_check(&d) == 198);
    double freq[3] = {1, 2.5, 0};
    hdfe_boot_design f = {3, 1, NULL, freq, 0};
    CHECK(hdfe_boot_design_check(&f) == 198); /* not an integer */
    freq[1] = 2;
    CHECK(hdfe_boot_design_check(&f) == 0);
    /* Frequency totals above INT_MAX are valid (draws(fast) draws them);
     * draws(stata) cannot (draw_units -1), and a unit stands for at most
     * INT_MAX / 2 observations */
    double big[3] = {1e9, 1e9, 1e9};
    hdfe_boot_design b = {3, 1, NULL, big, 77};
    CHECK(hdfe_boot_design_check(&b) == 0);
    CHECK(hdfe_boot_draw_units(&b) == -1 && hdfe_boot_draw_total(&b) == 3e9);
    ST_int mb[3];
    hdfe_boot_draw_fast(&b, 0, mb);
    CHECK(mb[0] > 0 && mb[1] > 0 && mb[2] > 0 &&
          (int64_t)mb[0] + (int64_t)mb[1] + (int64_t)mb[2] == (int64_t)3000000000LL);
    big[0] = 1073741824.0;
    CHECK(hdfe_boot_design_check(&b) == 198);
}

static void test_fast_draws(void)
{
    enum { U = 50, R = 4000 };
    static const ST_int starts[] = {0, 20, 35, 50};
    hdfe_boot_design d = {U, 3, starts, NULL, 0x1234567ull};
    CHECK(hdfe_boot_design_check(&d) == 0);
    ST_int a[U], b[U], c[U];
    hdfe_boot_draw_fast(&d, 5, a);
    hdfe_boot_draw_fast(&d, 9, c);
    hdfe_boot_draw_fast(&d, 5, b);
    CHECK(memcmp(a, b, sizeof(a)) == 0); /* pure function of (seed, r) */
    CHECK(memcmp(a, c, sizeof(a)) != 0);
    double mean[U] = {0};
    for (ST_int r = 0; r < R; r++)
    {
        hdfe_boot_draw_fast(&d, r, a);
        for (ST_int s = 0; s < 3; s++)
        {
            ST_int total = 0;
            for (ST_int u = starts[s]; u < starts[s + 1]; u++)
                total += a[u];
            CHECK(total == starts[s + 1] - starts[s]);
        }
        for (ST_int u = 0; u < U; u++)
            mean[u] += a[u] / (double)R;
    }
    /* each multiplicity has mean 1 and variance (n - 1)/n */
    for (ST_int u = 0; u < U; u++)
        CHECK(fabs(mean[u] - 1.0) < 5.0 / sqrt((double)R));

    /* Frequency weights: Multinomial(F, freq / F) within strata */
    double freq[U];
    for (ST_int u = 0; u < U; u++)
        freq[u] = (double)((u * 7) % 6); /* includes zeros */
    hdfe_boot_design f = {U, 3, starts, freq, 99};
    CHECK(hdfe_boot_design_check(&f) == 0);
    double fmean[U] = {0}, F[3] = {0};
    for (ST_int s = 0; s < 3; s++)
        for (ST_int u = starts[s]; u < starts[s + 1]; u++)
            F[s] += freq[u];
    for (ST_int r = 0; r < R; r++)
    {
        hdfe_boot_draw_fast(&f, r, a);
        for (ST_int s = 0; s < 3; s++)
        {
            ST_int total = 0;
            for (ST_int u = starts[s]; u < starts[s + 1]; u++)
            {
                if (freq[u] == 0)
                    CHECK(a[u] == 0);
                total += a[u];
            }
            CHECK(total == (ST_int)F[s]);
        }
        for (ST_int u = 0; u < U; u++)
            fmean[u] += a[u] / (double)R;
    }
    for (ST_int s = 0; s < 3; s++)
    {
        for (ST_int u = starts[s]; u < starts[s + 1]; u++)
        {
            double p = freq[u] / F[s], sd = sqrt(F[s] * p * (1 - p) / R);
            CHECK(fabs(fmean[u] - freq[u]) <= 5.0 * sd + 1e-12);
        }
    }
}

static void test_binomial(void)
{
    static const struct
    {
        int64_t n;
        double p;
    } cases[] = {{5, 0.3}, {40, 0.2}, {60, 0.9}, {1000, 0.5}, {100000, 0.01},
                 {10000000, 0.3}, {200, 0.04}, {3, 0.999}};
    ctools_rng rng;
    ctools_rng_seed(&rng, 2024);
    CHECK(ctools_rng_binomial(&rng, 0, 0.5) == 0);
    CHECK(ctools_rng_binomial(&rng, 10, 0.0) == 0);
    CHECK(ctools_rng_binomial(&rng, 10, 1.0) == 10);
    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++)
    {
        const int M = 40000;
        double n = (double)cases[c].n, p = cases[c].p, s = 0, s2 = 0;
        for (int i = 0; i < M; i++)
        {
            int64_t k = ctools_rng_binomial(&rng, cases[c].n, p);
            CHECK(k >= 0 && k <= cases[c].n);
            s += (double)k;
            s2 += (double)k * (double)k;
        }
        double mean = s / M, var = (s2 - s * mean) / (M - 1), v = n * p * (1 - p);
        CHECK(fabs(mean - n * p) < 5.0 * sqrt(v / M));
        /* Exact fourth central moment of Bin(n,p) and sampling variance of
         * the unbiased sample variance. Use the same five-SE bound as the
         * mean, rather than a fixed percentage unrelated to M or p. */
        double mu4 = 3 * v * v + v * (1 - 6 * p * (1 - p));
        double var_se = sqrt((mu4 - ((double)(M - 3) / (M - 1)) * v * v) / M);
        CHECK(fabs(var - v) < 5.0 * var_se);
    }
}

/* Toy model: one row per unit; b = (weighted mean, weighted mean of
 * squares); three FE levels (level of unit u is u % 3), each the weighted
 * mean of its drawn rows; d per unit the deviation from the replicate mean
 * when the unit's level was drawn. Replicates in which unit 0 is drawn
 * zero times report a collinear status when fail_unit0 is set. */
typedef struct
{
    ST_int U;
    const double *y;
    int fail_unit0;
} toy_model;

static size_t toy_bytes(const void *m)
{
    (void)m;
    return 64;
}
static void *toy_new(const void *m)
{
    (void)m;
    return malloc(64);
}
static void toy_free(void *ws)
{
    free(ws);
}
static ST_retcode toy_fit(const void *mv, void *ws, const ST_int *mult, hdfe_boot_rep *out)
{
    const toy_model *m = mv;
    (void)ws;
    double w = 0, s = 0, s2 = 0, lw[3] = {0}, ls[3] = {0};
    for (ST_int u = 0; u < m->U; u++)
    {
        w += mult[u];
        s += mult[u] * m->y[u];
        s2 += mult[u] * m->y[u] * m->y[u];
        lw[u % 3] += mult[u];
        ls[u % 3] += mult[u] * m->y[u];
    }
    out->b[0] = s / w;
    out->b[1] = s2 / w;
    /* fe and d are NULL when the run has no fixed-effect target */
    for (int l = 0; out->fe && l < 3; l++)
        out->fe[l] = lw[l] > 0 ? ls[l] / lw[l] : NAN;
    for (ST_int u = 0; out->d && u < m->U; u++)
        out->d[u] = lw[u % 3] > 0 ? m->y[u] - out->b[0] : NAN;
    out->iterations = 1;
    if (m->fail_unit0 && mult[0] == 0)
        out->status = HDFE_BOOT_COLLINEAR;
    return 0;
}
static const hdfe_boot_model_ops toy_ops = {toy_bytes, toy_new, toy_free, toy_fit};

static int cmp_dbl(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static void run_toy(const hdfe_boot_design *d, const toy_model *m, ST_int R, const ST_int *chunks,
                    ST_int threads, hdfe_boot_acc *acc)
{
    CHECK(hdfe_boot_acc_init(acc, 2, 3, m->U, R) == 0);
    ST_int first = 0;
    for (ST_int i = 0; first < R; i++)
    {
        ST_int n = chunks[i] < R - first ? chunks[i] : R - first;
        CHECK(hdfe_boot_run(d, NULL, &toy_ops, m, first, n, threads, (size_t)1 << 30, acc) == 0);
        first += n;
    }
}

static void test_engine(void)
{
    enum { U = 24, R = 300 };
    double y[U];
    for (ST_int u = 0; u < U; u++)
        y[u] = sin(1.7 * u) + 0.1 * u;
    toy_model m = {U, y, 0};
    hdfe_boot_design d = {U, 1, NULL, NULL, 77};

    /* Chunking and threads do not change anything */
    static const ST_int all[] = {R}, ones[] = {1, 1, 1, 7, 7, 1000}, sevens[] = {7, 7, 7, 7, 1000};
    hdfe_boot_acc a, b, c;
    run_toy(&d, &m, R, all, 1, &a);
    run_toy(&d, &m, R, ones, 4, &b);
    run_toy(&d, &m, R, sevens, 3, &c);
    CHECK(memcmp(a.b, b.b, sizeof(double) * 2 * R) == 0);
    CHECK(memcmp(a.b, c.b, sizeof(double) * 2 * R) == 0);
    CHECK(memcmp(a.fe_m2, c.fe_m2, sizeof(double) * 3) == 0);
    CHECK(memcmp(a.d_m2, b.d_m2, sizeof(double) * U) == 0);
    CHECK(a.iterations == R && a.status_count[HDFE_BOOT_OK] == R);

    /* Independent recomputation of every replicate and of the summary */
    double *rep = malloc(sizeof(double) * 2 * R), *x = malloc(sizeof(double) * R);
    ST_int mult[U];
    for (ST_int r = 0; r < R; r++)
    {
        hdfe_boot_draw_fast(&d, r, mult);
        double w = 0, s = 0, s2 = 0;
        for (ST_int u = 0; u < U; u++)
        {
            w += mult[u];
            s += mult[u] * y[u];
            s2 += mult[u] * y[u] * y[u];
        }
        rep[2 * r] = s / w;
        rep[2 * r + 1] = s2 / w;
        CHECK(rep[2 * r] == a.b[2 * r] && rep[2 * r + 1] == a.b[2 * r + 1]);
    }
    double bhat[2] = {0, 0};
    for (ST_int u = 0; u < U; u++)
    {
        bhat[0] += y[u] / U;
        bhat[1] += y[u] * y[u] / U;
    }
    hdfe_boot_summary s, smse;
    CHECK(hdfe_boot_summarize(&a, bhat, 0, 95, &s) == 0);
    CHECK(hdfe_boot_summarize(&a, bhat, 1, 95, &smse) == 0);
    CHECK(s.N_reps == R && s.N_misreps == 0);
    for (int j = 0; j < 2; j++)
    {
        double mean = 0;
        for (ST_int r = 0; r < R; r++)
            mean += rep[2 * r + j] / R;
        CHECK(fabs(s.b_bs[j] - mean) < 1e-14 * (1 + fabs(mean)));
        for (int k = 0; k < 2; k++)
        {
            double mk = 0, v = 0, vm = 0;
            for (ST_int r = 0; r < R; r++)
                mk += rep[2 * r + k] / R;
            for (ST_int r = 0; r < R; r++)
            {
                v += (rep[2 * r + j] - mean) * (rep[2 * r + k] - mk);
                vm += (rep[2 * r + j] - bhat[j]) * (rep[2 * r + k] - bhat[k]);
            }
            CHECK(fabs(s.V[j * 2 + k] - v / (R - 1)) < 1e-12 * fabs(v / (R - 1)));
            CHECK(fabs(smse.V[j * 2 + k] - vm / R) < 1e-12 * fabs(vm / R));
        }
        /* se is the standard deviation about the mean, also under mse */
        CHECK(fabs(s.se[j] - sqrt(s.V[j * 3])) < 1e-12 * s.se[j]);
        CHECK(smse.se[j] == s.se[j]);
        double z = hdfe_boot_invnormal(0.975);
        CHECK(fabs(s.ci_normal[j] - (bhat[j] - z * s.se[j])) < 1e-12);
        for (ST_int r = 0; r < R; r++)
            x[r] = rep[2 * r + j];
        qsort(x, R, sizeof(double), cmp_dbl);
        /* R = 300: P = 7.5 and 292.5, so x(8) and x(293) */
        CHECK(s.ci_percentile[j] == x[7] && s.ci_percentile[2 + j] == x[292]);
        ST_int below = 0;
        for (ST_int r = 0; r < R; r++)
            below += x[r] <= bhat[j];
        CHECK(fabs(s.z0[j] - hdfe_boot_invnormal((double)below / R)) < 1e-15);
    }

    /* Fixed effects: the SE over estimable replicates, from the replicates */
    double fe_se[3], d_se[U];
    hdfe_boot_fe_se(&a, fe_se, d_se);
    for (int l = 0; l < 3; l++)
    {
        double vals[R], mean = 0, ss = 0;
        int n = 0;
        for (ST_int r = 0; r < R; r++)
        {
            hdfe_boot_draw_fast(&d, r, mult);
            double lw = 0, ls = 0;
            for (ST_int u = l; u < U; u += 3)
            {
                lw += mult[u];
                ls += mult[u] * y[u];
            }
            if (lw > 0)
                vals[n++] = ls / lw;
        }
        for (int i = 0; i < n; i++)
            mean += vals[i] / n;
        for (int i = 0; i < n; i++)
            ss += (vals[i] - mean) * (vals[i] - mean);
        CHECK(a.fe_n[l] == n);
        CHECK(fabs(fe_se[l] - sqrt(ss / (n - 1))) < 1e-12 * fe_se[l]);
    }

    /* Supplied multiplicities give the same replicates */
    ST_int *supplied = malloc(sizeof(ST_int) * U * 10);
    for (ST_int r = 0; r < 10; r++)
        hdfe_boot_draw_fast(&d, r, supplied + r * U);
    hdfe_boot_acc e;
    CHECK(hdfe_boot_acc_init(&e, 2, 3, U, 10) == 0);
    CHECK(hdfe_boot_run(&d, supplied, &toy_ops, &m, 0, 10, 2, (size_t)1 << 30, &e) == 0);
    CHECK(memcmp(e.b, a.b, sizeof(double) * 20) == 0);

    /* Contract errors: out-of-order chunk, memory budget */
    CHECK(hdfe_boot_run(&d, NULL, &toy_ops, &m, 5, 1, 1, (size_t)1 << 30, &e) == 198);
    hdfe_boot_acc g;
    CHECK(hdfe_boot_acc_init(&g, 2, 3, U, 5) == 0);
    CHECK(hdfe_boot_run(&d, NULL, &toy_ops, &m, 0, 5, 1, 10, &g) == 909);
    CHECK(g.R == 0);

    /* Failed replicates are dropped whole (casewise) */
    toy_model fm = {U, y, 1};
    hdfe_boot_acc f;
    CHECK(hdfe_boot_acc_init(&f, 2, 3, U, R) == 0);
    CHECK(hdfe_boot_run(&d, NULL, &toy_ops, &fm, 0, R, 2, (size_t)1 << 30, &f) == 0);
    ST_int failed = f.status_count[HDFE_BOOT_COLLINEAR];
    CHECK(failed > 0 && failed < R);
    hdfe_boot_summary sf;
    CHECK(hdfe_boot_summarize(&f, bhat, 0, 95, &sf) == 0);
    CHECK(sf.N_misreps == failed && sf.N_reps == R - failed);
    double mean0 = 0;
    for (ST_int r = 0; r < R; r++)
    {
        hdfe_boot_draw_fast(&d, r, mult);
        if (mult[0] > 0)
            mean0 += rep[2 * r] / (R - failed);
    }
    CHECK(fabs(sf.b_bs[0] - mean0) < 1e-13);

    /* Fewer than two complete replicates */
    hdfe_boot_acc h;
    hdfe_boot_summary sh;
    CHECK(hdfe_boot_acc_init(&h, 2, 0, 0, 1) == 0);
    CHECK(hdfe_boot_run(&d, NULL, &toy_ops, &m, 0, 1, 1, (size_t)1 << 30, &h) == 0);
    CHECK(hdfe_boot_summarize(&h, bhat, 0, 95, &sh) == 2000);

    hdfe_boot_summary_free(&s);
    hdfe_boot_summary_free(&smse);
    hdfe_boot_summary_free(&sf);
    hdfe_boot_acc_free(&a);
    hdfe_boot_acc_free(&b);
    hdfe_boot_acc_free(&c);
    hdfe_boot_acc_free(&e);
    hdfe_boot_acc_free(&f);
    hdfe_boot_acc_free(&g);
    hdfe_boot_acc_free(&h);
    hdfe_boot_acc_free(&h); /* idempotent */
    free(rep);
    free(x);
    free(supplied);
}

static void test_percentiles_and_normal(void)
{
    double x[10];
    for (int i = 0; i < 10; i++)
        x[i] = i + 1;
    CHECK(hdfe_boot_pctile(x, 10, 25) == 3);    /* P = 2.5: x(3) */
    CHECK(hdfe_boot_pctile(x, 10, 50) == 5.5);  /* P = 5: (x(5) + x(6)) / 2 */
    CHECK(hdfe_boot_pctile(x, 10, 2.5) == 1);   /* P = 0.25: x(1) */
    CHECK(hdfe_boot_pctile(x, 10, 0) == 1);     /* (x(0) + x(1)) / 2, x(0) = x(1) */
    CHECK(hdfe_boot_pctile(x, 10, 100) == 10);  /* (x(10) + x(11)) / 2, x(11) = x(10) */
    CHECK(fabs(hdfe_boot_invnormal(0.975) - 1.959963984540054) < 1e-14);
    CHECK(fabs(hdfe_boot_invnormal(0.5)) < 1e-15);
    /* Above z = 3, Phi(z) is too close to 1 for the round trip to be posed */
    for (double z = -8; z <= 3; z += 0.37)
        CHECK(fabs(hdfe_boot_invnormal(hdfe_boot_normal(z)) - z) < 1e-9 * (1 + fabs(z)));
}

int main(void)
{
    test_bsample_mapping();
    test_design_check();
    test_fast_draws();
    test_binomial();
    test_engine();
    test_percentiles_and_normal();
    printf("boot_engine: all checks passed\n");
    return 0;
}
#endif

/* Wild-bootstrap engine tests (src/hdfe/hdfe_wild.c): the five-number form
 * against a direct dense computation of every replicate (the restricted
 * regression, the perturbed data partialled or not, OLS, cluster-robust
 * V), enumeration, p-values, determinism, and test inversion. */
#if defined(CTOOLS_CASE_wild_engine)
#undef NDEBUG
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdfe/hdfe_wild.h"

#define CHECK(cond)                                                                            \
    do                                                                                         \
    {                                                                                          \
        if (!(cond))                                                                           \
        {                                                                                      \
            fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);                    \
            exit(1);                                                                           \
        }                                                                                      \
    } while (0)

#define NROW 60
#define KX 3
#define NREP 40
#define NLEV_MAX 12

static uint64_t rng_state = 88172645463325252ull;

static double unif(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return ((double)(rng_state >> 11) + 0.5) / 9007199254740992.0;
}

static double gauss(void)
{
    double u1 = unif(), u2 = unif();
    return sqrt(-2 * log(u1)) * cos(6.283185307179586 * u2);
}

/* Inverse of a small symmetric positive definite matrix (Gauss-Jordan) */
static void invert(int n, const double *M, double *inv)
{
    double a[KX * 2 * KX];
    for (int r = 0; r < n; r++)
        for (int c = 0; c < 2 * n; c++)
            a[r * 2 * n + c] = c < n ? M[r * n + c] : (double)(c - n == r);
    for (int p = 0; p < n; p++)
    {
        double piv = a[p * 2 * n + p];
        CHECK(fabs(piv) > 1e-12);
        for (int c = 0; c < 2 * n; c++)
            a[p * 2 * n + c] /= piv;
        for (int r = 0; r < n; r++)
        {
            if (r == p)
                continue;
            double f = a[r * 2 * n + p];
            for (int c = 0; c < 2 * n; c++)
                a[r * 2 * n + c] -= f * a[p * 2 * n + c];
        }
    }
    for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++)
            inv[r * n + c] = a[r * 2 * n + n + c];
}

typedef struct
{
    double X[KX * NROW]; /* partialled, column-major */
    double y[NROW];
    double w[NROW];
    int fe[NROW];
    int vce[NROW], boot[NROW];
    int G, C;
    double A[KX * KX], beta[KX], u[NROW];
} toy;

/* Weighted demeaning within the fixed effect, in place */
static void partial(const toy *t, double *col)
{
    double sw[NLEV_MAX] = {0}, sx[NLEV_MAX] = {0};
    for (int i = 0; i < NROW; i++)
    {
        sw[t->fe[i]] += t->w[i];
        sx[t->fe[i]] += t->w[i] * col[i];
    }
    for (int i = 0; i < NROW; i++)
        col[i] -= sx[t->fe[i]] / sw[t->fe[i]];
}

static ST_retcode project_cb(void *ctx, ST_double *cols, ST_int ncols)
{
    for (ST_int c = 0; c < ncols; c++)
        partial((const toy *)ctx, cols + (size_t)c * NROW);
    return 0;
}

/* Weighted OLS of y on the columns of X listed in use (n of them) */
static void ols(const toy *t, const double *y, const int *use, int n, double *coef)
{
    double M[KX * KX], inv[KX * KX], r[KX];
    for (int a = 0; a < n; a++)
    {
        r[a] = 0;
        for (int i = 0; i < NROW; i++)
            r[a] += t->w[i] * t->X[use[a] * NROW + i] * y[i];
        for (int b = 0; b < n; b++)
        {
            M[a * n + b] = 0;
            for (int i = 0; i < NROW; i++)
                M[a * n + b] += t->w[i] * t->X[use[a] * NROW + i] * t->X[use[b] * NROW + i];
        }
    }
    invert(n, M, inv);
    for (int a = 0; a < n; a++)
    {
        coef[a] = 0;
        for (int b = 0; b < n; b++)
            coef[a] += inv[a * n + b] * r[b];
    }
}

/* Cluster-robust variance of coefficient k from residuals e (no factors);
 * vce NULL: each row its own cluster */
static double vkk(const toy *t, const int *vce, int G, const double *e, int k)
{
    static double s[NROW][KX];
    memset(s, 0, sizeof(s));
    for (int i = 0; i < NROW; i++)
        for (int l = 0; l < KX; l++)
            s[vce ? vce[i] : i][l] += t->w[i] * t->X[l * NROW + i] * e[i];
    double v = 0;
    for (int g = 0; g < (vce ? G : NROW); g++)
    {
        double ak = 0;
        for (int l = 0; l < KX; l++)
            ak += t->A[l * KX + k] * s[g][l];
        v += ak * ak;
    }
    return v;
}

/* The toy model: an effect nested in the six clusters, or crossing them */
static void make_toy(toy *t, int crossed, int coarse_boot)
{
    memset(t, 0, sizeof(*t));
    t->G = 6;
    for (int i = 0; i < NROW; i++)
    {
        t->fe[i] = crossed ? i % 5 : i % 12;
        t->vce[i] = crossed ? i / 10 : t->fe[i] / 2;
        t->boot[i] = coarse_boot ? t->vce[i] / 2 : t->vce[i];
        t->w[i] = 0.5 + unif();
    }
    t->C = coarse_boot ? 3 : 6;
    double raw[KX][NROW];
    for (int i = 0; i < NROW; i++)
    {
        raw[0][i] = gauss() + 0.3 * t->fe[i];
        raw[1][i] = gauss() + 0.2 * t->vce[i];
        raw[2][i] = 0.5 * raw[0][i] + gauss();
        t->y[i] = 0.4 * raw[0][i] - 0.2 * raw[1][i] + 0.1 * t->fe[i] +
                  gauss() * (1 + t->vce[i] % 3);
    }
    for (int l = 0; l < KX; l++)
    {
        memcpy(t->X + l * NROW, raw[l], sizeof(raw[l]));
        partial(t, t->X + l * NROW);
    }
    partial(t, t->y);
    double M[KX * KX];
    for (int a = 0; a < KX; a++)
        for (int b = 0; b < KX; b++)
        {
            M[a * KX + b] = 0;
            for (int i = 0; i < NROW; i++)
                M[a * KX + b] += t->w[i] * t->X[a * NROW + i] * t->X[b * NROW + i];
        }
    invert(KX, M, t->A);
    const int all[KX] = {0, 1, 2};
    ols(t, t->y, all, KX, t->beta);
    for (int i = 0; i < NROW; i++)
    {
        t->u[i] = t->y[i];
        for (int l = 0; l < KX; l++)
            t->u[i] -= t->X[l * NROW + i] * t->beta[l];
    }
}

/* The replicate statistic computed directly */
static double direct_tstar(const toy *t, const int *vce, int G, const int *boot, int k,
                           double b0, const double *v, int projected, int wcr)
{
    double bt[KX], ut[NROW], e[NROW], ys[NROW], bs[KX], us[NROW];
    const int all[KX] = {0, 1, 2};
    if (wcr)
    {
        int other[KX - 1], n = 0;
        double yr[NROW], g[KX];
        for (int l = 0; l < KX; l++)
            if (l != k)
                other[n++] = l;
        for (int i = 0; i < NROW; i++)
            yr[i] = t->y[i] - t->X[k * NROW + i] * b0;
        ols(t, yr, other, n, g);
        for (int l = 0, a = 0; l < KX; l++)
            bt[l] = l == k ? b0 : g[a++];
        for (int i = 0; i < NROW; i++)
        {
            ut[i] = t->y[i];
            for (int l = 0; l < KX; l++)
                ut[i] -= t->X[l * NROW + i] * bt[l];
        }
    }
    else
    {
        memcpy(bt, t->beta, sizeof(bt));
        memcpy(ut, t->u, sizeof(ut));
    }
    for (int i = 0; i < NROW; i++)
        e[i] = ut[i] * v[boot ? boot[i] : i];
    if (projected)
        partial(t, e);
    for (int i = 0; i < NROW; i++)
    {
        ys[i] = e[i];
        for (int l = 0; l < KX; l++)
            ys[i] += t->X[l * NROW + i] * bt[l];
    }
    ols(t, ys, all, KX, bs);
    for (int i = 0; i < NROW; i++)
    {
        us[i] = ys[i];
        for (int l = 0; l < KX; l++)
            us[i] -= t->X[l * NROW + i] * bs[l];
    }
    double V = vkk(t, vce, G, us, k);
    if (!(V > 0))
        return NAN;
    return (bs[k] - (wcr ? b0 : t->beta[k])) / sqrt(V);
}

/* Weights of four kinds: Rademacher, normal, Webb, Mammen */
static void draw_weights(double *v, int B, int C)
{
    static const double webb[6] = {-1.224744871391589, -1, -0.7071067811865476,
                                   0.7071067811865476, 1, 1.224744871391589};
    for (int b = 0; b < B; b++)
        for (int c = 0; c < C; c++)
        {
            double u = unif(), x;
            switch (b % 4)
            {
            case 0:
                x = u < 0.5 ? -1 : 1;
                break;
            case 1:
                x = gauss();
                break;
            case 2:
                x = webb[(int)(u * 6)];
                break;
            default:
                x = u < 0.7236067977499790 ? -0.6180339887498949 : 1.618033988749895;
            }
            v[(size_t)b * C + c] = x;
        }
}

static int close_to(double a, double b)
{
    if (isnan(a) || isnan(b))
        return isnan(a) && isnan(b);
    return fabs(a - b) <= 1e-9 * (1 + fabs(b));
}

static void setup_model(const toy *t, int robust, int projected, int wcr, const int *test,
                        int n_test, int B, hdfe_wild_model *m)
{
    hdfe_wild_input in;
    memset(&in, 0, sizeof(in));
    in.N = NROW;
    in.K = KX;
    in.X = t->X;
    in.w = t->w;
    in.resid = t->u;
    in.A = t->A;
    in.beta = t->beta;
    in.vce_id = robust ? NULL : t->vce;
    in.G = t->G;
    in.boot_id = robust ? NULL : t->boot;
    in.C = t->C;
    in.n_test = n_test;
    in.test = test;
    in.null_imposed = wcr;
    in.projected = projected;
    CHECK(hdfe_wild_setup(&in, B, m) == 0);
}

/* Every replicate, test and b0 against the direct computation */
static void check_design(int crossed, int coarse_boot, int robust)
{
    toy t;
    make_toy(&t, crossed, coarse_boot);
    const int C = robust ? NROW : t.C;
    double *v = malloc(sizeof(double) * NREP * C);
    CHECK(v);
    draw_weights(v, NREP, C);
    const int test[2] = {0, 2};
    for (int wcr = 0; wcr <= 1; wcr++)
        for (int projected = 0; projected <= 1; projected++)
        {
            hdfe_wild_model m;
            setup_model(&t, robust, projected, wcr, test, 2, NREP, &m);
            CHECK(hdfe_wild_run(&m, v, 0, NREP, 2, project_cb, &t, (size_t)1 << 20) == 0);
            for (int j = 0; j < 2; j++)
            {
                const int k = test[j];
                const double b0s[4] = {0.0, t.beta[k], t.beta[k] + 0.3, -0.7};
                const double V = vkk(&t, robust ? NULL : t.vce, t.G, t.u, k);
                for (int q = 0; q < 4; q++)
                {
                    CHECK(close_to(hdfe_wild_t(&m, j, b0s[q]), (t.beta[k] - b0s[q]) / sqrt(V)));
                    for (int b = 0; b < NREP; b++)
                    {
                        double d = direct_tstar(&t, robust ? NULL : t.vce, t.G,
                                                robust ? NULL : t.boot, k, b0s[q],
                                                v + (size_t)b * C, projected, wcr);
                        CHECK(close_to(hdfe_wild_tstar(&m, j, b, b0s[q]), d));
                    }
                }
                /* p-values: boottest's counts (the statistics equal the direct
                 * ones to rounding, checked above) */
                const double t0 = hdfe_wild_t(&m, j, 0.0);
                int lt = 0, gt = 0, absgt = 0;
                for (int b = 0; b < NREP; b++)
                {
                    const double d = hdfe_wild_tstar(&m, j, b, 0.0);
                    lt += d < t0;
                    gt += d > t0;
                    absgt += fabs(d) > fabs(t0);
                }
                ST_int F = -1;
                CHECK(hdfe_wild_pvalue(&m, j, 0.0, HDFE_WILD_SYMMETRIC, &F) ==
                      (double)absgt / NREP);
                CHECK(F == NREP);
                CHECK(hdfe_wild_pvalue(&m, j, 0.0, HDFE_WILD_EQUALTAIL, NULL) ==
                      2.0 * (lt < gt ? lt : gt) / NREP);
                CHECK(hdfe_wild_pvalue(&m, j, 0.0, HDFE_WILD_LOWER, NULL) == (double)lt / NREP);
                CHECK(hdfe_wild_pvalue(&m, j, 0.0, HDFE_WILD_UPPER, NULL) == (double)gt / NREP);
            }
            hdfe_wild_free(&m);
        }
    free(v);
}

/* Nested effects: the aggregate path is exact, so both paths agree */
static void check_nested_paths(void)
{
    toy t;
    make_toy(&t, 0, 0);
    double v[NREP * 6];
    draw_weights(v, NREP, 6);
    const int test[1] = {1};
    hdfe_wild_model a, p;
    setup_model(&t, 0, 0, 1, test, 1, NREP, &a);
    setup_model(&t, 0, 1, 1, test, 1, NREP, &p);
    CHECK(hdfe_wild_run(&a, v, 0, NREP, 1, NULL, NULL, 0) == 0);
    CHECK(hdfe_wild_run(&p, v, 0, NREP, 1, project_cb, &t, (size_t)1 << 20) == 0);
    for (int b = 0; b < NREP; b++)
        for (double b0 = -1; b0 <= 1; b0 += 0.25)
        {
            const double x = hdfe_wild_tstar(&a, 0, b, b0), y = hdfe_wild_tstar(&p, 0, b, b0);
            CHECK(fabs(x - y) < 1e-10 * (1 + fabs(x)));
        }
    /* exact ties on the projected path too, once a run has set the observed
     * statistic */
    hdfe_wild_model e;
    setup_model(&t, 0, 1, 1, test, 1, 64, &e);
    CHECK(isnan(hdfe_wild_t(&e, 0, 0.0)));
    CHECK(hdfe_wild_run(&e, NULL, 0, 64, 2, project_cb, &t, (size_t)1 << 20) == 0);
    CHECK(hdfe_wild_tstar(&e, 0, 63, 0.1) == hdfe_wild_t(&e, 0, 0.1));
    CHECK(hdfe_wild_tstar(&e, 0, 0, 0.1) == -hdfe_wild_t(&e, 0, 0.1));
    hdfe_wild_free(&e);
    hdfe_wild_free(&a);
    hdfe_wild_free(&p);
}

/* Chunks, threads and projection batches leave the results unchanged */
static void check_determinism(void)
{
    toy t;
    make_toy(&t, 1, 0);
    double v[NREP * 6];
    draw_weights(v, NREP, 6);
    const int test[2] = {0, 1};
    const size_t small = 2 * 3 * NROW * sizeof(double); /* two replicates' columns */
    for (int projected = 0; projected <= 1; projected++)
    {
        hdfe_wild_model a, b;
        setup_model(&t, 0, projected, 1, test, 2, NREP, &a);
        setup_model(&t, 0, projected, 1, test, 2, NREP, &b);
        CHECK(hdfe_wild_run(&a, v, 0, NREP, 1, project_cb, &t, (size_t)1 << 24) == 0);
        CHECK(hdfe_wild_run(&b, v, 0, 17, 4, project_cb, &t, small) == 0);
        CHECK(hdfe_wild_run(&b, v + 17 * 6, 17, 1, 3, project_cb, &t, small) == 0);
        CHECK(hdfe_wild_run(&b, v + 18 * 6, 18, NREP - 18, 8, project_cb, &t, small) == 0);
        CHECK(memcmp(a.five, b.five, sizeof(double) * 5 * 2 * NREP) == 0);
        CHECK(hdfe_wild_run(&b, v, 0, 1, 1, project_cb, &t, small) == 198); /* out of order */
        hdfe_wild_free(&a);
        hdfe_wild_free(&b);
    }
    hdfe_wild_model m;
    setup_model(&t, 0, 1, 1, test, 2, NREP, &m);
    CHECK(hdfe_wild_run(&m, v, 0, 1, 1, project_cb, &t, 3 * NROW * sizeof(double) - 1) == 909);
    CHECK(hdfe_wild_run(&m, v, 0, 1, 1, NULL, NULL, (size_t)1 << 20) == 198);
    CHECK(hdfe_wild_run(&m, v, 0, NREP + 1, 1, project_cb, &t, (size_t)1 << 20) == 198);
    CHECK(isnan(hdfe_wild_pvalue(&m, 0, 0.0, HDFE_WILD_SYMMETRIC, NULL))); /* incomplete */
    hdfe_wild_free(&m);
    hdfe_wild_free(&m); /* idempotent */

    hdfe_wild_input in;
    memset(&in, 0, sizeof(in));
    int bad_vce[NROW];
    memcpy(bad_vce, t.vce, sizeof(bad_vce));
    bad_vce[7] = 6;
    in.N = NROW;
    in.K = KX;
    in.X = t.X;
    in.resid = t.u;
    in.A = t.A;
    in.beta = t.beta;
    in.vce_id = bad_vce;
    in.G = 6;
    in.n_test = 1;
    in.test = test;
    CHECK(hdfe_wild_setup(&in, 10, &m) == 198);
    in.vce_id = t.vce;
    const int bad_test[1] = {KX};
    in.test = bad_test;
    CHECK(hdfe_wild_setup(&in, 10, &m) == 198);
}

/* boottest's count_binary(N, -1, 1), row-major N x 2^N */
static void count_binary(int N, double *out)
{
    const int cols = 1 << N;
    if (N == 1)
    {
        out[0] = -1;
        out[1] = 1;
        return;
    }
    double *tmp = malloc(sizeof(double) * (size_t)(N - 1) * (size_t)(cols / 2));
    CHECK(tmp);
    count_binary(N - 1, tmp);
    for (int c = 0; c < cols; c++)
        out[c] = c < cols / 2 ? -1 : 1;
    for (int r = 1; r < N; r++)
        for (int c = 0; c < cols; c++)
            out[r * cols + c] = tmp[(r - 1) * (cols / 2) + c % (cols / 2)];
    free(tmp);
}

static void check_enumeration(void)
{
    CHECK(hdfe_wild_enumerable(4, 16) == 16);
    CHECK(hdfe_wild_enumerable(4, 15) == 0);
    CHECK(hdfe_wild_enumerable(5, 999) == 32);
    CHECK(hdfe_wild_enumerable(9, 999) == 512);
    CHECK(hdfe_wild_enumerable(10, 999) == 0);
    double ref[6 * 64], vec[6];
    count_binary(6, ref);
    for (int j = 0; j < 64; j++)
    {
        hdfe_wild_enum_weights(6, (uint64_t)j, vec);
        for (int c = 0; c < 6; c++)
            CHECK(vec[c] == ref[c * 64 + j]);
    }

    /* The engine's enumeration equals the same vectors passed in; the
     * all-ones vector (the last) reproduces the observed statistic */
    toy t;
    make_toy(&t, 0, 0);
    const int test[1] = {0};
    double v[64 * 6];
    for (int j = 0; j < 64; j++)
        hdfe_wild_enum_weights(6, (uint64_t)j, v + j * 6);
    hdfe_wild_model a, b;
    setup_model(&t, 0, 0, 1, test, 1, 64, &a);
    setup_model(&t, 0, 0, 1, test, 1, 64, &b);
    CHECK(hdfe_wild_run(&a, NULL, 0, 30, 2, NULL, NULL, 0) == 0);
    CHECK(hdfe_wild_run(&a, NULL, 30, 34, 3, NULL, NULL, 0) == 0);
    CHECK(hdfe_wild_run(&b, v, 0, 64, 1, NULL, NULL, 0) == 0);
    CHECK(memcmp(a.five, b.five, sizeof(double) * 5 * 64) == 0);
    /* v = 1 and v = -1 tie with the observed statistic exactly, so the
     * symmetric p-value never counts them (boottest's strict inequality) */
    int ties = 0;
    for (double b0 = -0.5; b0 <= 0.5; b0 += 0.125)
    {
        const double t0 = hdfe_wild_t(&a, 0, b0);
        CHECK(hdfe_wild_tstar(&a, 0, 63, b0) == t0);
        CHECK(hdfe_wild_tstar(&a, 0, 0, b0) == -t0);
        int absgt = 0;
        for (int j = 0; j < 64; j++)
        {
            absgt += fabs(hdfe_wild_tstar(&a, 0, j, b0)) > fabs(t0);
            ties += fabs(hdfe_wild_tstar(&a, 0, j, b0)) == fabs(t0);
        }
        CHECK(hdfe_wild_pvalue(&a, 0, b0, HDFE_WILD_SYMMETRIC, NULL) == absgt / 64.0);
    }
    CHECK(ties >= 2 * 9);
    /* Under WCU, v and -v (columns j and 63 - j) give opposite statistics */
    hdfe_wild_model u;
    setup_model(&t, 0, 0, 0, test, 1, 64, &u);
    CHECK(hdfe_wild_run(&u, NULL, 0, 64, 1, NULL, NULL, 0) == 0);
    for (int j = 0; j < 64; j++)
        CHECK(fabs(hdfe_wild_tstar(&u, 0, j, 0) + hdfe_wild_tstar(&u, 0, 63 - j, 0)) < 1e-12);
    /* Enumeration needs B = 2^C */
    hdfe_wild_model c;
    setup_model(&t, 0, 0, 1, test, 1, 63, &c);
    CHECK(hdfe_wild_run(&c, NULL, 0, 63, 1, NULL, NULL, 0) == 198);
    hdfe_wild_free(&a);
    hdfe_wild_free(&b);
    hdfe_wild_free(&u);
    hdfe_wild_free(&c);
}

/* Test inversion on p-functions with known confidence sets */
static double p_tent(void *ctx, double x)
{
    (void)ctx;
    double p = 1 - fabs(x - 2) / 5;
    return p > 0 ? p : 0;
}

static double p_halfline(void *ctx, double x)
{
    (void)ctx;
    return x < 0.25 ? 0.01 : 0.5;
}

static double p_two(void *ctx, double x)
{
    (void)ctx;
    return fabs(x) > 0.5 && fabs(x) < 2 ? 0.5 : 0.01;
}

static double p_none(void *ctx, double x)
{
    (void)ctx;
    (void)x;
    return 0.0;
}

static double p_engine(void *ctx, double b0)
{
    return hdfe_wild_pvalue((const hdfe_wild_model *)ctx, 0, b0, HDFE_WILD_SYMMETRIC, NULL);
}

/* Accepts everything, counting its calls (the outward search must stop at
 * 1000 scale even where center + limit - center rounds below limit) */
static double p_all(void *ctx, double x)
{
    (void)x;
    if (ctx)
        ++*(long *)ctx;
    return 0.5;
}

static void p_all_batch(void *ctx, const double *b0, const int *from, const int *to, int n,
                        double *pv)
{
    (void)from;
    (void)to;
    for (int k = 0; k < n; k++)
        pv[k] = p_all(ctx, b0[k]);
}

/* p_tent, but its evaluation fails (NaN) above 3 */
static double p_tent_fails(void *ctx, double x)
{
    return x > 3 ? NAN : p_tent(ctx, x);
}

static void p_tent_fails_batch(void *ctx, const double *b0, const int *from, const int *to, int n,
                               double *pv)
{
    (void)from;
    (void)to;
    for (int k = 0; k < n; k++)
        pv[k] = p_tent_fails(ctx, b0[k]);
}

static void check_inversion(void)
{
    hdfe_wild_ci ci;
    hdfe_wild_invert(p_tent, NULL, 2, 1, 0.05, &ci);
    CHECK(!ci.empty && !ci.disjoint && !ci.open_lo && !ci.open_hi);
    CHECK(fabs(ci.lo - (2 - 4.75)) < 2e-6 && fabs(ci.hi - (2 + 4.75)) < 2e-6);
    hdfe_wild_invert(p_halfline, NULL, 1, 1, 0.05, &ci);
    CHECK(!ci.empty && ci.open_hi && !ci.open_lo && isnan(ci.hi));
    CHECK(fabs(ci.lo - 0.25) < 2e-6);
    hdfe_wild_invert(p_two, NULL, 0, 1, 0.05, &ci);
    CHECK(ci.disjoint && !ci.open_lo && !ci.open_hi);
    CHECK(fabs(ci.lo + 2) < 2e-6 && fabs(ci.hi - 2) < 2e-6);
    hdfe_wild_invert(p_none, NULL, 0, 1, 0.05, &ci);
    CHECK(ci.empty && isnan(ci.lo) && isnan(ci.hi));
    hdfe_wild_invert(p_tent, NULL, 2, 0, 0.05, &ci);
    CHECK(ci.empty);
    CHECK(!ci.failed);

    /* An unbounded set stops at 1000 scale, also where the capped point
     * center + limit lies a rounding error short of limit from center
     * (0.3808 + 64 - 0.3808 < 64; 0.598104 + 127.8 - 0.598104 < 127.8),
     * which used to loop forever */
    const double centers[2] = {0.3808, 0.598104}, scales[2] = {0.064, 0.1278};
    for (int c = 0; c < 2; c++)
    {
        CHECK(centers[c] + 1000.0 * scales[c] - centers[c] < 1000.0 * scales[c]);
        long calls = 0;
        hdfe_wild_invert(p_all, &calls, centers[c], scales[c], 0.05, &ci);
        CHECK(ci.open_lo && ci.open_hi && isnan(ci.lo) && isnan(ci.hi) && !ci.empty);
        CHECK(calls > 0 && calls < 200);
        calls = 0;
        hdfe_wild_invert_batch(p_all_batch, &calls, centers[c], scales[c], 0.05, 16, 1e-6, &ci);
        CHECK(ci.open_lo && ci.open_hi && isnan(ci.lo) && isnan(ci.hi) && !ci.empty);
        CHECK(calls > 0 && calls < 200);
    }

    /* A failed evaluation (NaN p) is not a rejection: the set above 3 is
     * left open instead of closing at 3, and the failure is flagged */
    hdfe_wild_invert(p_tent_fails, NULL, 2, 1, 0.05, &ci);
    CHECK(ci.failed > 0 && ci.open_hi && isnan(ci.hi) && !ci.open_lo);
    CHECK(fabs(ci.lo - (2 - 4.75)) < 2e-6);
    hdfe_wild_invert_batch(p_tent_fails_batch, NULL, 2, 1, 0.05, 16, 1e-6, &ci);
    CHECK(ci.failed > 0 && ci.open_hi && isnan(ci.hi) && !ci.open_lo);
    CHECK(fabs(ci.lo - (2 - 4.75)) < 2e-6);

    /* With the engine's p: just inside each bound p > alpha, just outside not */
    toy t;
    make_toy(&t, 1, 0);
    double v[NREP * 6];
    draw_weights(v, NREP, 6);
    const int test[1] = {0};
    hdfe_wild_model m;
    setup_model(&t, 0, 1, 1, test, 1, NREP, &m);
    CHECK(hdfe_wild_run(&m, v, 0, NREP, 1, project_cb, &t, (size_t)1 << 20) == 0);
    const double se = sqrt(t.A[0]);
    hdfe_wild_invert(p_engine, &m, t.beta[0], se, 0.1, &ci);
    CHECK(!ci.empty && !ci.open_lo && !ci.open_hi);
    CHECK(ci.lo < t.beta[0] && t.beta[0] < ci.hi);
    CHECK(p_engine(&m, ci.lo + 1e-5 * se) > 0.1 && p_engine(&m, ci.lo - 1e-5 * se) <= 0.1);
    CHECK(p_engine(&m, ci.hi - 1e-5 * se) > 0.1 && p_engine(&m, ci.hi + 1e-5 * se) <= 0.1);
    hdfe_wild_free(&m);
}


/* The score bootstrap's statistic against its definition */
static double direct_score(const double *e, const int *vce, int G, const int *boot,
                           const double *sw, const double *vc, int N)
{
    double num = 0, total = 0, sg[NROW] = {0}, shg[NROW] = {0};
    for (int i = 0; i < N; i++)
    {
        const double x = vc[boot ? boot[i] : i] * e[i];
        num += x;
        sg[vce ? vce[i] : i] += x;
        shg[vce ? vce[i] : i] += sw ? sw[i] : 1;
        total += sw ? sw[i] : 1;
    }
    double den = 0;
    for (int g = 0; g < (vce ? G : N); g++)
    {
        const double d = sg[g] - shg[g] / total * num;
        den += d * d;
    }
    return num / sqrt(den);
}

static void check_score(void)
{
    double e[NROW], sw[NROW];
    int vce[NROW], boot[NROW];
    for (int i = 0; i < NROW; i++)
    {
        e[i] = gauss() + 0.3;
        sw[i] = 0.5 + unif();
        vce[i] = i / 6;      /* 10 clusters */
        boot[i] = vce[i] / 2; /* 5 coarser bootstrap clusters */
    }
    for (int design = 0; design < 3; design++)
    {
        const int *vp = design == 2 ? NULL : vce, *bp = design == 2 ? NULL : (design ? boot : vce);
        const int G = 10, C = design == 2 ? NROW : (design ? 5 : 10);
        double *v = malloc(sizeof(double) * NREP * C), t[NREP + 1];
        CHECK(v);
        draw_weights(v, NREP, C);
        hdfe_wild_score_model m;
        CHECK(hdfe_wild_score_setup(NROW, vp, G, bp, C, design == 1 ? sw : NULL, NREP, 0, &m) == 0);
        CHECK(hdfe_wild_score_weights(&m, v, 0, 15) == 0);
        CHECK(hdfe_wild_score_weights(&m, v + 15 * C, 15, NREP - 15) == 0);
        CHECK(hdfe_wild_score_t(&m, e, 3, t) == 0);
        double ones[NROW];
        for (int c = 0; c < C; c++)
            ones[c] = 1;
        CHECK(close_to(t[0], direct_score(e, vp, G, bp, design == 1 ? sw : NULL, ones, NROW)));
        int lt = 0, gt = 0, absgt = 0;
        for (int b = 0; b < NREP; b++)
        {
            CHECK(close_to(t[1 + b], direct_score(e, vp, G, bp, design == 1 ? sw : NULL,
                                                 v + (size_t)b * C, NROW)));
            lt += t[1 + b] < t[0];
            gt += t[1 + b] > t[0];
            absgt += fabs(t[1 + b]) > fabs(t[0]);
        }
        ST_int F = -1;
        CHECK(hdfe_wild_pvalue_t(t, NREP, HDFE_WILD_SYMMETRIC, &F) == (double)absgt / NREP && F == NREP);
        CHECK(hdfe_wild_pvalue_t(t, NREP, HDFE_WILD_EQUALTAIL, NULL) == 2.0 * (lt < gt ? lt : gt) / NREP);
        CHECK(hdfe_wild_pvalue_t(t, NREP, HDFE_WILD_LOWER, NULL) == (double)lt / NREP);
        CHECK(hdfe_wild_pvalue_t(t, NREP, HDFE_WILD_UPPER, NULL) == (double)gt / NREP);
        hdfe_wild_score_free(&m);
        hdfe_wild_score_free(&m);
        free(v);
    }
    /* Enumeration: the same vectors as passed in, and exact ties at v = 1 and
     * v = -1 (the last and the first) */
    double v[32 * 5], t1[33], t2[33];
    for (int j = 0; j < 32; j++)
        hdfe_wild_enum_weights(5, (uint64_t)j, v + j * 5);
    hdfe_wild_score_model a, b;
    CHECK(hdfe_wild_score_setup(NROW, vce, 10, boot, 5, NULL, 32, 1, &a) == 0);
    CHECK(hdfe_wild_score_setup(NROW, vce, 10, boot, 5, NULL, 32, 0, &b) == 0);
    CHECK(hdfe_wild_score_weights(&b, v, 0, 32) == 0);
    CHECK(hdfe_wild_score_t(&a, e, 2, t1) == 0);
    CHECK(hdfe_wild_score_t(&b, e, 1, t2) == 0);
    CHECK(memcmp(t1, t2, sizeof(t1)) == 0);
    CHECK(t1[32] == t1[0] && t1[1] == -t1[0]);
    hdfe_wild_score_free(&a);
    hdfe_wild_score_free(&b);
    /* Errors: enumeration needs B = 2^C; ids in range */
    CHECK(hdfe_wild_score_setup(NROW, vce, 10, boot, 5, NULL, 31, 1, &a) == 198);
    CHECK(hdfe_wild_score_setup(NROW, vce, 9, boot, 5, NULL, 32, 0, &a) == 198);
    /* The weights arrive in order, all of them before the statistic */
    {
        double vv[NREP * 5], tt[NREP + 1];
        draw_weights(vv, NREP, 5);
        hdfe_wild_score_model c;
        CHECK(hdfe_wild_score_setup(NROW, vce, 10, boot, 5, NULL, NREP, 0, &c) == 0);
        CHECK(hdfe_wild_score_weights(&c, vv + 5 * 5, 5, 3) == 198);
        CHECK(hdfe_wild_score_weights(&c, vv, 0, 5) == 0);
        CHECK(hdfe_wild_score_t(&c, e, 1, tt) == 198);
        CHECK(hdfe_wild_score_weights(&c, vv + 5 * 5, 5, NREP - 5) == 0);
        CHECK(hdfe_wild_score_t(&c, e, 1, tt) == 0);
        hdfe_wild_score_free(&c);
    }
}

int main(void)
{
    check_design(0, 0, 0); /* nested effect, bootstrap = vce() clusters */
    check_design(1, 0, 0); /* crossed effect */
    check_design(1, 1, 0); /* coarser bootstrap clusters */
    check_design(1, 0, 1); /* robust: every row its own cluster */
    check_nested_paths();
    check_determinism();
    check_enumeration();
    check_inversion();
    check_score();
    printf("wild_engine: all checks passed\n");
    return 0;
}
#endif

/*
 * Fixed-effect projection: the iteration on the rows against the
 * coefficient-space solver (src/hdfe/hdfe_coef.h) on synthetic designs.
 * No Stata: the solver sources are compiled into this one program.
 *
 *   cc -O3 -fopenmp -fno-fast-math -ffp-contract=off -DSYSTEM=STUNIX \
 *      -DSD_FASTMODE -Isrc validation/performance/benchmark_hdfe_coef.c -lm
 *   (-DSYSTEM=STWIN32 and -static on Windows; APPLEMAC on macOS)
 *
 *   benchmark_hdfe_coef <design> <N> <K> <threads> <tolerance> <repeats>
 *                       <weights> <methods> [maxiter]
 *
 * design   one effect: one onefirm onefew (rows grouped by it) onecyc onernd;
 *          two: readme fy panel wy rr sparse akmchain akmrand akmreg akmhard
 *          akmsmall school; three: wfy rry fgy pyg trade
 * weights  0 none, 1 uniform(0.5, 1.5), 2 lognormal, 3 one or 2^26, 4 a
 *          third of the rows weigh nothing
 * methods  comma-separated: map (the rows), default (partial_out_columns),
 *          prepared (partial_out_prepare, then partial_out_columns_prepared:
 *          the same residuals as default; a second line gives the fastest
 *          time of each step), auto dense pcg lap (hdfe_coef_method; two or
 *          more effects)
 *
 * Every third column has level means a thousand times its noise. Each line
 * reports the status, the fastest and the median time, the iterations, the
 * largest correlation left with any level (orth), the distance from a
 * tight-tolerance solve relative to its norm (err) and a hash of the output
 * bits, which must not change with the thread count.
 *
 * -DGUARD puts every block the coefficient-space module allocates between
 * fences that free() checks.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

#include "stplugin.h"
static ST_plugin mock_plugin;
ST_plugin *_stata_ = &mock_plugin;
#include "ctools_cancel.c"
#include "ctools_hdfe_utils.c"
#include "hdfe/hdfe_solver.c"
#ifdef GUARD
#define FENCE 256
static void guard_fail(const char *what, void *p)
{
    fprintf(stderr, "GUARD: %s fence of block %p damaged\n", what, p);
    abort();
}
static void *guard_malloc(size_t n)
{
    unsigned char *raw = (unsigned char *)(malloc)(n + 2 * FENCE + 16);
    if (!raw)
        return NULL;
    memcpy(raw, &n, sizeof(n));
    memset(raw + 16, 0xA5, FENCE);
    memset(raw + 16 + FENCE + n, 0x5A, FENCE);
    memset(raw + 16 + FENCE, 0xCD, n);
    return raw + 16 + FENCE;
}
static void *guard_calloc(size_t a, size_t b)
{
    void *p = guard_malloc(a * b);
    if (p)
        memset(p, 0, a * b);
    return p;
}
static void guard_free(void *p)
{
    if (!p)
        return;
    unsigned char *raw = (unsigned char *)p - FENCE - 16;
    size_t n;
    memcpy(&n, raw, sizeof(n));
    for (int i = 0; i < FENCE; i++)
    {
        if (raw[16 + i] != 0xA5)
            guard_fail("front", p);
        if (raw[16 + FENCE + n + i] != 0x5A)
            guard_fail("back", p);
    }
    (free)(raw);
}
#define malloc guard_malloc
#define calloc guard_calloc
#define free guard_free
#endif
#include "hdfe/hdfe_coef_core.c"
#include "hdfe/hdfe_coef.c"
#include "hdfe/hdfe_coef_multi.c"
#include "hdfe/hdfe_coef_lap.c"
#ifdef GUARD
#undef malloc
#undef calloc
#undef free
#endif

static double now(void)
{
#ifdef _WIN32
    LARGE_INTEGER f, t;
    QueryPerformanceFrequency(&f);
    QueryPerformanceCounter(&t);
    return (double)t.QuadPart / (double)f.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
#endif
}

/* ------------------------------------------------------------------------
 * Designs
 * ---------------------------------------------------------------------- */

static uint64_t rng = 88172645463325252ULL;
static uint64_t draw(void)
{
    rng ^= rng << 13;
    rng ^= rng >> 7;
    rng ^= rng << 17;
    return rng;
}
static double uniform(void)
{
    return (double)(draw() >> 11) * (1.0 / 9007199254740992.0);
}
static int below(int n)
{
    return (int)(uniform() * n);
}
static double normal(void)
{
    return (uniform() + uniform() + uniform() + uniform() - 2.0) * 1.7320508075688772;
}
static int clamp(int x, int lo, int hi)
{
    return x < lo ? lo : x > hi ? hi : x;
}

typedef struct
{
    int N, G;
    int *id[4]; /* 0-based levels */
    int nlev[4];
} design;

/* Levels renumbered 0..n-1 in ascending order of their codes */
static void densify(int *id, int N, int *nlev)
{
    int mx = 0;
    for (int i = 0; i < N; i++)
    {
        if (id[i] > mx)
            mx = id[i];
    }
    int *map = (int *)calloc((size_t)mx + 1, sizeof(int));
    for (int i = 0; i < N; i++)
        map[id[i]] = 1;
    int n = 0;
    for (int c = 0; c <= mx; c++)
    {
        if (map[c])
            map[c] = ++n;
    }
    for (int i = 0; i < N; i++)
        id[i] = map[id[i]] - 1;
    free(map);
    *nlev = n;
}

/* Workers with obs_per rows each in firms grouped in regions of `region`
 * firms. A worker moves with probability pmove, within the region or, with
 * probability pfar, anywhere. chain: moves go to the neighbouring firm, so
 * the firms form one long path. */
static void workers_firms(design *P, int obs_per, int nfirm, double pmove, int region, double pfar,
                          int chain)
{
    const int N = P->N, W = (N + obs_per - 1) / obs_per;
    int firm = 0, dir = 1;
    for (int i = 0; i < N; i++)
    {
        const int wk = i / obs_per;
        if (i % obs_per == 0)
        {
            firm = chain ? (int)((long long)wk * nfirm / W) : below(nfirm);
            dir = (wk & 1) ? 1 : -1;
        }
        else if (uniform() < pmove)
        {
            if (chain)
                firm = clamp(firm + dir, 0, nfirm - 1);
            else if (uniform() < pfar)
                firm = below(nfirm);
            else
                firm = clamp(firm / region * region + below(region), 0, nfirm - 1);
        }
        P->id[0][i] = wk;
        P->id[1][i] = firm;
    }
}

static int design_make(design *P, const char *name, int N)
{
    memset(P, 0, sizeof(*P));
    P->N = N;
    P->G = 2;
    for (int g = 0; g < 4; g++)
        P->id[g] = (int *)malloc((size_t)N * sizeof(int));
    int *a = P->id[0], *b = P->id[1], *c = P->id[2];
    const int firms = N / 200 > 1 ? N / 200 : 2;
    if (!strcmp(name, "readme"))
    {
        /* the README benchmark: N / 200 firms in turn, 100 random groups */
        for (int i = 0; i < N; i++)
        {
            a[i] = i % firms;
            b[i] = below(100);
        }
    }
    else if (!strcmp(name, "fy"))
    {
        for (int i = 0; i < N; i++)
        {
            a[i] = i / 200;
            b[i] = below(20);
        }
    }
    else if (!strcmp(name, "panel"))
    {
        /* firms over 20 years, 15% of the firm-years missing */
        int i = 0;
        for (int f = 0; i < N; f++)
        {
            for (int t = 0; t < 20 && i < N; t++)
            {
                if (uniform() < 0.85)
                {
                    a[i] = f;
                    b[i] = t;
                    i++;
                }
            }
        }
    }
    else if (!strcmp(name, "wy"))
    {
        /* workers over ten consecutive years */
        int s = 0;
        for (int i = 0; i < N; i++)
        {
            if (i % 10 == 0)
                s = below(11);
            a[i] = i / 10;
            b[i] = s + i % 10;
        }
    }
    else if (!strcmp(name, "rr"))
    {
        for (int i = 0; i < N; i++)
        {
            a[i] = below(firms);
            b[i] = below(firms);
        }
    }
    else if (!strcmp(name, "sparse"))
    {
        /* two rows to a level of one effect, three to a level of the other */
        for (int i = 0; i < N; i++)
        {
            a[i] = below(N / 2);
            b[i] = below(N / 3);
        }
    }
    else if (!strcmp(name, "akmchain"))
        workers_firms(P, 10, firms, 0.15, 1, 0.0, 1);
    else if (!strcmp(name, "akmrand"))
        workers_firms(P, 10, firms, 0.10, 1, 1.0, 0);
    else if (!strcmp(name, "akmreg"))
        workers_firms(P, 10, firms, 0.10, 20, 0.05, 0);
    else if (!strcmp(name, "akmhard"))
        workers_firms(P, 10, firms, 0.05, 20, 0.01, 0);
    else if (!strcmp(name, "akmsmall"))
        workers_firms(P, 10, N / 20 > 1 ? N / 20 : 2, 0.10, 50, 0.02, 0);
    else if (!strcmp(name, "school"))
    {
        /* students with three rows, teachers in schools of 30 */
        const int nt = N / 60 > 30 ? N / 60 : 30;
        int school = 0;
        for (int i = 0; i < N; i++)
        {
            if (i % 3 == 0)
                school = below(nt / 30);
            const int s = uniform() < 0.02 ? below(nt / 30) : school;
            a[i] = i / 3;
            b[i] = clamp(s * 30 + below(30), 0, nt - 1);
        }
    }
    else if (!strcmp(name, "wfy"))
    {
        workers_firms(P, 10, firms, 0.10, 20, 0.05, 0);
        int s = 0;
        for (int i = 0; i < N; i++)
        {
            if (i % 10 == 0)
                s = below(11);
            c[i] = s + i % 10;
        }
        P->G = 3;
    }
    else if (!strcmp(name, "rry"))
    {
        for (int i = 0; i < N; i++)
        {
            a[i] = below(firms);
            b[i] = below(firms);
            c[i] = below(20);
        }
        P->G = 3;
    }
    else if (!strcmp(name, "fgy"))
    {
        for (int i = 0; i < N; i++)
        {
            a[i] = i % firms;
            b[i] = below(100);
            c[i] = below(20);
        }
        P->G = 3;
    }
    else if (!strcmp(name, "pyg"))
    {
        /* a panel of ten years per unit, and a group drawn for every row:
         * a level of the eliminated effect meets a few levels of the others */
        for (int i = 0; i < N; i++)
        {
            a[i] = i / 10;
            b[i] = i % 10;
            c[i] = below(100);
        }
        P->G = 3;
    }
    else if (!strncmp(name, "one", 3))
    {
        /* one effect: units of ten rows (one) or 200 (onefirm) in sorted
         * order, twenty sorted levels (onefew), firms in turn (onecyc) or
         * at random (onernd) */
        const int few = N / 20 + 1;
        for (int i = 0; i < N; i++)
        {
            a[i] = !strcmp(name, "one")       ? i / 10
                   : !strcmp(name, "onefirm") ? i / 200
                   : !strcmp(name, "onefew")  ? i / few
                   : !strcmp(name, "onecyc")  ? i % firms
                                              : below(firms);
            b[i] = 0;
        }
        if (strcmp(name, "one") && strcmp(name, "onefirm") && strcmp(name, "onefew") &&
            strcmp(name, "onecyc") && strcmp(name, "onernd"))
            return -1;
        P->G = 1;
    }
    else if (!strcmp(name, "trade"))
    {
        /* exporter-year, importer-year and pair effects of 160 countries */
        const int C = 160;
        int T = N / (C * (C - 1) * 7 / 10);
        if (T < 2)
            T = 2;
        unsigned char *alive = (unsigned char *)malloc((size_t)C * C);
        for (int x = 0; x < C * C; x++)
            alive[x] = uniform() < 0.7;
        int i = 0;
        for (int t = 0; t < T && i < N; t++)
        {
            for (int e = 0; e < C && i < N; e++)
            {
                for (int m = 0; m < C && i < N; m++)
                {
                    if (e == m || !alive[e * C + m])
                        continue;
                    a[i] = e * T + t;
                    b[i] = m * T + t;
                    c[i] = e * C + m;
                    i++;
                }
            }
        }
        P->N = i;
        free(alive);
        P->G = 3;
    }
    else
        return -1;
    for (int g = 0; g < P->G; g++)
        densify(P->id[g], P->N, &P->nlev[g]);
    return 0;
}

/* ------------------------------------------------------------------------
 * Harness
 * ---------------------------------------------------------------------- */

static void state_open(HDFE_State *S, const design *P, double *w, int K, double tolerance,
                       int maxiter, int threads)
{
    memset(S, 0, sizeof(*S));
    S->G = P->G;
    S->N = P->N;
    S->K = K;
    S->maxiter = maxiter;
    S->tolerance = tolerance;
    S->num_threads = threads < K ? threads : K;
    if (S->num_threads < 1)
        S->num_threads = 1;
    S->has_weights = w != NULL;
    S->weight_type = w ? 1 : 0;
    S->weights = w;
    S->factors = (FE_Factor *)calloc((size_t)P->G, sizeof(FE_Factor));
    S->factors_initialized = 1;
    for (int g = 0; g < P->G; g++)
    {
        FE_Factor *f = &S->factors[g];
        f->num_levels = f->max_level = P->nlev[g];
        f->has_intercept = 1;
        f->levels = (ST_int *)malloc((size_t)P->N * sizeof(ST_int));
        f->counts = (ST_double *)calloc((size_t)P->nlev[g], sizeof(ST_double));
        if (w)
            f->weighted_counts = (ST_double *)calloc((size_t)P->nlev[g], sizeof(ST_double));
        for (int i = 0; i < P->N; i++)
        {
            f->levels[i] = P->id[g][i] + 1;
            f->counts[P->id[g][i]] += 1.0;
            if (w)
                f->weighted_counts[P->id[g][i]] += w[i];
        }
    }
    if (ctools_hdfe_alloc_buffers(S, 0, K))
    {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
}
static void state_close(HDFE_State *S)
{
    S->weights = NULL;
    ctools_hdfe_state_cleanup(S);
}

/* Largest, over columns and effects, of the share of a column's (weighted)
 * norm that its level sums still account for: 0 for a projection residual */
static double orthogonality(const design *P, const double *X, int K, const double *w)
{
    const int N = P->N;
    double worst = 0.0;
    for (int k = 0; k < K; k++)
    {
        const double *x = X + (size_t)k * N;
        double ss = 0.0;
        for (int i = 0; i < N; i++)
            ss += (w ? w[i] : 1.0) * x[i] * x[i];
        for (int g = 0; g < P->G; g++)
        {
            double *s = (double *)calloc((size_t)P->nlev[g], sizeof(double));
            double *n = (double *)calloc((size_t)P->nlev[g], sizeof(double));
            for (int i = 0; i < N; i++)
            {
                const double wi = w ? w[i] : 1.0;
                s[P->id[g][i]] += wi * x[i];
                n[P->id[g][i]] += wi;
            }
            double q = 0.0;
            for (int l = 0; l < P->nlev[g]; l++)
            {
                if (n[l] > 0)
                    q += s[l] * s[l] / n[l];
            }
            free(s);
            free(n);
            const double rel = ss > 0 ? sqrt(q / ss) : 0.0;
            if (rel > worst)
                worst = rel;
        }
    }
    return worst;
}

static double distance(const double *X, const double *R, size_t len, const double *w, int N)
{
    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < len; i++)
    {
        const double wi = w ? w[i % (size_t)N] : 1.0;
        const double d = X[i] - R[i];
        num += wi * d * d;
        den += wi * R[i] * R[i];
    }
    return den > 0 ? sqrt(num / den) : 0.0;
}

static int no_break(void)
{
    return 0;
}

static int by_value(const void *a, const void *b)
{
    const double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

int main(int argc, char **argv)
{
    if (argc < 9)
    {
        fprintf(stderr,
                "usage: benchmark_hdfe_coef design N K threads tolerance repeats weights methods "
                "[maxiter]\n");
        return 2;
    }
    const char *name = argv[1];
    int N = atoi(argv[2]);
    const int K = atoi(argv[3]);
    const int threads = atoi(argv[4]);
    const double tolerance = atof(argv[5]);
    const int repeats = atoi(argv[6]) < 64 ? atoi(argv[6]) : 64;
    const int weights = atoi(argv[7]);
    char methods[256];
    snprintf(methods, sizeof(methods), "%s", argv[8]);
    const int maxiter = argc > 9 ? atoi(argv[9]) : 16000;
    if (N < 1 || K < 1 || threads < 1 || repeats < 1)
        return 2;
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
#endif
    /* as in the plugin: this thread polls the host for Break while a team
     * of the row-space paths works */
    ctools_cancel_begin(no_break);
    design P;
    if (design_make(&P, name, N))
    {
        fprintf(stderr, "unknown design %s\n", name);
        return 2;
    }
    N = P.N;
    double *w = NULL;
    if (weights)
    {
        w = (double *)malloc((size_t)N * sizeof(double));
        for (int i = 0; i < N; i++)
        {
            w[i] = weights == 1   ? 0.5 + uniform()
                   : weights == 2 ? exp(1.5 * normal())
                   : weights == 3 ? (uniform() < 0.6 ? 67108864.0 : 1.0)
                                  : (uniform() < 0.3 ? 0.0 : 1.0 + below(3));
        }
    }
    const size_t len = (size_t)N * K;
    double *X0 = (double *)malloc(len * sizeof(double));
    double *X = (double *)malloc(len * sizeof(double));
    double *R = (double *)malloc(len * sizeof(double));
    if (!X0 || !X || !R)
    {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
    for (int k = 0; k < K; k++)
    {
        const double c0 = 0.37 * (k + 1), c1 = 0.91 * (k + 1), c2 = 0.53 * (k + 1);
        const double shift = (k % 3 == 2) ? 1000.0 : 0.0;
        for (int i = 0; i < N; i++)
        {
            double v = normal() + 0.7 * sin(c0 * P.id[0][i]) + 0.7 * cos(c1 * P.id[1][i]);
            if (P.G > 2)
                v += 0.5 * sin(c2 * P.id[2][i]);
            X0[(size_t)k * N + i] = v + shift * (1 + P.id[0][i] % 7);
        }
    }
    printf("# %s N=%d K=%d G=%d levels=%d,%d,%d threads=%d tol=%g weights=%d\n", name, N, K, P.G,
           P.nlev[0], P.nlev[1], P.G > 2 ? P.nlev[2] : 0, threads, tolerance, weights);

    HDFE_State S;
    if (P.G == 1)
    {
        /* reference: the level means in extended precision */
        long double *sum = (long double *)malloc((size_t)P.nlev[0] * sizeof(long double));
        long double *den = (long double *)malloc((size_t)P.nlev[0] * sizeof(long double));
        for (int k = 0; k < K; k++)
        {
            for (int l = 0; l < P.nlev[0]; l++)
                sum[l] = den[l] = 0.0L;
            for (int i = 0; i < N; i++)
            {
                const long double wi = w ? (long double)w[i] : 1.0L;
                sum[P.id[0][i]] += wi * (long double)X0[(size_t)k * N + i];
                den[P.id[0][i]] += wi;
            }
            for (int i = 0; i < N; i++)
            {
                const int l = P.id[0][i];
                R[(size_t)k * N + i] = (double)((long double)X0[(size_t)k * N + i] -
                                                (den[l] > 0 ? sum[l] / den[l] : 0.0L));
            }
        }
        free(sum);
        free(den);
        printf("# reference: level means, orth=%.2e\n", orthogonality(&P, R, K, w));
    }
    else
    {
        /* reference: the strongest coefficient-space method at a tight tolerance */
        state_open(&S, &P, w, K, 1e-13, 200000, threads);
        memcpy(R, X0, len * sizeof(double));
        int kept = 0, most = 0;
        for (int g = 0; g < P.G; g++)
        {
            kept += P.nlev[g];
            if (P.nlev[g] > most)
                most = P.nlev[g];
        }
        kept -= most;
        const hdfe_coef_method m = kept <= 1000 ? HDFE_COEF_DENSE
                                   : P.G == 2   ? HDFE_COEF_PCG_LAP
                                                : HDFE_COEF_PCG;
        const HDFE_SolveResult r = hdfe_coef_partial_out(&S, R, N, K, m);
        printf("# reference: status=%d method=%d iters=%d orth=%.2e\n", (int)r.status, (int)m,
               (int)r.iterations, orthogonality(&P, R, K, w));
        state_close(&S);
    }
    state_open(&S, &P, w, K, tolerance, maxiter, threads);

    for (char *m = strtok(methods, ","); m; m = strtok(NULL, ","))
    {
        double times[64];
        double ahead = 1e30, columns = 1e30;
        HDFE_SolveResult r = {0, 0, 0};
        for (int rep = 0; rep < repeats; rep++)
        {
            memcpy(X, X0, len * sizeof(double));
            const double t0 = now();
            if (!strcmp(m, "map"))
                r = partial_out_columns_options(&S, X, N, K, 0, 1);
            else if (!strcmp(m, "default"))
                r = partial_out_columns(&S, X, N, K, threads);
            else if (!strcmp(m, "prepared"))
            {
                struct hdfe_coef_prepared *prepared = partial_out_prepare(&S, K);
                const double t1 = now();
                r = partial_out_columns_prepared(&S, X, N, K, threads, prepared);
                const double t2 = now();
                partial_out_prepared_free(prepared);
                if (t1 - t0 < ahead)
                    ahead = t1 - t0;
                if (t2 - t1 < columns)
                    columns = t2 - t1;
            }
            else if (!strcmp(m, "auto"))
                r = hdfe_coef_partial_out(&S, X, N, K, HDFE_COEF_AUTO);
            else if (!strcmp(m, "dense"))
                r = hdfe_coef_partial_out(&S, X, N, K, HDFE_COEF_DENSE);
            else if (!strcmp(m, "pcg"))
                r = hdfe_coef_partial_out(&S, X, N, K, HDFE_COEF_PCG);
            else if (!strcmp(m, "lap"))
                r = hdfe_coef_partial_out(&S, X, N, K, HDFE_COEF_PCG_LAP);
            else
            {
                fprintf(stderr, "unknown method %s\n", m);
                return 2;
            }
            times[rep] = now() - t0;
        }
        uint64_t hash = 1469598103934665603ULL;
        for (size_t x = 0; x < len; x++)
        {
            uint64_t bits;
            memcpy(&bits, &X[x], 8);
            hash = (hash ^ bits) * 1099511628211ULL;
        }
        const double orth = orthogonality(&P, X, K, w);
        const double err = distance(X, R, len, w, N);
        double fastest = times[0];
        for (int rep = 1; rep < repeats; rep++)
        {
            if (times[rep] < fastest)
                fastest = times[rep];
        }
        qsort(times, (size_t)repeats, sizeof(double), by_value);
        printf("%-7s rc=%-3d min=%8.4f med=%8.4f iters=%-6d orth=%.1e err=%.1e hash=%016llx\n", m,
               (int)r.status, fastest, times[repeats / 2], (int)r.iterations, orth, err,
               (unsigned long long)hash);
        if (ahead < 1e30)
            printf("        ahead=%8.4f columns=%8.4f\n", ahead, columns);
        fflush(stdout);
    }
    state_close(&S);
    free(w);
    free(X0);
    free(X);
    free(R);
    for (int g = 0; g < 4; g++)
        free(P.id[g]);
    return 0;
}

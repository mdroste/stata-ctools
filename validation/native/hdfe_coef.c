/* Coefficient-space projection: every method against an orthogonalization of
 * the dummies on small designs, the dispatch of partial_out_columns(), the
 * iteration limit, Break, thread-count invariance on large designs, and every
 * allocation failure without a leak. */
#include "newcommands_mock.h"
void ctools_msg(const char *m, const char *f, ...)
{
    (void)m;
    (void)f;
}
#include "ctools_hdfe_utils.c"
#include "hdfe/hdfe_solver.c"

enum
{
    MAXG = 4
};
typedef struct
{
    int N, G, nlev[MAXG];
    int *lv[MAXG]; /* 0-based */
    double *w;     /* NULL: unweighted */
} design;

static uint64_t rng = UINT64_C(0x9E3779B97F4A7C15);
static double uniform(void)
{
    rng ^= rng << 13;
    rng ^= rng >> 7;
    rng ^= rng << 17;
    return (double)(rng >> 11) / 9007199254740992.0;
}
static int below(int n)
{
    return (int)(uniform() * n);
}
static double normal(void)
{
    double u = uniform(), v = uniform();
    return sqrt(-2 * log(u + 1e-300)) * cos(6.283185307179586 * v);
}

static void design_free(design *d)
{
    for (int g = 0; g < MAXG; g++)
        free(d->lv[g]);
    free(d->w);
    memset(d, 0, sizeof(*d));
}

/* layout 0: every effect drawn at random (repeated cells, rows in no order);
 * 1: effects 0 and 1 in two blocks that share no level (two components);
 * 2: rows grouped by effect 0, the levels of effect 1 increasing within it
 *    (no repeated cell), one to four rows per level;
 * 3: as 0 with the last effect nested in the one before it.
 * weights 0: none, 1: positive, 2: a third of the rows weigh nothing. */
static void design_make(design *d, int N, int G, const int *nlev, int layout, int weights)
{
    memset(d, 0, sizeof(*d));
    d->G = G;
    for (int g = 0; g < G; g++)
        d->nlev[g] = nlev[g];
    if (layout == 2)
    {
        N = 0;
        int *cells = malloc((size_t)nlev[0] * 4 * sizeof(int));
        int *counts = malloc((size_t)nlev[0] * sizeof(int));
        assert(cells && counts);
        for (int f = 0; f < nlev[0]; f++)
        {
            int want = 1 + below(4), have = 0;
            for (int j = 0; j < nlev[1] && have < want; j++)
            {
                if (below(nlev[1] - j) < want - have)
                    cells[f * 4 + have++] = j;
            }
            counts[f] = have;
            N += have;
        }
        for (int g = 0; g < G; g++)
        {
            d->lv[g] = malloc((size_t)N * sizeof(int));
            assert(d->lv[g]);
        }
        for (int f = 0, i = 0; f < nlev[0]; f++)
        {
            for (int c = 0; c < counts[f]; c++, i++)
            {
                d->lv[0][i] = f;
                d->lv[1][i] = cells[f * 4 + c];
                for (int g = 2; g < G; g++)
                    d->lv[g][i] = below(nlev[g]);
            }
        }
        free(cells);
        free(counts);
    }
    else
    {
        for (int g = 0; g < G; g++)
        {
            d->lv[g] = malloc((size_t)N * sizeof(int));
            assert(d->lv[g]);
        }
        for (int i = 0; i < N; i++)
        {
            const int block = layout == 1 ? below(2) : 0;
            for (int g = 0; g < G; g++)
            {
                int l = below(nlev[g]);
                if (layout == 1 && g < 2)
                    l = block * (nlev[g] / 2) + below(nlev[g] / 2);
                if (layout == 3 && g == G - 1)
                    l = d->lv[g - 1][i] * nlev[g] / nlev[g - 1];
                d->lv[g][i] = l;
            }
        }
    }
    d->N = N;
    if (weights)
    {
        d->w = malloc((size_t)N * sizeof(double));
        assert(d->w);
        for (int i = 0; i < N; i++)
            d->w[i] = weights == 2 && below(3) == 0 ? 0.0 : 0.25 + 3 * uniform();
    }
}

/* K columns: noise, effects of every level, and in every third column level
 * means a million times the noise */
static double *columns_make(const design *d, int K)
{
    double *X = malloc((size_t)d->N * K * sizeof(double));
    assert(X);
    for (int k = 0; k < K; k++)
    {
        for (int i = 0; i < d->N; i++)
        {
            double v = normal();
            for (int g = 0; g < d->G; g++)
                v += sin(0.7 * (k + 1) * (d->lv[g][i] + 1) + g);
            if (k % 3 == 2)
                v += 1e6 * (1 + d->lv[0][i] % 5) + 3e5 * (d->lv[1][i] % 3);
            X[(size_t)k * d->N + i] = v;
        }
    }
    return X;
}

static void state_open(HDFE_State *S, const design *d, int K, double tolerance, int maxiter)
{
    memset(S, 0, sizeof(*S));
    S->G = d->G;
    S->N = d->N;
    S->K = K;
    S->maxiter = maxiter;
    S->tolerance = tolerance;
    S->num_threads = threads < K ? threads : K;
    S->has_weights = d->w != NULL;
    S->weight_type = d->w ? 1 : 0;
    S->weights = d->w;
    S->factors = calloc((size_t)d->G, sizeof(FE_Factor));
    S->factors_initialized = 1;
    assert(S->factors);
    for (int g = 0; g < d->G; g++)
    {
        FE_Factor *f = &S->factors[g];
        f->num_levels = f->max_level = d->nlev[g];
        f->has_intercept = 1;
        f->levels = malloc((size_t)d->N * sizeof(ST_int));
        f->counts = calloc((size_t)d->nlev[g], sizeof(ST_double));
        f->weighted_counts = d->w ? calloc((size_t)d->nlev[g], sizeof(ST_double)) : NULL;
        assert(f->levels && f->counts && (!d->w || f->weighted_counts));
        for (int i = 0; i < d->N; i++)
        {
            f->levels[i] = d->lv[g][i] + 1;
            f->counts[d->lv[g][i]] += 1;
            if (d->w)
                f->weighted_counts[d->lv[g][i]] += d->w[i];
        }
    }
    assert(!ctools_hdfe_alloc_buffers(S, 0, K));
}
static void state_close(HDFE_State *S)
{
    S->weights = NULL; /* the design's */
    ctools_hdfe_state_cleanup(S);
}

static double wdot(const design *d, const double *a, const double *b)
{
    double s = 0;
    for (int i = 0; i < d->N; i++)
        s += (d->w ? d->w[i] : 1.0) * a[i] * b[i];
    return s;
}

/* out = the K columns less their weighted least-squares fit on the dummies
 * of every effect: Gram-Schmidt on the dummies, each orthogonalized twice and
 * dropped when nothing of it is left */
static void reference(const design *d, const double *X, int K, double *out)
{
    const int N = d->N;
    int total = 0, nq = 0;
    for (int g = 0; g < d->G; g++)
        total += d->nlev[g];
    double *Q = malloc((size_t)total * N * sizeof(double));
    double *c = malloc((size_t)N * sizeof(double));
    assert(Q && c);
    for (int g = 0; g < d->G; g++)
    {
        for (int l = 0; l < d->nlev[g]; l++)
        {
            for (int i = 0; i < N; i++)
                c[i] = d->lv[g][i] == l;
            const double before = wdot(d, c, c);
            if (!(before > 0))
                continue;
            for (int pass = 0; pass < 2; pass++)
            {
                for (int q = 0; q < nq; q++)
                {
                    const double *col = Q + (size_t)q * N;
                    const double p = wdot(d, col, c);
                    for (int i = 0; i < N; i++)
                        c[i] -= p * col[i];
                }
            }
            const double after = wdot(d, c, c);
            if (!(after > 1e-12 * before))
                continue;
            const double scale = 1 / sqrt(after);
            for (int i = 0; i < N; i++)
                Q[(size_t)nq * N + i] = c[i] * scale;
            nq++;
        }
    }
    memcpy(out, X, (size_t)N * K * sizeof(double));
    for (int k = 0; k < K; k++)
    {
        double *x = out + (size_t)k * N;
        for (int pass = 0; pass < 2; pass++)
        {
            for (int q = 0; q < nq; q++)
            {
                const double *col = Q + (size_t)q * N;
                const double p = wdot(d, col, x);
                for (int i = 0; i < N; i++)
                    x[i] -= p * col[i];
            }
        }
    }
    free(Q);
    free(c);
}

/* Largest over the columns of the weighted distance from the reference
 * relative to the reference's weighted norm (rows that weigh nothing have no
 * unique residual). A design that fits a column exactly leaves a residual of
 * rounding only, so the norm is not taken below 1e-6 of the column's. */
static double distance(const design *d, const double *got, const double *want, const double *X,
                       int K)
{
    double worst = 0;
    for (int k = 0; k < K; k++)
    {
        double num = 0, den = 0, raw = 0;
        for (int i = 0; i < d->N; i++)
        {
            const double wi = d->w ? d->w[i] : 1.0;
            const double e = got[(size_t)k * d->N + i] - want[(size_t)k * d->N + i];
            num += wi * e * e;
            den += wi * want[(size_t)k * d->N + i] * want[(size_t)k * d->N + i];
            raw += wi * X[(size_t)k * d->N + i] * X[(size_t)k * d->N + i];
        }
        den += 1e-12 * raw;
        const double rel = den > 0 ? sqrt(num / den) : sqrt(num);
        if (rel > worst)
            worst = rel;
    }
    return worst;
}

static const hdfe_coef_method methods[] = {HDFE_COEF_AUTO, HDFE_COEF_DENSE, HDFE_COEF_PCG,
                                           HDFE_COEF_PCG_LAP};

static void small_designs(void)
{
    static const int shapes[][1 + MAXG] = {
        {2, 30, 6, 0, 0}, {2, 6, 30, 0, 0},  {2, 40, 24, 0, 0},
        {3, 25, 6, 4, 0}, {3, 6, 25, 12, 0}, {4, 20, 6, 4, 2},
    };
    for (size_t s = 0; s < sizeof(shapes) / sizeof(shapes[0]); s++)
        for (int layout = 0; layout < 4; layout++)
            for (int weights = 0; weights < 3; weights++)
                for (int K = 3; K <= 19; K += 16)
                {
                    const int G = shapes[s][0];
                    design d;
                    design_make(&d, 260, G, shapes[s] + 1, layout, weights);
                    double *X = columns_make(&d, K);
                    double *want = malloc((size_t)d.N * K * sizeof(double));
                    double *got = malloc((size_t)d.N * K * sizeof(double));
                    assert(want && got);
                    reference(&d, X, K, want);
                    HDFE_State S;
                    state_open(&S, &d, K, 1e-11, 10000);
                    assert(hdfe_coef_applies(&S));
                    for (int m = 0; m < (G == 2 ? 4 : 3); m++)
                    {
                        memcpy(got, X, (size_t)d.N * K * sizeof(double));
                        const HDFE_SolveResult r =
                            hdfe_coef_partial_out(&S, got, d.N, K, methods[m]);
                        const double off = distance(&d, got, want, X, K);
                        if (r.status || !r.converged || r.iterations < 1 || !(off < 1e-8))
                            fprintf(stderr,
                                    "shape %d layout %d weights %d K %d method %d: status %d "
                                    "iterations %d distance %g\n",
                                    (int)s, layout, weights, K, m, (int)r.status, (int)r.iterations,
                                    off);
                        assert(!r.status && r.converged && r.iterations >= 1 && off < 1e-8);
                    }
                    /* The dispatch: the coefficients by default, the rows when the limit
                     * is one on their iteration; both to the same projection */
                    memcpy(got, X, (size_t)d.N * K * sizeof(double));
                    HDFE_SolveResult r = hdfe_coef_partial_out(&S, got, d.N, K, HDFE_COEF_AUTO);
                    double *again = malloc((size_t)d.N * K * sizeof(double));
                    assert(again && !r.status);
                    memcpy(again, X, (size_t)d.N * K * sizeof(double));
                    r = partial_out_columns(&S, again, d.N, K, threads);
                    assert(!r.status && !memcmp(got, again, (size_t)d.N * K * sizeof(double)));
                    S.maxiter = HDFE_COEF_MIN_MAXITER - 1;
                    S.tolerance = 1e-9;
                    memcpy(got, X, (size_t)d.N * K * sizeof(double));
                    r = partial_out_columns_options(&S, got, d.N, K, 0, 1);
                    memcpy(again, X, (size_t)d.N * K * sizeof(double));
                    const HDFE_SolveResult rows = partial_out_columns(&S, again, d.N, K, threads);
                    assert(r.status == rows.status && r.iterations == rows.iterations);
                    assert(!memcmp(got, again, (size_t)d.N * K * sizeof(double)));
                    free(again);
                    state_close(&S);
                    free(X);
                    free(want);
                    free(got);
                    design_free(&d);
                }
}

/* What is not a set of categorical effects stays with the rows; states that
 * cannot be solved are refused before anything is read. */
static void refusals(void)
{
    const int nlev[] = {20, 5};
    design d;
    design_make(&d, 200, 2, nlev, 0, 0);
    double *X = columns_make(&d, 2);
    HDFE_State S;
    state_open(&S, &d, 2, 1e-8, 10000);
    const int base = live;
    assert(hdfe_coef_partial_out(&S, X, d.N + 1, 2, HDFE_COEF_AUTO).status == 198);
    assert(hdfe_coef_partial_out(&S, X, d.N, 0, HDFE_COEF_AUTO).status == 198);
    assert(hdfe_coef_partial_out(&S, NULL, d.N, 2, HDFE_COEF_AUTO).status == 198);
    S.tolerance = 0;
    assert(hdfe_coef_partial_out(&S, X, d.N, 2, HDFE_COEF_AUTO).status == 198);
    S.tolerance = 1e-8;
    S.maxiter = 0;
    assert(hdfe_coef_partial_out(&S, X, d.N, 2, HDFE_COEF_AUTO).status == 198);
    S.maxiter = 10000;
    for (int g = 0; g < 2; g++)
    {
        /* a level beyond the factor's count, in the eliminated and in the kept effect */
        const ST_int kept = S.factors[g].levels[17];
        S.factors[g].levels[17] = nlev[g] + 1;
        for (int m = 0; m < 4; m++)
            assert(hdfe_coef_partial_out(&S, X, d.N, 2, methods[m]).status == 198 && live == base);
        S.factors[g].levels[17] = 0;
        assert(hdfe_coef_partial_out(&S, X, d.N, 2, HDFE_COEF_AUTO).status == 198 && live == base);
        S.factors[g].levels[17] = kept;
    }
    S.G = 1;
    assert(!hdfe_coef_applies(&S));
    S.G = 2;
    S.factors[1].slope = calloc((size_t)d.N, sizeof(ST_double));
    assert(S.factors[1].slope && !hdfe_coef_applies(&S));
    assert(hdfe_coef_partial_out(&S, X, d.N, 2, HDFE_COEF_AUTO).status == 198);
    state_close(&S);
    free(X);
    design_free(&d);
}

/* The limit is honoured with r(430), finite output and no leak; a chain of
 * levels linked one to the next, which the plain iteration cannot finish in
 * its trial, is passed to the Laplacian preconditioner. */
static void iteration_limit(void)
{
    design d;
    {
        /* three cells to a level of the larger effect: neither iteration
         * is done in three steps */
        const int nlev[] = {2000, 1500};
        design_make(&d, 6000, 2, nlev, 0, 1);
        double *X = columns_make(&d, 2);
        double *got = malloc((size_t)d.N * 2 * sizeof(double));
        assert(got);
        HDFE_State S;
        state_open(&S, &d, 2, 1e-10, 3);
        const int base = live;
        for (int m = 2; m < 4; m++)
        {
            memcpy(got, X, (size_t)d.N * 2 * sizeof(double));
            const HDFE_SolveResult r = hdfe_coef_partial_out(&S, got, d.N, 2, methods[m]);
            assert(r.status == 430 && !r.converged && r.iterations == 3 && live == base);
            for (int i = 0; i < 2 * d.N; i++)
                assert(isfinite(got[i]));
        }
        state_close(&S);
        free(X);
        free(got);
        design_free(&d);
    }
    memset(&d, 0, sizeof(d));
    d.G = 2;
    d.N = 4000;
    d.nlev[0] = 2000;
    d.nlev[1] = 2001;
    for (int g = 0; g < 2; g++)
    {
        d.lv[g] = malloc((size_t)d.N * sizeof(int));
        assert(d.lv[g]);
    }
    for (int i = 0; i < d.N; i++)
    {
        d.lv[0][i] = i / 2;
        d.lv[1][i] = i / 2 + i % 2; /* level f of effect 0 meets levels f and f + 1 */
    }
    double *X = columns_make(&d, 2);
    double *got = malloc((size_t)d.N * 2 * sizeof(double));
    assert(got);
    HDFE_State S;
    state_open(&S, &d, 2, 1e-10, 100000);
    const int base = live;
    memcpy(got, X, (size_t)d.N * 2 * sizeof(double));
    const HDFE_SolveResult r = hdfe_coef_partial_out(&S, got, d.N, 2, HDFE_COEF_AUTO);
    assert(!r.status && r.converged && r.iterations > HDFE_COEF_PROBE_ITERATIONS && live == base);
    double *want = malloc((size_t)d.N * 2 * sizeof(double));
    assert(want);
    S.tolerance = 1e-13;
    memcpy(want, X, (size_t)d.N * 2 * sizeof(double));
    assert(!hdfe_coef_partial_out(&S, want, d.N, 2, HDFE_COEF_PCG_LAP).status);
    assert(distance(&d, got, want, X, 2) < 1e-7);
    state_close(&S);
    free(X);
    free(got);
    free(want);
    design_free(&d);
}

static int polls, cancel_after;
static int poll_host(void)
{
    return ++polls >= cancel_after;
}

/* Break at every point where a solve asks: status 1, nothing asked
 * afterwards, nothing left allocated. */
static void cancellation(void)
{
    static const int shapes[][1 + MAXG] = {
        {2, 60, 9, 0, 0}, {2, 1500, 1200, 0, 0}, {3, 50, 9, 5, 0}, {3, 1500, 900, 400, 0}};
    for (size_t s = 0; s < sizeof(shapes) / sizeof(shapes[0]); s++)
    {
        const int G = shapes[s][0];
        design d;
        design_make(&d, 6000, G, shapes[s] + 1, 0, 1);
        double *X = columns_make(&d, 3);
        double *got = malloc((size_t)d.N * 3 * sizeof(double));
        assert(got);
        HDFE_State S;
        state_open(&S, &d, 3, 1e-12, 100000);
        const int base = live;
        for (int m = 0; m < (G == 2 ? 4 : 3); m++)
        {
            int finished = 0;
            for (int after = 1; !finished; after++)
            {
                memcpy(got, X, (size_t)d.N * 3 * sizeof(double));
                polls = 0;
                cancel_after = after;
                ctools_cancel_begin(poll_host);
                const HDFE_SolveResult r = hdfe_coef_partial_out(&S, got, d.N, 3, methods[m]);
                ctools_cancel_end();
                assert(live == base);
                if (polls < after)
                {
                    /* the solve ended before the Break */
                    assert(!r.status && r.converged && after > 3);
                    finished = 1;
                }
                else
                    assert(r.status == 1 && !r.converged && polls == after);
            }
        }
        state_close(&S);
        free(X);
        free(got);
        design_free(&d);
    }
}

/* Every allocation fails in turn: r(920) or, where the code has a fallback,
 * the same residuals; never a leak. */
static void allocation_failures(const design *d, int K, hdfe_coef_method method)
{
    double *X = columns_make(d, K);
    double *want = malloc((size_t)d->N * K * sizeof(double));
    double *got = malloc((size_t)d->N * K * sizeof(double));
    assert(want && got);
    HDFE_State S;
    state_open(&S, d, K, 1e-9, 10000);
    memcpy(want, X, (size_t)d->N * K * sizeof(double));
    assert(!hdfe_coef_partial_out(&S, want, d->N, K, method).status);
    const int base = live;
    for (int n = 1;; n++)
    {
        memcpy(got, X, (size_t)d->N * K * sizeof(double));
        allocations = 0;
        fail_alloc = n;
        const HDFE_SolveResult r = hdfe_coef_partial_out(&S, got, d->N, K, method);
        const int reached = allocations >= n;
        fail_alloc = 0;
        assert(live == base);
        if (!reached)
        {
            assert(!r.status && n > 8);
            break;
        }
        assert(r.status == 920 || !r.status);
        if (!r.status)
            assert(!memcmp(got, want, (size_t)d->N * K * sizeof(double)));
        else
            assert(!r.converged);
    }
    state_close(&S);
    free(X);
    free(want);
    free(got);
}

static void allocation_sweeps(void)
{
    static const int shapes[][1 + MAXG] = {{2, 30, 6, 0, 0}, {3, 25, 6, 4, 0}, {4, 20, 6, 4, 2}};
    for (size_t s = 0; s < sizeof(shapes) / sizeof(shapes[0]); s++)
        for (int layout = 0; layout <= 2; layout += 2)
            for (int weights = 0; weights < 2; weights++)
            {
                const int G = shapes[s][0];
                design d;
                design_make(&d, 260, G, shapes[s] + 1, layout, weights);
                for (int m = 0; m < (G == 2 ? 4 : 3); m++)
                    allocation_failures(&d, 18, methods[m]);
                design_free(&d);
            }
}

#ifdef _OPENMP
/* Large enough for every parallel loop: the residuals of each method are the
 * same bits on one thread and on four, and agree with the iteration on the
 * rows. */
static void thread_counts(void)
{
    static const struct
    {
        int N, G, nlev[MAXG], layout, weights, K;
        hdfe_coef_method method;
    } cases[] = {
        {300000, 2, {30000, 2000, 0, 0}, 0, 1, 5, HDFE_COEF_AUTO},
        {300000, 2, {30000, 2000, 0, 0}, 0, 0, 3, HDFE_COEF_PCG_LAP},
        {300000, 2, {3000, 40, 0, 0}, 0, 1, 18, HDFE_COEF_AUTO},
        {300000, 2, {150000, 100000, 0, 0}, 0, 0, 3, HDFE_COEF_PCG},
        {300000, 2, {150000, 600, 0, 0}, 2, 1, 4, HDFE_COEF_AUTO},
        {300000, 3, {20000, 1500, 30, 0}, 0, 1, 5, HDFE_COEF_AUTO},
        {300000, 3, {3000, 300, 30, 0}, 3, 0, 4, HDFE_COEF_AUTO},
        {300000, 4, {20000, 900, 30, 7}, 0, 0, 3, HDFE_COEF_PCG},
    };
    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++)
    {
        design d;
        design_make(&d, cases[c].N, cases[c].G, cases[c].nlev, cases[c].layout, cases[c].weights);
        const int K = cases[c].K;
        const size_t bytes = (size_t)d.N * K * sizeof(double);
        double *X = columns_make(&d, K);
        double *one = malloc(bytes), *four = malloc(bytes);
        assert(one && four);
        threads = 4;
        HDFE_State S;
        state_open(&S, &d, K, 1e-9, 10000);
        HDFE_SolveResult r1, r4;
        omp_set_num_threads(1);
        memcpy(one, X, bytes);
        r1 = hdfe_coef_partial_out(&S, one, d.N, K, cases[c].method);
        omp_set_num_threads(4);
        memcpy(four, X, bytes);
        r4 = hdfe_coef_partial_out(&S, four, d.N, K, cases[c].method);
        assert(!r1.status && !r4.status && r1.iterations == r4.iterations);
        if (memcmp(one, four, bytes))
            fprintf(stderr, "case %d: residuals depend on the thread count\n", (int)c);
        assert(!memcmp(one, four, bytes));
        /* The iteration on the rows is the yardstick for the columns it
         * solves well. It measures its progress against the raw column, so
         * it stops short on those whose level means are a million times the
         * noise (1e-4 off here at this tolerance); the small designs above
         * hold such columns to the orthogonalized reference instead. */
        memcpy(one, X, bytes);
        assert(!partial_out_columns_options(&S, one, d.N, K, 0, 1).status);
        for (int k = 0; k < K; k++)
        {
            if (k % 3 == 2)
                continue;
            const size_t at = (size_t)k * d.N;
            const double off = distance(&d, four + at, one + at, X + at, 1);
            if (!(off < 1e-6))
                fprintf(stderr, "case %d column %d: distance from the row iteration %g\n", (int)c,
                        k, off);
            assert(off < 1e-6);
        }
        state_close(&S);
        free(X);
        free(one);
        free(four);
        design_free(&d);
    }
    /* allocations inside the parallel loops fail cleanly too */
    design d;
    const int nlev[] = {30000, 2000};
    design_make(&d, 300000, 2, nlev, 0, 1);
    allocation_failures(&d, 2, HDFE_COEF_PCG_LAP);
    design_free(&d);
    const int wide[] = {150000, 100000};
    design_make(&d, 300000, 2, wide, 0, 0);
    allocation_failures(&d, 2, HDFE_COEF_PCG);
    design_free(&d);
    const int three[] = {20000, 1500, 30};
    design_make(&d, 300000, 3, three, 0, 1);
    allocation_failures(&d, 2, HDFE_COEF_PCG);
    design_free(&d);
}
#endif

int main(void)
{
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    for (threads = 1; threads <= 4; threads *= 4)
    {
#ifdef _OPENMP
        omp_set_num_threads(threads);
#endif
        small_designs();
        refusals();
        iteration_limit();
        cancellation();
        allocation_sweeps();
    }
#ifdef _OPENMP
    thread_counts();
#endif
    assert(!live);
    puts("PASS HDFE coefficient-space methods, dispatch, limits, Break, threads and allocation "
         "failures");
    return 0;
}

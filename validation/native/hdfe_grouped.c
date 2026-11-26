/* Rows that arrive grouped by a fixed effect take one pass where the general
 * projections take two; both must give the same bits.
 *
 * One effect: partial_out_columns() on levels that never decrease (runs of a
 * few rows, runs longer than a block of the grouped pass, weights, more
 * columns than are walked together, a lone decrease after the rows that are
 * probed) against the scatter and gather of project_and_subtract_fe written
 * out here.
 *
 * Two or more effects: hdfe_coef_rows_rhs() on rows grouped by the eliminated
 * effect (levels without rows, a level longer than a stretch of the fused
 * pass, every arrangement of the row plan) against the two separate passes,
 * the level means of the eliminated effect and then the sums over the kept
 * effects, written out here. */
#undef NDEBUG
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include "ctools_cancel.h"
#include "ctools_hdfe_utils.c"
#include "hdfe/hdfe_solver.c"

ST_plugin *_stata_;

static uint64_t seed = UINT64_C(0x2545F4914F6CDD1D);
static double uniform(void)
{
    uint64_t z = (seed += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    z ^= z >> 31;
    return ((double)(z >> 11) + 0.5) / 9007199254740992.0;
}
static int below(int n)
{
    return (int)(uniform() * n);
}

static void set_threads(int n)
{
#ifdef _OPENMP
    omp_set_num_threads(n);
#else
    (void)n;
#endif
}

/* ------------------------------------------------------------------------
 * One effect
 * ---------------------------------------------------------------------- */

/* Levels 1, 1, .., 2, 2, ..: runs of 1 to most rows, except that the run
 * that reaches row long_at is long_rows long. Returns the number of levels. */
static int grouped_levels(ST_int *levels, int N, int most, int long_at, int long_rows)
{
    int level = 0, i = 0;
    while (i < N)
    {
        int run = 1 + below(most);
        if (long_rows > 0 && i <= long_at && long_at < i + run)
            run = long_rows;
        level++;
        for (int r = 0; r < run && i < N; r++)
            levels[i++] = level;
    }
    return level;
}

/* weights 0: none; 1: positive; 2: a third of the rows weigh nothing.
 * divide: the factor has no reciprocal counts, so the sums are divided. */
static void one_effect(const char *what, const ST_int *levels, int nlev, int N, int K, int weights,
                       int divide)
{
    HDFE_State S;
    memset(&S, 0, sizeof(S));
    S.N = N;
    S.G = 1;
    S.K = K;
    S.num_threads = 4 < K ? 4 : K;
    S.maxiter = 10000;
    S.tolerance = 1e-8;
    S.factors = calloc(1, sizeof(FE_Factor));
    assert(S.factors);
    FE_Factor *f = &S.factors[0];
    f->num_levels = f->max_level = nlev;
    f->has_intercept = 1;
    f->levels = malloc((size_t)N * sizeof(ST_int));
    f->counts = calloc((size_t)nlev, sizeof(double));
    assert(f->levels && f->counts);
    memcpy(f->levels, levels, (size_t)N * sizeof(ST_int));
    if (weights)
    {
        S.has_weights = 1;
        S.weights = malloc((size_t)N * sizeof(double));
        f->weighted_counts = calloc((size_t)nlev, sizeof(double));
        assert(S.weights && f->weighted_counts);
        for (int i = 0; i < N; i++)
            S.weights[i] = (weights == 2 && below(3) == 0) ? 0.0 : 0.25 + 2.0 * uniform();
    }
    for (int i = 0; i < N; i++)
    {
        f->counts[levels[i] - 1] += 1.0;
        if (weights)
            f->weighted_counts[levels[i] - 1] += S.weights[i];
    }
    assert(ctools_hdfe_alloc_buffers(&S, 0, K) == 0);
    if (divide)
    {
        free(f->inv_counts);
        f->inv_counts = NULL;
        free(f->inv_weighted_counts);
        f->inv_weighted_counts = NULL;
    }

    const size_t bytes = (size_t)N * K * sizeof(double);
    double *X = malloc(bytes), *want = malloc(bytes), *got = malloc(bytes);
    double *means = malloc((size_t)nlev * sizeof(double));
    assert(X && want && got && means);
    for (int k = 0; k < K; k++)
    {
        for (int i = 0; i < N; i++)
            X[(size_t)k * N + i] = (k % 3 == 2 ? 1000.0 * (levels[i] % 7) : 0.0) + uniform() - 0.5 +
                                   0.3 * (levels[i] % 11);
    }

    /* the scatter and gather of the general projection */
    const double *w = S.weights;
    const double *inv = w && f->inv_weighted_counts ? f->inv_weighted_counts : f->inv_counts;
    const double *cnt = w && f->weighted_counts ? f->weighted_counts : f->counts;
    memcpy(want, X, bytes);
    for (int k = 0; k < K; k++)
    {
        double *x = want + (size_t)k * N;
        memset(means, 0, (size_t)nlev * sizeof(double));
        for (int i = 0; i < N; i++)
        {
            if (w)
                means[levels[i] - 1] += x[i] * w[i];
            else
                means[levels[i] - 1] += x[i];
        }
        for (int l = 0; l < nlev; l++)
        {
            if (inv)
                means[l] *= inv[l];
            else if (cnt[l] > 0)
                means[l] /= cnt[l];
        }
        for (int i = 0; i < N; i++)
            x[i] -= means[levels[i] - 1];
    }

    for (int threads = 1; threads <= 4; threads += 3)
    {
        set_threads(threads);
        memcpy(got, X, bytes);
        const HDFE_SolveResult r = partial_out_columns(&S, got, N, K, S.num_threads);
        assert(r.status == 0 && r.converged && r.iterations == 1);
        if (memcmp(got, want, bytes) != 0)
        {
            fprintf(stderr, "one effect, %s, %d threads: residuals differ\n", what, threads);
            assert(0);
        }
    }
    free(X);
    free(want);
    free(got);
    free(means);
    free(S.weights);
    S.weights = NULL;
    ctools_hdfe_state_cleanup(&S);
}

static void one_effect_cases(void)
{
    enum { N = 300001 };
    ST_int *levels = malloc((size_t)N * sizeof(ST_int));
    assert(levels);
    int nlev;

    /* short runs: the blocks of the grouped pass */
    nlev = grouped_levels(levels, N, 12, -1, 0);
    one_effect("short runs", levels, nlev, N, 3, 0, 0);
    one_effect("short runs, weights", levels, nlev, N, 3, 1, 0);
    one_effect("short runs, zero weights", levels, nlev, N, 2, 2, 0);
    one_effect("short runs, eleven columns", levels, nlev, N, 11, 0, 0);
    one_effect("short runs, one column", levels, nlev, N, 1, 1, 0);
    one_effect("short runs, divided", levels, nlev, N, 3, 0, 1);
    one_effect("short runs, weights, divided", levels, nlev, N, 3, 1, 1);

    /* a run that crosses a block boundary by less than a block: uneven blocks */
    nlev = grouped_levels(levels, N, 12, 60000, 40000);
    one_effect("a run of 40000", levels, nlev, N, 3, 0, 0);
    one_effect("a run of 40000, weights", levels, nlev, N, 4, 1, 0);

    /* runs longer than that: a column is one task */
    nlev = grouped_levels(levels, N, 12, 60000, 110000);
    one_effect("a run of 110000", levels, nlev, N, 3, 0, 0);
    one_effect("a run of 110000, weights", levels, nlev, N, 3, 1, 0);
    nlev = grouped_levels(levels, N, 12, N - 1000, 90000);
    one_effect("a long last run", levels, nlev, N, 3, 0, 0);
    for (int i = 0; i < N; i++)
        levels[i] = 1 + i / (N / 7 + 1);
    one_effect("seven levels", levels, 7, N, 5, 0, 0);
    one_effect("seven levels, weights", levels, 7, N, 2, 1, 0);
    for (int i = 0; i < N; i++)
        levels[i] = 1;
    one_effect("one level", levels, 1, N, 3, 1, 0);

    /* not grouped: among the first rows, and only after them */
    nlev = 5000;
    for (int i = 0; i < N; i++)
        levels[i] = 1 + i % nlev;
    one_effect("levels in turn", levels, nlev, N, 3, 0, 0);
    one_effect("levels in turn, weights", levels, nlev, N, 3, 1, 0);
    nlev = grouped_levels(levels, N, 12, -1, 0);
    levels[200000] = levels[199999] - 1;
    one_effect("one decrease", levels, nlev, N, 3, 0, 0);
    one_effect("one decrease, weights", levels, nlev, N, 3, 2, 0);

    /* few rows */
    for (int n = 1; n <= 9; n++)
    {
        nlev = grouped_levels(levels, n, 3, -1, 0);
        one_effect("few rows", levels, nlev, n, 2, n % 2, 0);
    }
    nlev = grouped_levels(levels, 4096, 5, -1, 0);
    one_effect("4096 rows", levels, nlev, 4096, 3, 0, 0);
    nlev = grouped_levels(levels, 4097, 5, -1, 0);
    one_effect("4097 rows", levels, nlev, 4097, 3, 1, 0);
    free(levels);
}

/* ------------------------------------------------------------------------
 * The right-hand side on rows grouped by the eliminated effect
 * ---------------------------------------------------------------------- */

/* N rows grouped by B (levels 1..nB, some without rows, one level long rows
 * long when long_rows > 0), nk kept effects of m levels in all */
static void rhs_case(const char *what, int N, int nk, int m, int K, int weighted, int long_rows)
{
    hdfe_coef_rows D;
    memset(&D, 0, sizeof(D));
    ST_int *lvB = malloc((size_t)N * sizeof(ST_int));
    ST_int *lvk[HDFE_COEF_MAX_EFFECTS] = {NULL};
    assert(lvB && nk >= 1 && nk <= 3);
    int level = 0, i = 0;
    while (i < N)
    {
        int run = 1 + below(9);
        if (long_rows > 0 && i <= N / 3 && N / 3 < i + run)
            run = long_rows;
        level += 1 + (below(5) == 0); /* now and then a level without rows */
        for (int r = 0; r < run && i < N; r++)
            lvB[i++] = level;
    }
    const int nB = level + 3; /* and three after the last row */
    const int each = m / nk;
    for (int h = 0; h < nk; h++)
    {
        lvk[h] = malloc((size_t)N * sizeof(ST_int));
        assert(lvk[h]);
        for (i = 0; i < N; i++)
            lvk[h][i] = 1 + below(each);
        D.lvk[h] = lvk[h];
        D.off[h] = h * each;
    }
    D.off[nk] = nk * each;
    m = nk * each;
    double *w = weighted ? malloc((size_t)N * sizeof(double)) : NULL;
    double *invB = calloc((size_t)nB, sizeof(double));
    assert(invB && (!weighted || w));
    for (i = 0; i < N; i++)
    {
        if (w)
            w[i] = 0.25 + 2.0 * uniform();
        invB[lvB[i] - 1] += w ? w[i] : 1.0;
    }
    for (int f = 0; f < nB; f++)
        invB[f] = invB[f] > 0.0 ? 1.0 / invB[f] : 0.0;
    D.N = N;
    D.nB = nB;
    D.m = m;
    D.nk = nk;
    D.sorted = 1;
    D.lvB = lvB;
    D.w = w;
    D.invB = invB;

    double *X = malloc((size_t)N * K * sizeof(double));
    assert(X);
    for (int k = 0; k < K; k++)
    {
        for (i = 0; i < N; i++)
            X[(size_t)k * N + i] = (k % 3 == 2 ? 1000.0 * (lvB[i] % 7) : 0.0) + uniform() - 0.5 +
                                   0.2 * (lvk[0][i] % 5);
    }

    hdfe_coef_rowplan R;
    memset(&R, 0, sizeof(R));
    assert(hdfe_coef_rowplan_init(&R, &D, K) == 0);

    /* The two passes, with the row plan's own arrangement of the sums: one
     * column at a time (no chunks), or chunks of rows whose sums over the
     * kept levels are added up in order. */
    const size_t na = (size_t)nB * K, ng = (size_t)m * K;
    double *a0 = calloc(na, sizeof(double)), *g = calloc(ng, sizeof(double));
    double *f0 = calloc((size_t)K, sizeof(double)), *scale = calloc((size_t)K, sizeof(double));
    assert(a0 && g && f0 && scale);
    const int chunks = R.nrc > 0 ? R.nrc : 1;
    for (int k = 0; k < K; k++)
    {
        const double *x = X + (size_t)k * N;
        for (i = 0; i < N; i++)
            a0[(size_t)(lvB[i] - 1) * K + k] += w ? w[i] * x[i] : x[i];
        for (int f = 0; f < nB; f++)
            a0[(size_t)f * K + k] *= invB[f];
        double *gc = calloc((size_t)m, sizeof(double));
        assert(gc);
        double centred = 0.0, raw = 0.0;
        for (int c = 0; c < chunks; c++)
        {
            const int r0 = R.nrc > 0 ? R.bound[c] : 0, r1 = R.nrc > 0 ? R.bound[c + 1] : N;
            double ss = 0.0, ss_raw = 0.0;
            memset(gc, 0, (size_t)m * sizeof(double));
            for (i = r0; i < r1; i++)
            {
                const double wi = w ? w[i] : 1.0;
                const double v = x[i] - a0[(size_t)(lvB[i] - 1) * K + k];
                const double wv = wi * v;
                ss += wv * v;
                ss_raw += wi * x[i] * x[i];
                for (int h = 0; h < nk; h++)
                    gc[D.off[h] + lvk[h][i] - 1] += wv;
            }
            if (R.nrc > 0)
            {
                for (int j = 0; j < m; j++)
                    g[(size_t)j * K + k] += gc[j];
                centred += ss;
                raw += ss_raw;
            }
            else
            {
                for (int j = 0; j < m; j++)
                    g[(size_t)j * K + k] = gc[j];
                centred = ss;
                raw = ss_raw;
            }
        }
        f0[k] = centred;
        scale[k] = raw;
        free(gc);
    }

    double *a0_got = malloc(na * sizeof(double)), *g_got = malloc(ng * sizeof(double));
    double f0_got[HDFE_COEF_KMAX], scale_got[HDFE_COEF_KMAX];
    assert(a0_got && g_got && K <= HDFE_COEF_KMAX);
    for (int threads = 1; threads <= 4; threads += 3)
    {
        set_threads(threads);
        /* stale contents must not show through */
        for (size_t x = 0; x < na; x++)
            a0_got[x] = 7.0;
        for (size_t x = 0; x < ng; x++)
            g_got[x] = -3.0;
        hdfe_coef_rows_rhs(&D, &R, X, K, a0_got, g_got, f0_got, scale_got);
        if (memcmp(a0_got, a0, na * sizeof(double)) || memcmp(g_got, g, ng * sizeof(double)) ||
            memcmp(f0_got, f0, (size_t)K * sizeof(double)) ||
            memcmp(scale_got, scale, (size_t)K * sizeof(double)))
        {
            fprintf(stderr, "right-hand side, %s (%d chunks), %d threads: sums differ\n", what,
                    R.nrc, threads);
            assert(0);
        }
    }
    hdfe_coef_rowplan_free(&R);
    free(a0);
    free(g);
    free(f0);
    free(scale);
    free(a0_got);
    free(g_got);
    free(X);
    free(invB);
    free(w);
    for (int h = 0; h < nk; h++)
        free(lvk[h]);
    free(lvB);
}

static void rhs_cases(void)
{
    /* one chunk */
    rhs_case("few rows", 9000, 1, 40, 3, 0, 0);
    rhs_case("few rows, weights", 9000, 1, 40, 1, 1, 0);
    /* chunks of rows */
    rhs_case("chunks", 200000, 1, 300, 3, 0, 0);
    rhs_case("chunks, weights", 200000, 1, 300, 3, 1, 0);
    rhs_case("chunks, two kept effects", 200000, 2, 500, 5, 1, 0);
    rhs_case("chunks, three kept effects", 200000, 3, 600, 2, 0, 0);
    rhs_case("chunks, sixteen columns", 200000, 1, 300, 16, 0, 0);
    rhs_case("chunks, a level of 30000 rows", 200000, 1, 300, 3, 1, 30000);
    rhs_case("chunks, a level of 30000 rows, two kept", 200000, 2, 300, 4, 0, 30000);
    /* kept effects too large for a buffer per chunk: one column per thread */
    rhs_case("by column", 70000, 1, 600000, 2, 0, 0);
    rhs_case("by column, weights", 70000, 2, 1200000, 3, 1, 20000);
}

int main(void)
{
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    one_effect_cases();
    rhs_cases();
    puts("PASS grouped rows: one effect in one pass and the fused right-hand side, bit for bit");
    return 0;
}

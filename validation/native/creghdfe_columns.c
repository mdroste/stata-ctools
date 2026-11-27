/* The regression's columns from the loaded variables (creghdfe_column_sums
 * and creghdfe_columns_build) against the passes they replace, written out
 * here: copy the kept rows and add them up, the total sum of squares of the
 * copy, the copy divided by its standard deviation where it lies, and the
 * unrolled dot product of the result with itself. Means, sums of squares,
 * standard deviations, norms and every value of the matrix must have the
 * same bits: with and without dropped rows, factor-variable indicators,
 * weights of each type, the double-double sums, the uncentered total, a
 * matrix large enough for the stores that bypass the cache at either
 * alignment of a column, and one or four threads. */
#include "newcommands_mock.h"
#include <float.h>
#include "ctools_matrix.c"
#include "ctools_ols.c"
#include "creghdfe/creghdfe_core.c"

static uint64_t rng = UINT64_C(0x2545F4914F6CDD1D);
static double uniform(void)
{
    rng ^= rng << 13;
    rng ^= rng >> 7;
    rng ^= rng << 17;
    return (double)(rng >> 11) / 9007199254740992.0;
}

typedef struct
{
    int N_in, K;
    int drop;        /* one loaded row in this many is not kept (0: all kept) */
    int levels;      /* columns 1.. hold integers and column 1 is an indicator */
    int weight_type; /* 0 none, 1 aweight, 2 fweight, 3 pweight */
    int standardize, quad, uncentered, raw;
    int shift;       /* the matrix starts 8 bytes off a 16-byte boundary */
} spec;

static void set_threads(int n)
{
#ifdef _OPENMP
    omp_set_num_threads(n);
#else
    (void)n;
#endif
}

static void run(const char *what, spec s)
{
    const int N_in = s.N_in, K = s.K;
    double **cols = malloc((size_t)K * sizeof(*cols));
    double *level = s.levels ? malloc((size_t)K * sizeof(double)) : NULL;
    assert(cols && (!s.levels || level));
    for (int k = 0; k < K; k++)
    {
        cols[k] = malloc(((size_t)N_in + 1) * sizeof(double));
        assert(cols[k]);
        if (level)
            level[k] = k == 1 ? 3.0 : -1.0;
        for (int i = 0; i < N_in; i++)
        {
            if (s.levels && k >= 1)
                cols[k][i] = (double)(int)(uniform() * 6);
            else
                cols[k][i] = (k % 3 == 2 ? 5000.0 : 0.0) + (k + 1) * (uniform() - 0.3);
        }
    }
    int *rows = NULL, N = N_in;
    if (s.drop)
    {
        rows = malloc(((size_t)N_in + 1) * sizeof(int));
        assert(rows);
        N = 0;
        for (int i = 0; i < N_in; i++)
        {
            if ((int)(uniform() * s.drop) != 0)
                rows[N++] = i;
        }
        /* at least two rows: a single one has no standard deviation */
        if (N < 2)
        {
            rows[0] = 0;
            rows[1] = 1;
            N = 2;
        }
    }
    double *w_raw = NULL;
    if (s.weight_type)
    {
        w_raw = malloc((size_t)N * sizeof(double));
        assert(w_raw);
        for (int i = 0; i < N; i++)
            w_raw[i] = s.weight_type == 2 ? 1.0 + (int)(uniform() * 4) : 0.2 + 3.0 * uniform();
    }
    const size_t cells = (size_t)N * K;

    /* ---- the passes written out ---- */
    double *want = malloc((cells + 1) * sizeof(double));
    double *want_raw = malloc((cells + 1) * sizeof(double));
    double *w_ref = w_raw ? malloc((size_t)N * sizeof(double)) : NULL;
    double sums[16], means[16], tss[16], sd[16], norm2[16], sum_w = 0.0;
    assert(want && want_raw && K <= 16 && (!w_raw || w_ref));
    for (int k = 0; k < K; k++)
    {
        double total = 0.0;
        double *dst = want + (size_t)k * N;
        for (int ii = 0; ii < N; ii++)
        {
            double val = cols[k][rows ? rows[ii] : ii];
            if (level && level[k] >= 0)
                val = val == level[k] ? 1.0 : 0.0;
            dst[ii] = val;
            const double w = w_raw ? w_raw[ii] : 1.0;
            total += w * val;
        }
        sums[k] = total;
    }
    for (int ii = 0; ii < N && w_raw; ii++)
        sum_w += w_raw[ii];
    for (int k = 0; k < K; k++)
        means[k] = sums[k] / (w_raw ? sum_w : (double)N);
    memcpy(want_raw, want, cells * sizeof(double));
    if (w_raw)
    {
        const double scale = (double)N / sum_w;
        for (int ii = 0; ii < N; ii++)
            w_ref[ii] = (s.weight_type == 1 || s.weight_type == 3) ? w_raw[ii] * scale : w_raw[ii];
    }
    for (int k = 0; k < K; k++)
    {
        double *col = want + (size_t)k * N;
        double ss = 0.0;
        if (w_ref)
        {
            if (s.quad)
                ss = dd_sum_sq_dev_weighted(col, means[k], w_ref, N);
            else
            {
                double kc = 0.0;
                for (int i = 0; i < N; i++)
                {
                    double d = col[i] - means[k];
                    double sq = d * d;
                    double wsq = w_ref[i] * sq;
                    double t = ss + wsq;
                    kc += (ss - t) + wsq;
                    ss = t;
                }
                ss += kc;
            }
            const double df = s.weight_type == 2 ? sum_w - 1.0 : (double)(N - 1);
            tss[k] = ss;
            sd[k] = sqrt(ss / df);
        }
        else
        {
            if (s.quad)
                ss = dd_sum_sq_dev(col, means[k], N);
            else
            {
                double kc = 0.0;
                for (int i = 0; i < N; i++)
                {
                    double d = col[i] - means[k];
                    double sq = d * d;
                    double t = ss + sq;
                    kc += (ss - t) + sq;
                    ss = t;
                }
                ss += kc;
            }
            tss[k] = ss;
            sd[k] = sqrt(ss / (N - 1));
        }
        if (sd[k] < 1e-30)
            sd[k] = 1.0;
        if (s.uncentered)
        {
            double kc = 0.0;
            ss = 0.0;
            for (int i = 0; i < N; i++)
            {
                double sq = col[i] * col[i];
                double wsq = w_ref ? w_ref[i] * sq : sq;
                double t = ss + wsq;
                kc += (ss - t) + wsq;
                ss = t;
            }
            tss[k] = ss + kc;
        }
    }
    if (s.standardize)
    {
        if (isfinite(sd[0]) && sd[0] > 0.0)
        {
            int exponent;
            frexp(sd[0], &exponent);
            sd[0] = ldexp(1.0, exponent - 1);
        }
        for (int k = 0; k < K; k++)
        {
            double *col = want + (size_t)k * N;
            const double inv = 1.0 / sd[k];
            for (int i = 0; i < N; i++)
                col[i] *= inv;
        }
    }
    for (int k = 0; k < K; k++)
        norm2[k] = ctools_dot_unrolled(want + (size_t)k * N, want + (size_t)k * N, N);

    /* ---- the functions ---- */
    const double fence = 12345.678;
    double *block = malloc((cells + 4) * sizeof(double));
    double *got_raw = s.raw ? malloc((cells + 1) * sizeof(double)) : NULL;
    double *w_got = w_raw ? malloc((size_t)N * sizeof(double)) : NULL;
    assert(block && (!s.raw || got_raw) && (!w_raw || w_got));
    /* one double before and after the matrix, at either alignment */
    double *got = block + 1 + ((((uintptr_t)(block + 1) & 15) != 0) != (s.shift != 0));
    for (int threads = 1; threads <= 4; threads += 3)
    {
        set_threads(threads);
        double sums2[16], means2[16], tss2[16], sd2[16], norm22[16], sum_w2 = -1.0;
        got[-1] = got[cells] = fence;
        for (size_t x = 0; x < cells; x++)
            got[x] = -7.0;
        if (w_raw)
            memcpy(w_got, w_raw, (size_t)N * sizeof(double));
        const creghdfe_columns C = {
            .cols = (const double *const *)cols, .level = level, .rows = rows, .N = N, .K = K,
            .weights = w_got,
        };
        creghdfe_column_sums(&C, sums2, &sum_w2);
        if (w_raw)
            assert(memcmp(&sum_w2, &sum_w, sizeof(double)) == 0);
        if (memcmp(sums2, sums, (size_t)K * sizeof(double)) != 0)
        {
            fprintf(stderr, "%s, %d threads: column sums differ\n", what, threads);
            assert(0);
        }
        for (int k = 0; k < K; k++)
            means2[k] = sums2[k] / (w_raw ? sum_w2 : (double)N);
        if (w_raw && (s.weight_type == 1 || s.weight_type == 3))
        {
            const double scale = (double)N / sum_w2;
            for (int ii = 0; ii < N; ii++)
                w_got[ii] *= scale;
        }
        const creghdfe_columns_opts opts = {
            .means = means2, .weight_type = s.weight_type, .sum_weights = w_raw ? sum_w2 : 0.0,
            .standardize = s.standardize, .use_quad = s.quad, .uncentered = s.uncentered,
        };
        creghdfe_columns_build(&C, &opts, got, got_raw, tss2, sd2, norm22);
        const int same = !memcmp(got, want, cells * sizeof(double)) &&
                         !memcmp(tss2, tss, (size_t)K * sizeof(double)) &&
                         !memcmp(sd2, sd, (size_t)K * sizeof(double)) &&
                         !memcmp(norm22, norm2, (size_t)K * sizeof(double)) &&
                         (!got_raw || !memcmp(got_raw, want_raw, cells * sizeof(double)));
        if (!same || got[-1] != fence || got[cells] != fence)
        {
            fprintf(stderr, "%s, %d threads: %s\n", what, threads,
                    same ? "wrote outside the matrix" : "results differ");
            assert(0);
        }
    }
    /* ---- a single absorbed effect, projected as the column is written ----
     * (creghdfe_column_scale and creghdfe_column_absorb) against the matrix
     * above with the projection of hdfe_solver.c written out: the level sums
     * of every column in row order, times the inverse counts or divided by
     * the counts, subtracted. Levels in runs and in no order; a few levels
     * have no row. */
    for (int pattern = 0; pattern < 2 && !s.quad; pattern++)
    {
        int *lev = malloc((size_t)N * sizeof(int));
        assert(lev);
        int nlev;
        if (pattern == 0)
        {
            int l = 0, i = 0;
            while (i < N)
            {
                /* now and then a run longer than the blocks a column is read in */
                int len = 1 + (int)(uniform() * 40) + (uniform() < 0.002 ? 3000 : 0);
                l += 1 + (uniform() < 0.05);
                for (; len > 0 && i < N; len--)
                    lev[i++] = l;
            }
            nlev = l + 2;
        }
        else
        {
            nlev = N / 5 + 3;
            for (int i = 0; i < N; i++)
                lev[i] = 1 + (int)(uniform() * (nlev - 2));
        }
        double *cnt = calloc((size_t)nlev, sizeof(double));
        double *inv = malloc((size_t)nlev * sizeof(double));
        double *lmean = malloc((size_t)nlev * sizeof(double));
        double *ref = malloc((cells + 1) * sizeof(double));
        assert(cnt && inv && lmean && ref);
        for (int i = 0; i < N; i++)
            cnt[lev[i] - 1] += w_ref ? w_ref[i] : 1.0;
        for (int l = 0; l < nlev; l++)
            inv[l] = cnt[l] > 0 ? 1.0 / cnt[l] : 0.0;
        for (int use_inv = 1; use_inv >= 0; use_inv--)
        {
            for (int k = 0; k < K; k++)
            {
                const double *col = want + (size_t)k * N;
                memset(lmean, 0, (size_t)nlev * sizeof(double));
                for (int i = 0; i < N; i++)
                    lmean[lev[i] - 1] += w_ref ? col[i] * w_ref[i] : col[i];
                for (int l = 0; l < nlev; l++)
                {
                    if (use_inv)
                        lmean[l] *= inv[l];
                    else if (cnt[l] > 0)
                        lmean[l] /= cnt[l];
                }
                for (int i = 0; i < N; i++)
                    ref[(size_t)k * N + i] = col[i] - lmean[lev[i] - 1];
            }
            for (int threads = 1; threads <= 4; threads += 3)
            {
                set_threads(threads);
                double tss2[16], sd2[16], norm22[16];
                got[-1] = got[cells] = fence;
                for (size_t x = 0; x < cells; x++)
                    got[x] = -7.0;
                const creghdfe_columns C = {
                    .cols = (const double *const *)cols, .level = level, .rows = rows, .N = N,
                    .K = K, .weights = w_ref,
                };
                const creghdfe_columns_opts opts = {
                    .means = means, .weight_type = s.weight_type, .sum_weights = w_raw ? sum_w : 0.0,
                    .standardize = s.standardize, .use_quad = 0, .uncentered = s.uncentered,
                };
                const creghdfe_absorbed E = {
                    .levels = lev, .num_levels = nlev, .inv = use_inv ? inv : NULL, .counts = cnt,
                };
                int k;
#pragma omp parallel for schedule(static)
                for (k = 0; k < K; k++)
                {
                    double *scratch = malloc((size_t)nlev * sizeof(double));
                    assert(scratch);
                    creghdfe_column_scale(&C, &opts, k, tss2, sd2);
                    creghdfe_column_absorb(&C, &opts, k, &E, sd2, norm22, scratch, got);
                    free(scratch);
                }
                const int same = !memcmp(got, ref, cells * sizeof(double)) &&
                                 !memcmp(tss2, tss, (size_t)K * sizeof(double)) &&
                                 !memcmp(sd2, sd, (size_t)K * sizeof(double)) &&
                                 !memcmp(norm22, norm2, (size_t)K * sizeof(double));
                if (!same || got[-1] != fence || got[cells] != fence)
                {
                    fprintf(stderr, "%s, one effect (pattern %d, %s), %d threads: %s\n", what,
                            pattern, use_inv ? "inverse counts" : "counts", threads,
                            same ? "wrote outside the matrix" : "results differ");
                    assert(0);
                }
            }
        }
        free(lev);
        free(cnt);
        free(inv);
        free(lmean);
        free(ref);
    }
    free(block);
    free(got_raw);
    free(w_got);
    free(want);
    free(want_raw);
    free(w_ref);
    free(w_raw);
    free(rows);
    free(level);
    for (int k = 0; k < K; k++)
        free(cols[k]);
    free(cols);
}

int main(void)
{
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    /* small: ordinary stores, every block boundary and tail */
    static const int small[] = {2, 3, 7, 8, 9, 15, 16, 17, 1000, 1023, 1024, 1025, 2047, 2048, 2055, 5000};
    for (size_t x = 0; x < sizeof(small) / sizeof(small[0]); x++)
    {
        const int n = small[x];
        run("small", (spec){n, 3, 0, 0, 0, 1, 0, 0, 0, 0});
        run("small, weights", (spec){n, 2, 0, 0, 1, 1, 0, 0, 1, 0});
        run("small, dropped rows", (spec){n, 3, 5, 0, 0, 1, 0, 0, 1, 1});
        run("small, indicator", (spec){n, 4, 3, 1, 2, 1, 0, 0, 0, 0});
        run("small, not standardized", (spec){n, 2, 0, 0, 3, 0, 0, 1, 1, 0});
        run("small, quad", (spec){n, 3, 4, 1, 1, 1, 1, 1, 1, 1});
    }
    /* above the parallel threshold */
    run("medium", (spec){150001, 5, 0, 0, 0, 1, 0, 0, 0, 0});
    run("medium, fweights, dropped", (spec){150001, 5, 7, 1, 2, 1, 0, 0, 1, 0});
    run("medium, quad, pweights", (spec){150001, 3, 0, 0, 3, 1, 1, 0, 0, 0});
    run("medium, quad, uncentered", (spec){150000, 3, 9, 0, 0, 1, 1, 1, 1, 0});
    /* large: the matrix exceeds CREGHDFE_STREAM_BYTES. An odd row count
     * puts every other column 8 bytes off a 16-byte boundary. */
    run("large, odd rows", (spec){1500001, 3, 0, 0, 0, 1, 0, 0, 0, 0});
    run("large, odd rows, shifted", (spec){1500001, 3, 0, 0, 0, 1, 0, 0, 0, 1});
    run("large, even rows", (spec){1500000, 3, 0, 0, 1, 1, 0, 0, 1, 0});
    run("large, dropped rows, indicator", (spec){1700003, 3, 9, 1, 0, 1, 0, 0, 0, 1});
    run("large, weights, not standardized", (spec){1500007, 3, 0, 0, 3, 0, 0, 1, 0, 0});
    run("large, quad", (spec){1500001, 3, 0, 0, 2, 1, 1, 0, 1, 1});
    run("large, one column", (spec){4300001, 1, 0, 0, 0, 1, 0, 0, 0, 1});
    puts("PASS regression columns: sums, sums of squares, scales, norms and the matrix, bit for bit");
    return 0;
}

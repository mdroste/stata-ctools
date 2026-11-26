/* Cross products (ctools_matmul_ata/atb/atdb), compensated cross products
 * (ctools_kahan_cross and the sums beside it) and the sandwich meats built on
 * them: every result against an independent sequential reference, the two
 * loop arrangements of each kernel against each other (a failed scratch
 * allocation selects the second), thread-count invariance, and the
 * allocation failures of the variance estimators. */
#include "newcommands_mock.h"
#include <float.h>
#include "ctools_matrix.c"
#include "ctools_ols.c"

static uint64_t rng = UINT64_C(0x2545F4914F6CDD1D);
static double uniform(void)
{
    rng ^= rng << 13;
    rng ^= rng >> 7;
    rng ^= rng << 17;
    return (double)(rng >> 11) / 9007199254740992.0;
}

/* Columns of mixed scale; every third one rides on a large level, so that
 * its products cancel and a plain sum loses digits a compensated one keeps */
static double *columns(int N, int K)
{
    double *X = malloc(((size_t)N * K + 1) * sizeof(double));
    assert(X);
    for (int k = 0; k < K; k++)
    {
        const double scale = (k % 4 == 1) ? 1e5 : (k % 4 == 3 ? 1e-3 : 1.0);
        const double level = (k % 3 == 2) ? 1e4 : 0.0;
        for (int i = 0; i < N; i++)
            X[(size_t)k * N + i] = level + scale * (uniform() - 0.5);
    }
    return X;
}
static double *positive(int N)
{
    double *w = malloc(((size_t)N + 1) * sizeof(double));
    assert(w);
    for (int i = 0; i < N; i++)
        w[i] = 0.25 + 2 * uniform();
    return w;
}

/* sum of (w * x) * y in double-double with the split product, row by row,
 * and the sum of the magnitudes */
static double ref_dd(const double *x, const double *y, const double *w, int N, double *magnitude)
{
    dd_real s = {0, 0};
    double m = 0;
    for (int i = 0; i < N; i++)
    {
        const double wx = w ? w[i] * x[i] : x[i];
        const dd_real p = two_prod(wx, y[i]);
        s = dd_add_d(dd_add_d(s, p.hi), p.lo);
        m += fabs(p.hi);
    }
    if (magnitude)
        *magnitude = m;
    return s.hi + s.lo;
}

/* The compensated sum as the estimators always took it: the rows in order,
 * the product error by the split. The pair (s, c). */
static void ref_kahan(const double *x, const double *y, const double *w, int N, double *s_out,
                      double *c_out)
{
    double s = 0, c = 0;
    for (int i = 0; i < N; i++)
    {
        const double a = w ? w[i] * x[i] : x[i];
        const double p = a * y[i];
        const double pe = CTOOLS_PROD_ERR(a, y[i], p);
        const double t = s + p;
        const double v = t - s;
        const double te = (s - (t - v)) + (p - v);
        c += te + pe;
        s = t;
    }
    *s_out = s;
    *c_out = c;
}

static int same(const void *a, const void *b, size_t bytes)
{
    return bytes == 0 || memcmp(a, b, bytes) == 0;
}

/* f() once more with its first (and then second) allocation refused: the
 * same bits, and nothing left allocated */
#define SAME_WITHOUT_SCRATCH(call, out, expected, bytes)                                           \
    do {                                                                                           \
        const int base_ = live;                                                                    \
        for (int refuse_ = 1; refuse_ <= 2; refuse_++) {                                           \
            memset((out), 0x5A, (bytes));                                                          \
            allocations = 0;                                                                       \
            fail_alloc = refuse_;                                                                  \
            call;                                                                                  \
            fail_alloc = 0;                                                                        \
            assert(same((out), (expected), (bytes)) && live == base_);                             \
        }                                                                                          \
    } while (0)

#ifdef _OPENMP
#define SAME_ON_ONE_THREAD(call, out, expected, bytes)                                             \
    do {                                                                                           \
        const int threads_ = omp_get_max_threads();                                                \
        memset((out), 0x5A, (bytes));                                                              \
        omp_set_num_threads(1);                                                                    \
        call;                                                                                      \
        omp_set_num_threads(threads_);                                                             \
        assert(same((out), (expected), (bytes)));                                                  \
    } while (0)
#else
#define SAME_ON_ONE_THREAD(call, out, expected, bytes) ((void)0)
#endif

static void cross(int N, int K1, int K2, int weighted)
{
    double *A = columns(N, K1), *B = columns(N, K2), *w = weighted ? positive(N) : NULL;
    const size_t cells = (size_t)K1 * K2, bytes = cells * sizeof(double);
    double *C = malloc(bytes + 8), *C2 = malloc(bytes + 8);
    assert(C && C2);

    ctools_matmul_atdb(A, B, w, N, K1, K2, C);
    for (int j = 0; j < K2; j++)
        for (int i = 0; i < K1; i++)
        {
            double magnitude;
            const double ref = ref_dd(A + (size_t)i * N, B + (size_t)j * N, w, N, &magnitude);
            assert(fabs(C[(size_t)j * K1 + i] - ref) <= 1e-13 * magnitude);
        }
    SAME_WITHOUT_SCRATCH(ctools_matmul_atdb(A, B, w, N, K1, K2, C2), C2, C, bytes);
    SAME_ON_ONE_THREAD(ctools_matmul_atdb(A, B, w, N, K1, K2, C2), C2, C, bytes);
    if (!weighted)
    {
        ctools_matmul_atb(A, B, N, K1, K2, C2);
        assert(same(C2, C, bytes));
    }
    /* an element depends only on its two columns: the last column of A
     * against the first of B alone */
    {
        double one;
        ctools_matmul_atdb(A + (size_t)(K1 - 1) * N, B, w, N, 1, 1, &one);
        assert(same(&one, C + (K1 - 1), sizeof(double)));
    }

    /* Gram matrices: one triangle mirrored, the three entry points alike,
     * and the triangle equal to the product of A with a copy of itself */
    {
        const size_t gbytes = (size_t)K1 * K1 * sizeof(double);
        double *G = malloc(gbytes + 8), *G2 = malloc(gbytes + 8);
        double *copy = malloc(((size_t)N * K1 + 1) * sizeof(double));
        assert(G && G2 && copy);
        memcpy(copy, A, (size_t)N * K1 * sizeof(double));
        ctools_matmul_atdb(A, A, w, N, K1, K1, G);
        ctools_matmul_atdb(A, copy, w, N, K1, K1, G2);
        for (int j = 0; j < K1; j++)
            for (int i = 0; i <= j; i++)
            {
                assert(same(G + (size_t)j * K1 + i, G + (size_t)i * K1 + j, sizeof(double)));
                assert(same(G + (size_t)j * K1 + i, G2 + (size_t)j * K1 + i, sizeof(double)));
            }
        SAME_WITHOUT_SCRATCH(ctools_matmul_atdb(A, A, w, N, K1, K1, G2), G2, G, gbytes);
        SAME_ON_ONE_THREAD(ctools_matmul_atdb(A, A, w, N, K1, K1, G2), G2, G, gbytes);
        if (!weighted)
        {
            ctools_matmul_ata(A, N, K1, G2);
            assert(same(G2, G, gbytes));
            ctools_matmul_atb(A, A, N, K1, K1, G2);
            assert(same(G2, G, gbytes));
        }
        free(G);
        free(G2);
        free(copy);
    }
    free(A);
    free(B);
    free(w);
    free(C);
    free(C2);
}

static void compensated(int N, int K, int with_y)
{
    double *data = columns(N, K + 1), *w = positive(N); /* y, then X */
    const double *y = data, *X = data + N;
    const size_t bytes = (size_t)K * K * sizeof(double);
    double *xtx = malloc(bytes + 8), *xtx2 = malloc(bytes + 8);
    double *xty = malloc((size_t)K * sizeof(double) + 8), *xty2 = malloc((size_t)K * sizeof(double) + 8);
    long double *wide = malloc((size_t)K * K * sizeof(long double)), *wide2 = malloc((size_t)K * K * sizeof(long double));
    assert(xtx && xtx2 && xty && xty2 && wide && wide2);
    const int one_block = N <= CTOOLS_KAHAN_BLOCK;

    memset(xty, 0, (size_t)K * sizeof(double));
    memset(xty2, 0, (size_t)K * sizeof(double));
    ctools_kahan_cross(X, with_y ? y : NULL, N, K, xtx, with_y ? xty : NULL);
    for (int j = 0; j < K; j++)
        for (int k = -1; k < K; k++)
        {
            if (k < 0 && !with_y)
                continue;
            const double *right = k < 0 ? y : X + (size_t)k * N;
            const double got = k < 0 ? xty[j] : xtx[(size_t)j * K + k];
            double s, c, magnitude;
            ref_kahan(X + (size_t)j * N, right, NULL, N, &s, &c);
            const double sequential = s + c;
            /* the sequential sum itself: the instruction and the split agree */
            const double dot = kahan_dot(X + (size_t)j * N, right, N);
            assert(same(&dot, &sequential, sizeof(double)));
            if (one_block)
                assert(same(&got, &sequential, sizeof(double)));
            const double ref = ref_dd(X + (size_t)j * N, right, NULL, N, &magnitude);
            assert(fabs(got - ref) <= 2 * DBL_EPSILON * fabs(ref) + 1e-28 * magnitude);
            if (k >= 0)
                assert(same(&got, xtx + (size_t)k * K + j, sizeof(double)));
        }
    SAME_WITHOUT_SCRATCH(ctools_kahan_cross(X, with_y ? y : NULL, N, K, xtx2, with_y ? xty2 : NULL),
                         xtx2, xtx, bytes);
    assert(same(xty2, xty, (size_t)K * sizeof(double)));
    SAME_ON_ONE_THREAD(ctools_kahan_cross(X, with_y ? y : NULL, N, K, xtx2, with_y ? xty2 : NULL),
                       xtx2, xtx, bytes);
    assert(same(xty2, xty, (size_t)K * sizeof(double)));
    if (with_y)
    {
        compute_xtx_xty(data, N, K + 1, xtx2, xty2);
        assert(same(xtx2, xtx, bytes) && same(xty2, xty, (size_t)K * sizeof(double)));
    }

    /* the weighted form, in long double */
    for (size_t e = 0; e < (size_t)K * K; e++)
        wide[e] = wide2[e] = -1.0L;
    ctools_kahan_cross_weighted(X, w, N, K, wide);
    for (int j = 0; j < K; j++)
        for (int k = 0; k < K; k++)
        {
            const long double got = wide[(size_t)j * K + k];
            if (k < j)
            {
                assert(got == -1.0L); /* the lower triangle is the caller's */
                continue;
            }
            double s, c, magnitude;
            ref_kahan(X + (size_t)j * N, X + (size_t)k * N, w, N, &s, &c);
            if (one_block)
                assert(got == (long double)s + (long double)c);
            const double ref = ref_dd(X + (size_t)j * N, X + (size_t)k * N, w, N, &magnitude);
            assert(fabs((double)got - ref) <= 2 * DBL_EPSILON * fabs(ref) + 1e-28 * magnitude);
        }
    {
        const int base = live;
        for (int refuse = 1; refuse <= 2; refuse++)
        {
            allocations = 0;
            fail_alloc = refuse;
            ctools_kahan_cross_weighted(X, w, N, K, wide2);
            fail_alloc = 0;
            assert(live == base);
            for (size_t e = 0; e < (size_t)K * K; e++)
                assert(wide2[e] == wide[e]);
        }
    }
    free(data);
    free(w);
    free(xtx);
    free(xtx2);
    free(xty);
    free(xty2);
    free(wide);
    free(wide2);
}

/* The double-double sums against the same loops with the split product */
static void double_double(int N)
{
    double *data = columns(N, 7), *w = positive(N);
    const double *x = data + (size_t)2 * N;
    const double mean = 1e4;
    dd_real a = {0, 0}, b = {0, 0}, c = {0, 0}, d = {0, 0};
    for (int i = 0; i < N; i++)
    {
        dd_real p = two_prod(x[i], x[i]);
        a = dd_add_d(dd_add_d(a, p.hi), p.lo);
        dd_real hi = two_prod(w[i], p.hi), lo = two_prod(w[i], p.lo);
        b = dd_add_d(dd_add_d(dd_add_d(dd_add_d(b, hi.hi), hi.lo), lo.hi), lo.lo);
        const double dev = x[i] - mean;
        p = two_prod(dev, dev);
        c = dd_add_d(dd_add_d(c, p.hi), p.lo);
        hi = two_prod(w[i], p.hi);
        lo = two_prod(w[i], p.lo);
        d = dd_add_d(dd_add_d(dd_add_d(dd_add_d(d, hi.hi), hi.lo), lo.hi), lo.lo);
    }
    double want = a.hi + a.lo, got = dd_sum_sq(x, N);
    assert(same(&got, &want, sizeof(double)));
    want = b.hi + b.lo;
    got = dd_sum_sq_weighted(x, w, N);
    assert(same(&got, &want, sizeof(double)));
    want = c.hi + c.lo;
    got = dd_sum_sq_dev(x, mean, N);
    assert(same(&got, &want, sizeof(double)));
    want = d.hi + d.lo;
    got = dd_sum_sq_dev_weighted(x, mean, w, N);
    assert(same(&got, &want, sizeof(double)));

    const double *ys[6];
    dd_real many[6], blocked[6], blocked1[6];
    for (int t = 0; t < 6; t++)
        ys[t] = data + (size_t)t * N;
    for (int weighted = 0; weighted < 2; weighted++)
    {
        const double *ww = weighted ? w : NULL;
        dd_dot_weighted_multi(x, ys, 6, ww, N, many);
        assert(dd_dot_weighted_multi_blocked(x, ys, 6, ww, N, blocked) == 0);
        for (int t = 0; t < 6; t++)
        {
            dd_real s = {0, 0};
            for (int i = 0; i < N; i++)
            {
                const double wx = ww ? ww[i] * x[i] : x[i];
                const dd_real p = two_prod(wx, ys[t][i]);
                s = dd_add_d(dd_add_d(s, p.hi), p.lo);
            }
            const dd_real single = dd_dot_weighted(x, ys[t], ww, N);
            assert(same(&single, &s, sizeof(s)) && same(many + t, &s, sizeof(s)));
            assert(fabs((blocked[t].hi + blocked[t].lo) - (s.hi + s.lo)) <=
                   2 * DBL_EPSILON * fabs(s.hi));
        }
#ifdef _OPENMP
        const int threads_now = omp_get_max_threads();
        omp_set_num_threads(1);
        assert(dd_dot_weighted_multi_blocked(x, ys, 6, ww, N, blocked1) == 0);
        omp_set_num_threads(threads_now);
        assert(same(blocked1, blocked, sizeof(blocked)));
#else
        (void)blocked1;
#endif
    }
    free(data);
    free(w);
}

/* The sandwich estimators against the definition: the meat of the robust one
 * from sequential compensated sums, that of the clustered one from each
 * cluster's e'X summed in row order and the outer products added in cluster
 * order; then D * meat * D * dof in long double. */
static void sandwich(const double *D, const long double *upper, int K, double dof, double *V)
{
    long double *meat = malloc((size_t)K * K * sizeof(long double));
    long double *temp = malloc((size_t)K * K * sizeof(long double));
    assert(meat && temp);
    for (int j = 0; j < K; j++)
        for (int k = j; k < K; k++)
            meat[(size_t)j * K + k] = meat[(size_t)k * K + j] = upper[(size_t)j * K + k];
    for (int j = 0; j < K; j++)
        for (int i = 0; i < K; i++)
        {
            long double acc = 0.0L;
            for (int k = 0; k < K; k++)
                acc += (long double)D[(size_t)k * K + i] * meat[(size_t)j * K + k];
            temp[(size_t)j * K + i] = acc;
        }
    for (int j = 0; j < K; j++)
        for (int i = 0; i < K; i++)
        {
            long double acc = 0.0L;
            for (int k = 0; k < K; k++)
                acc += temp[(size_t)k * K + i] * (long double)D[(size_t)j * K + k];
            V[(size_t)j * K + i] = (double)(acc * (long double)dof);
        }
    free(meat);
    free(temp);
}

/* weight_type 0: none; 1: aweight (normalized); 2: fweight; 3: raw weights
 * (the IV convention). layout 0: clusters at random; 1: rows sorted by
 * cluster; 2: at random with some clusters empty. */
static void variance(int N, int K, int G, int weight_type, int layout)
{
    double *X = columns(N, K), *e = columns(N, 1), *w = weight_type ? positive(N) : NULL;
    double *D = malloc((size_t)K * K * sizeof(double));
    const size_t bytes = (size_t)K * K * sizeof(double);
    double *V = malloc(bytes), *want = malloc(bytes);
    int *ids = malloc(((size_t)N + 1) * sizeof(int));
    long double *upper = malloc((size_t)K * K * sizeof(long double));
    double *obs = malloc(((size_t)N + 1) * sizeof(double));
    assert(D && V && want && ids && upper && obs);
    for (int j = 0; j < K; j++)
        for (int k = 0; k <= j; k++)
            D[(size_t)j * K + k] = D[(size_t)k * K + j] = (j == k ? 2.0 : 0.25) / (1 + j + k) / N;
    for (int i = 0; i < N; i++)
    {
        if (layout == 1)
            ids[i] = (int)((long long)i * G / N);
        else
            ids[i] = (int)(uniform() * (layout == 2 && G > 2 ? G - 2 : G));
        if (weight_type == 2)
            w[i] = 1 + floor(3 * uniform());
    }
    ctools_vce_data d = {.X_eff = X, .D = D, .resid = e, .weights = w, .weight_type = weight_type == 3 ? 1 : weight_type,
                         .N = N, .K = K, .normalize_weights = weight_type == 1};
    const double dof = 1.25;
    double scale = 1.0;
    if (weight_type == 1)
    {
        double sum = 0;
        for (int i = 0; i < N; i++)
            sum += w[i];
        scale = (double)N / sum;
    }

    /* robust */
    for (int i = 0; i < N; i++)
    {
        if (!w)
            obs[i] = e[i] * e[i];
        else if (weight_type == 2)
            obs[i] = e[i] * e[i] * w[i];
        else if (weight_type == 1)
        {
            const double ew = e[i] * (w[i] * scale);
            obs[i] = ew * ew;
        }
        else
            obs[i] = w[i] * e[i] * e[i];
    }
    for (int j = 0; j < K; j++)
        for (int k = j; k < K; k++)
        {
            double s, c;
            ref_kahan(X + (size_t)j * N, X + (size_t)k * N, obs, N, &s, &c);
            upper[(size_t)j * K + k] = (long double)s + (long double)c;
        }
    sandwich(D, upper, K, dof, want);
    memset(V, 0, bytes);
    assert(ctools_vce_robust(&d, dof, V) == 0);
    if (N <= CTOOLS_KAHAN_BLOCK)
        assert(same(V, want, bytes));
    else
        for (size_t entry = 0; entry < (size_t)K * K; entry++)
            assert(fabs(V[entry] - want[entry]) <= 1e-12 * fabs(want[entry]) + 1e-300);

    /* clustered */
    for (int i = 0; i < N; i++)
        obs[i] = !w ? e[i] : (weight_type == 1 ? e[i] * (w[i] * scale) : e[i] * w[i]);
    {
        long double *sums = calloc((size_t)G * K, sizeof(long double));
        assert(sums);
        for (int i = 0; i < N; i++)
            for (int k = 0; k < K; k++)
                sums[(size_t)ids[i] * K + k] += (long double)obs[i] * (long double)X[(size_t)k * N + i];
        for (size_t entry = 0; entry < (size_t)K * K; entry++)
            upper[entry] = 0.0L;
        for (int g = 0; g < G; g++)
            for (int j = 0; j < K; j++)
                for (int k = j; k < K; k++)
                    upper[(size_t)j * K + k] += sums[(size_t)g * K + j] * sums[(size_t)g * K + k];
        free(sums);
    }
    sandwich(D, upper, K, dof, want);
    memset(V, 0, bytes);
    assert(ctools_vce_cluster(&d, ids, G, dof, V) == 0);
    assert(same(V, want, bytes));

    /* every allocation of the clustered estimator fails as 920, V untouched;
     * the robust one may instead complete with the same result when the
     * refused allocation was scratch of the cross product */
    for (int cluster = 0; cluster < 2 && N <= 5000; cluster++)
    {
        const int base = live;
        memset(want, 0, bytes);
        allocations = 0;
        assert((cluster ? ctools_vce_cluster(&d, ids, G, dof, want) : ctools_vce_robust(&d, dof, want)) == 0);
        const int count = allocations;
        assert(count > 0 && live == base);
        for (int n = 1; n <= count; n++)
        {
            for (size_t entry = 0; entry < (size_t)K * K; entry++)
                V[entry] = 123.0;
            allocations = 0;
            fail_alloc = n;
            const ST_retcode rc = cluster ? ctools_vce_cluster(&d, ids, G, dof, V) : ctools_vce_robust(&d, dof, V);
            fail_alloc = 0;
            assert(live == base);
            if (rc == 0)
                assert(!cluster && same(V, want, bytes));
            else
            {
                assert(rc == 920);
                for (size_t entry = 0; entry < (size_t)K * K; entry++)
                    assert(V[entry] == 123.0);
            }
        }
    }
    free(X);
    free(e);
    free(w);
    free(D);
    free(V);
    free(want);
    free(ids);
    free(upper);
    free(obs);
}

int main(void)
{
    mock.missval = 8.98846567431158e307;
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(4);
#endif
    /* block and chunk edges of the plain products */
    static const int rows[] = {0, 1, 7, 511, 512, 513, 4096, 7680, 7681, 20001, 70000};
    for (size_t r = 0; r < sizeof(rows) / sizeof(rows[0]); r++)
        for (int weighted = 0; weighted < 2; weighted++)
        {
            cross(rows[r], 1, 1, weighted);
            cross(rows[r], 3, 2, weighted);
            cross(rows[r], 5, 7, weighted);
        }
    /* more columns of A than one pass of the column arrangement takes */
    cross(8200, 300, 2, 0);
    cross(8200, 259, 1, 1);

    /* sub-block, block and wave edges of the compensated ones */
    static const int krows[] = {0, 1, 5, 2047, 2048, 2049, 32768, 32769, 100000};
    for (size_t r = 0; r < sizeof(krows) / sizeof(krows[0]); r++)
        for (int with_y = 0; with_y < 2; with_y++)
        {
            compensated(krows[r], 1, with_y);
            compensated(krows[r], 4, with_y);
            compensated(krows[r], 6, with_y);
        }
    compensated(65 * CTOOLS_KAHAN_BLOCK + 17, 2, 1);
    double_double(0);
    double_double(3);
    double_double(1000);
    double_double(100000);

    for (int weight_type = 0; weight_type <= 3; weight_type++)
        for (int layout = 0; layout < 3; layout++)
        {
            variance(400, 1, 2, weight_type, layout);
            variance(3000, 3, 37, weight_type, layout);
            variance(3000, 5, 2999, weight_type, layout);
            variance(60000, 4, 5000, weight_type, layout);
        }
    assert(!live);
    puts("cross products: ok");
    return 0;
}

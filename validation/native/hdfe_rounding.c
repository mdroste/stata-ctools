/* The row-space conjugate gradient at the rounding of the column it is given.
 * A column far from zero (mean 2000 or 1e6, unit within variation) is
 * projected as closely as its centered copy at tolerances from 1e-12 down to
 * one that double precision cannot meet; a column the effects explain goes to
 * zero; a projected column stays put; an exhausted iteration limit is still
 * reported; results do not depend on the thread count; and a solve that does
 * not come near the rounding is the plain Hestenes-Stiefel iteration, bit for
 * bit. Without the rounding floor the first of these fails: at tolerance
 * 1e-13 the mean-2000 column comes back a million times its size away. */
#undef NDEBUG
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include "ctools_cancel.h"
#include "ctools_hdfe_utils.c"
#include "hdfe/hdfe_solver.c"

ST_plugin *_stata_;

enum { N = 40000, K = 6, THREADS = 4 };

static uint64_t seed;
static double uniform(void)
{
    uint64_t z = (seed += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    z ^= z >> 31;
    return ((double)(z >> 11) + 0.5) / 9007199254740992.0;
}
static double normal(void)
{
    double u = uniform(), v = uniform();
    return sqrt(-2.0 * log(u)) * cos(6.283185307179586 * v);
}

/* Two crossed effects: 4000 x 400 levels drawn at random (every level used),
 * or a balanced 2000 x 20 panel, which one sweep projects exactly */
static void design(HDFE_State *s, int balanced, int weighted)
{
    memset(s, 0, sizeof(*s));
    s->N = N; s->G = 2; s->num_threads = THREADS;
    s->maxiter = 10000; s->tolerance = 1e-8;
    s->has_weights = weighted;
    if (weighted) {
        s->weights = malloc(N * sizeof(double));
        assert(s->weights);
        for (int i = 0; i < N; i++) s->weights[i] = 1 + (int)(uniform() * 5);
    }
    s->factors = calloc(2, sizeof(FE_Factor));
    assert(s->factors);
    for (int g = 0; g < 2; g++) {
        FE_Factor *f = &s->factors[g];
        f->num_levels = balanced ? (g ? 20 : 2000) : (g ? 400 : 4000);
        f->has_intercept = 1;
        f->levels = calloc(N, sizeof(ST_int));
        f->counts = calloc(f->num_levels, sizeof(double));
        f->weighted_counts = calloc(f->num_levels, sizeof(double));
        assert(f->levels && f->counts && f->weighted_counts);
        for (int i = 0; i < N; i++) {
            int l = balanced ? (g ? i / 2000 : i % 2000)
                             : i < f->num_levels ? i : (int)(uniform() * f->num_levels);
            f->levels[i] = 1 + l;
            f->counts[l]++;
            f->weighted_counts[l] += weighted ? s->weights[i] : 1;
        }
    }
    assert(!ctools_hdfe_alloc_buffers(s, 0, K));
}

static double dot(const HDFE_State *s, const double *a, const double *b)
{
    return s->weights ? weighted_dot_product(a, b, s->weights, N) : dot_product(a, b, N);
}

/* The Hestenes-Stiefel iteration with no rounding floor, on the solver's own
 * kernels: what cg_solve_column_threaded() must be while it stays above it */
static int plain_cg(HDFE_State *s, double *y)
{
    double *r = s->thread_cg_r[0], *u = s->thread_cg_u[0], *v = s->thread_cg_v[0];
    double potential = dot(s, y, y), tol_sq = s->tolerance * s->tolerance;
    transform_sym_kaczmarz_threaded(s, y, r, s->thread_fe_means, 0);
    double ssr = dot(s, r, r);
    memcpy(u, r, N * sizeof(double));
    for (int iter = 1; iter <= s->maxiter; iter++) {
        transform_sym_kaczmarz_threaded(s, u, v, s->thread_fe_means, 0);
        double uv = dot(s, u, v), alpha = fabs(uv) < 1e-30 ? 0.0 : ssr / uv;
        double recent = alpha * ssr, before = ssr;
        potential -= recent;
        ssr = s->weights ? ctools_simd_fused_axpy2_wdot(y, u, r, v, s->weights, alpha, (size_t)N)
                         : ctools_simd_fused_axpy2_dot(y, u, r, v, alpha, (size_t)N);
        if (potential <= 1e-30 || recent <= tol_sq * potential) return iter;
        ctools_simd_axpby(u, r, 1.0, fabs(before) < 1e-30 ? 0.0 : ssr / before, (size_t)N);
    }
    return -s->maxiter;
}

static double rms(const double *a)
{
    double sum = 0;
    for (int i = 0; i < N; i++) sum += a[i] * a[i];
    return sqrt(sum / N);
}
/* largest difference between two columns in units of the second one's size */
static double apart(const double *a, const double *b)
{
    double most = 0, size = rms(b);
    for (int i = 0; i < N; i++) most = fmax(most, fabs(a[i] - b[i]));
    return most / size;
}
static int solve(HDFE_State *s, double tolerance, const double *in, double *out)
{
    memcpy(out, in, N * sizeof(double));
    s->tolerance = tolerance;
    return cg_solve_column_threaded(s, out, 0);
}

static void rounding_case(int balanced, int weighted)
{
    HDFE_State s;
    seed = 4711 + 2 * balanced + weighted;
    design(&s, balanced, weighted);
    double *col = malloc((size_t)N * 8 * sizeof(double)), *data = malloc((size_t)N * K * 2 * sizeof(double));
    assert(col && data);
    double *centered = col, *far = col + N, *farther = col + 2 * N, *near = col + 3 * N;
    double *explained = col + 4 * N, *base = col + 5 * N, *got = col + 6 * N, *zero = col + 7 * N;
    double mean = 0;
    for (int i = 0; i < N; i++) {
        int a = s.factors[0].levels[i] - 1, b = s.factors[1].levels[i] - 1;
        centered[i] = normal() + (a % 7) / 3.0;
        explained[i] = 2000 + (a % 11) / 4.0 + (b % 13) / 5.0;
        zero[i] = 0;
        mean += centered[i] / N;
    }
    for (int i = 0; i < N; i++) {
        centered[i] -= mean;
        near[i] = centered[i] + 40;
        far[i] = centered[i] + 2000;
        farther[i] = centered[i] + 1e6;
    }

    /* The column with a level against its centered copy, tolerance by
     * tolerance; 1e-16 is below what double precision can deliver */
    const double tolerances[] = {1e-12, 1e-13, 1e-14, 1e-16};
    for (int t = 0; t < 4; t++) {
        assert(solve(&s, tolerances[t], centered, base) > 0);
        assert(solve(&s, tolerances[t], near, got) > 0);
        assert(apart(got, base) < 1e-9);
        assert(solve(&s, tolerances[t], far, got) > 0);
        assert(apart(got, base) < 1e-9);
        /* 1e6 + e holds e itself to 6e-11 only */
        assert(solve(&s, tolerances[t], farther, got) > 0);
        assert(apart(got, base) < 1e-7);
        /* nothing is left of a column the effects explain */
        assert(solve(&s, tolerances[t], explained, got) > 0);
        double most = 0;
        for (int i = 0; i < N; i++) most = fmax(most, fabs(got[i]));
        assert(most < 1e-12 * rms(explained));
        /* and a projected column is left where it is */
        assert(solve(&s, tolerances[t], base, got) > 0);
        assert(apart(got, base) < 1e-13);
    }

    /* Away from the rounding the solver is the plain iteration, bit for bit */
    if (!balanced)
        for (int t = 0; t < 3; t++) {
            const double *in = t == 2 ? near : centered;
            const double tolerance = t == 1 ? 1e-10 : 1e-8;
            int iterations = solve(&s, tolerance, in, got);
            memcpy(base, in, N * sizeof(double));
            assert(plain_cg(&s, base) == iterations && iterations > 2);
            assert(!memcmp(got, base, N * sizeof(double)));
        }

    /* All columns at once, on several threads and on one: the same bits, and
     * those of the single-column solves above */
    const double *inputs[K] = {centered, near, far, farther, explained, zero};
    s.tolerance = 1e-13;
    for (int pass = 0; pass < 2; pass++) {
        double *block = data + (size_t)pass * N * K;
        for (int k = 0; k < K; k++) memcpy(block + (size_t)k * N, inputs[k], N * sizeof(double));
        s.num_threads = pass ? 1 : THREADS;
        HDFE_SolveResult result = partial_out_columns_options(&s, block, N, K, 0, 1);
        assert(!result.status && result.converged && result.iterations > 0);
    }
    assert(!memcmp(data, data + (size_t)N * K, (size_t)N * K * sizeof(double)));
    assert(solve(&s, 1e-13, far, got) > 0 && !memcmp(got, data + (size_t)2 * N, N * sizeof(double)));

    /* One iteration is not enough for the random design, floor or no floor */
    if (!balanced) {
        s.maxiter = 1;
        s.tolerance = 1e-8;
        for (int k = 0; k < K; k++) memcpy(data + (size_t)k * N, inputs[k], N * sizeof(double));
        HDFE_SolveResult result = partial_out_columns_options(&s, data, N, K, 0, 1);
        assert(result.status == 430 && !result.converged);
    }
    free(col);
    free(data);
    free(s.weights);
    s.weights = NULL;
    ctools_hdfe_state_cleanup(&s);
}

int main(void)
{
    for (int balanced = 0; balanced < 2; balanced++)
        for (int weighted = 0; weighted < 2; weighted++)
            rounding_case(balanced, weighted);
    puts("PASS HDFE row-space rounding floor: levels 40, 2000 and 1e6, explained and projected columns, "
         "plain-iteration parity, threads and iteration limit");
    return 0;
}

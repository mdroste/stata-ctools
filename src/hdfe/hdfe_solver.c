/*
 * hdfe_solver.c
 *
 * CG solver, Kaczmarz transforms, FE projections
 * Highly optimized with OpenMP parallelization and SIMD hints
 * Part of the ctools Stata plugin suite
 */

#include "hdfe_solver.h"
#include "../ctools_unroll.h"
#include "../ctools_simd.h"

/* ========================================================================
 * Dot product
 *
 * On AVX2: explicit SIMD with FMA (4 doubles/vector, 2 accumulators = 8/iter)
 * On NEON: scalar 8-way unrolled auto-vectorizes to 4 NEON FMAs (8/iter),
 *          which matches or beats explicit 2-vector SIMD (4/iter)
 * ======================================================================== */

static ST_double dot_product(const ST_double *RESTRICT x, const ST_double *RESTRICT y, ST_int N)
{
#if CTOOLS_HAS_AVX2
    return ctools_simd_dot(x, y, (size_t)N);
#else
    return ctools_dot_unrolled(x, y, N);
#endif
}

/* ========================================================================
 * Weighted dot product: sum(w[i] * x[i] * y[i])
 * Uses explicit AVX2/NEON FMA when available, otherwise falls back to
 * compiler-vectorized scalar loop.
 * ======================================================================== */

static ST_double weighted_dot_product(const ST_double *RESTRICT x, const ST_double *RESTRICT y,
                                      const ST_double *RESTRICT w, ST_int N)
{
    return ctools_simd_dot_weighted(w, x, y, (size_t)N);
}

/* Weighted projection on one continuous regressor within each category.
 * PPML represents intercept-plus-slopes as separate projection
 * factors; ordinary categorical factors retain their optimized kernels. */
static void project_slope(const ST_double *source, ST_double *dest, const FE_Factor *f, ST_int N,
                          const ST_double *weights, ST_double *means)
{
    memset(means, 0, (size_t)f->num_levels * sizeof(*means));
    for (ST_int i = 0; i < N; i++)
    {
        ST_int level = f->levels[i] - 1;
        ST_double loading = f->slope[i] - (f->slope_center ? f->slope_center[level] : 0.0);
        means[level] += source[i] * loading * (weights ? weights[i] : 1.0);
    }
    for (ST_int g = 0; g < f->num_levels; g++)
        means[g] *= weights ? f->inv_weighted_counts[g] : f->inv_counts[g];
    for (ST_int i = 0; i < N; i++)
    {
        ST_int level = f->levels[i] - 1;
        ST_double loading = f->slope[i] - (f->slope_center ? f->slope_center[level] : 0.0);
        dest[i] = source[i] - loading * means[level];
    }
}

/* ========================================================================
 * FE Projection: project-and-subtract (scatter-gather)
 * ======================================================================== */

static void project_and_subtract_fe(ST_double *RESTRICT ans, const FE_Factor *RESTRICT f, ST_int N,
                                    const ST_double *RESTRICT weights, ST_double *RESTRICT means)
{
    if (f->slope)
    {
        project_slope(ans, ans, f, N, weights, means);
        return;
    }
    ST_int i, level;
    const ST_int num_levels = f->num_levels;
    const ST_int *RESTRICT levels = f->levels;

    const ST_double *RESTRICT inv_counts = (weights != NULL && f->inv_weighted_counts != NULL)
                                               ? f->inv_weighted_counts
                                               : f->inv_counts;
    const ST_double *RESTRICT counts =
        (weights != NULL && f->weighted_counts != NULL) ? f->weighted_counts : f->counts;

    memset(means, 0, num_levels * sizeof(ST_double));

    /* Scatter: accumulate ans into means by FE level.
     * Sequential scan of ans[] is critical — hardware prefetcher handles N-sized
     * arrays perfectly. Random access to means[] is fine since it's num_levels-sized
     * and typically fits in L1/L2 cache. */
    if (weights == NULL)
    {
        for (i = 0; i < N; i++)
        {
            means[levels[i] - 1] += ans[i];
        }
    }
    else
    {
        for (i = 0; i < N; i++)
        {
            means[levels[i] - 1] += ans[i] * weights[i];
        }
    }

    if (inv_counts != NULL)
    {
#pragma omp simd
        for (level = 0; level < num_levels; level++)
        {
            means[level] *= inv_counts[level];
        }
    }
    else
    {
#pragma omp simd
        for (level = 0; level < num_levels; level++)
        {
            if (counts[level] > 0)
            {
                means[level] /= counts[level];
            }
        }
    }

    /* Gather: subtract FE means from ans */
    for (i = 0; i < N; i++)
    {
        ans[i] -= means[levels[i] - 1];
    }
}

/* ========================================================================
 * FE Projection: first projection reads from source, writes to dest
 *
 * Replaces memcpy(ans, y) + project_and_subtract_fe(ans, ...) with a single
 * function that reads from source and writes dest = source - means(source).
 * Saves one full memcpy(N) per Kaczmarz call.
 * ======================================================================== */

static void project_from_source_fe(const ST_double *RESTRICT source, ST_double *RESTRICT dest,
                                   const FE_Factor *RESTRICT f, ST_int N,
                                   const ST_double *RESTRICT weights, ST_double *RESTRICT means)
{
    if (f->slope)
    {
        project_slope(source, dest, f, N, weights, means);
        return;
    }
    ST_int i, level;
    const ST_int num_levels = f->num_levels;
    const ST_int *RESTRICT levels = f->levels;

    const ST_double *RESTRICT inv_counts = (weights != NULL && f->inv_weighted_counts != NULL)
                                               ? f->inv_weighted_counts
                                               : f->inv_counts;
    const ST_double *RESTRICT counts =
        (weights != NULL && f->weighted_counts != NULL) ? f->weighted_counts : f->counts;

    memset(means, 0, num_levels * sizeof(ST_double));

    /* Scatter: accumulate source into means by FE level */
    if (weights == NULL)
    {
        for (i = 0; i < N; i++)
        {
            means[levels[i] - 1] += source[i];
        }
    }
    else
    {
        for (i = 0; i < N; i++)
        {
            means[levels[i] - 1] += source[i] * weights[i];
        }
    }

    if (inv_counts != NULL)
    {
#pragma omp simd
        for (level = 0; level < num_levels; level++)
        {
            means[level] *= inv_counts[level];
        }
    }
    else
    {
#pragma omp simd
        for (level = 0; level < num_levels; level++)
        {
            if (counts[level] > 0)
            {
                means[level] /= counts[level];
            }
        }
    }

    /* Gather: dest = source - means[level] */
    for (i = 0; i < N; i++)
    {
        dest[i] = source[i] - means[levels[i] - 1];
    }
}

/* ========================================================================
 * FE Projection: last backward projection with fused complement
 *
 * Computes the normal scatter-gather on ans, then fuses the complement
 * (ans = source - ans) into the gather step.
 * Result: ans = source - (ans - means(ans))
 * Saves one full N-sized pass per Kaczmarz call.
 * ======================================================================== */

static void project_complement_fe(const ST_double *RESTRICT source, ST_double *RESTRICT ans,
                                  const FE_Factor *RESTRICT f, ST_int N,
                                  const ST_double *RESTRICT weights, ST_double *RESTRICT means)
{
    if (f->slope)
    {
        project_slope(ans, ans, f, N, weights, means);
        for (ST_int i = 0; i < N; i++)
            ans[i] = source[i] - ans[i];
        return;
    }
    ST_int i, level;
    const ST_int num_levels = f->num_levels;
    const ST_int *RESTRICT levels = f->levels;

    const ST_double *RESTRICT inv_counts = (weights != NULL && f->inv_weighted_counts != NULL)
                                               ? f->inv_weighted_counts
                                               : f->inv_counts;
    const ST_double *RESTRICT counts =
        (weights != NULL && f->weighted_counts != NULL) ? f->weighted_counts : f->counts;

    memset(means, 0, num_levels * sizeof(ST_double));

    /* Scatter: accumulate ans into means by FE level */
    if (weights == NULL)
    {
        for (i = 0; i < N; i++)
        {
            means[levels[i] - 1] += ans[i];
        }
    }
    else
    {
        for (i = 0; i < N; i++)
        {
            means[levels[i] - 1] += ans[i] * weights[i];
        }
    }

    if (inv_counts != NULL)
    {
#pragma omp simd
        for (level = 0; level < num_levels; level++)
        {
            means[level] *= inv_counts[level];
        }
    }
    else
    {
#pragma omp simd
        for (level = 0; level < num_levels; level++)
        {
            if (counts[level] > 0)
            {
                means[level] /= counts[level];
            }
        }
    }

    /* Fused gather + complement: ans = source - (ans - means[level]) */
    for (i = 0; i < N; i++)
    {
        ans[i] = source[i] - ans[i] + means[levels[i] - 1];
    }
}

/* ========================================================================
 * Symmetric Kaczmarz transformation with thread-local buffers
 *
 * Optimized: first forward projection reads from y directly (no memcpy),
 * last backward projection fuses the complement. Saves 2 N-sized passes.
 * ======================================================================== */

static void transform_sym_kaczmarz_threaded(const HDFE_State *RESTRICT S,
                                            const ST_double *RESTRICT y, ST_double *RESTRICT ans,
                                            ST_double **RESTRICT fe_means, ST_int thread_id)
{
    ST_int g;
    const ST_int N = S->N;
    const ST_int G = S->G;

    /* Forward sweep — first projection reads from y, writes to ans (no memcpy) */
    project_from_source_fe(y, ans, &S->factors[0], N, S->weights, fe_means[thread_id * G + 0]);

    for (g = 1; g < G; g++)
    {
        project_and_subtract_fe(ans, &S->factors[g], N, S->weights, fe_means[thread_id * G + g]);
    }

    /* Backward sweep — all but last use standard projection */
    for (g = G - 2; g > 0; g--)
    {
        project_and_subtract_fe(ans, &S->factors[g], N, S->weights, fe_means[thread_id * G + g]);
    }

    /* Last backward projection on factor 0 with complement fused:
     * ans = y - (ans - P0(ans)) = y - (I-P0)(I-P_{G-1})...(I-P0)y */
    project_complement_fe(y, ans, &S->factors[0], N, S->weights, fe_means[thread_id * G + 0]);
}

/* ========================================================================
 * Direct demean for G=1 (bypasses CG solver entirely)
 *
 * With a single FE factor, the projection is a simple demean.
 * No iteration needed — one pass produces the exact result.
 * ======================================================================== */

static void demean_column_single_fe(HDFE_State *S, ST_double *y, ST_int thread_id)
{
    const ST_int N = S->N;
    const FE_Factor *RESTRICT f = &S->factors[0];
    ST_double *RESTRICT means = S->thread_fe_means[thread_id * S->G];

    project_and_subtract_fe(y, f, N, S->weights, means);
}

/* ========================================================================
 * CG solver for a single column with thread-local buffers
 *
 * Uses fused axpy2+dot kernel to save 2 passes per iteration.
 * Convergence check avoids sqrt by comparing squared quantities.
 * ======================================================================== */

ST_int cg_solve_column_threaded(HDFE_State *S, ST_double *y, ST_int thread_id)
{
    ST_int iter;
    ST_double *RESTRICT r = S->thread_cg_r[thread_id];
    ST_double *RESTRICT u = S->thread_cg_u[thread_id];
    ST_double *RESTRICT v = S->thread_cg_v[thread_id];
    ST_double ssr, ssr_old, alpha, beta;
    ST_double improvement_potential, recent_ssr, uv;
    const ST_int N = S->N;
    const ST_double *RESTRICT weights = S->weights;
    const ST_int has_weights = S->has_weights;
    const ST_double tol_sq = S->tolerance * S->tolerance;

    if (has_weights && weights != NULL)
    {
        improvement_potential = weighted_dot_product(y, y, weights, N);

        transform_sym_kaczmarz_threaded(S, y, r, S->thread_fe_means, thread_id);
        ssr = weighted_dot_product(r, r, weights, N);
        memcpy(u, r, N * sizeof(ST_double));

        for (iter = 1; iter <= S->maxiter; iter++)
        {
            if (ctools_cancel_requested()) return -iter;
            transform_sym_kaczmarz_threaded(S, u, v, S->thread_fe_means, thread_id);

            uv = weighted_dot_product(u, v, weights, N);
            alpha = (fabs(uv) < 1e-30) ? 0.0 : ssr / uv;

            recent_ssr = alpha * ssr;
            improvement_potential -= recent_ssr;

            /* Fused: y -= alpha*u, r -= alpha*v, ssr = w·(r·r) */
            ssr_old = ssr;
            ssr = ctools_simd_fused_axpy2_wdot(y, u, r, v, weights, alpha, (size_t)N);
            beta = (fabs(ssr_old) < 1e-30) ? 0.0 : ssr / ssr_old;

            /* CG update: u = r + beta*u (SIMD-accelerated) */
            ctools_simd_axpby(u, r, 1.0, beta, (size_t)N);

            /* Convergence: recent_ssr / improvement_potential <= tol^2 */
            if (improvement_potential > 1e-30 && recent_ssr <= tol_sq * improvement_potential)
            {
                return iter;
            }
            if (improvement_potential <= 1e-30)
            {
                return iter;
            }
        }
    }
    else
    {
        improvement_potential = dot_product(y, y, N);

        transform_sym_kaczmarz_threaded(S, y, r, S->thread_fe_means, thread_id);
        ssr = dot_product(r, r, N);
        memcpy(u, r, N * sizeof(ST_double));

        for (iter = 1; iter <= S->maxiter; iter++)
        {
            if (ctools_cancel_requested()) return -iter;
            transform_sym_kaczmarz_threaded(S, u, v, S->thread_fe_means, thread_id);

            uv = dot_product(u, v, N);
            alpha = (fabs(uv) < 1e-30) ? 0.0 : ssr / uv;

            recent_ssr = alpha * ssr;
            improvement_potential -= recent_ssr;

            /* Fused: y -= alpha*u, r -= alpha*v, ssr = r·r */
            ssr_old = ssr;
            ssr = ctools_simd_fused_axpy2_dot(y, u, r, v, alpha, (size_t)N);
            beta = (fabs(ssr_old) < 1e-30) ? 0.0 : ssr / ssr_old;

            /* CG update: u = r + beta*u (SIMD-accelerated) */
            ctools_simd_axpby(u, r, 1.0, beta, (size_t)N);

            /* Convergence: recent_ssr / improvement_potential <= tol^2 */
            if (improvement_potential > 1e-30 && recent_ssr <= tol_sq * improvement_potential)
            {
                return iter;
            }
            if (improvement_potential <= 1e-30)
            {
                return iter;
            }
        }
    }

    return -S->maxiter;
}

/* ========================================================================
 * Partial Out Multiple Columns (Shared Helper)
 *
 * Partials out fixed effects from K columns of data in parallel.
 * For G=1, uses direct demean (no CG iteration needed).
 * ======================================================================== */

HDFE_SolveResult partial_out_columns(HDFE_State *S, ST_double *data, ST_int N, ST_int K,
                                     ST_int num_threads)
{
    ST_int k;
    ST_int max_iters = 0;
    ST_int any_failed = 0;

    if (!S || !data || N != S->N || K < 1 || S->num_threads < 1 || S->maxiter < 1 ||
        !isfinite(S->tolerance) || S->tolerance <= 0)
    {
        return (HDFE_SolveResult){198, 0, 0};
    }
    (void)num_threads; /* Use S->num_threads or OpenMP default */
    if (ctools_cancel_requested())
        return (HDFE_SolveResult){1, 0, 0};
    ctools_cancel_team team = {0};

    /* G=0 short-circuit: no FE factors, nothing to partial out */
    if (S->G == 0)
    {
        return (HDFE_SolveResult){0, 1, 0};
    }

    /* G=1 short-circuit: direct demean, no CG iteration */
    if (S->G == 1)
    {
#pragma omp parallel num_threads(S->num_threads)
        {
#pragma omp for schedule(dynamic) nowait
            for (k = 0; k < K; k++)
            {
                if (ctools_cancel_requested()) continue;
                int tid = 0;
#ifdef _OPENMP
                tid = omp_get_thread_num();
#endif
                demean_column_single_fe(S, data + (size_t)k * S->N, tid);
            }
            ctools_cancel_wait(&team);
        }
        if (ctools_cancel_requested()) return (HDFE_SolveResult){1, 0, 1};
        return (HDFE_SolveResult){0, 1, 1};
    }

/* General CG solve for G >= 2 */
#pragma omp parallel num_threads(S->num_threads) reduction(max : max_iters)  \
    reduction(| : any_failed)
    {
#pragma omp for schedule(dynamic) nowait
        for (k = 0; k < K; k++)
        {
            if (ctools_cancel_requested()) continue;
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            ST_int iters = cg_solve_column_threaded(S, data + (size_t)k * S->N, tid);

            if (iters < 0)
            {
                any_failed = 1;
                if (-iters > max_iters)
                    max_iters = -iters;
            }
            else
            {
                if (iters > max_iters)
                    max_iters = iters;
            }
        }
        ctools_cancel_wait(&team);
    }
    if (ctools_cancel_requested()) return (HDFE_SolveResult){1, 0, max_iters};
    return (HDFE_SolveResult){any_failed ? 430 : 0, !any_failed, max_iters};
}

/* hdfe's experimental random Kaczmarz: a symmetric sweep (forward, then
 * back) over the factors in a random order. The order is drawn once per
 * column from a xorshift generator seeded by the column (so that every sweep
 * of a column applies the same symmetric operator, as conjugate gradient
 * needs; results do not depend on the thread schedule, and Stata's
 * random-number state is not touched). */
static int hdfe_random_sym_kaczmarz(HDFE_State *S, const double *y, double *out, int tid,
                                    unsigned long long seed)
{
    double **means = S->thread_fe_means;
    int order_small[64];
    int G = S->G;
    int *order = G <= 64 ? order_small : malloc((size_t)G * sizeof(*order));
    if (!order)
        return -1;
    for (int g = 0; g < G; g++)
        order[g] = g;
    unsigned long long x = seed;
    for (int g = G - 1; g > 0; g--)
    {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        int j = (int)(x % (unsigned long long)(g + 1));
        int t = order[g];
        order[g] = order[j];
        order[j] = t;
    }
    project_from_source_fe(y, out, &S->factors[order[0]], S->N, S->weights,
                           means[tid * G + order[0]]);
    for (int g = 1; g < G; g++)
        project_and_subtract_fe(out, &S->factors[order[g]], S->N, S->weights,
                                means[tid * G + order[g]]);
    for (int g = G - 2; g >= 0; g--)
        project_and_subtract_fe(out, &S->factors[order[g]], S->N, S->weights,
                                means[tid * G + order[g]]);
    for (int i = 0; i < S->N; i++)
        out[i] = y[i] - out[i];
    if (order != order_small)
        free(order);
    return 0;
}

/* Optional hdfe transforms. out is (I-T)y; scratch may not alias either.
 * 0: symmetric Kaczmarz, 1: forward Kaczmarz, 2: Cimmino, 3: random
 * (symmetric) Kaczmarz in the order that seed draws. Returns -1 when
 * memory runs out, else 0. */
static int hdfe_map_complement(HDFE_State *S, const double *y, double *out, double *scratch,
                               int tid, int transform, unsigned long long seed)
{
    double **means = S->thread_fe_means;
    if (transform == 0)
    {
        transform_sym_kaczmarz_threaded(S, y, out, means, tid);
        return 0;
    }
    if (transform == 3)
        return hdfe_random_sym_kaczmarz(S, y, out, tid, seed);
    if (transform == 1)
    {
        project_from_source_fe(y, out, &S->factors[0], S->N, S->weights, means[tid * S->G]);
        for (int g = 1; g < S->G; g++)
            project_and_subtract_fe(out, &S->factors[g], S->N, S->weights, means[tid * S->G + g]);
        for (int i = 0; i < S->N; i++)
            out[i] = y[i] - out[i];
    }
    else
    {
        memset(out, 0, (size_t)S->N * sizeof(*out));
        for (int g = 0; g < S->G; g++)
        {
            project_from_source_fe(y, scratch, &S->factors[g], S->N, S->weights,
                                   means[tid * S->G + g]);
            for (int i = 0; i < S->N; i++)
                out[i] += (y[i] - scratch[i]) / S->G;
        }
    }
    return 0;
}

static double hdfe_map_dot(HDFE_State *S, const double *a, const double *b)
{
    return S->weights ? weighted_dot_product(a, b, S->weights, S->N) : dot_product(a, b, S->N);
}

HDFE_SolveResult partial_out_columns_options(HDFE_State *S, double *data, int N, int K,
                                             int transform, int acceleration)
{
    if (!S || !data || N != S->N || K < 1 || S->num_threads < 1 || S->maxiter < 1 ||
        !isfinite(S->tolerance) || S->tolerance <= 0 || transform < 0 || transform > 3 ||
        acceleration < 0 || acceleration > 4)
        return (HDFE_SolveResult){198, 0, 0};
    if (ctools_cancel_requested())
        return (HDFE_SolveResult){1, 0, 0};
    if (S->G == 0)
        return (HDFE_SolveResult){0, 1, 0};
    if (transform == 0 && acceleration == 1)
        return partial_out_columns(S, data, N, K, S->num_threads);
    int failed = 0, max_iterations = 0;
    const double tol_sq = S->tolerance * S->tolerance;
    ctools_cancel_team team = {0};
#pragma omp parallel num_threads(S->num_threads) reduction(| : failed)       \
    reduction(max : max_iterations)
    {
#pragma omp for schedule(dynamic) nowait
        for (int k = 0; k < K; k++)
        {
            if (ctools_cancel_requested()) continue;
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            double *y = data + (size_t)k * N;
            double *r = S->thread_cg_r[tid], *u = S->thread_cg_u[tid], *v = S->thread_cg_v[tid];
            unsigned long long seed = 0x9E3779B97F4A7C15ULL * (unsigned long long)(k + 1);
            /* Cimmino needs a distinct projection buffer. Allocate per active
             * column, never through SPI on a worker thread. */
            double *scratch = malloc((size_t)N * sizeof(*scratch));
            if (!scratch)
            {
                failed |= 2;
                continue;
            }
            double scale = hdfe_map_dot(S, y, y);
            double threshold = tol_sq * scale;
            double roundoff = 2.220446049250313e-16 * 2.220446049250313e-16 * scale;
            /* Conjugate gradient stops by the Hestenes rule of reghdfe and of
             * the default solver: the last step's decrease of the squared norm
             * is at most tol^2 times what is left to decrease (improvement) */
            double improvement = scale;
            int oom = hdfe_map_complement(S, y, r, scratch, tid, transform, seed);
            memcpy(u, r, (size_t)N * sizeof(*u));
            double rr = hdfe_map_dot(S, r, r);
            int iter = 0;
            int done = rr <= threshold;
            while (!oom && iter < S->maxiter && !done)
            {
                if (ctools_cancel_requested()) break;
                iter++;
                /* hybrid (hdfe's experimental acceleration): three plain sweeps,
                 * then conjugate gradient from the current residual */
                int accel = acceleration;
                if (acceleration == 4)
                {
                    accel = iter <= 3 ? 0 : 1;
                    if (iter == 4)
                    {
                        memcpy(u, r, (size_t)N * sizeof(*u));
                        improvement = hdfe_map_dot(S, y, y);
                    }
                }
                double step = 1;
                if (accel == 1 || accel == 2)
                {
                    oom |= hdfe_map_complement(S, accel == 1 ? u : r, v, scratch, tid, transform, seed);
                    double denom = hdfe_map_dot(S, accel == 1 ? u : r, v);
                    if (accel == 1)
                        step = fabs(denom) < 1e-30 ? 0.0 : rr / denom;
                    else
                    {
                        if (denom > 0)
                            step = rr / denom;
                        if (!isfinite(step) || step <= 0)
                            step = 1;
                    }
                }
                if (accel == 3)
                {
                    /* Vector Aitken: every third step after five plain sweeps,
                     * extrapolate along the last change using its second difference.
                     * u retains the previous iterate. */
                    if (iter >= 6 && iter % 3 == 0)
                    {
                        for (int i = 0; i < N; i++)
                            v[i] = u[i] - y[i] - r[i];
                        double denom = hdfe_map_dot(S, v, v);
                        if (denom > 0)
                            step = 1 + hdfe_map_dot(S, r, v) / denom;
                        if (!isfinite(step))
                            step = 1;
                    }
                    memcpy(u, y, (size_t)N * sizeof(*u));
                }
                for (int i = 0; i < N; i++)
                    y[i] -= step * (accel == 1 ? u[i] : r[i]);
                if (accel == 1)
                {
                    double recent = step * rr;
                    improvement -= recent;
                    /* as reghdfe (edittozerotol(, 1e-15) on its standardized
                     * data, mean square one): a negligible step has converged */
                    if (fabs(recent) < 1e-15 * scale / N)
                        recent = 0;
                    for (int i = 0; i < N; i++)
                        r[i] -= step * v[i];
                    double next = hdfe_map_dot(S, r, r), beta = rr > 0 ? next / rr : 0;
                    for (int i = 0; i < N; i++)
                        u[i] = r[i] + beta * u[i];
                    rr = next;
                    if (transform == 1)
                    {
                        /* forward Kaczmarz is not symmetric, so the recurrence's
                         * residual need not shrink: stop by reghdfe's Hestenes rule */
                        done = improvement <= 1e-30 * scale || recent <= tol_sq * improvement;
                    }
                    else
                    {
                        threshold = fmax(roundoff, tol_sq * hdfe_map_dot(S, y, y));
                        done = rr <= threshold;
                    }
                }
                else
                {
                    oom |= hdfe_map_complement(S, y, r, scratch, tid, transform, seed);
                    rr = hdfe_map_dot(S, r, r);
                    threshold = fmax(roundoff, tol_sq * hdfe_map_dot(S, y, y));
                    done = rr <= threshold;
                }
                if (!isfinite(rr))
                    break;
            }
            if (oom)
                failed |= 2;
            else if (!done || !isfinite(rr))
                failed |= 1;
            if (iter > max_iterations)
                max_iterations = iter;
            free(scratch);
        }
        ctools_cancel_wait(&team);
    }
    if (ctools_cancel_requested()) return (HDFE_SolveResult){1, 0, max_iterations};
    return (HDFE_SolveResult){(failed & 2) ? 920 : failed ? 430 : 0, !failed, max_iterations};
}

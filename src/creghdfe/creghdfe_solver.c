/*
 * creghdfe_solver.c
 *
 * CG solver, Kaczmarz transforms, FE projections
 * Highly optimized with OpenMP parallelization and SIMD hints
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_solver.h"
#include "../ctools_unroll.h"
#include "../ctools_simd.h"

/* ========================================================================
 * Dot product
 *
 * On AVX2: explicit SIMD with FMA (4 doubles/vector, 2 accumulators = 8/iter)
 * On NEON: scalar 8-way unrolled auto-vectorizes to 4 NEON FMAs (8/iter),
 *          which matches or beats explicit 2-vector SIMD (4/iter)
 * ======================================================================== */

static ST_double dot_product(const ST_double * RESTRICT x,
                             const ST_double * RESTRICT y,
                             ST_int N)
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

static ST_double weighted_dot_product(const ST_double * RESTRICT x,
                                      const ST_double * RESTRICT y,
                                      const ST_double * RESTRICT w,
                                      ST_int N)
{
    return ctools_simd_dot_weighted(w, x, y, (size_t)N);
}

/* ========================================================================
 * FE Projection: project-and-subtract (scatter-gather)
 * ======================================================================== */

static void project_and_subtract_fe(ST_double * RESTRICT ans,
                                    const FE_Factor * RESTRICT f,
                                    ST_int N,
                                    const ST_double * RESTRICT weights,
                                    ST_double * RESTRICT means)
{
    ST_int i, level;
    const ST_int num_levels = f->num_levels;
    const ST_int * RESTRICT levels = f->levels;

    const ST_double * RESTRICT inv_counts = (weights != NULL && f->inv_weighted_counts != NULL)
        ? f->inv_weighted_counts : f->inv_counts;
    const ST_double * RESTRICT counts = (weights != NULL && f->weighted_counts != NULL)
        ? f->weighted_counts : f->counts;

    memset(means, 0, num_levels * sizeof(ST_double));

    /* Scatter: accumulate ans into means by FE level.
     * Sequential scan of ans[] is critical — hardware prefetcher handles N-sized
     * arrays perfectly. Random access to means[] is fine since it's num_levels-sized
     * and typically fits in L1/L2 cache. */
    if (weights == NULL) {
        for (i = 0; i < N; i++) {
            means[levels[i] - 1] += ans[i];
        }
    } else {
        for (i = 0; i < N; i++) {
            means[levels[i] - 1] += ans[i] * weights[i];
        }
    }

    if (inv_counts != NULL) {
        #pragma omp simd
        for (level = 0; level < num_levels; level++) {
            means[level] *= inv_counts[level];
        }
    } else {
        #pragma omp simd
        for (level = 0; level < num_levels; level++) {
            if (counts[level] > 0) {
                means[level] /= counts[level];
            }
        }
    }

    /* Gather: subtract FE means from ans */
    for (i = 0; i < N; i++) {
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

static void project_from_source_fe(const ST_double * RESTRICT source,
                                    ST_double * RESTRICT dest,
                                    const FE_Factor * RESTRICT f,
                                    ST_int N,
                                    const ST_double * RESTRICT weights,
                                    ST_double * RESTRICT means)
{
    ST_int i, level;
    const ST_int num_levels = f->num_levels;
    const ST_int * RESTRICT levels = f->levels;

    const ST_double * RESTRICT inv_counts = (weights != NULL && f->inv_weighted_counts != NULL)
        ? f->inv_weighted_counts : f->inv_counts;
    const ST_double * RESTRICT counts = (weights != NULL && f->weighted_counts != NULL)
        ? f->weighted_counts : f->counts;

    memset(means, 0, num_levels * sizeof(ST_double));

    /* Scatter: accumulate source into means by FE level */
    if (weights == NULL) {
        for (i = 0; i < N; i++) {
            means[levels[i] - 1] += source[i];
        }
    } else {
        for (i = 0; i < N; i++) {
            means[levels[i] - 1] += source[i] * weights[i];
        }
    }

    if (inv_counts != NULL) {
        #pragma omp simd
        for (level = 0; level < num_levels; level++) {
            means[level] *= inv_counts[level];
        }
    } else {
        #pragma omp simd
        for (level = 0; level < num_levels; level++) {
            if (counts[level] > 0) {
                means[level] /= counts[level];
            }
        }
    }

    /* Gather: dest = source - means[level] */
    for (i = 0; i < N; i++) {
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

static void project_complement_fe(const ST_double * RESTRICT source,
                                   ST_double * RESTRICT ans,
                                   const FE_Factor * RESTRICT f,
                                   ST_int N,
                                   const ST_double * RESTRICT weights,
                                   ST_double * RESTRICT means)
{
    ST_int i, level;
    const ST_int num_levels = f->num_levels;
    const ST_int * RESTRICT levels = f->levels;

    const ST_double * RESTRICT inv_counts = (weights != NULL && f->inv_weighted_counts != NULL)
        ? f->inv_weighted_counts : f->inv_counts;
    const ST_double * RESTRICT counts = (weights != NULL && f->weighted_counts != NULL)
        ? f->weighted_counts : f->counts;

    memset(means, 0, num_levels * sizeof(ST_double));

    /* Scatter: accumulate ans into means by FE level */
    if (weights == NULL) {
        for (i = 0; i < N; i++) {
            means[levels[i] - 1] += ans[i];
        }
    } else {
        for (i = 0; i < N; i++) {
            means[levels[i] - 1] += ans[i] * weights[i];
        }
    }

    if (inv_counts != NULL) {
        #pragma omp simd
        for (level = 0; level < num_levels; level++) {
            means[level] *= inv_counts[level];
        }
    } else {
        #pragma omp simd
        for (level = 0; level < num_levels; level++) {
            if (counts[level] > 0) {
                means[level] /= counts[level];
            }
        }
    }

    /* Fused gather + complement: ans = source - (ans - means[level]) */
    for (i = 0; i < N; i++) {
        ans[i] = source[i] - ans[i] + means[levels[i] - 1];
    }
}

/* ========================================================================
 * Symmetric Kaczmarz transformation with thread-local buffers
 *
 * Optimized: first forward projection reads from y directly (no memcpy),
 * last backward projection fuses the complement. Saves 2 N-sized passes.
 * ======================================================================== */

static void transform_sym_kaczmarz_threaded(const HDFE_State * RESTRICT S,
                                            const ST_double * RESTRICT y,
                                            ST_double * RESTRICT ans,
                                            ST_double ** RESTRICT fe_means,
                                            ST_int thread_id)
{
    ST_int g;
    const ST_int N = S->N;
    const ST_int G = S->G;

    /* Forward sweep — first projection reads from y, writes to ans (no memcpy) */
    project_from_source_fe(y, ans, &S->factors[0], N, S->weights,
                           fe_means[thread_id * G + 0]);

    for (g = 1; g < G; g++) {
        project_and_subtract_fe(ans, &S->factors[g], N, S->weights,
                                fe_means[thread_id * G + g]);
    }

    /* Backward sweep — all but last use standard projection */
    for (g = G - 2; g > 0; g--) {
        project_and_subtract_fe(ans, &S->factors[g], N, S->weights,
                                fe_means[thread_id * G + g]);
    }

    /* Last backward projection on factor 0 with complement fused:
     * ans = y - (ans - P0(ans)) = y - (I-P0)(I-P_{G-1})...(I-P0)y */
    project_complement_fe(y, ans, &S->factors[0], N, S->weights,
                          fe_means[thread_id * G + 0]);
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
    const FE_Factor * RESTRICT f = &S->factors[0];
    ST_double * RESTRICT means = S->thread_fe_means[thread_id * S->G];

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
    ST_double * RESTRICT r = S->thread_cg_r[thread_id];
    ST_double * RESTRICT u = S->thread_cg_u[thread_id];
    ST_double * RESTRICT v = S->thread_cg_v[thread_id];
    ST_double ssr, ssr_old, alpha, beta;
    ST_double improvement_potential, recent_ssr, uv;
    const ST_int N = S->N;
    const ST_double * RESTRICT weights = S->weights;
    const ST_int has_weights = S->has_weights;
    const ST_double tol_sq = S->tolerance * S->tolerance;

    if (has_weights && weights != NULL) {
        improvement_potential = weighted_dot_product(y, y, weights, N);

        transform_sym_kaczmarz_threaded(S, y, r, S->thread_fe_means, thread_id);
        ssr = weighted_dot_product(r, r, weights, N);
        memcpy(u, r, N * sizeof(ST_double));

        for (iter = 1; iter <= S->maxiter; iter++) {
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
            if (improvement_potential > 1e-30 &&
                recent_ssr <= tol_sq * improvement_potential) {
                return iter;
            }
            if (improvement_potential <= 1e-30) {
                return iter;
            }
        }
    } else {
        improvement_potential = dot_product(y, y, N);

        transform_sym_kaczmarz_threaded(S, y, r, S->thread_fe_means, thread_id);
        ssr = dot_product(r, r, N);
        memcpy(u, r, N * sizeof(ST_double));

        for (iter = 1; iter <= S->maxiter; iter++) {
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
            if (improvement_potential > 1e-30 &&
                recent_ssr <= tol_sq * improvement_potential) {
                return iter;
            }
            if (improvement_potential <= 1e-30) {
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

ST_int partial_out_columns(HDFE_State *S, ST_double *data, ST_int N, ST_int K, ST_int num_threads)
{
    ST_int k;
    ST_int max_iters = 0;
    ST_int any_failed = 0;

    (void)N;  /* N is stored in S->N */
    (void)num_threads;  /* Use S->num_threads or OpenMP default */

    /* G=0 short-circuit: no FE factors, nothing to partial out */
    if (S->G == 0) {
        return 0;
    }

    /* G=1 short-circuit: direct demean, no CG iteration */
    if (S->G == 1) {
        #pragma omp parallel for num_threads(S->num_threads) schedule(dynamic)
        for (k = 0; k < K; k++) {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            demean_column_single_fe(S, data + k * S->N, tid);
        }
        return 1;  /* "1 iteration" */
    }

    /* General CG solve for G >= 2 */
    #pragma omp parallel for num_threads(S->num_threads) schedule(dynamic) \
        reduction(max:max_iters) reduction(|:any_failed)
    for (k = 0; k < K; k++) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        ST_int iters = cg_solve_column_threaded(S, data + k * S->N, tid);

        if (iters < 0) {
            any_failed = 1;
            if (-iters > max_iters) max_iters = -iters;
        } else {
            if (iters > max_iters) max_iters = iters;
        }
    }

    return any_failed ? -max_iters : max_iters;
}

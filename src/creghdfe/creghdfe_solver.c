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
 * Dot product - uses ctools K-way unrolling abstraction
 * ======================================================================== */

ST_double dot_product(const ST_double * RESTRICT x,
                      const ST_double * RESTRICT y,
                      ST_int N)
{
    return ctools_dot_unrolled(x, y, N);
}

/* ========================================================================
 * CSR Format Construction for Fast Projection
 *
 * CSR format stores observations grouped by level for cache-friendly access.
 * Instead of scatter-gather over N observations with indirect indexing,
 * we can iterate over levels and access contiguous observation indices.
 * ======================================================================== */

/*
 * Build sorted observation indices using radix LSD sort on levels.
 * This enables cache-friendly projection by accessing means sequentially.
 *
 * After sorting:
 * - sorted_indices[i] = original observation index
 * - sorted_levels[i] = level of that observation (for sequential means access)
 *
 * Observations with the same level are grouped together.
 */
int build_sorted_indices(FE_Factor *f, ST_int N)
{
    if (!f || !f->levels || f->num_levels <= 0) return -1;

    const ST_int * RESTRICT levels = f->levels;
    const ST_int num_levels = f->num_levels;

    /* Allocate output arrays */
    f->sorted_indices = (ST_int *)malloc(N * sizeof(ST_int));
    f->sorted_levels = (ST_int *)malloc(N * sizeof(ST_int));

    if (!f->sorted_indices || !f->sorted_levels) {
        if (f->sorted_indices) { free(f->sorted_indices); f->sorted_indices = NULL; }
        if (f->sorted_levels) { free(f->sorted_levels); f->sorted_levels = NULL; }
        f->sorted_initialized = 0;
        return -1;
    }

    /* Use counting sort since levels are in range [1, num_levels] */
    /* This is O(N + L) which is faster than radix sort for small L */
    ST_int *counts = (ST_int *)calloc(num_levels, sizeof(ST_int));
    ST_int *offsets = (ST_int *)malloc((num_levels + 1) * sizeof(ST_int));

    if (!counts || !offsets) {
        if (counts) free(counts);
        if (offsets) free(offsets);
        free(f->sorted_indices); f->sorted_indices = NULL;
        free(f->sorted_levels); f->sorted_levels = NULL;
        f->sorted_initialized = 0;
        return -1;
    }

    /* Count observations per level */
    for (ST_int i = 0; i < N; i++) {
        counts[levels[i] - 1]++;
    }

    /* Compute prefix sums for starting positions */
    offsets[0] = 0;
    for (ST_int l = 0; l < num_levels; l++) {
        offsets[l + 1] = offsets[l] + counts[l];
    }

    /* Reset counts to use as insertion pointers */
    memset(counts, 0, num_levels * sizeof(ST_int));

    /* Fill sorted arrays */
    for (ST_int i = 0; i < N; i++) {
        ST_int l = levels[i] - 1;  /* 0-indexed level */
        ST_int pos = offsets[l] + counts[l];
        f->sorted_indices[pos] = i;
        f->sorted_levels[pos] = levels[i];  /* Keep 1-indexed for compatibility */
        counts[l]++;
    }

    free(counts);
    free(offsets);
    f->sorted_initialized = 1;
    return 0;
}

void free_sorted_indices(FE_Factor *f)
{
    if (f) {
        if (f->sorted_indices) { free(f->sorted_indices); f->sorted_indices = NULL; }
        if (f->sorted_levels) { free(f->sorted_levels); f->sorted_levels = NULL; }
        f->sorted_initialized = 0;
    }
}

/* ========================================================================
 * Weighted dot product: sum(w[i] * x[i] * y[i])
 * Uses explicit AVX2/NEON FMA when available, otherwise falls back to
 * compiler-vectorized scalar loop.
 * ======================================================================== */

ST_double weighted_dot_product(const ST_double * RESTRICT x,
                               const ST_double * RESTRICT y,
                               const ST_double * RESTRICT w,
                               ST_int N)
{
    return ctools_simd_dot_weighted(w, x, y, (size_t)N);
}

/* ========================================================================
 * Symmetric Kaczmarz transformation with thread-local buffers
 * ======================================================================== */

/*
 * Fused project-and-subtract: computes means and subtracts in-place.
 */
static void project_and_subtract_fe(ST_double * RESTRICT ans,
                                    const FE_Factor * RESTRICT f,
                                    ST_int N,
                                    const ST_double * RESTRICT weights,
                                    ST_double * RESTRICT means)
{
    ST_int i, level;
    const ST_int num_levels = f->num_levels;
    const ST_int * RESTRICT levels = f->levels;

    /* Use precomputed inverse counts if available, otherwise fall back to counts */
    const ST_double * RESTRICT inv_counts = (weights != NULL && f->inv_weighted_counts != NULL)
        ? f->inv_weighted_counts : f->inv_counts;
    const ST_double * RESTRICT counts = (weights != NULL && f->weighted_counts != NULL)
        ? f->weighted_counts : f->counts;

    memset(means, 0, num_levels * sizeof(ST_double));

    if (weights == NULL) {
        for (i = 0; i < N; i++) {
            means[levels[i] - 1] += ans[i];
        }
    } else {
        for (i = 0; i < N; i++) {
            means[levels[i] - 1] += ans[i] * weights[i];
        }
    }

    /* Use inverse counts (multiply) if available, otherwise divide by counts */
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

    for (i = 0; i < N; i++) {
        ans[i] -= means[levels[i] - 1];
    }
}

void transform_sym_kaczmarz_threaded(const HDFE_State * RESTRICT S,
                                     const ST_double * RESTRICT y,
                                     ST_double * RESTRICT ans,
                                     ST_double ** RESTRICT fe_means,
                                     ST_int thread_id)
{
    ST_int g, i;
    const ST_int N = S->N;
    const ST_int G = S->G;

    memcpy(ans, y, N * sizeof(ST_double));

    /* Forward sweep */
    for (g = 0; g < G; g++) {
        project_and_subtract_fe(ans, &S->factors[g], N, S->weights,
                                fe_means[thread_id * G + g]);
    }

    /* Backward sweep */
    for (g = G - 2; g >= 0; g--) {
        project_and_subtract_fe(ans, &S->factors[g], N, S->weights,
                                fe_means[thread_id * G + g]);
    }

    /* Return projection: y - ans */
    #pragma omp simd
    for (i = 0; i < N; i++) {
        ans[i] = y[i] - ans[i];
    }
}

/* ========================================================================
 * CG solver for a single column with thread-local buffers
 * ======================================================================== */

ST_int cg_solve_column_threaded(HDFE_State *S, ST_double *y, ST_int thread_id)
{
    ST_int iter;
    ST_double * RESTRICT r = S->thread_cg_r[thread_id];
    ST_double * RESTRICT u = S->thread_cg_u[thread_id];
    ST_double * RESTRICT v = S->thread_cg_v[thread_id];
    ST_double ssr, ssr_old, alpha, beta;
    ST_double improvement_potential, recent_ssr, update_error, uv;
    const ST_int N = S->N;
    const ST_double * RESTRICT weights = S->weights;
    const ST_int has_weights = S->has_weights;

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

            /* CG update: y -= alpha*u, r -= alpha*v (SIMD-accelerated) */
            ctools_simd_axpy(y, u, -alpha, (size_t)N);
            ctools_simd_axpy(r, v, -alpha, (size_t)N);

            ssr_old = ssr;
            ssr = weighted_dot_product(r, r, weights, N);
            beta = (fabs(ssr_old) < 1e-30) ? 0.0 : ssr / ssr_old;

            /* CG update: u = r + beta*u (SIMD-accelerated) */
            ctools_simd_axpby(u, r, 1.0, beta, (size_t)N);

            update_error = (improvement_potential > 1e-30) ? sqrt(recent_ssr / improvement_potential) : 0.0;
            if (update_error <= S->tolerance) {
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

            /* CG update: y -= alpha*u, r -= alpha*v (SIMD-accelerated) */
            ctools_simd_axpy(y, u, -alpha, (size_t)N);
            ctools_simd_axpy(r, v, -alpha, (size_t)N);

            ssr_old = ssr;
            ssr = dot_product(r, r, N);
            beta = (fabs(ssr_old) < 1e-30) ? 0.0 : ssr / ssr_old;

            /* CG update: u = r + beta*u (SIMD-accelerated) */
            ctools_simd_axpby(u, r, 1.0, beta, (size_t)N);

            update_error = (improvement_potential > 1e-30) ? sqrt(recent_ssr / improvement_potential) : 0.0;
            if (update_error <= S->tolerance) {
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
 * ======================================================================== */

ST_int partial_out_columns(HDFE_State *S, ST_double *data, ST_int N, ST_int K, ST_int num_threads)
{
    ST_int k;
    ST_int max_iters = 0;
    ST_int any_failed = 0;

    (void)N;  /* N is stored in S->N */
    (void)num_threads;  /* Use S->num_threads or OpenMP default */

    #pragma omp parallel for num_threads(S->num_threads) schedule(dynamic)
    for (k = 0; k < K; k++) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        ST_int iters = cg_solve_column_threaded(S, data + k * S->N, tid);

        #pragma omp critical
        {
            if (iters < 0) {
                any_failed = 1;
            }
            if (iters > max_iters) {
                max_iters = iters;
            }
        }
    }

    return any_failed ? -max_iters : max_iters;
}

/*
 * creghdfe_solver.c
 *
 * CG solver, Kaczmarz transforms, FE projections
 * Highly optimized with OpenMP parallelization and SIMD hints
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_solver.h"
#include "../ctools_unroll.h"

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
 * Project onto a single fixed effect using thread-local means buffer
 * ======================================================================== */

void project_one_fe_threaded(const ST_double * RESTRICT y,
                             const FE_Factor * RESTRICT f,
                             ST_int N,
                             const ST_double * RESTRICT weights,
                             ST_double * RESTRICT proj,
                             ST_double * RESTRICT means)
{
    ST_int i, level;
    const ST_int num_levels = f->num_levels;
    const ST_int * RESTRICT levels = f->levels;

    /* Use weighted_counts when weights are present, otherwise use unweighted counts */
    const ST_double * RESTRICT counts = (weights != NULL && f->weighted_counts != NULL)
        ? f->weighted_counts : f->counts;

    memset(means, 0, num_levels * sizeof(ST_double));

    if (weights == NULL) {
        /* Unweighted: accumulate y values */
        for (i = 0; i < N; i++) {
            means[levels[i] - 1] += y[i];
        }
    } else {
        /* Weighted: accumulate y * weight values */
        for (i = 0; i < N; i++) {
            means[levels[i] - 1] += y[i] * weights[i];
        }
    }

    /* Divide by counts (weighted or unweighted as appropriate) */
    #pragma omp simd
    for (level = 0; level < num_levels; level++) {
        if (counts[level] > 0) {
            means[level] /= counts[level];
        }
    }

    /* Broadcast means to projection (gather operation) */
    for (i = 0; i < N; i++) {
        proj[i] = means[levels[i] - 1];
    }
}

/* ========================================================================
 * Symmetric Kaczmarz transformation with thread-local buffers
 * ======================================================================== */

void transform_sym_kaczmarz_threaded(const HDFE_State * RESTRICT S,
                                     const ST_double * RESTRICT y,
                                     ST_double * RESTRICT ans,
                                     ST_double * RESTRICT proj,
                                     ST_double ** RESTRICT fe_means,
                                     ST_int thread_id)
{
    ST_int g, i;
    const ST_int N = S->N;
    const ST_int G = S->G;

    memcpy(ans, y, N * sizeof(ST_double));

    /* Forward sweep */
    for (g = 0; g < G; g++) {
        project_one_fe_threaded(ans, &S->factors[g], N, S->weights, proj,
                                fe_means[thread_id * G + g]);
        #pragma omp simd
        for (i = 0; i < N; i++) {
            ans[i] -= proj[i];
        }
    }

    /* Backward sweep */
    for (g = G - 2; g >= 0; g--) {
        project_one_fe_threaded(ans, &S->factors[g], N, S->weights, proj,
                                fe_means[thread_id * G + g]);
        #pragma omp simd
        for (i = 0; i < N; i++) {
            ans[i] -= proj[i];
        }
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
    ST_int iter, i;
    ST_double * RESTRICT r = S->thread_cg_r[thread_id];
    ST_double * RESTRICT u = S->thread_cg_u[thread_id];
    ST_double * RESTRICT v = S->thread_cg_v[thread_id];
    ST_double * RESTRICT proj = S->thread_proj[thread_id];
    ST_double ssr, ssr_old, alpha, beta;
    ST_double improvement_potential, recent_ssr, update_error, uv;
    const ST_int N = S->N;

    improvement_potential = dot_product(y, y, N);

    transform_sym_kaczmarz_threaded(S, y, r, proj, S->thread_fe_means, thread_id);
    ssr = dot_product(r, r, N);
    memcpy(u, r, N * sizeof(ST_double));

    for (iter = 1; iter <= S->maxiter; iter++) {
        transform_sym_kaczmarz_threaded(S, u, v, proj, S->thread_fe_means, thread_id);

        uv = dot_product(u, v, N);
        alpha = (fabs(uv) < 1e-30) ? 0.0 : ssr / uv;

        recent_ssr = alpha * ssr;
        improvement_potential -= recent_ssr;

        /* Fused update loop with SIMD */
        #pragma omp simd
        for (i = 0; i < N; i++) {
            y[i] -= alpha * u[i];
            r[i] -= alpha * v[i];
        }

        ssr_old = ssr;
        ssr = dot_product(r, r, N);
        beta = (fabs(ssr_old) < 1e-30) ? 0.0 : ssr / ssr_old;

        #pragma omp simd
        for (i = 0; i < N; i++) {
            u[i] = r[i] + beta * u[i];
        }

        update_error = (improvement_potential > 1e-30) ? sqrt(recent_ssr / improvement_potential) : 0.0;
        if (update_error <= S->tolerance) {
            return iter;
        }
    }

    return -S->maxiter;
}

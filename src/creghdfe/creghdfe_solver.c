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
 * CSR Format Construction for Fast Projection
 *
 * CSR format stores observations grouped by level for cache-friendly access.
 * Instead of scatter-gather over N observations with indirect indexing,
 * we can iterate over levels and access contiguous observation indices.
 * ======================================================================== */

int build_csr_format(FE_Factor *f, ST_int N)
{
    if (!f || !f->levels || f->num_levels <= 0) return -1;

    const ST_int num_levels = f->num_levels;
    const ST_int * RESTRICT levels = f->levels;

    /* Allocate offsets and indices */
    f->csr_offsets = (ST_int *)malloc((num_levels + 1) * sizeof(ST_int));
    f->csr_indices = (ST_int *)malloc(N * sizeof(ST_int));

    if (!f->csr_offsets || !f->csr_indices) {
        if (f->csr_offsets) { free(f->csr_offsets); f->csr_offsets = NULL; }
        if (f->csr_indices) { free(f->csr_indices); f->csr_indices = NULL; }
        f->csr_initialized = 0;
        return -1;
    }

    /* Count observations per level (use existing counts as double, convert to int) */
    ST_int *counts = (ST_int *)calloc(num_levels, sizeof(ST_int));
    if (!counts) {
        free(f->csr_offsets); f->csr_offsets = NULL;
        free(f->csr_indices); f->csr_indices = NULL;
        f->csr_initialized = 0;
        return -1;
    }

    /* Count observations per level */
    for (ST_int i = 0; i < N; i++) {
        counts[levels[i] - 1]++;  /* levels are 1-indexed */
    }

    /* Compute prefix sums for offsets */
    f->csr_offsets[0] = 0;
    for (ST_int l = 0; l < num_levels; l++) {
        f->csr_offsets[l + 1] = f->csr_offsets[l] + counts[l];
    }

    /* Reset counts to use as insertion pointers */
    memset(counts, 0, num_levels * sizeof(ST_int));

    /* Fill indices array: group observation indices by level */
    for (ST_int i = 0; i < N; i++) {
        ST_int l = levels[i] - 1;  /* 0-indexed level */
        ST_int pos = f->csr_offsets[l] + counts[l];
        f->csr_indices[pos] = i;
        counts[l]++;
    }

    free(counts);
    f->csr_initialized = 1;
    return 0;
}

void free_csr_format(FE_Factor *f)
{
    if (f) {
        if (f->csr_offsets) { free(f->csr_offsets); f->csr_offsets = NULL; }
        if (f->csr_indices) { free(f->csr_indices); f->csr_indices = NULL; }
        f->csr_initialized = 0;
    }
}

/* ========================================================================
 * CSR-accelerated FE projection
 *
 * Instead of:
 *   for i in 0..N:
 *     means[levels[i]-1] += y[i]  // scatter (poor locality)
 *     proj[i] = means[levels[i]-1]  // gather (poor locality)
 *
 * We do:
 *   for level in 0..num_levels:
 *     for pos in offsets[level]..offsets[level+1]:
 *       means[level] += y[indices[pos]]  // still gather, but we cache level
 *     for pos in offsets[level]..offsets[level+1]:
 *       proj[indices[pos]] = means[level]  // scatter with cached mean
 *
 * Benefits:
 * - Level mean is in register during inner loop
 * - Access patterns are more predictable for prefetcher
 * - Can parallelize over levels if num_levels is large
 * ======================================================================== */

void project_one_fe_csr(const ST_double * RESTRICT y,
                        const FE_Factor * RESTRICT f,
                        ST_int N,
                        const ST_double * RESTRICT weights,
                        ST_double * RESTRICT proj,
                        ST_double * RESTRICT means)
{
    const ST_int num_levels = f->num_levels;
    const ST_int * RESTRICT offsets = f->csr_offsets;
    const ST_int * RESTRICT indices = f->csr_indices;
    const ST_double * RESTRICT counts = (weights != NULL && f->weighted_counts != NULL)
        ? f->weighted_counts : f->counts;

    (void)N;  /* N is implicit in offsets */

    if (weights == NULL) {
        /* Unweighted case - parallelize over levels for large FEs */
        #ifdef _OPENMP
        #pragma omp parallel for if(num_levels > 1000) schedule(static)
        #endif
        for (ST_int level = 0; level < num_levels; level++) {
            ST_double sum = 0.0;
            const ST_int start = offsets[level];
            const ST_int end = offsets[level + 1];

            /* Accumulate y values for this level */
            for (ST_int pos = start; pos < end; pos++) {
                sum += y[indices[pos]];
            }

            /* Compute mean */
            ST_double mean = (counts[level] > 0) ? sum / counts[level] : 0.0;
            means[level] = mean;

            /* Broadcast mean to proj array */
            for (ST_int pos = start; pos < end; pos++) {
                proj[indices[pos]] = mean;
            }
        }
    } else {
        /* Weighted case */
        #ifdef _OPENMP
        #pragma omp parallel for if(num_levels > 1000) schedule(static)
        #endif
        for (ST_int level = 0; level < num_levels; level++) {
            ST_double sum = 0.0;
            const ST_int start = offsets[level];
            const ST_int end = offsets[level + 1];

            /* Accumulate weighted y values for this level */
            for (ST_int pos = start; pos < end; pos++) {
                ST_int i = indices[pos];
                sum += y[i] * weights[i];
            }

            /* Compute mean */
            ST_double mean = (counts[level] > 0) ? sum / counts[level] : 0.0;
            means[level] = mean;

            /* Broadcast mean to proj array */
            for (ST_int pos = start; pos < end; pos++) {
                proj[indices[pos]] = mean;
            }
        }
    }
}

/* ========================================================================
 * Weighted dot product: sum(w[i] * x[i] * y[i])
 * ======================================================================== */

ST_double weighted_dot_product(const ST_double * RESTRICT x,
                               const ST_double * RESTRICT y,
                               const ST_double * RESTRICT w,
                               ST_int N)
{
    ST_double sum = 0.0;
    ST_int i;

    #pragma omp simd reduction(+:sum)
    for (i = 0; i < N; i++) {
        sum += w[i] * x[i] * y[i];
    }
    return sum;
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

    /* Forward sweep - use CSR projection when available */
    for (g = 0; g < G; g++) {
        const FE_Factor *f = &S->factors[g];
        if (f->csr_initialized) {
            project_one_fe_csr(ans, f, N, S->weights, proj,
                               fe_means[thread_id * G + g]);
        } else {
            project_one_fe_threaded(ans, f, N, S->weights, proj,
                                    fe_means[thread_id * G + g]);
        }
        #pragma omp simd
        for (i = 0; i < N; i++) {
            ans[i] -= proj[i];
        }
    }

    /* Backward sweep - use CSR projection when available */
    for (g = G - 2; g >= 0; g--) {
        const FE_Factor *f = &S->factors[g];
        if (f->csr_initialized) {
            project_one_fe_csr(ans, f, N, S->weights, proj,
                               fe_means[thread_id * G + g]);
        } else {
            project_one_fe_threaded(ans, f, N, S->weights, proj,
                                    fe_means[thread_id * G + g]);
        }
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
    const ST_double * RESTRICT weights = S->weights;
    const ST_int has_weights = S->has_weights;

    /* Use weighted dot products when weights are present */
    if (has_weights && weights != NULL) {
        improvement_potential = weighted_dot_product(y, y, weights, N);

        transform_sym_kaczmarz_threaded(S, y, r, proj, S->thread_fe_means, thread_id);
        ssr = weighted_dot_product(r, r, weights, N);
        memcpy(u, r, N * sizeof(ST_double));

        for (iter = 1; iter <= S->maxiter; iter++) {
            transform_sym_kaczmarz_threaded(S, u, v, proj, S->thread_fe_means, thread_id);

            uv = weighted_dot_product(u, v, weights, N);
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
            ssr = weighted_dot_product(r, r, weights, N);
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
    } else {
        /* Unweighted case - original code */
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
    }

    return -S->maxiter;
}

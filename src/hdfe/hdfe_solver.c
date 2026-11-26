/*
 * hdfe_solver.c
 *
 * CG solver, Kaczmarz transforms, FE projections
 * Highly optimized with OpenMP parallelization and SIMD hints
 * Part of the ctools Stata plugin suite
 */

#include "hdfe_solver.h"
#include "hdfe_coef.h"
#include "../ctools_unroll.h"
#include "../ctools_simd.h"
#include <float.h>

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

/* A tiny denominator alone is not convergence. Recursive residuals can reach
 * the arithmetic floor, however, especially for an already partialled column.
 * Accept that floor only with well-scaled weights. Tiny positive weights can
 * hide a large unresolved FE component behind a small weighted residual; that
 * ambiguous case must fail. The weight scan is cold and replaces the old
 * zero-step update of four observation vectors, never adding a normal pass.
 *
 * The exception is S->floor_any_weight_span (cppmlhdfe). There the lightest
 * rows are zero outcomes whose IRLS weight, their mean, heads to zero: the
 * norm the projection is wanted in, not weights a user chose. Their range
 * passes 1/(128 eps) after enough iterations on such panels, and the floor
 * the loop reached with well-scaled weights one iteration before (and that
 * ppmlhdfe accepts) then became r(498). cppmlhdfe still leaves the flag
 * unset when its user's weights alone span that range. */
int hdfe_weight_span_too_wide(const ST_double *weights, ST_int N)
{
    ST_double smallest = DBL_MAX, largest = 0.0;
    for (ST_int i = 0; i < N; i++) {
        ST_double w = weights[i];
        if (w > largest) largest = w;
        if (w > 0.0 && w < smallest) smallest = w;
    }
    return smallest <= 128.0 * DBL_EPSILON * largest;
}

static ST_int cg_breakdown_result(const HDFE_State *S, ST_double residual_norm2,
                                  ST_double initial_norm2, ST_int iter)
{
    const ST_double roundoff = 64.0 * DBL_EPSILON * DBL_EPSILON;
    if (!isfinite(initial_norm2) || !isfinite(residual_norm2) ||
        residual_norm2 > roundoff * initial_norm2)
        return HDFE_CG_BREAKDOWN;
    if (S->has_weights && S->weights && !S->floor_any_weight_span &&
        hdfe_weight_span_too_wide(S->weights, S->N))
        return HDFE_CG_BREAKDOWN;
    return iter;
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

/* Level sums in means to level means, as the projections above do */
static inline void fe_means_finish(const FE_Factor *RESTRICT f, const ST_double *RESTRICT weights,
                                   ST_double *RESTRICT means)
{
    ST_int level;
    const ST_int num_levels = f->num_levels;
    const ST_double *RESTRICT inv_counts = (weights != NULL && f->inv_weighted_counts != NULL)
                                               ? f->inv_weighted_counts
                                               : f->inv_counts;
    const ST_double *RESTRICT counts =
        (weights != NULL && f->weighted_counts != NULL) ? f->weighted_counts : f->counts;
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
}

/* The sweep of transform_sym_kaczmarz_threaded for G >= 2 categorical
 * factors with each projection's gather fused into the next projection's
 * scatter: the value subtracted and stored is the value the next factor's
 * level sums accumulate, in the same row order, so every sum and every
 * stored value equals the unfused sweep's. 2G passes over the rows instead
 * of 4G - 2. */
static void transform_sym_kaczmarz_fused(const HDFE_State *RESTRICT S, const ST_double *RESTRICT y,
                                         ST_double *RESTRICT ans, ST_double **RESTRICT fe_means,
                                         ST_int thread_id)
{
    ST_int i;
    const ST_int N = S->N;
    const ST_int G = S->G;
    const ST_double *RESTRICT weights = S->weights;
    const FE_Factor *prev = &S->factors[0];
    const ST_int *prev_levels = prev->levels;
    ST_double *prev_means = fe_means[thread_id * G + 0];

    /* Level means of y over the first factor */
    memset(prev_means, 0, prev->num_levels * sizeof(ST_double));
    if (weights == NULL)
    {
        for (i = 0; i < N; i++)
        {
            prev_means[prev_levels[i] - 1] += y[i];
        }
    }
    else
    {
        for (i = 0; i < N; i++)
        {
            prev_means[prev_levels[i] - 1] += y[i] * weights[i];
        }
    }
    fe_means_finish(prev, weights, prev_means);

    /* Factors 1 .. G-1, then G-2 .. 0: subtract the previous factor's means
     * and accumulate this factor's level sums in one pass */
    for (ST_int step = 1; step < 2 * G - 1; step++)
    {
        const ST_int g = step < G ? step : 2 * G - 2 - step;
        const FE_Factor *cur = &S->factors[g];
        const ST_int *cur_levels = cur->levels;
        ST_double *cur_means = fe_means[thread_id * G + g];
        memset(cur_means, 0, cur->num_levels * sizeof(ST_double));
        if (step == 1)
        {
            if (weights == NULL)
            {
                for (i = 0; i < N; i++)
                {
                    const ST_double value = y[i] - prev_means[prev_levels[i] - 1];
                    ans[i] = value;
                    cur_means[cur_levels[i] - 1] += value;
                }
            }
            else
            {
                for (i = 0; i < N; i++)
                {
                    const ST_double value = y[i] - prev_means[prev_levels[i] - 1];
                    ans[i] = value;
                    cur_means[cur_levels[i] - 1] += value * weights[i];
                }
            }
        }
        else if (weights == NULL)
        {
            for (i = 0; i < N; i++)
            {
                const ST_double value = ans[i] - prev_means[prev_levels[i] - 1];
                ans[i] = value;
                cur_means[cur_levels[i] - 1] += value;
            }
        }
        else
        {
            for (i = 0; i < N; i++)
            {
                const ST_double value = ans[i] - prev_means[prev_levels[i] - 1];
                ans[i] = value;
                cur_means[cur_levels[i] - 1] += value * weights[i];
            }
        }
        fe_means_finish(cur, weights, cur_means);
        prev = cur;
        prev_levels = cur_levels;
        prev_means = cur_means;
    }

    /* Last gather (factor 0) with the complement:
     * ans = y - (ans - P0(ans)) */
    for (i = 0; i < N; i++)
    {
        ans[i] = y[i] - ans[i] + prev_means[prev_levels[i] - 1];
    }
}

static void transform_sym_kaczmarz_threaded(const HDFE_State *RESTRICT S,
                                            const ST_double *RESTRICT y, ST_double *RESTRICT ans,
                                            ST_double **RESTRICT fe_means, ST_int thread_id)
{
    ST_int g;
    const ST_int N = S->N;
    const ST_int G = S->G;

    if (G >= 2)
    {
        ST_int categorical = 1;
        for (g = 0; g < G; g++)
            if (S->factors[g].slope) categorical = 0;
        if (categorical)
        {
            transform_sym_kaczmarz_fused(S, y, ans, fe_means, thread_id);
            return;
        }
    }

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
 * One effect whose rows arrive grouped
 *
 * When the levels never decrease along the rows, the rows of a level are
 * consecutive (a panel sorted by its unit, say): the level's sum is complete
 * where its run ends, while those rows are still in the cache. Each run is
 * summed and its mean subtracted at once, so a column comes from memory once
 * where project_and_subtract_fe() reads it twice (the scatter into the level
 * sums, then the gather), and no add waits for the store of the row before
 * it, which is what holds the scatter back when equal levels follow one
 * another. A level's sum adds the same values in the same order as the
 * scatter does: the residuals are the same to the last bit.
 *
 * The rows are worked in blocks that hold whole runs, any number at a time:
 * nothing is shared between runs, so neither the blocks nor the thread count
 * can change a result.
 * ======================================================================== */

/* Rows per block, and columns walked together through one */
#define HDFE_GROUPED_BLOCK 32768
#define HDFE_GROUPED_COLUMNS 8

/* 1 when levels[] never decreases */
static int fe_levels_grouped(const ST_int *RESTRICT levels, ST_int N)
{
    /* rows in no particular order fail within the first few */
    const ST_int probe = N < 4096 ? N : 4096;
    for (ST_int i = 1; i < probe; i++)
    {
        if (levels[i] < levels[i - 1]) return 0;
    }
    int decreases = 0;
#pragma omp parallel for schedule(static) reduction(| : decreases)                                 \
    if (N > 200000 && omp_get_level() == 0)
    for (ST_int i = probe; i < N; i++)
    {
        decreases |= levels[i] < levels[i - 1];
    }
    return !decreases;
}

/* The rows i0 .. i1-1, which hold whole runs, of the K columns x, x + N, ...:
 * every run less its (weighted) level mean, as project_and_subtract_fe() */
static void demean_grouped_rows(ST_double *x, size_t N, ST_int K, const FE_Factor *RESTRICT f,
                                const ST_double *RESTRICT weights, ST_int i0, ST_int i1)
{
    const ST_int *RESTRICT levels = f->levels;
    const ST_double *RESTRICT inv_counts = (weights != NULL && f->inv_weighted_counts != NULL)
                                               ? f->inv_weighted_counts
                                               : f->inv_counts;
    const ST_double *RESTRICT counts =
        (weights != NULL && f->weighted_counts != NULL) ? f->weighted_counts : f->counts;
    ST_int i = i0;
    while (i < i1)
    {
        const ST_int level = levels[i] - 1;
        ST_int j = i + 1;
        while (j < i1 && levels[j] - 1 == level) j++;
        const ST_int n = j - i;
        for (ST_int k = 0; k < K; k++)
        {
            ST_double *RESTRICT col = x + (size_t)k * N + i;
            ST_double mean = 0.0;
            if (weights == NULL)
            {
                for (ST_int r = 0; r < n; r++)
                {
                    mean += col[r];
                }
            }
            else
            {
                const ST_double *RESTRICT w = weights + i;
                for (ST_int r = 0; r < n; r++)
                {
                    mean += col[r] * w[r];
                }
            }
            if (inv_counts != NULL)
            {
                mean *= inv_counts[level];
            }
            else if (counts[level] > 0)
            {
                mean /= counts[level];
            }
            for (ST_int r = 0; r < n; r++)
            {
                col[r] -= mean;
            }
        }
        i = j;
    }
}

/* The projection on one effect whose rows are grouped, for K columns.
 * Returns 0 when it does not apply (a slope, or levels that decrease
 * somewhere), and the caller projects the general way. */
static int demean_columns_grouped(const HDFE_State *S, ST_double *data, ST_int K)
{
    const ST_int N = S->N;
    const FE_Factor *f = &S->factors[0];
    if (f->slope || N < 2 || !fe_levels_grouped(f->levels, N)) return 0;

    /* Block b takes the rows start[b] .. start[b+1]-1: a run that crosses a
     * nominal boundary stays with the block it began in. A run that runs on
     * for more than a block past one (few levels, many rows each) would
     * leave whole blocks to a single thread: then, as without the memory for
     * the boundaries, each column is one task instead. */
    const ST_int *levels = f->levels;
    const ST_int nblk = (N + HDFE_GROUPED_BLOCK - 1) / HDFE_GROUPED_BLOCK;
    ST_int *start = (ST_int *)malloc(((size_t)nblk + 1) * sizeof(ST_int));
    int by_block = start != NULL && nblk > 1;
    if (by_block)
    {
        start[0] = 0;
        for (ST_int b = 1; b < nblk; b++)
        {
            ST_int i = b * HDFE_GROUPED_BLOCK;
            const ST_int limit = N - i > HDFE_GROUPED_BLOCK ? i + HDFE_GROUPED_BLOCK : N;
            while (i < limit && levels[i] == levels[i - 1]) i++;
            if (i == limit && i < N && levels[i] == levels[i - 1])
            {
                by_block = 0;
                break;
            }
            start[b] = i;
        }
        start[nblk] = N;
    }
    if (by_block)
    {
#pragma omp parallel for schedule(dynamic, 1)                                                      \
    if ((size_t)N * (size_t)K > 200000 && omp_get_level() == 0)
        for (ST_int b = 0; b < nblk; b++)
        {
            for (ST_int k0 = 0; k0 < K; k0 += HDFE_GROUPED_COLUMNS)
            {
                const ST_int kn = K - k0 < HDFE_GROUPED_COLUMNS ? K - k0 : HDFE_GROUPED_COLUMNS;
                demean_grouped_rows(data + (size_t)k0 * N, (size_t)N, kn, f, S->weights, start[b],
                                    start[b + 1]);
            }
        }
    }
    else
    {
#pragma omp parallel for schedule(dynamic, 1)                                                      \
    if (K > 1 && (size_t)N * (size_t)K > 200000 && omp_get_level() == 0)
        for (ST_int k = 0; k < K; k++)
        {
            demean_grouped_rows(data + (size_t)k * N, (size_t)N, 1, f, S->weights, 0, N);
        }
    }
    free(start);
    return 1;
}

/* ========================================================================
 * Rounding floor of the CG solver
 *
 * The first sweep and the first step of a solve are rounded at the size of
 * the column it starts from, and the operator is singular: the part of that
 * rounding outside its range stays in the residual for good. It measures
 * 0.25 to 0.6 units of 2.2e-16 of the column's root sum of squares (random,
 * balanced and weakly connected designs, weighted or not, column means from
 * 40 to 1e6 standard deviations). Once the residual has come down to it, the
 * step length r'r / u'Au divides noise by noise, and a few such steps destroy
 * the column: x = 2000 + e at tolerance 1e-13 ended almost a million times
 * the size of its projection away from it, with convergence reported.
 *
 * So a round ends when the residual sum of squares is within
 * HDFE_CG_ROUNDING^2 of the sum of squares the round started from, and no
 * step is taken from a residual that small. If the column has by then lost
 * most of its sum of squares, the next round starts from what is left of
 * it, whose rounding is that much smaller; otherwise double precision has no
 * more to give. Three rounds are what a column needs that the effects
 * explain entirely: the column, its remainder, the rounding of that.
 *
 * 32 units keep the last step three orders of magnitude (in the sum of
 * squares) above the noise. A wider margin ends rounds that had digits left
 * to gain: with 1024 units, columns whose mean is a few thousand standard
 * deviations started a second round at the default tolerance, and ordinary
 * columns lost one to two digits at tolerances of 1e-12 and below.
 * ======================================================================== */

#define HDFE_CG_ROUNDING (32.0 * 2.220446049250313e-16)
#define HDFE_CG_ROUNDS 3

/* Residual sum of squares at which a round that started from a column with
 * sum of squares scale ends; an overflowed scale sets no floor */
static inline ST_double cg_rounding_floor(ST_double scale)
{
    return isfinite(scale) ? HDFE_CG_ROUNDING * HDFE_CG_ROUNDING * scale : -1.0;
}

/* A round has ended: its residual is at the floor, or its direction has no
 * curvature (recent, the decrease the step not taken would have made, is
 * then negative). Returns 1 when the column is finished:
 *   - it has had its rounds, or it still has a good part of the sum of
 *     squares the round started from, so that starting again gains nothing;
 *   - or the tolerance asks for no more: neither the step not taken nor the
 *     rounding of this round is beyond what it allows of the column.
 * Otherwise the next round starts from the column as it now is: residual and
 * direction are recomputed from it, and the potential, the scale and the
 * floor are those of the new start; if that residual is at its floor
 * already, the column is finished too. */
static int cg_round_end(HDFE_State *S, ST_double *y, ST_int thread_id, ST_double recent,
                        ST_int *round, ST_double *potential, ST_double *scale,
                        ST_double *rounding, ST_double *ssr)
{
    const ST_int N = S->N;
    const ST_double *weights = (S->has_weights && S->weights != NULL) ? S->weights : NULL;
    ST_double *r = S->thread_cg_r[thread_id];
    const ST_double left =
        weights ? weighted_dot_product(y, y, weights, N) : dot_product(y, y, N);
    const ST_double allowed = S->tolerance * S->tolerance * left;

    if (*round >= HDFE_CG_ROUNDS || 4.0 * left >= *scale) return 1;
    if (recent >= 0.0 && recent <= allowed && *rounding <= allowed) return 1;
    (*round)++;
    *potential = *scale = left;
    *rounding = cg_rounding_floor(left);
    transform_sym_kaczmarz_threaded(S, y, r, S->thread_fe_means, thread_id);
    *ssr = weights ? weighted_dot_product(r, r, weights, N) : dot_product(r, r, N);
    memcpy(S->thread_cg_u[thread_id], r, (size_t)N * sizeof(ST_double));
    return *ssr <= *rounding;
}

/* ========================================================================
 * CG solver for a single column with thread-local buffers
 *
 * Uses fused axpy2+dot kernel to save 2 passes per iteration.
 * Convergence check avoids sqrt by comparing squared quantities.
 * The rounding floor is only consulted before a step: a solve that stays
 * above it is the plain Hestenes-Stiefel iteration.
 * ======================================================================== */

ST_int cg_solve_column_threaded(HDFE_State *S, ST_double *y, ST_int thread_id)
{
    ST_int iter;
    ST_double *RESTRICT r = S->thread_cg_r[thread_id];
    ST_double *RESTRICT u = S->thread_cg_u[thread_id];
    ST_double *RESTRICT v = S->thread_cg_v[thread_id];
    ST_double ssr, ssr_old, alpha, beta;
    ST_double improvement_potential, initial_norm2, recent_ssr, uv;
    const ST_int N = S->N;
    const ST_double *RESTRICT weights = S->weights;
    const ST_int has_weights = S->has_weights;
    const ST_double tol_sq = S->tolerance * S->tolerance;
    ST_double scale, rounding;
    ST_int round = 1;

    if (has_weights && weights != NULL)
    {
        improvement_potential = initial_norm2 = weighted_dot_product(y, y, weights, N);

        transform_sym_kaczmarz_threaded(S, y, r, S->thread_fe_means, thread_id);
        ssr = weighted_dot_product(r, r, weights, N);
        memcpy(u, r, N * sizeof(ST_double));
        scale = improvement_potential;
        rounding = cg_rounding_floor(scale);

        for (iter = 1; iter <= S->maxiter; iter++)
        {
            if (ctools_cancel_requested()) return -iter;
            if (ssr == 0.0) return iter;
            transform_sym_kaczmarz_threaded(S, u, v, S->thread_fe_means, thread_id);

            uv = weighted_dot_product(u, v, weights, N);
            /* I - symmetric MAP is a positive semidefinite contraction.
             * A direction whose curvature is at roundoff relative to r'r
             * cannot certify the projection. Never turn breakdown into a
             * zero step: zero decrease would pass the convergence test below.
             * This guard is invariant to column and common weight scaling;
             * it adds no dot products or data passes. */
            if (fabs(uv) < 1e-30 || !(uv > 128.0 * DBL_EPSILON * ssr))
            {
                if (ssr <= 64.0 * DBL_EPSILON * DBL_EPSILON * initial_norm2)
                    return cg_breakdown_result(S, ssr, initial_norm2, iter);
                if (!(uv > 128.0 * DBL_EPSILON * ssr)) return HDFE_CG_BREAKDOWN;
            }
            if (!isfinite(uv) || !isfinite(ssr)) return HDFE_CG_BREAKDOWN;
            alpha = ssr / uv;

            /* No step from a residual at the rounding floor, nor along a
             * direction without curvature, which is rounding as well: the
             * column is finished, or its next round starts */
            if (ssr <= rounding || alpha < 0.0)
            {
                if (cg_round_end(S, y, thread_id, alpha * ssr, &round, &improvement_potential, &scale,
                                 &rounding, &ssr))
                    return iter;
                continue;
            }
            recent_ssr = alpha * ssr;
            improvement_potential -= recent_ssr;

            /* Fused: y -= alpha*u, r -= alpha*v, ssr = w·(r·r) */
            ssr_old = ssr;
            ssr = ctools_simd_fused_axpy2_wdot(y, u, r, v, weights, alpha, (size_t)N);

            /* Convergence: recent_ssr / improvement_potential <= tol^2. Checked
             * before the direction update, which a converged column discards. */
            if (improvement_potential > 1e-30 && recent_ssr <= tol_sq * improvement_potential)
            {
                return iter;
            }
            if (improvement_potential <= 1e-30)
            {
                return iter;
            }

            /* CG update: u = r + beta*u (SIMD-accelerated) */
            beta = (fabs(ssr_old) < 1e-30) ? 0.0 : ssr / ssr_old;
            ctools_simd_axpby(u, r, 1.0, beta, (size_t)N);
        }
    }
    else
    {
        improvement_potential = initial_norm2 = dot_product(y, y, N);

        transform_sym_kaczmarz_threaded(S, y, r, S->thread_fe_means, thread_id);
        ssr = dot_product(r, r, N);
        memcpy(u, r, N * sizeof(ST_double));
        scale = improvement_potential;
        rounding = cg_rounding_floor(scale);

        for (iter = 1; iter <= S->maxiter; iter++)
        {
            if (ctools_cancel_requested()) return -iter;
            if (ssr == 0.0) return iter;
            transform_sym_kaczmarz_threaded(S, u, v, S->thread_fe_means, thread_id);

            uv = dot_product(u, v, N);
            if (fabs(uv) < 1e-30 || !(uv > 128.0 * DBL_EPSILON * ssr))
            {
                if (ssr <= 64.0 * DBL_EPSILON * DBL_EPSILON * initial_norm2)
                    return cg_breakdown_result(S, ssr, initial_norm2, iter);
                if (!(uv > 128.0 * DBL_EPSILON * ssr)) return HDFE_CG_BREAKDOWN;
            }
            if (!isfinite(uv) || !isfinite(ssr)) return HDFE_CG_BREAKDOWN;
            alpha = ssr / uv;

            /* Rounding floor or no curvature: no step (see above) */
            if (ssr <= rounding || alpha < 0.0)
            {
                if (cg_round_end(S, y, thread_id, alpha * ssr, &round, &improvement_potential, &scale,
                                 &rounding, &ssr))
                    return iter;
                continue;
            }
            recent_ssr = alpha * ssr;
            improvement_potential -= recent_ssr;

            /* Fused: y -= alpha*u, r -= alpha*v, ssr = r·r */
            ssr_old = ssr;
            ssr = ctools_simd_fused_axpy2_dot(y, u, r, v, alpha, (size_t)N);

            /* Convergence: recent_ssr / improvement_potential <= tol^2. Checked
             * before the direction update, which a converged column discards. */
            if (improvement_potential > 1e-30 && recent_ssr <= tol_sq * improvement_potential)
            {
                return iter;
            }
            if (improvement_potential <= 1e-30)
            {
                return iter;
            }

            /* CG update: u = r + beta*u (SIMD-accelerated) */
            beta = (fabs(ssr_old) < 1e-30) ? 0.0 : ssr / ssr_old;
            ctools_simd_axpby(u, r, 1.0, beta, (size_t)N);
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

/* coefficients != 0: categorical effects go to the coefficient-space solver
 * (hdfe_coef.h), which keeps its structure in plan when there is one;
 * otherwise every G >= 2 runs the CG iteration below. before_rows, when not
 * NULL, is called before any path that works on the rows (see
 * partial_out_columns_plan). */
static HDFE_SolveResult partial_out_columns_impl(HDFE_State *S, ST_double *data, ST_int N, ST_int K,
                                                 int coefficients, hdfe_coef_plan *plan,
                                                 void (*before_rows)(void *), void *before_ctx)
{
    ST_int k;
    ST_int max_iters = 0;
    ST_int any_failed = 0;

    if (!S || !data || N != S->N || K < 1 || S->num_threads < 1 || S->maxiter < 1 ||
        !isfinite(S->tolerance) || S->tolerance <= 0)
    {
        return (HDFE_SolveResult){198, 0, 0};
    }
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
        if (before_rows)
            before_rows(before_ctx);
        /* rows grouped by the effect: one pass, the same residuals */
        if (demean_columns_grouped(S, data, K))
        {
            if (ctools_cancel_requested()) return (HDFE_SolveResult){1, 0, 1};
            return (HDFE_SolveResult){0, 1, 1};
        }
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

    if (coefficients && S->maxiter >= HDFE_COEF_MIN_MAXITER)
    {
        if (plan)
        {
            int declined = 0;
            const HDFE_SolveResult solved =
                hdfe_coef_partial_out_plan(S, plan, data, N, K, HDFE_COEF_AUTO, &declined);
            if (!declined)
                return solved;
        }
        else if (hdfe_coef_applies(S))
            return hdfe_coef_partial_out(S, data, N, K, HDFE_COEF_AUTO);
    }
    if (before_rows)
        before_rows(before_ctx);

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

            if (iters == HDFE_CG_BREAKDOWN)
            {
                any_failed |= 2;
            }
            else if (iters < 0)
            {
                any_failed |= 1;
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
    return (HDFE_SolveResult){(any_failed & 2) ? 498 : any_failed ? 430 : 0,
                             !any_failed, max_iters};
}

HDFE_SolveResult partial_out_columns(HDFE_State *S, ST_double *data, ST_int N, ST_int K,
                                     ST_int num_threads)
{
    (void)num_threads; /* Use S->num_threads or OpenMP default */
    return partial_out_columns_impl(S, data, N, K, 1, NULL, NULL, NULL);
}

HDFE_SolveResult partial_out_columns_plan(HDFE_State *S, ST_double *data, ST_int N, ST_int K,
                                          ST_int num_threads, hdfe_coef_plan *plan,
                                          void (*before_rows)(void *ctx), void *ctx)
{
    (void)num_threads;
    return partial_out_columns_impl(S, data, N, K, 1, plan, before_rows, ctx);
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
    /* reghdfe's default iteration, asked for by name */
    if (transform == 0 && acceleration == 1)
        return partial_out_columns_impl(S, data, N, K, 0, NULL, NULL, NULL);
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
                    {
                        /* As in the default CG path, breakdown is not a zero
                         * step. Keep the signed steps of the nonsymmetric
                         * forward Kaczmarz path outside this tiny-step case. */
                        if (fabs(denom) < 1e-30 &&
                            (rr <= 64.0 * DBL_EPSILON * DBL_EPSILON * scale ||
                             !(denom > 128.0 * DBL_EPSILON * rr)))
                        {
                            if (cg_breakdown_result(S, rr, scale, iter) == HDFE_CG_BREAKDOWN)
                                failed |= 4;
                            else
                                done = 1;
                            break;
                        }
                        step = rr / denom;
                    }
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
    return (HDFE_SolveResult){(failed & 2) ? 920 : (failed & 4) ? 498 : failed ? 430 : 0,
                             !failed, max_iterations};
}

#include "cppmlhdfe_internal.h"
#include <math.h>
#include <stdio.h>

/* Iterated rectifier separation uses equality-constrained weighted least squares
 * (Correia, Guimaraes and Zylkin). Positive-y rows receive the Van Loan penalty
 * ceil(eps^(-1/2)). Partialling out FEs keeps memory O(N K). Pivoted, twice
 * orthogonalized QR avoids squaring the condition number. Small positive fitted
 * means alone are never evidence of separation. */

/* The regressors less the fixed effects, for the basis the ReLU residual is
 * orthogonalized against. Solved for the coefficients of the effects when
 * that applies (every thread on the rows), with the absorbed-column rule of
 * ppml_reference_partial(); otherwise, and when that solve runs out of
 * iterations, the iteration on the rows. (The residual of each ReLU step
 * itself always comes from the row iteration: see ppml_reference_partial.) */
static HDFE_SolveResult ppml_relu_partial_basis(HDFE_State *S, double *q, int N, int K,
                                                int nthreads)
{
    if (S->G >= 2 && S->maxiter >= HDFE_COEF_MIN_MAXITER && hdfe_coef_applies(S))
    {
        double *before = malloc((size_t)K * sizeof(double));
        if (before)
        {
            for (int k = 0; k < K; k++)
                before[k] = ppml_quad_dot(q + (size_t)k * N, q + (size_t)k * N, NULL, N);
            HDFE_SolveResult direct = hdfe_coef_partial_out(S, q, N, K, HDFE_COEF_AUTO);
            if (!direct.status)
                for (int k = 0; k < K; k++)
                {
                    double *column = q + (size_t)k * N;
                    if (ppml_quad_dot(column, column, NULL, N) <=
                        before[k] * fmin(1e-6, S->tolerance / 10))
                        memset(column, 0, (size_t)N * sizeof(double));
                }
            free(before);
            /* 430: the columns hold the last iterate, which differs from the
             * regressors by fixed effects only; the row iteration finishes */
            if (direct.status != 430)
                return direct;
        }
    }
    return ppml_reference_partial(S, q, N, K, nthreads);
}

int ppml_relu_basis(HDFE_State *S, const double *X, int N, int K, double *q, int nthreads)
{
    if (!K)
        return 0;
    int rank = 0;

    memcpy(q, X, (size_t)N * K * sizeof(*q));
    HDFE_SolveResult projection = ppml_relu_partial_basis(S, q, N, K, nthreads);
    if (projection.status)
        return -projection.status;
    ST_int col;
    /* Each column is scaled by sqrt(w) row by row; its norm and that of the
     * unprojected column are chains of hypot() over the rows in order. The
     * 2K chains are independent, so each has its own thread (one thread per
     * column ran both of its chains one after the other); the values are
     * the same. */
    ST_double *chain = (ST_double *)malloc((size_t)2 * K * sizeof(ST_double));
    if (!chain)
        return -920;
    for (col = 0; col < K; col++)
    {
        ST_double *column = q + (size_t)col * N;
        ST_int i;
#pragma omp parallel for schedule(static) if (N > 100000)
        for (i = 0; i < N; i++)
            column[i] *= sqrt(S->weights[i]);
    }
#pragma omp parallel for schedule(dynamic, 1) if ((size_t)N * K > 100000)
    for (col = 0; col < 2 * K; col++)
    {
        const ST_int k = col / 2;
        ST_double norm = 0;
        if (col % 2 == 0)
            for (ST_int i = 0; i < N; i++)
                norm = hypot(norm, q[(size_t)k * N + i]);
        else
            for (ST_int i = 0; i < N; i++)
                norm = hypot(norm, X[(size_t)k * N + i] * sqrt(S->weights[i]));
        chain[col] = norm;
    }
    for (col = 0; col < K; col++)
    {
        const ST_int k = col;
        const ST_double norm = chain[2 * k], original_norm = chain[2 * k + 1];
        ST_double *column = q + (size_t)k * N;
        ST_int i;
        /* Do not amplify projection noise in an FE-collinear regressor
         * into an artificial separating direction. */
        if (norm < fmax(1e-9, 1e-8 * original_norm))
            memset(column, 0, (size_t)N * sizeof(*q));
        else
        {
#pragma omp parallel for schedule(static) if (N > 100000)
            for (i = 0; i < N; i++)
                column[i] /= norm;
        }
    }
    free(chain);
    ST_double *norms = (ST_double *)malloc((size_t)K * sizeof(ST_double));
    if (!norms)
        return -920;
    for (ST_int k = 0; k < K; k++)
    {
        ST_int pivot = k;
        ST_double best = 0;
        /* column norms in parallel, then the first largest as before */
#pragma omp parallel for schedule(static) if (K - k > 2 && (size_t)N * (K - k) > 200000)
        for (col = k; col < K; col++)
            norms[col] = fast_dot(q + (size_t)col * N, q + (size_t)col * N, N);
        for (ST_int j = k; j < K; j++)
        {
            if (norms[j] > best)
            {
                best = norms[j];
                pivot = j;
            }
        }
        if (best < 1e-20)
            break;
        for (ST_int i = 0; i < N; i++)
        {
            ST_double tmp = q[(size_t)k * N + i];
            q[(size_t)k * N + i] = q[(size_t)pivot * N + i];
            q[(size_t)pivot * N + i] = tmp;
            q[(size_t)k * N + i] /= sqrt(best);
        }
        rank++;
        /* the later columns are orthogonalized against column k independently */
#pragma omp parallel for schedule(static) if (K - k > 2 && (size_t)N * (K - k) > 200000)
        for (col = k + 1; col < K; col++)
        {
            for (ST_int pass = 0; pass < 2; pass++)
            {
                ST_double dot = fast_dot(q + (size_t)k * N, q + (size_t)col * N, N);
                for (ST_int i = 0; i < N; i++)
                    q[(size_t)col * N + i] -= dot * q[(size_t)k * N + i];
            }
        }
    }
    free(norms);
    return rank;
}

/* How far inside the band that counts as zero every fitted value of the
 * screen's step must be for the check to end there */
#define PPML_RELU_SCREEN_MARGIN 4.0
/* Smallest pivot, on a unit diagonal, that the screen solves the regressors'
 * cross products with */
#define PPML_RELU_SCREEN_PIVOT 1e-9

/* Screen. On data without separation the first step of the iteration ends
 * the check: the fitted values of the zero-outcome rows all fall inside the
 * band that counts as zero. This takes that step ahead of the basis: the
 * regressors and u less the fixed effects from one solve for the effects'
 * coefficients (every thread on the rows; the iteration uses one), u's
 * regression on the regressors from their cross products, and the fitted
 * values in one pass over the rows. A regressor counts as it does in the
 * basis: not when the effects absorb it or what is left of it is projection
 * noise.
 *
 * Returns 1 when every fitted value is inside the band by
 * PPML_RELU_SCREEN_MARGIN: nothing is separated and the basis is never
 * built. Otherwise nothing computed here is used. 2: a fitted value is near
 * or outside the band (as on quasi-separated panels, where what the
 * iteration on the rows concludes follows its own inexactness on rows far
 * lighter than the others: ppml_reference_partial), or the solve does not
 * apply or did not converge; the iteration decides. 0: the cross products
 * cannot tell which regressors count (one sits at a threshold of the basis,
 * or is nearly a combination of the others); the basis can, and the step is
 * taken again with it. Negative: minus the status of a failure.
 * q is scratch of N (K + 1) doubles. */
static ST_int ppml_relu_screen(const PPML_Options *options, HDFE_State *S, const ST_double *y,
                               const ST_double *X, const ST_double *u, ST_int N, ST_int K,
                               ST_int zeros, ST_double ratio, ST_double penalty, ST_double *q)
{
    if (S->G < 2 || S->maxiter < HDFE_COEF_MIN_MAXITER || !hdfe_coef_applies(S))
        return 2;
    const ST_int K1 = K + 1;
    ST_int result = 2, used = 0, i, k;
    /* size: each regressor's sum of squares, the same with the weights, and
     * the coefficients; cross: the cross products, their scaled copy for the
     * factorization, and the right-hand side */
    ST_double *size = (ST_double *)malloc((size_t)K1 * 3 * sizeof(ST_double));
    long double *cross = (long double *)malloc(((size_t)K1 * K1 * 2 + K1) * sizeof(long double));
    ST_int *order = (ST_int *)malloc((size_t)K1 * sizeof(ST_int));
    if (!size || !cross || !order)
    {
        result = -920;
        goto out;
    }
    ST_double *beta = size + (size_t)K1 * 2;
    long double *factor = cross + (size_t)K1 * K1;
    long double *rhs = factor + (size_t)K1 * K1;
    for (k = 0; k < K; k++)
    {
        const ST_double *source = X + (size_t)k * N;
        ST_double *column = q + (size_t)k * N;
#pragma omp parallel for schedule(static) if (N > 100000)
        for (i = 0; i < N; i++)
            column[i] = source[i];
        size[k] = ppml_quad_dot(source, source, NULL, N);
        size[K1 + k] = ppml_quad_dot(source, source, S->weights, N);
    }
    {
        ST_double *column = q + (size_t)K * N;
#pragma omp parallel for schedule(static) if (N > 100000)
        for (i = 0; i < N; i++)
            column[i] = u[i];
    }
    const HDFE_SolveResult direct = hdfe_coef_partial_out(S, q, N, K1, HDFE_COEF_AUTO);
    if (direct.status == 1 || direct.status == 920)
    {
        result = -direct.status;
        goto out;
    }
    if (direct.status)
        goto out;
    ctools_kahan_cross_weighted(q, S->weights, N, K1, cross);

    /* Which regressors count: the two rules of the basis (ppml_relu_partial_basis
     * and ppml_relu_basis) on these sums. Within a small factor of either
     * threshold the basis decides. */
    const ST_double absorbed = fmin(1e-6, S->tolerance / 10);
    for (k = 0; k < K; k++)
    {
        const ST_double *column = q + (size_t)k * N;
        const ST_double left = ppml_quad_dot(column, column, NULL, N);
        const ST_double norm = sqrt((ST_double)cross[(size_t)k * K1 + k]);
        const ST_double floor_absorbed = size[k] * absorbed;
        const ST_double floor_noise = fmax(1e-9, 1e-8 * sqrt(size[K1 + k]));
        if (!isfinite(left) || !isfinite(norm) || !isfinite(floor_absorbed) ||
            !isfinite(floor_noise))
            goto out;
        if (left <= floor_absorbed / 4 || norm < floor_noise / 2)
            continue;
        if (left <= floor_absorbed * 4 || norm < floor_noise * 2)
        {
            result = 0;
            goto out;
        }
        order[used++] = k;
    }

    /* u on the regressors that count: Cholesky with pivoting of their cross
     * products scaled to a unit diagonal */
    for (ST_int a = 0; a < used; a++)
    {
        const ST_int ka = order[a];
        const long double scale_a = sqrtl(cross[(size_t)ka * K1 + ka]);
        for (ST_int b = 0; b < used; b++)
        {
            const ST_int kb = order[b];
            const long double scale_b = sqrtl(cross[(size_t)kb * K1 + kb]);
            const long double product =
                ka <= kb ? cross[(size_t)ka * K1 + kb] : cross[(size_t)kb * K1 + ka];
            factor[(size_t)a * K1 + b] = product / (scale_a * scale_b);
        }
        rhs[a] = cross[(size_t)ka * K1 + K] / scale_a;
    }
    for (ST_int j = 0; j < used; j++)
    {
        ST_int pivot = j;
        for (ST_int a = j + 1; a < used; a++)
            if (factor[(size_t)a * K1 + a] > factor[(size_t)pivot * K1 + pivot])
                pivot = a;
        if (!(factor[(size_t)pivot * K1 + pivot] >= PPML_RELU_SCREEN_PIVOT))
        {
            result = 0;
            goto out;
        }
        if (pivot != j)
        {
            for (ST_int b = 0; b < used; b++)
            {
                const long double held = factor[(size_t)j * K1 + b];
                factor[(size_t)j * K1 + b] = factor[(size_t)pivot * K1 + b];
                factor[(size_t)pivot * K1 + b] = held;
            }
            for (ST_int a = 0; a < used; a++)
            {
                const long double held = factor[(size_t)a * K1 + j];
                factor[(size_t)a * K1 + j] = factor[(size_t)a * K1 + pivot];
                factor[(size_t)a * K1 + pivot] = held;
            }
            const long double held = rhs[j];
            rhs[j] = rhs[pivot];
            rhs[pivot] = held;
            const ST_int moved = order[j];
            order[j] = order[pivot];
            order[pivot] = moved;
        }
        const long double diagonal = sqrtl(factor[(size_t)j * K1 + j]);
        factor[(size_t)j * K1 + j] = diagonal;
        for (ST_int a = j + 1; a < used; a++)
            factor[(size_t)a * K1 + j] /= diagonal;
        for (ST_int a = j + 1; a < used; a++)
            for (ST_int b = j + 1; b <= a; b++)
            {
                factor[(size_t)a * K1 + b] -=
                    factor[(size_t)a * K1 + j] * factor[(size_t)b * K1 + j];
                factor[(size_t)b * K1 + a] = factor[(size_t)a * K1 + b];
            }
    }
    for (ST_int a = 0; a < used; a++)
    {
        for (ST_int b = 0; b < a; b++)
            rhs[a] -= factor[(size_t)a * K1 + b] * rhs[b];
        rhs[a] /= factor[(size_t)a * K1 + a];
    }
    for (ST_int a = used - 1; a >= 0; a--)
    {
        for (ST_int b = a + 1; b < used; b++)
            rhs[a] -= factor[(size_t)b * K1 + a] * rhs[b];
        rhs[a] /= factor[(size_t)a * K1 + a];
    }
    for (ST_int a = 0; a < used; a++)
        beta[a] = (ST_double)(rhs[a] / sqrtl(cross[(size_t)order[a] * K1 + order[a]]));

    const ST_double delta =
        fmax((ST_double)zeros / penalty * ratio * ratio, 1e-8) + options->relu_tol;
    const ST_double *projected = q + (size_t)K * N;
    ST_double lowest = 0.0, highest = 0.0;
    int not_finite = 0;
#pragma omp parallel for schedule(static) reduction(min : lowest) reduction(max : highest)       \
    reduction(| : not_finite) if (N > 100000)
    for (i = 0; i < N; i++)
        if (y[i] == 0)
        {
            ST_double residual = projected[i];
            for (ST_int a = 0; a < used; a++)
                residual -= beta[a] * q[(size_t)order[a] * N + i];
            const ST_double fitted = u[i] - residual;
            if (!isfinite(fitted))
                not_finite |= 1;
            else
            {
                if (fitted < lowest)
                    lowest = fitted;
                if (fitted > highest)
                    highest = fitted;
            }
        }
    if (!not_finite && lowest >= -0.1 * delta / PPML_RELU_SCREEN_MARGIN &&
        highest <= delta / PPML_RELU_SCREEN_MARGIN)
        result = 1;
out:
    free(size);
    free(cross);
    free(order);
    return result;
}

ST_int ppml_relu_ex(const PPML_Options *options, HDFE_State *S, const ST_double *y,
                    const ST_double *X, ST_int N, ST_int K, ST_int *separated, ST_int nthreads,
                    ST_double *certificate)
{
    return ppml_relu_ex2(options, S, y, X, N, K, separated, nthreads, certificate, NULL);
}

ST_int ppml_relu_ex2(const PPML_Options *options, HDFE_State *S, const ST_double *y,
                     const ST_double *X, ST_int N, ST_int K, ST_int *separated, ST_int nthreads,
                     ST_double *certificate, ST_int *iterations)
{
    ST_int zeros = 0, rank = 0, rc = 0;
    if (iterations)
        *iterations = 0;
    ST_double saved_tol = S->tolerance;
    const ST_double penalty = 67108864.0;
    memset(separated, 0, (size_t)N * sizeof(*separated));
#pragma omp parallel for schedule(static) reduction(+ : zeros) if (N > 100000)
    for (ST_int i = 0; i < N; i++)
        zeros += y[i] == 0;
    if (!zeros)
        return 0;
    ST_double *q = malloc((size_t)N * (K + 1) * sizeof(*q));
    ST_double *u = malloc((size_t)N * sizeof(*u));
    ST_double *r = malloc((size_t)N * sizeof(*r));
    ST_double *fit = malloc((size_t)N * sizeof(*fit));
    ST_double *last = calloc((size_t)N * 4, sizeof(*last));
    if (!q || !u || !r || !fit || !last)
    {
        rc = -920;
        goto done;
    }
    /* The loops over the rows below that are marked parallel compute each
     * row from that row alone (minima and counts are exact in any order),
     * so the threads change no result. */
#pragma omp parallel for schedule(static) if (N > 100000)
    for (ST_int i = 0; i < N; i++)
    {
        S->weights[i] = y[i] == 0 ? 1.0 : penalty;
        u[i] = y[i] == 0;
    }
    ppml_update_fe_weights(S, S->weights, N);
    S->tolerance = fmax(1e-13, options->relu_tol * options->relu_tol);
    ST_double norm_zero = 0.0, norm_positive = 0.0;
    /* (a column's two sums run over its rows in order; the columns are
     * independent and their maxima exact in any order) */
#pragma omp parallel for schedule(dynamic, 1) reduction(max : norm_zero, norm_positive)           \
    if (K > 1 && (size_t)N * K > 200000)
    for (ST_int k = 0; k < K; k++)
    {
        ST_double zero_norm = 0.0, positive_norm = 0.0;
        for (ST_int i = 0; i < N; i++)
        {
            if (y[i] == 0)
                zero_norm += fabs(X[(size_t)k * N + i]);
            else
                positive_norm += fabs(X[(size_t)k * N + i]);
        }
        norm_zero = fmax(norm_zero, zero_norm);
        norm_positive = fmax(norm_positive, positive_norm);
    }
    ST_double ratio = norm_positive > 0 ? fmax(norm_zero / norm_positive, 1.0) : 1.0;
    const ST_int screened = ppml_relu_screen(options, S, y, X, u, N, K, zeros, ratio, penalty, q);
    if (screened < 0 || screened == 1)
    {
        if (screened < 0)
            rc = screened;
        else if (iterations)
            *iterations = 1;
        goto done;
    }
    rank = ppml_relu_basis(S, X, N, K, q, nthreads);
    if (rank < 0)
    {
        rc = rank;
        goto done;
    }
    /* Where the screen could not tell which regressors count, the same step
     * with the basis, which can: u less the fixed effects from the solve for
     * their coefficients, orthogonalized against the basis. Inside the band
     * by PPML_RELU_SCREEN_MARGIN there is nothing to do; otherwise nothing
     * of it is used and the iteration on the rows decides as before. */
    if (!screened && S->G >= 2 && S->maxiter >= HDFE_COEF_MIN_MAXITER && hdfe_coef_applies(S))
    {
        ST_int i;
#pragma omp parallel for schedule(static) if (N > 100000)
        for (i = 0; i < N; i++)
            r[i] = u[i];
        const HDFE_SolveResult direct = hdfe_coef_partial_out(S, r, N, 1, HDFE_COEF_AUTO);
        if (direct.status == 1 || direct.status == 920)
        {
            rc = -direct.status;
            goto done;
        }
        if (!direct.status)
        {
#pragma omp parallel for schedule(static) if (N > 100000)
            for (i = 0; i < N; i++)
                r[i] *= sqrt(S->weights[i]);
            for (ST_int pass = 0; pass < 2; pass++)
                for (ST_int k = 0; k < rank; k++)
                {
                    const ST_double dot = fast_dot(q + (size_t)k * N, r, N);
#pragma omp parallel for schedule(static) if (N > 100000)
                    for (i = 0; i < N; i++)
                        r[i] -= dot * q[(size_t)k * N + i];
                }
            const ST_double delta =
                fmax(fast_dot(u, u, N) / penalty * ratio * ratio, 1e-8) + options->relu_tol;
            ST_double lowest = 0.0, highest = 0.0;
            int not_finite = 0;
#pragma omp parallel for schedule(static) reduction(min : lowest) reduction(max : highest)       \
    reduction(| : not_finite) if (N > 100000)
            for (i = 0; i < N; i++)
                if (y[i] == 0)
                {
                    const ST_double fitted = u[i] - r[i] / sqrt(S->weights[i]);
                    if (!isfinite(fitted))
                        not_finite |= 1;
                    else
                    {
                        if (fitted < lowest)
                            lowest = fitted;
                        if (fitted > highest)
                            highest = fitted;
                    }
                }
            if (!not_finite && lowest >= -0.1 * delta / PPML_RELU_SCREEN_MARGIN &&
                highest <= delta / PPML_RELU_SCREEN_MARGIN)
            {
                if (iterations)
                    *iterations = 1;
                goto done;
            }
        }
    }
    double ee_cumulative = 0, ee_boundary = zeros, progress1 = 0, progress2 = 0,
           acceleration_value = 1;
    int candidates1 = -1, candidates2 = -1, stuck = 0;
    for (ST_int iter = 0; iter < options->relu_maxiter; iter++)
    {
        if (ctools_cancel_requested()) { rc = -1; goto done; }
        if (iterations)
            *iterations = iter + 1;
#pragma omp parallel for schedule(static) if (N > 100000)
        for (int i = 0; i < N; i++)
        {
            r[i] = u[i] + last[i] - last[(size_t)N + i];
            last[(size_t)N + i] = u[i];
        }
        HDFE_SolveResult projection = ppml_reference_partial(S, r, N, 1, nthreads);
        memcpy(last, r, (size_t)N * sizeof(*r));
        if (projection.status)
        {
            rc = -projection.status;
            goto done;
        }
#pragma omp parallel for schedule(static) if (N > 100000)
        for (ST_int i = 0; i < N; i++)
            r[i] *= sqrt(S->weights[i]);
        for (ST_int pass = 0; pass < 2; pass++)
        {
            for (ST_int k = 0; k < rank; k++)
            {
                ST_double dot = fast_dot(q + (size_t)k * N, r, N);
#pragma omp parallel for schedule(static) if (N > 100000)
                for (ST_int i = 0; i < N; i++)
                    r[i] -= dot * q[(size_t)k * N + i];
            }
        }
        ST_double delta =
            fmax(fast_dot(u, u, N) / penalty * ratio * ratio, 1e-8) + options->relu_tol;
        ST_double min_fit = 0.0, min_resid = 0.0;
#pragma omp parallel for schedule(static) reduction(min : min_fit, min_resid) if (N > 100000)
        for (ST_int i = 0; i < N; i++)
        {
            r[i] /= sqrt(S->weights[i]);
            fit[i] = y[i] == 0 ? u[i] - r[i] : 0.0;
            if (fit[i] >= -0.1 * delta && fit[i] <= delta)
                fit[i] = 0;
            if (y[i] == 0)
            {
                min_fit = fmin(min_fit, fit[i]);
                min_resid = fmin(min_resid, fabs(r[i]) < options->relu_zero_tol ? 0 : r[i]);
            }
        }
        int candidates = 0;
#pragma omp parallel for schedule(static) reduction(+ : candidates) if (N > 100000)
        for (int i = 0; i < N; i++)
            candidates += y[i] == 0 && fit[i] > delta;
        double ee = ppml_quad_dot(r, r, NULL, N);
        ee_cumulative += iter ? 2 * ee : ee;
        double progress = 100 * ee_cumulative / ee_boundary;
        if (options->relu_accelerate && iter >= 3 && progress - progress2 < 1 &&
            candidates == candidates2)
            stuck = 1;
        progress2 = progress1;
        progress1 = progress;
        candidates2 = candidates1;
        candidates1 = candidates;
        if (min_fit >= 0 || min_resid >= 0)
        {
            for (ST_int i = 0; i < N; i++)
            {
                separated[i] = y[i] == 0 && fit[i] > 0 && (min_fit >= 0 || r[i] <= delta);
                if (certificate)
                    certificate[i] = separated[i] ? fit[i] : 0;
                rc += separated[i];
            }
            goto done;
        }
        int changed_weights = 0, accelerated = 0;
        if (stuck)
            for (int i = 0; i < N; i++)
                accelerated += y[i] == 0 &&
                               last[(size_t)3 * N + i] < 1.01 * last[(size_t)2 * N + i] &&
                               last[(size_t)2 * N + i] < 1.01 * fit[i] && fit[i] < -.1 * delta;
        acceleration_value = accelerated ? fmin(256, 4 * acceleration_value) : 1;
        for (int i = 0; i < N; i++)
        {
            double w = y[i] > 0 ? penalty : 1;
            if (stuck && y[i] == 0 && last[(size_t)3 * N + i] < 1.01 * last[(size_t)2 * N + i] &&
                last[(size_t)2 * N + i] < 1.01 * fit[i] && fit[i] < -.1 * delta)
                w = acceleration_value;
            changed_weights |= S->weights[i] != w;
            S->weights[i] = w;
            last[(size_t)3 * N + i] = last[(size_t)2 * N + i];
            last[(size_t)2 * N + i] = fit[i];
            u[i] = fmax(fit[i], 0);
        }
        if (changed_weights)
        {
            ppml_update_fe_weights(S, S->weights, N);
            rank = ppml_relu_basis(S, X, N, K, q, nthreads);
            if (rank < 0)
            {
                rc = rank;
                goto done;
            }
        }
    }
    /* As in ppmlhdfe's default non-strict ReLU check, a limit does not
     * certify any observation as separated. IRLS still must converge. */
    if (options->relu_strict)
        rc = -9010;
done:
    S->tolerance = saved_tol;
    free(q);
    free(u);
    free(r);
    free(fit);
    free(last);
    return rc;
}

ST_int ppml_relu(const PPML_Options *options, HDFE_State *S, const double *y, const double *X,
                 int N, int K, int *mask, int threads)
{
    return ppml_relu_ex(options, S, y, X, N, K, mask, threads, NULL);
}

void ppml_format_small_g(char *buffer, size_t size, double value)
{
    /* Mata's %g of the small positive values ppmlhdfe reports: four
     * decimals and a two-digit exponent (three decimals with three digits) */
    if (!(value > 0) || !isfinite(value)) {
        snprintf(buffer, size, "%g", value);
        return;
    }
    int exponent = (int)floor(log10(value));
    double mantissa = value / pow(10.0, exponent);
    int digits = (exponent <= -100 || exponent >= 100) ? 3 : 4;
    double scale = pow(10.0, digits);
    mantissa = floor(mantissa * scale + 0.5) / scale;
    if (mantissa >= 10.0) {
        mantissa /= 10.0;
        exponent++;
    }
    else if (mantissa < 1.0) {
        mantissa *= 10.0;
        exponent--;
    }
    snprintf(buffer, size, "%.*fe%c%0*d", digits, mantissa, exponent < 0 ? '-' : '+',
             (exponent <= -100 || exponent >= 100) ? 3 : 2, exponent < 0 ? -exponent : exponent);
}

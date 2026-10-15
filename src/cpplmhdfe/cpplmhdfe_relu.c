#include "cpplmhdfe_internal.h"

/* Iterated rectifier separation uses equality-constrained weighted least squares
 * (Correia, Guimaraes and Zylkin). Positive-y rows receive the Van Loan penalty
 * ceil(eps^(-1/2)). Partialling out FEs keeps memory O(N K). Pivoted, twice
 * orthogonalized QR avoids squaring the condition number. Small positive fitted
 * means alone are never evidence of separation. */

int ppml_relu_basis(HDFE_State *S, const double *X, int N, int K, double *q, int nthreads)
{
    if (!K)
        return 0;
    int rank = 0;

    memcpy(q, X, (size_t)N * K * sizeof(*q));
    HDFE_SolveResult projection = ppml_reference_partial(S, q, N, K, nthreads);
    if (projection.status)
        return -projection.status;
    ST_int col;
#pragma omp parallel for schedule(static) if (K > 1 && (size_t)N * K > 200000)
    for (col = 0; col < K; col++)
    {
        const ST_int k = col;
        ST_double norm = 0, original_norm = 0;
        for (ST_int i = 0; i < N; i++)
        {
            q[(size_t)k * N + i] *= sqrt(S->weights[i]);
            norm = hypot(norm, q[(size_t)k * N + i]);
            original_norm = hypot(original_norm, X[(size_t)k * N + i] * sqrt(S->weights[i]));
        }
        /* Do not amplify projection noise in an FE-collinear regressor
         * into an artificial separating direction. */
        if (norm < fmax(1e-9, 1e-8 * original_norm))
            memset(q + (size_t)k * N, 0, (size_t)N * sizeof(*q));
        else
            for (ST_int i = 0; i < N; i++)
                q[(size_t)k * N + i] /= norm;
    }
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

ST_int ppml_relu_ex(const PPML_Options *options, HDFE_State *S, const ST_double *y,
                    const ST_double *X, ST_int N, ST_int K, ST_int *separated, ST_int nthreads,
                    ST_double *certificate)
{
    ST_int zeros = 0, rank = 0, rc = 0;
    ST_double saved_tol = S->tolerance;
    const ST_double penalty = 67108864.0;
    memset(separated, 0, (size_t)N * sizeof(*separated));
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
    for (ST_int i = 0; i < N; i++)
    {
        S->weights[i] = y[i] == 0 ? 1.0 : penalty;
        u[i] = y[i] == 0;
    }
    ppml_update_fe_weights(S, S->weights, N);
    S->tolerance = fmax(1e-13, options->relu_tol * options->relu_tol);
    ST_double norm_zero = 0.0, norm_positive = 0.0;
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
    rank = ppml_relu_basis(S, X, N, K, q, nthreads);
    if (rank < 0)
    {
        rc = rank;
        goto done;
    }
    double ee_cumulative = 0, ee_boundary = zeros, progress1 = 0, progress2 = 0,
           acceleration_value = 1;
    int candidates1 = -1, candidates2 = -1, stuck = 0;
    for (ST_int iter = 0; iter < options->relu_maxiter; iter++)
    {
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
        for (ST_int i = 0; i < N; i++)
            r[i] *= sqrt(S->weights[i]);
        for (ST_int pass = 0; pass < 2; pass++)
        {
            for (ST_int k = 0; k < rank; k++)
            {
                ST_double dot = fast_dot(q + (size_t)k * N, r, N);
                for (ST_int i = 0; i < N; i++)
                    r[i] -= dot * q[(size_t)k * N + i];
            }
        }
        ST_double delta =
            fmax(fast_dot(u, u, N) / penalty * ratio * ratio, 1e-8) + options->relu_tol;
        ST_double min_fit = 0.0, min_resid = 0.0;
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

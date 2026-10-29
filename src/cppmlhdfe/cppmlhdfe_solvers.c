#include "cppmlhdfe_internal.h"

double ppml_quad_dot(const double *x, const double *y, const double *w, int n)
{
    dd_real sum = {0, 0};
    for (int i = 0; i < n; i++)
    {
        double xy = x[i] * y[i];
        dd_real prod = two_prod(xy, w ? w[i] : 1.0);
        sum = dd_add_d(dd_add_d(sum, prod.hi), prod.lo);
    }
    return sum.hi + sum.lo;
}

int ppml_qrsolve(const double *X, const double *y, const double *w, int N, int K, double *b,
                 double *resid)
{
    double *q = malloc((size_t)N * (K ? K : 1) * sizeof(double));
    double *R = calloc((size_t)(K ? K * K : 1), sizeof(double));
    double *qy = calloc(K ? K : 1, sizeof(double));
    int *order = malloc((K ? K : 1) * sizeof(int));
    if (!q || !R || !qy || !order)
    {
        free(q);
        free(R);
        free(qy);
        free(order);
        return 920;
    }
    double maximum = 0;
    int rank = 0;
    for (int k = 0; k < K; k++)
    {
        order[k] = k;
        if (b)
            b[k] = 0;
        for (int i = 0; i < N; i++)
            q[(size_t)k * N + i] = X[(size_t)k * N + i] * sqrt(w ? w[i] : 1);
        maximum = fmax(maximum, ppml_quad_dot(q + (size_t)k * N, q + (size_t)k * N, NULL, N));
    }
    for (int k = 0; k < K; k++)
    {
        int pivot = k;
        double norm = 0;
        for (int j = k; j < K; j++)
        {
            double z = ppml_quad_dot(q + (size_t)j * N, q + (size_t)j * N, NULL, N);
            if (z > norm)
            {
                norm = z;
                pivot = j;
            }
        }
        if (norm <= fmax(DBL_MIN, maximum * 1e-24))
            break;
        int index = order[k];
        order[k] = order[pivot];
        order[pivot] = index;
        for (int i = 0; i < N; i++)
        {
            double z = q[(size_t)k * N + i];
            q[(size_t)k * N + i] = q[(size_t)pivot * N + i];
            q[(size_t)pivot * N + i] = z;
        }
        for (int i = 0; i < k; i++)
        {
            double z = R[i * K + k];
            R[i * K + k] = R[i * K + pivot];
            R[i * K + pivot] = z;
        }
        R[k * K + k] = sqrt(norm);
        for (int i = 0; i < N; i++)
            q[(size_t)k * N + i] /= R[k * K + k];
        dd_real sum = {0, 0};
        for (int i = 0; i < N; i++)
            sum = dd_add_d(sum, q[(size_t)k * N + i] * sqrt(w ? w[i] : 1) * y[i]);
        qy[k] = sum.hi + sum.lo;
        for (int j = k + 1; j < K; j++)
            for (int pass = 0; pass < 2; pass++)
            {
                double dot = ppml_quad_dot(q + (size_t)k * N, q + (size_t)j * N, NULL, N);
                R[k * K + j] += dot;
                for (int i = 0; i < N; i++)
                    q[(size_t)j * N + i] -= dot * q[(size_t)k * N + i];
            }
        rank++;
    }
    if (resid)
        for (int i = 0; i < N; i++)
        {
            double z = 0;
            for (int k = 0; k < rank; k++)
                z += q[(size_t)k * N + i] * qy[k];
            resid[i] = y[i] - z / sqrt(w ? w[i] : 1);
        }
    if (b)
        for (int k = rank - 1; k >= 0; k--)
        {
            double z = qy[k];
            for (int j = k + 1; j < rank; j++)
                z -= R[k * K + j] * b[order[j]];
            b[order[k]] = z / R[k * K + k];
        }
    free(q);
    free(R);
    free(qy);
    free(order);
    return 0;
}

int ppml_joint_slope(const HDFE_State *S, int g)
{
    if (S->factors[g].slope || g + 1 >= S->G)
        return 0;
    const FE_Factor *f = &S->factors[g], *h = &S->factors[g + 1];
    return h->slope && h->slope_center && f->num_levels == h->num_levels &&
           !memcmp(f->levels, h->levels, (size_t)S->N * sizeof(int)) &&
           !(g + 2 < S->G && S->factors[g + 2].slope &&
             S->factors[g + 2].num_levels == f->num_levels &&
             !memcmp(f->levels, S->factors[g + 2].levels, (size_t)S->N * sizeof(int)));
}

/* Projection on factor g using the level buffers of thread slot tid; the
 * coefficients are left in S->thread_fe_means[tid * G + g] (and g + 1 for a
 * joint slope). */
static void ppml_project_slot(HDFE_State *S, double *y, int g, int tid)
{
    FE_Factor *f = &S->factors[g];
    double *a = S->thread_fe_means[tid * S->G + g];
    int joint = ppml_joint_slope(S, g);
    FE_Factor *h = joint ? &S->factors[g + 1] : NULL;
    double *b = joint ? S->thread_fe_means[tid * S->G + g + 1] : NULL;
    memset(a, 0, f->num_levels * sizeof(double));
    for (int i = 0; i < S->N; i++)
    {
        int l = f->levels[i] - 1;
        double x = f->slope ? f->slope[i] - (f->slope_center ? f->slope_center[l] : 0) : 1;
        a[l] += y[i] * x * S->weights[i];
    }
    for (int l = 0; l < f->num_levels; l++)
        a[l] *= f->inv_weighted_counts[l];
    if (joint)
    {
        memset(b, 0, f->num_levels * sizeof(double));
        for (int i = 0; i < S->N; i++)
        {
            int l = f->levels[i] - 1;
            b[l] += (h->slope[i] - h->slope_center[l]) * S->weights[i] * (y[i] - a[l]);
        }
        for (int l = 0; l < f->num_levels; l++)
            b[l] *= h->inv_weighted_counts[l];
    }
    for (int i = 0; i < S->N; i++)
    {
        int l = f->levels[i] - 1;
        double x = f->slope ? f->slope[i] - (f->slope_center ? f->slope_center[l] : 0) : 1;
        double fitted = x * a[l];
        if (joint)
            fitted += (h->slope[i] - h->slope_center[l]) * b[l];
        y[i] -= fitted;
    }
}

void ppml_project(HDFE_State *S, double *y, int g)
{
    ppml_project_slot(S, y, g, 0);
}

static void ppml_transform_slot(HDFE_State *S, const double *y, double *ans, int tid)
{
    memcpy(ans, y, (size_t)S->N * sizeof(double));
    int groups[10], n = 0;
    for (int g = 0; g < S->G; g++)
    {
        groups[n++] = g;
        if (ppml_joint_slope(S, g))
            g++;
    }
    for (int g = 0; g < n; g++)
        ppml_project_slot(S, ans, groups[g], tid);
    for (int g = n - 2; g >= 0; g--)
        ppml_project_slot(S, ans, groups[g], tid);
    for (int i = 0; i < S->N; i++)
        ans[i] = y[i] - ans[i];
}

void ppml_transform(HDFE_State *S, const double *y, double *ans)
{
    ppml_transform_slot(S, y, ans, 0);
}

/* Reference projection: CG on the symmetric Kaczmarz operator with
 * double-double inner products, one column at a time. Columns are independent
 * and each uses the scratch buffers of its own thread slot, so they are solved
 * in parallel with results that do not depend on the thread count. */
HDFE_SolveResult ppml_reference_partial(HDFE_State *S, double *data, int N, int K, int threads)
{
    (void)threads;
    if (S->G == 1)
        return partial_out_columns(S, data, N, K, 1);
    int slots = S->num_threads < K ? S->num_threads : K;
    if (slots < 1)
        slots = 1;
    int maxiter = 0, failed = 0;
    int k;
#pragma omp parallel for num_threads(slots) schedule(dynamic) reduction(max : maxiter) reduction(| : failed) if (slots > 1)
    for (k = 0; k < K; k++)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        double *r = S->thread_cg_r[tid], *u = S->thread_cg_u[tid], *v = S->thread_cg_v[tid];
        double *y = data + (size_t)k * N;
        double initial = ppml_quad_dot(y, y, NULL, N);
        double potential = ppml_quad_dot(y, y, S->weights, N);
        ppml_transform_slot(S, y, r, tid);
        double ssr = ppml_quad_dot(r, r, S->weights, N);
        memcpy(u, r, N * sizeof(double));
        int iter;
        for (iter = 1; iter <= S->maxiter; iter++)
        {
            ppml_transform_slot(S, u, v, tid);
            double alpha = ssr / fmax(ppml_quad_dot(u, v, S->weights, N), DBL_EPSILON);
            double recent = alpha * ssr;
            potential -= recent;
            for (int i = 0; i < N; i++)
            {
                y[i] -= alpha * u[i];
                r[i] -= alpha * v[i];
            }
            double old = ssr;
            ssr = ppml_quad_dot(r, r, S->weights, N);
            double beta = ssr / fmax(old, DBL_EPSILON);
            for (int i = 0; i < N; i++)
                u[i] = r[i] + beta * u[i];
            double error = fabs(recent) < 1e-15 ? 0 : recent / fmax(potential, DBL_EPSILON);
            if (error <= S->tolerance * S->tolerance)
                break;
        }
        if (iter > S->maxiter)
        {
            failed |= 1;
            continue;
        }
        if (iter > maxiter)
            maxiter = iter;
        if (ppml_quad_dot(y, y, NULL, N) <= initial * fmin(1e-6, S->tolerance / 10))
            memset(y, 0, N * sizeof(double));
    }
    if (failed)
        return (HDFE_SolveResult){430, 0, S->maxiter + 1};
    return (HDFE_SolveResult){0, 1, maxiter};
}

double ppml_norm(const double *x, int n)
{
    return sqrt(fmax(0, ppml_quad_dot(x, x, NULL, n)));
}

void ppml_fe_multiply(HDFE_State *S, const int *offset, const double *x, double *out, int transpose)
{
    int M = offset[S->G];
    memset(out, 0, (size_t)(transpose ? M : S->N) * sizeof(double));
    for (int g = 0; g < S->G; g++)
    {
        FE_Factor *f = &S->factors[g];
        for (int i = 0; i < S->N; i++)
        {
            int l = f->levels[i] - 1;
            double loading =
                f->slope ? f->slope[i] - (f->slope_center ? f->slope_center[l] : 0) : 1;
            double a = sqrt(S->weights[i]) * loading * sqrt(f->inv_weighted_counts[l]);
            if (transpose)
                out[offset[g] + l] += a * x[i];
            else
                out[i] += a * x[offset[g] + l];
        }
    }
}

HDFE_SolveResult ppml_lsmr(const PPML_Options *options, HDFE_State *S, double *data, int N, int K)
{
    int offset[11] = {0}, M = 0, max_iterations = 0, rc = 0;
    for (int g = 0; g < S->G; g++)
    {
        if (S->factors[g].num_levels > INT_MAX - M)
            return (HDFE_SolveResult){920, 0, 0};
        M += S->factors[g].num_levels;
        offset[g + 1] = M;
    }
    double *obs = calloc((size_t)N * 3, sizeof(double));
    double *coef = calloc((size_t)M * 6, sizeof(double));
    if (!obs || !coef)
    {
        free(obs);
        free(coef);
        return (HDFE_SolveResult){920, 0, 0};
    }
    double *u = obs, *av = obs + N, *b = obs + (size_t)2 * N;
    double *v = coef, *atu = coef + M, *h = coef + (size_t)2 * M, *hbar = coef + (size_t)3 * M,
           *x = coef + (size_t)4 * M;
    for (int col = 0; col < K; col++)
    {
        double *y = data + (size_t)col * N;
        for (int i = 0; i < N; i++)
            u[i] = b[i] = y[i] * sqrt(S->weights[i]);
        double beta = ppml_norm(u, N), normb = beta;
        if (beta == 0)
            continue;
        for (int i = 0; i < N; i++)
            u[i] /= beta;
        ppml_fe_multiply(S, offset, u, v, 1);
        double alpha = ppml_norm(v, M);
        if (alpha == 0)
            continue;
        for (int j = 0; j < M; j++)
        {
            v[j] /= alpha;
            h[j] = v[j];
            hbar[j] = x[j] = 0;
        }
        double zetabar = alpha * beta, alphabar = alpha, rho = 1, rhobar = 1, cbar = 1, sbar = 0;
        double betadd = beta, betad = 0, rhodold = 1, tautildeold = 0, thetatilde = 0, zeta = 0;
        double normA2 = alpha * alpha, maxrbar = 0, minrbar = DBL_MAX;
        int iter, converged = 0;
        for (iter = 1; iter <= S->maxiter; iter++)
        {
            ppml_fe_multiply(S, offset, v, av, 0);
            for (int i = 0; i < N; i++)
                u[i] = av[i] - alpha * u[i];
            beta = ppml_norm(u, N);
            if (beta > DBL_EPSILON)
                for (int i = 0; i < N; i++)
                    u[i] /= beta;
            ppml_fe_multiply(S, offset, u, atu, 1);
            for (int j = 0; j < M; j++)
                v[j] = atu[j] - beta * v[j];
            alpha = ppml_norm(v, M);
            if (alpha > DBL_EPSILON)
                for (int j = 0; j < M; j++)
                    v[j] /= alpha;
            double rhoold = rho;
            rho = hypot(alphabar, beta);
            if (rho == 0)
            {
                converged = 1;
                break;
            }
            double c = alphabar / rho, s = beta / rho, thetanew = s * alpha;
            alphabar = c * alpha;
            double rhobarold = rhobar, zetaold = zeta, thetabar = sbar * rho, rhotemp = cbar * rho;
            rhobar = hypot(rhotemp, thetanew);
            if (rhobar == 0)
            {
                converged = 1;
                break;
            }
            cbar = rhotemp / rhobar;
            sbar = thetanew / rhobar;
            zeta = cbar * zetabar;
            zetabar = -sbar * zetabar;
            for (int j = 0; j < M; j++)
            {
                hbar[j] = h[j] - (thetabar * rho / (rhoold * rhobarold)) * hbar[j];
                x[j] += zeta / (rho * rhobar) * hbar[j];
                h[j] = v[j] - (thetanew / rho) * h[j];
            }
            double betahat = c * betadd;
            betadd = -s * betadd;
            double oldtheta = thetatilde, rhotildeold = hypot(rhodold, thetabar);
            double ct = rhodold / rhotildeold, st = thetabar / rhotildeold;
            thetatilde = st * rhobar;
            rhodold = ct * rhobar;
            betad = -st * betad + ct * betahat;
            tautildeold = (zetaold - oldtheta * tautildeold) / rhotildeold;
            double taud = (zeta - thetatilde * tautildeold) / rhodold;
            double normr = hypot(betad - taud, betadd);
            normA2 += beta * beta;
            double normA = sqrt(normA2);
            normA2 += alpha * alpha;
            maxrbar = fmax(maxrbar, rhobarold);
            if (iter > 1)
                minrbar = fmin(minrbar, rhobarold);
            double condA =
                fmax(maxrbar, fabs(rhotemp)) / fmax(DBL_MIN, fmin(minrbar, fabs(rhotemp)));
            double rtol = options->btol + S->tolerance * normA * ppml_norm(x, M) / normb;
            if (normr <= rtol * normb || fabs(zetabar) <= S->tolerance * normA * normr ||
                condA >= options->conlim)
            {
                converged = 1;
                break;
            }
        }
        if (!converged)
        {
            rc = 430;
            break;
        }
        if (iter > max_iterations)
            max_iterations = iter;
        ppml_fe_multiply(S, offset, x, av, 0);
        for (int i = 0; i < N; i++)
        {
            if (S->weights[i] > 0)
                y[i] = (b[i] - av[i]) / sqrt(S->weights[i]);
            if (!isfinite(y[i]))
            {
                rc = 430;
                break;
            }
        }
        if (rc)
            break;
    }
    free(obs);
    free(coef);
    return (HDFE_SolveResult){rc, !rc, max_iterations};
}

void ppml_transform_option(const PPML_Options *options, HDFE_State *S, const double *y, double *out,
                           double *work)
{
    int groups[10], n = 0;
    for (int g = 0; g < S->G; g++)
    {
        groups[n++] = g;
        if (ppml_joint_slope(S, g))
            g++;
    }
    if (options->transform == 2)
    {
        memset(out, 0, (size_t)S->N * sizeof(double));
        for (int g = 0; g < n; g++)
        {
            memcpy(work, y, (size_t)S->N * sizeof(double));
            ppml_project(S, work, groups[g]);
            for (int i = 0; i < S->N; i++)
                out[i] += (y[i] - work[i]) / n;
        }
    }
    else
    {
        memcpy(out, y, (size_t)S->N * sizeof(double));
        for (int g = 0; g < n; g++)
            ppml_project(S, out, groups[g]);
        if (!options->transform)
            for (int g = n - 2; g >= 0; g--)
                ppml_project(S, out, groups[g]);
        for (int i = 0; i < S->N; i++)
            out[i] = y[i] - out[i];
    }
}

HDFE_SolveResult ppml_map_options(const PPML_Options *options, HDFE_State *S, double *data, int N,
                                  int K)
{
    double *buffer = calloc((size_t)N * 5, sizeof(double));
    if (!buffer)
        return (HDFE_SolveResult){920, 0, 0};
    double *r = buffer, *u = buffer + N, *v = buffer + (size_t)2 * N,
           *work = buffer + (size_t)3 * N, *last = buffer + (size_t)4 * N;
    int max_iterations = 0, rc = 0;
    for (int k = 0; k < K; k++)
    {
        double *y = data + (size_t)k * N;
        double initial = ppml_quad_dot(y, y, NULL, N),
               potential = ppml_quad_dot(y, y, S->weights, N), ssr = 0;
        int method = options->acceleration, iter;
        if (method == 0 || method == 4)
        {
            ppml_transform_option(options, S, y, r, work);
            ssr = ppml_quad_dot(r, r, S->weights, N);
            memcpy(u, r, (size_t)N * sizeof(double));
        }
        for (iter = 1; iter <= S->maxiter; iter++)
        {
            double error = 0;
            if (method == 0 || (method == 4 && iter > 6))
            {
                if (method == 4 && iter == 7)
                {
                    ppml_transform_option(options, S, y, r, work);
                    ssr = ppml_quad_dot(r, r, S->weights, N);
                    potential = ppml_quad_dot(y, y, S->weights, N);
                    memcpy(u, r, (size_t)N * sizeof(double));
                }
                ppml_transform_option(options, S, u, v, work);
                double alpha = ssr / fmax(ppml_quad_dot(u, v, S->weights, N), DBL_EPSILON);
                double recent = alpha * ssr;
                potential -= recent;
                for (int i = 0; i < N; i++)
                {
                    y[i] -= alpha * u[i];
                    r[i] -= alpha * v[i];
                }
                double old = ssr;
                ssr = ppml_quad_dot(r, r, S->weights, N);
                double beta = ssr / fmax(old, DBL_EPSILON);
                for (int i = 0; i < N; i++)
                    u[i] = r[i] + beta * u[i];
                error = sqrt(fabs(recent) < 1e-15 ? 0 : recent / fmax(potential, DBL_EPSILON));
            }
            else
            {
                ppml_transform_option(options, S, y, r, work);
                double step = 1;
                if (method == 2 || method == 4)
                {
                    double rr = ppml_quad_dot(r, r, S->weights, N);
                    /* Near the solution y'Wr is cancellation dominated. A
                     * unit projection remains stable; a tiny accelerated
                     * step can falsely report convergence before demeaning. */
                    if (rr > DBL_EPSILON * fmax(1, ppml_quad_dot(y, y, S->weights, N)) && iter % 10)
                        step = ppml_quad_dot(y, r, S->weights, N) / rr;
                    if (!isfinite(step) || step <= 0)
                        step = 1;
                }
                for (int i = 0; i < N; i++)
                    v[i] = y[i] - step * r[i];
                if (method == 3 && iter >= 6 && iter % 3 == 0)
                {
                    for (int i = 0; i < N; i++)
                        u[i] = v[i] - 2 * y[i] + last[i];
                    double factor = -ppml_quad_dot(r, u, S->weights, N) /
                                    fmax(ppml_quad_dot(u, u, S->weights, N), DBL_EPSILON);
                    for (int i = 0; i < N; i++)
                        v[i] += factor * r[i];
                }
                double ws = 0;
                for (int i = 0; i < N; i++)
                {
                    error += S->weights[i] * fabs(v[i] - y[i]) / (1 + fabs(y[i]));
                    ws += S->weights[i];
                    last[i] = y[i];
                    y[i] = v[i];
                }
                error /= ws;
            }
            if (error <= S->tolerance)
                break;
        }
        if (iter > S->maxiter)
        {
            rc = 430;
            break;
        }
        if (iter > max_iterations)
            max_iterations = iter;
        if (ppml_quad_dot(y, y, NULL, N) <= initial * fmin(1e-6, S->tolerance / 10))
            memset(y, 0, N * sizeof(double));
    }
    free(buffer);
    return (HDFE_SolveResult){rc, !rc, max_iterations};
}

HDFE_SolveResult ppml_partial(const PPML_Options *options, HDFE_State *S, double *data, int N,
                              int K, int threads)
{
    if (options->acceleration == 5)
        return ppml_lsmr(options, S, data, N, K);
    if (options->acceleration || options->transform)
        return ppml_map_options(options, S, data, N, K);
    int pool = options->pool > 0 ? options->pool : K;
    HDFE_SolveResult out = {0, 1, 0};
    for (int k = 0; k < K; k += pool)
    {
        HDFE_SolveResult part =
            partial_out_columns(S, data + (size_t)k * N, N, K - k < pool ? K - k : pool, threads);
        if (part.status)
            return part;
        if (part.iterations > out.iterations)
            out.iterations = part.iterations;
    }
    return out;
}

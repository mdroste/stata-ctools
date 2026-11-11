/*
 * hdfe_wild.c
 *
 * Wild-bootstrap engine shared by the HDFE estimators (see hdfe_wild.h and
 * docs/development/BOOTSTRAP_PLAN.md, sections 5.9 and 6.10). No SPI calls.
 *
 * Why the aggregates suffice. With z = sum_h v_h (S_u,h + delta S_r,h), the
 * replicate's estimates are b* = b~ + A z (X'W M_D = X'W, so the perturbed
 * residuals need no partialling for b*), and the score of vce() cluster g
 * is X_g'W_g u*_g = sum_{h in g} v_h S_h - X_g'W_g X_g A z, or, projected,
 * X_g'W_g M_D(u~ (.) v)_g - X_g'W_g X_g A z. The statistic needs a'score_g
 * only, and a'X_g'W_g X_g A z = d_g'z with d_g = A X_g'W_g rho_g.
 * Part of the ctools Stata plugin suite
 */

#include "hdfe_wild.h"
#include "../ctools_cancel.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define FIVE(m, j, b) ((m)->five + (((size_t)(j) * (size_t)(m)->B + (size_t)(b)) * 5))

void hdfe_wild_free(hdfe_wild_model *m)
{
    if (!m)
        return;
    free(m->g_start);
    free(m->h_boot);
    free(m->S_u);
    free(m->S_r);
    free(m->q_u);
    free(m->q_r);
    free(m->d);
    free(m->a);
    free(m->b_hat);
    free(m->se);
    free(m->row_g);
    free(m->row_c);
    free(m->w);
    free(m->u);
    free(m->rho);
    free(m->r);
    free(m->five);
    free(m->obs);
    memset(m, 0, sizeof(*m));
}

static ST_retcode observed(hdfe_wild_model *m, hdfe_wild_project_fn project, void *ctx);

static void *alloc_zero(size_t n, size_t size)
{
    return calloc(n > 0 ? n : 1, size);
}

/* ========================================================================
 * Setup
 * ======================================================================== */

ST_retcode hdfe_wild_setup(const hdfe_wild_input *in, ST_int B, hdfe_wild_model *m)
{
    memset(m, 0, sizeof(*m));
    if (!in || in->N < 1 || in->K < 1 || in->n_test < 1 || B < 1 || !in->X || !in->resid ||
        !in->A || !in->beta || !in->test)
        return 198;
    const ST_int N = in->N, K = in->K, T = in->n_test;
    const ST_int G = in->vce_id ? in->G : N, C = in->boot_id ? in->C : N;
    if (G < 1 || C < 1)
        return 198;
    for (ST_int j = 0; j < T; j++)
    {
        if (in->test[j] < 0 || in->test[j] >= K)
            return 198;
    }
    for (ST_int i = 0; i < N; i++)
    {
        ST_int g = in->vce_id ? in->vce_id[i] : i, c = in->boot_id ? in->boot_id[i] : i;
        if (g < 0 || g >= G || c < 0 || c >= C)
            return 198;
    }

    m->N = N;
    m->K = K;
    m->G = G;
    m->C = C;
    m->n_test = T;
    m->B = B;
    m->null_imposed = in->null_imposed != 0;
    m->projected = in->projected != 0;

    ST_retcode rc = 920;
    ST_int *row_g = malloc((size_t)N * sizeof(ST_int));
    ST_int *row_c = malloc((size_t)N * sizeof(ST_int));
    ST_int *order = malloc((size_t)N * sizeof(ST_int));
    ST_int *row_h = malloc((size_t)N * sizeof(ST_int));
    ST_int *pos = alloc_zero((size_t)G + 1, sizeof(ST_int));
    ST_int *stamp = malloc((size_t)C * sizeof(ST_int));
    ST_int *slot = malloc((size_t)C * sizeof(ST_int));
    ST_double *w = malloc((size_t)N * sizeof(ST_double));
    ST_double *rho = malloc((size_t)T * (size_t)N * sizeof(ST_double));
    ST_double *tmp = alloc_zero((size_t)G * (size_t)K, sizeof(ST_double));
    m->g_start = malloc(((size_t)G + 1) * sizeof(ST_int));
    m->h_boot = malloc((size_t)N * sizeof(ST_int));
    m->a = malloc((size_t)T * (size_t)K * sizeof(ST_double));
    m->b_hat = malloc((size_t)T * sizeof(ST_double));
    m->se = malloc((size_t)T * sizeof(ST_double));
    if (!row_g || !row_c || !order || !row_h || !pos || !stamp || !slot || !w || !rho || !tmp ||
        !m->g_start || !m->h_boot || !m->a || !m->b_hat || !m->se)
        goto done;

    for (ST_int i = 0; i < N; i++)
    {
        row_g[i] = in->vce_id ? in->vce_id[i] : i;
        row_c[i] = in->boot_id ? in->boot_id[i] : i;
        w[i] = in->w ? in->w[i] : 1.0;
    }

    /* Rows by vce() cluster (counting sort, stable), then the intersection
     * groups of each cluster in order of first appearance */
    for (ST_int i = 0; i < N; i++)
        pos[row_g[i] + 1]++;
    for (ST_int g = 0; g < G; g++)
        pos[g + 1] += pos[g];
    for (ST_int i = 0; i < N; i++)
        order[pos[row_g[i]]++] = i;
    for (ST_int c = 0; c < C; c++)
        stamp[c] = -1;
    ST_int H = 0, at = 0;
    for (ST_int g = 0; g < G; g++)
    {
        m->g_start[g] = H;
        for (; at < N && row_g[order[at]] == g; at++)
        {
            ST_int i = order[at], c = row_c[i];
            if (stamp[c] != g)
            {
                stamp[c] = g;
                slot[c] = H;
                m->h_boot[H++] = c;
            }
            row_h[i] = slot[c];
        }
    }
    m->g_start[G] = H;
    m->H = H;

    m->S_u = alloc_zero((size_t)H * (size_t)K, sizeof(ST_double));
    m->q_u = alloc_zero((size_t)T * (size_t)H, sizeof(ST_double));
    m->d = alloc_zero((size_t)T * (size_t)G * (size_t)K, sizeof(ST_double));
    m->five = alloc_zero((size_t)T * (size_t)B * 5, sizeof(ST_double));
    m->obs = alloc_zero((size_t)T * 5, sizeof(ST_double));
    if (!m->S_u || !m->q_u || !m->d || !m->five || !m->obs)
        goto done;
    if (m->null_imposed)
    {
        m->S_r = alloc_zero((size_t)T * (size_t)H * (size_t)K, sizeof(ST_double));
        m->q_r = alloc_zero((size_t)T * (size_t)H, sizeof(ST_double));
        if (!m->S_r || !m->q_r)
            goto done;
    }

    /* Tested coefficients: a = A e_k, rho = X a, r = rho / a_k */
    for (ST_int j = 0; j < T; j++)
    {
        ST_int k = in->test[j];
        ST_double *a = m->a + (size_t)j * K;
        for (ST_int l = 0; l < K; l++)
            a[l] = in->A[(size_t)k * K + l];
        m->b_hat[j] = in->beta[k];
        if (!(a[k] > 0) || !isfinite(a[k]))
        {
            rc = 198;
            goto done;
        }
        ST_double *rj = rho + (size_t)j * N;
        for (ST_int i = 0; i < N; i++)
            rj[i] = 0;
        for (ST_int l = 0; l < K; l++)
        {
            const ST_double *x = in->X + (size_t)l * N;
            for (ST_int i = 0; i < N; i++)
                rj[i] += x[i] * a[l];
        }
    }

    /* S_u,h = X_h'W_h u_h; S_r,h = X_h'W_h r_h */
    for (ST_int l = 0; l < K; l++)
    {
        const ST_double *x = in->X + (size_t)l * N;
        for (ST_int i = 0; i < N; i++)
            m->S_u[(size_t)row_h[i] * K + l] += w[i] * in->resid[i] * x[i];
        if (m->null_imposed)
        {
            for (ST_int j = 0; j < T; j++)
            {
                const ST_double akk = m->a[(size_t)j * K + in->test[j]];
                const ST_double *rj = rho + (size_t)j * N;
                ST_double *S = m->S_r + (size_t)j * H * K;
                for (ST_int i = 0; i < N; i++)
                    S[(size_t)row_h[i] * K + l] += w[i] * (rj[i] / akk) * x[i];
            }
        }
    }

    for (ST_int j = 0; j < T; j++)
    {
        const ST_double *a = m->a + (size_t)j * K;
        const ST_double *rj = rho + (size_t)j * N;
        for (ST_int h = 0; h < H; h++)
        {
            ST_double qu = 0, qr = 0;
            for (ST_int l = 0; l < K; l++)
            {
                qu += a[l] * m->S_u[(size_t)h * K + l];
                if (m->null_imposed)
                    qr += a[l] * m->S_r[((size_t)j * H + h) * K + l];
            }
            m->q_u[(size_t)j * H + h] = qu;
            if (m->null_imposed)
                m->q_r[(size_t)j * H + h] = qr;
        }

        /* d_g = A X_g'W_g rho_g */
        memset(tmp, 0, (size_t)G * (size_t)K * sizeof(ST_double));
        for (ST_int l = 0; l < K; l++)
        {
            const ST_double *x = in->X + (size_t)l * N;
            for (ST_int i = 0; i < N; i++)
                tmp[(size_t)row_g[i] * K + l] += w[i] * rj[i] * x[i];
        }
        for (ST_int g = 0; g < G; g++)
        {
            ST_double *dg = m->d + ((size_t)j * G + g) * K;
            const ST_double *tg = tmp + (size_t)g * K;
            for (ST_int l = 0; l < K; l++)
            {
                ST_double s = 0;
                for (ST_int l2 = 0; l2 < K; l2++)
                    s += in->A[(size_t)l2 * K + l] * tg[l2];
                dg[l] = s;
            }
        }

        /* The raw cluster-robust standard error of the full fit */
        ST_double ss = 0;
        for (ST_int g = 0; g < G; g++)
        {
            ST_double s = 0;
            for (ST_int h = m->g_start[g]; h < m->g_start[g + 1]; h++)
                s += m->q_u[(size_t)j * H + h];
            ss += s * s;
        }
        m->se[j] = sqrt(ss);
    }

    if (m->projected)
    {
        m->row_g = row_g;
        m->row_c = row_c;
        m->w = w;
        m->rho = rho;
        row_g = row_c = NULL;
        w = rho = NULL;
        m->u = malloc((size_t)N * sizeof(ST_double));
        if (!m->u)
            goto done;
        memcpy(m->u, in->resid, (size_t)N * sizeof(ST_double));
        if (m->null_imposed)
        {
            m->r = malloc((size_t)T * (size_t)N * sizeof(ST_double));
            if (!m->r)
                goto done;
            for (ST_int j = 0; j < T; j++)
            {
                const ST_double akk = m->a[(size_t)j * K + in->test[j]];
                for (ST_int i = 0; i < N; i++)
                    m->r[(size_t)j * N + i] = m->rho[(size_t)j * N + i] / akk;
            }
        }
    }
    rc = m->projected ? 0 : observed(m, NULL, NULL);

done:
    free(row_g);
    free(row_c);
    free(order);
    free(row_h);
    free(pos);
    free(stamp);
    free(slot);
    free(w);
    free(rho);
    free(tmp);
    if (rc)
        hdfe_wild_free(m);
    return rc;
}

/* ========================================================================
 * Replicates
 * ======================================================================== */

ST_int hdfe_wild_enumerable(ST_int C, ST_int B)
{
    if (C < 1 || C > 30 || B < 1)
        return 0;
    return (double)C * log(2.0) < log((double)B) + 1e-6 ? (ST_int)1 << C : 0;
}

void hdfe_wild_enum_weights(ST_int C, uint64_t j, ST_double *v)
{
    for (ST_int c = 0; c < C; c++)
        v[c] = ((j >> (C - 1 - c)) & 1u) ? 1.0 : -1.0;
}

/* Scratch per replicate: v_h (H), z_u and z_r (K each), the vce() cluster
 * sums of the projected path (2 G), and enumerated weights (C). */
static size_t scratch_len(const hdfe_wild_model *m)
{
    return (size_t)m->H + 2 * (size_t)m->K + 2 * (size_t)m->G + (size_t)m->C;
}

/*
 * One replicate from its bootstrap-cluster weights, written to out + j
 * stride for each test j. cols: NULL on the
 * aggregate path; on the projected path the replicate's partialled columns
 * M_D(u (.) v) and, under WCR, M_D(r_j (.) v) for each test.
 */
static void replicate(const hdfe_wild_model *m, ST_double *out, size_t stride,
                      const ST_double *vrow, const ST_double *cols, ST_double *scr)
{
    const ST_int K = m->K, H = m->H, G = m->G, N = m->N, T = m->n_test;
    ST_double *vh = scr, *zu = vh + H, *zr = zu + K, *acc0 = zr + K, *acc1 = acc0 + G;

    for (ST_int h = 0; h < H; h++)
        vh[h] = vrow[m->h_boot[h]];
    for (ST_int l = 0; l < K; l++)
        zu[l] = 0;
    for (ST_int h = 0; h < H; h++)
    {
        const ST_double *S = m->S_u + (size_t)h * K;
        for (ST_int l = 0; l < K; l++)
            zu[l] += vh[h] * S[l];
    }

    for (ST_int j = 0; j < T; j++)
    {
        const ST_double *a = m->a + (size_t)j * K;
        ST_double n0 = 0, n1 = 0;
        for (ST_int l = 0; l < K; l++)
            n0 += a[l] * zu[l];
        for (ST_int l = 0; l < K; l++)
            zr[l] = 0;
        if (m->null_imposed)
        {
            const ST_double *Sr = m->S_r + (size_t)j * H * K;
            for (ST_int h = 0; h < H; h++)
                for (ST_int l = 0; l < K; l++)
                    zr[l] += vh[h] * Sr[(size_t)h * K + l];
            for (ST_int l = 0; l < K; l++)
                n1 += a[l] * zr[l];
        }

        if (cols)
        {
            const ST_double *P = cols;
            const ST_double *R = m->null_imposed ? cols + (size_t)(1 + j) * N : NULL;
            const ST_double *rj = m->rho + (size_t)j * N;
            for (ST_int g = 0; g < G; g++)
                acc0[g] = acc1[g] = 0;
            for (ST_int i = 0; i < N; i++)
            {
                const ST_double wr = m->w[i] * rj[i];
                acc0[m->row_g[i]] += wr * P[i];
                if (R)
                    acc1[m->row_g[i]] += wr * R[i];
            }
        }

        ST_double D00 = 0, D01 = 0, D11 = 0;
        for (ST_int g = 0; g < G; g++)
        {
            const ST_double *dg = m->d + ((size_t)j * G + g) * K;
            ST_double e0 = 0, e1 = 0;
            if (cols)
            {
                e0 = acc0[g];
                e1 = acc1[g];
            }
            else
            {
                for (ST_int h = m->g_start[g]; h < m->g_start[g + 1]; h++)
                {
                    e0 += vh[h] * m->q_u[(size_t)j * H + h];
                    if (m->null_imposed)
                        e1 += vh[h] * m->q_r[(size_t)j * H + h];
                }
            }
            for (ST_int l = 0; l < K; l++)
            {
                e0 -= dg[l] * zu[l];
                e1 -= dg[l] * zr[l];
            }
            D00 += e0 * e0;
            D01 += e0 * e1;
            D11 += e1 * e1;
        }
        ST_double *f = out + (size_t)j * stride;
        f[0] = n0;
        f[1] = n1;
        f[2] = D00;
        f[3] = D01;
        f[4] = D11;
    }
}

/* The observed statistic: the replicate with v = 1 (cols as for replicates
 * on the projected path, NULL otherwise). */
static ST_retcode observed(hdfe_wild_model *m, hdfe_wild_project_fn project, void *ctx)
{
    const size_t slen = scratch_len(m);
    const ST_int per = m->null_imposed ? 1 + m->n_test : 1;
    ST_double *scr = malloc(slen * sizeof(ST_double));
    ST_double *ones = malloc((size_t)m->C * sizeof(ST_double));
    ST_double *cols = m->projected ? malloc((size_t)m->N * per * sizeof(ST_double)) : NULL;
    ST_retcode rc = 0;
    if (!scr || !ones || (m->projected && !cols))
        rc = 920;
    if (!rc)
    {
        for (ST_int c = 0; c < m->C; c++)
            ones[c] = 1.0;
        if (cols)
        {
            memcpy(cols, m->u, (size_t)m->N * sizeof(ST_double));
            if (m->null_imposed)
                memcpy(cols + m->N, m->r, (size_t)m->N * m->n_test * sizeof(ST_double));
            rc = project(ctx, cols, per);
        }
    }
    if (!rc)
    {
        replicate(m, m->obs, 5, ones, cols, scr);
        m->obs_ready = 1;
    }
    free(scr);
    free(ones);
    free(cols);
    return rc;
}

/* The weights of replicate b: row b - first of v, or enumerated. */
static const ST_double *weights_of(const hdfe_wild_model *m, const ST_double *v, ST_int first,
                                   ST_int b, ST_double *buf)
{
    if (v)
        return v + (size_t)(b - first) * (size_t)m->C;
    hdfe_wild_enum_weights(m->C, (uint64_t)b, buf);
    return buf;
}

ST_retcode hdfe_wild_run(hdfe_wild_model *m, const ST_double *v, ST_int first, ST_int count,
                         ST_int threads, hdfe_wild_project_fn project, void *ctx,
                         size_t budget)
{
    if (!m || !m->five || first != m->done || count < 0 || count > m->B - first)
        return 198;
    if (count == 0)
        return 0;
    if (!v && hdfe_wild_enumerable(m->C, m->B) != m->B)
        return 198;
    if (m->projected && !project)
        return 198;
    const ST_int T_thr = threads < 1 ? 1 : threads;
    const size_t slen = scratch_len(m);
    volatile int fail = 0;

    if (!m->projected)
    {
        ctools_cancel_team team = {0};
        const ST_int nthr = T_thr < count ? T_thr : count;
#ifdef _OPENMP
#pragma omp parallel num_threads(nthr) if (nthr > 1)
#endif
        {
            ST_double *scr = malloc(slen * sizeof(ST_double));
            if (!scr)
            {
#ifdef _OPENMP
#pragma omp critical(hdfe_wild_fail)
#endif
                fail = 1;
            }
#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif
            for (ST_int b = first; b < first + count; b++)
            {
                if (ctools_cancel_requested()) continue;
                if (scr)
                {
                    ST_double *vbuf = scr + (slen - (size_t)m->C);
                    replicate(m, FIVE(m, 0, b), (size_t)m->B * 5,
                              weights_of(m, v, first, b, vbuf), NULL, scr);
                }
            }
            free(scr);
            ctools_cancel_wait(&team);
        }
        (void)nthr;
        if (ctools_cancel_requested()) return 1;
        if (fail)
            return 920;
        m->done += count;
        return 0;
    }

    /* Projected path: batches of replicates, about HDFE_WILD_BATCH_BYTES of
     * columns at a time; one replicate's columns may take more, up to the
     * budget (the memory the caller allows) */
    if (!m->obs_ready)
    {
        ST_retcode orc = observed(m, project, ctx);
        if (orc)
            return orc;
    }
    const ST_int per = m->null_imposed ? 1 + m->n_test : 1;
    const size_t col_bytes = (size_t)m->N * (size_t)per * sizeof(ST_double);
    if (budget < col_bytes)
        return 909;
    const size_t room = budget < HDFE_WILD_BATCH_BYTES ? budget : HDFE_WILD_BATCH_BYTES;
    size_t fit = room / col_bytes;
    if (fit < 1)
        fit = 1;
    ST_int nb = fit < (size_t)count ? (ST_int)fit : count;
    ST_double *cols = malloc((size_t)nb * col_bytes);
    ST_double *vbufs = v ? NULL : malloc((size_t)nb * (size_t)m->C * sizeof(ST_double));
    if (!cols || (!v && !vbufs))
    {
        free(cols);
        free(vbufs);
        return 920;
    }
    ST_retcode rc = 0;
    for (ST_int start = first; start < first + count && !rc; start += nb)
    {
        if (ctools_cancel_requested()) { rc = 1; break; }
        ctools_cancel_team team = {0};
        const ST_int n = first + count - start < nb ? first + count - start : nb;
        for (ST_int bb = 0; bb < n; bb++)
        {
            const ST_double *vr =
                weights_of(m, v, first, start + bb, vbufs ? vbufs + (size_t)bb * m->C : NULL);
            ST_double *c0 = cols + (size_t)bb * per * m->N;
            for (ST_int i = 0; i < m->N; i++)
                c0[i] = m->u[i] * vr[m->row_c[i]];
            for (ST_int j = 0; m->null_imposed && j < m->n_test; j++)
            {
                ST_double *cj = c0 + (size_t)(1 + j) * m->N;
                const ST_double *rj = m->r + (size_t)j * m->N;
                for (ST_int i = 0; i < m->N; i++)
                    cj[i] = rj[i] * vr[m->row_c[i]];
            }
        }
        rc = project(ctx, cols, n * per);
        if (rc)
            break;
        const ST_int nthr = T_thr < n ? T_thr : n;
#ifdef _OPENMP
#pragma omp parallel num_threads(nthr) if (nthr > 1)
#endif
        {
            ST_double *scr = malloc(slen * sizeof(ST_double));
            if (!scr)
            {
#ifdef _OPENMP
#pragma omp critical(hdfe_wild_fail)
#endif
                fail = 1;
            }
#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif
            for (ST_int bb = 0; bb < n; bb++)
            {
                if (ctools_cancel_requested()) continue;
                if (scr)
                {
                    ST_double *vbuf = scr + (slen - (size_t)m->C);
                    const ST_double *vr = weights_of(m, v, first, start + bb, vbuf);
                    replicate(m, FIVE(m, 0, start + bb), (size_t)m->B * 5, vr,
                              cols + (size_t)bb * per * m->N, scr);
                }
            }
            free(scr);
            ctools_cancel_wait(&team);
        }
        (void)nthr;
        if (ctools_cancel_requested())
            rc = 1;
        else if (fail)
            rc = 920;
    }
    free(cols);
    free(vbufs);
    if (!rc)
        m->done += count;
    return rc;
}

/* ========================================================================
 * Statistics and p-values
 * ======================================================================== */

double hdfe_wild_t(const hdfe_wild_model *m, ST_int j, double b0)
{
    if (!m->obs_ready)
        return NAN;
    const ST_double *f = m->obs + (size_t)j * 5;
    const double delta = m->b_hat[j] - b0;
    if (!m->null_imposed)
        return f[2] > 0 ? delta / sqrt(f[2]) : NAN;
    const double num = f[0] + delta * f[1];
    const double den2 = f[2] + delta * (2.0 * f[3] + delta * f[4]);
    if (!(den2 > 0) || !isfinite(den2) || !isfinite(num))
        return NAN;
    return num / sqrt(den2);
}

double hdfe_wild_tstar(const hdfe_wild_model *m, ST_int j, ST_int b, double b0)
{
    const ST_double *f = FIVE(m, j, b);
    const double delta = m->null_imposed ? m->b_hat[j] - b0 : 0.0;
    const double num = f[0] + delta * f[1];
    const double den2 = f[2] + delta * (2.0 * f[3] + delta * f[4]);
    if (!(den2 > 0) || !isfinite(den2) || !isfinite(num))
        return NAN;
    return num / sqrt(den2);
}

static double p_from_counts(ST_int F, ST_int lt, ST_int gt, ST_int absgt, int ptype)
{
    if (F == 0)
        return NAN;
    switch (ptype)
    {
    case HDFE_WILD_EQUALTAIL:
        return 2.0 * (double)(lt < gt ? lt : gt) / F;
    case HDFE_WILD_LOWER:
        return (double)lt / F;
    case HDFE_WILD_UPPER:
        return (double)gt / F;
    default:
        return (double)absgt / F;
    }
}

double hdfe_wild_pvalue_t(const ST_double *t, ST_int B, int ptype, ST_int *feasible)
{
    if (feasible)
        *feasible = 0;
    if (isnan(t[0]))
        return NAN;
    ST_int F = 0, lt = 0, gt = 0, absgt = 0;
    for (ST_int b = 1; b <= B; b++)
    {
        const double ts = t[b];
        if (isnan(ts))
            continue;
        F++;
        lt += ts < t[0];
        gt += ts > t[0];
        absgt += fabs(ts) > fabs(t[0]);
    }
    if (feasible)
        *feasible = F;
    return p_from_counts(F, lt, gt, absgt, ptype);
}

double hdfe_wild_pvalue(const hdfe_wild_model *m, ST_int j, double b0, int ptype,
                        ST_int *feasible)
{
    if (feasible)
        *feasible = 0;
    if (m->done != m->B || j < 0 || j >= m->n_test)
        return NAN;
    const double t = hdfe_wild_t(m, j, b0);
    if (isnan(t))
        return NAN;
    ST_int F = 0, lt = 0, gt = 0, absgt = 0;
    for (ST_int b = 0; b < m->B; b++)
    {
        const double ts = hdfe_wild_tstar(m, j, b, b0);
        if (isnan(ts))
            continue;
        F++;
        lt += ts < t;
        gt += ts > t;
        absgt += fabs(ts) > fabs(t);
    }
    if (feasible)
        *feasible = F;
    return p_from_counts(F, lt, gt, absgt, ptype);
}

/* ========================================================================
 * Confidence sets by test inversion
 * ======================================================================== */

/* x with P(Z > x) = a, by bisection (enough for a grid width) */
static double upper_normal_quantile(double a)
{
    double lo = -40, hi = 40;
    for (int it = 0; it < 200; it++)
    {
        double mid = 0.5 * (lo + hi);
        if (0.5 * erfc(mid / sqrt(2.0)) > a)
            lo = mid;
        else
            hi = mid;
    }
    return 0.5 * (lo + hi);
}

/* p(x) > alpha. A NaN p (the evaluation failed) counts as not rejected, so
 * that a failure can only widen the set, and is counted in *failed. */
static int accept_p(double px, double alpha, int *failed)
{
    if (isnan(px))
    {
        ++*failed;
        return 1;
    }
    return px > alpha;
}

static int accepted(hdfe_wild_pfn p, void *ctx, double x, double alpha, int *failed)
{
    return accept_p(p(ctx, x), alpha, failed);
}

/* The boundary between an accepted and a rejected point */
static double bisect(hdfe_wild_pfn p, void *ctx, double acc, double rej, double alpha,
                     double tol, int *failed)
{
    for (int it = 0; it < 200 && fabs(rej - acc) > tol; it++)
    {
        const double mid = 0.5 * (acc + rej);
        if (accepted(p, ctx, mid, alpha, failed))
            acc = mid;
        else
            rej = mid;
    }
    return 0.5 * (acc + rej);
}

/* Doublings of the outward search: from scale, 1000 scale is reached in 10;
 * the cap only guards against a non-finite scale or center */
#define WILD_OUTWARD_MAX 64

/* Outward from an accepted point x_acc (side +1 or -1): the bound, or NaN
 * when none lies within the limit. The stop test is on the step itself,
 * which is set to exactly limit when capped: center + limit - center can
 * round below limit, and testing that would query the capped point forever. */
static double outward(hdfe_wild_pfn p, void *ctx, double center, double x_acc, int side,
                      double scale, double alpha, double tol, int *failed)
{
    const double limit = 1000.0 * scale;
    double dist = fabs(x_acc - center);
    if (dist < scale)
        dist = scale;
    for (int it = 0; it < WILD_OUTWARD_MAX; it++)
    {
        if (dist >= limit)
            return NAN;
        double next = 2.0 * dist;
        if (next > limit)
            next = limit;
        const double x = center + side * next;
        if (!accepted(p, ctx, x, alpha, failed))
            return bisect(p, ctx, x_acc, x, alpha, tol, failed);
        x_acc = x;
        dist = next;
    }
    return NAN;
}

void hdfe_wild_invert(hdfe_wild_pfn p, void *ctx, double center, double scale, double alpha,
                      hdfe_wild_ci *ci)
{
    hdfe_wild_invert_ex(p, ctx, center, scale, alpha, 32, 1e-6, ci);
}

void hdfe_wild_invert_ex(hdfe_wild_pfn p, void *ctx, double center, double scale, double alpha,
                         int half, double rel_tol, hdfe_wild_ci *ci)
{
    enum
    {
        HALF_MAX = 32,
        NPT_MAX = 2 * HALF_MAX + 1
    };
    const int HALF = half < 1 ? 1 : (half > HALF_MAX ? HALF_MAX : half);
    const int NPT = 2 * HALF + 1;
    memset(ci, 0, sizeof(*ci));
    ci->lo = ci->hi = NAN;
    if (!isfinite(center) || !(scale > 0) || !isfinite(scale) || !(alpha > 0) || !(alpha < 1))
    {
        ci->empty = 1;
        return;
    }
    const double z = upper_normal_quantile(alpha / 2);
    const double span = 2.0 * (z > 1 ? z : 1) * scale;
    const double tol = rel_tol * scale;
    double xs[NPT_MAX];
    int acc[NPT_MAX];
    int first = -1, last = -1;
    for (int i = 0; i < NPT; i++)
    {
        xs[i] = center + span * (double)(i - HALF) / HALF;
        acc[i] = accepted(p, ctx, xs[i], alpha, &ci->failed);
        if (acc[i])
        {
            if (first < 0)
                first = i;
            last = i;
        }
    }
    if (first < 0)
    {
        ci->empty = 1;
        return;
    }
    for (int i = first; i <= last; i++)
    {
        if (!acc[i])
            ci->disjoint = 1;
    }
    ci->hi = last < NPT - 1
                 ? bisect(p, ctx, xs[last], xs[last + 1], alpha, tol, &ci->failed)
                 : outward(p, ctx, center, xs[last], 1, scale, alpha, tol, &ci->failed);
    ci->lo = first > 0
                 ? bisect(p, ctx, xs[first], xs[first - 1], alpha, tol, &ci->failed)
                 : outward(p, ctx, center, xs[first], -1, scale, alpha, tol, &ci->failed);
    ci->open_hi = isnan(ci->hi);
    ci->open_lo = isnan(ci->lo);
}

/* One bound's search for hdfe_wild_invert_batch: outward() and bisect() as
 * a state machine that yields one query at a time */
typedef struct
{
    int active, bisecting, side, it, from, to;
    double acc, rej;    /* bisection bracket */
    double x_acc, dist; /* outward: the last accepted point, its distance */
    double query, next, result;
} bound_search;

static void bound_start(bound_search *b, int side, double center, double scale, double x_acc,
                        const double *x_rej, int from, int to)
{
    memset(b, 0, sizeof(*b));
    b->active = 1;
    b->side = side;
    b->from = from;
    b->to = to;
    if (x_rej)
    {
        b->bisecting = 1;
        b->acc = x_acc;
        b->rej = *x_rej;
    }
    else
    {
        b->x_acc = x_acc;
        b->dist = fabs(x_acc - center);
        if (b->dist < scale)
            b->dist = scale;
    }
}

/* Sets b->query and returns 1, or ends the search (b->result) and returns 0 */
static int bound_query(bound_search *b, double center, double scale, double tol)
{
    if (b->bisecting)
    {
        if (b->it >= 200 || !(fabs(b->rej - b->acc) > tol))
        {
            b->result = 0.5 * (b->acc + b->rej);
            b->active = 0;
            return 0;
        }
        b->query = 0.5 * (b->acc + b->rej);
        return 1;
    }
    /* The stop test is on the step, as in outward() */
    const double limit = 1000.0 * scale;
    if (b->dist >= limit || b->it >= WILD_OUTWARD_MAX)
    {
        b->result = NAN;
        b->active = 0;
        return 0;
    }
    b->next = 2.0 * b->dist;
    if (b->next > limit)
        b->next = limit;
    b->query = center + b->side * b->next;
    return 1;
}

static void bound_update(bound_search *b, int accepted)
{
    b->from = b->to;
    if (b->bisecting)
    {
        if (accepted)
            b->acc = b->query;
        else
            b->rej = b->query;
        b->it++;
    }
    else if (!accepted)
    {
        b->bisecting = 1;
        b->acc = b->x_acc;
        b->rej = b->query;
        b->it = 0;
    }
    else
    {
        b->x_acc = b->query;
        b->dist = b->next;
        b->it++;
    }
}

void hdfe_wild_invert_batch(hdfe_wild_pbatch_fn p, void *ctx, double center, double scale,
                            double alpha, int half, double rel_tol, hdfe_wild_ci *ci)
{
    enum
    {
        HALF_MAX = 32,
        NPT_MAX = 2 * HALF_MAX + 1
    };
    const int HALF = half < 1 ? 1 : (half > HALF_MAX ? HALF_MAX : half);
    const int NPT = 2 * HALF + 1;
    memset(ci, 0, sizeof(*ci));
    ci->lo = ci->hi = NAN;
    if (!isfinite(center) || !(scale > 0) || !isfinite(scale) || !(alpha > 0) || !(alpha < 1))
    {
        ci->empty = 1;
        return;
    }
    const double z = upper_normal_quantile(alpha / 2);
    const double span = 2.0 * (z > 1 ? z : 1) * scale;
    const double tol = rel_tol * scale;
    double xs[NPT_MAX], pv[NPT_MAX];
    int acc[NPT_MAX], from[NPT_MAX], to[NPT_MAX];
    int first = -1, last = -1;
    for (int i = 0; i < NPT; i++)
    {
        xs[i] = center + span * (double)(i - HALF) / HALF;
        from[i] = -1;
        to[i] = i;
    }
    p(ctx, xs, from, to, NPT, pv);
    for (int i = 0; i < NPT; i++)
    {
        acc[i] = accept_p(pv[i], alpha, &ci->failed);
        if (acc[i])
        {
            if (first < 0)
                first = i;
            last = i;
        }
    }
    if (first < 0)
    {
        ci->empty = 1;
        return;
    }
    for (int i = first; i <= last; i++)
    {
        if (!acc[i])
            ci->disjoint = 1;
    }
    bound_search b[2]; /* 0: upper, 1: lower */
    bound_start(&b[0], 1, center, scale, xs[last], last < NPT - 1 ? &xs[last + 1] : NULL, last,
                NPT);
    bound_start(&b[1], -1, center, scale, xs[first], first > 0 ? &xs[first - 1] : NULL, first,
                NPT + 1);
    for (;;)
    {
        double q[2], pq[2];
        int qf[2], qt[2], who[2], n = 0;
        for (int s = 0; s < 2; s++)
        {
            if (b[s].active && bound_query(&b[s], center, scale, tol))
            {
                q[n] = b[s].query;
                qf[n] = b[s].from;
                qt[n] = b[s].to;
                who[n++] = s;
            }
        }
        if (n == 0)
            break;
        p(ctx, q, qf, qt, n, pq);
        for (int k = 0; k < n; k++)
            bound_update(&b[who[k]], accept_p(pq[k], alpha, &ci->failed));
    }
    ci->hi = b[0].result;
    ci->lo = b[1].result;
    ci->open_hi = isnan(ci->hi);
    ci->open_lo = isnan(ci->lo);
}

/* ========================================================================
 * Score bootstrap
 * ======================================================================== */

void hdfe_wild_score_free(hdfe_wild_score_model *m)
{
    if (!m)
        return;
    free(m->g_start);
    free(m->h_boot);
    free(m->row_h);
    free(m->share);
    free(m->v);
    memset(m, 0, sizeof(*m));
}

ST_retcode hdfe_wild_score_setup(ST_int N, const ST_int *vce_id, ST_int G, const ST_int *boot_id,
                                 ST_int C, const ST_double *share_w, ST_int B, int enumerate,
                                 hdfe_wild_score_model *m)
{
    memset(m, 0, sizeof(*m));
    if (N < 1 || B < 1)
        return 198;
    if (!vce_id)
        G = N;
    if (!boot_id)
        C = N;
    if (G < 1 || C < 1 || (enumerate && hdfe_wild_enumerable(C, B) != B))
        return 198;
    for (ST_int i = 0; i < N; i++)
    {
        const ST_int g = vce_id ? vce_id[i] : i, c = boot_id ? boot_id[i] : i;
        if (g < 0 || g >= G || c < 0 || c >= C)
            return 198;
    }
    m->N = N;
    m->G = G;
    m->C = C;
    m->B = B;
    ST_int *pos = calloc((size_t)G + 1, sizeof(ST_int));
    ST_int *order = malloc((size_t)N * sizeof(ST_int));
    ST_int *stamp = malloc((size_t)C * sizeof(ST_int));
    ST_int *slot = malloc((size_t)C * sizeof(ST_int));
    m->g_start = malloc(((size_t)G + 1) * sizeof(ST_int));
    m->h_boot = malloc((size_t)N * sizeof(ST_int));
    m->row_h = malloc((size_t)N * sizeof(ST_int));
    m->share = calloc((size_t)G, sizeof(ST_double));
    m->v = enumerate ? NULL : malloc((size_t)B * (size_t)C * sizeof(ST_double));
    m->got = enumerate ? B : 0;
    ST_retcode rc = 0;
    if (!pos || !order || !stamp || !slot || !m->g_start || !m->h_boot || !m->row_h ||
        !m->share || (!enumerate && !m->v))
        rc = 920;
    if (!rc)
    {
        double total = 0.0;
        for (ST_int i = 0; i < N; i++)
        {
            const ST_int g = vce_id ? vce_id[i] : i;
            const double sw = share_w ? share_w[i] : 1.0;
            pos[g + 1]++;
            m->share[g] += sw;
            total += sw;
        }
        for (ST_int g = 0; g < G; g++)
        {
            pos[g + 1] += pos[g];
            m->share[g] /= total;
        }
        for (ST_int i = 0; i < N; i++)
            order[pos[vce_id ? vce_id[i] : i]++] = i;
        for (ST_int c = 0; c < C; c++)
            stamp[c] = -1;
        ST_int H = 0, at = 0;
        for (ST_int g = 0; g < G; g++)
        {
            m->g_start[g] = H;
            for (; at < N && (vce_id ? vce_id[order[at]] : order[at]) == g; at++)
            {
                const ST_int i = order[at], c = boot_id ? boot_id[i] : i;
                if (stamp[c] != g)
                {
                    stamp[c] = g;
                    slot[c] = H;
                    m->h_boot[H++] = c;
                }
                m->row_h[i] = slot[c];
            }
        }
        m->g_start[G] = H;
        m->H = H;
    }
    free(pos);
    free(order);
    free(stamp);
    free(slot);
    if (rc)
        hdfe_wild_score_free(m);
    return rc;
}

ST_retcode hdfe_wild_score_weights(hdfe_wild_score_model *m, const ST_double *v, ST_int first,
                                   ST_int count)
{
    if (!m->v || first != m->got || count < 0 || first + count > m->B)
        return 198;
    memcpy(m->v + (size_t)first * m->C, v, (size_t)count * (size_t)m->C * sizeof(ST_double));
    m->got += count;
    return 0;
}

/* One replicate's statistic from group sums E (H) and weights vrow (C) */
static double score_stat(const hdfe_wild_score_model *m, const ST_double *E,
                         const ST_double *vrow, ST_double *vE)
{
    double num = 0.0;
    for (ST_int h = 0; h < m->H; h++)
    {
        vE[h] = vrow[m->h_boot[h]] * E[h];
        num += vE[h];
    }
    double den = 0.0;
    for (ST_int g = 0; g < m->G; g++)
    {
        double s = 0.0;
        for (ST_int h = m->g_start[g]; h < m->g_start[g + 1]; h++)
            s += vE[h];
        s -= m->share[g] * num;
        den += s * s;
    }
    if (!(den > 0) || !isfinite(den) || !isfinite(num))
        return NAN;
    return num / sqrt(den);
}

ST_retcode hdfe_wild_score_t(const hdfe_wild_score_model *m, const ST_double *e, ST_int threads,
                             ST_double *t)
{
    /* every replicate's weights must have arrived */
    if (m->got != m->B)
        return 198;
    ST_double *E = calloc((size_t)m->H, sizeof(ST_double));
    ST_double *ones = malloc((size_t)m->C * sizeof(ST_double));
    ST_double *vE = malloc((size_t)m->H * sizeof(ST_double));
    if (!E || !ones || !vE)
    {
        free(E);
        free(ones);
        free(vE);
        return 920;
    }
    for (ST_int i = 0; i < m->N; i++)
        E[m->row_h[i]] += e[i];
    for (ST_int c = 0; c < m->C; c++)
        ones[c] = 1.0;
    t[0] = score_stat(m, E, ones, vE);
    free(ones);
    free(vE);
    const ST_int T_thr = threads < 1 ? 1 : (threads < m->B ? threads : m->B);
    volatile int fail = 0;
    ctools_cancel_team team = {0};
#ifdef _OPENMP
#pragma omp parallel num_threads(T_thr) if (T_thr > 1)
#endif
    {
        ST_double *buf = malloc(((size_t)m->H + (size_t)m->C) * sizeof(ST_double));
        if (!buf)
        {
#ifdef _OPENMP
#pragma omp critical(hdfe_wild_fail)
#endif
            fail = 1;
        }
#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif
        for (ST_int b = 0; b < m->B; b++)
        {
            if (ctools_cancel_requested()) continue;
            if (!buf)
                continue;
            const ST_double *vrow;
            if (m->v)
                vrow = m->v + (size_t)b * m->C;
            else
            {
                hdfe_wild_enum_weights(m->C, (uint64_t)b, buf + m->H);
                vrow = buf + m->H;
            }
            t[1 + b] = score_stat(m, E, vrow, buf);
        }
        free(buf);
        ctools_cancel_wait(&team);
    }
    (void)T_thr;
    free(E);
    return ctools_cancel_requested() ? 1 : fail ? 920 : 0;
}

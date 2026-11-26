/*
 * hdfe_coef_multi.c
 *
 * Coefficient-space projection on three or more fixed effects (hdfe_coef.h)
 * Part of the ctools Stata plugin suite
 *
 * B, the effect with the most levels, is eliminated; the coefficients of
 * the other effects are stacked in one vector beta and solve
 *
 *     S beta = Dk' W M_B x,     S = Dk' W M_B Dk.
 *
 * S is applied through the rows grouped by B: within a level of B, every
 * row's sum of coefficients less the level's weighted mean of those sums,
 * scattered back to the row's kept levels. A level of B with one row
 * contributes nothing and is skipped.
 *
 * The preconditioner is one symmetric block Gauss-Seidel sweep over the kept
 * effects on Dk' W Dk (its diagonal blocks are the level counts, the others
 * the cross-tabulations of pairs of kept effects). In exact arithmetic that
 * is the conjugate-gradient iteration of the row-space solver for the sweep
 * that starts with B. When the kept effects have few levels in total, S is
 * assembled and factored instead.
 */

#include "hdfe_coef_internal.h"

/* As HDFE_COEF_DENSE_MAX and HDFE_COEF_DENSE_SWEEPS in hdfe_coef.c */
#define HDFE_COEF_MULTI_DENSE_MAX 1000
#define HDFE_COEF_MULTI_DENSE_SWEEPS 20.0
/* Rows of a level of B whose sums of coefficients are kept between the two
 * loops of the sweep rather than gathered twice */
#define HDFE_COEF_MULTI_ROWBUF 64
#define HDFE_COEF_MULTI_CHUNKS 32
#define HDFE_COEF_MULTI_CHUNK_ROWS 16384
/* As HDFE_COEF_PIVOT_TOL in hdfe_coef_core.c: the rank deficiency of S is
 * found by the factorization here */

typedef struct
{
    hdfe_coef_rows D; /* the data rows */
    int nk;
    ST_int m;
    ST_int nlev[HDFE_COEF_MAX_EFFECTS];
    ST_int *rowptr;                          /* nB + 1: rows grouped by B */
    const ST_int *gk[HDFE_COEF_MAX_EFFECTS]; /* kept levels in grouped order */
    ST_int *gk_own[HDFE_COEF_MAX_EFFECTS];   /* storage unless the rows arrive grouped */
    const ST_double *gw;                     /* weights in grouped order, NULL: one */
    ST_double *gw_own;
    ST_double *invB;      /* nB */
    ST_double *invk;      /* m: 1 / weighted count of each kept level */
    ST_int nact, nnz_act; /* levels of B with two or more rows */
    ST_int *act;
    int nch;
    ST_int *ch_act;
    /* cross-tabulations of the kept effects: rows a, columns b (a != b) */
    hdfe_crosstab X[HDFE_COEF_MAX_EFFECTS][HDFE_COEF_MAX_EFFECTS];
} hcm_system;

typedef struct
{
    const hcm_system *M;
    int K;
    ST_double *buf; /* nch * m * K */
} hcm_sweep;

static void hcm_system_free(hcm_system *M)
{
    free(M->rowptr);
    free(M->gw_own);
    free(M->invB);
    free(M->invk);
    free(M->act);
    free(M->ch_act);
    for (int a = 0; a < HDFE_COEF_MAX_EFFECTS; a++)
    {
        free(M->gk_own[a]);
        for (int b = 0; b < HDFE_COEF_MAX_EFFECTS; b++)
            hdfe_crosstab_free(&M->X[a][b]);
    }
}

/* T = A with rows and columns exchanged */
static ST_retcode hcm_transpose(const hdfe_crosstab *A, hdfe_crosstab *T)
{
    memset(T, 0, sizeof(*T));
    T->nrow = A->ncol;
    T->ncol = A->nrow;
    T->nnz = A->nnz;
    ST_int *dest = NULL;
    int sorted = 0;
    ST_retcode rc =
        hdfe_coef_group_rows(A->nnz, A->col, A->ncol, NULL, 0, &T->rowptr, &dest, &sorted, NULL);
    if (rc)
        return rc;
    T->col_own = (ST_int *)malloc((size_t)A->nnz * sizeof(ST_int));
    T->val_own = A->val ? (ST_double *)malloc((size_t)A->nnz * sizeof(ST_double)) : NULL;
    if (!T->col_own || (A->val && !T->val_own))
    {
        free(dest);
        hdfe_crosstab_free(T);
        return 920;
    }
    ST_int *tcol = T->col_own;
    ST_double *tval = T->val_own;
#pragma omp parallel for schedule(dynamic, 1024) if (A->nnz > HDFE_COEF_PARALLEL_MIN)
    for (ST_int r = 0; r < A->nrow; r++)
    {
        for (ST_int e = A->rowptr[r]; e < A->rowptr[r + 1]; e++)
        {
            const ST_int at = dest ? dest[e] : e;
            tcol[at] = r + 1;
            if (tval)
                tval[at] = A->val[e];
        }
    }
    free(dest);
    T->col = T->col_own;
    T->val = T->val_own;
    return 0;
}

static ST_retcode hcm_system_build(hcm_system *M, const HDFE_State *S, int need_crosstabs)
{
    memset(M, 0, sizeof(*M));
    const ST_int N = S->N;
    const int G = (int)S->G;
    int b = 0;
    for (int g = 1; g < G; g++)
    {
        if (S->factors[g].num_levels > S->factors[b].num_levels)
            b = g;
    }
    hdfe_coef_rows *D = &M->D;
    D->N = N;
    D->nB = S->factors[b].num_levels;
    D->lvB = S->factors[b].levels;
    D->w = S->weights;
    int nk = 0;
    ST_int m = 0;
    for (int g = 0; g < G; g++)
    {
        if (g == b)
            continue;
        D->lvk[nk] = S->factors[g].levels;
        D->off[nk] = m;
        M->nlev[nk] = S->factors[g].num_levels;
        if (M->nlev[nk] > 2147483647 - m)
            return 920;
        m += M->nlev[nk];
        nk++;
    }
    D->off[nk] = m;
    D->nk = nk;
    D->m = m;
    M->nk = nk;
    M->m = m;
    const ST_int nB = D->nB;
    const ST_double *w = D->w;
    const int par = N > HDFE_COEF_PARALLEL_MIN;

    ST_int *dest = NULL;
    int sorted = 0;
    ST_retcode rc = hdfe_coef_group_rows(N, D->lvB, nB, NULL, 0, &M->rowptr, &dest, &sorted, NULL);
    if (rc)
        return rc;
    D->sorted = sorted;
    int bad = 0;
#pragma omp parallel for schedule(static) reduction(| : bad) if (par)
    for (ST_int i = 0; i < N; i++)
    {
        for (int h = 0; h < nk; h++)
        {
            const ST_int l = D->lvk[h][i];
            if (l < 1 || l > M->nlev[h])
                bad |= 1;
        }
    }
    if (bad)
    {
        free(dest);
        return 198;
    }

    /* the rows in the order of B */
    if (sorted)
    {
        for (int h = 0; h < nk; h++)
            M->gk[h] = D->lvk[h];
        M->gw = w;
    }
    else
    {
        int failed = 0;
        for (int h = 0; h < nk; h++)
        {
            M->gk_own[h] = (ST_int *)malloc((size_t)N * sizeof(ST_int));
            if (!M->gk_own[h])
                failed = 1;
            M->gk[h] = M->gk_own[h];
        }
        if (w)
        {
            M->gw_own = (ST_double *)malloc((size_t)N * sizeof(ST_double));
            if (!M->gw_own)
                failed = 1;
        }
        if (failed)
        {
            free(dest);
            return 920;
        }
        M->gw = M->gw_own;
#pragma omp parallel for schedule(static) if (par)
        for (ST_int i = 0; i < N; i++)
        {
            const ST_int at = dest[i];
            for (int h = 0; h < nk; h++)
                M->gk_own[h][at] = D->lvk[h][i];
            if (w)
                M->gw_own[at] = w[i];
        }
        free(dest);
        dest = NULL;
    }

    M->invB = (ST_double *)malloc((size_t)nB * sizeof(ST_double));
    M->invk = (ST_double *)calloc((size_t)m, sizeof(ST_double));
    M->act = (ST_int *)malloc((size_t)nB * sizeof(ST_int));
    if (!M->invB || !M->invk || !M->act)
        return 920;
    D->invB = M->invB;
    const ST_int *rowptr = M->rowptr;
    const ST_double *gw = M->gw;
    ST_double *invB = M->invB, *invk = M->invk;
#pragma omp parallel for schedule(static) if (par)
    for (ST_int f = 0; f < nB; f++)
    {
        ST_double total = 0.0;
        if (gw)
        {
            for (ST_int e = rowptr[f]; e < rowptr[f + 1]; e++)
                total += gw[e];
        }
        else
            total = (ST_double)(rowptr[f + 1] - rowptr[f]);
        invB[f] = total > 0.0 ? 1.0 / total : 0.0;
    }
    /* weighted counts of the kept levels: one effect per thread, in row order */
#pragma omp parallel for schedule(dynamic, 1) if (par)
    for (int h = 0; h < nk; h++)
    {
        ST_double *count = invk + D->off[h];
        const ST_int *lv = D->lvk[h];
        for (ST_int i = 0; i < N; i++)
            count[lv[i] - 1] += w ? w[i] : 1.0;
        for (ST_int j = 0; j < M->nlev[h]; j++)
            count[j] = count[j] > 0.0 ? 1.0 / count[j] : 0.0;
    }
    ST_int nact = 0, nnz_act = 0;
    for (ST_int f = 0; f < nB; f++)
    {
        const ST_int rows = rowptr[f + 1] - rowptr[f];
        if (rows >= 2)
        {
            M->act[nact++] = f;
            nnz_act += rows;
        }
    }
    M->nact = nact;
    M->nnz_act = nnz_act;
    int nch = nnz_act / HDFE_COEF_MULTI_CHUNK_ROWS;
    if (nch > HDFE_COEF_MULTI_CHUNKS)
        nch = HDFE_COEF_MULTI_CHUNKS;
    if (nch < 1)
        nch = 1;
    while (nch > 1 && (long long)m * nch * 2 > (long long)nnz_act * nk)
        nch /= 2;
    M->nch = nch;
    M->ch_act = hdfe_coef_chunks(rowptr, M->act, nact, nch);
    if (!M->ch_act)
        return 920;

    if (need_crosstabs)
    {
        for (int a = 0; a < nk; a++)
        {
            for (int c = a + 1; c < nk; c++)
            {
                rc = hdfe_crosstab_build(&M->X[a][c], N, D->lvk[a], M->nlev[a], D->lvk[c],
                                         M->nlev[c], w);
                if (!rc)
                    rc = hcm_transpose(&M->X[a][c], &M->X[c][a]);
                if (rc)
                    return rc;
            }
        }
    }
    return 0;
}

/* ========================================================================
 * Sweep and preconditioner
 * ======================================================================== */

/* buf (m x K, zeroed here) = S p restricted to the levels a0..a1-1 of act */
HDFE_COEF_KERNEL void hcm_chunk_scatter(const hcm_system *M, ST_int a0, ST_int a1,
                                        const ST_double *RESTRICT p, ST_double *RESTRICT buf,
                                        const int nk, const int K)
{
    const ST_int *RESTRICT rowptr = M->rowptr;
    const ST_double *RESTRICT gw = M->gw, *RESTRICT invB = M->invB;
    const ST_int *off = M->D.off;
    ST_double s[HDFE_COEF_KMAX], v[HDFE_COEF_KMAX];
    ST_double kept[HDFE_COEF_MULTI_ROWBUF * HDFE_COEF_KMAX];
    memset(buf, 0, (size_t)M->m * K * sizeof(ST_double));
    for (ST_int a = a0; a < a1; a++)
    {
        const ST_int f = M->act[a];
        const ST_int e0 = rowptr[f], e1 = rowptr[f + 1];
        const int small = e1 - e0 <= HDFE_COEF_MULTI_ROWBUF;
        for (int k = 0; k < K; k++)
            s[k] = 0.0;
        for (ST_int e = e0; e < e1; e++)
        {
            const ST_double we = gw ? gw[e] : 1.0;
            const ST_double *p0 = p + (size_t)(off[0] + M->gk[0][e] - 1) * K;
            for (int k = 0; k < K; k++)
                v[k] = p0[k];
            for (int h = 1; h < nk; h++)
            {
                const ST_double *ph = p + (size_t)(off[h] + M->gk[h][e] - 1) * K;
                for (int k = 0; k < K; k++)
                    v[k] += ph[k];
            }
            for (int k = 0; k < K; k++)
                s[k] += we * v[k];
            if (small)
            {
                for (int k = 0; k < K; k++)
                    kept[(size_t)(e - e0) * K + k] = v[k];
            }
        }
        const ST_double ib = invB[f];
        for (int k = 0; k < K; k++)
            s[k] *= ib;
        for (ST_int e = e0; e < e1; e++)
        {
            const ST_double we = gw ? gw[e] : 1.0;
            if (small)
            {
                for (int k = 0; k < K; k++)
                    v[k] = kept[(size_t)(e - e0) * K + k];
            }
            else
            {
                const ST_double *p0 = p + (size_t)(off[0] + M->gk[0][e] - 1) * K;
                for (int k = 0; k < K; k++)
                    v[k] = p0[k];
                for (int h = 1; h < nk; h++)
                {
                    const ST_double *ph = p + (size_t)(off[h] + M->gk[h][e] - 1) * K;
                    for (int k = 0; k < K; k++)
                        v[k] += ph[k];
                }
            }
            for (int k = 0; k < K; k++)
                v[k] = we * (v[k] - s[k]);
            for (int h = 0; h < nk; h++)
            {
                ST_double *bj = buf + (size_t)(off[h] + M->gk[h][e] - 1) * K;
                for (int k = 0; k < K; k++)
                    bj[k] += v[k];
            }
        }
    }
}

/* q = S p (both m x K) */
static void hcm_apply(void *ctx, const ST_double *p, ST_double *q)
{
    const hcm_sweep *W = (const hcm_sweep *)ctx;
    const hcm_system *M = W->M;
    const int K = W->K, nk = M->nk, nch = M->nch;
    const size_t stride = (size_t)M->m * (size_t)K;
    const int par = (size_t)M->nnz_act * (size_t)K * (size_t)nk > HDFE_COEF_PARALLEL_MIN;
#pragma omp parallel for schedule(dynamic, 1) if (par && nch > 1)
    for (int c = 0; c < nch; c++)
    {
        HDFE_COEF_K_DISPATCH(K, hcm_chunk_scatter(M, M->ch_act[c], M->ch_act[c + 1], p,
                                                  W->buf + (size_t)c * stride, nk, KK))
    }
#pragma omp parallel for schedule(static) if (par && stride > 4096)
    for (ST_int j = 0; j < M->m; j++)
    {
        for (int k = 0; k < K; k++)
        {
            const size_t ix = (size_t)j * K + k;
            ST_double acc = 0.0;
            for (int c = 0; c < nch; c++)
                acc += W->buf[(size_t)c * stride + ix];
            q[ix] = acc;
        }
    }
}

/* One block of the Gauss-Seidel sweep: forward, z_a = N_a^-1 (r_a - sum
 * over earlier effects of C_ab z_b); backward, z_a -= N_a^-1 (sum over
 * later effects of C_ab z_b). Each level is computed by one thread. */
HDFE_COEF_KERNEL void hcm_sgs_block(const hcm_system *M, int a, int backward,
                                    const ST_double *RESTRICT r, ST_double *z, const int K)
{
    const ST_int na = M->nlev[a], off = M->D.off[a];
    const int nk = M->nk;
    const int b0 = backward ? a + 1 : 0, b1 = backward ? nk : a;
#pragma omp parallel for schedule(static) if ((size_t)na * (size_t)K > 20000)
    for (ST_int j = 0; j < na; j++)
    {
        ST_double acc[HDFE_COEF_KMAX];
        for (int k = 0; k < K; k++)
            acc[k] = 0.0;
        for (int b = b0; b < b1; b++)
        {
            const hdfe_crosstab *X = &M->X[a][b];
            const ST_int boff = M->D.off[b];
            for (ST_int e = X->rowptr[j]; e < X->rowptr[j + 1]; e++)
            {
                const ST_double *zc = z + (size_t)(boff + X->col[e] - 1) * K;
                const ST_double v = X->val ? X->val[e] : 1.0;
                for (int k = 0; k < K; k++)
                    acc[k] += v * zc[k];
            }
        }
        const ST_double inv = M->invk[off + j];
        ST_double *zj = z + (size_t)(off + j) * K;
        if (backward)
        {
            for (int k = 0; k < K; k++)
                zj[k] -= inv * acc[k];
        }
        else
        {
            const ST_double *rj = r + (size_t)(off + j) * K;
            for (int k = 0; k < K; k++)
                zj[k] = (rj[k] - acc[k]) * inv;
        }
    }
}

/* z = (N + L)^-T N (N + L)^-1 r over the kept effects */
static void hcm_sgs(void *ctx, const ST_double *r, ST_double *z)
{
    const hcm_sweep *W = (const hcm_sweep *)ctx;
    const hcm_system *M = W->M;
    const int K = W->K;
    for (int a = 0; a < M->nk; a++)
    {
        HDFE_COEF_K_DISPATCH(K, hcm_sgs_block(M, a, 0, r, z, KK))
    }
    for (int a = M->nk - 2; a >= 0; a--)
    {
        HDFE_COEF_K_DISPATCH(K, hcm_sgs_block(M, a, 1, r, z, KK))
    }
}

/* a[f] -= (sum over the rows of level f of w * the row's coefficients) / n_f */
HDFE_COEF_KERNEL void hcm_rows_eliminated(const hcm_system *M, ST_int f0, ST_int f1,
                                          const ST_double *RESTRICT beta, ST_double *RESTRICT a,
                                          const int nk, const int K)
{
    const ST_int *RESTRICT rowptr = M->rowptr;
    const ST_double *RESTRICT gw = M->gw;
    const ST_int *off = M->D.off;
    ST_double s[HDFE_COEF_KMAX];
    for (ST_int f = f0; f < f1; f++)
    {
        for (int k = 0; k < K; k++)
            s[k] = 0.0;
        for (ST_int e = rowptr[f]; e < rowptr[f + 1]; e++)
        {
            const ST_double we = gw ? gw[e] : 1.0;
            for (int h = 0; h < nk; h++)
            {
                const ST_double *bj = beta + (size_t)(off[h] + M->gk[h][e] - 1) * K;
                for (int k = 0; k < K; k++)
                    s[k] += we * bj[k];
            }
        }
        ST_double *af = a + (size_t)f * K;
        const ST_double ib = M->invB[f];
        for (int k = 0; k < K; k++)
            af[k] -= s[k] * ib;
    }
}

static void hcm_eliminated(const hcm_system *M, int K, ST_double *a0, const ST_double *beta)
{
    const ST_int nB = M->D.nB;
    const int nk = M->nk;
    const ST_int nbB = (nB + 4095) / 4096;
#pragma omp parallel for schedule(dynamic,                                                         \
                                      4) if ((size_t)M->D.N * (size_t)K > HDFE_COEF_PARALLEL_MIN)
    for (ST_int blk = 0; blk < nbB; blk++)
    {
        const ST_int f0 = blk * 4096, f1 = nB - f0 > 4096 ? f0 + 4096 : nB;
        HDFE_COEF_K_DISPATCH(K, hcm_rows_eliminated(M, f0, f1, beta, a0, nk, KK))
    }
}

/* ========================================================================
 * Dense S
 * ======================================================================== */

/* Lower triangle of S (row-major m x m) */
static ST_retcode hcm_dense(const hcm_system *M, ST_double *Sd)
{
    const ST_int m = M->m;
    const int nk = M->nk;
    const size_t m2 = (size_t)m * (size_t)m;
    const ST_int *rowptr = M->rowptr, *off = M->D.off;
    const ST_double *gw = M->gw;
    int nch = M->nch;
    while (nch > 1 && m2 * (size_t)nch * sizeof(ST_double) > 2 * HDFE_COEF_BUFFER_BYTES)
        nch /= 2;
    ST_int *bound = hdfe_coef_chunks(rowptr, M->act, M->nact, nch);
    ST_double *buf = (ST_double *)malloc(m2 * (size_t)nch * sizeof(ST_double));
    ST_double *sums = (ST_double *)calloc((size_t)m * (size_t)nch, sizeof(ST_double));
    ST_int *lists = (ST_int *)malloc((size_t)m * (size_t)nch * sizeof(ST_int));
    unsigned char *marks = (unsigned char *)calloc((size_t)m * (size_t)nch, 1);
    if (!bound || !buf || !sums || !lists || !marks)
    {
        free(bound);
        free(buf);
        free(sums);
        free(lists);
        free(marks);
        return 920;
    }
#pragma omp parallel for schedule(dynamic, 1) if (nch > 1)
    for (int c = 0; c < nch; c++)
    {
        ST_double *RESTRICT B = buf + (size_t)c * m2;
        ST_double *RESTRICT u = sums + (size_t)c * m;
        ST_int *RESTRICT list = lists + (size_t)c * m;
        unsigned char *RESTRICT seen = marks + (size_t)c * m;
        memset(B, 0, m2 * sizeof(ST_double));
        for (ST_int a = bound[c]; a < bound[c + 1]; a++)
        {
            const ST_int f = M->act[a];
            const ST_double ib = M->invB[f];
            ST_int distinct = 0;
            for (ST_int e = rowptr[f]; e < rowptr[f + 1]; e++)
            {
                const ST_double we = gw ? gw[e] : 1.0;
                ST_int at[HDFE_COEF_MAX_EFFECTS];
                for (int h = 0; h < nk; h++)
                {
                    const ST_int j = off[h] + M->gk[h][e] - 1;
                    at[h] = j;
                    if (!seen[j])
                    {
                        seen[j] = 1;
                        list[distinct++] = j;
                    }
                    u[j] += we;
                }
                /* Dk' W Dk */
                for (int h = 0; h < nk; h++)
                {
                    for (int g = 0; g <= h; g++)
                    {
                        const ST_int hi = at[h] > at[g] ? at[h] : at[g];
                        const ST_int lo = at[h] > at[g] ? at[g] : at[h];
                        B[(size_t)hi * m + lo] += we;
                    }
                }
            }
            /* less the level's rank-one term */
            if (2 * distinct > m)
            {
                for (ST_int j1 = 0; j1 < m; j1++)
                {
                    const ST_double v1 = u[j1] * ib;
                    if (v1 == 0.0)
                        continue;
                    ST_double *RESTRICT Bj = B + (size_t)j1 * m;
                    for (ST_int j2 = 0; j2 <= j1; j2++)
                        Bj[j2] -= v1 * u[j2];
                }
            }
            else
            {
                for (ST_int x = 0; x < distinct; x++)
                {
                    const ST_int j1 = list[x];
                    const ST_double v1 = u[j1] * ib;
                    ST_double *RESTRICT Bj = B + (size_t)j1 * m;
                    for (ST_int y = 0; y < distinct; y++)
                    {
                        const ST_int j2 = list[y];
                        if (j2 <= j1)
                            Bj[j2] -= v1 * u[j2];
                    }
                }
            }
            for (ST_int x = 0; x < distinct; x++)
            {
                u[list[x]] = 0.0;
                seen[list[x]] = 0;
            }
        }
    }
#pragma omp parallel for schedule(static) if (m2 > HDFE_COEF_PARALLEL_MIN)
    for (ST_int j1 = 0; j1 < m; j1++)
    {
        for (ST_int j2 = 0; j2 <= j1; j2++)
        {
            ST_double total = 0.0;
            for (int c = 0; c < nch; c++)
                total += buf[(size_t)c * m2 + (size_t)j1 * m + j2];
            Sd[(size_t)j1 * m + j2] = total;
        }
    }
    free(bound);
    free(buf);
    free(sums);
    free(lists);
    free(marks);
    return 0;
}

/* ========================================================================
 * Driver
 * ======================================================================== */

static hdfe_coef_method hcm_choose(const hcm_system *M, ST_int K, hdfe_coef_method method)
{
    if (method == HDFE_COEF_DENSE)
        return HDFE_COEF_DENSE;
    if (method != HDFE_COEF_AUTO || M->m > HDFE_COEF_MULTI_DENSE_MAX)
        return HDFE_COEF_PCG;
    const ST_double m = (ST_double)M->m;
    ST_double work = m * m * m / 3.0;
    for (ST_int a = 0; a < M->nact; a++)
    {
        ST_double cells = (ST_double)(M->rowptr[M->act[a] + 1] - M->rowptr[M->act[a]]) * M->nk;
        if (cells > m)
            cells = m;
        work += cells * cells;
    }
    return work <= HDFE_COEF_MULTI_DENSE_SWEEPS * (ST_double)M->nnz_act * M->nk * (ST_double)K + 2e6
               ? HDFE_COEF_DENSE
               : HDFE_COEF_PCG;
}

HDFE_SolveResult hdfe_coef_multi_partial_out(const HDFE_State *S, ST_double *data, ST_int N,
                                             ST_int K, hdfe_coef_method method)
{
    HDFE_SolveResult result = {0, 1, 0};
    hcm_system M;
    hdfe_coef_rowplan R;
    hcm_sweep W = {&M, 0, NULL};
    ST_double *a0 = NULL, *g = NULL, *beta = NULL, *pcg = NULL, *Sd = NULL;
    ST_double weakest = 1.0;
    unsigned char *dead = NULL, *pin = NULL;
    memset(&R, 0, sizeof(R));

    /* the cross-tabulations serve the iteration only; the choice needs the
     * system, so they are built for every method but the forced dense one */
    ST_retcode rc = hcm_system_build(&M, S, method != HDFE_COEF_DENSE);
    if (rc)
        goto done;
    const ST_int nB = M.D.nB, m = M.m;
    if (method == HDFE_COEF_DENSE && m > 46340)
    {
        rc = 198;
        goto done;
    }
    const hdfe_coef_method chosen = hcm_choose(&M, K, method);
    if (ctools_cancel_requested())
    {
        rc = 1;
        goto done;
    }

    const int KB = K < HDFE_COEF_KMAX ? (int)K : HDFE_COEF_KMAX;
    a0 = (ST_double *)malloc((size_t)nB * (size_t)KB * sizeof(ST_double));
    g = (ST_double *)malloc((size_t)m * (size_t)KB * sizeof(ST_double));
    beta = (ST_double *)malloc((size_t)m * (size_t)KB * sizeof(ST_double));
    if (!a0 || !g || !beta || (rc = hdfe_coef_rowplan_init(&R, &M.D, KB)))
    {
        if (!rc)
            rc = 920;
        goto done;
    }
    if (chosen == HDFE_COEF_DENSE)
    {
        Sd = (ST_double *)malloc((size_t)m * (size_t)m * sizeof(ST_double));
        dead = (unsigned char *)malloc((size_t)m);
        pin = (unsigned char *)calloc((size_t)m, 1);
        if (!Sd || !dead || !pin)
        {
            rc = 920;
            goto done;
        }
        if ((rc = hcm_dense(&M, Sd)))
            goto done;
        /* A constant added to every level of one kept effect is absorbed by
         * B, so S is singular at least once per kept effect whatever the
         * data. Holding one level of each at zero by name leaves the pivot
         * threshold to find only what the data add (nested or disconnected
         * effects). The last level with a nonzero diagonal is the one the
         * elimination would reach with nothing left. */
        for (int h = 0; h < M.nk; h++)
        {
            for (ST_int j = M.D.off[h + 1] - 1; j >= M.D.off[h]; j--)
            {
                if (Sd[(size_t)j * m + j] > 0.0)
                {
                    pin[j] = 1;
                    break;
                }
            }
        }
        if ((rc = hdfe_coef_chol_factor(Sd, m, pin, dead, &weakest)))
            goto done;
    }
    else
    {
        W.buf = (ST_double *)malloc((size_t)M.nch * (size_t)m * (size_t)KB * sizeof(ST_double));
        pcg = (ST_double *)malloc(HDFE_COEF_PCG_WORK(m, KB) * sizeof(ST_double));
        if (!W.buf || !pcg)
        {
            rc = 920;
            goto done;
        }
    }
    if (ctools_cancel_requested())
    {
        rc = 1;
        goto done;
    }

    for (ST_int k0 = 0; k0 < K; k0 += KB)
    {
        /* The columns first..first+kb-1 of the batch; after the first round,
         * those that stopped at their rounding (hdfe_coef_internal.h) */
        int first = 0, kb = K - k0 < KB ? (int)(K - k0) : KB;
        ST_int total = 0;
        for (int round = 0; round < HDFE_COEF_ROUNDS; round++)
        {
            ST_double *X = data + (size_t)(k0 + first) * (size_t)N;
            ST_double f0[HDFE_COEF_KMAX], scale[HDFE_COEF_KMAX];
            unsigned char state[HDFE_COEF_KMAX] = {0};
            ST_int iterations = 1;
            W.K = kb;
            hdfe_coef_rows_rhs(&M.D, &R, X, kb, a0, g, f0, scale);
            if (ctools_cancel_requested())
            {
                rc = 1;
                goto done;
            }
            if (chosen == HDFE_COEF_DENSE)
            {
                memcpy(beta, g, (size_t)m * (size_t)kb * sizeof(ST_double));
                hdfe_coef_chol_solve(Sd, m, dead, beta, kb);
                hdfe_coef_direct_state(m, kb, g, beta, f0, scale, S->tolerance, weakest, state);
            }
            else
            {
                const hdfe_coef_pcg_problem P = {m,    kb,           hcm_apply,  &W, hcm_sgs, &W,
                                                 NULL, S->tolerance, S->maxiter, 0,  0};
                const ST_retcode solve =
                    hdfe_coef_pcg(&P, g, beta, f0, scale, pcg, state, &iterations);
                if (solve == 1)
                {
                    rc = 1;
                    goto done;
                }
                if (solve == 430)
                    result.converged = 0;
            }
            /* a round counts as one iteration even when nothing was left
             * for the kept effects to explain */
            total += iterations > 1 ? iterations : 1;
            if (ctools_cancel_requested())
            {
                rc = 1;
                goto done;
            }
            hcm_eliminated(&M, kb, a0, beta);
            hdfe_coef_rows_residual(&M.D, X, kb, a0, beta);
            if (!hdfe_coef_refine_span(state, kb, &first, &kb))
                break;
        }
        if (total > result.iterations)
            result.iterations = total;
    }

done:
    hdfe_coef_rowplan_free(&R);
    free(W.buf);
    free(a0);
    free(g);
    free(beta);
    free(pcg);
    free(Sd);
    free(dead);
    free(pin);
    hcm_system_free(&M);
    if (rc)
    {
        result.status = rc;
        result.converged = 0;
    }
    else if (!result.converged)
        result.status = 430;
    return result;
}

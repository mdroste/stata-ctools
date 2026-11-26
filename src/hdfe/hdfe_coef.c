/*
 * hdfe_coef.c
 *
 * Coefficient-space projection on two fixed effects, and the entry point
 * of the module (hdfe_coef.h)
 * Part of the ctools Stata plugin suite
 *
 * B is the effect with more levels and is eliminated; S is kept. With C the
 * cross-tabulation (rows B, columns S), N_B and N_S the weighted level
 * counts,
 *
 *     S = N_S - C' N_B^-1 C
 *
 * is applied through the cells. A row of C with a single cell (a worker
 * who never changes firm) links nothing: it only lowers the diagonal, so
 * the sweeps skip it.
 *
 * For each batch of columns:
 *   1. a0 = level means of B                       (one pass over the rows)
 *   2. g = D_S' W (x - D_B a0), f0 = |x - D_B a0|^2 (one pass over the rows)
 *   3. S b = g                                      (cells only)
 *   4. a = a0 - N_B^-1 C b, x <- x - a - b          (one pass over the rows)
 * Step 2 works on the B-demeaned column so that level means much larger
 * than the within variation cost no digits.
 */

#include "hdfe_coef_internal.h"

/* Largest kept effect solved by dense Cholesky: the factorization is n^3/3
 * operations on one core's worth of dot products (0.03 s at 1,000). */
#define HDFE_COEF_DENSE_MAX 1000

/* The dense path is taken when assembling and factoring S costs no more
 * than about a dozen sweeps of the iteration: one unit of dense work costs
 * 0.6 of one cell-column of a sweep (validation/performance/
 * benchmark_hdfe_coef.c, 10M rows: 2.5e8 units in 8 ms against 4 sweeps of
 * 4.3M cells and 3 columns in 4 ms). */
#define HDFE_COEF_DENSE_SWEEPS 20.0

/* Chunks of rows of B that scatter into private buffers in the fused sweep */
#define HDFE_COEF_SWEEP_CHUNKS 32
/* Cells per chunk below which more chunks only add overhead */
#define HDFE_COEF_CHUNK_CELLS 16384

typedef struct
{
    ST_int N, nB, nS;
    const ST_int *lvB, *lvS; /* 1-based levels of every row */
    const ST_double *w;      /* row weights, NULL when unweighted */
    hdfe_crosstab C;         /* rows B, columns S */
    ST_double *invB;         /* nB: 1 / weighted count (0 for an empty level) */
    ST_double *nSw;          /* nS: weighted counts */
    ST_double *dS;           /* nS: nSw less the single-cell rows (fused sweep) */
    ST_double *diag;         /* nS: diagonal of S */
    ST_int nact, nnz_act;    /* rows of B with two or more cells */
    ST_int *act;
    int twopass; /* gather-only sweep through the cells by column */
    int nch;     /* chunks of the fused sweep */
    ST_int *ch_act, *ch_all;
    ST_int *cptr, *crow; /* two-pass: cells grouped by level of S */
    ST_double *cval;
    hdfe_coef_rows D; /* the data rows, for the accumulating passes */
} hc2_system;

typedef struct
{
    const hc2_system *M;
    int K;          /* columns in the batch */
    ST_double *buf; /* fused sweep: nch * nS * K */
    ST_double *tB;  /* two-pass sweep: nB * K */
} hc2_sweep;

static void hc2_system_free(hc2_system *M)
{
    hdfe_crosstab_free(&M->C);
    free(M->invB);
    free(M->nSw);
    free(M->dS);
    free(M->diag);
    free(M->act);
    free(M->ch_act);
    free(M->ch_all);
    free(M->cptr);
    free(M->crow);
    free(M->cval);
}

/* ========================================================================
 * Sweeps over the cells
 * ======================================================================== */

/* buf (nS x K, zeroed here) = sum over the rows a0..a1-1 of rows[] of
 * c_f (c_f' p) / n_f, with c_f the row's cells */
HDFE_COEF_KERNEL void hc2_chunk_scatter(const hc2_system *M, const ST_int *rows, ST_int a0,
                                        ST_int a1, const ST_double *RESTRICT p,
                                        ST_double *RESTRICT buf, const int K)
{
    const ST_int *RESTRICT rowptr = M->C.rowptr, *RESTRICT col = M->C.col;
    const ST_double *RESTRICT val = M->C.val, *RESTRICT invB = M->invB;
    ST_double s[HDFE_COEF_KMAX];
    memset(buf, 0, (size_t)M->nS * K * sizeof(ST_double));
    for (ST_int a = a0; a < a1; a++)
    {
        const ST_int f = rows[a];
        const ST_int e0 = rowptr[f], e1 = rowptr[f + 1];
        for (int k = 0; k < K; k++)
            s[k] = 0.0;
        if (val)
        {
            for (ST_int e = e0; e < e1; e++)
            {
                const ST_double *pj = p + (size_t)(col[e] - 1) * K;
                const ST_double v = val[e];
                for (int k = 0; k < K; k++)
                    s[k] += v * pj[k];
            }
        }
        else
        {
            for (ST_int e = e0; e < e1; e++)
            {
                const ST_double *pj = p + (size_t)(col[e] - 1) * K;
                for (int k = 0; k < K; k++)
                    s[k] += pj[k];
            }
        }
        const ST_double ib = invB[f];
        for (int k = 0; k < K; k++)
            s[k] *= ib;
        if (val)
        {
            for (ST_int e = e0; e < e1; e++)
            {
                ST_double *bj = buf + (size_t)(col[e] - 1) * K;
                const ST_double v = val[e];
                for (int k = 0; k < K; k++)
                    bj[k] += v * s[k];
            }
        }
        else
        {
            for (ST_int e = e0; e < e1; e++)
            {
                ST_double *bj = buf + (size_t)(col[e] - 1) * K;
                for (int k = 0; k < K; k++)
                    bj[k] += s[k];
            }
        }
    }
}

/* t[f] = (c_f' p) / n_f for the rows f0..f1-1 of B; with subtract,
 * t[f] -= (c_f' p) / n_f */
HDFE_COEF_KERNEL void hc2_rows_gather(const hc2_system *M, ST_int f0, ST_int f1,
                                      const ST_double *RESTRICT p, ST_double *RESTRICT t,
                                      const int subtract, const int K)
{
    const ST_int *RESTRICT rowptr = M->C.rowptr, *RESTRICT col = M->C.col;
    const ST_double *RESTRICT val = M->C.val, *RESTRICT invB = M->invB;
    ST_double s[HDFE_COEF_KMAX];
    for (ST_int f = f0; f < f1; f++)
    {
        for (int k = 0; k < K; k++)
            s[k] = 0.0;
        if (val)
        {
            for (ST_int e = rowptr[f]; e < rowptr[f + 1]; e++)
            {
                const ST_double *pj = p + (size_t)(col[e] - 1) * K;
                const ST_double v = val[e];
                for (int k = 0; k < K; k++)
                    s[k] += v * pj[k];
            }
        }
        else
        {
            for (ST_int e = rowptr[f]; e < rowptr[f + 1]; e++)
            {
                const ST_double *pj = p + (size_t)(col[e] - 1) * K;
                for (int k = 0; k < K; k++)
                    s[k] += pj[k];
            }
        }
        ST_double *tf = t + (size_t)f * K;
        const ST_double ib = invB[f];
        if (subtract)
        {
            for (int k = 0; k < K; k++)
                tf[k] -= s[k] * ib;
        }
        else
        {
            for (int k = 0; k < K; k++)
                tf[k] = s[k] * ib;
        }
    }
}

/* q[j] = n_j p[j] - sum over the cells of level j of c t[row] */
HDFE_COEF_KERNEL void hc2_cols_gather(const hc2_system *M, ST_int j0, ST_int j1,
                                      const ST_double *RESTRICT p, const ST_double *RESTRICT t,
                                      ST_double *RESTRICT q, const int K)
{
    const ST_int *RESTRICT cptr = M->cptr, *RESTRICT crow = M->crow;
    const ST_double *RESTRICT cval = M->cval;
    ST_double s[HDFE_COEF_KMAX];
    for (ST_int j = j0; j < j1; j++)
    {
        for (int k = 0; k < K; k++)
            s[k] = 0.0;
        if (cval)
        {
            for (ST_int e = cptr[j]; e < cptr[j + 1]; e++)
            {
                const ST_double *tf = t + (size_t)crow[e] * K;
                const ST_double v = cval[e];
                for (int k = 0; k < K; k++)
                    s[k] += v * tf[k];
            }
        }
        else
        {
            for (ST_int e = cptr[j]; e < cptr[j + 1]; e++)
            {
                const ST_double *tf = t + (size_t)crow[e] * K;
                for (int k = 0; k < K; k++)
                    s[k] += tf[k];
            }
        }
        const ST_double nj = M->nSw[j];
        const ST_double *pj = p + (size_t)j * K;
        ST_double *qj = q + (size_t)j * K;
        for (int k = 0; k < K; k++)
            qj[k] = nj * pj[k] - s[k];
    }
}

/* q = S p (both nS x K) */
static void hc2_apply(void *ctx, const ST_double *p, ST_double *q)
{
    const hc2_sweep *W = (const hc2_sweep *)ctx;
    const hc2_system *M = W->M;
    const int K = W->K;
    const ST_int nS = M->nS, nB = M->nB;
    const int par = (size_t)M->nnz_act * (size_t)K > HDFE_COEF_PARALLEL_MIN;
    if (!M->twopass)
    {
        const int nch = M->nch;
        const size_t stride = (size_t)nS * (size_t)K;
#pragma omp parallel for schedule(dynamic, 1) if (par && nch > 1)
        for (int c = 0; c < nch; c++)
        {
            HDFE_COEF_K_DISPATCH(K, hc2_chunk_scatter(M, M->act, M->ch_act[c], M->ch_act[c + 1], p,
                                                      W->buf + (size_t)c * stride, KK))
        }
        /* chunks added in order */
#pragma omp parallel for schedule(static) if (par && stride > 4096)
        for (ST_int j = 0; j < nS; j++)
        {
            const ST_double dj = M->dS[j];
            for (int k = 0; k < K; k++)
            {
                const size_t ix = (size_t)j * K + k;
                ST_double acc = 0.0;
                for (int c = 0; c < nch; c++)
                    acc += W->buf[(size_t)c * stride + ix];
                q[ix] = dj * p[ix] - acc;
            }
        }
    }
    else
    {
        const ST_int nbB = (nB + 4095) / 4096, nbS = (nS + 1023) / 1024;
#pragma omp parallel for schedule(dynamic, 4) if (par)
        for (ST_int blk = 0; blk < nbB; blk++)
        {
            const ST_int f0 = blk * 4096, f1 = nB - f0 > 4096 ? f0 + 4096 : nB;
            HDFE_COEF_K_DISPATCH(K, hc2_rows_gather(M, f0, f1, p, W->tB, 0, KK))
        }
#pragma omp parallel for schedule(dynamic, 4) if (par)
        for (ST_int blk = 0; blk < nbS; blk++)
        {
            const ST_int j0 = blk * 1024, j1 = nS - j0 > 1024 ? j0 + 1024 : nS;
            HDFE_COEF_K_DISPATCH(K, hc2_cols_gather(M, j0, j1, p, W->tB, q, KK))
        }
    }
}

/* ========================================================================
 * Building the system
 * ======================================================================== */

/* Cells grouped by level of S, for the two-pass sweep */
static ST_retcode hc2_by_column(hc2_system *M)
{
    const hdfe_crosstab *C = &M->C;
    const ST_int nnz = C->nnz, nB = M->nB;
    ST_int *dest = NULL;
    int sorted = 0;
    ST_retcode rc =
        hdfe_coef_group_rows(nnz, C->col, M->nS, NULL, 0, &M->cptr, &dest, &sorted, NULL);
    if (rc)
        return rc;
    M->crow = (ST_int *)malloc((size_t)nnz * sizeof(ST_int));
    M->cval = C->val ? (ST_double *)malloc((size_t)nnz * sizeof(ST_double)) : NULL;
    if (!M->crow || (C->val && !M->cval))
    {
        free(dest);
        return 920;
    }
    ST_int *crow = M->crow;
    ST_double *cval = M->cval;
#pragma omp parallel for schedule(dynamic, 1024) if (nnz > HDFE_COEF_PARALLEL_MIN)
    for (ST_int f = 0; f < nB; f++)
    {
        for (ST_int e = C->rowptr[f]; e < C->rowptr[f + 1]; e++)
        {
            const ST_int at = dest ? dest[e] : e;
            crow[at] = f;
            if (cval)
                cval[at] = C->val[e];
        }
    }
    free(dest);
    return 0;
}

static ST_retcode hc2_system_build(hc2_system *M, const HDFE_State *S, int allow_twopass)
{
    memset(M, 0, sizeof(*M));
    /* eliminate the effect with more levels */
    const int b = S->factors[1].num_levels > S->factors[0].num_levels ? 1 : 0;
    const FE_Factor *fB = &S->factors[b], *fS = &S->factors[1 - b];
    const ST_int N = S->N, nB = fB->num_levels, nS = fS->num_levels;
    M->N = N;
    M->nB = nB;
    M->nS = nS;
    M->lvB = fB->levels;
    M->lvS = fS->levels;
    M->w = S->weights;

    ST_retcode rc = hdfe_crosstab_build(&M->C, N, M->lvB, nB, M->lvS, nS, M->w);
    if (rc)
        return rc;
    const hdfe_crosstab *C = &M->C;
    const ST_int *rowptr = C->rowptr, *col = C->col;
    const ST_double *val = C->val;
    const int par = C->nnz > HDFE_COEF_PARALLEL_MIN;

    M->invB = (ST_double *)malloc((size_t)nB * sizeof(ST_double));
    M->nSw = (ST_double *)malloc((size_t)nS * sizeof(ST_double));
    M->dS = (ST_double *)malloc((size_t)nS * sizeof(ST_double));
    M->diag = (ST_double *)malloc((size_t)nS * sizeof(ST_double));
    M->act = (ST_int *)malloc((size_t)nB * sizeof(ST_int));
    if (!M->invB || !M->nSw || !M->dS || !M->diag || !M->act)
        return 920;

    /* weighted counts of B, each over its cells in order */
    ST_double *invB = M->invB;
#pragma omp parallel for schedule(static) if (par)
    for (ST_int f = 0; f < nB; f++)
    {
        ST_double total = 0.0;
        if (val)
        {
            for (ST_int e = rowptr[f]; e < rowptr[f + 1]; e++)
                total += val[e];
        }
        else
            total = (ST_double)(rowptr[f + 1] - rowptr[f]);
        invB[f] = total > 0.0 ? 1.0 / total : 0.0;
    }
    ST_int nact = 0, nnz_act = 0;
    for (ST_int f = 0; f < nB; f++)
    {
        const ST_int cells = rowptr[f + 1] - rowptr[f];
        if (cells >= 2)
        {
            M->act[nact++] = f;
            nnz_act += cells;
        }
    }
    M->nact = nact;
    M->nnz_act = nnz_act;

    /* The fused sweep scatters into one private buffer per chunk; it is
     * used while those buffers stay small next to the cells they replace.
     * Otherwise the sweep gathers twice, by row and then by column. */
    int nch = nnz_act / HDFE_COEF_CHUNK_CELLS;
    if (nch > HDFE_COEF_SWEEP_CHUNKS)
        nch = HDFE_COEF_SWEEP_CHUNKS;
    if (nch < 1)
        nch = 1;
    const int wanted = nch;
    while (nch > 1 && (long long)nS * nch * 4 > (long long)nnz_act)
        nch /= 2;
    M->twopass = allow_twopass && nch < wanted && nch < 8 && nnz_act >= 200000;
    M->nch = nch;
    M->ch_act = hdfe_coef_chunks(rowptr, M->act, nact, nch);
    M->ch_all = hdfe_coef_chunks(rowptr, NULL, nB, nch);
    if (!M->ch_act || !M->ch_all)
        return 920;

    /* Weighted counts of S, the part of them in single-cell rows and the
     * diagonal of S: per chunk, then chunks in order. */
    {
        const size_t stride = (size_t)3 * (size_t)nS;
        ST_double *acc = (ST_double *)malloc((size_t)nch * stride * sizeof(ST_double));
        if (!acc)
            return 920;
#pragma omp parallel for schedule(dynamic, 1) if (par && nch > 1)
        for (int c = 0; c < nch; c++)
        {
            ST_double *RESTRICT count = acc + (size_t)c * stride;
            ST_double *RESTRICT single = count + nS, *RESTRICT square = count + 2 * (size_t)nS;
            memset(count, 0, stride * sizeof(ST_double));
            for (ST_int f = M->ch_all[c]; f < M->ch_all[c + 1]; f++)
            {
                const ST_int e0 = rowptr[f], e1 = rowptr[f + 1];
                if (e1 - e0 == 1)
                {
                    const ST_double v = val ? val[e0] : 1.0;
                    count[col[e0] - 1] += v;
                    single[col[e0] - 1] += v;
                    continue;
                }
                const ST_double ib = invB[f];
                for (ST_int e = e0; e < e1; e++)
                {
                    const ST_double v = val ? val[e] : 1.0;
                    count[col[e] - 1] += v;
                    square[col[e] - 1] += v * v * ib;
                }
            }
        }
        for (ST_int j = 0; j < nS; j++)
        {
            ST_double count = 0.0, single = 0.0, square = 0.0;
            for (int c = 0; c < nch; c++)
            {
                const ST_double *a = acc + (size_t)c * stride;
                count += a[j];
                single += a[(size_t)nS + j];
                square += a[2 * (size_t)nS + j];
            }
            M->nSw[j] = count;
            M->dS[j] = count - single;
            M->diag[j] = count - single - square;
        }
        free(acc);
    }
    M->D.N = N;
    M->D.nB = nB;
    M->D.m = nS;
    M->D.nk = 1;
    M->D.sorted = C->sorted;
    M->D.lvB = M->lvB;
    M->D.lvk[0] = M->lvS;
    M->D.off[0] = 0;
    M->D.off[1] = nS;
    M->D.w = M->w;
    M->D.invB = M->invB;
    if (M->twopass)
        return hc2_by_column(M);
    return 0;
}

/* ========================================================================
 * Dense S
 * ======================================================================== */

/* Lower triangle of sum over the active rows of c c' / n (row-major nS x nS) */
static ST_retcode hc2_dense_sum(const hc2_system *M, ST_double *Sd)
{
    const ST_int nS = M->nS;
    const size_t n2 = (size_t)nS * (size_t)nS;
    const ST_int *rowptr = M->C.rowptr, *col = M->C.col;
    const ST_double *val = M->C.val;
    int nch = M->nch;
    while (nch > 1 && n2 * (size_t)nch * sizeof(ST_double) > 2 * HDFE_COEF_BUFFER_BYTES)
        nch /= 2;
    ST_int *bound = hdfe_coef_chunks(rowptr, M->act, M->nact, nch);
    ST_double *buf = (ST_double *)malloc(n2 * (size_t)nch * sizeof(ST_double));
    ST_double *dense = (ST_double *)calloc((size_t)nS * (size_t)nch, sizeof(ST_double));
    if (!bound || !buf || !dense)
    {
        free(bound);
        free(buf);
        free(dense);
        return 920;
    }
#pragma omp parallel for schedule(dynamic, 1) if (nch > 1)
    for (int c = 0; c < nch; c++)
    {
        ST_double *RESTRICT B = buf + (size_t)c * n2;
        ST_double *RESTRICT row = dense + (size_t)c * nS;
        memset(B, 0, n2 * sizeof(ST_double));
        for (ST_int a = bound[c]; a < bound[c + 1]; a++)
        {
            const ST_int f = M->act[a];
            const ST_int e0 = rowptr[f], e1 = rowptr[f + 1], cells = e1 - e0;
            const ST_double ib = M->invB[f];
            if (2 * cells > nS)
            {
                /* a row that meets most levels: rank-one update on the row
                 * written out densely */
                for (ST_int e = e0; e < e1; e++)
                    row[col[e] - 1] = val ? val[e] : 1.0;
                for (ST_int j1 = 0; j1 < nS; j1++)
                {
                    const ST_double v1 = row[j1] * ib;
                    if (v1 == 0.0)
                        continue;
                    ST_double *RESTRICT Bj = B + (size_t)j1 * nS;
                    for (ST_int j2 = 0; j2 <= j1; j2++)
                        Bj[j2] += v1 * row[j2];
                }
                for (ST_int e = e0; e < e1; e++)
                    row[col[e] - 1] = 0.0;
            }
            else
            {
                for (ST_int x = e0; x < e1; x++)
                {
                    const ST_int j1 = col[x] - 1;
                    const ST_double v1 = (val ? val[x] : 1.0) * ib;
                    ST_double *RESTRICT Bj = B + (size_t)j1 * nS;
                    for (ST_int y = e0; y < e1; y++)
                    {
                        const ST_int j2 = col[y] - 1;
                        if (j2 <= j1)
                            Bj[j2] += v1 * (val ? val[y] : 1.0);
                    }
                }
            }
        }
    }
#pragma omp parallel for schedule(static) if (n2 > HDFE_COEF_PARALLEL_MIN)
    for (ST_int j1 = 0; j1 < nS; j1++)
    {
        for (ST_int j2 = 0; j2 <= j1; j2++)
        {
            ST_double total = 0.0;
            for (int c = 0; c < nch; c++)
                total += buf[(size_t)c * n2 + (size_t)j1 * nS + j2];
            Sd[(size_t)j1 * nS + j2] = total;
        }
    }
    free(bound);
    free(buf);
    free(dense);
    return 0;
}

static ST_int hc2_find(ST_int *parent, ST_int x)
{
    while (parent[x] != x)
    {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

/* Factor of S in Sd (nS x nS) with one level of each connected component
 * held at zero; dead marks those levels. */
static ST_retcode hc2_dense_factor(const hc2_system *M, ST_double *Sd, unsigned char *dead,
                                   ST_double *weakest)
{
    const ST_int nS = M->nS;
    ST_retcode rc = hc2_dense_sum(M, Sd);
    if (rc)
        return rc;
    ST_int *parent = (ST_int *)malloc((size_t)nS * sizeof(ST_int));
    unsigned char *pin = (unsigned char *)calloc((size_t)nS, 1);
    if (!parent || !pin)
    {
        free(parent);
        free(pin);
        return 920;
    }
    /* Two levels are linked when some row of B meets both. The level with
     * the largest index in each component is the one a Cholesky sweep would
     * find with a zero pivot; fixing it by name avoids a threshold. */
    for (ST_int j = 0; j < nS; j++)
        parent[j] = j;
    for (ST_int j1 = 1; j1 < nS; j1++)
    {
        for (ST_int j2 = 0; j2 < j1; j2++)
        {
            if (Sd[(size_t)j1 * nS + j2] != 0.0)
            {
                const ST_int r1 = hc2_find(parent, j1), r2 = hc2_find(parent, j2);
                if (r1 != r2)
                    parent[r1 < r2 ? r1 : r2] = r1 < r2 ? r2 : r1;
            }
        }
    }
    /* a root is the largest index of its component */
    for (ST_int j = 0; j < nS; j++)
        pin[j] = hc2_find(parent, j) == j;
    for (ST_int j1 = 0; j1 < nS; j1++)
    {
        for (ST_int j2 = 0; j2 < j1; j2++)
            Sd[(size_t)j1 * nS + j2] = -Sd[(size_t)j1 * nS + j2];
        Sd[(size_t)j1 * nS + j1] = M->dS[j1] - Sd[(size_t)j1 * nS + j1];
    }
    rc = hdfe_coef_chol_factor(Sd, nS, pin, dead, weakest);
    free(parent);
    free(pin);
    return rc;
}

/* ========================================================================
 * Approximate Cholesky of S
 * ======================================================================== */

/* Every row of B is a star centred on it. Eliminating the centre exactly
 * links all pairs of its cells; stars of up to four cells are eliminated
 * exactly and larger ones by the factorization's own sampling. NULL when
 * memory runs out, or with *oversize set when the graph has more edges than
 * the factorization can index (the plain iteration then carries on). */
static hdfe_lap *hc2_laplacian(const hc2_system *M, int *oversize)
{
    const ST_int *rowptr = M->C.rowptr, *col = M->C.col;
    const ST_double *val = M->C.val;
    long long m = 0;
    ST_int longest = 0;
    *oversize = 0;
    for (ST_int a = 0; a < M->nact; a++)
    {
        const ST_int cells = rowptr[M->act[a] + 1] - rowptr[M->act[a]];
        m += cells <= 4 ? (long long)cells * (cells - 1) / 2 : cells - 1;
        if (cells > longest)
            longest = cells;
    }
    if (m > 2147483647LL / 2)
    {
        *oversize = 1;
        return NULL;
    }
    const size_t mm = (size_t)(m > 0 ? m : 1);
    ST_int *eu = (ST_int *)malloc(mm * sizeof(ST_int));
    ST_int *ev = (ST_int *)malloc(mm * sizeof(ST_int));
    ST_double *ew = (ST_double *)malloc(mm * sizeof(ST_double));
    hdfe_lap_leaf *cs = (hdfe_lap_leaf *)malloc((size_t)(longest > 0 ? longest : 1) * sizeof(*cs));
    ST_double *cum = (ST_double *)malloc((size_t)(longest > 0 ? longest : 1) * sizeof(ST_double));
    hdfe_lap *L = NULL;
    if (eu && ev && ew && cs && cum)
    {
        uint64_t state = UINT64_C(0xD1B54A32D192ED03);
        ST_int ne = 0;
        for (ST_int a = 0; a < M->nact; a++)
        {
            const ST_int f = M->act[a], e0 = rowptr[f], cells = rowptr[f + 1] - e0;
            if (cells <= 4)
            {
                const ST_double ib = M->invB[f];
                for (ST_int x = 0; x < cells; x++)
                {
                    for (ST_int y = x + 1; y < cells; y++)
                    {
                        eu[ne] = col[e0 + x] - 1;
                        ev[ne] = col[e0 + y] - 1;
                        ew[ne] = (val ? val[e0 + x] * val[e0 + y] : 1.0) * ib;
                        ne++;
                    }
                }
                continue;
            }
            for (ST_int x = 0; x < cells; x++)
            {
                cs[x].entry = x;
                cs[x].row = col[e0 + x] - 1;
                cs[x].val = val ? val[e0 + x] : 1.0;
            }
            ne += hdfe_lap_sample_star(cs, cum, cells, &state, eu + ne, ev + ne, ew + ne);
        }
        L = hdfe_lap_factor(M->nS, ne, eu, ev, ew);
    }
    free(eu);
    free(ev);
    free(ew);
    free(cs);
    free(cum);
    return L;
}

typedef struct
{
    const hdfe_lap *L;
    ST_int n;
    int K;
} hc2_lap_ctx;

static void hc2_lap_apply(void *ctx, const ST_double *r, ST_double *z)
{
    const hc2_lap_ctx *P = (const hc2_lap_ctx *)ctx;
    memcpy(z, r, (size_t)P->n * (size_t)P->K * sizeof(ST_double));
    hdfe_lap_solve(P->L, z, P->K);
}

/* ========================================================================
 * The eliminated effect's coefficients
 * ======================================================================== */

/* a = a0 - N_B^-1 C b, in place in a0 */
static void hc2_eliminated(const hc2_system *M, int K, ST_double *a0, const ST_double *b)
{
    const ST_int nB = M->nB;
    const ST_int nbB = (nB + 4095) / 4096;
#pragma omp parallel for schedule(dynamic, 4) if ((size_t)M->C.nnz * (size_t)K >                   \
                                                      HDFE_COEF_PARALLEL_MIN)
    for (ST_int blk = 0; blk < nbB; blk++)
    {
        const ST_int f0 = blk * 4096, f1 = nB - f0 > 4096 ? f0 + 4096 : nB;
        HDFE_COEF_K_DISPATCH(K, hc2_rows_gather(M, f0, f1, b, a0, 1, KK))
    }
}

/* ========================================================================
 * Driver
 * ======================================================================== */

static hdfe_coef_method hc2_choose(const hc2_system *M, ST_int K, hdfe_coef_method method)
{
    if (method != HDFE_COEF_AUTO)
        return method;
    if (M->nS > HDFE_COEF_DENSE_MAX)
        return HDFE_COEF_PCG;
    const ST_double nS = (ST_double)M->nS;
    ST_double work = nS * nS * nS / 3.0;
    for (ST_int a = 0; a < M->nact; a++)
    {
        const ST_double cells = (ST_double)(M->C.rowptr[M->act[a] + 1] - M->C.rowptr[M->act[a]]);
        work += 2.0 * cells > nS ? 0.5 * nS * nS : cells * cells;
    }
    return work <= HDFE_COEF_DENSE_SWEEPS * (ST_double)M->nnz_act * (ST_double)K + 2e6
               ? HDFE_COEF_DENSE
               : HDFE_COEF_PCG;
}

static HDFE_SolveResult hc2_partial_out(const HDFE_State *S, ST_double *data, ST_int N, ST_int K,
                                        hdfe_coef_method method)
{
    HDFE_SolveResult result = {0, 1, 0};
    hc2_system M;
    hdfe_coef_rowplan R;
    hc2_sweep W = {&M, 0, NULL, NULL};
    ST_double *a0 = NULL, *g = NULL, *b = NULL, *pcg = NULL, *Sd = NULL, *minv = NULL;
    ST_double weakest = 1.0;
    unsigned char *dead = NULL;
    hdfe_lap *lap = NULL;
    memset(&R, 0, sizeof(R));

    /* the dense path needs the fused bookkeeping, which is all it uses */
    ST_retcode rc = hc2_system_build(&M, S, method != HDFE_COEF_DENSE);
    if (rc)
        goto done;
    const ST_int nB = M.nB, nS = M.nS;
    hdfe_coef_method chosen = hc2_choose(&M, K, method);
    if (chosen == HDFE_COEF_DENSE && M.twopass)
        chosen = HDFE_COEF_PCG;
    if (ctools_cancel_requested())
    {
        rc = 1;
        goto done;
    }

    const int KB = K < HDFE_COEF_KMAX ? (int)K : HDFE_COEF_KMAX;
    a0 = (ST_double *)malloc((size_t)nB * (size_t)KB * sizeof(ST_double));
    g = (ST_double *)malloc((size_t)nS * (size_t)KB * sizeof(ST_double));
    b = (ST_double *)malloc((size_t)nS * (size_t)KB * sizeof(ST_double));
    if (!a0 || !g || !b || (rc = hdfe_coef_rowplan_init(&R, &M.D, KB)))
    {
        if (!rc)
            rc = 920;
        goto done;
    }
    if (chosen == HDFE_COEF_DENSE)
    {
        Sd = (ST_double *)malloc((size_t)nS * (size_t)nS * sizeof(ST_double));
        dead = (unsigned char *)malloc((size_t)nS);
        if (!Sd || !dead)
        {
            rc = 920;
            goto done;
        }
        if ((rc = hc2_dense_factor(&M, Sd, dead, &weakest)))
            goto done;
    }
    else
    {
        if (M.twopass)
            W.tB = (ST_double *)malloc((size_t)nB * (size_t)KB * sizeof(ST_double));
        else
            W.buf =
                (ST_double *)malloc((size_t)M.nch * (size_t)nS * (size_t)KB * sizeof(ST_double));
        pcg = (ST_double *)malloc(HDFE_COEF_PCG_WORK(nS, KB) * sizeof(ST_double));
        minv = (ST_double *)malloc((size_t)nS * sizeof(ST_double));
        if ((!W.tB && !W.buf) || !pcg || !minv)
        {
            rc = 920;
            goto done;
        }
        /* The diagonal of S, not the level counts that the row-space sweep
         * uses implicitly: where most rows of B have a single cell the two
         * differ widely across levels, and the diagonal needs 25 to 40%
         * fewer iterations (benchmark_hdfe_coef.c, worker-firm graphs). A
         * level whose diagonal is rounding is not linked to any other. */
        const ST_double *scale = M.twopass ? M.nSw : M.dS;
        for (ST_int j = 0; j < nS; j++)
            minv[j] = (M.diag[j] > 1e-12 * scale[j] && M.diag[j] > 0.0) ? 1.0 / M.diag[j] : 0.0;
        if (chosen == HDFE_COEF_PCG_LAP)
        {
            int oversize = 0;
            lap = hc2_laplacian(&M, &oversize);
            if (!lap)
            {
                rc = oversize ? 198 : 920;
                goto done;
            }
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
            ST_int iterations = 0;
            W.K = kb;
            hdfe_coef_rows_rhs(&M.D, &R, X, kb, a0, g, f0, scale);
            if (ctools_cancel_requested())
            {
                rc = 1;
                goto done;
            }
            if (chosen == HDFE_COEF_DENSE)
            {
                memcpy(b, g, (size_t)nS * (size_t)kb * sizeof(ST_double));
                hdfe_coef_chol_solve(Sd, nS, dead, b, kb);
                hdfe_coef_direct_state(nS, kb, g, b, f0, scale, S->tolerance, weakest, state);
            }
            else
            {
                hc2_lap_ctx lap_ctx = {lap, nS, kb};
                hdfe_coef_pcg_problem P = {nS,   kb,           hc2_apply,  &W, NULL, NULL,
                                           minv, S->tolerance, S->maxiter, 0,  0};
                ST_retcode solve = 430;
                ST_int used = 0;
                int resume = lap != NULL;
                if (!lap)
                {
                    /* automatic: a short plain run settles well connected data */
                    const int probe =
                        method == HDFE_COEF_AUTO && S->maxiter > HDFE_COEF_PROBE_ITERATIONS;
                    if (probe)
                    {
                        /* half as many again when it is about to finish:
                         * the factor costs about that much */
                        P.maxiter = HDFE_COEF_PROBE_ITERATIONS;
                        P.extend = HDFE_COEF_PROBE_ITERATIONS / 2;
                        if (P.extend > S->maxiter - P.maxiter)
                            P.extend = S->maxiter - P.maxiter;
                    }
                    solve = hdfe_coef_pcg(&P, g, b, f0, scale, pcg, state, &used);
                    iterations = used;
                    P.extend = 0;
                    if (solve == 430 && probe)
                    {
                        int oversize = 0;
                        lap = hc2_laplacian(&M, &oversize);
                        if (!lap && !oversize)
                        {
                            rc = 920;
                            goto done;
                        }
                        lap_ctx.L = lap;
                        P.warm = 1;
                        resume = 1;
                    }
                }
                if (resume && solve == 430)
                {
                    /* from where the plain run stopped, with the factor when
                     * there is one; what that run used counts against the
                     * limit */
                    if (lap)
                    {
                        P.precond = hc2_lap_apply;
                        P.precond_ctx = &lap_ctx;
                        P.minv = NULL;
                    }
                    P.maxiter = S->maxiter - iterations > 0 ? S->maxiter - iterations : 1;
                    solve = hdfe_coef_pcg(&P, g, b, f0, scale, pcg, state, &used);
                    iterations += used;
                }
                if (solve == 1)
                {
                    rc = 1;
                    goto done;
                }
                if (solve == 430)
                    result.converged = 0;
            }
            /* a round counts as one iteration even when nothing was left
             * for the kept effect to explain */
            total += iterations > 1 ? iterations : 1;
            if (ctools_cancel_requested())
            {
                rc = 1;
                goto done;
            }
            hc2_eliminated(&M, kb, a0, b);
            hdfe_coef_rows_residual(&M.D, X, kb, a0, b);
            if (!hdfe_coef_refine_span(state, kb, &first, &kb))
                break;
        }
        if (total > result.iterations)
            result.iterations = total;
    }

done:
    hdfe_lap_free(lap);
    hdfe_coef_rowplan_free(&R);
    free(W.buf);
    free(W.tB);
    free(a0);
    free(g);
    free(b);
    free(pcg);
    free(Sd);
    free(dead);
    free(minv);
    hc2_system_free(&M);
    if (rc)
    {
        result.status = rc;
        result.converged = 0;
    }
    else if (!result.converged)
        result.status = 430;
    return result;
}

/* ========================================================================
 * Entry point
 * ======================================================================== */

int hdfe_coef_applies(const HDFE_State *S)
{
    if (!S || !S->factors || S->G < 2 || S->G > HDFE_COEF_MAX_EFFECTS || S->N < 1)
        return 0;
    for (ST_int g = 0; g < S->G; g++)
    {
        const FE_Factor *f = &S->factors[g];
        if (f->slope || !f->levels || f->num_levels < 1)
            return 0;
    }
    return 1;
}

HDFE_SolveResult hdfe_coef_partial_out(const HDFE_State *S, ST_double *data, ST_int N, ST_int K,
                                       hdfe_coef_method method)
{
    if (!hdfe_coef_applies(S) || !data || N != S->N || K < 1 || S->maxiter < 1 ||
        !isfinite(S->tolerance) || S->tolerance <= 0)
        return (HDFE_SolveResult){198, 0, 0};
    if (S->G > 2)
        return hdfe_coef_multi_partial_out(S, data, N, K, method);
    if (method == HDFE_COEF_DENSE && S->factors[0].num_levels > 46340 &&
        S->factors[1].num_levels > 46340)
        return (HDFE_SolveResult){198, 0, 0};
    return hc2_partial_out(S, data, N, K, method);
}

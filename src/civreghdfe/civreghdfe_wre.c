/*
 * civreghdfe_wre.c
 *
 * The WRE bootstrap of civreghdfe's k-class estimators (see
 * civreghdfe_wre.h). Matrices are column-major unless noted; the basis rows
 * are kept row-major ([y, Y2, X] of row i at i * nb) for the row loops.
 * Part of the ctools Stata plugin suite
 */

#include "civreghdfe_wre.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "civreghdfe_matrix.h"
#include "hdfe/hdfe_wild.h"

struct civ_wre_model_s
{
    ST_int N, K2, K1, Kz, L, K, F, E, nb, iX;
    int est;
    double kclass, fuller_adj;
    int vce, null_imposed, projected;
    ST_int G, C, H;
    int score_agg; /* robust or cluster: sums per (cluster x bootstrap cluster) group */
    int aggc;      /* nested: sums per bootstrap cluster for X'WQ and Q'WQ */
    ST_int n_test;
    ST_int *test;
    ST_double *Br; /* N x nb, row-major */
    ST_double *w;  /* N */
    ST_int *row_c; /* N bootstrap cluster */
    /* rows by vce() cluster g, and by group h (g x bootstrap cluster) */
    ST_int *rows;    /* N */
    ST_int *g_hs;    /* G + 1: groups of g are g_hs[g] .. g_hs[g+1]-1 */
    ST_int *h_rs;    /* H + 1: rows of h are rows[h_rs[h]] .. */
    ST_int *h_c;     /* H */
    ST_double *XBh;  /* H x (L x nb) */
    ST_double *XXg;  /* G x (L x L) */
    ST_double *XBc;  /* C x (L x nb) */
    ST_double *BBc;  /* C x (nb x nb) */
    ST_double *BB;   /* nb x nb */
    ST_double *XX, *XXinv, *XX11inv, *Pi0;
    /* weights of replicates 1 .. B (stored at r - 1) */
    ST_int B, got;
    uint8_t *vcode;
    double vtab[16];
    int ntab;
    double *vdbl;
    /* per replicate 0 .. B: X'WQ (L x nb) and Q'WQ (nb x nb) */
    ST_double *XQ, *QQ;
    int prepared, cached;
    ST_double *XQg; /* projected, score_agg: (B+1) x G x (L x nb) */
    ST_double *Qr;  /* projected, rows: (B+1) x N x nb, row-major */
};

/* ------------------------------------------------------------------------
 * Small dense algebra
 * ------------------------------------------------------------------------ */

/* Cholesky in place (lower triangle); 1 when not positive definite */
static int wre_chol(ST_int n, double *A)
{
    for (ST_int j = 0; j < n; j++)
    {
        const double d0 = A[j + (size_t)j * n];
        double s = d0;
        for (ST_int k = 0; k < j; k++)
            s -= A[j + (size_t)k * n] * A[j + (size_t)k * n];
        if (!(s > 1e-13 * fabs(d0)) || !(s > 0.0))
            return 1;
        const double l = sqrt(s);
        A[j + (size_t)j * n] = l;
        for (ST_int i = j + 1; i < n; i++)
        {
            double t = A[i + (size_t)j * n];
            for (ST_int k = 0; k < j; k++)
                t -= A[i + (size_t)k * n] * A[j + (size_t)k * n];
            A[i + (size_t)j * n] = t / l;
        }
    }
    return 0;
}

static void wre_chol_solve(ST_int n, const double *Lc, double *b)
{
    for (ST_int i = 0; i < n; i++)
    {
        double t = b[i];
        for (ST_int k = 0; k < i; k++)
            t -= Lc[i + (size_t)k * n] * b[k];
        b[i] = t / Lc[i + (size_t)i * n];
    }
    for (ST_int i = n - 1; i >= 0; i--)
    {
        double t = b[i];
        for (ST_int k = i + 1; k < n; k++)
            t -= Lc[k + (size_t)i * n] * b[k];
        b[i] = t / Lc[i + (size_t)i * n];
    }
}

/* Inverse of a symmetric positive definite matrix; 1 when not */
static int wre_spd_inverse(ST_int n, const double *A, double *inv, double *work)
{
    memcpy(work, A, (size_t)n * n * sizeof(double));
    if (wre_chol(n, work))
        return 1;
    for (ST_int j = 0; j < n; j++)
    {
        double *c = inv + (size_t)j * n;
        for (ST_int i = 0; i < n; i++)
            c[i] = i == j;
        wre_chol_solve(n, work, c);
    }
    for (ST_int j = 0; j < n; j++)
        for (ST_int i = 0; i < j; i++)
            inv[i + (size_t)j * n] = inv[j + (size_t)i * n] =
                0.5 * (inv[i + (size_t)j * n] + inv[j + (size_t)i * n]);
    return 0;
}

/*
 * The smallest root of |A - lambda Bm| = 0 for symmetric A and positive
 * definite Bm (n x n, overwritten): LIML's lambda. Equilibrated by
 * diag(Bm)^(-1/2) first, as civreghdfe computes lambda. NaN on failure.
 */
static double wre_min_geneig(ST_int n, double *A, double *Bm, double *work)
{
    for (ST_int i = 0; i < n; i++)
        if (!(Bm[i + (size_t)i * n] > 0.0))
            return NAN;
    double ds[64];
    if (n > 64)
        return NAN;
    for (ST_int i = 0; i < n; i++)
        ds[i] = 1.0 / sqrt(Bm[i + (size_t)i * n]);
    for (ST_int j = 0; j < n; j++)
        for (ST_int i = 0; i < n; i++)
        {
            A[i + (size_t)j * n] *= ds[i] * ds[j];
            Bm[i + (size_t)j * n] *= ds[i] * ds[j];
        }
    if (wre_chol(n, Bm))
        return NAN;
    /* C = L^-1 A L^-T: columns of L^-1 A, then rows */
    double *Cm = work;
    memcpy(Cm, A, (size_t)n * n * sizeof(double));
    for (ST_int j = 0; j < n; j++) /* forward solve each column */
    {
        double *c = Cm + (size_t)j * n;
        for (ST_int i = 0; i < n; i++)
        {
            double t = c[i];
            for (ST_int k = 0; k < i; k++)
                t -= Bm[i + (size_t)k * n] * c[k];
            c[i] = t / Bm[i + (size_t)i * n];
        }
    }
    for (ST_int i = 0; i < n; i++) /* then each row: (L^-1 (L^-1 A)')' */
    {
        for (ST_int j = 0; j < n; j++)
        {
            double t = Cm[i + (size_t)j * n];
            for (ST_int k = 0; k < j; k++)
                t -= Bm[j + (size_t)k * n] * Cm[i + (size_t)k * n];
            Cm[i + (size_t)j * n] = t / Bm[j + (size_t)j * n];
        }
    }
    if (n == 1)
        return Cm[0];
    if (n == 2)
    {
        const double a = Cm[0], d = Cm[3], b = 0.5 * (Cm[1] + Cm[2]);
        const double h = 0.5 * (a - d);
        return 0.5 * (a + d) - sqrt(h * h + b * b);
    }
    for (ST_int j = 0; j < n; j++)
        for (ST_int i = 0; i < j; i++)
            Cm[i + (size_t)j * n] = Cm[j + (size_t)i * n] =
                0.5 * (Cm[i + (size_t)j * n] + Cm[j + (size_t)i * n]);
    return civreghdfe_jacobi_eigenvalues(Cm, n, NULL);
}

/* ------------------------------------------------------------------------
 * Setup
 * ------------------------------------------------------------------------ */

void civ_wre_free(civ_wre_model *m)
{
    if (!m)
        return;
    free(m->test);
    free(m->Br);
    free(m->w);
    free(m->row_c);
    free(m->rows);
    free(m->g_hs);
    free(m->h_rs);
    free(m->h_c);
    free(m->XBh);
    free(m->XXg);
    free(m->XBc);
    free(m->BBc);
    free(m->BB);
    free(m->XX);
    free(m->XXinv);
    free(m->XX11inv);
    free(m->Pi0);
    free(m->vcode);
    free(m->vdbl);
    free(m->XQ);
    free(m->QQ);
    free(m->XQg);
    free(m->Qr);
    free(m);
}

ST_int civ_wre_B(const civ_wre_model *m) { return m ? m->B : 0; }
ST_int civ_wre_C(const civ_wre_model *m) { return m ? m->C : 0; }
int civ_wre_cached(const civ_wre_model *m) { return m && m->cached; }

typedef struct
{
    int64_t key;
    ST_int i;
} wre_key;

static int wre_key_cmp(const void *a, const void *b)
{
    const wre_key *x = a, *y = b;
    if (x->key != y->key)
        return x->key < y->key ? -1 : 1;
    return (x->i > y->i) - (x->i < y->i);
}

/* acc (r x c, column-major) += f * a b' for row vectors a (r) and b (c) */
static inline void wre_outer(double *acc, ST_int r, ST_int c, double f, const double *a,
                             const double *b)
{
    for (ST_int j = 0; j < c; j++)
    {
        const double fb = f * b[j];
        double *col = acc + (size_t)j * r;
        for (ST_int i = 0; i < r; i++)
            col[i] += a[i] * fb;
    }
}

ST_retcode civ_wre_setup(const civ_wre_input *in, ST_int B, civ_wre_model **out)
{
    *out = NULL;
    const ST_int N = in->N, K2 = in->K2, K1 = in->K1, Kz = in->Kz;
    if (N < 1 || K2 < 1 || K1 < 0 || Kz < 0 || K1 + Kz < K2 + K1 || B < 1 || in->n_test < 1 ||
        !in->y || !in->Y2 || (K1 && !in->X1) || (Kz && !in->X2) || in->vce < 0 || in->vce > 2 ||
        in->est < 0 || in->est > 3 || (in->vce == 2 && (!in->vce_id || in->G < 1)) ||
        (in->boot_id && in->C < 1) || K1 + Kz + K2 + 1 > 4096)
        return 198;
    civ_wre_model *m = calloc(1, sizeof(*m));
    if (!m)
        return 920;
    m->N = N;
    m->K2 = K2;
    m->K1 = K1;
    m->Kz = Kz;
    m->L = K1 + Kz;
    m->K = K2 + K1;
    m->F = 1 + m->K;
    m->E = 1 + K2;
    m->nb = 1 + K2 + m->L;
    m->iX = 1 + K2;
    m->est = in->est;
    m->kclass = in->kclass;
    m->fuller_adj = in->fuller_adj;
    m->vce = in->vce;
    m->null_imposed = in->null_imposed;
    m->projected = in->projected;
    m->B = B;
    m->n_test = in->n_test;
    const ST_int L = m->L, nb = m->nb, iX = m->iX;
    ST_retcode rc = 0;

    m->test = malloc((size_t)in->n_test * sizeof(ST_int));
    m->Br = malloc((size_t)N * nb * sizeof(ST_double));
    m->w = malloc((size_t)N * sizeof(ST_double));
    m->row_c = malloc((size_t)N * sizeof(ST_int));
    m->rows = malloc((size_t)N * sizeof(ST_int));
    m->BB = calloc((size_t)nb * nb, sizeof(ST_double));
    m->XX = malloc((size_t)L * L * sizeof(ST_double));
    m->XXinv = malloc((size_t)L * L * sizeof(ST_double));
    m->XX11inv = malloc((size_t)(K1 > 0 ? K1 * K1 : 1) * sizeof(ST_double));
    m->Pi0 = malloc((size_t)L * K2 * sizeof(ST_double));
    wre_key *keys = malloc((size_t)N * sizeof(wre_key));
    double *work = malloc((size_t)L * L * sizeof(double));
    if (!m->test || !m->Br || !m->w || !m->row_c || !m->rows || !m->BB || !m->XX || !m->XXinv ||
        !m->XX11inv || !m->Pi0 || !keys || !work)
    {
        rc = 920;
        goto fail;
    }
    for (ST_int j = 0; j < in->n_test; j++)
    {
        if (in->test[j] < 0 || in->test[j] >= m->K)
        {
            rc = 198;
            goto fail;
        }
        m->test[j] = in->test[j];
    }
    for (ST_int i = 0; i < N; i++)
    {
        double *b = m->Br + (size_t)i * nb;
        b[0] = in->y[i];
        for (ST_int j = 0; j < K2; j++)
            b[1 + j] = in->Y2[(size_t)j * N + i];
        for (ST_int j = 0; j < K1; j++)
            b[iX + j] = in->X1[(size_t)j * N + i];
        for (ST_int j = 0; j < Kz; j++)
            b[iX + K1 + j] = in->X2[(size_t)j * N + i];
        m->w[i] = in->w ? in->w[i] : 1.0;
        if (!(m->w[i] >= 0.0) || !isfinite(m->w[i]))
        {
            rc = 198;
            goto fail;
        }
    }

    /* Clusters: bootstrap (each row without ids) and vce() (one group for
     * the unadjusted variance, each row for robust) */
    m->C = in->boot_id ? in->C : N;
    m->G = in->vce == 2 ? in->G : in->vce == 1 ? N : 1;
    for (ST_int i = 0; i < N; i++)
    {
        const ST_int c = in->boot_id ? in->boot_id[i] : i;
        const ST_int g = in->vce == 2 ? in->vce_id[i] : in->vce == 1 ? i : 0;
        if (c < 0 || c >= m->C || g < 0 || g >= m->G)
        {
            rc = 198;
            goto fail;
        }
        m->row_c[i] = c;
        keys[i].key = (int64_t)g * (int64_t)m->C + c;
        keys[i].i = i;
    }
    qsort(keys, (size_t)N, sizeof(wre_key), wre_key_cmp);
    ST_int H = 0;
    for (ST_int r = 0; r < N; r++)
    {
        m->rows[r] = keys[r].i;
        H += (r == 0 || keys[r].key != keys[r - 1].key);
    }
    m->H = H;
    m->g_hs = calloc((size_t)m->G + 1, sizeof(ST_int));
    m->h_rs = malloc((size_t)(H + 1) * sizeof(ST_int));
    m->h_c = malloc((size_t)H * sizeof(ST_int));
    if (!m->g_hs || !m->h_rs || !m->h_c)
    {
        rc = 920;
        goto fail;
    }
    for (ST_int r = 0, h = -1; r < N; r++)
    {
        if (r == 0 || keys[r].key != keys[r - 1].key)
        {
            h++;
            m->h_rs[h] = r;
            m->h_c[h] = (ST_int)(keys[r].key % m->C);
            m->g_hs[(ST_int)(keys[r].key / m->C) + 1]++;
        }
    }
    m->h_rs[H] = N;
    for (ST_int g = 0; g < m->G; g++)
        m->g_hs[g + 1] += m->g_hs[g];

    /* Cross products of the basis */
    for (ST_int i = 0; i < N; i++)
    {
        const double *b = m->Br + (size_t)i * nb;
        wre_outer(m->BB, nb, nb, m->w[i], b, b);
    }
    for (ST_int j = 0; j < L; j++)
        for (ST_int i = 0; i < L; i++)
            m->XX[i + (size_t)j * L] = m->BB[(iX + i) + (size_t)(iX + j) * nb];
    if (wre_spd_inverse(L, m->XX, m->XXinv, work))
    {
        rc = 198;
        goto fail;
    }
    if (K1 > 0)
    {
        double *x11 = malloc((size_t)K1 * K1 * sizeof(double));
        if (!x11)
        {
            rc = 920;
            goto fail;
        }
        for (ST_int j = 0; j < K1; j++)
            for (ST_int i = 0; i < K1; i++)
                x11[i + (size_t)j * K1] = m->XX[i + (size_t)j * L];
        int bad = wre_spd_inverse(K1, x11, m->XX11inv, work);
        free(x11);
        if (bad)
        {
            rc = 198;
            goto fail;
        }
    }
    for (ST_int j = 0; j < K2; j++)
        for (ST_int l = 0; l < L; l++)
        {
            double s = 0;
            for (ST_int q = 0; q < L; q++)
                s += m->XXinv[l + (size_t)q * L] * m->BB[(iX + q) + (size_t)(1 + j) * nb];
            m->Pi0[l + (size_t)j * L] = s;
        }

    /* Sums per group for the sandwich, and per bootstrap cluster for the
     * replicates' cross products on the nested path */
    m->score_agg = in->vce >= 1 && (size_t)H * (size_t)L <= (size_t)N;
    m->aggc = !in->projected && (size_t)m->C * (size_t)(L + nb) <= (size_t)N;
    if (m->score_agg)
    {
        m->XXg = calloc((size_t)m->G * L * L, sizeof(ST_double));
        if (!in->projected)
            m->XBh = calloc((size_t)H * L * nb, sizeof(ST_double));
        if (!m->XXg || (!in->projected && !m->XBh))
        {
            rc = 920;
            goto fail;
        }
        for (ST_int g = 0; g < m->G; g++)
            for (ST_int h = m->g_hs[g]; h < m->g_hs[g + 1]; h++)
                for (ST_int r = m->h_rs[h]; r < m->h_rs[h + 1]; r++)
                {
                    const ST_int i = m->rows[r];
                    const double *b = m->Br + (size_t)i * nb;
                    wre_outer(m->XXg + (size_t)g * L * L, L, L, m->w[i], b + iX, b + iX);
                    if (m->XBh)
                        wre_outer(m->XBh + (size_t)h * L * nb, L, nb, m->w[i], b + iX, b);
                }
    }
    if (m->aggc)
    {
        m->XBc = calloc((size_t)m->C * L * nb, sizeof(ST_double));
        m->BBc = calloc((size_t)m->C * nb * nb, sizeof(ST_double));
        if (!m->XBc || !m->BBc)
        {
            rc = 920;
            goto fail;
        }
        for (ST_int i = 0; i < N; i++)
        {
            const double *b = m->Br + (size_t)i * nb;
            const ST_int c = m->row_c[i];
            wre_outer(m->XBc + (size_t)c * L * nb, L, nb, m->w[i], b + iX, b);
            wre_outer(m->BBc + (size_t)c * nb * nb, nb, nb, m->w[i], b, b);
        }
    }
    free(keys);
    free(work);
    *out = m;
    return 0;

fail:
    free(keys);
    free(work);
    civ_wre_free(m);
    return rc;
}

/* ------------------------------------------------------------------------
 * Weights
 * ------------------------------------------------------------------------ */

static inline double wre_v(const civ_wre_model *m, ST_int r, ST_int c)
{
    if (r == 0)
        return 1.0;
    const size_t at = (size_t)(r - 1) * (size_t)m->C + (size_t)c;
    return m->vdbl ? m->vdbl[at] : m->vtab[m->vcode[at]];
}

ST_retcode civ_wre_weights(civ_wre_model *m, const ST_double *v, ST_int first, ST_int count,
                           size_t budget)
{
    if (!m || first != m->got || count < 1 || first + count > m->B)
        return 198;
    /* enumeration (v NULL) generates sign vectors 0 .. 2^C - 1: B must be 2^C */
    if (!v && hdfe_wild_enumerable(m->C, m->B) != m->B)
        return 198;
    const ST_int C = m->C;
    const size_t total = (size_t)m->B * (size_t)C;
    if (!m->vcode && !m->vdbl)
    {
        if (total > budget)
            return 909;
        if (!(m->vcode = malloc(total)))
            return 920;
    }
    double *e = NULL;
    if (!v && !(e = malloc((size_t)C * sizeof(double))))
        return 920;
    for (ST_int r = 0; r < count; r++)
    {
        if (e)
            hdfe_wild_enum_weights(C, (uint64_t)(first + r), e);
        const double *row = e ? e : v + (size_t)r * C;
        for (ST_int c = 0; c < C; c++)
        {
            const size_t at = (size_t)(first + r) * (size_t)C + (size_t)c;
            const double x = row[c];
            if (!isfinite(x))
            {
                free(e);
                return 198;
            }
            if (m->vdbl)
            {
                m->vdbl[at] = x;
                continue;
            }
            int q = 0;
            while (q < m->ntab && m->vtab[q] != x)
                q++;
            if (q == m->ntab && m->ntab < 16)
                m->vtab[m->ntab++] = x;
            if (q < m->ntab)
            {
                m->vcode[at] = (uint8_t)q;
                continue;
            }
            /* continuous weights: keep doubles */
            if (total > budget / sizeof(double))
            {
                free(e);
                return 909;
            }
            double *d = malloc(total * sizeof(double));
            if (!d)
            {
                free(e);
                return 920;
            }
            for (size_t s = 0; s < at; s++)
                d[s] = m->vtab[m->vcode[s]];
            d[at] = x;
            free(m->vcode);
            m->vcode = NULL;
            m->vdbl = d;
        }
    }
    free(e);
    m->got += count;
    return 0;
}

/* ------------------------------------------------------------------------
 * Replicates' cross products
 * ------------------------------------------------------------------------ */

/* Nested path: X'WQ and Q'WQ of replicate r with Q = basis .* v */
static void wre_rep_nested(const civ_wre_model *m, ST_int r, double *XQ, double *QQ)
{
    const ST_int L = m->L, nb = m->nb, iX = m->iX;
    memset(XQ, 0, (size_t)L * nb * sizeof(double));
    memset(QQ, 0, (size_t)nb * nb * sizeof(double));
    if (m->aggc)
    {
        for (ST_int c = 0; c < m->C; c++)
        {
            const double v = wre_v(m, r, c);
            const double *xb = m->XBc + (size_t)c * L * nb, *bb = m->BBc + (size_t)c * nb * nb;
            for (size_t q = 0; q < (size_t)L * nb; q++)
                XQ[q] += v * xb[q];
            const double v2 = v * v;
            for (size_t q = 0; q < (size_t)nb * nb; q++)
                QQ[q] += v2 * bb[q];
        }
        return;
    }
    for (ST_int i = 0; i < m->N; i++)
    {
        const double *b = m->Br + (size_t)i * nb;
        const double v = wre_v(m, r, m->row_c[i]);
        wre_outer(XQ, L, nb, m->w[i] * v, b + iX, b);
        wre_outer(QQ, nb, nb, m->w[i] * v * v, b, b);
    }
}

/* The perturbed basis of replicates r0 .. r0 + nr - 1, projected:
 * buf (N x nr nb, column-major; replicate rr's columns at rr nb) */
static ST_retcode wre_project_batch(const civ_wre_model *m, ST_int r0, ST_int nr, double *buf,
                                    civ_wre_project_fn project, void *ctx)
{
    const ST_int N = m->N, nb = m->nb;
    for (ST_int rr = 0; rr < nr; rr++)
        for (ST_int j = 0; j < nb; j++)
        {
            double *col = buf + ((size_t)rr * nb + j) * N;
            for (ST_int i = 0; i < N; i++)
                col[i] = m->Br[(size_t)i * nb + j] * wre_v(m, r0 + rr, m->row_c[i]);
        }
    return project(ctx, buf, nr * nb);
}

/* X'WQ, Q'WQ (and the cluster sums X_g'W_g Q_g when kept) from projected
 * columns Q (N x nb at q, column-major) */
static void wre_rep_from_q(const civ_wre_model *m, const double *q, double *XQ, double *QQ,
                           double *XQg, double *Qrow)
{
    const ST_int N = m->N, L = m->L, nb = m->nb, iX = m->iX;
    double qi[4096];
    memset(XQ, 0, (size_t)L * nb * sizeof(double));
    memset(QQ, 0, (size_t)nb * nb * sizeof(double));
    if (XQg)
        memset(XQg, 0, (size_t)m->G * L * nb * sizeof(double));
    for (ST_int g = 0; g < m->G; g++)
        for (ST_int h = m->g_hs[g]; h < m->g_hs[g + 1]; h++)
            for (ST_int s = m->h_rs[h]; s < m->h_rs[h + 1]; s++)
            {
                const ST_int i = m->rows[s];
                const double *b = m->Br + (size_t)i * nb;
                for (ST_int j = 0; j < nb; j++)
                    qi[j] = q[(size_t)j * N + i];
                wre_outer(XQ, L, nb, m->w[i], b + iX, qi);
                wre_outer(QQ, nb, nb, m->w[i], qi, qi);
                if (XQg)
                    wre_outer(XQg + (size_t)g * L * nb, L, nb, m->w[i], b + iX, qi);
                if (Qrow)
                    memcpy(Qrow + (size_t)i * nb, qi, (size_t)nb * sizeof(double));
            }
}

static int wre_threads(ST_int threads)
{
    return threads > 1 ? (int)threads : 1;
}

ST_retcode civ_wre_prepare(civ_wre_model *m, ST_int threads, civ_wre_project_fn project,
                           void *ctx, size_t budget)
{
    if (!m || m->got != m->B)
        return 198;
    const ST_int N = m->N, L = m->L, nb = m->nb, R = m->B + 1;
    free(m->XQ);
    free(m->QQ);
    free(m->XQg);
    free(m->Qr);
    m->XQg = m->Qr = NULL;
    m->XQ = malloc((size_t)R * L * nb * sizeof(double));
    m->QQ = malloc((size_t)R * nb * nb * sizeof(double));
    if (!m->XQ || !m->QQ)
        return 920;
    const int P = wre_threads(threads);
    (void)P;
    if (!m->projected)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(P) if (P > 1)
#endif
        for (ST_int r = 0; r < R; r++)
            wre_rep_nested(m, r, m->XQ + (size_t)r * L * nb, m->QQ + (size_t)r * nb * nb);
        m->cached = 1;
        m->prepared = 1;
        return 0;
    }
    if (!project)
        return 198;
    /* The projected basis is kept, per cluster or per row, when it fits */
    const size_t per_rep = m->vce == 0 ? 0
                           : m->score_agg ? (size_t)m->G * L * nb
                                          : (size_t)N * nb;
    /* Caches and batches are sized on the soft size; budget only bounds
     * one replicate's projection */
    const size_t soft = budget < CIV_WRE_SOFT_BYTES ? budget : CIV_WRE_SOFT_BYTES;
    m->cached = per_rep == 0 || (double)per_rep * R * sizeof(double) <= 0.5 * (double)soft;
    if (m->cached && per_rep > 0)
    {
        if (m->score_agg)
            m->XQg = malloc((size_t)R * per_rep * sizeof(double));
        else
            m->Qr = malloc((size_t)R * per_rep * sizeof(double));
        if (!m->XQg && !m->Qr)
            m->cached = 0;
    }
    const size_t col_bytes = (size_t)N * nb * sizeof(double);
    if (col_bytes > budget)
        return 909;
    ST_int batch = (ST_int)((soft / 4) / col_bytes);
    if (batch < 1)
        batch = 1;
    if (batch > R)
        batch = R;
    double *buf = malloc((size_t)batch * col_bytes);
    if (!buf)
        return 920;
    ST_retcode rc = 0;
    for (ST_int r0 = 0; r0 < R && !rc; r0 += batch)
    {
        const ST_int nr = (r0 + batch <= R) ? batch : R - r0;
        if ((rc = wre_project_batch(m, r0, nr, buf, project, ctx)))
            break;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(P) if (P > 1)
#endif
        for (ST_int rr = 0; rr < nr; rr++)
        {
            const ST_int r = r0 + rr;
            wre_rep_from_q(m, buf + (size_t)rr * nb * N, m->XQ + (size_t)r * L * nb,
                           m->QQ + (size_t)r * nb * nb,
                           m->XQg ? m->XQg + (size_t)r * per_rep : NULL,
                           m->Qr ? m->Qr + (size_t)r * per_rep : NULL);
        }
    }
    free(buf);
    if (!rc)
        m->prepared = 1;
    return rc;
}

/* ------------------------------------------------------------------------
 * The bootstrap DGP at one point
 * ------------------------------------------------------------------------ */

typedef struct
{
    int ok;
    double *ThF; /* L x F: the bootstrap columns [y*, Z*] = X ThF + Q CF */
    double *CF;  /* nb x E (the remaining columns of CF are 0) */
    double *XXT; /* L x F: XX ThF */
    double *G0;  /* F x F: ThF' XX ThF */
} wre_dgp;

/* basis position of regressor r of Z = [Y2, X1] */
static inline ST_int wre_bas(const civ_wre_model *m, ST_int r)
{
    return r < m->K2 ? 1 + r : m->iX + (r - m->K2);
}

static ST_retcode wre_dgp_build(const civ_wre_model *m, ST_int k, double b0, wre_dgp *d)
{
    const ST_int L = m->L, K = m->K, K2 = m->K2, K1 = m->K1, nb = m->nb, iX = m->iX, F = m->F,
                 E = m->E;
    const int null = m->null_imposed;
    const ST_int Kr = K - (null ? 1 : 0), n1 = 1 + Kr;
    d->ok = 0;
    ST_int regs[4096], is_end[4096];
    for (ST_int r = 0, q = 0; r < K; r++)
    {
        if (null && r == k)
            continue;
        regs[q] = r;
        is_end[q] = r < K2;
        q++;
    }
    size_t need = (size_t)nb * n1 + (size_t)n1 * n1 * 3 + (size_t)L * n1 * 2 +
                  (size_t)n1 * 6 + (size_t)nb * 2 + (size_t)L * (K2 + 4) + (size_t)K * 2 +
                  (size_t)n1 * n1 * 4 + 64;
    double *ws = malloc(need * sizeof(double));
    if (!ws)
        return 920;
    double *p = ws;
    double *V = p;   p += (size_t)nb * n1;  /* coordinates of [dep, regs] */
    double *GV = p;  p += (size_t)n1 * n1;  /* V'BB V */
    double *XV = p;  p += (size_t)L * n1;   /* X'W V */
    double *IXV = p; p += (size_t)L * n1;   /* XXinv X'W V */
    double *PV = p;  p += (size_t)n1 * n1;  /* V' P_X V */
    double *H = p;   p += (size_t)n1 * n1;
    double *rhs = p; p += (size_t)n1;
    double *bt = p;  p += (size_t)K;        /* restricted coefficients, Z order */
    double *cu = p;  p += (size_t)nb;
    double *xbc = p; p += (size_t)L;
    double *gx = p;  p += (size_t)L;
    double *Pi = p;  p += (size_t)L * K2;
    double *c1 = p;  p += (size_t)nb;
    double *Aw = p;  p += (size_t)n1 * n1;
    double *Bw = p;  p += (size_t)n1 * n1;
    double *Wk = p;  p += (size_t)n1 * n1 * 2;
    ST_retcode rc = 0;

    memset(V, 0, (size_t)nb * n1 * sizeof(double));
    V[0] = 1.0;
    if (null)
        V[wre_bas(m, k)] -= b0;
    for (ST_int q = 0; q < Kr; q++)
        V[(size_t)(1 + q) * nb + wre_bas(m, regs[q])] = 1.0;
    for (ST_int b = 0; b < n1; b++)
    {
        const double *vb = V + (size_t)b * nb;
        double tmp[4096];
        for (ST_int i = 0; i < nb; i++)
        {
            double s = 0;
            for (ST_int j = 0; j < nb; j++)
                s += m->BB[i + (size_t)j * nb] * vb[j];
            tmp[i] = s;
        }
        for (ST_int a = 0; a < n1; a++)
        {
            double s = 0;
            const double *va = V + (size_t)a * nb;
            for (ST_int i = 0; i < nb; i++)
                s += va[i] * tmp[i];
            GV[a + (size_t)b * n1] = s;
        }
        for (ST_int l = 0; l < L; l++)
            XV[l + (size_t)b * L] = tmp[iX + l];
    }
    for (ST_int b = 0; b < n1; b++)
        for (ST_int l = 0; l < L; l++)
        {
            double s = 0;
            for (ST_int q = 0; q < L; q++)
                s += m->XXinv[l + (size_t)q * L] * XV[q + (size_t)b * L];
            IXV[l + (size_t)b * L] = s;
        }
    for (ST_int b = 0; b < n1; b++)
        for (ST_int a = 0; a < n1; a++)
        {
            double s = 0;
            for (ST_int l = 0; l < L; l++)
                s += XV[l + (size_t)a * L] * IXV[l + (size_t)b * L];
            PV[a + (size_t)b * n1] = s;
        }

    /* LIML's kappa for the DGP; 1 (OLS) without endogenous regressors */
    ST_int ne = 0, nx = 0, ie[4096], ix[4096];
    ie[ne++] = 0;
    for (ST_int q = 0; q < Kr; q++)
    {
        if (is_end[q])
            ie[ne++] = 1 + q;
        else
            ix[nx++] = 1 + q;
    }
    double kappa = 1.0;
    if (ne > 1)
    {
        /* A = Y'M_X1 Y: Y'Y less its projection on the exogenous columns */
        double *Gx = Wk, *Gxi = Wk + (size_t)nx * nx;
        for (ST_int b = 0; b < ne; b++)
            for (ST_int a = 0; a < ne; a++)
            {
                Aw[a + (size_t)b * ne] = GV[ie[a] + (size_t)ie[b] * n1];
                Bw[a + (size_t)b * ne] = GV[ie[a] + (size_t)ie[b] * n1] - PV[ie[a] + (size_t)ie[b] * n1];
            }
        if (nx > 0)
        {
            double *tmpw = malloc((size_t)nx * nx * 2 * sizeof(double) + (size_t)nx * ne * sizeof(double));
            if (!tmpw)
            {
                rc = 920;
                goto done;
            }
            double *Gxx = tmpw, *Gxe = tmpw + (size_t)nx * nx;
            (void)Gx;
            (void)Gxi;
            for (ST_int b = 0; b < nx; b++)
                for (ST_int a = 0; a < nx; a++)
                    Gxx[a + (size_t)b * nx] = GV[ix[a] + (size_t)ix[b] * n1];
            if (wre_chol(nx, Gxx))
            {
                free(tmpw);
                goto done;
            }
            for (ST_int b = 0; b < ne; b++)
            {
                double *col = Gxe + (size_t)b * nx;
                for (ST_int a = 0; a < nx; a++)
                    col[a] = GV[ix[a] + (size_t)ie[b] * n1];
            }
            for (ST_int b = 0; b < ne; b++)
            {
                double sol[4096];
                memcpy(sol, Gxe + (size_t)b * nx, (size_t)nx * sizeof(double));
                wre_chol_solve(nx, Gxx, sol);
                for (ST_int a = 0; a < ne; a++)
                {
                    double s = 0;
                    for (ST_int q = 0; q < nx; q++)
                        s += Gxe[q + (size_t)a * nx] * sol[q];
                    Aw[a + (size_t)b * ne] -= s;
                }
            }
            free(tmpw);
        }
        kappa = wre_min_geneig(ne, Aw, Bw, Wk);
        if (!isfinite(kappa))
            goto done;
    }
    for (ST_int b = 0; b < Kr; b++)
    {
        for (ST_int a = 0; a < Kr; a++)
            H[a + (size_t)b * Kr] = (1.0 - kappa) * GV[(1 + a) + (size_t)(1 + b) * n1] +
                                    kappa * PV[(1 + a) + (size_t)(1 + b) * n1];
        rhs[b] = (1.0 - kappa) * GV[1 + b] + kappa * PV[1 + b];
    }
    if (Kr > 0)
    {
        if (wre_chol(Kr, H))
            goto done;
        wre_chol_solve(Kr, H, rhs);
    }
    for (ST_int r = 0; r < K; r++)
        bt[r] = 0;
    if (null)
        bt[k] = b0;
    for (ST_int q = 0; q < Kr; q++)
        bt[regs[q]] = rhs[q];

    /* Structural residual u = y - Z bt; the efficient first stage: Y2 on
     * X and u (FWL), Pi = Pi0 - (XXinv X'Wu) d */
    memset(cu, 0, (size_t)nb * sizeof(double));
    cu[0] = 1.0;
    for (ST_int r = 0; r < K; r++)
        cu[wre_bas(m, r)] -= bt[r];
    double uu = 0;
    {
        double tmp[4096];
        for (ST_int i = 0; i < nb; i++)
        {
            double s = 0;
            for (ST_int j = 0; j < nb; j++)
                s += m->BB[i + (size_t)j * nb] * cu[j];
            tmp[i] = s;
        }
        for (ST_int i = 0; i < nb; i++)
            uu += cu[i] * tmp[i];
        for (ST_int l = 0; l < L; l++)
            xbc[l] = tmp[iX + l];
        for (ST_int l = 0; l < L; l++)
        {
            double s = 0;
            for (ST_int q = 0; q < L; q++)
                s += m->XXinv[l + (size_t)q * L] * xbc[q];
            gx[l] = s;
        }
        double sres = uu;
        for (ST_int l = 0; l < L; l++)
            sres -= xbc[l] * gx[l];
        memcpy(Pi, m->Pi0, (size_t)L * K2 * sizeof(double));
        if (uu > 0.0 && sres > 1e-12 * uu)
        {
            for (ST_int j = 0; j < K2; j++)
            {
                double dj = tmp[1 + j];
                for (ST_int l = 0; l < L; l++)
                    dj -= xbc[l] * m->Pi0[l + (size_t)j * L];
                dj /= sres;
                for (ST_int l = 0; l < L; l++)
                    Pi[l + (size_t)j * L] -= gx[l] * dj;
            }
        }
    }

    /* The bootstrap columns [y*, Y2*, X1] = X ThF + Q CF */
    memset(d->ThF, 0, (size_t)L * F * sizeof(double));
    memset(d->CF, 0, (size_t)nb * E * sizeof(double));
    for (ST_int j = 0; j < K2; j++)
    {
        double *cj = d->CF + (size_t)(1 + j) * nb;
        cj[1 + j] = 1.0;
        for (ST_int l = 0; l < L; l++)
        {
            cj[iX + l] = -Pi[l + (size_t)j * L];
            d->ThF[l + (size_t)(1 + j) * L] = Pi[l + (size_t)j * L];
        }
    }
    for (ST_int q = 0; q < K1; q++)
        d->ThF[q + (size_t)(1 + K2 + q) * L] = 1.0;
    memcpy(c1, cu, (size_t)nb * sizeof(double));
    for (ST_int j = 0; j < K2; j++)
    {
        const double *cj = d->CF + (size_t)(1 + j) * nb;
        for (ST_int i = 0; i < nb; i++)
            c1[i] += cj[i] * bt[j];
    }
    memcpy(d->CF, c1, (size_t)nb * sizeof(double));
    for (ST_int l = 0; l < L; l++)
    {
        double s = l < K1 ? bt[K2 + l] : 0.0;
        for (ST_int j = 0; j < K2; j++)
            s += Pi[l + (size_t)j * L] * bt[j];
        d->ThF[l] = s;
    }
    for (ST_int f = 0; f < F; f++)
        for (ST_int l = 0; l < L; l++)
        {
            double s = 0;
            for (ST_int q = 0; q < L; q++)
                s += m->XX[l + (size_t)q * L] * d->ThF[q + (size_t)f * L];
            d->XXT[l + (size_t)f * L] = s;
        }
    for (ST_int b = 0; b < F; b++)
        for (ST_int a = 0; a < F; a++)
        {
            double s = 0;
            for (ST_int l = 0; l < L; l++)
                s += d->ThF[l + (size_t)a * L] * d->XXT[l + (size_t)b * L];
            d->G0[a + (size_t)b * F] = s;
        }
    d->ok = 1;
done:
    free(ws);
    return rc;
}

/* ------------------------------------------------------------------------
 * One replicate at one point
 * ------------------------------------------------------------------------ */

typedef struct
{
    double *M1, *XF, *T, *QC, *CQC, *FF, *W, *FPF, *H, *beta, *ak, *a, *mv, *th, *cv, *Aw, *Bw,
        *Wk, *u;
} wre_scratch;

static size_t wre_scratch_doubles(const civ_wre_model *m)
{
    const size_t L = m->L, F = m->F, E = m->E, nb = m->nb, K = m->K;
    return L * E + L * F + F * E + nb * E + E * E + F * F + L * F + F * F + K * K + K + K + F + L +
           L + nb + E * E * 4 + E * E * 2 + nb + 16;
}

static void wre_scratch_init(const civ_wre_model *m, double *p, wre_scratch *s)
{
    const size_t L = m->L, F = m->F, E = m->E, nb = m->nb, K = m->K;
    s->M1 = p;  p += L * E;
    s->XF = p;  p += L * F;
    s->T = p;   p += F * E;
    s->QC = p;  p += nb * E;
    s->CQC = p; p += E * E;
    s->FF = p;  p += F * F;
    s->W = p;   p += L * F;
    s->FPF = p; p += F * F;
    s->H = p;   p += K * K;
    s->beta = p; p += K;
    s->ak = p;  p += K;
    s->a = p;   p += F;
    s->mv = p;  p += L;
    s->th = p;  p += L;
    s->cv = p;  p += nb;
    s->Aw = p;  p += E * E * 2;
    s->Bw = p;  p += E * E * 2;
    s->Wk = p;  p += E * E * 2;
    s->u = p;
}

/*
 * The user's estimator on the bootstrap data of replicate r (cross
 * products XQ, QQ) under DGP d: the tested coefficient into *bk, and for
 * the unadjusted variance its value into *var; for the sandwich, m'X_g'W_g
 * u_g is then sum over the scores of mv, th, cv (see the callers).
 */
static int wre_core(const civ_wre_model *m, const wre_dgp *d, ST_int k, const double *XQ,
                    const double *QQ, wre_scratch *s, double *bk, double *var)
{
    const ST_int L = m->L, K = m->K, K2 = m->K2, K1 = m->K1, nb = m->nb, F = m->F, E = m->E;
    *bk = NAN;
    *var = NAN;
    /* M1 = XQ CF; XF = XXT + [M1 0] */
    for (ST_int e = 0; e < E; e++)
        for (ST_int l = 0; l < L; l++)
        {
            double t = 0;
            for (ST_int j = 0; j < nb; j++)
                t += XQ[l + (size_t)j * L] * d->CF[j + (size_t)e * nb];
            s->M1[l + (size_t)e * L] = t;
        }
    memcpy(s->XF, d->XXT, (size_t)L * F * sizeof(double));
    for (ST_int e = 0; e < E; e++)
        for (ST_int l = 0; l < L; l++)
            s->XF[l + (size_t)e * L] += s->M1[l + (size_t)e * L];
    /* FF = G0 + ThF'M1 + M1'ThF + CF'QQ CF */
    for (ST_int e = 0; e < E; e++)
        for (ST_int a = 0; a < F; a++)
        {
            double t = 0;
            for (ST_int l = 0; l < L; l++)
                t += d->ThF[l + (size_t)a * L] * s->M1[l + (size_t)e * L];
            s->T[a + (size_t)e * F] = t;
        }
    for (ST_int e = 0; e < E; e++)
        for (ST_int i = 0; i < nb; i++)
        {
            double t = 0;
            for (ST_int j = 0; j < nb; j++)
                t += QQ[i + (size_t)j * nb] * d->CF[j + (size_t)e * nb];
            s->QC[i + (size_t)e * nb] = t;
        }
    for (ST_int e = 0; e < E; e++)
        for (ST_int a = 0; a < E; a++)
        {
            double t = 0;
            for (ST_int i = 0; i < nb; i++)
                t += d->CF[i + (size_t)a * nb] * s->QC[i + (size_t)e * nb];
            s->CQC[a + (size_t)e * E] = t;
        }
    for (ST_int b = 0; b < F; b++)
        for (ST_int a = 0; a < F; a++)
        {
            double t = d->G0[a + (size_t)b * F];
            if (b < E)
                t += s->T[a + (size_t)b * F];
            if (a < E)
                t += s->T[b + (size_t)a * F];
            if (a < E && b < E)
                t += s->CQC[a + (size_t)b * E];
            s->FF[a + (size_t)b * F] = t;
        }
    /* W = XXinv XF; FPF = XF'W */
    for (ST_int f = 0; f < F; f++)
        for (ST_int l = 0; l < L; l++)
        {
            double t = 0;
            for (ST_int q = 0; q < L; q++)
                t += m->XXinv[l + (size_t)q * L] * s->XF[q + (size_t)f * L];
            s->W[l + (size_t)f * L] = t;
        }
    for (ST_int b = 0; b < F; b++)
        for (ST_int a = 0; a <= b; a++)
        {
            double t = 0;
            for (ST_int l = 0; l < L; l++)
                t += s->XF[l + (size_t)a * L] * s->W[l + (size_t)b * L];
            s->FPF[a + (size_t)b * F] = s->FPF[b + (size_t)a * F] = t;
        }

    /* kappa */
    double kappa = 1.0;
    if (m->est == 3)
        kappa = m->kclass;
    else if (m->est == 1 || m->est == 2)
    {
        double lambda = 1.0;
        if (m->L != m->K)
        {
            /* A = Y'M_X1 Y, Bm = Y'M_X Y for Y = [y*, Y2*] */
            for (ST_int b = 0; b < E; b++)
                for (ST_int a = 0; a < E; a++)
                {
                    s->Aw[a + (size_t)b * E] = s->FF[a + (size_t)b * F];
                    s->Bw[a + (size_t)b * E] = s->FF[a + (size_t)b * F] - s->FPF[a + (size_t)b * F];
                }
            if (K1 > 0)
            {
                /* less FF_ex XX11inv FF_xe (FF_xx = X1'WX1) */
                for (ST_int b = 0; b < E; b++)
                {
                    for (ST_int q = 0; q < K1; q++)
                    {
                        double t = 0;
                        for (ST_int p = 0; p < K1; p++)
                            t += m->XX11inv[q + (size_t)p * K1] * s->FF[(E + p) + (size_t)b * F];
                        s->u[q] = t;
                    }
                    for (ST_int a = 0; a < E; a++)
                    {
                        double t = 0;
                        for (ST_int q = 0; q < K1; q++)
                            t += s->FF[(E + q) + (size_t)a * F] * s->u[q];
                        s->Aw[a + (size_t)b * E] -= t;
                    }
                }
            }
            lambda = wre_min_geneig(E, s->Aw, s->Bw, s->Wk);
            if (!isfinite(lambda))
                return 0;
        }
        kappa = m->est == 2 ? lambda - m->fuller_adj : lambda;
    }

    /* H = (1-k) FF_ZZ + k FPF_ZZ; beta = H^-1 ((1-k) FF_Zy + k FPF_Zy) */
    for (ST_int b = 0; b < K; b++)
    {
        for (ST_int a = 0; a < K; a++)
            s->H[a + (size_t)b * K] = (1.0 - kappa) * s->FF[(1 + a) + (size_t)(1 + b) * F] +
                                      kappa * s->FPF[(1 + a) + (size_t)(1 + b) * F];
        s->beta[b] = (1.0 - kappa) * s->FF[1 + b] + kappa * s->FPF[1 + b];
    }
    if (wre_chol(K, s->H))
        return 0;
    wre_chol_solve(K, s->H, s->beta);
    for (ST_int r = 0; r < K; r++)
        s->ak[r] = r == k;
    wre_chol_solve(K, s->H, s->ak);
    for (ST_int r = 0; r < K; r++)
        if (!isfinite(s->beta[r]))
            return 0;
    *bk = s->beta[k];
    s->a[0] = 1.0;
    for (ST_int r = 0; r < K; r++)
        s->a[1 + r] = -s->beta[r];
    if (m->vce == 0)
    {
        double s2 = 0;
        for (ST_int b = 0; b < F; b++)
            for (ST_int a = 0; a < F; a++)
                s2 += s->a[a] * s->FF[a + (size_t)b * F] * s->a[b];
        *var = s->ak[k] * s2;
        return 1;
    }
    /* m = W_Z ak, th = ThF a, cv = CF a */
    for (ST_int l = 0; l < L; l++)
    {
        double t = 0;
        for (ST_int r = 0; r < K; r++)
            t += s->W[l + (size_t)(1 + r) * L] * s->ak[r];
        s->mv[l] = t;
        double u = 0;
        for (ST_int f = 0; f < F; f++)
            u += d->ThF[l + (size_t)f * L] * s->a[f];
        s->th[l] = u;
    }
    for (ST_int i = 0; i < nb; i++)
    {
        double t = 0;
        for (ST_int e = 0; e < E; e++)
            t += d->CF[i + (size_t)e * nb] * s->a[e];
        s->cv[i] = t;
    }
    (void)K2;
    return 1;
}

/* The sandwich's middle from per-cluster sums: sum over g of
 * (m'XX_g th + m'X_g'W_g Q_g cv)^2, with X_g'W_g Q_g from the groups'
 * sums and the weights (nested) or kept (projected) */
static double wre_meat_agg(const civ_wre_model *m, ST_int r, const wre_scratch *s,
                           const double *XQg)
{
    const ST_int L = m->L, nb = m->nb;
    double tot = 0;
    for (ST_int g = 0; g < m->G; g++)
    {
        const double *xx = m->XXg + (size_t)g * L * L;
        double sg = 0;
        for (ST_int b = 0; b < L; b++)
        {
            double t = 0;
            for (ST_int a = 0; a < L; a++)
                t += s->mv[a] * xx[a + (size_t)b * L];
            sg += t * s->th[b];
        }
        if (XQg)
        {
            const double *xq = XQg + (size_t)g * L * nb;
            for (ST_int j = 0; j < nb; j++)
            {
                double t = 0;
                for (ST_int a = 0; a < L; a++)
                    t += s->mv[a] * xq[a + (size_t)j * L];
                sg += t * s->cv[j];
            }
        }
        else
        {
            for (ST_int h = m->g_hs[g]; h < m->g_hs[g + 1]; h++)
            {
                const double *xb = m->XBh + (size_t)h * L * nb;
                double sh = 0;
                for (ST_int j = 0; j < nb; j++)
                {
                    double t = 0;
                    for (ST_int a = 0; a < L; a++)
                        t += s->mv[a] * xb[a + (size_t)j * L];
                    sh += t * s->cv[j];
                }
                sg += wre_v(m, r, m->h_c[h]) * sh;
            }
        }
        tot += sg * sg;
    }
    return tot;
}

/* The same from the rows: Q given column-major (q, stride N), row-major
 * (qr), or the basis times the weights (neither) */
static double wre_meat_rows(const civ_wre_model *m, ST_int r, const wre_scratch *s,
                            const double *q, const double *qr)
{
    const ST_int N = m->N, L = m->L, nb = m->nb, iX = m->iX;
    double tot = 0;
    for (ST_int g = 0; g < m->G; g++)
    {
        double sg = 0;
        for (ST_int h = m->g_hs[g]; h < m->g_hs[g + 1]; h++)
        {
            const double v = (q || qr) ? 1.0 : wre_v(m, r, m->h_c[h]);
            for (ST_int t = m->h_rs[h]; t < m->h_rs[h + 1]; t++)
            {
                const ST_int i = m->rows[t];
                const double *b = m->Br + (size_t)i * nb;
                double mx = 0, xt = 0, qc = 0;
                for (ST_int a = 0; a < L; a++)
                {
                    mx += s->mv[a] * b[iX + a];
                    xt += s->th[a] * b[iX + a];
                }
                if (q)
                    for (ST_int j = 0; j < nb; j++)
                        qc += q[(size_t)j * N + i] * s->cv[j];
                else if (qr)
                    for (ST_int j = 0; j < nb; j++)
                        qc += qr[(size_t)i * nb + j] * s->cv[j];
                else
                {
                    for (ST_int j = 0; j < nb; j++)
                        qc += b[j] * s->cv[j];
                    qc *= v;
                }
                sg += m->w[i] * mx * (xt + qc);
            }
        }
        tot += sg * sg;
    }
    return tot;
}

ST_retcode civ_wre_eval(civ_wre_model *m, ST_int j, ST_int npts, const double *b0, double *bk,
                        double *var, ST_int threads, civ_wre_project_fn project, void *ctx,
                        size_t budget)
{
    if (!m || !m->prepared || j < 0 || j >= m->n_test || npts < 1)
        return 198;
    const ST_int N = m->N, L = m->L, nb = m->nb, F = m->F, E = m->E, R = m->B + 1;
    const ST_int k = m->test[j];
    ST_retcode rc = 0;
    wre_dgp *d = calloc((size_t)npts, sizeof(wre_dgp));
    if (!d)
        return 920;
    const size_t dsz = (size_t)L * F * 2 + (size_t)nb * E + (size_t)F * F;
    double *dbuf = malloc((size_t)npts * dsz * sizeof(double));
    if (!dbuf)
    {
        free(d);
        return 920;
    }
    for (ST_int p = 0; p < npts && !rc; p++)
    {
        double *q = dbuf + (size_t)p * dsz;
        d[p].ThF = q;
        d[p].XXT = q + (size_t)L * F;
        d[p].CF = q + (size_t)L * F * 2;
        d[p].G0 = q + (size_t)L * F * 2 + (size_t)nb * E;
        rc = wre_dgp_build(m, k, b0[p], &d[p]);
    }
    for (size_t q = 0; q < (size_t)npts * R; q++)
        bk[q] = var[q] = NAN;
    if (rc)
    {
        free(dbuf);
        free(d);
        return rc;
    }
    const int P = wre_threads(threads);
    (void)P;
    const size_t sdbl = wre_scratch_doubles(m);
    const int reproject = m->projected && m->vce >= 1 && !m->cached;
    double *buf = NULL;
    ST_int batch = R;
    if (reproject)
    {
        const size_t col_bytes = (size_t)N * nb * sizeof(double);
        const size_t soft = budget < CIV_WRE_SOFT_BYTES ? budget : CIV_WRE_SOFT_BYTES;
        if (col_bytes > budget)
            rc = 909;
        batch = (ST_int)((soft / 4) / (col_bytes ? col_bytes : 1));
        if (batch < 1)
            batch = 1;
        if (batch > R)
            batch = R;
        if (!rc && !(buf = malloc((size_t)batch * col_bytes)))
            rc = 920;
    }
    int alloc_fail = 0;
    for (ST_int r0 = 0; r0 < R && !rc; r0 += batch)
    {
        const ST_int nr = (r0 + batch <= R) ? batch : R - r0;
        if (reproject && (rc = wre_project_batch(m, r0, nr, buf, project, ctx)))
            break;
#ifdef _OPENMP
#pragma omp parallel num_threads(P) if (P > 1)
#endif
        {
            double *sp = malloc(sdbl * sizeof(double));
            wre_scratch s;
            if (!sp)
            {
#ifdef _OPENMP
#pragma omp atomic write
#endif
                alloc_fail = 1;
            }
            else
                wre_scratch_init(m, sp, &s);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 2)
#endif
            for (ST_int rr = 0; rr < nr; rr++)
            {
                if (!sp)
                    continue;
                const ST_int r = r0 + rr;
                const double *XQ = m->XQ + (size_t)r * L * nb, *QQ = m->QQ + (size_t)r * nb * nb;
                for (ST_int p = 0; p < npts; p++)
                {
                    if (!d[p].ok)
                        continue;
                    double b = NAN, v = NAN;
                    if (!wre_core(m, &d[p], k, XQ, QQ, &s, &b, &v))
                        continue;
                    if (m->vce >= 1)
                    {
                        if (reproject)
                            v = wre_meat_rows(m, r, &s, buf + (size_t)rr * nb * N, NULL);
                        else if (m->score_agg)
                            v = wre_meat_agg(m, r, &s,
                                             m->XQg ? m->XQg + (size_t)r * m->G * L * nb : NULL);
                        else
                            v = wre_meat_rows(m, r, &s, NULL,
                                              m->Qr ? m->Qr + (size_t)r * N * nb : NULL);
                    }
                    bk[(size_t)p * R + r] = b;
                    var[(size_t)p * R + r] = v;
                }
            }
            free(sp);
        }
        if (alloc_fail)
            rc = 920;
    }
    free(buf);
    free(dbuf);
    free(d);
    return rc;
}

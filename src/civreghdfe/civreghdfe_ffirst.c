/*
    civreghdfe_ffirst.c
    ivreg2's first-stage and weak-instrument statistics (its s_ffirst),
    computed from the partialled data the estimation uses.

    With y the dependent variable, X1 the K1 endogenous regressors, X2 the
    K2 included and Z1 the L1 excluded instruments (Z = [X2 Z1], L columns),
    N, d = dofminus(), s = sdofminus:
      - the system of reduced form and first stages [y X1] = Z B + E, its
        classical VCE (E'WE/(N-d)) (x) (Z'WZ)^-1 and its VCE
        (I (x) (Z'WZ)^-1) S (I (x) (Z'WZ)^-1) with S the moment covariance of
        the model's VCE (civ_omega_build) of the residuals E, times N/(N-d)
        without clusters;
      - per endogenous regressor: rmse E'WE/(N-L-d-s), Shea's partial R2
        [(X'WX)^-1]_ii / [(X'P_Z X)^-1]_ii, the partial R2 W/(N-d)/(1+W/(N-d))
        of the classical Wald statistic W of the excluded instruments, and
        F = W/(N-d) (N-d-L-s)/L1 (clusters: W/(N-1) (N-L-s) (G-1)/G/L1),
        the Wald statistic robust when S is not iid;
      - Angrist-Pischke and Sanderson-Windmeijer statistics (K1 > 1): the
        regressor's residual from its regression on the other fitted (AP)
        or actual (SW) regressors with the 2SLS coefficients of that
        regression, then the Wald test of the excluded instruments in the
        regression of that residual on Z, with L1-K1+1 numerator dof;
      - Anderson-Rubin: the Wald test of the excluded instruments in the
        reduced form; Stock-Wright S: N g' Omega^-1 g, g the moments of y
        residualized on X2.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "civreghdfe_ffirst.h"
#include "../ctools_matrix.h"
#include "../ctools_spi.h"
#include "../ctools_config.h"
#include "../ctools_alloc.h"
#include "../ctools_ols.h"

/* Stata's Ftail and chi2tail are not available to plugins; the ado
   computes the p-values from the stored statistics and dof. */

/* q = b' A^-1 b for symmetric A (n x n) by Cholesky, as ivreg2's
   cholsolve(): NaN (missing) when A is not positive definite */
static ST_double civ_quad(const ST_double *A, ST_int n, const ST_double *b)
{
    if (n <= 0) return NAN;
    ST_double *L = (ST_double *)malloc((size_t)n * n * sizeof(ST_double));
    ST_double *x = (ST_double *)malloc((size_t)n * sizeof(ST_double));
    ST_double q = NAN;
    if (L && x) {
        memcpy(L, A, (size_t)n * n * sizeof(ST_double));
        if (ctools_cholesky(L, n) == 0) {
            /* L is lower triangular, row-major in the ctools convention:
               solve L u = b, then q = u'u */
            for (ST_int i = 0; i < n; i++) {
                ST_double s = b[i];
                for (ST_int j = 0; j < i; j++) s -= L[i * n + j] * x[j];
                x[i] = s / L[i * n + i];
            }
            q = 0.0;
            for (ST_int i = 0; i < n; i++) q += x[i] * x[i];
        }
    }
    free(L); free(x);
    return isfinite(q) ? q : NAN;
}

/* M[ex, ex] block of f * Zi S Zi (all L x L), ex = rows K2..L-1; out is L1 x L1 */
static void civ_excl_sandwich(const ST_double *Zi, const ST_double *S, ST_int L, ST_int K2,
                              ST_double f, ST_double *T, ST_double *out)
{
    ST_int L1 = L - K2;
    /* T = S Zi[:, ex] (L x L1) */
    for (ST_int b = 0; b < L1; b++) {
        const ST_double *zb = Zi + (size_t)(K2 + b) * L;
        for (ST_int k = 0; k < L; k++) {
            ST_double s = 0.0;
            for (ST_int j = 0; j < L; j++) s += S[(size_t)j * L + k] * zb[j];
            T[(size_t)b * L + k] = s;
        }
    }
    for (ST_int b = 0; b < L1; b++)
        for (ST_int a = 0; a < L1; a++) {
            const ST_double *za = Zi + (size_t)(K2 + a) * L;
            ST_double s = 0.0;
            for (ST_int k = 0; k < L; k++) s += za[k] * T[(size_t)b * L + k];
            out[(size_t)b * L1 + a] = f * s;
        }
}

/* residual r = v - Z (Z'WZ)^-1 Z'W v; coef (L) returned */
static void civ_resid_on_z(const ST_double *v, const ST_double *Z, ST_int L, const ST_double *w,
                           ST_int N, const ST_double *Zi, ST_double *zv, ST_double *coef,
                           ST_double *r)
{
    ctools_matmul_atdb(Z, v, w, N, L, 1, zv);
    for (ST_int i = 0; i < L; i++) {
        ST_double s = 0.0;
        for (ST_int j = 0; j < L; j++) s += Zi[(size_t)j * L + i] * zv[j];
        coef[i] = s;
    }
    memcpy(r, v, (size_t)N * sizeof(ST_double));
    for (ST_int j = 0; j < L; j++) {
        const ST_double *zj = Z + (size_t)j * N;
        ST_double c = coef[j];
        if (c == 0.0) continue;
        for (ST_int i = 0; i < N; i++) r[i] -= zj[i] * c;
    }
}

static ST_double civ_wss1(const ST_double *e, const ST_double *w, ST_int N)
{
    ST_double s = 0.0;
    for (ST_int i = 0; i < N; i++) s += (w ? w[i] : 1.0) * e[i] * e[i];
    return s;
}

/* ivreg2's F scaling of a Wald statistic with df1 numerator dof */
static void civ_fstat(const civ_ffirst_dof *dof, ST_double wald, ST_int df1,
                      ST_double *F, ST_double *df_r)
{
    ST_double N = dof->N;
    if (!dof->cluster) {
        *df_r = N - dof->dofminus - dof->L_dof - dof->sdofminus;
        *F = wald / (N - dof->dofminus) * (*df_r) / (ST_double)df1;
    } else {
        ST_double G = (ST_double)dof->nclust;
        *df_r = G - 1.0;
        *F = wald / (N - 1.0) * (N - dof->L_dof - dof->sdofminus) * (G - 1.0) / G / (ST_double)df1;
    }
    if (!isfinite(*F)) *F = SV_missval;
}

ST_retcode civ_ffirst(const ST_double *y, const ST_double *X1, ST_int K1,
                      const ST_double *Z, ST_int L, ST_int K2,
                      const ST_double *w, ST_int N, const civ_omega *om,
                      const civ_ffirst_dof *dof)
{
    ST_int L1 = L - K2, Q = K1 + 1, D = Q * L, K = K1 + K2;
    ST_double Ne = dof->N, Nd = dof->N - dof->dofminus;
    /* S carries 1/(N-d) without clusters (m_omega's dofminus) */
    ST_double f = dof->cluster ? 1.0 : Ne / Nd;
    ST_retcode rc = 0;
    if (K1 <= 0 || L1 <= 0) return 0;

    ST_double *Yall = NULL, *E = NULL, *ZZ = NULL, *Zi = NULL, *ZY = NULL, *B = NULL;
    ST_double *ee = NULL, *S = NULL, *Vr = NULL, *T = NULL, *blk = NULL, *zv = NULL, *coef = NULL;
    ST_double *XX = NULL, *XXi = NULL, *ZX = NULL, *XPZX = NULL, *XPZXi = NULL;
    ST_double *Xact = NULL, *Xhat = NULL, *Hh = NULL, *Hx = NULL, *e1 = NULL, *e2 = NULL;
    ST_double *firstmat = NULL, *Rb = NULL, *S1 = NULL;

    Yall = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)Q, sizeof(ST_double));
    E = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)Q, sizeof(ST_double));
    ZZ = (ST_double *)calloc((size_t)L * L, sizeof(ST_double));
    Zi = (ST_double *)calloc((size_t)L * L, sizeof(ST_double));
    ZY = (ST_double *)calloc((size_t)L * Q, sizeof(ST_double));
    B = (ST_double *)calloc((size_t)L * Q, sizeof(ST_double));
    ee = (ST_double *)calloc((size_t)Q * Q, sizeof(ST_double));
    S = (ST_double *)ctools_safe_calloc3((size_t)D, (size_t)D, sizeof(ST_double));
    Vr = (ST_double *)ctools_safe_calloc3((size_t)D, (size_t)D, sizeof(ST_double));
    T = (ST_double *)calloc((size_t)L * L, sizeof(ST_double));
    blk = (ST_double *)calloc((size_t)L * L, sizeof(ST_double));
    zv = (ST_double *)calloc((size_t)L, sizeof(ST_double));
    coef = (ST_double *)calloc((size_t)(L > K ? L : K), sizeof(ST_double));
    Rb = (ST_double *)calloc((size_t)L, sizeof(ST_double));
    S1 = (ST_double *)calloc((size_t)L * L, sizeof(ST_double));
    firstmat = (ST_double *)calloc((size_t)CIV_FFIRST_ROWS * K1, sizeof(ST_double));
    e1 = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double));
    e2 = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double));
    if (!Yall || !E || !ZZ || !Zi || !ZY || !B || !ee || !S || !Vr || !T || !blk || !zv ||
        !coef || !Rb || !S1 || !firstmat || !e1 || !e2) { rc = 920; goto done; }

    memcpy(Yall, y, (size_t)N * sizeof(ST_double));
    memcpy(Yall + N, X1, (size_t)N * K1 * sizeof(ST_double));

    ctools_matmul_atdb(Z, Z, w, N, L, L, ZZ);
    if (civ_sym_inverse(ZZ, L, Zi) != 0) { rc = 920; goto done; }
    ctools_matmul_atdb(Z, Yall, w, N, L, Q, ZY);
    ctools_matmul_ab(Zi, ZY, L, L, Q, B);
    /* E = [y X1] - Z B */
    ctools_matmul_xb(Z, N, L, B, Q, E);
    for (size_t j = 0; j < (size_t)N * Q; j++) E[j] = Yall[j] - E[j];
    ctools_matmul_atdb(E, E, w, N, Q, Q, ee);

    /* VCE of the system: Vr = f (I (x) Zi) S (I (x) Zi) */
    if (civ_omega_build(om, E, Q, Z, L, S) != 0) { rc = 920; goto done; }
    for (ST_int c2 = 0; c2 < Q; c2++)
        for (ST_int c1 = 0; c1 < Q; c1++) {
            /* blk = S_{c1,c2} (L x L), then Zi blk Zi */
            for (ST_int b = 0; b < L; b++)
                for (ST_int a = 0; a < L; a++)
                    blk[(size_t)b * L + a] = S[(size_t)(c2 * L + b) * D + (c1 * L + a)];
            ctools_matmul_ab(Zi, blk, L, L, L, T);
            ctools_matmul_ab(T, Zi, L, L, L, blk);
            for (ST_int b = 0; b < L; b++)
                for (ST_int a = 0; a < L; a++)
                    Vr[(size_t)(c2 * L + b) * D + (c1 * L + a)] = f * blk[(size_t)b * L + a];
        }
    /* symmetrize (ivreg2's makesymmetric) */
    for (ST_int c = 0; c < D; c++)
        for (ST_int r = c + 1; r < D; r++) {
            ST_double m = Vr[(size_t)c * D + r];
            Vr[(size_t)r * D + c] = m;
        }

    /* Anderson-Rubin: excluded instruments in the reduced form */
    {
        for (ST_int a = 0; a < L1; a++) Rb[a] = B[K2 + a];
        for (ST_int b = 0; b < L1; b++)
            for (ST_int a = 0; a < L1; a++)
                T[(size_t)b * L1 + a] = Vr[(size_t)(K2 + b) * D + (K2 + a)];
        ST_double arw = civ_quad(T, L1, Rb);
        ST_double arf, ardf2;
        if (!dof->cluster) {
            ardf2 = Ne - dof->dofminus - dof->L_dof - dof->sdofminus;
            arf = arw / Nd * ardf2 / (ST_double)L1;
        } else {
            ardf2 = (ST_double)dof->nclust - 1.0;
            arf = arw / (Ne - 1.0) * (Ne - dof->L_dof - dof->sdofminus) / (ST_double)dof->nclust * ardf2 / (ST_double)L1;
        }
        if ((rc = ctools_scal_save("__civreghdfe_archi2", isfinite(arw) ? arw : SV_missval)) ||
            (rc = ctools_scal_save("__civreghdfe_arf", isfinite(arf) ? arf : SV_missval)) ||
            (rc = ctools_scal_save("__civreghdfe_ardf", (ST_double)L1)) ||
            (rc = ctools_scal_save("__civreghdfe_ardf_r", ardf2)))
            goto done;
    }

    /* Stock-Wright S: y residualized on X2, then N g' Omega^-1 g */
    {
        if (K2 > 0) {
            ST_double *XtX = (ST_double *)calloc((size_t)K2 * K2, sizeof(ST_double));
            ST_double *Xi = (ST_double *)calloc((size_t)K2 * K2, sizeof(ST_double));
            ST_double *Xy = (ST_double *)calloc((size_t)K2, sizeof(ST_double));
            if (!XtX || !Xi || !Xy) { free(XtX); free(Xi); free(Xy); rc = 920; goto done; }
            ctools_matmul_atdb(Z, Z, w, N, K2, K2, XtX);
            ctools_matmul_atdb(Z, y, w, N, K2, 1, Xy);
            ST_int irc = civ_sym_inverse(XtX, K2, Xi);
            memcpy(e1, y, (size_t)N * sizeof(ST_double));
            if (!irc) {
                for (ST_int j = 0; j < K2; j++) {
                    ST_double c = 0.0;
                    for (ST_int k = 0; k < K2; k++) c += Xi[(size_t)k * K2 + j] * Xy[k];
                    const ST_double *zj = Z + (size_t)j * N;
                    for (ST_int i = 0; i < N; i++) e1[i] -= zj[i] * c;
                }
            }
            free(XtX); free(Xi); free(Xy);
            if (irc) { rc = 920; goto done; }
        } else {
            memcpy(e1, y, (size_t)N * sizeof(ST_double));
        }
        if (civ_omega_build(om, e1, 1, Z, L, S1) != 0) { rc = 920; goto done; }
        for (size_t j = 0; j < (size_t)L * L; j++) S1[j] *= f;
        ctools_matmul_atdb(Z, e1, w, N, L, 1, zv);
        ST_double ss = civ_quad(S1, L, zv);
        if ((rc = ctools_scal_save("__civreghdfe_sstat", isfinite(ss) ? ss : SV_missval)) ||
            (rc = ctools_scal_save("__civreghdfe_sstatdf", (ST_double)L1)))
            goto done;
    }

    /* Shea's partial R2: X = [X1 X2] */
    XX = (ST_double *)calloc((size_t)K * K, sizeof(ST_double));
    XXi = (ST_double *)calloc((size_t)K * K, sizeof(ST_double));
    ZX = (ST_double *)calloc((size_t)L * K, sizeof(ST_double));
    XPZX = (ST_double *)calloc((size_t)K * K, sizeof(ST_double));
    XPZXi = (ST_double *)calloc((size_t)K * K, sizeof(ST_double));
    Xact = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K, sizeof(ST_double));
    Xhat = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K, sizeof(ST_double));
    Hh = (ST_double *)calloc((size_t)K * K, sizeof(ST_double));
    Hx = (ST_double *)calloc((size_t)K * K1, sizeof(ST_double));
    if (!XX || !XXi || !ZX || !XPZX || !XPZXi || !Xact || !Xhat || !Hh || !Hx) { rc = 920; goto done; }
    memcpy(Xact, X1, (size_t)N * K1 * sizeof(ST_double));
    if (K2 > 0) memcpy(Xact + (size_t)N * K1, Z, (size_t)N * K2 * sizeof(ST_double));
    /* fitted X1 = Z B[:, 1..K1], then X2 */
    ctools_matmul_xb(Z, N, L, B + L, K1, Xhat);
    if (K2 > 0) memcpy(Xhat + (size_t)N * K1, Z, (size_t)N * K2 * sizeof(ST_double));
    ctools_matmul_atdb(Xact, Xact, w, N, K, K, XX);
    ctools_matmul_atdb(Z, Xact, w, N, L, K, ZX);
    {
        ST_double *t = (ST_double *)calloc((size_t)L * K, sizeof(ST_double));
        if (!t) { rc = 920; goto done; }
        ctools_matmul_ab(Zi, ZX, L, L, K, t);
        for (ST_int c = 0; c < K; c++)
            for (ST_int r = 0; r < K; r++) {
                ST_double s = 0.0;
                for (ST_int k = 0; k < L; k++) s += ZX[(size_t)r * L + k] * t[(size_t)c * L + k];
                XPZX[(size_t)c * K + r] = s;
            }
        free(t);
    }
    if (civ_sym_inverse(XX, K, XXi) != 0 || civ_sym_inverse(XPZX, K, XPZXi) != 0) { rc = 920; goto done; }
    ctools_matmul_atdb(Xhat, Xhat, w, N, K, K, Hh);
    ctools_matmul_atdb(Xhat, X1, w, N, K, K1, Hx);

    for (ST_int i = 0; i < K1; i++) {
        ST_double *fm = firstmat + (size_t)i * CIV_FFIRST_ROWS;
        ST_double eii = ee[(size_t)(i + 1) * Q + (i + 1)];
        fm[0] = sqrt(eii / (Ne - dof->L_dof - dof->dofminus - dof->sdofminus));
        fm[1] = XXi[(size_t)i * K + i] / XPZXi[(size_t)i * K + i];
        const ST_double *bi = B + (size_t)(i + 1) * L;
        for (ST_int a = 0; a < L1; a++) Rb[a] = bi[K2 + a];
        /* classical Wald: (e'We/(N-d)) Zi[ex,ex] */
        ST_double sig = eii / Nd;
        for (ST_int b = 0; b < L1; b++)
            for (ST_int a = 0; a < L1; a++)
                T[(size_t)b * L1 + a] = sig * Zi[(size_t)(K2 + b) * L + (K2 + a)];
        ST_double wald = civ_quad(T, L1, Rb);
        ST_double pr2 = (wald / Nd) / (1.0 + wald / Nd);
        if (dof->robust) {
            ST_int off = (i + 1) * L;
            for (ST_int b = 0; b < L1; b++)
                for (ST_int a = 0; a < L1; a++)
                    T[(size_t)b * L1 + a] = Vr[(size_t)(off + K2 + b) * D + (off + K2 + a)];
            wald = civ_quad(T, L1, Rb);
        }
        ST_double F, df_r;
        civ_fstat(dof, wald, L1, &F, &df_r);
        fm[2] = pr2; fm[3] = F; fm[4] = L1; fm[5] = df_r; fm[6] = SV_missval;  /* p by the ado */

        ST_double swF, swdf1, swdf2, swchi2, swr2, apF, apchi2, apr2;
        if (K1 == 1) {
            swF = apF = F; swdf1 = L1; swdf2 = df_r;
            swchi2 = apchi2 = wald; swr2 = apr2 = pr2;
        } else {
            /* the other regressors, fitted (AP) and actual (SW) */
            ST_int Km = K - 1;
            ST_double *Hm = (ST_double *)calloc((size_t)Km * Km, sizeof(ST_double));
            ST_double *Hmi = (ST_double *)calloc((size_t)Km * Km, sizeof(ST_double));
            ST_double *hx = (ST_double *)calloc((size_t)Km, sizeof(ST_double));
            if (!Hm || !Hmi || !hx) { free(Hm); free(Hmi); free(hx); rc = 920; goto done; }
            for (ST_int c = 0, cc = 0; c < K; c++) {
                if (c == i) continue;
                hx[cc] = Hx[(size_t)i * K + c];
                for (ST_int r = 0, rr = 0; r < K; r++) {
                    if (r == i) continue;
                    Hm[(size_t)cc * Km + rr] = Hh[(size_t)c * K + r];
                    rr++;
                }
                cc++;
            }
            ST_int hrc = civ_sym_inverse(Hm, Km, Hmi);
            ST_double *b1 = coef;  /* Km <= L */
            for (ST_int r = 0; r < Km && !hrc; r++) {
                ST_double s = 0.0;
                for (ST_int c = 0; c < Km; c++) s += Hmi[(size_t)c * Km + r] * hx[c];
                b1[r] = s;
            }
            free(Hm); free(Hmi); free(hx);
            if (hrc) { rc = 920; goto done; }
            ST_double stat[2][3];  /* [AP, SW] x [Wald, r2, robust Wald] */
            for (int which = 0; which < 2; which++) {
                const ST_double *Xs = which == 0 ? Xhat : Xact;
                const ST_double *xi = X1 + (size_t)i * N;
                memcpy(e1, xi, (size_t)N * sizeof(ST_double));
                for (ST_int c = 0, cc = 0; c < K; c++) {
                    if (c == i) continue;
                    const ST_double *xc = Xs + (size_t)c * N;
                    ST_double bc = b1[cc++];
                    for (ST_int r = 0; r < N; r++) e1[r] -= xc[r] * bc;
                }
                ST_double *b2 = Rb;  /* reuse: L entries */
                ST_double *b2full = (ST_double *)calloc((size_t)L, sizeof(ST_double));
                if (!b2full) { rc = 920; goto done; }
                civ_resid_on_z(e1, Z, L, w, N, Zi, zv, b2full, e2);
                for (ST_int a = 0; a < L1; a++) b2[a] = b2full[K2 + a];
                free(b2full);
                ST_double s2 = civ_wss1(e2, w, N) / Nd;
                for (ST_int b = 0; b < L1; b++)
                    for (ST_int a = 0; a < L1; a++)
                        T[(size_t)b * L1 + a] = s2 * Zi[(size_t)(K2 + b) * L + (K2 + a)];
                ST_double wd = civ_quad(T, L1, b2);
                stat[which][0] = wd;
                stat[which][1] = (wd / Nd) / (1.0 + wd / Nd);
                stat[which][2] = wd;
                if (dof->robust) {
                    if (civ_omega_build(om, e2, 1, Z, L, S1) != 0) { rc = 920; goto done; }
                    ST_double *Tt = (ST_double *)calloc((size_t)L * L1, sizeof(ST_double));
                    ST_double *Vb = (ST_double *)calloc((size_t)L1 * L1, sizeof(ST_double));
                    if (!Tt || !Vb) { free(Tt); free(Vb); rc = 920; goto done; }
                    civ_excl_sandwich(Zi, S1, L, K2, f, Tt, Vb);
                    stat[which][2] = civ_quad(Vb, L1, b2);
                    free(Tt); free(Vb);
                }
            }
            ST_int fdf1 = L1 - K1 + 1;
            swdf1 = fdf1;
            civ_fstat(dof, stat[0][2], fdf1, &apF, &swdf2);
            civ_fstat(dof, stat[1][2], fdf1, &swF, &swdf2);
            apchi2 = stat[0][2]; apr2 = stat[0][1];
            swchi2 = stat[1][2]; swr2 = stat[1][1];
        }
        fm[7] = swF; fm[8] = swdf1; fm[9] = swdf2; fm[10] = SV_missval;
        fm[11] = swchi2; fm[12] = SV_missval; fm[13] = swr2;
        fm[14] = apF; fm[15] = swdf1; fm[16] = swdf2; fm[17] = SV_missval;
        fm[18] = apchi2; fm[19] = SV_missval; fm[20] = apr2;
    }

    /* Store: e(first) without the p-values (the ado fills them in) */
    for (ST_int i = 0; i < K1 && !rc; i++)
        for (ST_int r = 0; r < CIV_FFIRST_ROWS && !rc; r++) {
            ST_double v = firstmat[(size_t)i * CIV_FFIRST_ROWS + r];
            rc = ctools_mat_store("__civreghdfe_firstmat", r + 1, i + 1, isfinite(v) ? v : SV_missval);
        }
    if (rc) goto done;

    /* The system in ivreg2's order: equations [y X1], instruments
       [excluded, included]; S in ivreg2's scale (divided by N) */
    {
        ST_int *perm = (ST_int *)malloc((size_t)L * sizeof(ST_int));
        if (!perm) { rc = 920; goto done; }
        for (ST_int r = 0; r < L; r++) perm[r] = r < L1 ? K2 + r : r - L1;
        for (ST_int c = 0; c < Q && !rc; c++)
            for (ST_int r = 0; r < L && !rc; r++)
                rc = ctools_mat_store("__civreghdfe_sfb", 1, c * L + r + 1, B[(size_t)c * L + perm[r]]);
        for (ST_int c2 = 0; c2 < Q && !rc; c2++)
            for (ST_int r2 = 0; r2 < L && !rc; r2++)
                for (ST_int c1 = 0; c1 < Q && !rc; c1++)
                    for (ST_int r1 = 0; r1 < L && !rc; r1++) {
                        size_t src = (size_t)(c2 * L + perm[r2]) * D + (c1 * L + perm[r1]);
                        ST_int row = c1 * L + r1 + 1, col = c2 * L + r2 + 1;
                        rc = ctools_mat_store("__civreghdfe_sfV", row, col, Vr[src]);
                        if (!rc) rc = ctools_mat_store("__civreghdfe_sfS", row, col, S[src] * f / Ne);
                    }
        free(perm);
    }
    if (!rc) rc = ctools_scal_save("__civreghdfe_rmse_rf",
                                   sqrt(ee[0] / (Ne - dof->L_dof - dof->dofminus - dof->sdofminus)));
    if (!rc) rc = ctools_scal_save("__civreghdfe_ff_K1", (ST_double)K1);
    if (!rc) rc = ctools_scal_save("__civreghdfe_ff_L", (ST_double)L);

done:
    free(Yall); free(E); free(ZZ); free(Zi); free(ZY); free(B); free(ee); free(S); free(Vr);
    free(T); free(blk); free(zv); free(coef); free(Rb); free(S1); free(firstmat); free(e1); free(e2);
    free(XX); free(XXi); free(ZX); free(XPZX); free(XPZXi); free(Xact); free(Xhat); free(Hh); free(Hx);
    return rc;
}

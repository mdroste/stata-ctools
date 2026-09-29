/*
    civreghdfe_tests.c
    IV Diagnostic Tests Implementation

    Implements:
    - First-stage F statistics
    - Underidentification tests (Anderson LM, Kleibergen-Paap rk LM)
    - Weak instrument tests (Cragg-Donald F, Kleibergen-Paap rk Wald F)
    - Overidentification tests (Sargan, Hansen J)
    - Endogeneity tests (Durbin-Wu-Hausman)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "civreghdfe_tests.h"
#include "civreghdfe_matrix.h"
#include "civreghdfe_vce.h"
#include "../ctools_config.h"

/* Shared OLS functions */
#include "../ctools_ols.h"

/* ------------------------------------------------------------------------
   Linear-algebra helpers (column-major, small matrices)
   ------------------------------------------------------------------------ */

/* g = Z'We for the L columns of Z */
static void civ_zwe(const ST_double *Z, ST_int L, const ST_double *e,
                    const ST_double *w, ST_int N, ST_double *g)
{
    for (ST_int l = 0; l < L; l++) {
        const ST_double *z = Z + (size_t)l * N;
        ST_double s = 0.0;
        if (w) { for (ST_int i = 0; i < N; i++) s += z[i] * w[i] * e[i]; }
        else   { for (ST_int i = 0; i < N; i++) s += z[i] * e[i]; }
        g[l] = s;
    }
}

/* A'WB (ka x kb), W = diag(w) or I */
static void civ_cross(const ST_double *A, ST_int ka, const ST_double *B, ST_int kb,
                      const ST_double *w, ST_int N, ST_double *out)
{
    if (w) ctools_matmul_atdb(A, B, w, N, ka, kb, out);
    else   ctools_matmul_atb(A, B, N, ka, kb, out);
}

/* Residuals of the ncol columns of A on the P columns of X by weighted least
   squares (ivreg2's and ranktest's partialling), written to out */
static ST_int civ_partial_cols(const ST_double *A, ST_int ncol, const ST_double *X, ST_int P,
                               const ST_double *w, ST_int N, ST_double *out)
{
    if (out != A) memcpy(out, A, (size_t)N * ncol * sizeof(ST_double));
    if (P <= 0 || ncol <= 0) return 0;
    ST_double *XX = (ST_double *)calloc((size_t)P * P, sizeof(ST_double));
    ST_double *XA = (ST_double *)calloc((size_t)P * ncol, sizeof(ST_double));
    ST_int rc = (XX && XA) ? 0 : -1;
    if (!rc) {
        civ_cross(X, P, A, ncol, w, N, XA);
        civ_cross(X, P, X, P, w, N, XX);
        rc = civ_sym_solve(XX, P, XA, ncol, XA);
    }
    if (!rc) {
        for (ST_int c = 0; c < ncol; c++) {
            ST_double *o = out + (size_t)c * N;
            const ST_double *b = XA + (size_t)c * P;
            for (ST_int p = 0; p < P; p++) {
                const ST_double *x = X + (size_t)p * N;
                ST_double bp = b[p];
                if (bp != 0.0) for (ST_int i = 0; i < N; i++) o[i] -= x[i] * bp;
            }
        }
    }
    free(XX); free(XA);
    return rc;
}

/* ranktest divides the partialled variables by their standard deviations
   (unweighted, n - 1; 1 when zero). The rank statistics do not depend on it,
   except through invsym's generalized inverse when their S is singular
   (e.g. too few clusters), so it is replicated. */
static void civ_standardize_cols(ST_double *A, ST_int ncol, ST_int N)
{
    for (ST_int c = 0; c < ncol; c++) {
        ST_double *a = A + (size_t)c * N;
        ST_double mean = 0.0, ss = 0.0;
        for (ST_int i = 0; i < N; i++) mean += a[i];
        mean /= (ST_double)N;
        for (ST_int i = 0; i < N; i++) { ST_double d = a[i] - mean; ss += d * d; }
        ST_double sd = (N > 1) ? sqrt(ss / (ST_double)(N - 1)) : 0.0;
        if (!(sd > 0.0)) sd = 1.0;
        for (ST_int i = 0; i < N; i++) a[i] /= sd;
    }
}

/* Lower Cholesky factor of the SPD n x n A (column-major): A = L L' */
static ST_int civ_chol_lower(const ST_double *A, ST_int n, ST_double *Lf)
{
    memset(Lf, 0, (size_t)n * n * sizeof(ST_double));
    for (ST_int j = 0; j < n; j++) {
        ST_double s = A[(size_t)j * n + j];
        for (ST_int k = 0; k < j; k++) s -= Lf[(size_t)k * n + j] * Lf[(size_t)k * n + j];
        if (!(s > 0.0)) return -1;
        ST_double d = sqrt(s);
        Lf[(size_t)j * n + j] = d;
        for (ST_int i = j + 1; i < n; i++) {
            ST_double t = A[(size_t)j * n + i];
            for (ST_int k = 0; k < j; k++) t -= Lf[(size_t)k * n + i] * Lf[(size_t)k * n + j];
            Lf[(size_t)j * n + i] = t / d;
        }
    }
    return 0;
}

/* Inverse of the lower-triangular n x n Lf (column-major) */
static void civ_lower_inverse(const ST_double *Lf, ST_int n, ST_double *Li)
{
    memset(Li, 0, (size_t)n * n * sizeof(ST_double));
    for (ST_int j = 0; j < n; j++) {
        Li[(size_t)j * n + j] = 1.0 / Lf[(size_t)j * n + j];
        for (ST_int i = j + 1; i < n; i++) {
            ST_double s = 0.0;
            for (ST_int k = j; k < i; k++) s -= Lf[(size_t)k * n + i] * Li[(size_t)j * n + k];
            Li[(size_t)j * n + i] = s / Lf[(size_t)i * n + i];
        }
    }
}

/* Eigen-decomposition of the symmetric n x n A (column-major) by cyclic
   Jacobi sweeps: eigenvalues ascending in ev, eigenvectors in the columns of
   Q (optional). Returns 0, or -1 on allocation failure. */
static ST_int civ_sym_eigen(const ST_double *A, ST_int n, ST_double *ev, ST_double *Q)
{
    ST_double *D = (ST_double *)malloc((size_t)n * n * sizeof(ST_double));
    ST_double *V = (ST_double *)calloc((size_t)n * n, sizeof(ST_double));
    ST_int *ord = (ST_int *)malloc((size_t)(n > 0 ? n : 1) * sizeof(ST_int));
    if (!D || !V || !ord) { free(D); free(V); free(ord); return -1; }
    memcpy(D, A, (size_t)n * n * sizeof(ST_double));
    for (ST_int i = 0; i < n; i++) V[(size_t)i * n + i] = 1.0;

    for (ST_int sweep = 0; sweep < 100; sweep++) {
        ST_double off = 0.0, diag = 0.0;
        for (ST_int j = 0; j < n; j++)
            for (ST_int i = 0; i < n; i++) {
                ST_double a = D[(size_t)j * n + i];
                if (i == j) diag += a * a; else off += a * a;
            }
        if (off <= 1e-30 * diag || off == 0.0) break;
        for (ST_int p = 0; p < n - 1; p++) {
            for (ST_int q = p + 1; q < n; q++) {
                ST_double apq = D[(size_t)q * n + p];
                if (apq == 0.0) continue;
                ST_double app = D[(size_t)p * n + p], aqq = D[(size_t)q * n + q];
                ST_double theta = (aqq - app) / (2.0 * apq);
                ST_double t = (theta >= 0.0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
                ST_double c = 1.0 / sqrt(t * t + 1.0), s = t * c;
                for (ST_int k = 0; k < n; k++) {   /* columns p, q */
                    ST_double dkp = D[(size_t)p * n + k], dkq = D[(size_t)q * n + k];
                    D[(size_t)p * n + k] = c * dkp - s * dkq;
                    D[(size_t)q * n + k] = s * dkp + c * dkq;
                }
                for (ST_int k = 0; k < n; k++) {   /* rows p, q */
                    ST_double dpk = D[(size_t)k * n + p], dqk = D[(size_t)k * n + q];
                    D[(size_t)k * n + p] = c * dpk - s * dqk;
                    D[(size_t)k * n + q] = s * dpk + c * dqk;
                }
                for (ST_int k = 0; k < n; k++) {
                    ST_double vkp = V[(size_t)p * n + k], vkq = V[(size_t)q * n + k];
                    V[(size_t)p * n + k] = c * vkp - s * vkq;
                    V[(size_t)q * n + k] = s * vkp + c * vkq;
                }
            }
        }
    }
    for (ST_int i = 0; i < n; i++) ord[i] = i;
    for (ST_int i = 1; i < n; i++) {        /* sort ascending */
        ST_int o = ord[i], j = i;
        while (j > 0 && D[(size_t)ord[j - 1] * n + ord[j - 1]] > D[(size_t)o * n + o]) { ord[j] = ord[j - 1]; j--; }
        ord[j] = o;
    }
    for (ST_int i = 0; i < n; i++) {
        ev[i] = D[(size_t)ord[i] * n + ord[i]];
        if (Q) memcpy(Q + (size_t)i * n, V + (size_t)ord[i] * n, (size_t)n * sizeof(ST_double));
    }
    free(D); free(V); free(ord);
    return 0;
}

/* Solve the general square system A X = B (n x n, nrhs columns) by Gaussian
   elimination with partial pivoting; B is overwritten with X */
static ST_int civ_gauss_solve(const ST_double *A, ST_int n, ST_double *B, ST_int nrhs)
{
    ST_double *M = (ST_double *)malloc((size_t)n * n * sizeof(ST_double));
    if (!M) return -1;
    memcpy(M, A, (size_t)n * n * sizeof(ST_double));
    for (ST_int k = 0; k < n; k++) {
        ST_int piv = k;
        for (ST_int i = k + 1; i < n; i++)
            if (fabs(M[(size_t)k * n + i]) > fabs(M[(size_t)k * n + piv])) piv = i;
        if (M[(size_t)k * n + piv] == 0.0) { free(M); return -1; }
        if (piv != k) {
            for (ST_int j = 0; j < n; j++) {
                ST_double t = M[(size_t)j * n + k]; M[(size_t)j * n + k] = M[(size_t)j * n + piv]; M[(size_t)j * n + piv] = t;
            }
            for (ST_int j = 0; j < nrhs; j++) {
                ST_double t = B[(size_t)j * n + k]; B[(size_t)j * n + k] = B[(size_t)j * n + piv]; B[(size_t)j * n + piv] = t;
            }
        }
        for (ST_int i = k + 1; i < n; i++) {
            ST_double f = M[(size_t)k * n + i] / M[(size_t)k * n + k];
            if (f == 0.0) continue;
            for (ST_int j = k; j < n; j++) M[(size_t)j * n + i] -= f * M[(size_t)j * n + k];
            for (ST_int j = 0; j < nrhs; j++) B[(size_t)j * n + i] -= f * B[(size_t)j * n + k];
        }
    }
    for (ST_int j = 0; j < nrhs; j++) {
        ST_double *b = B + (size_t)j * n;
        for (ST_int i = n - 1; i >= 0; i--) {
            ST_double s = b[i];
            for (ST_int k = i + 1; k < n; k++) s -= M[(size_t)k * n + i] * b[k];
            b[i] = s / M[(size_t)i * n + i];
        }
    }
    free(M);
    return 0;
}

/* Symmetric square root of the PSD n x n A (negative eigenvalues set to 0) */
static ST_int civ_sym_sqrt(const ST_double *A, ST_int n, ST_double *R)
{
    ST_double *ev = (ST_double *)malloc((size_t)n * sizeof(ST_double));
    ST_double *Q = (ST_double *)malloc((size_t)n * n * sizeof(ST_double));
    ST_int rc = (ev && Q) ? civ_sym_eigen(A, n, ev, Q) : -1;
    if (!rc) {
        for (ST_int j = 0; j < n; j++)
            for (ST_int i = 0; i < n; i++) {
                ST_double s = 0.0;
                for (ST_int k = 0; k < n; k++)
                    s += Q[(size_t)k * n + i] * (ev[k] > 0.0 ? sqrt(ev[k]) : 0.0) * Q[(size_t)k * n + j];
                R[(size_t)j * n + i] = s;
            }
    }
    free(ev); free(Q);
    return rc;
}

/* ------------------------------------------------------------------------
   GMM J statistics
   ------------------------------------------------------------------------ */

/*
    Efficient GMM with the weighting matrix S^-1 (S raw, L x L):
    beta = (X'WZ S^-1 Z'WX)^-1 X'WZ S^-1 Z'Wy, and J = g(e)' S^-1 g(e) at its
    residuals e (ivreg2's s_gmm1s with a supplied S, and the J of s_iegmm).
    Returns 0, or -1 if a matrix is singular or memory runs out.
*/
static ST_int civ_gmm_eff_j(const ST_double *y, const ST_double *X, ST_int K,
                            const ST_double *Z, ST_int L, const ST_double *w, ST_int N,
                            const ST_double *S, ST_double *J)
{
    ST_double *ZX = (ST_double *)calloc((size_t)L * K, sizeof(ST_double));
    ST_double *Zy = (ST_double *)calloc((size_t)L, sizeof(ST_double));
    ST_double *SiZX = (ST_double *)calloc((size_t)L * K, sizeof(ST_double));
    ST_double *Sig = (ST_double *)calloc((size_t)L, sizeof(ST_double));
    ST_double *H = (ST_double *)calloc((size_t)K * K, sizeof(ST_double));
    ST_double *b = (ST_double *)calloc((size_t)K, sizeof(ST_double));
    ST_double *e = (ST_double *)malloc((size_t)N * sizeof(ST_double));
    ST_int rc = (ZX && Zy && SiZX && Sig && H && b && e) ? 0 : -1;
    if (!rc) {
        civ_cross(Z, L, X, K, w, N, ZX);
        civ_zwe(Z, L, y, w, N, Zy);
        memcpy(SiZX, ZX, (size_t)L * K * sizeof(ST_double));
        memcpy(Sig, Zy, (size_t)L * sizeof(ST_double));
        rc = civ_sym_solve(S, L, SiZX, K, SiZX);
        if (!rc) rc = civ_sym_solve(S, L, Sig, 1, Sig);
    }
    if (!rc) {
        /* H = ZX' S^-1 ZX, b = H^-1 ZX' S^-1 Zy */
        for (ST_int c = 0; c < K; c++) {
            for (ST_int r = 0; r < K; r++) {
                ST_double s = 0.0;
                for (ST_int l = 0; l < L; l++) s += ZX[(size_t)r * L + l] * SiZX[(size_t)c * L + l];
                H[(size_t)c * K + r] = s;
            }
            ST_double s = 0.0;
            for (ST_int l = 0; l < L; l++) s += ZX[(size_t)c * L + l] * Sig[l];
            b[c] = s;
        }
        rc = civ_sym_solve(H, K, b, 1, b);
    }
    if (!rc) {
        for (ST_int i = 0; i < N; i++) {
            ST_double f = 0.0;
            for (ST_int k = 0; k < K; k++) f += X[(size_t)k * N + i] * b[k];
            e[i] = y[i] - f;
        }
        civ_zwe(Z, L, e, w, N, Zy);
        memcpy(Sig, Zy, (size_t)L * sizeof(ST_double));
        rc = civ_sym_solve(S, L, Sig, 1, Sig);
    }
    if (!rc) {
        ST_double j = 0.0;
        for (ST_int l = 0; l < L; l++) j += Zy[l] * Sig[l];
        *J = j;
    }
    free(ZX); free(Zy); free(SiZX); free(Sig); free(H); free(b); free(e);
    return rc;
}

/* J = g(e)' S^-1 g(e) for given residuals e */
static ST_int civ_j_direct(const ST_double *e, const ST_double *Z, ST_int L,
                           const ST_double *w, ST_int N, const ST_double *S, ST_double *J)
{
    ST_double *g = (ST_double *)calloc((size_t)L * 2, sizeof(ST_double));
    if (!g) return -1;
    civ_zwe(Z, L, e, w, N, g);
    ST_int rc = civ_sym_solve(S, L, g, 1, g + L);
    if (!rc) {
        ST_double j = 0.0;
        for (ST_int l = 0; l < L; l++) j += g[l] * g[L + l];
        *J = j;
    }
    free(g);
    return rc;
}

/* ------------------------------------------------------------------------
   Kleibergen-Paap rk statistic (ranktest's s_jstat with m_svd)
   ------------------------------------------------------------------------ */

/*
    Kleibergen-Paap rk statistic for H0: rank(Pi) = K - 1, the full-rank test
    behind ivreghdfe's idstat (LM) and widstat (Wald). Y (N x K) and Zx
    (N x L, L >= K) have the included exogenous regressors partialled out;
    om gives the moment covariance (ranktest's m_omega). With
    Qyy = Y'WY/N = Ly Ly', Qzz = Zx'WZx/N = Lz Lz', Pi = Qzz^-1 Zx'WY/N and
    shat0 = S(vhat, Zx)/N (vhat = Y for LM, Y - Zx Pi for Wald):
      that  = Lz' Pi Ly'^-1,  kpvar = (Ly^-1 # Lz^-1) shat0 (Ly^-1 # Lz^-1)'
    and, from the SVD that = U diag(s) V' with U2 = the last L-K+1 columns of
    U, u22 its last L-K+1 rows and v = the last column of V,
      aq = U2 u22^-1 (u22 u22')^1/2,  bq = sign(v[K]) v'
      rk = N lab' (Q kpvar Q')^- lab,  lab = aq' that bq', Q = bq # aq'
    U2 only enters through the orthogonal complement of the first K-1 left
    singular vectors, and v up to sign, so both come from symmetric
    eigen-decompositions of that*that' and that'*that.
    Returns 0 with the statistic in *stat, or -1.
*/
static ST_int civ_kp_rk(const ST_double *Y, ST_int K, const ST_double *Zx, ST_int L,
                        const ST_double *w, ST_int N, const civ_omega *om,
                        ST_int wald, ST_double *stat)
{
    ST_double NN = om->N_eff;
    ST_int D = K * L, m = L - K + 1, rc = 0;
    if (K <= 0 || L < K) return -1;

    ST_double *Qyy = (ST_double *)calloc((size_t)K * K, sizeof(ST_double));
    ST_double *Qzz = (ST_double *)calloc((size_t)L * L, sizeof(ST_double));
    ST_double *Pi = (ST_double *)calloc((size_t)L * K, sizeof(ST_double));
    ST_double *Ly = (ST_double *)calloc((size_t)K * K, sizeof(ST_double));
    ST_double *Lz = (ST_double *)calloc((size_t)L * L, sizeof(ST_double));
    ST_double *Ay = (ST_double *)calloc((size_t)K * K, sizeof(ST_double));
    ST_double *Az = (ST_double *)calloc((size_t)L * L, sizeof(ST_double));
    ST_double *that = (ST_double *)calloc((size_t)L * K, sizeof(ST_double));
    ST_double *S = (ST_double *)calloc((size_t)D * D, sizeof(ST_double));
    ST_double *B = (ST_double *)calloc((size_t)D * D, sizeof(ST_double));
    ST_double *T = (ST_double *)calloc((size_t)D * D, sizeof(ST_double));
    ST_double *Mu = (ST_double *)calloc((size_t)L * L, sizeof(ST_double));
    ST_double *Uq = (ST_double *)calloc((size_t)L * L, sizeof(ST_double));
    ST_double *evu = (ST_double *)calloc((size_t)L, sizeof(ST_double));
    ST_double *Mv = (ST_double *)calloc((size_t)K * K, sizeof(ST_double));
    ST_double *Vq = (ST_double *)calloc((size_t)K * K, sizeof(ST_double));
    ST_double *evv = (ST_double *)calloc((size_t)K, sizeof(ST_double));
    ST_double *u22 = (ST_double *)calloc((size_t)m * m, sizeof(ST_double));
    ST_double *uu = (ST_double *)calloc((size_t)m * m, sizeof(ST_double));
    ST_double *u22h = (ST_double *)calloc((size_t)m * m, sizeof(ST_double));
    ST_double *aq = (ST_double *)calloc((size_t)L * m, sizeof(ST_double));
    ST_double *Qm = (ST_double *)calloc((size_t)m * D, sizeof(ST_double));
    ST_double *vlab = (ST_double *)calloc((size_t)m * m, sizeof(ST_double));
    ST_double *lab = (ST_double *)calloc((size_t)m * 2, sizeof(ST_double));
    ST_double *vhat = wald ? (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K, sizeof(ST_double)) : NULL;
    if (!Qyy || !Qzz || !Pi || !Ly || !Lz || !Ay || !Az || !that || !S || !B || !T || !Mu ||
        !Uq || !evu || !Mv || !Vq || !evv || !u22 || !uu || !u22h || !aq || !Qm || !vlab ||
        !lab || (wald && !vhat)) rc = -1;

    if (!rc) {
        civ_cross(Y, K, Y, K, w, N, Qyy);
        civ_cross(Zx, L, Zx, L, w, N, Qzz);
        civ_cross(Zx, L, Y, K, w, N, Pi);
        for (ST_int i = 0; i < K * K; i++) Qyy[i] /= NN;
        for (ST_int i = 0; i < L * L; i++) Qzz[i] /= NN;
        for (ST_int i = 0; i < L * K; i++) Pi[i] /= NN;
        rc = civ_sym_solve(Qzz, L, Pi, K, Pi);
    }
    if (!rc && wald) {
        /* reduced-form residuals */
        memcpy(vhat, Y, (size_t)N * K * sizeof(ST_double));
        for (ST_int k = 0; k < K; k++)
            for (ST_int l = 0; l < L; l++) {
                ST_double p = Pi[(size_t)k * L + l];
                const ST_double *z = Zx + (size_t)l * N;
                ST_double *v = vhat + (size_t)k * N;
                for (ST_int i = 0; i < N; i++) v[i] -= z[i] * p;
            }
    }
    if (!rc) rc = civ_omega_build(om, wald ? vhat : Y, K, Zx, L, S) ? -1 : 0;
    if (!rc && K == 1) {
        /* One endogenous regressor: the full-rank test is ranktest's test
           of rank 0, N g' shat0^- g with g = Zx'WY/N and shat0 = S/N
           (algebraically the statistic below; it differs only when shat0
           is singular, e.g. with too few clusters, where both use invsym) */
        ST_double *g = lab, *sg = lab + m;   /* m = L when K = 1 */
        civ_cross(Zx, L, Y, 1, w, N, g);
        for (ST_int l = 0; l < L; l++) g[l] /= NN;
        for (ST_int i = 0; i < D * D; i++) S[i] /= NN;
        rc = civ_sym_solve(S, L, g, 1, sg);
        if (!rc) {
            ST_double q = 0.0;
            for (ST_int l = 0; l < L; l++) q += g[l] * sg[l];
            *stat = NN * q;
        }
        goto kp_done;
    }
    if (!rc) rc = civ_chol_lower(Qyy, K, Ly);
    if (!rc) rc = civ_chol_lower(Qzz, L, Lz);
    if (!rc) {
        civ_lower_inverse(Ly, K, Ay);
        civ_lower_inverse(Lz, L, Az);
        /* that = Lz' Pi Ay' */
        for (ST_int k = 0; k < K; k++)
            for (ST_int l = 0; l < L; l++) {
                ST_double s = 0.0;
                for (ST_int j = 0; j < K; j++) {
                    ST_double lp = 0.0;   /* (Lz' Pi)[l, j] */
                    for (ST_int r = l; r < L; r++) lp += Lz[(size_t)l * L + r] * Pi[(size_t)j * L + r];
                    s += lp * Ay[(size_t)j * K + k];
                }
                that[(size_t)k * L + l] = s;
            }
        /* kpvar = B shat0 B', B = Ay # Az, shat0 = S/N */
        for (ST_int k2 = 0; k2 < K; k2++)
            for (ST_int l2 = 0; l2 < L; l2++)
                for (ST_int k1 = 0; k1 < K; k1++)
                    for (ST_int l1 = 0; l1 < L; l1++)
                        B[(size_t)(k2 * L + l2) * D + (k1 * L + l1)] =
                            Ay[(size_t)k2 * K + k1] * Az[(size_t)l2 * L + l1];
        for (ST_int i = 0; i < D * D; i++) S[i] /= NN;
        ctools_matmul_ab(B, S, D, D, D, T);           /* T = B shat0 */
        for (ST_int c = 0; c < D; c++)                 /* S <- T B' */
            for (ST_int r = 0; r < D; r++) {
                ST_double s = 0.0;
                for (ST_int j = 0; j < D; j++) s += T[(size_t)j * D + r] * B[(size_t)j * D + c];
                S[(size_t)c * D + r] = s;
            }

        /* Singular vectors: v for the smallest singular value, U2 spanning
           the complement of the first K-1 left singular vectors */
        for (ST_int c = 0; c < K; c++)
            for (ST_int r = 0; r < K; r++) {
                ST_double s = 0.0;
                for (ST_int l = 0; l < L; l++) s += that[(size_t)r * L + l] * that[(size_t)c * L + l];
                Mv[(size_t)c * K + r] = s;
            }
        for (ST_int c = 0; c < L; c++)
            for (ST_int r = 0; r < L; r++) {
                ST_double s = 0.0;
                for (ST_int k = 0; k < K; k++) s += that[(size_t)k * L + r] * that[(size_t)k * L + c];
                Mu[(size_t)c * L + r] = s;
            }
        rc = civ_sym_eigen(Mv, K, evv, Vq);
        if (!rc) rc = civ_sym_eigen(Mu, L, evu, Uq);
    }
    if (!rc) {
        /* U2 = Uq[:, 0..m-1]; u22 = its rows K-1..L-1 */
        for (ST_int c = 0; c < m; c++)
            for (ST_int r = 0; r < m; r++) u22[(size_t)c * m + r] = Uq[(size_t)c * L + (K - 1 + r)];
        for (ST_int c = 0; c < m; c++)
            for (ST_int r = 0; r < m; r++) {
                ST_double s = 0.0;
                for (ST_int j = 0; j < m; j++) s += u22[(size_t)j * m + r] * u22[(size_t)j * m + c];
                uu[(size_t)c * m + r] = s;
            }
        rc = civ_sym_sqrt(uu, m, u22h);
        /* u22h <- u22^-1 u22h, then aq = U2 u22h */
        if (!rc) rc = civ_gauss_solve(u22, m, u22h, m);
    }
    ST_double v22 = (!rc) ? Vq[K - 1] : 0.0;   /* last element of v (column 0) */
    if (!rc && v22 == 0.0) rc = -1;
    if (!rc) {
        ST_double sg = (v22 > 0.0) ? 1.0 : -1.0;
        for (ST_int c = 0; c < m; c++)
            for (ST_int l = 0; l < L; l++) {
                ST_double s = 0.0;
                for (ST_int j = 0; j < m; j++) s += Uq[(size_t)j * L + l] * u22h[(size_t)c * m + j];
                aq[(size_t)c * L + l] = s;
            }
        /* lab = aq' that bq', Q = bq # aq', bq = sg * v' */
        for (ST_int r = 0; r < m; r++) {
            ST_double s = 0.0;
            for (ST_int k = 0; k < K; k++) {
                ST_double bk = sg * Vq[k];
                for (ST_int l = 0; l < L; l++) {
                    s += aq[(size_t)r * L + l] * that[(size_t)k * L + l] * bk;
                    Qm[(size_t)(k * L + l) * m + r] = bk * aq[(size_t)r * L + l];
                }
            }
            lab[r] = s;
        }
        /* vlab = Q kpvar Q' */
        ST_double *QS = T;   /* m x D scratch */
        for (ST_int c = 0; c < D; c++)
            for (ST_int r = 0; r < m; r++) {
                ST_double s = 0.0;
                for (ST_int j = 0; j < D; j++) s += Qm[(size_t)j * m + r] * S[(size_t)c * D + j];
                QS[(size_t)c * m + r] = s;
            }
        for (ST_int c = 0; c < m; c++)
            for (ST_int r = 0; r < m; r++) {
                ST_double s = 0.0;
                for (ST_int j = 0; j < D; j++) s += QS[(size_t)j * m + r] * Qm[(size_t)j * m + c];
                vlab[(size_t)c * m + r] = s;
            }
        rc = civ_sym_solve(vlab, m, lab, 1, lab + m);
    }
    if (!rc) {
        ST_double q = 0.0;
        for (ST_int r = 0; r < m; r++) q += lab[r] * lab[m + r];
        *stat = NN * q;
    }

kp_done:
    free(Qyy); free(Qzz); free(Pi); free(Ly); free(Lz); free(Ay); free(Az); free(that);
    free(S); free(B); free(T); free(Mu); free(Uq); free(evu); free(Mv); free(Vq); free(evv);
    free(u22); free(uu); free(u22h); free(aq); free(Qm); free(vlab); free(lab); free(vhat);
    return rc;
}

/*
    Compute underidentification test (Anderson LM / Kleibergen-Paap rk LM)
    and weak identification statistics (Cragg-Donald F / KP rk Wald F).
*/
void civreghdfe_compute_underid_test(
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int N_eff,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int df_a,
    const civ_omega *rk_om,
    ST_int rk_cluster,
    ST_int rk_nclust,
    ST_double *underid_stat,
    ST_int *underid_df,
    ST_double *cd_f,
    ST_double *kp_f
)
{
    ST_int L = K_iv - K_exog;  /* Number of excluded instruments */
    const ST_double *w = (weights && weight_type != 0) ? weights : NULL;

    *underid_stat = 0.0;
    *underid_df = L - K_endog + 1;
    if (*underid_df < 1) *underid_df = 1;
    *cd_f = 0.0;
    *kp_f = 0.0;

    if (K_endog <= 0 || L < K_endog) return;

    ST_int df_resid = N_eff - K_iv - df_a;
    if (df_resid <= 0) df_resid = 1;

    /* ranktest partials the included exogenous regressors out of the
       endogenous regressors and the excluded instruments */
    ST_double *Y = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_endog, sizeof(ST_double));
    ST_double *Zx = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)L, sizeof(ST_double));
    if (!Y || !Zx ||
        civ_partial_cols(X_endog, K_endog, Z, K_exog, w, N, Y) != 0 ||
        civ_partial_cols(Z + (size_t)K_exog * N, L, Z, K_exog, w, N, Zx) != 0) {
        free(Y); free(Zx);
        return;
    }
    civ_standardize_cols(Y, K_endog, N);
    civ_standardize_cols(Zx, L, N);

    {
        /* Squared canonical correlations = eigenvalues of Qyy^-1/2 Qzy'
           Qzz^-1 Qzy Qyy^-1/2; Anderson LM = N * the smallest, Cragg-Donald
           F = ev/(1-ev) * df_resid / L (ivreghdfe's idstat and cdf). With one
           endogenous regressor ev is the partial R2 of the excluded
           instruments; the first-stage F is not used because its residual
           dof keeps nested fixed effects out, which ivreghdfe's cdf does
           not. */
        ST_double *Qyy = (ST_double *)calloc((size_t)K_endog * K_endog, sizeof(ST_double));
        ST_double *Qzz = (ST_double *)calloc((size_t)L * L, sizeof(ST_double));
        ST_double *Qzy = (ST_double *)calloc((size_t)L * K_endog, sizeof(ST_double));
        ST_double *P = (ST_double *)calloc((size_t)L * K_endog, sizeof(ST_double));
        ST_double *M = (ST_double *)calloc((size_t)K_endog * K_endog, sizeof(ST_double));
        ST_double *ev = (ST_double *)calloc((size_t)K_endog, sizeof(ST_double));
        ST_int ok = (Qyy && Qzz && Qzy && P && M && ev);
        if (ok) {
            civ_cross(Y, K_endog, Y, K_endog, w, N, Qyy);
            civ_cross(Zx, L, Zx, L, w, N, Qzz);
            civ_cross(Zx, L, Y, K_endog, w, N, Qzy);
            memcpy(P, Qzy, (size_t)L * K_endog * sizeof(ST_double));
            ok = (civ_sym_solve(Qzz, L, P, K_endog, P) == 0);
        }
        if (ok) {
            /* M = Qzy' Qzz^-1 Qzy; the eigenvalues of Qyy^-1/2 M Qyy^-1/2
               are those of Ay M Ay' with Ay = chol(Qyy)^-1 */
            for (ST_int c = 0; c < K_endog; c++)
                for (ST_int r = 0; r < K_endog; r++) {
                    ST_double s = 0.0;
                    for (ST_int l = 0; l < L; l++) s += Qzy[(size_t)r * L + l] * P[(size_t)c * L + l];
                    M[(size_t)c * K_endog + r] = s;
                }
            ST_double *Ly = (ST_double *)calloc((size_t)K_endog * K_endog, sizeof(ST_double));
            ST_double *Ay = (ST_double *)calloc((size_t)K_endog * K_endog, sizeof(ST_double));
            ST_double *T = (ST_double *)calloc((size_t)K_endog * K_endog, sizeof(ST_double));
            ok = Ly && Ay && T && civ_chol_lower(Qyy, K_endog, Ly) == 0;
            if (ok) {
                civ_lower_inverse(Ly, K_endog, Ay);
                ctools_matmul_ab(Ay, M, K_endog, K_endog, K_endog, T);
                for (ST_int c = 0; c < K_endog; c++)
                    for (ST_int r = 0; r < K_endog; r++) {
                        ST_double s = 0.0;
                        for (ST_int j = 0; j < K_endog; j++)
                            s += T[(size_t)j * K_endog + r] * Ay[(size_t)j * K_endog + c];
                        M[(size_t)c * K_endog + r] = s;
                    }
                ok = (civ_sym_eigen(M, K_endog, ev, NULL) == 0);
            }
            free(Ly); free(Ay); free(T);
        }
        if (ok) {
            ST_double min_eval = ev[0];
            *underid_stat = (ST_double)N_eff * min_eval;
            if (min_eval > 0.0 && min_eval < 1.0)
                *cd_f = min_eval / (1.0 - min_eval) * ((ST_double)df_resid / (ST_double)L);
        }
        free(Qyy); free(Qzz); free(Qzy); free(P); free(M); free(ev);
    }

    if (rk_om->kind == CIV_OMEGA_IID && rk_om->kernel <= 0) {
        /* ranktest without robust/cluster/bw (default VCE, and kiefer, whose
           bw ivreghdfe sets after parsing): its Wald statistic is the
           Cragg-Donald statistic */
        *kp_f = *cd_f;
    } else {
        /* Kleibergen-Paap rk statistics with ranktest's S (robust, cluster,
           two-way, HAC/AC/Driscoll-Kraay with the Bartlett kernel): LM for
           idstat, Wald for widstat, converted as ivreghdfe does */
        ST_double lm = 0.0, wald = 0.0;
        if (civ_kp_rk(Y, K_endog, Zx, L, w, N, rk_om, 0, &lm) == 0) *underid_stat = lm;
        if (civ_kp_rk(Y, K_endog, Zx, L, w, N, rk_om, 1, &wald) == 0) {
            if (rk_cluster) {
                ST_double G = (ST_double)rk_nclust;
                *kp_f = wald / (ST_double)(N_eff - 1) * (ST_double)df_resid *
                        (G > 0 ? (G - 1.0) / G : 1.0) / (ST_double)L;
            } else {
                *kp_f = wald / (ST_double)N_eff * (ST_double)df_resid / (ST_double)L;
            }
        }
    }
    free(Y); free(Zx);
}

/*
    Compute Sargan/Hansen J overidentification test: J = g(e)' S^-1 g(e),
    g = Z'We, with the model's residuals (resid != NULL: the Sargan statistic
    when S is iid, and the J of efficient GMM2S / CUE) or with the residuals
    of efficient GMM for S (resid == NULL: the Hansen J of an inefficient
    2SLS / k-class fit, ivreg2's s_iegmm and s_liml).
*/
void civreghdfe_compute_hansen_j(
    const ST_double *y,
    const ST_double *X_all,
    ST_int K_total,
    const ST_double *Z,
    ST_int K_iv,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    const ST_double *S,
    const ST_double *resid,
    ST_double *J,
    ST_int *overid_df
)
{
    const ST_double *w = (weights && weight_type != 0) ? weights : NULL;
    *overid_df = K_iv - K_total;
    *J = 0.0;
    if (*overid_df <= 0) return;
    ST_double j = 0.0;
    ST_int rc = resid ? civ_j_direct(resid, Z, K_iv, w, N, S, &j)
                      : civ_gmm_eff_j(y, X_all, K_total, Z, K_iv, w, N, S, &j);
    if (rc == 0) *J = j;
}

/*
    Compute Durbin-Wu-Hausman endogeneity test.
*/
void civreghdfe_compute_dwh_test(
    const ST_double *y,
    const ST_double *X_exog,
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *temp1,
    ST_int N,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int df_a,
    ST_double *endog_chi2,
    ST_double *endog_f,
    ST_int *endog_df
)
{
    ST_int K_total = K_exog + K_endog;
    ST_int i, j, k;

    *endog_chi2 = 0.0;
    *endog_f = 0.0;
    *endog_df = K_endog;

    if (K_endog <= 0) return;

    /* Compute first-stage residuals: v = X_endog - Z * pi */
    ST_double *v_resid = (ST_double *)calloc((size_t)N * K_endog, sizeof(ST_double));
    if (!v_resid) return;

    for (ST_int e = 0; e < K_endog; e++) {
        const ST_double *x_endog_col = X_endog + (size_t)e * N;
        ST_double *v_col = v_resid + (size_t)e * N;
        const ST_double *pi_col = temp1 + (K_exog + e) * K_iv;

        for (i = 0; i < N; i++) {
            ST_double pred = 0.0;
            for (k = 0; k < K_iv; k++) {
                pred += Z[(size_t)k * N + i] * pi_col[k];
            }
            v_col[i] = x_endog_col[i] - pred;
        }
    }

    /* Build augmented design matrix: [X_exog, X_endog, v] */
    ST_int K_aug = K_total + K_endog;
    ST_double *X_aug = (ST_double *)calloc((size_t)N * K_aug, sizeof(ST_double));
    ST_double *XaXa = (ST_double *)calloc(K_aug * K_aug, sizeof(ST_double));
    ST_double *Xay = (ST_double *)calloc(K_aug, sizeof(ST_double));

    if (!X_aug || !XaXa || !Xay) {
        free(v_resid);
        free(X_aug); free(XaXa); free(Xay);
        return;
    }

    /* Copy X_exog, X_endog, v_resid */
    if (X_exog) {
        for (j = 0; j < K_exog; j++) {
            for (i = 0; i < N; i++) {
                X_aug[(size_t)j * N + i] = X_exog[(size_t)j * N + i];
            }
        }
    }
    if (X_endog) {
        for (j = 0; j < K_endog; j++) {
            for (i = 0; i < N; i++) {
                X_aug[(size_t)(K_exog + j) * N + i] = X_endog[(size_t)j * N + i];
            }
        }
    }
    for (j = 0; j < K_endog; j++) {
        for (i = 0; i < N; i++) {
            X_aug[(size_t)(K_total + j) * N + i] = v_resid[(size_t)j * N + i];
        }
    }

    /* Compute X_aug'X_aug */
    ctools_matmul_atb(X_aug, X_aug, N, K_aug, K_aug, XaXa);

    /* Compute X_aug'y */
    for (j = 0; j < K_aug; j++) {
        ST_double sum = 0.0;
        const ST_double *xa_col = X_aug + (size_t)j * N;
        for (i = 0; i < N; i++) {
            sum += xa_col[i] * y[i];
        }
        Xay[j] = sum;
    }

    /* Solve for augmented OLS coefficients */
    ST_double *XaXa_L = (ST_double *)calloc(K_aug * K_aug, sizeof(ST_double));
    ST_double *XaXa_inv = (ST_double *)calloc(K_aug * K_aug, sizeof(ST_double));
    ST_double *beta_aug = (ST_double *)calloc(K_aug, sizeof(ST_double));

    if (XaXa_L && XaXa_inv && beta_aug) {
        memcpy(XaXa_L, XaXa, K_aug * K_aug * sizeof(ST_double));
        if (ctools_cholesky(XaXa_L, K_aug) == 0) {
            ctools_invert_from_cholesky(XaXa_L, K_aug, XaXa_inv);

            /* beta_aug = XaXa_inv * Xay */
            for (i = 0; i < K_aug; i++) {
                ST_double sum = 0.0;
                for (k = 0; k < K_aug; k++) {
                    sum += XaXa_inv[k * K_aug + i] * Xay[k];
                }
                beta_aug[i] = sum;
            }

            /* Compute residuals and sigma^2 */
            ST_double sse_aug = 0.0;
            for (i = 0; i < N; i++) {
                ST_double fitted = 0.0;
                for (j = 0; j < K_aug; j++) {
                    fitted += X_aug[(size_t)j * N + i] * beta_aug[j];
                }
                ST_double r = y[i] - fitted;
                sse_aug += r * r;
            }

            ST_int df_aug = N - K_aug - df_a;
            if (df_aug <= 0) df_aug = 1;
            ST_double sigma2_aug = sse_aug / df_aug;

            /* gamma = beta_aug[K_total:K_aug-1] */
            ST_double *gamma = beta_aug + K_total;

            /* Extract v'v block from XaXa */
            ST_double *XaXa_vv = (ST_double *)calloc(K_endog * K_endog, sizeof(ST_double));
            if (XaXa_vv) {
                for (j = 0; j < K_endog; j++) {
                    for (i = 0; i < K_endog; i++) {
                        XaXa_vv[j * K_endog + i] = XaXa[(K_total + j) * K_aug + (K_total + i)];
                    }
                }

                /* chi2 = gamma' * XaXa_vv * gamma / sigma2_aug */
                ST_double quad = 0.0;
                for (j = 0; j < K_endog; j++) {
                    ST_double sum = 0.0;
                    for (i = 0; i < K_endog; i++) {
                        sum += XaXa_vv[j * K_endog + i] * gamma[i];
                    }
                    quad += gamma[j] * sum;
                }

                *endog_chi2 = quad / sigma2_aug;
                *endog_f = (*endog_chi2) / K_endog;

                free(XaXa_vv);
            }
        }
    }

    free(v_resid);
    free(X_aug);
    free(XaXa);
    free(Xay);
    free(XaXa_L);
    free(XaXa_inv);
    free(beta_aug);
}

/*
    Compute C-statistic for testing instrument orthogonality (ivreg2's
    orthog()): C = J - J_r, where J_r is the efficient GMM J of the model
    without the tested instruments, weighted by the corresponding block of
    the full model's S, as ivreg2 passes smatrix() to its restricted fit.
    With the iid S = sigma^2 Z'WZ this is the difference of Sargan
    statistics that both use the full model's sigma^2.
*/
void civreghdfe_compute_cstat(
    const ST_double *y,
    const ST_double *X_all,
    ST_int K_total,
    const ST_double *Z,
    ST_int K_iv,
    ST_int K_exog,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    const ST_double *S,
    ST_double J_full,
    const ST_int *orthog_indices,
    ST_int n_orthog,
    ST_double *cstat,
    ST_int *cstat_df
)
{
    const ST_double *w = (weights && weight_type != 0) ? weights : NULL;
    ST_int K_rest = K_iv - n_orthog;  /* Instruments not being tested */

    *cstat = 0.0;
    *cstat_df = n_orthog;

    /* A restricted model that is not identified has no C statistic (ivreg2
       posts 0, as when its recursive call fails with r(481)) */
    if (n_orthog <= 0 || K_rest < K_total) return;

    /* Keep mask: the tested columns are distinct, in-range excluded
       instruments (mapped in civreghdfe_impl.c); refuse anything else */
    ST_int *keep = (ST_int *)malloc((size_t)K_iv * sizeof(ST_int));
    if (!keep) return;
    for (ST_int k = 0; k < K_iv; k++) keep[k] = 1;
    ST_int n_masked = 0;
    for (ST_int i = 0; i < n_orthog; i++) {
        ST_int idx = orthog_indices[i] - 1;
        if (idx >= 0 && idx < K_iv - K_exog && keep[K_exog + idx]) {
            keep[K_exog + idx] = 0;
            n_masked++;
        }
    }
    ST_double *Zr = NULL, *Srr = NULL;
    if (n_masked == n_orthog) {
        Zr = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_rest, sizeof(ST_double));
        Srr = (ST_double *)calloc((size_t)K_rest * K_rest, sizeof(ST_double));
    }
    if (Zr && Srr) {
        /* Z_r and S_rr: the retained instruments and their block of S */
        ST_int a = 0;
        for (ST_int k = 0; k < K_iv; k++) {
            if (!keep[k]) continue;
            memcpy(Zr + (size_t)a * N, Z + (size_t)k * N, (size_t)N * sizeof(ST_double));
            ST_int b = 0;
            for (ST_int k2 = 0; k2 < K_iv; k2++) {
                if (!keep[k2]) continue;
                Srr[(size_t)b * K_rest + a] = S[(size_t)k2 * K_iv + k];
                b++;
            }
            a++;
        }
        ST_double Jr = 0.0;
        if (civ_gmm_eff_j(y, X_all, K_total, Zr, K_rest, w, N, Srr, &Jr) == 0) {
            *cstat = J_full - Jr;
            if (*cstat < 0) *cstat = 0;
        }
    }
    free(keep); free(Zr); free(Srr);
}

/*
    Compute endogeneity test for a subset of endogenous regressors (ivreg2's
    endog()): the C statistic of the model that treats the tested regressors
    as exogenous (they join the instruments). That model is re-estimated
    with the same estimator (2SLS, LIML or GMM2S, as ivreg2 passes liml and
    gmm2s but not cue or kclass) and the S structure in om; its own S comes
    from its residuals (GMM2S: its first step), and
      stat = J_exog - J_r,
    J_r = efficient GMM J of the original model weighted by the block of
    the exogenous model's S for the original instruments. With the iid S
    both J's use the exogenous model's sigma^2.
*/
void civreghdfe_compute_endogtest_subset(
    const ST_double *y,
    const ST_double *X_exog,
    const ST_double *X_endog,
    const ST_double *Z,
    ST_int N,
    const ST_double *weights,
    ST_int weight_type,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int est_type,
    const civ_omega *om,
    const ST_int *endogtest_indices,
    ST_int n_endogtest,
    ST_double *endogtest_stat,
    ST_int *endogtest_df
)
{
    const ST_double *w = (weights && weight_type != 0) ? weights : NULL;
    ST_int K_total = K_exog + K_endog;
    ST_int K_aug = K_iv + n_endogtest;
    ST_int iid = (om->kind == CIV_OMEGA_IID && om->kernel <= 0);

    *endogtest_stat = 0.0;
    *endogtest_df = n_endogtest;

    if (n_endogtest <= 0 || n_endogtest > K_endog) return;

    ST_double *X_all = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_total, sizeof(ST_double));
    ST_double *Z_aug = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_aug, sizeof(ST_double));
    ST_double *e = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double));
    ST_double *S_aug = (ST_double *)calloc((size_t)K_aug * K_aug, sizeof(ST_double));
    ST_double *S_r = (ST_double *)calloc((size_t)K_iv * K_iv, sizeof(ST_double));
    ST_double *ZZ = (ST_double *)calloc((size_t)K_aug * K_aug, sizeof(ST_double));
    ST_double *ZX = (ST_double *)calloc((size_t)K_aug * K_total, sizeof(ST_double));
    ST_double *Zy = (ST_double *)calloc((size_t)K_aug, sizeof(ST_double));
    ST_double *A = (ST_double *)calloc((size_t)K_aug * K_total, sizeof(ST_double));
    ST_double *XkX = (ST_double *)calloc((size_t)K_total * K_total, sizeof(ST_double));
    ST_double *b = (ST_double *)calloc((size_t)K_total, sizeof(ST_double));
    ST_double *X_ex = NULL, *X_en = NULL;
    ST_double kclass = 1.0, J_aug = 0.0, J_r = 0.0;
    if (!X_all || !Z_aug || !e || !S_aug || !S_r || !ZZ || !ZX || !Zy || !A || !XkX || !b)
        goto endogtest_cleanup;

    if (K_exog > 0) memcpy(X_all, X_exog, (size_t)N * K_exog * sizeof(ST_double));
    memcpy(X_all + (size_t)N * K_exog, X_endog, (size_t)N * K_endog * sizeof(ST_double));

    /* Z_aug = [Z, tested endogenous regressors] (distinct, in range: a
       repeated column would make Z_aug singular) */
    memcpy(Z_aug, Z, (size_t)N * K_iv * sizeof(ST_double));
    for (ST_int t = 0; t < n_endogtest; t++) {
        ST_int idx = endogtest_indices[t] - 1;
        if (idx < 0 || idx >= K_endog) goto endogtest_cleanup;
        for (ST_int u = 0; u < t; u++)
            if (endogtest_indices[u] - 1 == idx) goto endogtest_cleanup;
        memcpy(Z_aug + (size_t)(K_iv + t) * N, X_endog + (size_t)idx * N, (size_t)N * sizeof(ST_double));
    }

    /* k-class fit of the exogenous model: k = 1 (2SLS, and GMM2S's first
       step) or its LIML lambda (tested regressors exogenous) */
    if (est_type == 1 && K_aug > K_total) {
        ST_int K_en = K_endog - n_endogtest, K_ex = K_exog + n_endogtest;
        X_ex = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_ex, sizeof(ST_double));
        X_en = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)(K_en > 0 ? K_en : 1), sizeof(ST_double));
        if (!X_ex || !X_en) goto endogtest_cleanup;
        if (K_exog > 0) memcpy(X_ex, X_exog, (size_t)N * K_exog * sizeof(ST_double));
        memcpy(X_ex + (size_t)N * K_exog, Z_aug + (size_t)N * K_iv, (size_t)N * n_endogtest * sizeof(ST_double));
        ST_int ne = 0;
        for (ST_int k = 0; k < K_endog; k++) {
            ST_int tested = 0;
            for (ST_int t = 0; t < n_endogtest; t++) if (endogtest_indices[t] - 1 == k) tested = 1;
            if (!tested) memcpy(X_en + (size_t)(ne++) * N, X_endog + (size_t)k * N, (size_t)N * sizeof(ST_double));
        }
        kclass = civreghdfe_compute_liml_lambda(y, X_en, X_ex, Z_aug, weights, weight_type,
                                                N, K_ex, K_en, K_aug);
    }
    civ_cross(Z_aug, K_aug, Z_aug, K_aug, w, N, ZZ);
    civ_cross(Z_aug, K_aug, X_all, K_total, w, N, ZX);
    civ_zwe(Z_aug, K_aug, y, w, N, Zy);
    memcpy(A, ZX, (size_t)K_aug * K_total * sizeof(ST_double));
    if (civ_sym_solve(ZZ, K_aug, A, K_total, A) != 0) goto endogtest_cleanup;   /* (Z'Z)^-1 Z'X */
    /* XkX = (1-k) X'WX + k X'Pz X, b = (1-k) X'Wy + k X'Pz y */
    if (kclass != 1.0) {
        civ_cross(X_all, K_total, X_all, K_total, w, N, XkX);
        civ_zwe(X_all, K_total, y, w, N, b);
    }
    for (ST_int c = 0; c < K_total; c++) {
        for (ST_int r = 0; r < K_total; r++) {
            ST_double s = 0.0;
            for (ST_int l = 0; l < K_aug; l++) s += ZX[(size_t)r * K_aug + l] * A[(size_t)c * K_aug + l];
            XkX[(size_t)c * K_total + r] = (1.0 - kclass) * XkX[(size_t)c * K_total + r] + kclass * s;
        }
        ST_double s = 0.0;
        for (ST_int l = 0; l < K_aug; l++) s += A[(size_t)c * K_aug + l] * Zy[l];
        b[c] = (1.0 - kclass) * b[c] + kclass * s;
    }
    if (civ_sym_solve(XkX, K_total, b, 1, b) != 0) goto endogtest_cleanup;
    for (ST_int i = 0; i < N; i++) {
        ST_double f = 0.0;
        for (ST_int k = 0; k < K_total; k++) f += X_all[(size_t)k * N + i] * b[k];
        e[i] = y[i] - f;
    }

    /* The exogenous model's S and J: iid Sargan with its sigma^2; otherwise
       the efficient GMM J (GMM2S's own J, or the 2SLS/LIML Hansen J) */
    if (civ_omega_build(om, e, 1, Z_aug, K_aug, S_aug) != 0) goto endogtest_cleanup;
    /* ivreg2 reports no C statistic when that S is rank deficient */
    if (!iid && civ_sym_rank(S_aug, K_aug) < K_aug) {
        *endogtest_stat = SV_missval;
        goto endogtest_cleanup;
    }
    if ((iid ?civ_j_direct(e, Z_aug, K_aug, w, N, S_aug, &J_aug)
             : civ_gmm_eff_j(y, X_all, K_total, Z_aug, K_aug, w, N, S_aug, &J_aug)) != 0)
        goto endogtest_cleanup;

    /* Original model, efficient GMM with the Z block of S_aug */
    for (ST_int c = 0; c < K_iv; c++)
        memcpy(S_r + (size_t)c * K_iv, S_aug + (size_t)c * K_aug, (size_t)K_iv * sizeof(ST_double));
    if (civ_gmm_eff_j(y, X_all, K_total, Z, K_iv, w, N, S_r, &J_r) != 0) goto endogtest_cleanup;

    *endogtest_stat = J_aug - J_r;
    if (*endogtest_stat < 0) *endogtest_stat = 0;

endogtest_cleanup:
    free(X_all); free(Z_aug); free(e); free(S_aug); free(S_r); free(ZZ);
    free(ZX); free(Zy); free(A); free(XkX); free(b); free(X_ex); free(X_en);
}

/*
    Compute instrument redundancy test (ivreghdfe's redundant()): ranktest's
    LM test of H0: rank = 0 for the tested instruments Z_t and the
    endogenous regressors Y, both with the included exogenous regressors and
    the other excluded instruments partialled out:
      stat = vec(Z_t'WY)' S^-1 vec(Z_t'WY), S = m_omega(e = Y, Z = Z_t)
    With the iid S this is N times the sum of the squared canonical
    correlations.
*/
void civreghdfe_compute_redundant(
    const ST_double *X_endog,
    const ST_double *Z,
    ST_int N,
    const ST_double *weights,
    ST_int weight_type,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    const civ_omega *om,
    const ST_int *redund_indices,
    ST_int n_redund,
    ST_double *redund_stat,
    ST_int *redund_df
)
{
    const ST_double *w = (weights && weight_type != 0) ? weights : NULL;
    ST_int K_rest = K_iv - n_redund;
    ST_int D = K_endog * n_redund;

    *redund_stat = 0.0;
    *redund_df = D;

    if (n_redund <= 0 || K_endog <= 0 || K_rest < K_exog) return;

    ST_int *test = (ST_int *)calloc((size_t)K_iv, sizeof(ST_int));
    ST_double *Z_rest = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)(K_rest > 0 ? K_rest : 1), sizeof(ST_double));
    ST_double *Z_test = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)n_redund, sizeof(ST_double));
    ST_double *Y = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_endog, sizeof(ST_double));
    ST_double *g = (ST_double *)calloc((size_t)D * 2, sizeof(ST_double));
    ST_double *S = (ST_double *)calloc((size_t)D * D, sizeof(ST_double));
    if (!test || !Z_rest || !Z_test || !Y || !g || !S) goto redund_cleanup;

    /* Tested columns: distinct, in-range excluded instruments (mapped in
       civreghdfe_impl.c); refuse anything else */
    ST_int n_masked = 0;
    for (ST_int i = 0; i < n_redund; i++) {
        ST_int idx = redund_indices[i] - 1;
        if (idx >= 0 && idx < K_iv - K_exog && !test[K_exog + idx]) {
            test[K_exog + idx] = 1;
            n_masked++;
        }
    }
    if (n_masked != n_redund) goto redund_cleanup;

    ST_int ir = 0, it = 0;
    for (ST_int k = 0; k < K_iv; k++) {
        ST_double *dst = test[k] ? Z_test + (size_t)(it++) * N : Z_rest + (size_t)(ir++) * N;
        memcpy(dst, Z + (size_t)k * N, (size_t)N * sizeof(ST_double));
    }
    if (civ_partial_cols(X_endog, K_endog, Z_rest, K_rest, w, N, Y) != 0 ||
        civ_partial_cols(Z_test, n_redund, Z_rest, K_rest, w, N, Z_test) != 0)
        goto redund_cleanup;
    civ_standardize_cols(Y, K_endog, N);
    civ_standardize_cols(Z_test, n_redund, N);

    /* g = vec(Z_t'WY), index k*n_redund + l, as the moments e_k z_l */
    for (ST_int k = 0; k < K_endog; k++)
        civ_zwe(Z_test, n_redund, Y + (size_t)k * N, w, N, g + (size_t)k * n_redund);
    if (civ_omega_build(om, Y, K_endog, Z_test, n_redund, S) != 0) goto redund_cleanup;
    if (civ_sym_solve(S, D, g, 1, g + D) != 0) goto redund_cleanup;
    ST_double q = 0.0;
    for (ST_int d = 0; d < D; d++) q += g[d] * g[D + d];
    *redund_stat = q;

redund_cleanup:
    free(test); free(Z_rest); free(Z_test); free(Y); free(g); free(S);
}

/*
    civreghdfe_omega.c
    Covariance matrix of the IV moment conditions (port of ivreg2's m_omega)

    S sums the moment contributions u_i = e_i (x) z_i of the estimation
    sample with the structure of the requested VCE:
      - homoskedastic:  (E'WE/N) (x) Z'WZ; with a kernel (AC) plus, for each
        lag tau, kw(tau) * (sig_tau (x) ZZ_tau + its transpose), where
        sig_tau = sum w_t w_s e_t e_s'/N and ZZ_tau = sum w_t w_s z_t z_s'
        over the pairs s = t - tau
      - heteroskedastic: sum c_i u_i u_i' (c = w^2 for aw/pw, w for fw); with
        a kernel (HAC) plus kw(tau) * sum w_t w_s (u_t u_s' + u_s u_t')
      - cluster: sum_g h_g h_g', h_g = sum_{i in g} w_i u_i; two-way:
        S(c1) + S(c2) - S(c1 x c2); two-way with a kernel (c1 = panel,
        c2 = time): S(c1) + Driscoll-Kraay(c2) - HAC, the HAC S taking the
        place of the intersection
      - Driscoll-Kraay: h_p summed per time period, sum_p h_p h_p' plus
        kw(tau) * sum_p (h_p h_{p-tau}' + h_{p-tau} h_p')
    Kernel lags are actual time differences (t - s)/tdelta within a panel,
    as ivreg2 pairs t with L(tau).t: a gap or a dropped observation breaks
    the pair and an empty period contributes nothing. Lag windows use lags
    1..bw, spectral windows every lag up to T/tdelta - 1 (T = time span).
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "civreghdfe_omega.h"
#include "civreghdfe_matrix.h"
#include "civreghdfe_vce.h"
#include "../ctools_config.h"
#include "../ctools_ols.h"

/* Lags summed by a kernel: bw for lag windows; every lag up to
   T/tdelta - 1 for the spectral Quadratic Spectral window (m_omega's TAU) */
static ST_int omega_max_lag(ST_int kernel, ST_int bw, ST_int tmax_idx, ST_double tdelta)
{
    if (kernel <= 0) return 0;
    if (kernel == CIVREGHDFE_KERNEL_QS) {
        if (tdelta <= 0.0) tdelta = 1.0;
        ST_double tau = floor((ST_double)tmax_idx - 1.0 + 1.0 / tdelta + 1e-9);
        return tau > 0.0 ? (ST_int)tau : 0;
    }
    return bw > 0 ? bw : 0;
}

typedef struct { ST_int a, b, i; } omega_key;

static int cmp_key(const void *x, const void *y)
{
    const omega_key *p = (const omega_key *)x, *q = (const omega_key *)y;
    if (p->a != q->a) return p->a < q->a ? -1 : 1;
    if (p->b != q->b) return p->b < q->b ? -1 : 1;
    return (p->i > q->i) - (p->i < q->i);
}

ST_retcode civ_omega_init(civ_omega *om, ST_int kind, ST_int kernel, ST_int bw,
                          ST_int N, ST_double N_eff, const ST_double *w, ST_int wtype,
                          const ST_int *cl1, ST_int ncl1,
                          const ST_int *cl2, ST_int ncl2,
                          const civreghdfe_ts *ts)
{
    memset(om, 0, sizeof(*om));
    om->N = N;
    om->N_eff = N_eff;
    om->w = (w && wtype != 0) ? w : NULL;
    om->wtype = om->w ? wtype : 0;
    om->kind = kind;
    om->kernel = (kind != CIV_OMEGA_CLUSTER) ? kernel : 0;
    om->bw = bw;
    om->cl1 = cl1; om->ncl1 = ncl1;
    om->cl2 = cl2; om->ncl2 = ncl2;
    om->owns = 1;

    if (om->kernel > 0 || kind == CIV_OMEGA_DKRAAY) {
        if (!ts || !ts->tidx) return 198;
        om->tidx = ts->tidx;
        om->tmax_idx = ts->tmax_idx;
        om->tdelta = ts->tdelta;
        om->max_lag = omega_max_lag(om->kernel, bw, ts->tmax_idx, ts->tdelta);
    }

    /* Two-way: intersection clusters by sorting the (c1, c2) pairs (with a
       kernel the HAC S takes their place) */
    if (kind == CIV_OMEGA_TWOWAY && om->kernel <= 0) {
        omega_key *key = (omega_key *)ctools_safe_malloc2((size_t)N, sizeof(omega_key));
        om->cl12 = (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int));
        if (!key || !om->cl12) { free(key); civ_omega_free(om); return 920; }
        for (ST_int i = 0; i < N; i++) { key[i].a = cl1[i]; key[i].b = cl2[i]; key[i].i = i; }
        qsort(key, (size_t)N, sizeof(omega_key), cmp_key);
        ST_int g = 0;
        for (ST_int j = 0; j < N; j++) {
            if (j == 0 || key[j].a != key[j - 1].a || key[j].b != key[j - 1].b) g++;
            om->cl12[key[j].i] = g;
        }
        om->ncl12 = g;
        free(key);
    }

    /* Kernel lags within panels: observations sorted by (panel, time) */
    if (om->kernel > 0 && (kind == CIV_OMEGA_IID || kind == CIV_OMEGA_ROBUST || kind == CIV_OMEGA_TWOWAY)) {
        ST_int P = (ts->panel && ts->num_panels > 0) ? ts->num_panels : 1;
        omega_key *key = (omega_key *)ctools_safe_malloc2((size_t)N, sizeof(omega_key));
        om->order = (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int));
        om->pstart = (ST_int *)calloc((size_t)P + 1, sizeof(ST_int));
        if (!key || !om->order || !om->pstart) { free(key); civ_omega_free(om); return 920; }
        for (ST_int i = 0; i < N; i++) {
            key[i].a = ts->panel ? ts->panel[i] - 1 : 0;
            key[i].b = ts->tidx[i];
            key[i].i = i;
        }
        qsort(key, (size_t)N, sizeof(omega_key), cmp_key);
        for (ST_int j = 0; j < N; j++) {
            om->order[j] = key[j].i;
            if (key[j].a >= 0 && key[j].a < P) om->pstart[key[j].a + 1]++;
        }
        for (ST_int p = 0; p < P; p++) om->pstart[p + 1] += om->pstart[p];
        om->npanels = P;
        free(key);
    }
    return 0;
}

civ_omega civ_omega_with_kernel(const civ_omega *om, ST_int kernel)
{
    civ_omega c = *om;
    c.owns = 0;
    if (c.kind != CIV_OMEGA_CLUSTER) {
        /* A copy can switch kernels, not add one: the lag structure is built
           by civ_omega_init only when the original has a kernel. */
        if (kernel > 0 && om->kernel <= 0) kernel = 0;
        c.kernel = kernel;
        c.max_lag = omega_max_lag(kernel, c.bw, c.tmax_idx, c.tdelta);
    }
    return c;
}

void civ_omega_free(civ_omega *om)
{
    if (!om || !om->owns) return;
    free(om->cl12); om->cl12 = NULL;
    free(om->order); om->order = NULL;
    free(om->pstart); om->pstart = NULL;
}

/* S += scale * (A (x) B) (trans = 0) or its transpose (trans = 1); A is K x K,
   B is L x L, S is (K*L) x (K*L), all column-major */
static void add_kron(ST_double *S, const ST_double *A, ST_int K, const ST_double *B, ST_int L,
                     ST_double scale, int trans)
{
    ST_int D = K * L;
    for (ST_int k2 = 0; k2 < K; k2++)
        for (ST_int l2 = 0; l2 < L; l2++)
            for (ST_int k1 = 0; k1 < K; k1++)
                for (ST_int l1 = 0; l1 < L; l1++) {
                    ST_double v = trans ? A[k1 * K + k2] * B[l1 * L + l2]
                                        : A[k2 * K + k1] * B[l2 * L + l1];
                    S[(size_t)(k2 * L + l2) * D + (k1 * L + l1)] += scale * v;
                }
}

/* u = e_i (x) z_i (length K*L) */
static inline void moment_row(const ST_double *E, ST_int K, const ST_double *Z, ST_int L,
                              ST_int N, ST_int i, ST_double *u)
{
    for (ST_int k = 0; k < K; k++) {
        ST_double e = E[(size_t)k * N + i];
        for (ST_int l = 0; l < L; l++) u[k * L + l] = e * Z[(size_t)l * N + i];
    }
}

/* Iterator over the pairs (now, lag) of observations in the same panel
   whose time indices differ by exactly tau (a two-pointer scan of each
   panel's time-sorted observations) */
typedef struct {
    const civ_omega *om;
    const ST_int *ob;
    ST_int tau, p, m, a, b;
} pair_iter;

static void pair_begin(pair_iter *it, const civ_omega *om, ST_int tau)
{
    it->om = om; it->ob = NULL; it->tau = tau;
    it->p = -1; it->m = 0; it->a = 0; it->b = 0;
}

static int pair_next(pair_iter *it, ST_int *now, ST_int *lagobs)
{
    const civ_omega *om = it->om;
    const ST_int *t = om->tidx;
    for (;;) {
        while (it->b >= it->m) {
            if (++it->p >= om->npanels) return 0;
            it->ob = om->order + om->pstart[it->p];
            it->m = om->pstart[it->p + 1] - om->pstart[it->p];
            it->a = 0; it->b = 1;
            /* skip panels spanning fewer than tau periods */
            if (it->m < 2 || t[it->ob[it->m - 1]] - t[it->ob[0]] < it->tau) it->m = 0;
        }
        ST_int b = it->b++;
        ST_int target = t[it->ob[b]] - it->tau;
        if (target < t[it->ob[0]]) continue;
        while (t[it->ob[it->a]] < target) it->a++;
        if (t[it->ob[it->a]] != target) continue;
        *now = it->ob[b];
        *lagobs = it->ob[it->a];
        return 1;
    }
}

/* Sum over groups of h_g h_g', h_g = sum_{i in g} w_i u_i (1-based ids);
   sign +1 or -1. Used for clusters and Driscoll-Kraay periods (with the
   kernel lags when kernel > 0). */
static ST_retcode omega_groups(const civ_omega *om, const ST_double *E, ST_int K,
                               const ST_double *Z, ST_int L, const ST_int *gid, ST_int G,
                               ST_int use_kernel, ST_double sign, ST_double *S)
{
    ST_int D = K * L, N = om->N;
    ST_double *h = (ST_double *)ctools_safe_calloc3((size_t)G, (size_t)D, sizeof(ST_double));
    ST_double *u = (ST_double *)malloc((size_t)D * sizeof(ST_double));
    if (!h || !u) { free(h); free(u); return 920; }

    for (ST_int i = 0; i < N; i++) {
        ST_int g = gid ? gid[i] - 1 : om->tidx[i];
        if (g < 0 || g >= G) continue;
        ST_double wi = om->w ? om->w[i] : 1.0;
        moment_row(E, K, Z, L, N, i, u);
        ST_double *hg = h + (size_t)g * D;
        for (ST_int d = 0; d < D; d++) hg[d] += wi * u[d];
    }
    free(u);

    ST_int maxlag = use_kernel ? om->max_lag : 0;
    if (maxlag >= G) maxlag = G - 1;
    for (ST_int tau = 0; tau <= maxlag; tau++) {
        ST_double kw = civreghdfe_kernel_weight(om->kernel, tau, om->bw);
        if (kw == 0.0) continue;
        ST_double f = sign * kw;
        for (ST_int g = tau; g < G; g++) {
            const ST_double *a = h + (size_t)g * D;
            const ST_double *b = h + (size_t)(g - tau) * D;
            for (ST_int c = 0; c < D; c++) {
                ST_double *Sc = S + (size_t)c * D;
                if (tau == 0) {
                    for (ST_int r = 0; r < D; r++) Sc[r] += f * a[r] * a[c];
                } else {
                    for (ST_int r = 0; r < D; r++) Sc[r] += f * (a[r] * b[c] + b[r] * a[c]);
                }
            }
        }
    }
    free(h);
    return 0;
}

/* Heteroskedastic S (HC; with a kernel HAC), added to S with the given
   sign: lag 0 = sum c_i u_i u_i', c = w^2 (aw/pw) or w (fw), and for each
   kernel lag kw(tau) * sum w_t w_s (u_t u_s' + u_s u_t') */
static ST_retcode omega_hetero(const civ_omega *om, const ST_double *E, ST_int K,
                               const ST_double *Z, ST_int L, ST_double sign, ST_double *S)
{
    ST_int N = om->N, D = K * L;
    const ST_double *w = om->w;
    size_t DD = (size_t)D * D;
    int nt = ctools_get_openmp_threads();
    if (nt < 1) nt = 1;
    ST_int maxlag = om->kernel > 0 ? om->max_lag : 0;

    ST_double *buf = (ST_double *)ctools_safe_calloc3((size_t)nt, DD + 2 * (size_t)D, sizeof(ST_double));
    if (!buf) return 920;
    size_t per = DD + 2 * (size_t)D;
    #pragma omp parallel num_threads(nt) if (N > 5000 && nt > 1)
    {
        int t = 0;
#ifdef _OPENMP
        t = omp_get_thread_num();
#endif
        ST_double *Sl = buf + (size_t)t * per, *u = Sl + DD;
        #pragma omp for schedule(static)
        for (ST_int i = 0; i < N; i++) {
            ST_double c = 1.0;
            if (w) c = (om->wtype == 2) ? w[i] : w[i] * w[i];
            moment_row(E, K, Z, L, N, i, u);
            for (ST_int b = 0; b < D; b++) {
                ST_double cb = c * u[b];
                ST_double *Sb = Sl + (size_t)b * D;
                for (ST_int a = 0; a < D; a++) Sb[a] += u[a] * cb;
            }
        }
    }

    /* HAC lags: kw(tau) * sum w_t w_s (u_t u_s' + u_s u_t') */
    if (maxlag > 0) {
        #pragma omp parallel for schedule(dynamic) num_threads(nt) if (maxlag > 1 && nt > 1)
        for (ST_int tau = 1; tau <= maxlag; tau++) {
            int t = 0;
#ifdef _OPENMP
            t = omp_get_thread_num();
#endif
            ST_double kw = civreghdfe_kernel_weight(om->kernel, tau, om->bw);
            if (kw == 0.0) continue;
            ST_double *Sl = buf + (size_t)t * per, *u = Sl + DD, *v = u + D;
            ST_int now, lagobs;
            pair_iter it;
            pair_begin(&it, om, tau);
            while (pair_next(&it, &now, &lagobs)) {
                ST_double wv = kw * (w ? w[now] * w[lagobs] : 1.0);
                moment_row(E, K, Z, L, N, now, u);
                moment_row(E, K, Z, L, N, lagobs, v);
                for (ST_int b = 0; b < D; b++) {
                    ST_double *Sb = Sl + (size_t)b * D;
                    ST_double vb = wv * v[b], ub = wv * u[b];
                    for (ST_int a = 0; a < D; a++) Sb[a] += u[a] * vb + v[a] * ub;
                }
            }
        }
    }
    for (int t = 0; t < nt; t++) {
        const ST_double *Sl = buf + (size_t)t * per;
        for (size_t j = 0; j < DD; j++) S[j] += sign * Sl[j];
    }
    free(buf);
    return 0;
}

ST_retcode civ_omega_build(const civ_omega *om, const ST_double *E, ST_int K,
                           const ST_double *Z, ST_int L, ST_double *S)
{
    ST_int N = om->N, D = K * L;
    const ST_double *w = om->w;
    size_t DD = (size_t)D * D;
    memset(S, 0, DD * sizeof(ST_double));

    if (om->kind == CIV_OMEGA_CLUSTER)
        return omega_groups(om, E, K, Z, L, om->cl1, om->ncl1, 0, 1.0, S);
    if (om->kind == CIV_OMEGA_TWOWAY) {
        ST_retcode rc = omega_groups(om, E, K, Z, L, om->cl1, om->ncl1, 0, 1.0, S);
        if (om->kernel > 0) {
            /* ivreg2's two-way clustering on panel and time with a kernel:
               panel clusters + Driscoll-Kraay over the periods - HAC */
            if (!rc) rc = omega_groups(om, E, K, Z, L, NULL, om->tmax_idx + 1, 1, 1.0, S);
            if (!rc) rc = omega_hetero(om, E, K, Z, L, -1.0, S);
        } else {
            if (!rc) rc = omega_groups(om, E, K, Z, L, om->cl2, om->ncl2, 0, 1.0, S);
            if (!rc) rc = omega_groups(om, E, K, Z, L, om->cl12, om->ncl12, 0, -1.0, S);
        }
        return rc;
    }
    if (om->kind == CIV_OMEGA_DKRAAY)
        return omega_groups(om, E, K, Z, L, NULL, om->tmax_idx + 1, om->kernel > 0, 1.0, S);

    int nt = ctools_get_openmp_threads();
    if (nt < 1) nt = 1;
    ST_int maxlag = om->kernel > 0 ? om->max_lag : 0;

    if (om->kind == CIV_OMEGA_IID) {
        /* sigma (x) Z'WZ, sigma = E'WE / N */
        ST_double *sig = (ST_double *)calloc((size_t)K * K, sizeof(ST_double));
        ST_double *ZZ = (ST_double *)calloc((size_t)L * L, sizeof(ST_double));
        if (!sig || !ZZ) { free(sig); free(ZZ); return 920; }
        if (w) {
            ctools_matmul_atdb(E, E, w, N, K, K, sig);
            ctools_matmul_atdb(Z, Z, w, N, L, L, ZZ);
        } else {
            ctools_matmul_atb(E, E, N, K, K, sig);
            ctools_matmul_atb(Z, Z, N, L, L, ZZ);
        }
        for (ST_int j = 0; j < K * K; j++) sig[j] /= om->N_eff;
        add_kron(S, sig, K, ZZ, L, 1.0, 0);
        free(sig); free(ZZ);
        if (maxlag <= 0) return 0;

        /* AC lags: per-lag sig_tau (x) ZZ_tau; one lag per iteration */
        size_t per = DD + (size_t)K * K + (size_t)L * L;
        ST_double *buf = (ST_double *)ctools_safe_calloc3((size_t)nt, per, sizeof(ST_double));
        if (!buf) return 920;
        #pragma omp parallel for schedule(dynamic) num_threads(nt) if (maxlag > 1 && nt > 1)
        for (ST_int tau = 1; tau <= maxlag; tau++) {
            int t = 0;
#ifdef _OPENMP
            t = omp_get_thread_num();
#endif
            ST_double kw = civreghdfe_kernel_weight(om->kernel, tau, om->bw);
            if (kw == 0.0) continue;
            ST_double *Sl = buf + (size_t)t * per;
            ST_double *st = Sl + DD, *zt = st + (size_t)K * K;
            memset(st, 0, ((size_t)K * K + (size_t)L * L) * sizeof(ST_double));
            ST_int npairs = 0, now, lagobs;
            pair_iter it;
            pair_begin(&it, om, tau);
            while (pair_next(&it, &now, &lagobs)) {
                ST_double wv = w ? w[now] * w[lagobs] : 1.0;
                for (ST_int k2 = 0; k2 < K; k2++)
                    for (ST_int k1 = 0; k1 < K; k1++)
                        st[k2 * K + k1] += wv * E[(size_t)k1 * N + now] * E[(size_t)k2 * N + lagobs];
                for (ST_int l2 = 0; l2 < L; l2++)
                    for (ST_int l1 = 0; l1 < L; l1++)
                        zt[l2 * L + l1] += wv * Z[(size_t)l1 * N + now] * Z[(size_t)l2 * N + lagobs];
                npairs++;
            }
            if (npairs == 0) continue;
            for (ST_int j = 0; j < K * K; j++) st[j] /= om->N_eff;
            add_kron(Sl, st, K, zt, L, kw, 0);
            add_kron(Sl, st, K, zt, L, kw, 1);
        }
        for (int t = 0; t < nt; t++) {
            const ST_double *Sl = buf + (size_t)t * per;
            for (size_t j = 0; j < DD; j++) S[j] += Sl[j];
        }
        free(buf);
        return 0;
    }

    return omega_hetero(om, E, K, Z, L, 1.0, S);
}

ST_int civ_sym_inverse(const ST_double *A, ST_int n, ST_double *inv)
{
    if (n <= 0) return 0;
    ST_double *L = (ST_double *)ctools_safe_malloc3((size_t)n, (size_t)n, sizeof(ST_double));
    if (!L) return -1;
    memcpy(L, A, (size_t)n * n * sizeof(ST_double));
    ST_int ok = (ctools_cholesky(L, n) == 0) && ctools_invert_from_cholesky(L, n, inv) == 0;
    /* The Cholesky inverse is Stata's invsym result when no column is
       nearly dependent on the others: its variance conditional on all of
       them, 1/inv[j,j], at least 1e-9 of its diagonal and at least 1e-19
       (invsym's tolerances; the conditional variance at any pivot step can
       only be larger). Otherwise use the sweep inverse, which zeroes the
       dependent columns as invsym does. */
    for (ST_int j = 0; j < n && ok; j++) {
        ST_double ijj = inv[(size_t)j * n + j];
        if (!(ijj > 0.0) || 1.0 / ijj < 1e-9 * A[(size_t)j * n + j] || 1.0 / ijj < 1e-19) ok = 0;
    }
    free(L);
    if (!ok) return ctools_invsym(A, n, inv, NULL) == 0 ? 0 : -1;
    return 0;
}

ST_int civ_sym_rank(const ST_double *A, ST_int n, ST_double scale)
{
    if (n <= 0) return 0;
    ST_double *B = (ST_double *)ctools_safe_malloc3((size_t)n, (size_t)n, sizeof(ST_double));
    ST_double *inv = (ST_double *)ctools_safe_malloc3((size_t)n, (size_t)n, sizeof(ST_double));
    ST_int rank = -1;
    if (B && inv) {
        for (size_t j = 0; j < (size_t)n * n; j++) B[j] = scale * A[j];
        if (ctools_invsym(B, n, inv, &rank) != 0) rank = -1;
    }
    free(B); free(inv);
    return rank;
}

ST_int civ_sym_solve(const ST_double *A, ST_int n, const ST_double *b, ST_int nrhs, ST_double *x)
{
    ST_double *inv = (ST_double *)ctools_safe_malloc3((size_t)n, (size_t)n, sizeof(ST_double));
    ST_double *t = (ST_double *)ctools_safe_malloc3((size_t)n, (size_t)(nrhs > 0 ? nrhs : 1), sizeof(ST_double));
    if (!inv || !t || civ_sym_inverse(A, n, inv) != 0) { free(inv); free(t); return -1; }
    for (ST_int c = 0; c < nrhs; c++)
        for (ST_int i = 0; i < n; i++) {
            ST_double s = 0.0;
            for (ST_int j = 0; j < n; j++) s += inv[(size_t)j * n + i] * b[(size_t)c * n + j];
            t[(size_t)c * n + i] = s;
        }
    memcpy(x, t, (size_t)n * nrhs * sizeof(ST_double));
    free(inv); free(t);
    return 0;
}

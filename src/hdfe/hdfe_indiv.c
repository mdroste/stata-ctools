/*
 * hdfe_indiv.c
 *
 * Individual fixed effects in a group setting (see hdfe_indiv.h)
 * Part of the ctools Stata plugin suite
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "hdfe_indiv.h"
#include "../ctools_cancel.h"
#include "../ctools_alloc.h"
#include "../ctools_data.h"
#include "../ctools_data_io.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void hdfe_indiv_free(HDFE_Indiv *I)
{
    if (!I) return;
    free(I->ptr);
    free(I->ind);
    I->ptr = I->ind = NULL;
    I->N_rows = I->NI = 0;
}

ST_retcode hdfe_indiv_load(HDFE_Indiv *I, const ST_double *row_group, ST_int N_rows, int gvar,
                           int ivar, int tvar, int mean)
{
    ST_retcode rc = 0;
    ctools_filtered_data e;
    ST_int *row_of = NULL, *count = NULL;
    ST_int max_g = 0, max_i = 0;
    memset(I, 0, sizeof(*I));
    I->mean = mean;
    I->N_rows = N_rows;
    ctools_filtered_data_init(&e);

    for (ST_int r = 0; r < N_rows; r++) {
        ST_double v = row_group[r];
        if (!(v >= 1.0) || v > 2147483646.0 || v != floor(v)) return 198;
        if ((ST_int)v > max_g) max_g = (ST_int)v;
    }
    int vars[3] = {gvar, ivar, tvar};
    stata_retcode lrc = ctools_data_load(&e, vars, 3, 0, 0, CTOOLS_LOAD_SKIP_IF | CTOOLS_LOAD_NO_SORT_ORDER);
    if (lrc != STATA_OK) {
        ctools_filtered_data_free(&e);
        return ctools_stata_rc(lrc);
    }
    const size_t n = e.data.nobs;
    const ST_double *gv = e.data.vars[0].data.dbl, *iv = e.data.vars[1].data.dbl,
                    *tv = e.data.vars[2].data.dbl;
    if (n > 0 && (!gv || !iv || !tv)) { rc = 920; goto done; }

    /* The regression row of each group number */
    row_of = (ST_int *)malloc((size_t)(max_g + 1) * sizeof(ST_int));
    count = (ST_int *)calloc((size_t)N_rows + 1, sizeof(ST_int));
    if (!row_of || !count) { rc = 920; goto done; }
    for (ST_int g = 0; g <= max_g; g++) row_of[g] = -1;
    for (ST_int r = 0; r < N_rows; r++) {
        ST_int g = (ST_int)row_group[r];
        if (row_of[g] >= 0) { rc = 198; goto done; }
        row_of[g] = r;
    }

    /* Edges: the marked rows whose group has a regression row */
    for (size_t j = 0; j < n; j++) {
        if (!(tv[j] != 0.0) || !(tv[j] < SV_missval)) continue;
        ST_double gd = gv[j], id = iv[j];
        if (!(gd >= 1.0) || !(id >= 1.0) || !(gd < SV_missval) || !(id < SV_missval)) continue;
        if (gd > (ST_double)max_g) continue;
        ST_int r = row_of[(ST_int)gd];
        if (r < 0) continue;
        if (id > 2147483646.0) { rc = 198; goto done; }
        if ((ST_int)id > max_i) max_i = (ST_int)id;
        count[r + 1]++;
    }
    I->ptr = (ST_int *)malloc(((size_t)N_rows + 1) * sizeof(ST_int));
    if (!I->ptr) { rc = 920; goto done; }
    I->ptr[0] = 0;
    for (ST_int r = 0; r < N_rows; r++) I->ptr[r + 1] = I->ptr[r] + count[r + 1];
    I->ind = (ST_int *)malloc((size_t)(I->ptr[N_rows] > 0 ? I->ptr[N_rows] : 1) * sizeof(ST_int));
    if (!I->ind) { rc = 920; goto done; }
    for (ST_int r = 0; r < N_rows; r++) count[r] = I->ptr[r];
    for (size_t j = 0; j < n; j++) {
        if (!(tv[j] != 0.0) || !(tv[j] < SV_missval)) continue;
        ST_double gd = gv[j], id = iv[j];
        if (!(gd >= 1.0) || !(id >= 1.0) || !(gd < SV_missval) || !(id < SV_missval)) continue;
        if (gd > (ST_double)max_g) continue;
        ST_int r = row_of[(ST_int)gd];
        if (r < 0) continue;
        I->ind[count[r]++] = (ST_int)id - 1;
    }
    I->NI = max_i;
done:
    ctools_filtered_data_free(&e);
    free(row_of);
    free(count);
    if (rc) hdfe_indiv_free(I);
    return rc;
}

/* A stack of rows to drop (a row may be pushed more than once) */
typedef struct {
    ST_int *v;
    size_t n, cap;
} row_stack;

static int stack_push(row_stack *s, ST_int r)
{
    if (s->n == s->cap) {
        size_t cap = s->cap ? 2 * s->cap : 1024;
        ST_int *v = (ST_int *)realloc(s->v, cap * sizeof(ST_int));
        if (!v) return -1;
        s->v = v;
        s->cap = cap;
    }
    s->v[s->n++] = r;
    return 0;
}

ST_int hdfe_indiv_singletons(ST_int *const *levels, ST_int G, const HDFE_Indiv *I, ST_int *mask)
{
    const ST_int N = I->N_rows, NI = I->NI;
    ST_int dropped = -1;
    ST_int **cnt = NULL, **rx = NULL, *deg = NULL, *ix = NULL;
    row_stack st = {0};
    for (ST_int r = 0; r < N; r++) mask[r] = 1;
    if (N <= 0) return 0;

    cnt = (ST_int **)calloc((size_t)(G > 0 ? G : 1), sizeof(ST_int *));
    rx = (ST_int **)calloc((size_t)(G > 0 ? G : 1), sizeof(ST_int *));
    deg = (ST_int *)calloc((size_t)(NI > 0 ? NI : 1), sizeof(ST_int));
    ix = (ST_int *)calloc((size_t)(NI > 0 ? NI : 1), sizeof(ST_int));
    if (!cnt || !rx || !deg || !ix) goto cleanup;
    for (ST_int g = 0; g < G; g++) {
        ST_int L = 0;
        for (ST_int r = 0; r < N; r++) {
            if (levels[g][r] < 1) goto cleanup;
            if (levels[g][r] > L) L = levels[g][r];
        }
        cnt[g] = (ST_int *)calloc((size_t)L + 1, sizeof(ST_int));
        rx[g] = (ST_int *)calloc((size_t)L + 1, sizeof(ST_int));
        if (!cnt[g] || !rx[g]) goto cleanup;
        for (ST_int r = 0; r < N; r++) {
            cnt[g][levels[g][r]]++;
            rx[g][levels[g][r]] ^= r;
        }
        for (ST_int l = 1; l <= L; l++)
            if (cnt[g][l] == 1 && stack_push(&st, rx[g][l])) goto cleanup;
    }
    /* Each individual's remaining groups: their count and the xor of their rows */
    for (ST_int r = 0; r < N; r++) {
        for (ST_int e = I->ptr[r]; e < I->ptr[r + 1]; e++) {
            deg[I->ind[e]]++;
            ix[I->ind[e]] ^= r;
        }
    }
    for (ST_int i = 0; i < NI; i++)
        if (deg[i] == 1 && stack_push(&st, ix[i])) goto cleanup;

    dropped = 0;
    while (st.n > 0) {
        ST_int r = st.v[--st.n];
        if (!mask[r]) continue;
        mask[r] = 0;
        dropped++;
        for (ST_int g = 0; g < G; g++) {
            ST_int l = levels[g][r];
            cnt[g][l]--;
            rx[g][l] ^= r;
            if (cnt[g][l] == 1 && stack_push(&st, rx[g][l])) { dropped = -1; goto cleanup; }
        }
        for (ST_int e = I->ptr[r]; e < I->ptr[r + 1]; e++) {
            ST_int i = I->ind[e];
            deg[i]--;
            ix[i] ^= r;
            if (deg[i] == 1 && stack_push(&st, ix[i])) { dropped = -1; goto cleanup; }
        }
    }
cleanup:
    for (ST_int g = 0; cnt && g < G; g++) free(cnt[g]);
    for (ST_int g = 0; rx && g < G; g++) free(rx[g]);
    free(cnt); free(rx); free(deg); free(ix); free(st.v);
    return dropped;
}

ST_retcode hdfe_indiv_compact(HDFE_Indiv *I, const ST_int *mask)
{
    const ST_int N = I->N_rows;
    ST_int *num = (ST_int *)malloc((size_t)(I->NI > 0 ? I->NI : 1) * sizeof(ST_int));
    if (!num) return 920;
    for (ST_int i = 0; i < I->NI; i++) num[i] = -1;
    ST_int nr = 0, ne = 0, ni = 0;
    for (ST_int r = 0; r < N; r++) {
        if (!mask[r]) continue;
        ST_int start = I->ptr[r], end = I->ptr[r + 1];
        I->ptr[nr++] = ne;
        for (ST_int e = start; e < end; e++) {
            ST_int i = I->ind[e];
            if (num[i] < 0) num[i] = ni++;
            I->ind[ne++] = num[i];
        }
    }
    I->ptr[nr] = ne;
    I->N_rows = nr;
    I->NI = ni;
    free(num);
    return 0;
}

/* The preconditioned, row-weighted design: A D with rows scaled by the
 * square roots of the weights. Columns: the levels of each factor, then the
 * individuals. */
typedef struct {
    const HDFE_State *S;
    const HDFE_Indiv *I;
    ST_int m, n, G;
    const ST_int *off;    /* G + 1 column offsets of the factors; individuals from off[G] */
    const ST_double *d;   /* n column scalings (0 for an empty column) */
    const ST_double *s;   /* m row scalings sqrt(w), or NULL */
    const ST_double *agg; /* m: 1/#individuals (mean) or 1 (sum) */
} indiv_op;

/* out = A D v */
static void op_mult(const indiv_op *op, const ST_double *v, ST_double *out)
{
    const ST_int m = op->m, G = op->G, base = op->off[G];
    for (ST_int r = 0; r < m; r++) {
        ST_double acc = 0.0;
        for (ST_int g = 0; g < G; g++) {
            ST_int c = op->off[g] + op->S->factors[g].levels[r] - 1;
            acc += op->d[c] * v[c];
        }
        ST_double t = 0.0;
        for (ST_int e = op->I->ptr[r]; e < op->I->ptr[r + 1]; e++) {
            ST_int c = base + op->I->ind[e];
            t += op->d[c] * v[c];
        }
        acc += op->agg[r] * t;
        out[r] = op->s ? op->s[r] * acc : acc;
    }
}

/* out = D A' u */
static void op_tmult(const indiv_op *op, const ST_double *u, ST_double *out)
{
    const ST_int m = op->m, G = op->G, base = op->off[G];
    memset(out, 0, (size_t)op->n * sizeof(ST_double));
    for (ST_int r = 0; r < m; r++) {
        ST_double ur = op->s ? op->s[r] * u[r] : u[r];
        for (ST_int g = 0; g < G; g++) out[op->off[g] + op->S->factors[g].levels[r] - 1] += ur;
        ST_double ua = op->agg[r] * ur;
        for (ST_int e = op->I->ptr[r]; e < op->I->ptr[r + 1]; e++) out[base + op->I->ind[e]] += ua;
    }
    for (ST_int c = 0; c < op->n; c++) out[c] *= op->d[c];
}

static ST_double vnorm(const ST_double *x, ST_int n)
{
    ST_double s = 0.0;
    for (ST_int i = 0; i < n; i++) s += x[i] * x[i];
    return sqrt(s);
}

/* Stable plane rotation (scipy's _sym_ortho): c, s, r with r = hypot(a, b) */
static void sym_ortho(ST_double a, ST_double b, ST_double *c, ST_double *s, ST_double *r)
{
    if (b == 0.0) {
        *c = (a > 0.0) - (a < 0.0);
        if (a == 0.0) *c = 1.0;
        *s = 0.0;
        *r = fabs(a);
    } else if (a == 0.0) {
        *c = 0.0;
        *s = (b > 0.0) - (b < 0.0);
        *r = fabs(b);
    } else if (fabs(b) > fabs(a)) {
        ST_double tau = a / b;
        *s = ((b > 0.0) - (b < 0.0)) / sqrt(1.0 + tau * tau);
        *c = *s * tau;
        *r = b / *s;
    } else {
        ST_double tau = b / a;
        *c = ((a > 0.0) - (a < 0.0)) / sqrt(1.0 + tau * tau);
        *s = *c * tau;
        *r = a / *c;
    }
}

/* LSMR (Fong and Saunders 2011, as scipy's lsmr without damping) on the
 * scaled column b (m), leaving the residual b - A D x in b. ws: 5n + m
 * doubles. Returns the iterations, or -1 without convergence. */
static ST_int lsmr_column(const indiv_op *op, ST_double *b, ST_double tol, ST_int maxiter, ST_double *ws)
{
    const ST_int m = op->m, n = op->n;
    ST_double *u = ws, *v = u + m, *h = v + n, *hbar = h + n, *x = hbar + n, *w = x + n;
    memcpy(u, b, (size_t)m * sizeof(ST_double));
    ST_double beta = vnorm(u, m);
    const ST_double normb = beta;
    if (beta == 0.0) return 0;
    for (ST_int i = 0; i < m; i++) u[i] /= beta;
    op_tmult(op, u, v);
    ST_double alpha = vnorm(v, n);
    if (alpha == 0.0) return 0; /* b is orthogonal to the design: its residual is b */
    for (ST_int j = 0; j < n; j++) v[j] /= alpha;

    ST_double zetabar = alpha * beta, alphabar = alpha, rho = 1.0, rhobar = 1.0, cbar = 1.0,
              sbar = 0.0;
    memcpy(h, v, (size_t)n * sizeof(ST_double));
    memset(hbar, 0, (size_t)n * sizeof(ST_double));
    memset(x, 0, (size_t)n * sizeof(ST_double));
    /* estimates of ||r||, ||A|| and cond(A) */
    ST_double betadd = beta, betad = 0.0, rhodold = 1.0, tautildeold = 0.0, thetatilde = 0.0,
              zeta = 0.0, dsum = 0.0;
    ST_double normA2 = alpha * alpha, maxrbar = 0.0, minrbar = 1e100;
    ST_int itn = 0, converged = 0;
    while (itn < maxiter) {
        if (ctools_cancel_requested()) return -1;
        itn++;
        op_mult(op, v, w);
        for (ST_int i = 0; i < m; i++) u[i] = w[i] - alpha * u[i];
        beta = vnorm(u, m);
        if (beta > 0.0) {
            for (ST_int i = 0; i < m; i++) u[i] /= beta;
            op_tmult(op, u, w);
            for (ST_int j = 0; j < n; j++) v[j] = w[j] - beta * v[j];
            alpha = vnorm(v, n);
            if (alpha > 0.0)
                for (ST_int j = 0; j < n; j++) v[j] /= alpha;
        }
        /* no damping: the rotation Qhat is the identity */
        ST_double chat = 1.0, shat = 0.0, alphahat = alphabar;
        ST_double rhoold = rho, c, s;
        sym_ortho(alphahat, beta, &c, &s, &rho);
        ST_double thetanew = s * alpha;
        alphabar = c * alpha;
        ST_double rhobarold = rhobar, zetaold = zeta;
        ST_double thetabar = sbar * rho, rhotemp = cbar * rho;
        sym_ortho(cbar * rho, thetanew, &cbar, &sbar, &rhobar);
        zeta = cbar * zetabar;
        zetabar = -sbar * zetabar;
        ST_double f1 = -(thetabar * rho / (rhoold * rhobarold)), f2 = zeta / (rho * rhobar),
                  f3 = -(thetanew / rho);
        for (ST_int j = 0; j < n; j++) {
            hbar[j] = h[j] + f1 * hbar[j];
            x[j] += f2 * hbar[j];
            h[j] = v[j] + f3 * h[j];
        }
        /* ||r|| */
        ST_double betaacute = chat * betadd, betacheck = -shat * betadd;
        ST_double betahat = c * betaacute;
        betadd = -s * betaacute;
        ST_double thetatildeold = thetatilde, ctildeold, stildeold, rhotildeold;
        sym_ortho(rhodold, thetabar, &ctildeold, &stildeold, &rhotildeold);
        thetatilde = stildeold * rhobar;
        rhodold = ctildeold * rhobar;
        betad = -stildeold * betad + ctildeold * betahat;
        tautildeold = (zetaold - thetatildeold * tautildeold) / rhotildeold;
        ST_double taud = (zeta - thetatilde * tautildeold) / rhodold;
        dsum += betacheck * betacheck;
        ST_double normr = sqrt(dsum + (betad - taud) * (betad - taud) + betadd * betadd);
        /* ||A|| and cond(A) */
        normA2 += beta * beta;
        ST_double normA = sqrt(normA2);
        normA2 += alpha * alpha;
        if (rhobarold > maxrbar) maxrbar = rhobarold;
        if (itn > 1 && rhobarold < minrbar) minrbar = rhobarold;
        (void)rhotemp;
        /* convergence */
        ST_double normar = fabs(zetabar), normx = vnorm(x, n);
        ST_double test1 = normr / normb;
        ST_double test2 = (normA * normr) != 0.0 ? normar / (normA * normr) : INFINITY;
        ST_double rtol = tol + tol * normA * normx / normb;
        if (test2 <= tol || test1 <= rtol || normar == 0.0) {
            converged = 1;
            break;
        }
    }
    /* the residual, computed rather than estimated */
    op_mult(op, x, w);
    for (ST_int i = 0; i < m; i++) b[i] -= w[i];
    return converged ? itn : -1;
}

ST_retcode hdfe_indiv_partial_out(const HDFE_State *S, const HDFE_Indiv *I, ST_double *data, ST_int N,
                                  ST_int K, ST_double tol, ST_int maxiter, ST_int *iterations)
{
    const ST_int G = S->G;
    ST_retcode rc = 0;
    ST_int *off = NULL, nonconv = 0, most = 0;
    ST_double *d = NULL, *sw = NULL, *agg = NULL;
    if (I->N_rows != N || (G > 0 && S->N != N)) return 198;
    off = (ST_int *)malloc((size_t)(G + 1) * sizeof(ST_int));
    if (!off) return 920;
    off[0] = 0;
    for (ST_int g = 0; g < G; g++) off[g + 1] = off[g] + S->factors[g].num_levels;
    const ST_int n = off[G] + I->NI;
    d = (ST_double *)calloc((size_t)(n > 0 ? n : 1), sizeof(ST_double));
    agg = (ST_double *)malloc((size_t)(N > 0 ? N : 1) * sizeof(ST_double));
    if (S->weights) sw = (ST_double *)malloc((size_t)(N > 0 ? N : 1) * sizeof(ST_double));
    if (!d || !agg || (S->weights && !sw)) { rc = 920; goto done; }
    for (ST_int r = 0; r < N; r++) {
        ST_int k = I->ptr[r + 1] - I->ptr[r];
        agg[r] = I->mean ? (k > 0 ? 1.0 / k : 0.0) : 1.0;
        if (sw) sw[r] = sqrt(S->weights[r]);
    }
    /* Diagonal preconditioner: the inverse norm of each column */
    for (ST_int r = 0; r < N; r++) {
        ST_double wr = S->weights ? S->weights[r] : 1.0;
        for (ST_int g = 0; g < G; g++) d[off[g] + S->factors[g].levels[r] - 1] += wr;
        for (ST_int e = I->ptr[r]; e < I->ptr[r + 1]; e++) d[off[G] + I->ind[e]] += wr * agg[r] * agg[r];
    }
    for (ST_int c = 0; c < n; c++) d[c] = d[c] > 0.0 ? 1.0 / sqrt(d[c]) : 0.0;

    indiv_op op = {S, I, N, n, G, off, d, sw, agg};
    ST_int k;
    ctools_cancel_team team = {0};
#pragma omp parallel reduction(+ : nonconv) reduction(max : most) if (K > 1)
    {
#pragma omp for schedule(dynamic) nowait
        for (k = 0; k < K; k++) {
            if (ctools_cancel_requested()) continue;
            ST_double *ws = (ST_double *)malloc(((size_t)5 * n + 2 * (size_t)N) * sizeof(ST_double));
            if (!ws) {
                nonconv += 1000000;
                continue;
            }
            ST_double *col = data + (size_t)k * N;
            if (sw)
                for (ST_int r = 0; r < N; r++) col[r] *= sw[r];
            ST_int it = lsmr_column(&op, col, tol, maxiter, ws);
            if (sw)
                for (ST_int r = 0; r < N; r++) col[r] = sw[r] > 0.0 ? col[r] / sw[r] : 0.0;
            if (it < 0) nonconv++;
            else if (it > most) most = it;
            free(ws);
        }
        ctools_cancel_wait(&team);
    }
    if (ctools_cancel_requested()) rc = 1;
    else if (nonconv >= 1000000) rc = 920;
    else if (nonconv) rc = 430;
    if (iterations) *iterations = most;
done:
    free(off); free(d); free(agg); free(sw);
    return rc;
}

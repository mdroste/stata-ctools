/*
 * hdfe_wild.h
 *
 * Wild-bootstrap engine shared by the HDFE estimators: the wild cluster
 * bootstrap of a linear model after its fixed effects are partialled out,
 * restricted (WCR) or unrestricted (WCU), by the algorithm of Roodman,
 * MacKinnon, Nielsen and Webb (2019); enumeration of Rademacher sign
 * vectors; bootstrap p-values; and confidence sets by test inversion.
 * Design: docs/development/BOOTSTRAP_PLAN.md (sections 5.9 and 6.10).
 *
 * Notation (all on the partialled data): X is N x K, w the weights, u the
 * residuals of the full fit, A = (X'WX)^-1, b the estimates. For a tested
 * coefficient k, a = A e_k, rho = X a, and r = rho / a_k is x_k's residual
 * on the other regressors (FWL). Under H0: b_k = b0, with delta = b_k - b0,
 * the restricted residuals are u + delta r. Replicate weights v are drawn
 * per bootstrap cluster; the statistic uses the cluster-robust variance on
 * the vce() clusters, without small-sample factors (they cancel in every
 * p-value).
 *
 * For each replicate and tested coefficient the engine keeps five numbers
 * (n0, n1, D00, D01, D11) such that for any b0
 *     t*(b0) = (n0 + delta n1) / sqrt(D00 + 2 delta D01 + delta^2 D11),
 * so p(b0) is exact for every b0 at O(B). WCU keeps n1 = D01 = D11 = 0.
 *
 * Fixed effects (plan 5.9): on the aggregate path, the perturbed residuals
 * are used as they are. That is exact when every absorbed effect is nested
 * in the bootstrap clusters, and is the fedgp(partialled) convention
 * otherwise. With projected set, the caller's project() callback partials
 * the perturbed residuals by the fixed effects in every replicate: the
 * fedgp(full) convention, which the bootstrap of the full model gives.
 *
 * Ids are 0-based. Nothing here calls the Stata SPI.
 * Part of the ctools Stata plugin suite
 */

#ifndef CTOOLS_HDFE_WILD_H
#define CTOOLS_HDFE_WILD_H

#include <stddef.h>
#include <stdint.h>

#include "stplugin.h"

/* p-value types, as boottest's ptype() */
enum
{
    HDFE_WILD_SYMMETRIC = 0,
    HDFE_WILD_EQUALTAIL = 1,
    HDFE_WILD_LOWER = 2,
    HDFE_WILD_UPPER = 3
};

/* The fitted model. Arrays are the caller's and are read during setup
 * only. */
typedef struct
{
    ST_int N, K;
    const ST_double *X;     /* N x K partialled regressors, column-major */
    const ST_double *w;     /* N weights, or NULL for unweighted */
    const ST_double *resid; /* N partialled residuals of the full fit */
    const ST_double *A;     /* K x K (X'WX)^-1, column-major */
    const ST_double *beta;  /* K estimates */
    const ST_int *vce_id;   /* N vce() clusters 0 .. G-1, or NULL: each row */
    ST_int G;
    const ST_int *boot_id; /* N bootstrap clusters 0 .. C-1, or NULL: each row */
    ST_int C;
    ST_int n_test;       /* tested coefficients */
    const ST_int *test;  /* their columns, 0 .. K-1 */
    int null_imposed;    /* 1: WCR; 0: WCU */
    int projected;       /* 1: project the perturbed residuals (fedgp(full)) */
} hdfe_wild_input;

/*
 * Partial ncols columns of N rows (column-major, in place) by the fixed
 * effects, with the fit's weights. Called on the caller's thread.
 */
typedef ST_retcode (*hdfe_wild_project_fn)(void *ctx, ST_double *cols, ST_int ncols);

typedef struct
{
    ST_int N, K, G, C, H, n_test, B;
    int null_imposed, projected;
    ST_int done; /* replicates computed so far (in order) */

    /* Intersection groups h of the vce() and bootstrap clusters, listed by
     * vce() cluster: the groups of cluster g are g_start[g] .. g_start[g+1]-1 */
    ST_int *g_start;  /* G + 1 */
    ST_int *h_boot;   /* H: bootstrap cluster of each group */
    ST_double *S_u;   /* H x K: X_h'W_h u_h (row h at h*K) */
    ST_double *S_r;   /* n_test x H x K: X_h'W_h r_h (WCR only) */
    ST_double *q_u;   /* n_test x H: a'S_u */
    ST_double *q_r;   /* n_test x H: a'S_r (WCR only) */
    ST_double *d;     /* n_test x G x K: A X_g'W_g rho_g */
    ST_double *a;     /* n_test x K: A e_k */
    ST_double *b_hat; /* n_test: tested estimates */
    ST_double *se;    /* n_test: raw cluster-robust standard errors */

    /* Projected path only: row data */
    ST_int *row_g;    /* N: vce() cluster of each row */
    ST_int *row_c;    /* N: bootstrap cluster of each row */
    ST_double *w;     /* N weights (1 when unweighted) */
    ST_double *u;     /* N residuals */
    ST_double *rho;   /* n_test x N */
    ST_double *r;     /* n_test x N (WCR only) */

    ST_double *five;  /* n_test x B x 5: n0, n1, D00, D01, D11 */
    ST_double *obs;   /* n_test x 5: the same for v = 1, the observed statistic */
    int obs_ready;    /* on the projected path, set by the first run */
} hdfe_wild_model;

/*
 * Build the model for B replicates. Returns 0, 198 for invalid input
 * (an id out of range, a tested column out of range, a nonpositive
 * a_k), or 920.
 */
ST_retcode hdfe_wild_setup(const hdfe_wild_input *in, ST_int B, hdfe_wild_model *m);
void hdfe_wild_free(hdfe_wild_model *m); /* idempotent */

/*
 * Compute replicates first .. first + count - 1 (first must equal m->done).
 * v: count x C weights, replicate-major (row b holds replicate first + b),
 * or NULL to use the Rademacher sign vectors of enumeration, replicate j
 * being hdfe_wild_enum_weights(C, j). Up to `threads` replicates run at once
 * on the aggregate path. The projected path needs project() and partials
 * batches of replicates, about HDFE_WILD_BATCH_BYTES of columns at a time
 * (at least one replicate, whatever its size, within budget bytes).
 * Results do not depend on threads, chunks or batches.
 * Returns 0, 198, 909 (one replicate's columns exceed budget), 920, or
 * project()'s error.
 */
#define HDFE_WILD_BATCH_BYTES ((size_t)1 << 31)
ST_retcode hdfe_wild_run(hdfe_wild_model *m, const ST_double *v, ST_int first, ST_int count,
                         ST_int threads, hdfe_wild_project_fn project, void *ctx,
                         size_t budget);

/* Observed statistic of tested coefficient j at b0 (raw scale, as t*):
 * computed as the replicate with every weight 1, so that the replicates
 * with v = 1 and v = -1 tie with it exactly, as in boottest. NaN before the
 * first run on the projected path. */
double hdfe_wild_t(const hdfe_wild_model *m, ST_int j, double b0);

/* Replicate b's statistic of tested coefficient j at b0: NaN when its
 * standard error is zero. */
double hdfe_wild_tstar(const hdfe_wild_model *m, ST_int j, ST_int b, double b0);

/*
 * Bootstrap p-value of H0: b_j = b0 over the replicates with a defined
 * statistic (boottest's getp: strict inequalities, denominator the
 * feasible replicates). *feasible receives their count; p is NaN when it
 * is 0. Needs all B replicates.
 */
double hdfe_wild_pvalue(const hdfe_wild_model *m, ST_int j, double b0, int ptype,
                        ST_int *feasible);

/* Enumeration (boottest): all 2^C Rademacher vectors replace B draws when
 * C ln 2 < ln B + 1e-6. Returns 2^C then, otherwise 0. */
ST_int hdfe_wild_enumerable(ST_int C, ST_int B);

/* Sign vector j (0 .. 2^C - 1) in boottest's order (count_binary): v_c is
 * +1 when bit C-1-c of j is set, else -1. */
void hdfe_wild_enum_weights(ST_int C, uint64_t j, ST_double *v);

/*
 * Confidence set {b0 : p(b0) > alpha} by test inversion: p is evaluated on
 * a grid around center +/- 2 z scale, widened outward until p <= alpha,
 * and each bound refined by bisection to 1e-6 scale. A side with no bound
 * within 1000 scale of center is open (its bound NaN). disjoint: the grid
 * found a rejected point between accepted ones; lo and hi then bound the
 * whole set. empty: no grid point was accepted (both bounds NaN). A trial
 * value whose p is NaN (its evaluation failed) is treated as not rejected,
 * so that a failure never narrows the set, and counted in failed.
 */
typedef double (*hdfe_wild_pfn)(void *ctx, double b0);

typedef struct
{
    double lo, hi;
    int open_lo, open_hi, disjoint, empty;
    int failed; /* trial values whose p could not be evaluated */
} hdfe_wild_ci;

void hdfe_wild_invert(hdfe_wild_pfn p, void *ctx, double center, double scale, double alpha,
                      hdfe_wild_ci *ci);

/*
 * Score bootstrap (Kline and Santos 2012) of an ML estimate, as boottest
 * computes it for one null restriction: e is each row's contribution to the
 * efficient score of the tested coefficient at the restricted estimate. With
 * E_h its sums over the groups h of vce() and bootstrap clusters, a
 * replicate's statistic is t* = N* / sqrt(D*) with N* = sum_h v_h E_h and
 * D* = sum_g (sum_{h in g} v_h E_h - share_g N*)^2 over the vce() clusters g,
 * share_g the cluster's share of the rows (of share_w when given); the
 * observed statistic is the one for v = 1, so that v = 1 and v = -1 tie with
 * it exactly. Weights v as for hdfe_wild_run: B x C replicate-major, or NULL
 * for enumeration (B = 2^C).
 */
typedef struct
{
    ST_int N, G, C, H, B;
    ST_int *g_start; /* G + 1: groups of each vce() cluster */
    ST_int *h_boot;  /* H: bootstrap cluster of each group */
    ST_int *row_h;   /* N: group of each row */
    ST_double *share; /* G */
    ST_double *v;    /* B x C weights, or NULL for enumeration */
    ST_int got;      /* weight rows received so far (B when enumerating) */
} hdfe_wild_score_model;

/* Returns 0, 198 (ids out of range; enumerate without B = 2^C) or 920. */
ST_retcode hdfe_wild_score_setup(ST_int N, const ST_int *vce_id, ST_int G, const ST_int *boot_id,
                                 ST_int C, const ST_double *share_w, ST_int B, int enumerate,
                                 hdfe_wild_score_model *m);
/* Rows first .. first + count - 1 of the weights (count x C, replicate-major),
 * in order: first must equal the rows received so far. Returns 0 or 198. */
ST_retcode hdfe_wild_score_weights(hdfe_wild_score_model *m, const ST_double *v, ST_int first,
                                   ST_int count);
/* t[0] the observed statistic, t[1..B] the replicates' (NaN when D* = 0).
 * Returns 0 or 920. */
ST_retcode hdfe_wild_score_t(const hdfe_wild_score_model *m, const ST_double *e, ST_int threads,
                             ST_double *t);
void hdfe_wild_score_free(hdfe_wild_score_model *m); /* idempotent */

/* boottest's p-value from t[0] (observed) and t[1..B], as hdfe_wild_pvalue. */
double hdfe_wild_pvalue_t(const ST_double *t, ST_int B, int ptype, ST_int *feasible);

/* hdfe_wild_invert with a grid of 2 half + 1 points and bisection to
 * rel_tol scale (hdfe_wild_invert: 32 and 1e-6), for costly p-functions. */
void hdfe_wild_invert_ex(hdfe_wild_pfn p, void *ctx, double center, double scale, double alpha,
                         int half, double rel_tol, hdfe_wild_ci *ci);

/*
 * p at n points b0[0 .. n) at once, into pv. Query k belongs to a chain of
 * queries whose state (a warm start) it may begin from: from[k] the chain
 * whose last point it starts from (-1: none), to[k] the chain it extends.
 * No query of one call extends a chain that another query of that call
 * starts from, so the n points may be evaluated concurrently.
 */
typedef void (*hdfe_wild_pbatch_fn)(void *ctx, const double *b0, const int *from, const int *to,
                                    int n, double *pv);

/*
 * hdfe_wild_invert_ex evaluating the same points, so the same set when p
 * depends on b0 alone, but every grid point in one call and the searches
 * for the two bounds side by side (a point of each per call). Chains: grid
 * point i starts chain i (0 .. 2 half); the upper bound's search starts
 * from the chain of its accepted grid point and extends chain 2 half + 1,
 * the lower bound's likewise into chain 2 half + 2.
 */
void hdfe_wild_invert_batch(hdfe_wild_pbatch_fn p, void *ctx, double center, double scale,
                            double alpha, int half, double rel_tol, hdfe_wild_ci *ci);

#endif /* CTOOLS_HDFE_WILD_H */

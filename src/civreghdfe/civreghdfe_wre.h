/*
 * civreghdfe_wre.h
 *
 * The wild restricted efficient (WRE) bootstrap of Davidson and MacKinnon
 * (2010) for civreghdfe's k-class estimators (2SLS, LIML, Fuller, k-class),
 * as boottest computes it (docs/development/BOOTSTRAP_PLAN.md, 5.9 and 9.4).
 *
 * The model is the fit's data with the fixed effects partialled out: y,
 * the endogenous regressors Y2, the exogenous regressors X1 (partial()
 * columns included, as regressors) and the excluded instruments X2; the
 * instruments are X = [X1, X2] and the regressors Z = [Y2, X1]. For H0 on
 * coefficient k at b0 (null imposed), each point builds:
 *   - the restricted structural fit by LIML (as boottest, whatever the
 *     user's estimator), residual u;
 *   - the efficient first stage: Y2 on X and u, giving Pi and the residuals
 *     U2 = Y2 - X Pi;
 *   - the bootstrap data Y2* = X Pi + U2 .* v, y* = X1 g + Y2* b + u .* v
 *     with the restricted coefficients (g, b), the same weight v for both;
 *   - the user's estimator on (y*, Y2*) and its analytic variance (the
 *     unadjusted, robust or one-way cluster sandwich, with P_X Z in the
 *     middle as civreghdfe and ivreg2 compute it; small-sample factors
 *     cancel in every p-value).
 * nonull: the DGP is the unrestricted LIML fit, and the statistic is
 * centered on the estimate (boottest's WUE).
 *
 * Every column involved is a combination of the basis [y, Y2, X], so a
 * replicate needs only cross products of X with the perturbed basis
 * Q = M_D(basis .* v): X'WQ, Q'WQ and, for the sandwich, each vce()
 * cluster's X_g'W_g Q_g. When every absorbed effect is nested in the
 * bootstrap clusters, M_D(basis .* v) = basis .* v and these come from
 * cluster sums kept at setup; otherwise (fedgp(full)) the caller's project()
 * partials the perturbed basis by the fixed effects. Any b0 is then a few
 * small dense solves per replicate, which makes test inversion affordable.
 *
 * The observed statistic is replicate 0, every weight 1, computed by the
 * same code as the replicates, so that a replicate with every weight 1 ties
 * with it exactly, as in boottest.
 *
 * Ids are 0-based. Nothing here calls the Stata SPI.
 * Part of the ctools Stata plugin suite
 */

#ifndef CIVREGHDFE_WRE_H
#define CIVREGHDFE_WRE_H

#include <stddef.h>
#include <stdint.h>

#include "../stplugin.h"

typedef struct
{
    ST_int N, K2, K1, Kz;        /* rows; endogenous, exogenous regressors; excluded instruments */
    const ST_double *y;          /* N */
    const ST_double *Y2;         /* N x K2, column-major */
    const ST_double *X1;         /* N x K1 */
    const ST_double *X2;         /* N x Kz */
    const ST_double *w;          /* N weights, or NULL */
    int est;                     /* 0 2SLS, 1 LIML, 2 Fuller, 3 k-class */
    double kclass;               /* est 3 */
    double fuller_adj;           /* est 2: alpha / (N_eff - K_iv), the fit's */
    int vce;                     /* 0 unadjusted, 1 robust, 2 one-way cluster */
    const ST_int *vce_id;        /* vce 2: N clusters 0 .. G-1 */
    ST_int G;
    const ST_int *boot_id;       /* N bootstrap clusters 0 .. C-1, or NULL: each row */
    ST_int C;
    ST_int n_test;
    const ST_int *test;          /* regressor positions in Z = [Y2, X1] */
    int null_imposed;            /* 1 WRE, 0 WUE */
    int projected;               /* 1: project the perturbed basis (project() needed) */
} civ_wre_input;

/* Partial ncols columns of N rows (column-major, in place) by the fixed
 * effects with the fit's weights; called on the caller's thread. */
typedef ST_retcode (*civ_wre_project_fn)(void *ctx, ST_double *cols, ST_int ncols);

typedef struct civ_wre_model_s civ_wre_model;

/* Build the model for B replicates. Returns 0, 198 (bad dimensions or ids,
 * an underidentified model, a singular X'WX), or 920. */
ST_retcode civ_wre_setup(const civ_wre_input *in, ST_int B, civ_wre_model **out);
void civ_wre_free(civ_wre_model *m);

ST_int civ_wre_B(const civ_wre_model *m);
ST_int civ_wre_C(const civ_wre_model *m);

/*
 * Replicates first .. first + count - 1 (0-based, in order) of the weights:
 * count x C values, replicate-major, or NULL for the Rademacher sign
 * vectors of enumeration (hdfe_wild_enum_weights). Discrete weights are
 * kept as one byte each. Returns 0, 198, 909 (they would exceed budget
 * bytes) or 920.
 */
ST_retcode civ_wre_weights(civ_wre_model *m, const ST_double *v, ST_int first, ST_int count,
                           size_t budget);

/*
 * After the last weights: each replicate's cross products, and on the
 * projected path its projected basis (kept when it fits in half of
 * CIV_WRE_SOFT_BYTES, or of budget when smaller, otherwise projected again at
 * each evaluation, in batches of about a quarter of that). One replicate's
 * projection may take more, up to budget bytes. Returns 0, 198, 909 (one
 * replicate's projection exceeds budget), 920 or project()'s error.
 */
#define CIV_WRE_SOFT_BYTES ((size_t)1 << 32)
ST_retcode civ_wre_prepare(civ_wre_model *m, ST_int threads, civ_wre_project_fn project,
                           void *ctx, size_t budget);

/* 1 when prepare() kept the projected basis (or none was needed) */
int civ_wre_cached(const civ_wre_model *m);

/*
 * Tested coefficient j at npts points b0[0 .. npts): for replicate r = 0
 * (observed) .. B, the replicate's estimate bk[p (B+1) + r] of the tested
 * coefficient and its raw variance var[p (B+1) + r] (NaN when the fit or
 * the DGP fails). With the null imposed, t = (bk - b0) / sqrt(var); without
 * it the points do not matter (the DGP is the unrestricted fit). Results do
 * not depend on threads. Returns 0, 198, 909, 920 or project()'s error.
 */
ST_retcode civ_wre_eval(civ_wre_model *m, ST_int j, ST_int npts, const double *b0, double *bk,
                        double *var, ST_int threads, civ_wre_project_fn project, void *ctx,
                        size_t budget);

#endif /* CIVREGHDFE_WRE_H */

/*
 * creghdfe_core.h
 *
 * The numerical core of creghdfe's fit on partialled data, without the Stata
 * SPI: the full-sample fit and bootstrap replicates share it.
 * Part of the ctools Stata plugin suite
 */

#ifndef CREGHDFE_CORE_H
#define CREGHDFE_CORE_H

#include "creghdfe_types.h"

/*
 * Collinearity of the partialled regressors, as reghdfe judges it.
 *
 * data: N x K column-major, y then the K - 1 regressors, partialled, in the
 * units that were partialled; norm2: each column's (unweighted) sum of
 * squares before partialling, both counting each row count_w times (NULL:
 * once; bootstrap replicates pass their multiplicities so that a row drawn
 * twice counts as the expanded data would). A column whose partialled sum of squares is
 * at most collinear_tol * norm2 is collinear with the fixed effects; for y
 * that means the effects explain it, and its column is zeroed.
 * weights: the N weights or NULL; weight_type 1 aw, 2 fw, 3 pw.
 *
 * Outputs (caller-allocated, K - 1 entries unless noted): fe_flagged (the
 * fixed-effect rule), is_collinear (that rule or reghdfe_rmcoll's invsym
 * rule), *depvar_explained, and xtx ((K-1)^2) and xty_all: the (weighted)
 * cross products of the partialled regressors with themselves and with y,
 * which the caller may reuse. scratch: (K-1)^2 doubles.
 * Returns 0 or 920.
 */
ST_retcode creghdfe_collinear(ST_double *data, ST_int N, ST_int K, const ST_double *norm2,
                              const ST_double *count_w, const ST_double *weights, ST_int weight_type,
                              ST_double collinear_tol, ST_int *fe_flagged,
                              ST_int *is_collinear, ST_int *depvar_explained, ST_double *xtx,
                              ST_double *xty_all, ST_double *scratch);

typedef struct
{
    ST_double *data;            /* N x K partialled columns, y first */
    ST_int N, K;
    const ST_int *is_collinear; /* K - 1 flags from creghdfe_collinear */
    const ST_double *xtx;       /* cross products from creghdfe_collinear */
    const ST_double *xty_all;
    const ST_double *weights;   /* N or NULL */
    ST_int weight_type;
    ST_double sum_weights;      /* sum of the fweights (the 1'W1 term) */
    const ST_double *means;     /* K weighted means of the unpartialled columns */
    const ST_double *stdevs;    /* K scales the columns were divided by */
    ST_int standardize;
    ST_int use_quad;
} creghdfe_ols_in;

/* Every array is allocated by creghdfe_ols and owned by the caller, also
 * after a failure (free each one; data_keep may alias in->data). */
typedef struct
{
    ST_int K_keep;
    ST_int *keep_idx;       /* K_keep columns of data kept (1-based: y is 0) */
    ST_double *data_keep;   /* N x (K_keep + 1): y and the kept regressors */
    ST_double *xtx_keep;    /* K_keep^2 */
    ST_double *xty_keep;    /* K_keep */
    ST_double *beta_keep;   /* K_keep + 1, _cons last, in the partialled units */
    ST_double *inv_xx_keep; /* (K_keep + 1)^2, with the constant */
    ST_double *means_x;     /* K_keep means of the kept regressors, in those units */
    ST_double tss_within;
    ST_double rss;
} creghdfe_ols_out;

/* OLS of y on the kept regressors and the constant. Returns 0, 920, 198
 * when X'X is not positive definite, or the inversion's error. */
ST_retcode creghdfe_ols(const creghdfe_ols_in *in, creghdfe_ols_out *out);

#endif /* CREGHDFE_CORE_H */

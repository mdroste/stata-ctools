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
 * The regression's columns, read from the loaded variables.
 *
 * The matrix that is partialled holds the kept rows of every column divided
 * by its standard deviation. Each column needs its mean, then its sum of
 * squared deviations (which needs the mean), then its scale (which needs
 * that sum): three passes over the column at least. They are made here on
 * the loaded variable, and the matrix is written once, already scaled, by
 * the pass that also adds up its sum of squares. At 10M rows the columns do
 * not fit a cache, so the time of this step is the memory it moves: copying
 * the rows first and then standardizing the copy in place moves nearly twice
 * as much.
 *
 * Every sum runs over the kept rows in row order, one column per thread, so
 * the results do not depend on the thread count.
 */
typedef struct
{
    const ST_double *const *cols; /* K loaded columns, y first */
    const ST_double *level;       /* NULL, or K entries: level[k] >= 0 makes column k the
                                     indicator of its loaded value being level[k] */
    const ST_int *rows;           /* the loaded row of each kept row, ascending; NULL when
                                     every loaded row is kept */
    ST_int N;                     /* kept rows */
    ST_int K;
    const ST_double *weights;     /* N weights of the kept rows, or NULL */
} creghdfe_columns;

/* sums[k] (K): the sum over the kept rows of weight * column k, and with
 * weights *sum_w: their sum. */
void creghdfe_column_sums(const creghdfe_columns *C, ST_double *sums, ST_double *sum_w);

typedef struct
{
    const ST_double *means; /* K (weighted) means of the columns */
    ST_int weight_type;     /* 2: frequency weights, whose sum counts the observations */
    ST_double sum_weights;
    ST_int standardize;     /* divide each column by its standard deviation */
    ST_int use_quad;        /* sums of squared deviations in double-double arithmetic */
    ST_int uncentered;      /* tss about zero (a regression without a constant) */
} creghdfe_columns_opts;

/*
 * data (N x K, column-major): the kept rows of every column, divided by
 * stdevs[k] when standardize is set. raw (N x K, or NULL): the same rows
 * undivided. tss[k]: the (weighted) sum of squared deviations from means[k],
 * or the sum of squares with uncentered. stdevs[k]: its root over N - 1 (the
 * sum of the frequency weights less one), 1 when that is below 1e-30; with
 * standardize the dependent variable's (column 0) is replaced by the power
 * of two at or below it, which changes no bit of its mantissas. norm2[k]:
 * the unweighted sum of squares of column k of data. C->weights are the
 * weights the regression uses (normalized where it normalizes them).
 */
void creghdfe_columns_build(const creghdfe_columns *C, const creghdfe_columns_opts *opts,
                            ST_double *data, ST_double *raw, ST_double *tss, ST_double *stdevs,
                            ST_double *norm2);

/*
 * Column k alone: its entry of creghdfe_column_sums(), and its part of
 * creghdfe_columns_build() (column k of data and raw, entry k of the other
 * arrays; opts->means[k] must hold its mean). For a caller that runs the
 * columns as tasks of its own, any number of them at once: a column's passes
 * read the loaded variables and write nothing but what is its own, on the
 * calling thread.
 */
ST_double creghdfe_column_sum(const creghdfe_columns *C, ST_int k);
void creghdfe_column_build(const creghdfe_columns *C, const creghdfe_columns_opts *opts, ST_int k,
                           ST_double *data, ST_double *raw, ST_double *tss, ST_double *stdevs,
                           ST_double *norm2);

/*
 * A regression with a single absorbed effect (no slope): the projection on
 * it done while the column is written.
 *
 * The projection subtracts from every value the (weighted) mean of its
 * level. Done on the finished matrix it reads every column twice more (the
 * level sums, then the subtraction) after the pass that wrote it. Here the
 * values pass through once less: the level sums are taken from the scaled
 * values as they are computed, and the column is written already projected.
 * At 10M rows in no order that takes creghdfe's partial-out phase from 52
 * to 39 ms with three columns and from 124 to 86 ms with nine. (Rows grouped
 * by the effect are better left to partial_out_columns(), which takes them
 * run by run on every thread: the phase was 44 ms that way and 46 ms with
 * the runs handled here, a column per thread.)
 * Every number is the one creghdfe_column_build() followed by
 * partial_out_columns() gives, to the last bit: norm2[k] is the sum of
 * squares of the column before the projection, added in the same eight
 * chains, and each level's sum adds the same scaled values in row order.
 *
 * Two calls per column, since the effect's counts may still be on their way
 * while the first runs (it does not read them):
 *
 *   creghdfe_column_scale()   opts->means[k] must hold the column's mean.
 *       Sets tss[k] and stdevs[k], as creghdfe_column_build() does.
 *   creghdfe_column_absorb()  sets norm2[k] and writes column k of data.
 *       sums (E->num_levels entries) is scratch for the level sums.
 *
 * Not for the double-double sums (use_quad) nor for a caller that wants the
 * unscaled rows as well: those take creghdfe_column_build().
 */
typedef struct
{
    const ST_int *levels;    /* the 1-based level of every kept row */
    ST_int num_levels;
    const ST_double *inv;    /* 1 / the (weighted) count of each level; NULL: divide by counts */
    const ST_double *counts; /* the (weighted) counts, read when inv is NULL */
} creghdfe_absorbed;

void creghdfe_column_scale(const creghdfe_columns *C, const creghdfe_columns_opts *opts, ST_int k,
                           ST_double *tss, ST_double *stdevs);
void creghdfe_column_absorb(const creghdfe_columns *C, const creghdfe_columns_opts *opts, ST_int k,
                            const creghdfe_absorbed *E, const ST_double *stdevs, ST_double *norm2,
                            ST_double *sums, ST_double *data);

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

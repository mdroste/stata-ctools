/*
 * cbinscatter_fit.h
 *
 * Line fitting functions for cbinscatter
 * Part of the ctools Stata plugin suite
 */

#ifndef CBINSCATTER_FIT_H
#define CBINSCATTER_FIT_H

#include "cbinscatter_types.h"

/* ========================================================================
 * Line Fitting
 * ======================================================================== */

/*
 * Fit lines for all by-groups in results
 *
 * Parameters:
 *   y           - y values (full dataset)
 *   x           - x values (full dataset)
 *   group_start - first row of each by-group and the end of the last
 *                 (rows ordered by group); NULL if no by
 *   weights     - weights (NULL if unweighted)
 *   N           - total observations
 *   config      - binscatter configuration
 *   results     - results structure (groups will be updated with fit coefs)
 *
 * Returns:
 *   0 on success, error code on failure
 */
ST_retcode fit_all_groups(
    const ST_double *y,
    const ST_double *x,
    const ST_int *group_start,   /* num_by_groups + 1 row offsets, NULL without by() */
    const ST_double *weights,
    ST_int N,
    const BinscatterConfig *config,
    BinscatterResults *results
);

/*
 * Fit one group's polynomial of order config->linetype, jointly with K
 * controls and G absorbed effects when given (method(binsreg), binsreg's
 * polyreg() convention: the line is evaluated at the controls' weighted
 * means). Sets group->fit_coefs (powers of x), fit_order, fit_r2, fit_n.
 *
 * Parameters:
 *   controls    - control variables (column-major: K x N), NULL if K == 0
 *   fe_vars     - FE levels 1..L (column-major: G x N), NULL if G == 0
 */
ST_retcode fit_group_polynomial(
    const ST_double *y,
    const ST_double *x,
    const ST_double *weights,
    ST_int N,
    const ST_double *controls,
    ST_int K,
    const ST_int *fe_vars,
    ST_int G,
    const BinscatterConfig *config,
    ByGroupResult *group
);

/*
 * rd(): fit each group's polynomial separately on each side of the break
 * points config->rd[0..n_rd-1], as binscatter's rd() does: interval 0 is
 * x <= rd[0], interval r is rd[r-1] < x <= rd[r], the last x > rd[n_rd-1].
 * results->rd_coefs and rd_stats (allocated here) get one row per group and
 * interval; an interval without observations keeps missing values. An
 * interval with one observation is an error (2001), as regress makes it.
 */
ST_retcode fit_all_groups_rd(
    const ST_double *y,
    const ST_double *x,
    const ST_int *group_start,   /* num_by_groups + 1 row offsets, NULL without by() */
    const ST_double *weights,
    ST_int N,
    const BinscatterConfig *config,
    BinscatterResults *results
);

#endif /* CBINSCATTER_FIT_H */

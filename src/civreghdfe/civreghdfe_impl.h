/*
    civreghdfe_impl.h
    Instrumental Variables Regression with High-Dimensional Fixed Effects

    This module implements 2SLS/IV regression with HDFE absorption,
    reusing the creghdfe infrastructure for HDFE demeaning.
*/

#ifndef CIVREGHDFE_IMPL_H
#define CIVREGHDFE_IMPL_H

#include "../stplugin.h"
#include "../ctools_types.h"
#include "hdfe/hdfe_types.h"
#include "hdfe/hdfe_solver.h"
#include "../ctools_ols.h"
#include "hdfe/hdfe_remap.h"

/*
    Main entry point for civreghdfe plugin calls.

    Subcommands:
    - "iv_regression" - Full IV/2SLS regression with HDFE

    Returns STATA_OK on success, error code on failure.
*/
ST_retcode civreghdfe_main(const char *args);

/* The full fit (iv_regression); civreghdfe_main dispatches to it */
ST_retcode civreghdfe_fit(void);

/* Shared with the bootstrap replicates (civreghdfe_boot.c) */

/* A column is absorbed by the fixed effects (or by partial()) when its
   partialled (weighted) SS is at most CIV_FE_NOISE times its uncentered SS
   before partialling (civreghdfe_impl.c) */
#define CIV_FE_NOISE 1e-20

/* Weighted sum of squares sum(w x^2), and about the weighted mean; w may be NULL */
ST_double civ_wss(const ST_double *x, const ST_double *w, ST_int N);
ST_double civ_wss_dev(const ST_double *x, const ST_double *w, ST_int N);

/* partial(): FWL residualization of every column on the partial() columns
   (see civreghdfe_impl.c); messages only with report */
ST_retcode civ_fwl_partial_out(ST_double *data, ST_int N, ST_int total_cols,
                               ST_int exog_start, ST_int zexog_start, ST_int K_exog,
                               ST_int *is_partial, const ST_double *w,
                               const ST_double *ss_exog, ST_double absorb_tol, int report);

/* Columns Stata's invsym() drops from the symmetric K x K A (see
   civreghdfe_impl.c); returns the number newly dropped, or -1 */
ST_int civ_invsym_drops(const ST_double *A, ST_int K, ST_int n_ordered, ST_int *drop);

/*
 * Cleanup function for civreghdfe persistent state.
 * Frees the global HDFE state if allocated.
 * Safe to call multiple times (idempotent).
 * Note: civreghdfe normally cleans up after itself, but this
 * handles cases where execution was interrupted.
 */
void civreghdfe_cleanup_state(void);

#endif /* CIVREGHDFE_IMPL_H */

/*
 * hdfe_solver.h
 *
 * CG solver, Kaczmarz transforms, FE projections
 * Part of the ctools Stata plugin suite
 */

#ifndef CTOOLS_HDFE_SOLVER_H
#define CTOOLS_HDFE_SOLVER_H

#include "hdfe_types.h"
#include "../ctools_cancel.h"
#include <limits.h>

/* ========================================================================
 * Conjugate Gradient Solver
 * ======================================================================== */

/* Positive iteration count on convergence, negative count on exhaustion,
 * or this distinct sentinel on numerical breakdown (not overridable by abort(0)). */
#define HDFE_CG_BREAKDOWN INT_MIN
ST_int cg_solve_column_threaded(HDFE_State *S, ST_double *y, ST_int thread_id);

/* 1 when the positive values of weights[0..N-1] span more than the iteration
 * on the rows can certify at its rounding floor: there a residual that small
 * is a breakdown, since rows that light can hide an unresolved effect, unless
 * S->floor_any_weight_span. A caller whose own weights are only part of the
 * solver's (cppmlhdfe) applies the same test to them. */
int hdfe_weight_span_too_wide(const ST_double *weights, ST_int N);

/* ========================================================================
 * Partial Out Multiple Columns (Shared Helper)
 *
 * Partials out fixed effects from K columns of data in parallel.
 * Data should be in column-major format with N rows.
 *
 * Parameters:
 *   S: HDFE state (must be fully initialized)
 *   data: Column-major matrix (N x K), modified in place
 *   N: Number of observations
 *   K: Number of columns to partial out
 *   num_threads: Number of threads to use
 *
 * Returns: status (0, 1 for cancellation, 198 for invalid configuration,
 * 920 for allocation failure, 498 for numerical breakdown, or 430 for nonconvergence),
 * convergence flag, and maximum iterations used across columns.
 *
 * Two or more categorical effects (no slopes) are solved for their
 * coefficients (hdfe_coef.h), unless S->maxiter is below
 * HDFE_COEF_MIN_MAXITER. Everything else, and every call of
 * partial_out_columns_options(), runs the iteration on the rows: conjugate
 * gradients on the symmetric Kaczmarz sweep, reghdfe's default.
 * ======================================================================== */
typedef struct
{
    ST_retcode status;
    ST_int converged;
    ST_int iterations;
} HDFE_SolveResult;

HDFE_SolveResult partial_out_columns(HDFE_State *S, ST_double *data, ST_int N, ST_int K,
                                     ST_int num_threads);

/* partial_out_columns() for a caller that projects on the same effects again
 * and again while only S->weights changes (the iterations of a weighted
 * fit). plan (hdfe_coef.h; may be NULL) keeps what the solve for the
 * coefficients finds from the levels alone, so that later calls skip it; the
 * residuals are those of partial_out_columns(), to the last bit.
 *
 * The solve for the coefficients reads the weights but not the factors'
 * level counts (weighted_counts and their inverses); every other path does.
 * A caller that brings the counts up to date only when they are needed
 * passes before_rows: it is called, with ctx, before any of those paths
 * touches the data, and not at all when the coefficients are solved for. */
struct hdfe_coef_plan;
HDFE_SolveResult partial_out_columns_plan(HDFE_State *S, ST_double *data, ST_int N, ST_int K,
                                          ST_int num_threads, struct hdfe_coef_plan *plan,
                                          void (*before_rows)(void *ctx), void *ctx);

/* Optional transforms: 0 symmetric Kaczmarz, 1 Kaczmarz, 2 Cimmino, 3 random
 * (symmetric) Kaczmarz, hdfe's experimental transform: a new factor order
 * for every sweep, from a fixed per-column generator.
 * Acceleration: 0 none, 1 CG, 2 steepest descent, 3 Aitken, 4 hybrid (hdfe's
 * experimental acceleration: three plain sweeps, then CG). Every pair runs,
 * Kaczmarz with CG included (CG may then fail to converge, as in hdfe).
 * Always iterates on the rows: (0, 1) is the iteration partial_out_columns()
 * falls back on, for callers whose user named a transform or acceleration.
 * Uses the same state/buffer and column-major ownership contract above. */
HDFE_SolveResult partial_out_columns_options(HDFE_State *S, ST_double *data, ST_int N, ST_int K,
                                            int transform, int acceleration);

#endif /* CTOOLS_HDFE_SOLVER_H */

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

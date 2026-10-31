/*
 * hdfe_solver.h
 *
 * CG solver, Kaczmarz transforms, FE projections
 * Part of the ctools Stata plugin suite
 */

#ifndef CTOOLS_HDFE_SOLVER_H
#define CTOOLS_HDFE_SOLVER_H

#include "hdfe_types.h"

/* ========================================================================
 * Conjugate Gradient Solver
 * ======================================================================== */

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
 * Returns: status (0, 198 for invalid configuration, or 430 for nonconvergence),
 * convergence flag, and maximum iterations used across columns.
 * ======================================================================== */
typedef struct
{
    ST_retcode status;
    ST_int converged;
    ST_int iterations;
} HDFE_SolveResult;

HDFE_SolveResult partial_out_columns(HDFE_State *S, ST_double *data, ST_int N, ST_int K,
                                     ST_int num_threads);

/* Optional transforms: 0 symmetric Kaczmarz, 1 Kaczmarz, 2 Cimmino.
 * Acceleration: 0 none, 1 CG, 2 steepest descent, 3 Aitken.
 * Uses the same state/buffer and column-major ownership contract above. */
HDFE_SolveResult partial_out_columns_options(HDFE_State *S, ST_double *data, ST_int N, ST_int K,
                                            int transform, int acceleration);

#endif /* CTOOLS_HDFE_SOLVER_H */

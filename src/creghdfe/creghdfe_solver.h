/*
 * creghdfe_solver.h
 *
 * CG solver, Kaczmarz transforms, FE projections
 * Part of the ctools Stata plugin suite
 */

#ifndef CREGHDFE_SOLVER_H
#define CREGHDFE_SOLVER_H

#include "creghdfe_types.h"

/* ========================================================================
 * Dot product computation
 * ======================================================================== */

ST_double dot_product(const ST_double * RESTRICT x,
                      const ST_double * RESTRICT y,
                      ST_int N);

ST_double weighted_dot_product(const ST_double * RESTRICT x,
                               const ST_double * RESTRICT y,
                               const ST_double * RESTRICT w,
                               ST_int N);

/* ========================================================================
 * Sorted Indices for Fast Projection
 * ======================================================================== */

/* Build sorted observation indices using counting sort on levels.
 * This enables cache-friendly projection by accessing means sequentially.
 * Returns 0 on success, non-zero on error. */
int build_sorted_indices(FE_Factor *f, ST_int N);

/* Free sorted indices data */
void free_sorted_indices(FE_Factor *f);

/* ========================================================================
 * Symmetric Kaczmarz Transformation
 * ======================================================================== */

void transform_sym_kaczmarz_threaded(const HDFE_State * RESTRICT S,
                                     const ST_double * RESTRICT y,
                                     ST_double * RESTRICT ans,
                                     ST_double ** RESTRICT fe_means,
                                     ST_int thread_id);

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
 * Returns: Maximum iterations used across all columns (negative if any failed)
 * ======================================================================== */
ST_int partial_out_columns(HDFE_State *S, ST_double *data, ST_int N, ST_int K, ST_int num_threads);

#endif /* CREGHDFE_SOLVER_H */

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

#endif /* CREGHDFE_SOLVER_H */

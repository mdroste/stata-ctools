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
 * CSR Format for Fast Projection
 * ======================================================================== */

/* Build CSR (Compressed Sparse Row) format from levels array.
 * This enables cache-friendly projection by grouping observations by level.
 * Returns 0 on success, non-zero on error. */
int build_csr_format(FE_Factor *f, ST_int N);

/* Free CSR format data */
void free_csr_format(FE_Factor *f);

/* ========================================================================
 * Fixed Effect Projection
 * ======================================================================== */

void project_one_fe_threaded(const ST_double * RESTRICT y,
                             const FE_Factor * RESTRICT f,
                             ST_int N,
                             const ST_double * RESTRICT weights,
                             ST_double * RESTRICT proj,
                             ST_double * RESTRICT means);

/* CSR-accelerated projection - faster for large FEs */
void project_one_fe_csr(const ST_double * RESTRICT y,
                        const FE_Factor * RESTRICT f,
                        ST_int N,
                        const ST_double * RESTRICT weights,
                        ST_double * RESTRICT proj,
                        ST_double * RESTRICT means);

/* ========================================================================
 * Symmetric Kaczmarz Transformation
 * ======================================================================== */

void transform_sym_kaczmarz_threaded(const HDFE_State * RESTRICT S,
                                     const ST_double * RESTRICT y,
                                     ST_double * RESTRICT ans,
                                     ST_double * RESTRICT proj,
                                     ST_double ** RESTRICT fe_means,
                                     ST_int thread_id);

/* ========================================================================
 * Conjugate Gradient Solver
 * ======================================================================== */

ST_int cg_solve_column_threaded(HDFE_State *S, ST_double *y, ST_int thread_id);

#endif /* CREGHDFE_SOLVER_H */

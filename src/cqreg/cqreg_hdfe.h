/*
 * cqreg_hdfe.h
 *
 * HDFE (High-Dimensional Fixed Effects) integration for quantile regression.
 * Provides wrappers around creghdfe's CG solver for partialling out FEs.
 * Part of the ctools suite.
 */

#ifndef CQREG_HDFE_H
#define CQREG_HDFE_H

#include "cqreg_types.h"

/* Forward declaration to avoid circular dependency */
struct HDFE_State;

/* ============================================================================
 * Main Interface
 * ============================================================================ */

/*
 * Initialize HDFE state for quantile regression.
 *
 * This sets up the factor structures and CG solver state needed for
 * partialling out fixed effects. If fe_data is non-NULL, uses pre-loaded
 * arrays directly (no SPI reads). Otherwise reads FE variables from Stata.
 *
 * Parameters:
 *   state     - cqreg state (will have hdfe_state set)
 *   fe_vars    - Array of Stata variable indices for FE variables (1-indexed),
 *                or NULL if fe_data is provided
 *   G          - Number of FE variables
 *   N_filtered - Number of observations after if-condition filtering
 *   in1, in2   - Observation range (from Stata's _N())
 *   maxiter    - Maximum CG iterations (default: 10000)
 *   tolerance  - CG convergence tolerance (default: 1e-8)
 *   fe_data    - Pre-loaded FE data arrays [G][N_filtered], or NULL to
 *                read from Stata via SF_vdata/SF_ifobs
 *
 * Returns:
 *   0 on success, error code otherwise
 */
ST_int cqreg_hdfe_init(cqreg_state *state,
                       const ST_int *fe_vars,
                       ST_int G,
                       ST_int N_filtered,
                       ST_int in1, ST_int in2,
                       ST_int maxiter,
                       ST_double tolerance,
                       ST_double **fe_data);

/*
 * Partial out fixed effects from y and X.
 *
 * Uses the CG solver from creghdfe to project y and each column of X
 * onto the orthogonal complement of the FE space.
 *
 * This modifies y and X in place.
 *
 * Parameters:
 *   state - cqreg state with initialized hdfe_state
 *   y     - Dependent variable (N), modified in place
 *   X     - Design matrix (N x K, column-major), modified in place
 *   N     - Number of observations
 *   K     - Number of regressors
 *
 * Returns:
 *   0 on success, error code otherwise
 */
ST_int cqreg_hdfe_partial_out(cqreg_state *state,
                              ST_double *y,
                              ST_double *X,
                              ST_int N, ST_int K);

/*
 * Partial out a single column using the CG solver.
 *
 * Parameters:
 *   state     - cqreg state with initialized hdfe_state
 *   col       - Column vector (N), modified in place
 *   N         - Number of observations
 *   thread_id - Thread ID for parallel execution
 *
 * Returns:
 *   Number of CG iterations (positive), or negative on error
 */
ST_int cqreg_hdfe_partial_out_column(cqreg_state *state,
                                     ST_double *col,
                                     ST_int N,
                                     ST_int thread_id);

/*
 * Get degrees of freedom absorbed by fixed effects.
 *
 * Parameters:
 *   state - cqreg state with initialized hdfe_state
 *
 * Returns:
 *   Degrees of freedom absorbed (sum of (num_levels - 1) for each FE,
 *   minus mobility groups adjustment if G >= 2)
 */
ST_int cqreg_hdfe_get_df_absorbed(const cqreg_state *state);

/*
 * Get number of singletons dropped.
 *
 * Parameters:
 *   state - cqreg state with initialized hdfe_state
 *
 * Returns:
 *   Number of observations dropped due to singleton FE levels
 */
ST_int cqreg_hdfe_get_num_singletons(const cqreg_state *state);

/*
 * Clean up HDFE state.
 * Frees all memory associated with the HDFE solver.
 *
 * Parameters:
 *   state - cqreg state with hdfe_state to clean up
 */
void cqreg_hdfe_cleanup(cqreg_state *state);

/* ============================================================================
 * Internal Functions
 * ============================================================================ */

/*
 * Create factor from Stata variable.
 * Reads values and creates level assignments, filtering by if-condition.
 *
 * Parameters:
 *   var_idx    - Stata variable index (1-indexed)
 *   in1        - First observation
 *   in2        - Last observation
 *   N_filtered - Number of observations after if-condition filtering
 *   levels     - Output: level assignments (N_filtered)
 *   counts     - Output: count per level
 *   num_levels - Output: number of unique levels
 *
 * Returns:
 *   0 on success, error code otherwise
 */
ST_int cqreg_hdfe_create_factor(ST_int var_idx,
                                ST_int in1, ST_int in2,
                                ST_int N_filtered,
                                ST_int *levels,
                                ST_double *counts,
                                ST_int *num_levels);

/*
 * Create factor from pre-loaded data array.
 * Same as cqreg_hdfe_create_factor but reads from a pre-loaded
 * ST_double array instead of Stata's SPI, avoiding redundant I/O.
 */
ST_int cqreg_hdfe_create_factor_from_data(const ST_double *data,
                                           ST_int N_filtered,
                                           ST_int *levels,
                                           ST_double *counts,
                                           ST_int *num_levels);

/*
 * Detect and mark singleton observations.
 * Iteratively removes observations that belong to singleton FE levels.
 *
 * Parameters:
 *   state - cqreg state with initialized factors
 *   mask  - Output: 1 = keep, 0 = drop (N)
 *
 * Returns:
 *   Number of singletons dropped
 */
ST_int cqreg_hdfe_detect_singletons(cqreg_state *state, ST_int *mask);

#endif /* CQREG_HDFE_H */

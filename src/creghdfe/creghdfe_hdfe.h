/*
 * creghdfe_hdfe.h
 *
 * Global HDFE state (g_state) and its cleanup
 * Part of the ctools Stata plugin suite
 */

#ifndef CREGHDFE_HDFE_H
#define CREGHDFE_HDFE_H

#include "creghdfe_types.h"

/* Private to creghdfe; IV and PPML own independent state. */
extern HDFE_State *g_state;

/* Clean up global state */
void cleanup_state(void);

/* vce(bootstrap) hooks, installed by creghdfe_boot_install() when the plugin
 * dispatches creghdfe, so that the fit and its cleanup do not link the
 * bootstrap module (native harnesses compile them alone). NULL: absent. */
struct creghdfe_boot_source_s;
extern ST_retcode (*creghdfe_boot_retain_hook)(const struct creghdfe_boot_source_s *);
extern void (*creghdfe_boot_release_hook)(void);

/* ========================================================================
 * Heterogeneous slopes (absorb(id##c.x) and the like)
 *
 * Each intercept and each slope of an absorbed effect is a factor of the
 * HDFE state (reghdfe's extended absorbed effects). A slope factor projects
 * on its variable within each level, centered within the level when its
 * effect has an intercept, so that the intercept and slope projections of
 * a level are orthogonal. The solver alternates over these factors; at
 * convergence the partialled columns are the projection reghdfe computes
 * with all of an effect's coefficients at once.
 * ======================================================================== */
#include "creghdfe_input.h"

/* Slope factors of S (built by hdfe_build_factors on the kept rows), from
 * the slope columns of the N_orig loaded rows (slope_cols[g], NULL for an
 * intercept) and the loaded weights (NULL unweighted). Call before
 * ctools_hdfe_alloc_buffers: counts and weighted_counts of a slope factor
 * become the (weighted) sums of squares of its centered slope per level,
 * zero where reghdfe's invsym() drops the slope (constant within a level of
 * an effect with an intercept, zero without one). */
ST_retcode creghdfe_setup_slopes(HDFE_State *S, const creghdfe_input *in,
                                 const ST_double *const *slope_cols, const ST_int *mask,
                                 ST_int N_orig, const ST_double *weights);

/* Redundant coefficients of every factor, as reghdfe's estimate_dof()
 * counts them: intercepts as ctools_hdfe_redundant_intercepts() counts them
 * among the intercept factors (nested: NULL or G flags), and, with
 * dofadjustments(continuous), the levels whose standardized slope is
 * constant (with an intercept) or zero (without one). exact (G flags, may
 * be NULL) marks the counts reghdfe reports as exact: nested effects, the
 * first remaining intercept and, with mobility-group adjustments, the
 * second. Returns 0 or 920. */
ST_retcode creghdfe_redundant(const HDFE_State *S, const creghdfe_input *in, const ST_int *nested,
                              ST_int dof_pairs, ST_int *redundant, ST_int *exact);

/* hdfe_fe_recover() with slopes: the same sweeps and steepest-descent steps,
 * with each absorbed effect projected on all of its intercept and slopes at
 * once within each level, as reghdfe's store_alphas does. alpha[g] receives
 * factor g's coefficient per level: the intercept, or the slope (in the
 * units of the slope variable). Returns 0, 920 or 430. */
ST_retcode creghdfe_fe_recover_slopes(const HDFE_State *S, const creghdfe_input *in,
                                      ST_double *error, ST_double **alpha);

#endif /* CREGHDFE_HDFE_H */

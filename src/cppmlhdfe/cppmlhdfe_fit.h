/*
 * cppmlhdfe_fit.h
 *
 * The PPML fit's numerical phases without the Stata SPI, shared by the full
 * fit (do_ppml_regression) and the bootstrap replicates
 * (docs/development/BOOTSTRAP_PLAN.md, 8.1): the IRLS loop, the final
 * projection and collinearity check, and the final solve. Arrays are the
 * caller's; the functions hand back the ones they allocate in the io
 * struct, also on failure, for the caller to free.
 * Part of the ctools Stata plugin suite
 */

#ifndef CPPMLHDFE_FIT_H
#define CPPMLHDFE_FIT_H

#include "cppmlhdfe_internal.h"

typedef struct
{
    /* Inputs */
    int report;               /* SF_error / SF_display messages (full fit only) */
    int iteration_log;        /* ppmlhdfe's iteration log */
    const PPML_Options *options;
    HDFE_State *S;            /* factors on the sample, IRLS weights */
    ST_double *y, *X;         /* N standardized outcomes, N x K_x regressors */
    const ST_double *offset;  /* N, or NULL */
    const ST_double *w_user;  /* N user weights, or NULL */
    ST_double user_weight_scale, tolerance, irls_tol, stdev_y, inner_tol_override;
    ST_int N, K_x, K_aug, G, irls_maxiter, num_threads, verbose;
    ST_int has_weights, weight_type;
    ST_int K_irls;            /* regressors in IRLS after the collinearity pre-pass */
    ST_int *irls_keep_idx;    /* their columns, or NULL for all */
    /* Workspace (caller's) */
    ST_double *mu, *eta, *irls_w, *z_orig, *aug_copy, *previous_working;
    ST_double *xtx, *xty, *irls_beta, *inv_xx, *xb_work, *fe_work;
    ST_double *xtx_k, *xty_k, *beta_k, *inv_xx_k;
    ST_int *sep_m;
    /* ppml_irls_loop */
    ST_double *beta_prev, *eta_refit, *eta_refit_prev;
    ST_int iterations, converged;
    ST_double deviance, subiterations, target_inner_tol;
    /* ppml_final_collinear */
    ST_int *is_collinear;
    ST_int num_collinear, K_keep;
    /* ppml_final_solve */
    ST_int *keep_idx;
    ST_double *beta_final, *w_reg, *xtx_keep, *xty_keep, *inv_xx_final, *V_keep;
    dd_real *xty_dd;
    ST_double N_ref, rss;
} ppml_irls_io;

ST_retcode ppml_irls_loop(ppml_irls_io *io);
ST_retcode ppml_final_collinear(ppml_irls_io *io);
ST_retcode ppml_final_solve(ppml_irls_io *io);

/* Minimum-norm fixed-effect coefficients of d (N, overwritten) on the
 * factors of S with its weights, by CGLS from zero: alpha[g] holds
 * S->factors[g].num_levels zeros on entry. Returns 0, 430 (not converged
 * within S->maxiter) or 920. */
int ppml_fe_minnorm(HDFE_State *S, double *d, double **alpha, int report);

#endif /* CPPMLHDFE_FIT_H */

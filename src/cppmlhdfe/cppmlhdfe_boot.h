/*
 * cppmlhdfe_boot.h
 *
 * vce(bootstrap) for cppmlhdfe: the model retained after the full fit and
 * the plugin calls that run pairs-bootstrap replicates on it
 * (docs/development/BOOTSTRAP_PLAN.md, section 8). A replicate is the PPML
 * fit on the drawn rows, weighted by their multiplicities: the sample
 * screening on the multiplicities, the separation checks the full fit ran,
 * then the same IRLS and final solve (cppmlhdfe_fit.h), warm-started from
 * the full fit's linear predictor, in the full fit's standardized units.
 * Part of the ctools Stata plugin suite
 */

#ifndef CPPMLHDFE_BOOT_H
#define CPPMLHDFE_BOOT_H

#include "cppmlhdfe_options.h"
#include "cppmlhdfe_types.h"
#include "../ctools_types.h"

/* What the full fit hands over (borrowed; retain copies it). */
typedef struct
{
    ST_int N, K_x, G;             /* compact rows, regressors, effects */
    const ST_double *y, *X;       /* standardized outcome and regressors */
    const ST_double *offset;      /* N, or NULL */
    const ST_double *w_user;      /* N user weights, or NULL */
    ST_int weight_type;           /* 0 none, 1 aw, 2 fw, 3 pw */
    ST_int *const *levels;        /* G x N 1-based levels on the compact rows */
    const ST_int *num_levels;     /* G */
    const ST_double *eta;         /* N linear predictor without the offset */
    const ST_double *stdev_x;     /* K_x scales */
    ST_double stdev_y;
    const PPML_Options *options;  /* with the separation checks to run */
    ST_double tolerance, irls_tol, inner_tol_override;
    ST_int maxiter, irls_maxiter, keep_singletons, check_fe, check_relu;
    const ST_int *is_collinear;   /* K_x */
    ST_int K_keep;
    const ST_double *b_hat;       /* K_keep + 1 in the data's units, _cons last */
    const perm_idx_t *obs_map;    /* N 1-based Stata observations */
    /* Heterogeneous slopes: factor g's slope variable on the compact rows,
     * or NULL for an intercept; slope_centered[g] when it is centered
     * within its levels (an intercept on the same levels is absorbed) */
    const ST_double *const *slope;
    const ST_int *slope_centered;
} cppmlhdfe_boot_source;

/* Hooks the fit calls (NULL unless cppmlhdfe_boot_install ran), so that
 * harnesses that compile the fit alone still link. */
extern ST_retcode (*cppmlhdfe_boot_retain_hook)(const cppmlhdfe_boot_source *src);
extern void (*cppmlhdfe_boot_release_hook)(void);
/* full_regression retain: the fit calls the retain hook when it succeeds */
extern int cppmlhdfe_boot_retain_requested;

void cppmlhdfe_boot_install(void);

/* Keep the model; reads the resampling unit and strata variables (plugin
 * indices in scalars __cppmlhdfe_bs_unit_idx and _strata_idx) and sets
 * scalars __cppmlhdfe_bs_U and __cppmlhdfe_bs_S. Returns 0 or 920. */
ST_retcode cppmlhdfe_boot_retain(const cppmlhdfe_boot_source *src);
void cppmlhdfe_boot_release(void);

/* boot_chunk first=# count=# total=# (seedhi=# seedlo=# | stata | supplied)
 * [fe] [dse] [ci] [rows]: as creghdfe's; fe and dse (read at first=0) also
 * recover each replicate's fixed effects and their sum d. Statuses in
 * matrix __cppmlhdfe_bs_status (1 x count). */
ST_retcode cppmlhdfe_boot_chunk(const char *args);

/* boot_finish level=# [coef] [mse] [nreps]: with fe/dse, the fixed
 * effects' and d's SEs into the call's variables (see finish_effects); with
 * coef, the summary into __cppmlhdfe_bs_V, _b, _se, _bias, _z0, _ci_normal,
 * _ci_percentile, _ci_bc, _counts, and scalars __cppmlhdfe_bs_N_reps,
 * __cppmlhdfe_bs_N_misreps. */
ST_retcode cppmlhdfe_boot_finish(const char *args);

/* wildbootstrap(): the score bootstrap (plan 8.4) on the retained fit.
 *
 * wild_setup [vce] [boot]: the varlist holds the vce() cluster variable and
 * the bootstrap cluster variable, numbered here in ascending order (neither:
 * each row); matrices __cppmlhdfe_wb_test (tested regressor positions,
 * 1-based) and __cppmlhdfe_wb_scale (their analytic standard errors, the
 * scale of the interval search). Sets scalar __cppmlhdfe_wb_C.
 *
 * wild_chunk first=# count=# total=# [enum]: the weights of replicates
 * first .. first+count-1 (C per replicate, packed by `_ctools_hdfe_boot
 * wildweights`), or none with enum (Rademacher enumeration).
 *
 * wild_finish level=# ptype=# [ci] [keep]: per tested coefficient, the fit
 * restricted by H0: b = 0 and its p-value into __cppmlhdfe_wb_p (1 x T),
 * __cppmlhdfe_wb_fail (1 x T: the restricted fit failed); with ci, the
 * confidence set by inverting the test (a restricted fit per trial value)
 * into __cppmlhdfe_wb_ci (T x 2) and __cppmlhdfe_wb_cistat (T x 4);
 * scalars __cppmlhdfe_wb_reps and __cppmlhdfe_wb_misreps. keep holds the
 * fit for the pairs bootstrap. */
ST_retcode cppmlhdfe_wild_setup(const char *args);
ST_retcode cppmlhdfe_wild_chunk(const char *args);
ST_retcode cppmlhdfe_wild_finish(const char *args);

#endif /* CPPMLHDFE_BOOT_H */

/*
 * creghdfe_boot.h
 *
 * vce(bootstrap) for creghdfe: the model retained after the full-sample fit
 * and the plugin calls that run replicates on it
 * (docs/development/BOOTSTRAP_PLAN.md, section 7).
 * Part of the ctools Stata plugin suite
 */

#ifndef CREGHDFE_BOOT_H
#define CREGHDFE_BOOT_H

#include "creghdfe_types.h"
#include "creghdfe_input.h"
#include "../ctools_types.h"

/* What the full-sample fit hands over (borrowed; retain copies it). */
typedef struct creghdfe_boot_source_s
{
    ST_int N, K, G;                 /* compact rows; y + regressors; effects */
    const ST_double *raw;           /* N x K unstandardized columns */
    const ST_double *tilde;         /* N x K partialled columns, in the units partialled */
    const ST_double *stdevs;        /* K scales (all 1 without standardize) */
    const ST_double *weights;       /* N weights as the fit used them, or NULL */
    ST_int weight_type;             /* 1 aw, 2 fw, 3 pw */
    const FE_Factor *factors;       /* G factors on the compact rows */
    const ST_int *is_collinear;     /* K - 1 full-sample flags */
    const ST_double *b_hat;         /* K_keep + 1 coefficients, _cons last */
    ST_int K_keep;
    ST_int maxiter, drop_singletons, use_quad;
    ST_double tolerance, collinear_tol;
    const ST_int *mask;             /* N_orig loaded rows, 1 for a compact row */
    const perm_idx_t *obs_map;      /* N_orig 1-based Stata observations */
    ST_int N_orig;
    ST_int unit_var, strata_var;    /* plugin variable indices, or 0 */
    const creghdfe_input *in;       /* the fit's input: its slope factors
                                       (slope_idx, fe_absvar, num_slopes) */
} creghdfe_boot_source;

/* Keep the model for boot_chunk; reads the resampling unit and strata
 * variables (a cluster is a unit within each stratum, as in bsample).
 * Returns 0 or 920. Sets the
 * scalars __creghdfe_bs_U (units) and __creghdfe_bs_S (strata). */
ST_retcode creghdfe_boot_retain(const creghdfe_boot_source *src);

/* Free the retained model (idempotent). */
void creghdfe_boot_release(void);

/* Install creghdfe_boot_retain and _release as the hooks the fit and its
 * cleanup call (creghdfe_hdfe.h). */
void creghdfe_boot_install(void);

/* boot_chunk first=# count=# total=# (seedhi=# seedlo=# | stata | supplied)
 * [fe] [dse] [ci] [rows]: run replicates first .. first+count-1. With
 * stata, the varlist holds bsample's draw and tie-break uniforms, packed
 * by `_ctools_hdfe_boot uniforms` (hdfe_boot_read_draws); with supplied,
 * count _bsdraws() multiplicity variables. fe, dse and ci (read at
 * first=0) recover each replicate's effects, d, and keep the effects for
 * percentile bounds. rows stores the chunk's coefficients in matrix
 * __creghdfe_bs_rows (count x K_keep + 1) for saving(). Leaves the
 * replicates' statuses in matrix __creghdfe_bs_status (1 x count). */
ST_retcode creghdfe_boot_chunk(const char *args);

/* boot_finish level=# [coef] [fe] [nreps] [ci] [dse] [mse]: with fe, the
 * effects' SEs (then counts with nreps, bounds with ci, d's SE with dse)
 * into the call's variables; with coef, the summary into matrices
 * __creghdfe_bs_V,
 * _b, _se, _bias, _z0, _ci_normal, _ci_percentile, _ci_bc, _status, and
 * scalars __creghdfe_bs_N_reps, __creghdfe_bs_N_misreps; the caller sizes
 * the matrices (K_keep + 1 columns). */
ST_retcode creghdfe_boot_finish(const char *args);

/* wildbootstrap() (plan 7.5), on the retained fit:
 *
 * wild_setup [vce] [boot] [nonull] [partialled] [projected]: the varlist
 * holds the vce() cluster variable (vce) and the bootstrap cluster variable
 * (boot), numbered here in ascending order; the tested regressors are the
 * 1-based positions in matrix __creghdfe_wb_test. Without either, each row
 * is a cluster. Sets scalars __creghdfe_wb_C (bootstrap clusters),
 * __creghdfe_wb_G and __creghdfe_wb_path (0 nested, 1 projected,
 * 2 partialled, 3 projected although nested).
 *
 * wild_chunk first=# count=# total=# [enum]: replicates first ..
 * first+count-1 of total; the varlist holds their weights, C per replicate
 * packed by `_ctools_hdfe_boot wildweights`, or none with enum (Rademacher
 * enumeration).
 *
 * wild_finish level=# ptype=# [ci] [keep]: p-values of H0: b = 0 into
 * matrix __creghdfe_wb_p (1 x T); with ci, confidence sets into
 * __creghdfe_wb_ci (T x 2) and __creghdfe_wb_cistat (T x 4: open lower,
 * open upper, disjoint, empty); scalars __creghdfe_wb_reps and
 * __creghdfe_wb_misreps. keep holds the fit for fevce(bootstrap). */
ST_retcode creghdfe_wild_setup(const char *args);
ST_retcode creghdfe_wild_chunk(const char *args);
ST_retcode creghdfe_wild_finish(const char *args);

#endif /* CREGHDFE_BOOT_H */

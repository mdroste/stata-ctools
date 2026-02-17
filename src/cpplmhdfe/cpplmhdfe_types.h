/*
 * cpplmhdfe_types.h
 *
 * Type definitions for cpplmhdfe (Poisson PML with HDFE)
 * Part of the ctools Stata plugin suite
 */

#ifndef CPPLMHDFE_TYPES_H
#define CPPLMHDFE_TYPES_H

#include "stplugin.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Shared FE_Factor and HDFE_State types */
#include "../ctools_hdfe_utils.h"

/* Backward compatibility macros */
#define cholesky ctools_cholesky
#define invert_from_cholesky ctools_invert_from_cholesky

/* Compiler hints */
#if defined(__GNUC__) || defined(__clang__)
#define RESTRICT __restrict__
#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define RESTRICT
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/* ========================================================================
 * PPML State structure
 * ======================================================================== */

typedef struct {
    ST_double *mu;              /* Fitted values exp(eta), N */
    ST_double *eta;             /* Linear predictor Xb + FE + offset, N */
    ST_double *irls_weights;    /* mu * w_user, N */
    ST_double *working_dep;     /* z = eta + (y - mu)/mu, N */
    ST_double *offset;          /* Offset/exposure, N (or NULL) */
    ST_double *w_user;          /* User weights, N (or NULL) */
    ST_double *score_resid;     /* Score residuals (y - mu), N */
    ST_double deviance;
    ST_double ll;               /* Log-likelihood */
    ST_int max_irls_iter;
    ST_double irls_tol;
    ST_double sep_tol;          /* Separation tolerance */
    ST_int num_separated;
    ST_int verbose;
    ST_int converged;
    ST_int iterations;
} PPML_State;

/* ========================================================================
 * Factor data structure for initialization (mirrors creghdfe)
 * ======================================================================== */

typedef struct {
    ST_int *levels;      /* Level assignment for each obs (1-indexed) */
    ST_int *counts;      /* Count of obs per level */
    ST_int num_levels;   /* Number of unique levels */
    ST_int num_obs;      /* Number of observations */
} PPMLFactorData;

#endif /* CPPLMHDFE_TYPES_H */

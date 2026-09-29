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

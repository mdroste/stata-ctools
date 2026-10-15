/*
 * creghdfe_types.h
 *
 * Type definitions, structures, and macros for creghdfe plugin
 * Part of the ctools Stata plugin suite
 *
 * FE_Factor and HDFE_State are defined in the shared ctools_hdfe_utils.h
 */

#ifndef CREGHDFE_TYPES_H
#define CREGHDFE_TYPES_H

#include "stplugin.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* Shared FE_Factor and HDFE_State types */
#include "../ctools_hdfe_utils.h"
#include "../hdfe/hdfe_sample.h"

/* Thread limit is managed centrally via ctools_get_max_threads() in ctools_config.h */

/* Compiler hints for optimization */
#if defined(__GNUC__) || defined(__clang__)
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

/* Global state pointer - defined in creghdfe_hdfe.c */


/* ========================================================================
 * Factor data structure for initialization
 * ======================================================================== */

typedef HDFE_LoadedFactor FactorData;

#endif /* CREGHDFE_TYPES_H */

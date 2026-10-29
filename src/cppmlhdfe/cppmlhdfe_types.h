/*
 * cppmlhdfe_types.h
 *
 * Type definitions for cppmlhdfe (Poisson PML with HDFE)
 * Part of the ctools Stata plugin suite
 */

#ifndef CPPMLHDFE_TYPES_H
#define CPPMLHDFE_TYPES_H

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
#include "../hdfe/hdfe_sample.h"

/* ========================================================================
 * Factor data structure for initialization (shared with creghdfe)
 * ======================================================================== */

typedef HDFE_LoadedFactor PPMLFactorData;

#endif /* CPPMLHDFE_TYPES_H */

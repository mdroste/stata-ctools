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

/* Thread limit is managed centrally via ctools_get_max_threads() in ctools_config.h */

/* Compiler hints for optimization */
#if defined(__GNUC__) || defined(__clang__)
#define RESTRICT __restrict__
#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define RESTRICT
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/* Cache-friendly block size for matrix operations */
#define BLOCK_SIZE 512

/* Global state pointer - defined in creghdfe_hdfe.c */
extern HDFE_State *g_state;

/* ========================================================================
 * Hash table for factor creation
 * Simple open-addressing hash table for integer keys
 * ======================================================================== */

#define HASH_EMPTY -1
#define HASH_LOAD_FACTOR 0.7

typedef struct {
    ST_int *keys;      /* Original values */
    ST_int *values;    /* Level assignments (1-indexed) */
    ST_int capacity;
    ST_int size;
} IntHashTable;

/* ========================================================================
 * Factor data structure for initialization
 * ======================================================================== */

typedef struct {
    ST_int *levels;      /* Level assignment for each obs (1-indexed) */
    ST_int *counts;      /* Count of obs per level */
    ST_int num_levels;   /* Number of unique levels */
    ST_int num_obs;      /* Number of observations */
} FactorData;

#endif /* CREGHDFE_TYPES_H */

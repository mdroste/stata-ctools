/*
 * creghdfe_types.h
 *
 * Type definitions, structures, and macros for creghdfe plugin
 * Part of the ctools Stata plugin suite
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
#define OMP_NUM_THREADS 8
#endif

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

/* ========================================================================
 * Fixed Effect Factor structure
 * ======================================================================== */

typedef struct {
    ST_int num_levels;
    ST_int has_intercept;
    ST_int num_slopes;
    ST_int *levels;
    ST_double *counts;
    ST_double *means;
} FE_Factor;

/* ========================================================================
 * Global HDFE State structure
 * Holds all state needed for the CG solver and factor operations
 * ======================================================================== */

typedef struct {
    ST_int G;
    ST_int N;
    ST_int K;
    ST_int in1;
    ST_int in2;
    ST_int has_weights;
    ST_double *weights;
    FE_Factor *factors;
    ST_int maxiter;
    ST_double tolerance;
    ST_int verbose;
    /* Per-thread working buffers for parallel column processing */
    ST_double **thread_cg_r;
    ST_double **thread_cg_u;
    ST_double **thread_cg_v;
    ST_double **thread_proj;
    ST_double **thread_fe_means;  /* Per-thread means buffers for each FE */
    ST_int num_threads;
    /* Cached data from HDFE init for reuse in partial_out */
    ST_int factors_initialized;  /* Flag indicating factors are ready */
    ST_int df_a;                 /* Degrees of freedom absorbed */
    ST_int mobility_groups;      /* Number of mobility groups (connected components) */
} HDFE_State;

/* Global state pointer - defined in creghdfe_hdfe.c */
extern HDFE_State *g_state;

/* ========================================================================
 * Union-Find (Disjoint Set Union) for connected components
 * Used to efficiently count mobility groups in bipartite graphs
 * ======================================================================== */

typedef struct {
    ST_int *parent;
    ST_int *rank;
    ST_int size;
} UnionFind;

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

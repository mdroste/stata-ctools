#ifndef CPSMATCH_INTERNAL_H
#define CPSMATCH_INTERNAL_H
/*
 * cpsmatch_impl.c
 *
 * C-accelerated propensity score matching for Stata
 * Part of the ctools suite
 *
 * High-performance implementation of propensity score matching algorithms:
 * - Parallel nearest neighbor search over controls sorted by propensity score
 * - Radius/caliper matching
 * - Kernel matching with various kernel functions
 *
 * Algorithm Overview:
 * 1. Load treatment indicator and propensity scores from Stata using ctools_data_load()
 * 2. Sort control observations by propensity score
 * 3. For each treated observation, find matches using specified method
 * 4. Compute matching weights
 * 5. Store results back to Stata (_weight, _id, _support, etc.)
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Checked SPI writes resolve positions in the plugin varlist, even after
 * estimation has introduced internal Stata variables. */
#define SD_SAFEMODE
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_runtime.h"
#include "ctools_config.h"
#include "ctools_threads.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CPSMATCH_MAX_NEIGHBORS 100
#define CPSMATCH_DEFAULT_CALIPER 0.25 /* Default: 0.25 * SD of pscore */

/* Matching methods */
typedef enum
{
    MATCH_NEAREST = 0, /* Nearest neighbor matching */
    MATCH_RADIUS = 1,  /* Radius/caliper matching */
    MATCH_KERNEL = 2   /* Kernel matching */
} match_method_t;

/* Kernel types for kernel matching */
typedef enum
{
    KERNEL_EPAN = 0,     /* Epanechnikov kernel */
    KERNEL_NORMAL = 1,   /* Gaussian/normal kernel */
    KERNEL_BIWEIGHT = 2, /* Biweight kernel */
    KERNEL_UNIFORM = 3,  /* Uniform kernel */
    KERNEL_TRICUBE = 4   /* Tricube kernel */
} kernel_type_t;

/* Matching options structure */
typedef struct
{
    int treat_idx;         /* 1-based index of treatment variable */
    int pscore_idx;        /* 1-based index of propensity score variable */
    int outcome_idx;       /* 1-based index of outcome variable (0 if none) */
    int out_weight_idx;    /* 1-based index for output weight variable */
    int out_nn_idx;        /* actual eligible controls used for each treated row */
    int out_match_idx;     /* 1-based index for output match ID variable */
    int out_support_idx;   /* 1-based index for common support indicator */
    size_t nobs;           /* Number of observations */
    match_method_t method; /* Matching method */
    int n_neighbors;       /* Number of neighbors for NN matching */
    double caliper;        /* Caliper width (in SD units or absolute) */
    int common_support;    /* Impose common support */
    kernel_type_t kernel;  /* Kernel type for kernel matching */
    double bandwidth;      /* Bandwidth for kernel matching */
    int with_replace;      /* Match with replacement */
    int ties;              /* Handle ties (use all tied matches) */
    int descending;        /* Process in descending pscore order */
    int verbose;           /* Verbose output */
    int touse;             /* Last plugin variable is the sample marker */
    int statistics;        /* Stored _treated follows _nn, before the marker */
} cpsmatch_options;

/* Match result for a single observation */
typedef struct
{
    size_t *match_ids;     /* Array of matched control indices (0-based) */
    double *match_weights; /* Array of weights for matches */
    int n_matches;         /* Number of matches found */
    int on_support;        /* 1 if on common support, 0 otherwise */
} match_result;

/* ============================================================================
 * Kernel Functions
 * ============================================================================ */

/* ============================================================================
 * Comparison functions for sorting
 * ============================================================================ */

typedef struct
{
    size_t idx;
    double dist;
} idx_dist_pair;

/* ============================================================================
 * Per-thread reusable workspace to eliminate per-match malloc/free
 * ============================================================================ */

#define MATCH_WORKSPACE_INITIAL_CAP 256

typedef struct
{
    idx_dist_pair *candidates; /* NN: candidate pairs during search */
    size_t candidates_cap;
    size_t *match_ids;     /* Output: matched control indices */
    double *match_weights; /* Output: match weights */
    size_t match_cap;
} match_workspace;

/* Structure for sorted control index */
typedef struct
{
    size_t orig_idx; /* Original index in control arrays */
    size_t obs_idx;  /* Original Stata observation index (for stable sort) */
    double pscore;   /* Propensity score for binary search */
} sorted_control_t;

/* Sort by pscore, then by observation index for stability (matches psmatch2 behavior) */

/* Binary search to find leftmost index where pscore >= target */

/* Binary search to find rightmost index where pscore <= target */

/* ============================================================================
 * Argument Parsing
 * ============================================================================ */

/* ============================================================================
 * Nearest Neighbor Matching (Original O(n) version for noreplacement)
 * Note: Currently unused - replaced by binary search version.
 *       Retained for potential future use with different matching algorithms.
 * ============================================================================ */

/* ============================================================================
 * FAST Nearest Neighbor Matching using sorted index (O(log n) per treated)
 * ============================================================================ */

/* ============================================================================
 * FAST Radius Matching - returns range info only, no memory allocation
 * ============================================================================ */

typedef struct
{
    size_t start;   /* Start index in sorted_controls */
    size_t end;     /* End index (exclusive) in sorted_controls */
    int n_matches;  /* Number of matches = end - start */
    int on_support; /* 1 if any matches found */
} radius_match_info;

/* ============================================================================
 * FAST Kernel Matching using sorted index
 * ============================================================================ */

/* ============================================================================
 * Kernel Matching - currently unused, retained for future kernel matching
 * ============================================================================ */

/* ============================================================================
 * Radius Matching - currently unused, replaced by optimized binary search version
 * ============================================================================ */

/* ============================================================================
 * Main Entry Point
 * ============================================================================ */

/* Private kernels: caller owns samples, workspaces and output buffers.
 * Matching/window traversal preserve stable order and existing definitions. */
double cpsmatch_kernel_epan(double u);
double cpsmatch_kernel_normal(double u);
double cpsmatch_kernel_biweight(double u);
double cpsmatch_kernel_uniform(double u);
double cpsmatch_kernel_tricube(double u);
double cpsmatch_apply_kernel(double u, kernel_type_t kernel);
void cpsmatch_match_workspace_init(match_workspace *ws, size_t initial_cap);
void cpsmatch_match_workspace_free(match_workspace *ws);
int cpsmatch_match_workspace_ensure_candidates(match_workspace *ws, size_t needed);
int cpsmatch_match_workspace_ensure_match(match_workspace *ws, size_t needed);
int cpsmatch_compare_by_pscore(const void *a, const void *b);
int cpsmatch_compare_by_pscore_desc(const void *a, const void *b);
size_t cpsmatch_binary_search_left(const sorted_control_t *sorted, size_t n, double target);
size_t cpsmatch_binary_search_right(const sorted_control_t *sorted, size_t n, double target);
int cpsmatch_parse_options(const char *args, cpsmatch_options *opts);
/* Radius outputs use filtered-row indices; controls and endpoints are sorted.
 * No host callbacks occur here. Allocation failures return 920 before writes. */
int cpsmatch_match_radius(const sorted_control_t *controls, size_t n_controls,
                         const double *treated_pscore, const size_t *treated_idx, size_t n_treated,
                         const size_t *control_idx, const perm_idx_t *obs_map,
                         double radius, int common, double common_min, double common_max,
                         double *out_weight, double *out_match, double *out_nn, double *out_support,
                         int *n_matched, int *n_off);
int cpsmatch_find_nearest_neighbors_fast(double pscore_treated,
                                         const sorted_control_t *sorted_controls,
                                         size_t *control_indices, size_t n_controls,
                                         int n_neighbors, double caliper, int handle_ties,
                                         match_workspace *ws, match_result *result);
void cpsmatch_find_radius_range_fast(double pscore_treated, const sorted_control_t *sorted_controls,
                                     size_t n_controls, double radius, radius_match_info *info);
int cpsmatch_find_kernel_matches_fast(double pscore_treated,
                                      const sorted_control_t *sorted_controls,
                                      size_t *control_indices, size_t n_controls, double bandwidth,
                                      kernel_type_t kernel, double caliper, match_workspace *ws,
                                      match_result *result);
#endif

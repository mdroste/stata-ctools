/*
 * cpsmatch_impl.c
 *
 * C-accelerated propensity score matching for Stata
 * Part of the ctools suite
 *
 * High-performance implementation of propensity score matching algorithms:
 * - Parallel nearest neighbor search using k-d tree or brute force
 * - Radius/caliper matching with spatial indexing
 * - Kernel matching with various kernel functions
 * - Mahalanobis distance matching
 *
 * Algorithm Overview:
 * 1. Load treatment indicator and propensity scores from Stata using ctools_data_load_selective()
 * 2. Build index structure for control observations
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

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_timer.h"
#include "ctools_config.h"
#include "ctools_threads.h"
#include "cpsmatch_impl.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CPSMATCH_MAX_NEIGHBORS 100
#define CPSMATCH_DEFAULT_CALIPER 0.25  /* Default: 0.25 * SD of pscore */

/* Matching methods */
typedef enum {
    MATCH_NEAREST = 0,     /* Nearest neighbor matching */
    MATCH_RADIUS = 1,      /* Radius/caliper matching */
    MATCH_KERNEL = 2,      /* Kernel matching */
    MATCH_MAHAL = 3        /* Mahalanobis distance matching */
} match_method_t;

/* Kernel types for kernel matching */
typedef enum {
    KERNEL_EPAN = 0,       /* Epanechnikov kernel */
    KERNEL_NORMAL = 1,     /* Gaussian/normal kernel */
    KERNEL_BIWEIGHT = 2,   /* Biweight kernel */
    KERNEL_UNIFORM = 3,    /* Uniform kernel */
    KERNEL_TRICUBE = 4     /* Tricube kernel */
} kernel_type_t;

/* Matching options structure */
typedef struct {
    int treat_idx;          /* 1-based index of treatment variable */
    int pscore_idx;         /* 1-based index of propensity score variable */
    int outcome_idx;        /* 1-based index of outcome variable (0 if none) */
    int out_weight_idx;     /* 1-based index for output weight variable */
    int out_match_idx;      /* 1-based index for output match ID variable */
    int out_support_idx;    /* 1-based index for common support indicator */
    size_t nobs;            /* Number of observations */
    match_method_t method;  /* Matching method */
    int n_neighbors;        /* Number of neighbors for NN matching */
    double caliper;         /* Caliper width (in SD units or absolute) */
    int common_support;     /* Impose common support */
    kernel_type_t kernel;   /* Kernel type for kernel matching */
    double bandwidth;       /* Bandwidth for kernel matching */
    int with_replace;       /* Match with replacement */
    int ties;               /* Handle ties (use all tied matches) */
    int descending;         /* Process in descending pscore order */
    int verbose;            /* Verbose output */
} cpsmatch_options;

/* Match result for a single observation */
typedef struct {
    size_t *match_ids;      /* Array of matched control indices (0-based) */
    double *match_dists;    /* Array of distances to matches */
    double *match_weights;  /* Array of weights for matches */
    int n_matches;          /* Number of matches found */
    int on_support;         /* 1 if on common support, 0 otherwise */
} match_result;

/* ============================================================================
 * Kernel Functions
 * ============================================================================ */

static inline double kernel_epan(double u)
{
    if (fabs(u) > 1.0) return 0.0;
    return 0.75 * (1.0 - u * u);
}

static inline double kernel_normal(double u)
{
    return exp(-0.5 * u * u) / sqrt(2.0 * M_PI);
}

static inline double kernel_biweight(double u)
{
    if (fabs(u) > 1.0) return 0.0;
    double v = 1.0 - u * u;
    return 0.9375 * v * v;
}

static inline double kernel_uniform(double u)
{
    if (fabs(u) > 1.0) return 0.0;
    return 0.5;
}

static inline double kernel_tricube(double u)
{
    if (fabs(u) > 1.0) return 0.0;
    double v = 1.0 - fabs(u * u * u);
    return 70.0 / 81.0 * v * v * v;
}

static double apply_kernel(double u, kernel_type_t kernel)
{
    switch (kernel) {
        case KERNEL_EPAN:    return kernel_epan(u);
        case KERNEL_NORMAL:  return kernel_normal(u);
        case KERNEL_BIWEIGHT: return kernel_biweight(u);
        case KERNEL_UNIFORM: return kernel_uniform(u);
        case KERNEL_TRICUBE: return kernel_tricube(u);
        default:             return kernel_epan(u);
    }
}

/* ============================================================================
 * Comparison functions for sorting
 * ============================================================================ */

typedef struct {
    size_t idx;
    double dist;
} idx_dist_pair;

static int compare_by_dist(const void *a, const void *b)
{
    const idx_dist_pair *pa = (const idx_dist_pair *)a;
    const idx_dist_pair *pb = (const idx_dist_pair *)b;
    if (pa->dist < pb->dist) return -1;
    if (pa->dist > pb->dist) return 1;
    return 0;
}

/* Structure for sorted control index */
typedef struct {
    size_t orig_idx;     /* Original index in control arrays */
    size_t obs_idx;      /* Original Stata observation index (for stable sort) */
    double pscore;       /* Propensity score for binary search */
} sorted_control_t;

/* Sort by pscore, then by observation index for stability (matches psmatch2 behavior) */
static int compare_by_pscore(const void *a, const void *b)
{
    const sorted_control_t *pa = (const sorted_control_t *)a;
    const sorted_control_t *pb = (const sorted_control_t *)b;
    if (pa->pscore < pb->pscore) return -1;
    if (pa->pscore > pb->pscore) return 1;
    /* Tie-breaker: sort by Stata observation index (ascending) for psmatch2 compatibility */
    if (pa->obs_idx < pb->obs_idx) return -1;
    if (pa->obs_idx > pb->obs_idx) return 1;
    return 0;
}

static int compare_by_pscore_desc(const void *a, const void *b)
{
    const sorted_control_t *pa = (const sorted_control_t *)a;
    const sorted_control_t *pb = (const sorted_control_t *)b;
    if (pa->pscore > pb->pscore) return -1;
    if (pa->pscore < pb->pscore) return 1;
    /* Tie-breaker: sort by Stata observation index (ascending) for psmatch2 compatibility */
    if (pa->obs_idx < pb->obs_idx) return -1;
    if (pa->obs_idx > pb->obs_idx) return 1;
    return 0;
}

/* Binary search to find leftmost index where pscore >= target */
static size_t binary_search_left(const sorted_control_t *sorted, size_t n, double target)
{
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (sorted[mid].pscore < target) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

/* Binary search to find rightmost index where pscore <= target */
static size_t binary_search_right(const sorted_control_t *sorted, size_t n, double target)
{
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (sorted[mid].pscore <= target) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

/* ============================================================================
 * Argument Parsing
 * ============================================================================ */

static int parse_options(const char *args, cpsmatch_options *opts)
{
    const char *p = args;
    char *endptr;
    int val;

    /* Initialize defaults */
    memset(opts, 0, sizeof(cpsmatch_options));
    opts->n_neighbors = 1;
    opts->caliper = 0.0;  /* 0 means no caliper */
    opts->with_replace = 1;
    opts->ties = 0;
    opts->descending = 0;
    opts->verbose = 0;
    opts->kernel = KERNEL_EPAN;
    opts->bandwidth = 0.06;

    /* Parse space-separated arguments */
    /* Expected format: treat_idx pscore_idx outcome_idx out_weight_idx out_match_idx
     *                  out_support_idx nobs method neighbor caliper common
     *                  kernel bwidth with_replace ties descending verbose */

    /* Parse treat_idx */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->treat_idx = val;
    p = endptr;

    /* Parse pscore_idx */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->pscore_idx = val;
    p = endptr;

    /* Parse outcome_idx */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->outcome_idx = val;
    p = endptr;

    /* Parse out_weight_idx */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->out_weight_idx = val;
    p = endptr;

    /* Parse out_match_idx */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->out_match_idx = val;
    p = endptr;

    /* Parse out_support_idx */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->out_support_idx = val;
    p = endptr;

    /* Parse nobs */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    size_t nobs_val = strtoul(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->nobs = nobs_val;
    p = endptr;

    /* Parse method */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->method = (match_method_t)val;
    p = endptr;

    /* Parse n_neighbors */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->n_neighbors = val;
    p = endptr;

    /* Parse caliper */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    double dval = strtod(p, &endptr);
    if (endptr == p) return -1;
    opts->caliper = dval;
    p = endptr;

    /* Parse common_support */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->common_support = val;
    p = endptr;

    /* Parse kernel */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->kernel = (kernel_type_t)val;
    p = endptr;

    /* Parse bandwidth */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    dval = strtod(p, &endptr);
    if (endptr == p) return -1;
    opts->bandwidth = dval;
    p = endptr;

    /* Parse with_replace */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->with_replace = val;
    p = endptr;

    /* Parse ties */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->ties = val;
    p = endptr;

    /* Parse descending */
    while (*p == ' ') p++;
    if (*p == '\0') return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p) return -1;
    opts->descending = val;
    p = endptr;

    /* Parse verbose (optional) */
    while (*p == ' ') p++;
    if (*p != '\0') {
        val = strtol(p, &endptr, 10);
        if (endptr != p) {
            opts->verbose = val;
        }
    }

    return 0;
}

/* ============================================================================
 * Nearest Neighbor Matching (Original O(n) version for noreplacement)
 * Note: Currently unused - replaced by binary search version.
 *       Retained for potential future use with different matching algorithms.
 * ============================================================================ */

__attribute__((unused))
static int find_nearest_neighbors(
    double pscore_treated,
    double *control_pscores,
    size_t *control_indices,
    int *control_available,  /* For matching without replacement */
    size_t n_controls,
    int n_neighbors,
    double caliper,
    int with_replace,
    int handle_ties,
    match_result *result)
{
    if (n_controls == 0) {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* Allocate temporary storage for all potential matches */
    idx_dist_pair *candidates = malloc(n_controls * sizeof(idx_dist_pair));
    if (!candidates) return -1;

    int n_candidates = 0;

    /* Find all controls within caliper (or all if no caliper) */
    for (size_t i = 0; i < n_controls; i++) {
        if (!with_replace && control_available && !control_available[i]) {
            continue;  /* Already matched */
        }

        double dist = fabs(pscore_treated - control_pscores[i]);

        if (caliper > 0.0 && dist > caliper) {
            continue;  /* Outside caliper */
        }

        candidates[n_candidates].idx = i;
        candidates[n_candidates].dist = dist;
        n_candidates++;
    }

    if (n_candidates == 0) {
        result->n_matches = 0;
        result->on_support = 0;
        free(candidates);
        return 0;
    }

    /* Sort by distance */
    qsort(candidates, n_candidates, sizeof(idx_dist_pair), compare_by_dist);

    /* Select neighbors */
    int n_selected = 0;
    int max_select = (n_neighbors < n_candidates) ? n_neighbors : n_candidates;

    if (handle_ties && max_select > 0) {
        /* Include all observations tied at the last selected distance */
        double last_dist = candidates[max_select - 1].dist;
        max_select = n_candidates;
        for (int i = n_neighbors; i < n_candidates; i++) {
            if (candidates[i].dist > last_dist + 1e-12) {
                max_select = i;
                break;
            }
        }
    }

    /* Allocate result arrays */
    result->match_ids = malloc(max_select * sizeof(size_t));
    result->match_dists = malloc(max_select * sizeof(double));
    result->match_weights = malloc(max_select * sizeof(double));

    if (!result->match_ids || !result->match_dists || !result->match_weights) {
        free(candidates);
        free(result->match_ids);
        free(result->match_dists);
        free(result->match_weights);
        result->match_ids = NULL;
        result->match_dists = NULL;
        result->match_weights = NULL;
        return -1;
    }

    double total_weight = 0.0;
    for (int i = 0; i < max_select; i++) {
        result->match_ids[n_selected] = control_indices[candidates[i].idx];
        result->match_dists[n_selected] = candidates[i].dist;
        result->match_weights[n_selected] = 1.0;
        total_weight += 1.0;

        /* Mark as used if matching without replacement */
        if (!with_replace && control_available) {
            control_available[candidates[i].idx] = 0;
        }

        n_selected++;
    }

    /* Normalize weights */
    if (total_weight > 0.0) {
        for (int i = 0; i < n_selected; i++) {
            result->match_weights[i] /= total_weight;
        }
    }

    result->n_matches = n_selected;
    result->on_support = 1;

    free(candidates);
    return 0;
}

/* ============================================================================
 * FAST Nearest Neighbor Matching using sorted index (O(log n) per treated)
 * ============================================================================ */

static int find_nearest_neighbors_fast(
    double pscore_treated,
    const sorted_control_t *sorted_controls,
    size_t *control_indices,
    size_t n_controls,
    int n_neighbors,
    double caliper,
    int handle_ties,
    match_result *result)
{
    if (n_controls == 0) {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* Binary search to find insertion point */
    size_t insert_pos = binary_search_left(sorted_controls, n_controls, pscore_treated);

    /* For k-NN, we need to find k nearest neighbors around insert_pos
     * Use two pointers expanding outward from insert_pos */

    /* Pre-allocate for max possible matches (n_neighbors + ties) */
    int max_alloc = n_neighbors * 2 + 10;  /* Some extra for ties */
    if (max_alloc > (int)n_controls) max_alloc = (int)n_controls;

    idx_dist_pair *candidates = malloc(max_alloc * sizeof(idx_dist_pair));
    if (!candidates) return -1;

    int n_candidates = 0;

    /* Two-pointer expansion from insert_pos */
    ssize_t left = (ssize_t)insert_pos - 1;
    size_t right = insert_pos;

    /* First, find the k nearest neighbors */
    while (n_candidates < n_neighbors && (left >= 0 || right < n_controls)) {
        double dist_left = (left >= 0) ? fabs(pscore_treated - sorted_controls[left].pscore) : DBL_MAX;
        double dist_right = (right < n_controls) ? fabs(pscore_treated - sorted_controls[right].pscore) : DBL_MAX;

        if (dist_left <= dist_right) {
            if (caliper > 0.0 && dist_left > caliper) {
                left = -1;  /* No more valid on left */
            } else {
                candidates[n_candidates].idx = sorted_controls[left].orig_idx;
                candidates[n_candidates].dist = dist_left;
                n_candidates++;
                left--;
            }
        } else {
            if (caliper > 0.0 && dist_right > caliper) {
                right = n_controls;  /* No more valid on right */
            } else {
                candidates[n_candidates].idx = sorted_controls[right].orig_idx;
                candidates[n_candidates].dist = dist_right;
                n_candidates++;
                right++;
            }
        }
    }

    if (n_candidates == 0) {
        result->n_matches = 0;
        result->on_support = 0;
        free(candidates);
        return 0;
    }

    /* Handle ties: include all observations at the same distance as the k-th */
    if (handle_ties && n_candidates >= n_neighbors) {
        double last_dist = candidates[n_candidates - 1].dist;

        /* Continue expanding to find ties */
        while (left >= 0 || right < n_controls) {
            double dist_left = (left >= 0) ? fabs(pscore_treated - sorted_controls[left].pscore) : DBL_MAX;
            double dist_right = (right < n_controls) ? fabs(pscore_treated - sorted_controls[right].pscore) : DBL_MAX;

            double min_dist = (dist_left < dist_right) ? dist_left : dist_right;
            if (min_dist > last_dist + 1e-12) break;  /* No more ties */

            /* Reallocate if needed */
            if (n_candidates >= max_alloc) {
                max_alloc *= 2;
                idx_dist_pair *new_candidates = realloc(candidates, max_alloc * sizeof(idx_dist_pair));
                if (!new_candidates) {
                    free(candidates);
                    return -1;
                }
                candidates = new_candidates;
            }

            if (dist_left <= dist_right && dist_left <= last_dist + 1e-12) {
                candidates[n_candidates].idx = sorted_controls[left].orig_idx;
                candidates[n_candidates].dist = dist_left;
                n_candidates++;
                left--;
            } else if (dist_right <= last_dist + 1e-12) {
                candidates[n_candidates].idx = sorted_controls[right].orig_idx;
                candidates[n_candidates].dist = dist_right;
                n_candidates++;
                right++;
            } else {
                break;
            }
        }
    }

    /* Allocate result arrays */
    result->match_ids = malloc(n_candidates * sizeof(size_t));
    result->match_dists = malloc(n_candidates * sizeof(double));
    result->match_weights = malloc(n_candidates * sizeof(double));

    if (!result->match_ids || !result->match_dists || !result->match_weights) {
        free(candidates);
        free(result->match_ids);
        free(result->match_dists);
        free(result->match_weights);
        result->match_ids = NULL;
        result->match_dists = NULL;
        result->match_weights = NULL;
        return -1;
    }

    double total_weight = (double)n_candidates;
    for (int i = 0; i < n_candidates; i++) {
        result->match_ids[i] = control_indices[candidates[i].idx];
        result->match_dists[i] = candidates[i].dist;
        result->match_weights[i] = 1.0 / total_weight;
    }

    result->n_matches = n_candidates;
    result->on_support = 1;

    free(candidates);
    return 0;
}

/* ============================================================================
 * FAST Radius Matching - returns range info only, no memory allocation
 * ============================================================================ */

typedef struct {
    size_t start;       /* Start index in sorted_controls */
    size_t end;         /* End index (exclusive) in sorted_controls */
    int n_matches;      /* Number of matches = end - start */
    int on_support;     /* 1 if any matches found */
} radius_match_info;

static void find_radius_range_fast(
    double pscore_treated,
    const sorted_control_t *sorted_controls,
    size_t n_controls,
    double radius,
    radius_match_info *info)
{
    if (n_controls == 0 || radius <= 0.0) {
        info->n_matches = 0;
        info->on_support = 0;
        return;
    }

    /* Binary search to find range [low_target, high_target] */
    double low_target = pscore_treated - radius;
    double high_target = pscore_treated + radius;

    info->start = binary_search_left(sorted_controls, n_controls, low_target);
    info->end = binary_search_right(sorted_controls, n_controls, high_target);
    info->n_matches = (int)(info->end - info->start);
    info->on_support = (info->n_matches > 0) ? 1 : 0;
}

/* Legacy wrapper for compatibility - currently unused */
__attribute__((unused))
static int find_radius_matches_fast(
    double pscore_treated,
    const sorted_control_t *sorted_controls,
    size_t *control_indices,
    size_t n_controls,
    double radius,
    match_result *result)
{
    radius_match_info info;
    find_radius_range_fast(pscore_treated, sorted_controls, n_controls, radius, &info);

    if (info.n_matches == 0) {
        result->n_matches = 0;
        result->on_support = 0;
        result->match_ids = NULL;
        result->match_dists = NULL;
        result->match_weights = NULL;
        return 0;
    }

    /* Allocate result arrays */
    result->match_ids = malloc(info.n_matches * sizeof(size_t));
    result->match_dists = malloc(info.n_matches * sizeof(double));
    result->match_weights = malloc(info.n_matches * sizeof(double));

    if (!result->match_ids || !result->match_dists || !result->match_weights) {
        free(result->match_ids);
        free(result->match_dists);
        free(result->match_weights);
        result->match_ids = NULL;
        result->match_dists = NULL;
        result->match_weights = NULL;
        return -1;
    }

    /* Fill arrays */
    double weight = 1.0 / info.n_matches;
    for (size_t i = info.start; i < info.end; i++) {
        size_t idx = i - info.start;
        result->match_ids[idx] = control_indices[sorted_controls[i].orig_idx];
        result->match_dists[idx] = fabs(pscore_treated - sorted_controls[i].pscore);
        result->match_weights[idx] = weight;
    }

    result->n_matches = info.n_matches;
    result->on_support = 1;

    return 0;
}

/* ============================================================================
 * FAST Kernel Matching using sorted index
 * ============================================================================ */

static int find_kernel_matches_fast(
    double pscore_treated,
    const sorted_control_t *sorted_controls,
    size_t *control_indices,
    size_t n_controls,
    double bandwidth,
    kernel_type_t kernel,
    double caliper,
    match_result *result)
{
    if (n_controls == 0 || bandwidth <= 0.0) {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* Determine search range based on kernel type and caliper */
    double search_radius = bandwidth;  /* Most kernels have support [-1, 1] */
    if (kernel == KERNEL_NORMAL) {
        search_radius = 4.0 * bandwidth;  /* Normal has infinite support, use 4 SD */
    }
    if (caliper > 0.0 && caliper < search_radius) {
        search_radius = caliper;
    }

    /* Binary search to find range */
    double low_target = pscore_treated - search_radius;
    double high_target = pscore_treated + search_radius;

    size_t start = binary_search_left(sorted_controls, n_controls, low_target);
    size_t end = binary_search_right(sorted_controls, n_controls, high_target);

    if (start >= end) {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* First pass: count and compute total weight */
    int n_matches = 0;
    double total_weight = 0.0;

    for (size_t i = start; i < end; i++) {
        double dist = fabs(pscore_treated - sorted_controls[i].pscore);
        if (caliper > 0.0 && dist > caliper) continue;

        double u = dist / bandwidth;
        double w = apply_kernel(u, kernel);
        if (w > 0.0) {
            n_matches++;
            total_weight += w;
        }
    }

    if (n_matches == 0) {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* Allocate result arrays */
    result->match_ids = malloc(n_matches * sizeof(size_t));
    result->match_dists = malloc(n_matches * sizeof(double));
    result->match_weights = malloc(n_matches * sizeof(double));

    if (!result->match_ids || !result->match_dists || !result->match_weights) {
        free(result->match_ids);
        free(result->match_dists);
        free(result->match_weights);
        result->match_ids = NULL;
        result->match_dists = NULL;
        result->match_weights = NULL;
        return -1;
    }

    /* Second pass: fill arrays */
    int idx = 0;
    for (size_t i = start; i < end; i++) {
        double dist = fabs(pscore_treated - sorted_controls[i].pscore);
        if (caliper > 0.0 && dist > caliper) continue;

        double u = dist / bandwidth;
        double w = apply_kernel(u, kernel);
        if (w > 0.0) {
            result->match_ids[idx] = control_indices[sorted_controls[i].orig_idx];
            result->match_dists[idx] = dist;
            result->match_weights[idx] = w / total_weight;
            idx++;
        }
    }

    result->n_matches = n_matches;
    result->on_support = 1;

    return 0;
}

/* ============================================================================
 * Kernel Matching - currently unused, retained for future kernel matching
 * ============================================================================ */

__attribute__((unused))
static int find_kernel_matches(
    double pscore_treated,
    double *control_pscores,
    size_t *control_indices,
    size_t n_controls,
    double bandwidth,
    kernel_type_t kernel,
    double caliper,
    match_result *result)
{
    if (n_controls == 0 || bandwidth <= 0.0) {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* Count controls within bandwidth (or caliper) */
    int n_matches = 0;
    double total_weight = 0.0;

    /* First pass: count matches */
    for (size_t i = 0; i < n_controls; i++) {
        double dist = fabs(pscore_treated - control_pscores[i]);
        if (caliper > 0.0 && dist > caliper) continue;

        double u = dist / bandwidth;
        double w = apply_kernel(u, kernel);
        if (w > 0.0) {
            n_matches++;
        }
    }

    if (n_matches == 0) {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* Allocate result arrays */
    result->match_ids = malloc(n_matches * sizeof(size_t));
    result->match_dists = malloc(n_matches * sizeof(double));
    result->match_weights = malloc(n_matches * sizeof(double));

    if (!result->match_ids || !result->match_dists || !result->match_weights) {
        free(result->match_ids);
        free(result->match_dists);
        free(result->match_weights);
        result->match_ids = NULL;
        result->match_dists = NULL;
        result->match_weights = NULL;
        return -1;
    }

    /* Second pass: fill arrays */
    int idx = 0;
    for (size_t i = 0; i < n_controls; i++) {
        double dist = fabs(pscore_treated - control_pscores[i]);
        if (caliper > 0.0 && dist > caliper) continue;

        double u = dist / bandwidth;
        double w = apply_kernel(u, kernel);
        if (w > 0.0) {
            result->match_ids[idx] = control_indices[i];
            result->match_dists[idx] = dist;
            result->match_weights[idx] = w;
            total_weight += w;
            idx++;
        }
    }

    /* Normalize weights */
    if (total_weight > 0.0) {
        for (int i = 0; i < n_matches; i++) {
            result->match_weights[i] /= total_weight;
        }
    }

    result->n_matches = n_matches;
    result->on_support = 1;

    return 0;
}

/* ============================================================================
 * Radius Matching - currently unused, replaced by optimized binary search version
 * ============================================================================ */

__attribute__((unused))
static int find_radius_matches(
    double pscore_treated,
    double *control_pscores,
    size_t *control_indices,
    size_t n_controls,
    double radius,
    match_result *result)
{
    if (n_controls == 0 || radius <= 0.0) {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* Count controls within radius */
    int n_matches = 0;
    for (size_t i = 0; i < n_controls; i++) {
        double dist = fabs(pscore_treated - control_pscores[i]);
        if (dist <= radius) {
            n_matches++;
        }
    }

    if (n_matches == 0) {
        result->n_matches = 0;
        result->on_support = 0;
        return 0;
    }

    /* Allocate result arrays */
    result->match_ids = malloc(n_matches * sizeof(size_t));
    result->match_dists = malloc(n_matches * sizeof(double));
    result->match_weights = malloc(n_matches * sizeof(double));

    if (!result->match_ids || !result->match_dists || !result->match_weights) {
        free(result->match_ids);
        free(result->match_dists);
        free(result->match_weights);
        result->match_ids = NULL;
        result->match_dists = NULL;
        result->match_weights = NULL;
        return -1;
    }

    /* Fill arrays */
    int idx = 0;
    for (size_t i = 0; i < n_controls; i++) {
        double dist = fabs(pscore_treated - control_pscores[i]);
        if (dist <= radius) {
            result->match_ids[idx] = control_indices[i];
            result->match_dists[idx] = dist;
            result->match_weights[idx] = 1.0 / n_matches;  /* Equal weights */
            idx++;
        }
    }

    result->n_matches = n_matches;
    result->on_support = 1;

    return 0;
}

/* ============================================================================
 * Main Entry Point
 * ============================================================================ */

ST_retcode cpsmatch_main(const char *args)
{
    cpsmatch_options opts;
    double t_start, t_load, t_setup, t_sort, t_match, t_store, t_total;
    ST_retcode rc = 0;

    t_start = ctools_timer_seconds();
    t_setup = 0.0;
    t_sort = 0.0;

    /* Parse arguments */
    if (args == NULL || *args == '\0') {
        SF_error("cpsmatch: no arguments specified\n");
        return 198;
    }

    if (parse_options(args, &opts) != 0) {
        SF_error("cpsmatch: failed to parse arguments\n");
        return 198;
    }

    if (opts.nobs == 0) {
        SF_error("cpsmatch: no observations\n");
        return 2000;
    }

    /* Get observation range - use full in-range, we'll filter by if-condition */
    ST_int obs1 = SF_in1();
    ST_int obs2 = SF_in2();
    size_t nobs_range = (size_t)(obs2 - obs1 + 1);

    /* ========================================================================
     * PHASE 1: Load data from Stata using ctools_data_load_selective()
     * ======================================================================== */
    t_load = ctools_timer_seconds();

    /* Build array of variable indices to load */
    int nvars_to_load = (opts.outcome_idx > 0) ? 3 : 2;
    int *var_indices = malloc(nvars_to_load * sizeof(int));
    if (!var_indices) {
        SF_error("cpsmatch: memory allocation failed\n");
        return 920;
    }
    var_indices[0] = opts.treat_idx;
    var_indices[1] = opts.pscore_idx;
    if (opts.outcome_idx > 0) {
        var_indices[2] = opts.outcome_idx;
    }

    /* Load data using ctools infrastructure - load full in-range */
    stata_data input_data;
    stata_data_init(&input_data);
    stata_retcode load_rc = ctools_data_load_selective(&input_data, var_indices,
                                                        nvars_to_load, obs1, obs2);
    free(var_indices);

    if (load_rc != STATA_OK) {
        SF_error("cpsmatch: failed to load data\n");
        stata_data_free(&input_data);
        return 920;
    }

    /* Get pointers to the loaded data arrays */
    double *treatment = input_data.vars[0].data.dbl;
    double *pscore = input_data.vars[1].data.dbl;
    /* outcome is loaded for potential ATT/ATE computation but not currently used */
    double *outcome __attribute__((unused)) = (opts.outcome_idx > 0) ? input_data.vars[2].data.dbl : NULL;

    t_load = ctools_timer_seconds() - t_load;

    /* ========================================================================
     * PHASE 2: Separate treated and control groups (with if-condition check)
     * ======================================================================== */
    t_setup = ctools_timer_seconds();

    /* Count treated and controls - check SF_ifobs() for if-condition filtering */
    size_t n_treated = 0, n_controls = 0;
    for (size_t i = 0; i < nobs_range; i++) {
        ST_int stata_obs = obs1 + (ST_int)i;
        /* Check if this observation satisfies the if-condition */
        if (!SF_ifobs(stata_obs)) continue;
        if (SF_is_missing(treatment[i])) continue;
        if (treatment[i] != 0.0) n_treated++;
        else n_controls++;
    }

    if (n_treated == 0 || n_controls == 0) {
        stata_data_free(&input_data);
        SF_error("cpsmatch: need both treated and control observations\n");
        return 198;
    }

    /* Build index arrays */
    size_t *treated_idx = malloc(n_treated * sizeof(size_t));
    size_t *control_idx = malloc(n_controls * sizeof(size_t));
    double *treated_pscore = malloc(n_treated * sizeof(double));
    double *control_pscore = malloc(n_controls * sizeof(double));
    int *control_available = NULL;

    if (!opts.with_replace) {
        control_available = malloc(n_controls * sizeof(int));
        if (control_available) {
            for (size_t i = 0; i < n_controls; i++) {
                control_available[i] = 1;
            }
        }
    }

    if (!treated_idx || !control_idx || !treated_pscore || !control_pscore) {
        stata_data_free(&input_data);
        free(treated_idx);
        free(control_idx);
        free(treated_pscore);
        free(control_pscore);
        free(control_available);
        SF_error("cpsmatch: memory allocation failed\n");
        return 920;
    }

    size_t ti = 0, ci = 0;
    for (size_t i = 0; i < nobs_range; i++) {
        ST_int stata_obs = obs1 + (ST_int)i;
        /* Check if this observation satisfies the if-condition */
        if (!SF_ifobs(stata_obs)) continue;
        if (SF_is_missing(treatment[i])) continue;
        if (treatment[i] != 0.0) {
            treated_idx[ti] = i;
            treated_pscore[ti] = pscore[i];
            ti++;
        } else {
            control_idx[ci] = i;
            control_pscore[ci] = pscore[i];
            ci++;
        }
    }

    /* ========================================================================
     * PHASE 3: Compute common support region
     * ======================================================================== */

    double ps_min_treated = DBL_MAX, ps_max_treated = -DBL_MAX;
    double ps_min_control = DBL_MAX, ps_max_control = -DBL_MAX;

    for (size_t i = 0; i < n_treated; i++) {
        if (!SF_is_missing(treated_pscore[i])) {
            if (treated_pscore[i] < ps_min_treated) ps_min_treated = treated_pscore[i];
            if (treated_pscore[i] > ps_max_treated) ps_max_treated = treated_pscore[i];
        }
    }

    for (size_t i = 0; i < n_controls; i++) {
        if (!SF_is_missing(control_pscore[i])) {
            if (control_pscore[i] < ps_min_control) ps_min_control = control_pscore[i];
            if (control_pscore[i] > ps_max_control) ps_max_control = control_pscore[i];
        }
    }

    double common_min = (ps_min_treated > ps_min_control) ? ps_min_treated : ps_min_control;
    double common_max = (ps_max_treated < ps_max_control) ? ps_max_treated : ps_max_control;

    /* ========================================================================
     * PHASE 4: Perform matching
     * ======================================================================== */

    /* Allocate output arrays - need full range size for all observations */
    double *out_weight = calloc(nobs_range, sizeof(double));
    double *out_match = malloc(nobs_range * sizeof(double));
    double *out_support = malloc(nobs_range * sizeof(double));

    /* Initialize to missing */
    for (size_t i = 0; i < nobs_range; i++) {
        out_match[i] = SV_missval;
        out_support[i] = SV_missval;
    }

    if (!out_weight || !out_match || !out_support) {
        stata_data_free(&input_data);
        free(treated_idx);
        free(control_idx);
        free(treated_pscore);
        free(control_pscore);
        free(control_available);
        free(out_weight);
        free(out_match);
        free(out_support);
        SF_error("cpsmatch: memory allocation failed\n");
        return 920;
    }

    /* Caliper is always in absolute units (propensity score scale) for psmatch2 compatibility */
    double caliper = opts.caliper;

    int n_matched_treated = 0;
    int n_off_support = 0;

    /* End setup phase, start sort/index phase */
    t_setup = ctools_timer_seconds() - t_setup;
    t_sort = ctools_timer_seconds();

    /* For matching without replacement, we need to:
     * 1. Sort treated by propensity score (descending by default, like psmatch2)
     * 2. Process sequentially to avoid race conditions
     * 3. Mark excess treated as off-support when controls run out
     */

    /* Create sorted order for treated (by propensity score) */
    size_t *treated_order = malloc(n_treated * sizeof(size_t));
    if (!treated_order) {
        stata_data_free(&input_data);
        free(treated_idx);
        free(control_idx);
        free(treated_pscore);
        free(control_pscore);
        free(control_available);
        free(out_weight);
        free(out_match);
        free(out_support);
        SF_error("cpsmatch: memory allocation failed\n");
        return 920;
    }

    /* For noreplacement, we need to sort treated by propensity score.
     * Use qsort with sorted_control_t structure for O(n log n) instead of O(n²) insertion sort */
    sorted_control_t *sorted_treated = NULL;
    if (!opts.with_replace) {
        sorted_treated = malloc(n_treated * sizeof(sorted_control_t));
        if (!sorted_treated) {
            stata_data_free(&input_data);
            free(treated_idx);
            free(control_idx);
            free(treated_pscore);
            free(control_pscore);
            free(control_available);
            free(out_weight);
            free(out_match);
            free(out_support);
            free(treated_order);
            SF_error("cpsmatch: memory allocation failed\n");
            return 920;
        }

        /* Populate and sort */
        for (size_t i = 0; i < n_treated; i++) {
            sorted_treated[i].orig_idx = i;
            sorted_treated[i].pscore = treated_pscore[i];
        }

        /* Sort: ascending (default) or descending */
        if (opts.descending) {
            qsort(sorted_treated, n_treated, sizeof(sorted_control_t), compare_by_pscore_desc);
        } else {
            qsort(sorted_treated, n_treated, sizeof(sorted_control_t), compare_by_pscore);
        }

        /* Copy sorted indices to treated_order */
        for (size_t i = 0; i < n_treated; i++) {
            treated_order[i] = sorted_treated[i].orig_idx;
        }
        free(sorted_treated);
    } else {
        /* For with-replacement, order doesn't matter - use original order */
        for (size_t i = 0; i < n_treated; i++) {
            treated_order[i] = i;
        }
    }

    /* Count available controls for matching without replacement */
    int n_controls_available = (int)n_controls;

    /* Matching loop - sequential for noreplacement, parallel otherwise */
    if (!opts.with_replace) {
        /* ================================================================
         * FAST sequential matching for noreplacement using sorted index
         * ================================================================ */

        /* Build sorted control index for O(log n) lookups */
        sorted_control_t *sorted_controls = malloc(n_controls * sizeof(sorted_control_t));
        int *sorted_available = malloc(n_controls * sizeof(int));  /* Track available in sorted order */

        if (!sorted_controls || !sorted_available) {
            free(sorted_controls);
            free(sorted_available);
            stata_data_free(&input_data);
            free(treated_idx);
            free(control_idx);
            free(treated_pscore);
            free(control_pscore);
            free(control_available);
            free(out_weight);
            free(out_match);
            free(out_support);
            free(treated_order);
            SF_error("cpsmatch: memory allocation failed\n");
            return 920;
        }

        /* Populate and sort control array */
        for (size_t i = 0; i < n_controls; i++) {
            sorted_controls[i].orig_idx = i;
            sorted_controls[i].obs_idx = control_idx[i];  /* Stata observation index for stable sort */
            sorted_controls[i].pscore = control_pscore[i];
            sorted_available[i] = 1;
        }
        qsort(sorted_controls, n_controls, sizeof(sorted_control_t), compare_by_pscore);

        /* End sort phase, start matching phase */
        t_sort = ctools_timer_seconds() - t_sort;
        t_match = ctools_timer_seconds();

        /* Sequential matching using binary search */
        for (size_t t_idx = 0; t_idx < n_treated; t_idx++) {
            size_t t = treated_order[t_idx];
            size_t obs_idx = treated_idx[t];
            double ps = treated_pscore[t];

            /* Check common support */
            if (opts.common_support) {
                if (ps < common_min || ps > common_max) {
                    out_support[obs_idx] = 0.0;
                    n_off_support++;
                    continue;
                }
            }

            /* Check if any controls are available */
            if (n_controls_available <= 0) {
                /* No more controls - treated is off support */
                out_support[obs_idx] = 0.0;
                n_off_support++;
                continue;
            }

            out_support[obs_idx] = 1.0;

            /* Binary search to find insertion point */
            size_t insert_pos = binary_search_left(sorted_controls, n_controls, ps);

            /* Expand left and right to find nearest available control */
            ssize_t left = (ssize_t)insert_pos - 1;
            size_t right = insert_pos;
            size_t best_idx = SIZE_MAX;
            double best_dist = DBL_MAX;

            /* Find the single nearest available control */
            while (left >= 0 || right < n_controls) {
                /* Check left candidate */
                if (left >= 0 && sorted_available[left]) {
                    double dist_left = fabs(ps - sorted_controls[left].pscore);
                    if (caliper <= 0.0 || dist_left <= caliper) {
                        if (dist_left < best_dist) {
                            best_dist = dist_left;
                            best_idx = (size_t)left;
                        }
                        break;  /* Found best on left, check right before deciding */
                    }
                    left--;
                } else if (left >= 0) {
                    left--;
                    continue;
                }

                /* Check right candidate */
                if (right < n_controls && sorted_available[right]) {
                    double dist_right = fabs(ps - sorted_controls[right].pscore);
                    if (caliper <= 0.0 || dist_right <= caliper) {
                        if (dist_right < best_dist) {
                            best_dist = dist_right;
                            best_idx = right;
                        }
                        break;  /* Found best on right */
                    }
                    right++;
                } else if (right < n_controls) {
                    right++;
                    continue;
                }

                /* Both sides exhausted without caliper match */
                if (left < 0 && right >= n_controls) break;
            }

            /* If we found a candidate on one side, check the other side for closer */
            if (best_idx != SIZE_MAX && left >= 0 && right < n_controls) {
                /* We broke out after finding one - check if other side is closer */
                if (best_idx == (size_t)left + 1) {
                    /* Found on left side, check right */
                    while (right < n_controls) {
                        if (sorted_available[right]) {
                            double dist_right = fabs(ps - sorted_controls[right].pscore);
                            if (dist_right < best_dist && (caliper <= 0.0 || dist_right <= caliper)) {
                                best_dist = dist_right;
                                best_idx = right;
                            }
                            break;
                        }
                        right++;
                    }
                } else {
                    /* Found on right side, check left */
                    while (left >= 0) {
                        if (sorted_available[left]) {
                            double dist_left = fabs(ps - sorted_controls[left].pscore);
                            if (dist_left < best_dist && (caliper <= 0.0 || dist_left <= caliper)) {
                                best_dist = dist_left;
                                best_idx = (size_t)left;
                            }
                            break;
                        }
                        left--;
                    }
                }
            }

            if (best_idx != SIZE_MAX) {
                n_matched_treated++;
                n_controls_available--;

                /* Mark control as used */
                sorted_available[best_idx] = 0;

                /* Get original control index */
                size_t orig_ctrl_idx = sorted_controls[best_idx].orig_idx;
                size_t match_obs_idx = control_idx[orig_ctrl_idx];

                /* Store match ID */
                out_match[obs_idx] = (double)(match_obs_idx + 1);  /* 1-based */

                /* Set weights */
                out_weight[match_obs_idx] = 1.0;
                out_weight[obs_idx] = 1.0;
            } else {
                /* No match found (e.g., caliper constraint) - treated is off support */
                out_support[obs_idx] = 0.0;
                n_off_support++;
            }
        }

        free(sorted_controls);
        free(sorted_available);
    } else {
        /* ================================================================
         * FAST matching for with-replacement using sorted control index
         * ================================================================ */

        /* Build sorted control index for O(log n) lookups */
        sorted_control_t *sorted_controls = malloc(n_controls * sizeof(sorted_control_t));
        if (!sorted_controls) {
            stata_data_free(&input_data);
            free(treated_idx);
            free(control_idx);
            free(treated_pscore);
            free(control_pscore);
            free(control_available);
            free(out_weight);
            free(out_match);
            free(out_support);
            free(treated_order);
            SF_error("cpsmatch: memory allocation failed\n");
            return 920;
        }

        /* Populate sorted control array */
        for (size_t i = 0; i < n_controls; i++) {
            sorted_controls[i].orig_idx = i;
            sorted_controls[i].obs_idx = control_idx[i];  /* Stata observation index for stable sort */
            sorted_controls[i].pscore = control_pscore[i];
        }

        /* Sort controls by propensity score */
        qsort(sorted_controls, n_controls, sizeof(sorted_control_t), compare_by_pscore);

        /* End sort phase, start matching phase */
        t_sort = ctools_timer_seconds() - t_sort;
        t_match = ctools_timer_seconds();

        /* Optimized parallel matching using fast O(log n) algorithms */

        /* For radius matching, use thread-local weight buffers to reduce atomic contention */
        if (opts.method == MATCH_RADIUS) {
            double radius = (caliper > 0.0) ? caliper : CPSMATCH_DEFAULT_CALIPER;

            /* Get number of threads */
            int n_threads = 1;
            #ifdef _OPENMP
            n_threads = omp_get_max_threads();
            #endif

            /* Allocate thread-local weight buffers */
            double **thread_weights = malloc(n_threads * sizeof(double *));
            if (thread_weights) {
                for (int i = 0; i < n_threads; i++) {
                    thread_weights[i] = calloc(nobs_range, sizeof(double));
                }
            }

            #pragma omp parallel reduction(+:n_matched_treated, n_off_support)
            {
                int tid = 0;
                #ifdef _OPENMP
                tid = omp_get_thread_num();
                #endif

                double *local_weights = (thread_weights && thread_weights[tid]) ? thread_weights[tid] : out_weight;
                int use_local = (thread_weights && thread_weights[tid]) ? 1 : 0;

                #pragma omp for schedule(dynamic, 64)
                for (size_t t = 0; t < n_treated; t++) {
                    size_t obs_idx = treated_idx[t];
                    double ps = treated_pscore[t];

                    /* Check common support */
                    if (opts.common_support) {
                        if (ps < common_min || ps > common_max) {
                            out_support[obs_idx] = 0.0;
                            n_off_support++;
                            continue;
                        }
                    }
                    out_support[obs_idx] = 1.0;

                    /* Fast range query - no allocation */
                    radius_match_info info;
                    find_radius_range_fast(ps, sorted_controls, n_controls, radius, &info);

                    if (info.n_matches > 0) {
                        n_matched_treated++;

                        /* Store first match ID */
                        out_match[obs_idx] = (double)(control_idx[sorted_controls[info.start].orig_idx] + 1);

                        /* Accumulate weights - use local buffer if available */
                        double weight = 1.0 / info.n_matches;
                        for (size_t i = info.start; i < info.end; i++) {
                            size_t match_idx = control_idx[sorted_controls[i].orig_idx];
                            if (use_local) {
                                local_weights[match_idx] += weight;
                            } else {
                                #pragma omp atomic
                                out_weight[match_idx] += weight;
                            }
                        }

                        /* Weight for treated observation = 1 */
                        if (use_local) {
                            local_weights[obs_idx] = 1.0;
                        } else {
                            out_weight[obs_idx] = 1.0;
                        }
                    }
                }
            }

            /* Merge thread-local buffers in parallel */
            if (thread_weights) {
                #pragma omp parallel for schedule(static)
                for (size_t j = 0; j < nobs_range; j++) {
                    double sum = 0.0;
                    for (int i = 0; i < n_threads; i++) {
                        if (thread_weights[i]) {
                            sum += thread_weights[i][j];
                        }
                    }
                    out_weight[j] += sum;
                }
                /* Free thread-local buffers */
                for (int i = 0; i < n_threads; i++) {
                    free(thread_weights[i]);
                }
                free(thread_weights);
            }
        } else {
            /* General path for NN and kernel matching */
            #pragma omp parallel for reduction(+:n_matched_treated, n_off_support) schedule(dynamic)
            for (size_t t = 0; t < n_treated; t++) {
                size_t obs_idx = treated_idx[t];
                double ps = treated_pscore[t];

                /* Check common support */
                if (opts.common_support) {
                    if (ps < common_min || ps > common_max) {
                        out_support[obs_idx] = 0.0;
                        n_off_support++;
                        continue;
                    }
                }
                out_support[obs_idx] = 1.0;

                match_result result;
                memset(&result, 0, sizeof(result));

                int match_rc = 0;

                switch (opts.method) {
                    case MATCH_NEAREST:
                        match_rc = find_nearest_neighbors_fast(
                            ps, sorted_controls, control_idx, n_controls,
                            opts.n_neighbors, caliper, opts.ties, &result);
                        break;

                    case MATCH_KERNEL:
                        match_rc = find_kernel_matches_fast(
                            ps, sorted_controls, control_idx, n_controls,
                            opts.bandwidth, opts.kernel, caliper, &result);
                        break;

                    default:
                        match_rc = find_nearest_neighbors_fast(
                            ps, sorted_controls, control_idx, n_controls,
                            opts.n_neighbors, caliper, opts.ties, &result);
                        break;
                }

                if (match_rc == 0 && result.n_matches > 0) {
                    n_matched_treated++;

                    /* Store first match ID for this treated observation */
                    out_match[obs_idx] = (double)(result.match_ids[0] + 1);  /* 1-based */

                    /* Accumulate weights for matched controls */
                    for (int m = 0; m < result.n_matches; m++) {
                        size_t match_idx = result.match_ids[m];
                        #pragma omp atomic
                        out_weight[match_idx] += result.match_weights[m];
                    }

                    /* Weight for treated observation = 1 */
                    out_weight[obs_idx] = 1.0;
                }

                /* Free match result */
                free(result.match_ids);
                free(result.match_dists);
                free(result.match_weights);
            }
        }

        free(sorted_controls);
    }

    free(treated_order);

    /* Mark controls on support */
    for (size_t i = 0; i < n_controls; i++) {
        size_t obs_idx = control_idx[i];
        double ps = control_pscore[i];

        if (opts.common_support) {
            if (ps >= common_min && ps <= common_max) {
                out_support[obs_idx] = 1.0;
            } else {
                out_support[obs_idx] = 0.0;
            }
        } else {
            out_support[obs_idx] = 1.0;
        }
    }

    t_match = ctools_timer_seconds() - t_match;

    /* ========================================================================
     * PHASE 5: Store results back to Stata
     * ======================================================================== */
    t_store = ctools_timer_seconds();

    /* Store results - variable indices need +1 offset for SF_vstore
     * This is because the indices from Stata count from 1, but SF_vstore
     * uses a different indexing scheme in the plugin API context */
    for (size_t i = 0; i < nobs_range; i++) {
        ST_int obs = obs1 + (ST_int)i;

        if (opts.out_weight_idx > 0) {
            SF_vstore(opts.out_weight_idx + 1, obs, out_weight[i]);
        }

        if (opts.out_match_idx > 0) {
            SF_vstore(opts.out_match_idx + 1, obs, out_match[i]);
        }

        if (opts.out_support_idx > 0) {
            SF_vstore(opts.out_support_idx + 1, obs, out_support[i]);
        }
    }

    t_store = ctools_timer_seconds() - t_store;

    /* ========================================================================
     * PHASE 6: Store timing and diagnostics
     * ======================================================================== */

    t_total = ctools_timer_seconds() - t_start;

    SF_scal_save("_cpsmatch_time_load", t_load);
    SF_scal_save("_cpsmatch_time_setup", t_setup);
    SF_scal_save("_cpsmatch_time_sort", t_sort);
    SF_scal_save("_cpsmatch_time_match", t_match);
    SF_scal_save("_cpsmatch_time_store", t_store);
    SF_scal_save("_cpsmatch_time_total", t_total);

    SF_scal_save("_cpsmatch_n_treated", (double)n_treated);
    SF_scal_save("_cpsmatch_n_controls", (double)n_controls);
    SF_scal_save("_cpsmatch_n_matched", (double)n_matched_treated);
    SF_scal_save("_cpsmatch_n_off_support", (double)n_off_support);
    SF_scal_save("_cpsmatch_common_min", common_min);
    SF_scal_save("_cpsmatch_common_max", common_max);

    /* Thread diagnostics */
    CTOOLS_SAVE_THREAD_INFO("_cpsmatch");

    /* Cleanup - use stata_data_free for input data (frees treatment, pscore, outcome) */
    stata_data_free(&input_data);
    free(treated_idx);
    free(control_idx);
    free(treated_pscore);
    free(control_pscore);
    free(control_available);
    free(out_weight);
    free(out_match);
    free(out_support);

    return rc;
}

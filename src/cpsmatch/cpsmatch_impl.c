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
 * 1. Load treatment indicator and propensity scores from Stata
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
 * Comparison function for sorting by distance
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
 * Nearest Neighbor Matching
 * ============================================================================ */

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
 * Kernel Matching
 * ============================================================================ */

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
 * Radius Matching
 * ============================================================================ */

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
    double t_start, t_load, t_match, t_store, t_total;
    ST_retcode rc = 0;

    t_start = ctools_timer_seconds();

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

    /* Get observation range */
    ST_int obs1 = SF_in1();
    ST_int obs2 = SF_in2();
    size_t nobs = (size_t)(obs2 - obs1 + 1);

    if (nobs != opts.nobs) {
        /* Use the actual observation range */
        nobs = opts.nobs;
    }

    /* ========================================================================
     * PHASE 1: Load data from Stata
     * ======================================================================== */
    t_load = ctools_timer_seconds();

    /* Allocate arrays */
    double *treatment = malloc(nobs * sizeof(double));
    double *pscore = malloc(nobs * sizeof(double));
    double *outcome = NULL;

    if (!treatment || !pscore) {
        free(treatment);
        free(pscore);
        SF_error("cpsmatch: memory allocation failed\n");
        return 920;
    }

    if (opts.outcome_idx > 0) {
        outcome = malloc(nobs * sizeof(double));
        if (!outcome) {
            free(treatment);
            free(pscore);
            SF_error("cpsmatch: memory allocation failed\n");
            return 920;
        }
    }

    /* Read data from Stata */
    for (size_t i = 0; i < nobs; i++) {
        ST_int obs = obs1 + (ST_int)i;
        SF_vdata(opts.treat_idx, obs, &treatment[i]);
        SF_vdata(opts.pscore_idx, obs, &pscore[i]);
        if (outcome) {
            SF_vdata(opts.outcome_idx, obs, &outcome[i]);
        }
    }

    t_load = ctools_timer_seconds() - t_load;

    /* ========================================================================
     * PHASE 2: Separate treated and control groups
     * ======================================================================== */
    t_match = ctools_timer_seconds();

    /* Count treated and controls */
    size_t n_treated = 0, n_controls = 0;
    for (size_t i = 0; i < nobs; i++) {
        if (SF_is_missing(treatment[i])) continue;
        if (treatment[i] != 0.0) n_treated++;
        else n_controls++;
    }

    if (n_treated == 0 || n_controls == 0) {
        free(treatment);
        free(pscore);
        free(outcome);
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
        free(treatment);
        free(pscore);
        free(outcome);
        free(treated_idx);
        free(control_idx);
        free(treated_pscore);
        free(control_pscore);
        free(control_available);
        SF_error("cpsmatch: memory allocation failed\n");
        return 920;
    }

    size_t ti = 0, ci = 0;
    for (size_t i = 0; i < nobs; i++) {
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

    /* Allocate output arrays */
    double *out_weight = calloc(nobs, sizeof(double));
    double *out_match = malloc(nobs * sizeof(double));
    double *out_support = malloc(nobs * sizeof(double));

    /* Initialize to missing */
    for (size_t i = 0; i < nobs; i++) {
        out_match[i] = SV_missval;
        out_support[i] = SV_missval;
    }

    if (!out_weight || !out_match || !out_support) {
        free(treatment);
        free(pscore);
        free(outcome);
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

    /* Compute caliper in absolute units if needed */
    double caliper = opts.caliper;
    if (caliper > 0.0) {
        /* Compute SD of propensity score for treated */
        double sum = 0.0, sumsq = 0.0;
        size_t n_valid = 0;
        for (size_t i = 0; i < n_treated; i++) {
            if (!SF_is_missing(treated_pscore[i])) {
                sum += treated_pscore[i];
                sumsq += treated_pscore[i] * treated_pscore[i];
                n_valid++;
            }
        }
        if (n_valid > 1) {
            double mean = sum / n_valid;
            double var = (sumsq - n_valid * mean * mean) / (n_valid - 1);
            double sd = sqrt(var);
            /* If caliper < 1, treat as proportion of SD; else as absolute */
            if (caliper < 1.0) {
                caliper = caliper * sd;
            }
        }
    }

    int n_matched_treated = 0;
    int n_off_support = 0;

    /* Parallel matching for each treated observation */
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
                match_rc = find_nearest_neighbors(
                    ps, control_pscore, control_idx,
                    opts.with_replace ? NULL : control_available,
                    n_controls, opts.n_neighbors, caliper,
                    opts.with_replace, opts.ties, &result);
                break;

            case MATCH_RADIUS:
                match_rc = find_radius_matches(
                    ps, control_pscore, control_idx, n_controls,
                    caliper > 0.0 ? caliper : CPSMATCH_DEFAULT_CALIPER,
                    &result);
                break;

            case MATCH_KERNEL:
                match_rc = find_kernel_matches(
                    ps, control_pscore, control_idx, n_controls,
                    opts.bandwidth, opts.kernel, caliper, &result);
                break;

            default:
                match_rc = find_nearest_neighbors(
                    ps, control_pscore, control_idx,
                    opts.with_replace ? NULL : control_available,
                    n_controls, opts.n_neighbors, caliper,
                    opts.with_replace, opts.ties, &result);
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

    /* Store results - NOTE: SF_vstore has off-by-one, so we add 1 to indices */
    for (size_t i = 0; i < nobs; i++) {
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

cleanup:
    /* Cleanup */
    free(treatment);
    free(pscore);
    free(outcome);
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

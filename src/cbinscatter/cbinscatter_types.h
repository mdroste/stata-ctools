/*
 * cbinscatter_types.h
 *
 * Type definitions for cbinscatter module
 * Part of the ctools Stata plugin suite
 */

#ifndef CBINSCATTER_TYPES_H
#define CBINSCATTER_TYPES_H

#include "stplugin.h"
#include "ctools_types.h"

/* ========================================================================
 * Configuration structure - populated from Stata scalars
 * ======================================================================== */

typedef struct {
    ST_int nquantiles;       /* Number of bins (default 20) */
    ST_int linetype;         /* 0=none, 1=linear, 2=quadratic, 3=cubic */
    ST_int compute_se;       /* Whether to compute standard errors */
    ST_int discrete;         /* Treat x as discrete (unique values = bins) */
    ST_int has_controls;     /* Has control variables */
    ST_int num_controls;     /* Number of control variables */
    ST_int has_absorb;       /* Has absorb() option */
    ST_int num_absorb;       /* Number of absorb variables */
    ST_int has_by;           /* Has by() option */
    ST_int has_weights;      /* Has weights */
    ST_int weight_type;      /* 0=none, 1=aweight, 2=fweight, 3=pweight, 4=iweight */
    ST_int verbose;          /* Verbose output */
    ST_int reportreg;        /* Report regression underlying binscatter */
    /*
     * Residualization method (following Cattaneo et al. "On Binscatter")
     * 0 = classic: residualize both Y and X, bin on residualized X
     * 1 = binsreg: bin on raw X, residualize Y only (conditional means)
     */
    ST_int method;
    /* HDFE parameters (if absorb specified) */
    ST_int maxiter;          /* CG solver max iterations */
    ST_double tolerance;     /* CG solver tolerance */
} BinscatterConfig;

/* ========================================================================
 * Bin statistics for a single bin
 * ======================================================================== */

typedef struct {
    ST_int bin_id;           /* Bin number (1 to nquantiles) */
    ST_double x_mean;        /* Mean of x within bin */
    ST_double y_mean;        /* Mean of y within bin */
    ST_double x_se;          /* SE of x mean (optional) */
    ST_double y_se;          /* SE of y mean (optional) */
    ST_int n_obs;            /* Observations in bin */
    ST_double sum_weights;   /* Sum of weights in bin (if weighted) */
} BinStats;

/* ========================================================================
 * Complete results for one by-group
 * ======================================================================== */

typedef struct {
    ST_int by_group_id;      /* Group identifier (1-based) */
    ST_int num_bins;         /* Number of bins with data */
    BinStats *bins;          /* Array of bin statistics */
    /* Fit coefficients (if line requested) */
    ST_int fit_order;        /* Polynomial order (1=linear, 2=quad, etc.) */
    ST_double *fit_coefs;    /* Coefficients [fit_order+1] */
    ST_double fit_r2;        /* R-squared of fit */
    ST_int fit_n;            /* Observations used in fit */
} ByGroupResult;

/* ========================================================================
 * Master results structure
 * ======================================================================== */

typedef struct {
    ST_int num_by_groups;    /* Number of by-groups (1 if no by()) */
    ByGroupResult *groups;   /* Array of by-group results */
    ST_int nquantiles;       /* Number of bins requested */
    ST_int total_obs;        /* Total observations used */
    ST_int obs_dropped;      /* Observations dropped (missing, etc.) */
} BinscatterResults;

/* ========================================================================
 * Working data for bin computation
 * ======================================================================== */

typedef struct {
    ST_int N;                /* Number of valid observations */
    ST_double *y;            /* y values (residualized if controls/absorb) */
    ST_double *x;            /* x values (residualized if controls/absorb) */
    ST_int y_owned;          /* Whether y is owned (should be freed) */
    ST_int x_owned;          /* Whether x is owned (should be freed) */
} BinscatterWorkData;

/* ========================================================================
 * Return codes
 * ======================================================================== */

#define CBINSCATTER_OK              0
#define CBINSCATTER_ERR_MEMORY      920   /* Memory allocation failed */
#define CBINSCATTER_ERR_NOOBS       2001  /* No valid observations */
#define CBINSCATTER_ERR_FEW_OBS     2002  /* Fewer obs than bins */
#define CBINSCATTER_ERR_SYNTAX      198   /* Syntax error */
#define CBINSCATTER_ERR_SINGULAR    199   /* Singular matrix */
#define CBINSCATTER_ERR_HDFE        459   /* HDFE solver failed */
#define CBINSCATTER_ERR_INVALID     2003  /* Invalid parameter value */

/* ========================================================================
 * Helper functions
 * ======================================================================== */

/* Initialize config to defaults */
void cbinscatter_config_init(BinscatterConfig *config);

/* Initialize results structure */
void cbinscatter_results_init(BinscatterResults *results);

/* Free results structure */
void cbinscatter_results_free(BinscatterResults *results);

/* Initialize work data */
void cbinscatter_workdata_init(BinscatterWorkData *work);

/* Free work data */
void cbinscatter_workdata_free(BinscatterWorkData *work);

#endif /* CBINSCATTER_TYPES_H */

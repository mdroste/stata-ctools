/*
 * cqreg_types.h
 *
 * Core data structures for C-accelerated quantile regression.
 * Part of the ctools suite.
 */

#ifndef CQREG_TYPES_H
#define CQREG_TYPES_H

#include "stplugin.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Default configuration */
#define CQREG_DEFAULT_MAXITER     50
#define CQREG_DEFAULT_TOL         1e-8
#define CQREG_DEFAULT_MU_INIT     0.1
#define CQREG_DEFAULT_SIGMA       0.1
#define CQREG_NUM_THREADS         8
#define CQREG_CACHE_LINE          64

/* Compiler hints */
#if defined(__GNUC__) || defined(__clang__)
#define CQREG_RESTRICT   __restrict__
#define CQREG_LIKELY(x)   __builtin_expect(!!(x), 1)
#define CQREG_UNLIKELY(x) __builtin_expect(!!(x), 0)
#define CQREG_ALIGNED(n)  __attribute__((aligned(n)))
#else
#define CQREG_RESTRICT
#define CQREG_LIKELY(x)   (x)
#define CQREG_UNLIKELY(x) (x)
#define CQREG_ALIGNED(n)
#endif

/* VCE types */
typedef enum {
    CQREG_VCE_IID     = 0,   /* IID with sparsity estimation */
    CQREG_VCE_ROBUST  = 1,   /* Robust (sandwich) */
    CQREG_VCE_CLUSTER = 2    /* Clustered */
} cqreg_vce_type;

/* Bandwidth selection methods for sparsity estimation */
typedef enum {
    CQREG_BW_HSHEATHER   = 0,   /* Hall-Sheather (default) */
    CQREG_BW_BOFINGER    = 1,   /* Bofinger */
    CQREG_BW_CHAMBERLAIN = 2    /* Chamberlain */
} cqreg_bw_method;

/* Density estimation methods for VCE */
typedef enum {
    CQREG_DENSITY_RESIDUAL = 0, /* Difference quotient on residuals (Stata: vce(iid, residual)) */
    CQREG_DENSITY_FITTED   = 1, /* Kernel density on fitted values (Stata default) */
    CQREG_DENSITY_KERNEL   = 2  /* Kernel density on residuals */
} cqreg_density_method;

/* IPM solver configuration */
typedef struct {
    ST_int maxiter;           /* Maximum iterations (default: 50) */
    ST_double tol_primal;     /* Primal feasibility tolerance */
    ST_double tol_dual;       /* Dual feasibility tolerance */
    ST_double tol_gap;        /* Duality gap tolerance */
    ST_int verbose;           /* Verbosity level (0=none, 1=summary, 2=detailed) */
    ST_int use_mehrotra;      /* Use predictor-corrector (default: 1) */
    ST_double mu_init;        /* Initial barrier parameter */
    ST_double sigma;          /* Centering parameter (0.1-0.3) */
} cqreg_ipm_config;

/* IPM solver state - primal-dual interior point method */
typedef struct {
    /* Problem dimensions */
    ST_int N;                 /* Number of observations */
    ST_int K;                 /* Number of regressors (including constant) */

    /* Coefficient vector */
    ST_double *beta;          /* Regression coefficients (K) */

    /* Primal variables (size N each) */
    ST_double *u;             /* Positive deviations: y - X*beta = u - v */
    ST_double *v;             /* Negative deviations */

    /* Dual variables (size N each) */
    ST_double *lambda_u;      /* Dual for u >= 0 */
    ST_double *lambda_v;      /* Dual for v >= 0 */

    /* Diagonal scaling for normal equations */
    ST_double *D;             /* D[i] = 1/(lambda_u[i]/u[i] + lambda_v[i]/v[i]) */

    /* Normal equations storage */
    ST_double *XDX;           /* X' * D * X (K x K) */
    ST_double *XDXcopy;       /* Copy for Cholesky (K x K) */
    ST_double *L;             /* Cholesky factor of XDX (K x K, lower triangular) */
    ST_double *rhs;           /* Right-hand side for normal equations (K) */

    /* Search directions */
    ST_double *delta_beta;    /* Search direction for beta (K) */
    ST_double *delta_u;       /* Search direction for u (N) */
    ST_double *delta_v;       /* Search direction for v (N) */
    ST_double *delta_lambda_u;/* Search direction for lambda_u (N) */
    ST_double *delta_lambda_v;/* Search direction for lambda_v (N) */

    /* Mehrotra predictor-corrector affine directions */
    ST_double *delta_u_aff;       /* Affine direction for u (N) */
    ST_double *delta_v_aff;       /* Affine direction for v (N) */
    ST_double *delta_lambda_u_aff;/* Affine direction for lambda_u (N) */
    ST_double *delta_lambda_v_aff;/* Affine direction for lambda_v (N) */

    /* Residuals */
    ST_double *r_primal;      /* Primal residual: y - X*beta - u + v (N) */
    ST_double *r_dual;        /* Dual residual: X'*(lambda_u - lambda_v) - c (K) */

    /* Working arrays */
    ST_double *work_N;        /* General workspace (N) */
    ST_double *work_K;        /* General workspace (K) */

    /* Per-thread buffers for parallel operations */
    ST_double **thread_buf;   /* num_threads x N scratch buffers */
    ST_int num_threads;

    /* Frisch-Newton solver workspace (pre-allocated to avoid malloc in hot loop) */
    ST_double *fn_xp;         /* Primal x (N) */
    ST_double *fn_s;          /* Primal s (N) */
    ST_double *fn_yd;         /* Dual y = -beta (K) */
    ST_double *fn_z;          /* Dual slack z (N) */
    ST_double *fn_w;          /* Dual slack w (N) */
    ST_double *fn_dx;         /* Direction for x (N) */
    ST_double *fn_ds;         /* Direction for s (N) */
    ST_double *fn_dy;         /* Direction for y (K) */
    ST_double *fn_dz;         /* Direction for z (N) */
    ST_double *fn_dw;         /* Direction for w (N) */
    ST_double *fn_fx;         /* Bound for x (N) */
    ST_double *fn_fs;         /* Bound for s (N) */
    ST_double *fn_fz;         /* Bound for z (N) */
    ST_double *fn_fw;         /* Bound for w (N) */
    ST_double *fn_q;          /* Diagonal scaling (N) */
    ST_double *fn_r;          /* Residual (N) */
    ST_double *fn_tmp;        /* Temp array (N) */
    ST_double *fn_sinv;       /* 1/s for corrector (N) */
    ST_double *fn_Xq;         /* Scaled X, Q*X (N*K) */

    /* Barrier and convergence tracking */
    ST_double mu;             /* Current barrier parameter */
    ST_double primal_obj;     /* Primal objective value */
    ST_double dual_obj;       /* Dual objective value */
    ST_double gap;            /* Duality gap */
    ST_double primal_infeas;  /* Primal infeasibility */
    ST_double dual_infeas;    /* Dual infeasibility */

    /* Configuration */
    cqreg_ipm_config config;

    /* Status */
    ST_int converged;         /* 1 if converged, 0 otherwise */
    ST_int iterations;        /* Number of iterations used */
} cqreg_ipm_state;

/* Sparsity estimation state */
typedef struct {
    ST_double quantile;       /* Quantile being estimated */
    ST_double bandwidth;      /* Kernel bandwidth */
    ST_double sparsity;       /* Estimated sparsity (1/f(F^{-1}(q))) */
    cqreg_bw_method bw_method;/* Bandwidth selection method */
    ST_int N;                 /* Number of observations */
    ST_double *sorted_resid;  /* Sorted residuals for density estimation (N) */
} cqreg_sparsity_state;

/* Main cqreg state - orchestrates the full regression */
typedef struct {
    /* Dimensions */
    ST_int N;                 /* Observations (after singleton removal if HDFE) */
    ST_int N_original;        /* Original number of observations */
    ST_int K;                 /* Regressors (excluding absorbed FEs) */
    ST_int G;                 /* Number of FE groups (0 if no absorb) */

    /* Quantile parameter */
    ST_double quantile;       /* Target quantile (0 < q < 1) */

    /* Data (column-major layout: data[k*N + i] = obs i of var k) */
    ST_double *y;             /* Dependent variable (N) */
    ST_double *X;             /* Regressors (N x K), constant is last column */
    ST_double *weights;       /* Observation weights (N), NULL if unweighted */

    /* Observation mask for HDFE singleton removal */
    ST_int *obs_mask;         /* 1 = keep, 0 = dropped (N_original) */

    /* HDFE integration */
    ST_int has_hdfe;          /* Whether absorb() was specified */
    void *hdfe_state;         /* Pointer to HDFE_State from creghdfe */
    ST_int df_a;              /* Degrees of freedom absorbed by FEs */

    /* Results */
    ST_double *beta;          /* Coefficient estimates (K) */
    ST_double *V;             /* Variance-covariance matrix (K x K) */
    ST_double *residuals;     /* Residuals y - X*beta (N) */
    ST_double *fitted;        /* Fitted values X*beta (N) */
    ST_double *obs_density;   /* Per-observation density estimates f_i(0) (N) */
    ST_double sum_adev;       /* Sum of absolute deviations (objective value) */
    ST_double sum_rdev;       /* Sum of raw deviations */
    ST_double sparsity;       /* Estimated sparsity */
    ST_double bandwidth;      /* Bandwidth used for sparsity estimation */

    /* VCE configuration */
    cqreg_vce_type vce_type;
    cqreg_bw_method bw_method;
    cqreg_density_method density_method;  /* Fitted vs residual density estimation */
    ST_int *cluster_ids;      /* Cluster identifiers (N), NULL if not clustered */
    ST_int num_clusters;      /* Number of unique clusters */

    /* IPM solver state */
    cqreg_ipm_state *ipm;

    /* Sparsity estimation state */
    cqreg_sparsity_state *sparsity_state;

    /* Convergence info */
    ST_int iterations;        /* IPM iterations to convergence */
    ST_int converged;         /* 1 if converged, 0 otherwise */

    /* Timing (if verbose) */
    ST_double time_load;      /* Data loading time */
    ST_double time_hdfe;      /* HDFE projection time */
    ST_double time_ipm;       /* IPM solver time */
    ST_double time_vce;       /* VCE computation time */
    ST_double time_total;     /* Total time */

    /* Thread count */
    ST_int num_threads;
} cqreg_state;

/* ============================================================================
 * Memory management function declarations
 * ============================================================================ */

/* Initialize IPM config with defaults */
void cqreg_ipm_config_init(cqreg_ipm_config *config);

/* Create and initialize IPM solver state */
cqreg_ipm_state *cqreg_ipm_create(ST_int N, ST_int K, const cqreg_ipm_config *config);

/* Free IPM solver state */
void cqreg_ipm_free(cqreg_ipm_state *ipm);

/* Create sparsity estimation state */
cqreg_sparsity_state *cqreg_sparsity_create(ST_int N, ST_double quantile, cqreg_bw_method method);

/* Free sparsity estimation state */
void cqreg_sparsity_free(cqreg_sparsity_state *sp);

/* Create main cqreg state */
cqreg_state *cqreg_state_create(ST_int N, ST_int K);

/* Free main cqreg state */
void cqreg_state_free(cqreg_state *state);

/* Aligned memory allocation */
void *cqreg_aligned_alloc(size_t size, size_t alignment);

/* Aligned memory free */
void cqreg_aligned_free(void *ptr);

#endif /* CQREG_TYPES_H */

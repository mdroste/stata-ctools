/*
 * cqreg_types.c
 *
 * Memory management for cqreg data structures.
 * Part of the ctools suite.
 */

#include "cqreg_types.h"
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Aligned memory allocation
 * ============================================================================ */

void *cqreg_aligned_alloc(size_t size, size_t alignment)
{
    void *ptr = NULL;
#if defined(_WIN32)
    ptr = _aligned_malloc(size, alignment);
#elif defined(__APPLE__) || defined(__linux__)
    if (posix_memalign(&ptr, alignment, size) != 0) {
        ptr = NULL;
    }
#else
    ptr = malloc(size);
#endif
    return ptr;
}

void cqreg_aligned_free(void *ptr)
{
    if (ptr == NULL) return;
#if defined(_WIN32)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

/* ============================================================================
 * IPM Configuration
 * ============================================================================ */

void cqreg_ipm_config_init(cqreg_ipm_config *config)
{
    if (config == NULL) return;

    config->maxiter = CQREG_DEFAULT_MAXITER;
    config->tol_primal = CQREG_DEFAULT_TOL;
    config->tol_dual = CQREG_DEFAULT_TOL;
    config->tol_gap = CQREG_DEFAULT_TOL;
    config->verbose = 0;
    config->use_mehrotra = 1;
    config->mu_init = CQREG_DEFAULT_MU_INIT;
    config->sigma = CQREG_DEFAULT_SIGMA;
}

/* ============================================================================
 * IPM State Management
 * ============================================================================ */

cqreg_ipm_state *cqreg_ipm_create(ST_int N, ST_int K, const cqreg_ipm_config *config)
{
    cqreg_ipm_state *ipm = NULL;
    ST_int t;

    if (N <= 0 || K <= 0) return NULL;

    ipm = (cqreg_ipm_state *)calloc(1, sizeof(cqreg_ipm_state));
    if (ipm == NULL) return NULL;

    ipm->N = N;
    ipm->K = K;

    /* Copy configuration */
    if (config != NULL) {
        memcpy(&ipm->config, config, sizeof(cqreg_ipm_config));
    } else {
        cqreg_ipm_config_init(&ipm->config);
    }

    /* Determine thread count */
#ifdef _OPENMP
    ipm->num_threads = omp_get_max_threads();
    if (ipm->num_threads > CQREG_NUM_THREADS) {
        ipm->num_threads = CQREG_NUM_THREADS;
    }
#else
    ipm->num_threads = 1;
#endif

    /* Allocate coefficient vector */
    ipm->beta = (ST_double *)cqreg_aligned_alloc(K * sizeof(ST_double), CQREG_CACHE_LINE);
    if (ipm->beta == NULL) goto cleanup;

    /* Allocate primal variables */
    ipm->u = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->v = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    if (ipm->u == NULL || ipm->v == NULL) goto cleanup;

    /* Allocate dual variables */
    ipm->lambda_u = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->lambda_v = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    if (ipm->lambda_u == NULL || ipm->lambda_v == NULL) goto cleanup;

    /* Allocate diagonal scaling */
    ipm->D = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    if (ipm->D == NULL) goto cleanup;

    /* Allocate normal equations storage */
    ipm->XDX = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->XDXcopy = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->L = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->rhs = (ST_double *)cqreg_aligned_alloc(K * sizeof(ST_double), CQREG_CACHE_LINE);
    if (ipm->XDX == NULL || ipm->XDXcopy == NULL || ipm->L == NULL || ipm->rhs == NULL) goto cleanup;

    /* Allocate search directions */
    ipm->delta_beta = (ST_double *)cqreg_aligned_alloc(K * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->delta_u = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->delta_v = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->delta_lambda_u = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->delta_lambda_v = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    if (ipm->delta_beta == NULL || ipm->delta_u == NULL || ipm->delta_v == NULL ||
        ipm->delta_lambda_u == NULL || ipm->delta_lambda_v == NULL) goto cleanup;

    /* Allocate Mehrotra affine directions */
    ipm->delta_u_aff = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->delta_v_aff = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->delta_lambda_u_aff = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->delta_lambda_v_aff = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    if (ipm->delta_u_aff == NULL || ipm->delta_v_aff == NULL ||
        ipm->delta_lambda_u_aff == NULL || ipm->delta_lambda_v_aff == NULL) goto cleanup;

    /* Allocate residuals */
    ipm->r_primal = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->r_dual = (ST_double *)cqreg_aligned_alloc(K * sizeof(ST_double), CQREG_CACHE_LINE);
    if (ipm->r_primal == NULL || ipm->r_dual == NULL) goto cleanup;

    /* Allocate working arrays */
    ipm->work_N = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    ipm->work_K = (ST_double *)cqreg_aligned_alloc(K * sizeof(ST_double), CQREG_CACHE_LINE);
    if (ipm->work_N == NULL || ipm->work_K == NULL) goto cleanup;

    /* Allocate per-thread buffers */
    ipm->thread_buf = (ST_double **)calloc(ipm->num_threads, sizeof(ST_double *));
    if (ipm->thread_buf == NULL) goto cleanup;

    for (t = 0; t < ipm->num_threads; t++) {
        ipm->thread_buf[t] = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
        if (ipm->thread_buf[t] == NULL) goto cleanup;
    }

    /* Initialize convergence tracking */
    ipm->mu = ipm->config.mu_init;
    ipm->converged = 0;
    ipm->iterations = 0;

    return ipm;

cleanup:
    cqreg_ipm_free(ipm);
    return NULL;
}

void cqreg_ipm_free(cqreg_ipm_state *ipm)
{
    ST_int t;

    if (ipm == NULL) return;

    /* Free coefficient vector */
    cqreg_aligned_free(ipm->beta);

    /* Free primal variables */
    cqreg_aligned_free(ipm->u);
    cqreg_aligned_free(ipm->v);

    /* Free dual variables */
    cqreg_aligned_free(ipm->lambda_u);
    cqreg_aligned_free(ipm->lambda_v);

    /* Free diagonal scaling */
    cqreg_aligned_free(ipm->D);

    /* Free normal equations storage */
    cqreg_aligned_free(ipm->XDX);
    cqreg_aligned_free(ipm->XDXcopy);
    cqreg_aligned_free(ipm->L);
    cqreg_aligned_free(ipm->rhs);

    /* Free search directions */
    cqreg_aligned_free(ipm->delta_beta);
    cqreg_aligned_free(ipm->delta_u);
    cqreg_aligned_free(ipm->delta_v);
    cqreg_aligned_free(ipm->delta_lambda_u);
    cqreg_aligned_free(ipm->delta_lambda_v);

    /* Free Mehrotra affine directions */
    cqreg_aligned_free(ipm->delta_u_aff);
    cqreg_aligned_free(ipm->delta_v_aff);
    cqreg_aligned_free(ipm->delta_lambda_u_aff);
    cqreg_aligned_free(ipm->delta_lambda_v_aff);

    /* Free residuals */
    cqreg_aligned_free(ipm->r_primal);
    cqreg_aligned_free(ipm->r_dual);

    /* Free working arrays */
    cqreg_aligned_free(ipm->work_N);
    cqreg_aligned_free(ipm->work_K);

    /* Free per-thread buffers */
    if (ipm->thread_buf != NULL) {
        for (t = 0; t < ipm->num_threads; t++) {
            cqreg_aligned_free(ipm->thread_buf[t]);
        }
        free(ipm->thread_buf);
    }

    free(ipm);
}

/* ============================================================================
 * Sparsity State Management
 * ============================================================================ */

cqreg_sparsity_state *cqreg_sparsity_create(ST_int N, ST_double quantile, cqreg_bw_method method)
{
    cqreg_sparsity_state *sp = NULL;

    if (N <= 0 || quantile <= 0.0 || quantile >= 1.0) return NULL;

    sp = (cqreg_sparsity_state *)calloc(1, sizeof(cqreg_sparsity_state));
    if (sp == NULL) return NULL;

    sp->N = N;
    sp->quantile = quantile;
    sp->bw_method = method;
    sp->bandwidth = 0.0;
    sp->sparsity = 0.0;

    sp->sorted_resid = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    if (sp->sorted_resid == NULL) {
        free(sp);
        return NULL;
    }

    return sp;
}

void cqreg_sparsity_free(cqreg_sparsity_state *sp)
{
    if (sp == NULL) return;
    cqreg_aligned_free(sp->sorted_resid);
    free(sp);
}

/* ============================================================================
 * Main State Management
 * ============================================================================ */

cqreg_state *cqreg_state_create(ST_int N, ST_int K)
{
    cqreg_state *state = NULL;

    if (N <= 0 || K <= 0) return NULL;

    state = (cqreg_state *)calloc(1, sizeof(cqreg_state));
    if (state == NULL) return NULL;

    state->N = N;
    state->N_original = N;
    state->K = K;
    state->G = 0;
    state->quantile = 0.5;

    /* Determine thread count */
#ifdef _OPENMP
    state->num_threads = omp_get_max_threads();
    if (state->num_threads > CQREG_NUM_THREADS) {
        state->num_threads = CQREG_NUM_THREADS;
    }
#else
    state->num_threads = 1;
#endif

    /* Allocate data arrays */
    state->y = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    state->X = (ST_double *)cqreg_aligned_alloc(N * K * sizeof(ST_double), CQREG_CACHE_LINE);
    if (state->y == NULL || state->X == NULL) goto cleanup;

    /* Allocate results */
    state->beta = (ST_double *)cqreg_aligned_alloc(K * sizeof(ST_double), CQREG_CACHE_LINE);
    state->V = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    state->residuals = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    state->fitted = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    state->obs_density = (ST_double *)cqreg_aligned_alloc(N * sizeof(ST_double), CQREG_CACHE_LINE);
    if (state->beta == NULL || state->V == NULL || state->residuals == NULL ||
        state->fitted == NULL || state->obs_density == NULL) goto cleanup;

    /* Initialize to zero */
    memset(state->beta, 0, K * sizeof(ST_double));
    memset(state->V, 0, K * K * sizeof(ST_double));
    memset(state->residuals, 0, N * sizeof(ST_double));
    memset(state->fitted, 0, N * sizeof(ST_double));
    memset(state->obs_density, 0, N * sizeof(ST_double));

    /* Initialize other fields */
    state->weights = NULL;
    state->obs_mask = NULL;
    state->has_hdfe = 0;
    state->hdfe_state = NULL;
    state->df_a = 0;
    state->sum_adev = 0.0;
    state->sum_rdev = 0.0;
    state->sparsity = 0.0;
    state->bandwidth = 0.0;
    state->vce_type = CQREG_VCE_IID;
    state->bw_method = CQREG_BW_HSHEATHER;
    state->density_method = CQREG_DENSITY_FITTED;  /* Match Stata's default */
    state->cluster_ids = NULL;
    state->num_clusters = 0;
    state->ipm = NULL;
    state->sparsity_state = NULL;
    state->iterations = 0;
    state->converged = 0;

    return state;

cleanup:
    cqreg_state_free(state);
    return NULL;
}

void cqreg_state_free(cqreg_state *state)
{
    if (state == NULL) return;

    /* Free data arrays */
    cqreg_aligned_free(state->y);
    cqreg_aligned_free(state->X);
    cqreg_aligned_free(state->weights);
    cqreg_aligned_free(state->obs_mask);

    /* Free results */
    cqreg_aligned_free(state->beta);
    cqreg_aligned_free(state->V);
    cqreg_aligned_free(state->residuals);
    cqreg_aligned_free(state->fitted);
    cqreg_aligned_free(state->obs_density);

    /* Free cluster IDs */
    cqreg_aligned_free(state->cluster_ids);

    /* Free IPM state */
    cqreg_ipm_free(state->ipm);

    /* Free sparsity state */
    cqreg_sparsity_free(state->sparsity_state);

    /* Note: hdfe_state is owned by creghdfe and should be freed separately */

    free(state);
}

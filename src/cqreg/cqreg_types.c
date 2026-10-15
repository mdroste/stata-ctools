/*
 * cqreg_types.c
 *
 * Memory management for cqreg data structures.
 * Part of the ctools suite.
 */

#include "cqreg_types.h"
#include "../ctools_config.h"
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * IPM Configuration
 * ============================================================================ */

void cqreg_ipm_config_init(cqreg_ipm_config *config)
{
    if (config == NULL) return;

    config->maxiter = CQREG_DEFAULT_MAXITER;
    config->tol_gap = CQREG_DEFAULT_TOL;
    config->verbose = 0;
    config->skip_crossover = 0;
}

/* ============================================================================
 * IPM State Management
 * ============================================================================ */

/*
 * Frisch-Newton workspace of cqreg_fn_solve: 17 N-sized and 2 K-sized arrays.
 * Returns 0, or -1 on overflow or allocation failure; the caller releases
 * whatever was allocated with cqreg_ipm_free().
 */
static int cqreg_ipm_alloc_fn_workspace(cqreg_ipm_state *ipm)
{
    size_t n_size, k_size;

    if (ctools_safe_mul_size((size_t)ipm->N, sizeof(ST_double), &n_size) != 0 ||
        ctools_safe_mul_size((size_t)ipm->K, sizeof(ST_double), &k_size) != 0) {
        return -1;
    }

    ipm->fn_xp = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_s = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_yd = (ST_double *)ctools_cacheline_alloc(k_size);
    ipm->fn_z = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_w = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_dx = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_ds = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_dy = (ST_double *)ctools_cacheline_alloc(k_size);
    ipm->fn_dz = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_dw = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_fx = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_fs = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_fz = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_fw = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_q = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_r = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_tmp = (ST_double *)ctools_cacheline_alloc(n_size);
    ipm->fn_sinv = (ST_double *)ctools_cacheline_alloc(n_size);
    /* The weighted right-hand side, not a scaled design matrix */
    ipm->fn_Xq = (ST_double *)ctools_cacheline_alloc(n_size);

    if (!ipm->fn_xp || !ipm->fn_s || !ipm->fn_yd || !ipm->fn_z || !ipm->fn_w ||
        !ipm->fn_dx || !ipm->fn_ds || !ipm->fn_dy || !ipm->fn_dz || !ipm->fn_dw ||
        !ipm->fn_fx || !ipm->fn_fs || !ipm->fn_fz || !ipm->fn_fw ||
        !ipm->fn_q || !ipm->fn_r || !ipm->fn_tmp || !ipm->fn_sinv || !ipm->fn_Xq)
        return -1;
    return 0;
}

cqreg_ipm_state *cqreg_ipm_create_lite(ST_int N, ST_int K, const cqreg_ipm_config *config)
{
    cqreg_ipm_state *ipm = NULL;

    if (N <= 0 || K <= 0) return NULL;

    /* calloc zeros all pointers to NULL, so cqreg_ipm_free works after a
     * partial allocation and when the crossover arrays are skipped */
    ipm = (cqreg_ipm_state *)calloc(1, sizeof(cqreg_ipm_state));
    if (ipm == NULL) return NULL;

    ipm->N = N;
    ipm->K = K;

    if (config != NULL) {
        memcpy(&ipm->config, config, sizeof(cqreg_ipm_config));
    } else {
        cqreg_ipm_config_init(&ipm->config);
    }

    /* Compute sizes with overflow check */
    size_t kk_size;
    if (ctools_safe_alloc_size((size_t)K, (size_t)K, sizeof(ST_double), &kk_size) != 0) {
        goto cleanup;
    }

    /* beta (K) - stores result */
    ipm->beta = (ST_double *)ctools_cacheline_alloc(K * sizeof(ST_double));
    if (!ipm->beta) goto cleanup;

    /* Normal equations: XDX (K*K), rhs (K); FN reuses the Cholesky factor. */
    ipm->XDX = (ST_double *)ctools_cacheline_alloc(kk_size);
    ipm->rhs = (ST_double *)ctools_cacheline_alloc(K * sizeof(ST_double));
    if (!ipm->XDX || !ipm->rhs) goto cleanup;

    /* r_primal, u, v (N each) and the basis (K) - used by crossover + final residuals */
    if (!config || !config->skip_crossover) {
        ipm->r_primal = (ST_double *)ctools_cacheline_alloc((size_t)N * sizeof(ST_double));
        ipm->u = (ST_double *)ctools_cacheline_alloc((size_t)N * sizeof(ST_double));
        ipm->v = (ST_double *)ctools_cacheline_alloc((size_t)N * sizeof(ST_double));
        ipm->basis = (ST_int *)ctools_cacheline_alloc((size_t)K * sizeof(ST_int));
        if (!ipm->r_primal || !ipm->u || !ipm->v || !ipm->basis) goto cleanup;
    }

    if (cqreg_ipm_alloc_fn_workspace(ipm) != 0) goto cleanup;

    ipm->converged = 0;
    ipm->iterations = 0;

    return ipm;

cleanup:
    cqreg_ipm_free(ipm);
    return NULL;
}

void cqreg_ipm_free(cqreg_ipm_state *ipm)
{
    if (ipm == NULL) return;

    ctools_aligned_free(ipm->beta);
    ctools_aligned_free(ipm->XDX);
    ctools_aligned_free(ipm->rhs);

    /* Crossover results and basis */
    ctools_aligned_free(ipm->u);
    ctools_aligned_free(ipm->v);
    ctools_aligned_free(ipm->r_primal);
    ctools_aligned_free(ipm->basis);

    /* Free Frisch-Newton solver workspace */
    ctools_aligned_free(ipm->fn_xp);
    ctools_aligned_free(ipm->fn_s);
    ctools_aligned_free(ipm->fn_yd);
    ctools_aligned_free(ipm->fn_z);
    ctools_aligned_free(ipm->fn_w);
    ctools_aligned_free(ipm->fn_dx);
    ctools_aligned_free(ipm->fn_ds);
    ctools_aligned_free(ipm->fn_dy);
    ctools_aligned_free(ipm->fn_dz);
    ctools_aligned_free(ipm->fn_dw);
    ctools_aligned_free(ipm->fn_fx);
    ctools_aligned_free(ipm->fn_fs);
    ctools_aligned_free(ipm->fn_fz);
    ctools_aligned_free(ipm->fn_fw);
    ctools_aligned_free(ipm->fn_q);
    ctools_aligned_free(ipm->fn_r);
    ctools_aligned_free(ipm->fn_tmp);
    ctools_aligned_free(ipm->fn_sinv);
    ctools_aligned_free(ipm->fn_Xq);

    free(ipm);
}

ST_double cqreg_ipm_get_objective(const cqreg_ipm_state *ipm, ST_double q)
{
    ST_int N = ipm->N;
    ST_double obj = 0.0;

    for (ST_int i = 0; i < N; i++) {
        obj += q * ipm->u[i] + (1.0 - q) * ipm->v[i];
    }

    return obj;
}

void cqreg_ipm_get_residuals(const cqreg_ipm_state *ipm, ST_double *residuals)
{
    ST_int N = ipm->N;

    /* Residuals = u - v */
    for (ST_int i = 0; i < N; i++) {
        residuals[i] = ipm->u[i] - ipm->v[i];
    }
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

    sp->sorted_resid = (ST_double *)ctools_cacheline_alloc(N * sizeof(ST_double));
    if (sp->sorted_resid == NULL) {
        free(sp);
        return NULL;
    }

    return sp;
}

static void cqreg_sparsity_free(cqreg_sparsity_state *sp)
{
    if (sp == NULL) return;
    ctools_aligned_free(sp->sorted_resid);
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
    state->K = K;
    state->quantile = 0.5;

    /* Determine thread count using unified policy */
    state->num_threads = ctools_get_max_threads();

    /*
     * Data arrays (y, X) are NOT allocated here.
     * They will be set up by ctools_data_load in cqreg_main.c:
     * - state->y points to data->vars[0].data.dbl (view, not owned)
     * - state->X is allocated contiguously and populated from data->vars[1..K-1]
     */
    state->data = NULL;
    state->data_owned = 0;
    state->y = NULL;
    state->X = NULL;
    state->y_owned = 0;
    state->X_owned = 0;

    /* Allocate results */
    state->beta = (ST_double *)ctools_cacheline_alloc(K * sizeof(ST_double));
    state->V = (ST_double *)ctools_cacheline_alloc((size_t)K * K * sizeof(ST_double));
    state->residuals = (ST_double *)ctools_cacheline_alloc(N * sizeof(ST_double));
    state->fitted = (ST_double *)ctools_cacheline_alloc(N * sizeof(ST_double));
    state->obs_density = (ST_double *)ctools_cacheline_alloc(N * sizeof(ST_double));
    if (state->beta == NULL || state->V == NULL || state->residuals == NULL ||
        state->fitted == NULL || state->obs_density == NULL) goto cleanup;

    /* Initialize to zero */
    memset(state->beta, 0, K * sizeof(ST_double));
    memset(state->V, 0, (size_t)K * K * sizeof(ST_double));
    memset(state->residuals, 0, N * sizeof(ST_double));
    memset(state->fitted, 0, N * sizeof(ST_double));
    memset(state->obs_density, 0, N * sizeof(ST_double));

    /* Initialize other fields */
    state->weights = NULL;
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

    /*
     * Free data arrays based on ownership flags.
     * - y_owned: If set, y was separately allocated.
     *            If not set, y is a view into stata_data (don't free).
     * - X is always a contiguous copy we allocated, so always free it.
     */
    if (state->y_owned && state->y != NULL) {
        ctools_aligned_free(state->y);
    }
    /* y is NOT freed if !y_owned because it points into stata_data */

    /* X is always our contiguous allocation */
    if (state->X != NULL) {
        ctools_aligned_free(state->X);
    }

    ctools_aligned_free(state->weights);

    /* Free stata_data if we own it */
    if (state->data_owned && state->data != NULL) {
        stata_data_free(state->data);
        free(state->data);
    }

    /* Free results */
    ctools_aligned_free(state->beta);
    ctools_aligned_free(state->V);
    ctools_aligned_free(state->residuals);
    ctools_aligned_free(state->fitted);
    ctools_aligned_free(state->obs_density);

    /* Free cluster IDs */
    ctools_aligned_free(state->cluster_ids);

    /* Free IPM state */
    cqreg_ipm_free(state->ipm);

    /* Free sparsity state */
    cqreg_sparsity_free(state->sparsity_state);

    free(state);
}

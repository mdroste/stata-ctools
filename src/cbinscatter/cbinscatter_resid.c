/*
 * cbinscatter_resid.c
 *
 * Residualization implementation for cbinscatter
 * Part of the ctools Stata plugin suite
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "stplugin.h"
#include "cbinscatter_resid.h"
#include "../ctools_config.h"
#include "../ctools_ols.h"
#include "../hdfe/hdfe_sample.h"
#include "../hdfe/hdfe_solver.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* Threads for the loops over independent columns (1 without OpenMP) */
static int resid_threads(void)
{
#ifdef _OPENMP
    int nt = ctools_get_openmp_threads();
    return nt > 1 ? nt : 1;
#else
    return 1;
#endif
}

static int resid_thread_id(void)
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

/* Use SF_is_missing() from stplugin.h for missing value checks */

/* ========================================================================
 * Helper: Solve linear system with collinearity handling
 * Uses modified Cholesky to detect collinear columns, drops them,
 * solves reduced system. beta entries for collinear columns are set to 0.
 * Returns 0 on success, non-zero on failure.
 * ======================================================================== */

ST_int cbinscatter_solve_with_collinearity(
    const ST_double *XtX,
    const ST_double *rhs,
    ST_int K,
    ST_double *beta
) {
    /* Normalize the Gram matrix so rank decisions do not depend on units.
     * Modified Cholesky omits dependent columns and solves the retained span. */
    if (K < 1) return CBINSCATTER_ERR_SINGULAR;
    ST_double *L = ctools_safe_malloc3((size_t)K, (size_t)K, sizeof(*L));
    ST_double *scale = ctools_safe_malloc2((size_t)K, sizeof(*scale));
    ST_double *solution = ctools_safe_calloc2((size_t)K, sizeof(*solution));
    ST_int rc = CBINSCATTER_OK;
    if (!L || !scale || !solution) { rc = CBINSCATTER_ERR_MEMORY; goto done; }
    for (ST_int i = 0; i < K; i++) {
        ST_double diagonal = XtX[(size_t)i * K + i];
        if (!isfinite(diagonal) || diagonal < 0 || !isfinite(rhs[i])) {
            rc = CBINSCATTER_ERR_SINGULAR; goto done;
        }
        scale[i] = sqrt(diagonal);
    }
    for (ST_int i = 0; i < K; i++) {
        for (ST_int j = 0; j < K; j++) {
            L[(size_t)i * K + j] = scale[i] && scale[j] ?
                XtX[(size_t)i * K + j] / scale[i] / scale[j] : 0;
            if (!isfinite(L[(size_t)i * K + j])) { rc = CBINSCATTER_ERR_SINGULAR; goto done; }
        }
    }
    const ST_double tolerance = 1e-10;
    for (ST_int k = 0; k < K; k++) {
        ST_double sum = L[(size_t)k * K + k];
        for (ST_int j = 0; j < k; j++) sum -= L[(size_t)k * K + j] * L[(size_t)k * K + j];
        if (sum <= tolerance) {
            if (sum < -tolerance) { rc = CBINSCATTER_ERR_SINGULAR; goto done; }
            L[(size_t)k * K + k] = 0;
            for (ST_int i = k + 1; i < K; i++) L[(size_t)i * K + k] = 0;
            continue;
        }
        L[(size_t)k * K + k] = sqrt(sum);
        for (ST_int i = k + 1; i < K; i++) {
            sum = L[(size_t)i * K + k];
            for (ST_int j = 0; j < k; j++) sum -= L[(size_t)i * K + j] * L[(size_t)k * K + j];
            L[(size_t)i * K + k] = sum / L[(size_t)k * K + k];
        }
    }
    for (ST_int i = 0; i < K; i++) {
        if (!L[(size_t)i * K + i]) continue;
        ST_double sum = rhs[i] / scale[i];
        for (ST_int j = 0; j < i; j++) sum -= L[(size_t)i * K + j] * solution[j];
        solution[i] = sum / L[(size_t)i * K + i];
    }
    for (ST_int i = K - 1; i >= 0; i--) {
        if (!L[(size_t)i * K + i]) continue;
        ST_double sum = solution[i];
        for (ST_int j = i + 1; j < K; j++) sum -= L[(size_t)j * K + i] * solution[j];
        solution[i] = sum / L[(size_t)i * K + i];
        if (!isfinite(solution[i]) || !isfinite(solution[i] / scale[i])) {
            rc = CBINSCATTER_ERR_SINGULAR; goto done;
        }
    }
    for (ST_int i = 0; i < K; i++) beta[i] = scale[i] ? solution[i] / scale[i] : 0;
done:
    free(L); free(scale); free(solution);
    return rc;
}

/* ========================================================================
 * Weighted X'X Computation
 * ======================================================================== */

static ST_retcode compute_xtx_weighted(
    const ST_double *X,
    ST_int N,
    ST_int K,
    const ST_double *weights,
    ST_double *XtX
) {
    ST_int i, j, k;
    ST_double w, xj;

    /* Initialize to zero */
    memset(XtX, 0, (size_t)K * (size_t)K * sizeof(ST_double));

    /* Compute X'WX: element (j, k) sums (w_i x_ij) x_ik over the
     * observations in order; the rows of the upper triangle are independent */
    #pragma omp parallel for schedule(dynamic) private(i, k, w, xj) if (K > 1 && (size_t)N * (size_t)K > 200000)
    for (j = 0; j < K; j++) {
        const ST_double *col_j = X + (size_t)j * N;
        for (k = j; k < K; k++) {
            const ST_double *col_k = X + (size_t)k * N;
            ST_double acc = 0.0;
            for (i = 0; i < N; i++) {
                w = (weights != NULL) ? weights[i] : 1.0;
                xj = col_j[i];
                ST_double wxj = w * xj;
                acc += wxj * col_k[i];
            }
            XtX[(size_t)j * K + k] = acc;
        }
    }

    /* Mirror to lower triangle */
    for (j = 0; j < K; j++) {
        for (k = j + 1; k < K; k++) {
            XtX[(size_t)k * K + j] = XtX[(size_t)j * K + k];
        }
    }

    return CBINSCATTER_OK;
}

/* ========================================================================
 * Weighted normal equations by column pairs (see cbinscatter_resid.h)
 * ======================================================================== */

void cbinscatter_normal_equations(
    const ST_double *const *cols,
    ST_int R,
    const ST_double *y,
    const ST_double *weights,
    const ST_int *bin_ids,
    ST_int num_bins,
    ST_int N,
    ST_double *XtX,
    ST_double *Xty
) {
    ST_int j;

    #pragma omp parallel for schedule(dynamic) if (R > 1 && (size_t)N * (size_t)R > 200000)
    for (j = 0; j < R; j++) {
        const ST_double *xj = cols[j];
        /* right-hand columns: y, then x_j .. x_{R-1}, four per pass */
        const ST_int n_rhs = R - j + 1;
        for (ST_int t = 0; t < n_rhs; t += 4) {
            const ST_double *r[4];
            ST_double acc[4] = {0.0, 0.0, 0.0, 0.0};
            const ST_int n_now = (n_rhs - t < 4) ? n_rhs - t : 4;
            for (ST_int q = 0; q < 4; q++) {
                ST_int rhs = t + (q < n_now ? q : 0);
                r[q] = (rhs == 0) ? y : cols[j + rhs - 1];
            }
            for (ST_int i = 0; i < N; i++) {
                if (bin_ids) {
                    ST_int b = bin_ids[i] - 1;
                    if (b < 0 || b >= num_bins) continue;
                }
                ST_double w = (weights != NULL) ? weights[i] : 1.0;
                ST_double wx = w * xj[i];
                acc[0] += wx * r[0][i];
                acc[1] += wx * r[1][i];
                acc[2] += wx * r[2][i];
                acc[3] += wx * r[3][i];
            }
            for (ST_int q = 0; q < n_now; q++) {
                ST_int rhs = t + q;
                if (rhs == 0) Xty[j] = acc[q];
                else XtX[(size_t)j * R + (j + rhs - 1)] = acc[q];
            }
        }
    }
}

/* ========================================================================
 * OLS Residualization
 * ======================================================================== */

ST_retcode ols_residualize(
    ST_double *y,
    ST_double *x,
    const ST_double *controls,
    ST_int N,
    ST_int K,
    const ST_double *weights,
    ST_int weight_type
) {
    ST_double *X_full = NULL;  /* Controls + constant */
    ST_double *XtX = NULL;
    ST_double *Xty = NULL;
    ST_double *Xtx = NULL;
    ST_double *beta_y = NULL;
    ST_double *beta_x = NULL;
    ST_retcode rc = CBINSCATTER_OK;
    ST_int i, j;
    ST_double w, xj;
    ST_int K_full = K + 1;  /* K controls + 1 constant */

    (void)weight_type;  /* Currently unused, reserved for future */

    /* Allocate working arrays (with overflow-safe multiplication)
     * Include an extra column for the constant term */
    X_full = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_full, sizeof(ST_double));
    XtX = (ST_double *)ctools_safe_malloc3((size_t)K_full, (size_t)K_full, sizeof(ST_double));
    Xty = (ST_double *)ctools_safe_calloc2((size_t)K_full, sizeof(ST_double));
    Xtx = (ST_double *)ctools_safe_calloc2((size_t)K_full, sizeof(ST_double));
    beta_y = (ST_double *)ctools_safe_malloc2((size_t)K_full, sizeof(ST_double));
    beta_x = (ST_double *)ctools_safe_malloc2((size_t)K_full, sizeof(ST_double));

    if (!X_full || !XtX || !Xty || !Xtx || !beta_y || !beta_x) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }

    /* Build design matrix: [constant | controls] */
    /* Column 0: constant (all 1s) */
    for (i = 0; i < N; i++) {
        X_full[0 * N + i] = 1.0;
    }
    /* Columns 1 to K: controls */
    for (j = 0; j < K; j++) {
        for (i = 0; i < N; i++) {
            X_full[(size_t)(j + 1) * N + i] = controls[(size_t)j * N + i];
        }
    }

    /* Compute X'X */
    rc = compute_xtx_weighted(X_full, N, K_full, weights, XtX);
    if (rc != CBINSCATTER_OK) goto cleanup;

    /* Compute X'y and X'x: one independent pair of sums per column */
    #pragma omp parallel for schedule(static) private(i, w, xj) if (K_full > 1 && (size_t)N * (size_t)K_full > 200000)
    for (j = 0; j < K_full; j++) {
        const ST_double *col = X_full + (size_t)j * N;
        ST_double sum_y = 0.0, sum_x = 0.0;
        for (i = 0; i < N; i++) {
            w = (weights != NULL) ? weights[i] : 1.0;
            ST_double wy = w * y[i];
            ST_double wx = w * x[i];
            xj = col[i];
            sum_y += xj * wy;
            sum_x += xj * wx;
        }
        Xty[j] = sum_y;
        Xtx[j] = sum_x;
    }

    /* Solve for beta_y and beta_x, handling collinear controls */
    if (ctools_solve_cholesky(XtX, Xty, K_full, beta_y) != 0) {
        /* Singular matrix - detect and drop collinear columns */
        if (cbinscatter_solve_with_collinearity(XtX, Xty, K_full, beta_y) != 0 ||
            cbinscatter_solve_with_collinearity(XtX, Xtx, K_full, beta_x) != 0) {
            rc = CBINSCATTER_ERR_SINGULAR;
            goto cleanup;
        }
    } else if (ctools_solve_cholesky(XtX, Xtx, K_full, beta_x) != 0) {
        rc = CBINSCATTER_ERR_SINGULAR;
        goto cleanup;
    }

    /* Compute residuals: y - X*beta_y, x - X*beta_x. Each fitted value
     * accumulates over the columns in order; rows are independent, so they
     * run in blocks (a block of the design stays in cache), in parallel. */
    {
        enum { BLOCK = 2048 };
        const ST_int num_blocks = N / BLOCK + (N % BLOCK != 0); /* no N + BLOCK overflow */
        ST_int blk;
        #pragma omp parallel for schedule(static) private(i, j, xj) if ((size_t)N * (size_t)K_full > 200000)
        for (blk = 0; blk < num_blocks; blk++) {
            const ST_int i0 = blk * BLOCK;
            const ST_int i1 = (N - i0 < BLOCK) ? N : i0 + BLOCK;
            ST_double fitted_y[BLOCK], fitted_x[BLOCK];
            for (i = i0; i < i1; i++) {
                fitted_y[i - i0] = 0.0;
                fitted_x[i - i0] = 0.0;
            }
            for (j = 0; j < K_full; j++) {
                const ST_double *col = X_full + (size_t)j * N;
                const ST_double by = beta_y[j], bx = beta_x[j];
                for (i = i0; i < i1; i++) {
                    xj = col[i];
                    fitted_y[i - i0] += xj * by;
                    fitted_x[i - i0] += xj * bx;
                }
            }
            for (i = i0; i < i1; i++) {
                y[i] -= fitted_y[i - i0];
                x[i] -= fitted_x[i - i0];
            }
        }
    }

cleanup:
    free(X_full);
    free(XtX);
    free(Xty);
    free(Xtx);
    free(beta_y);
    free(beta_x);

    return rc;
}

/* ========================================================================
 * HDFE Residualization
 *
 * One absorbed factor is a single exact projection (blocks of rows in
 * parallel); several use the shared conjugate-gradient solver of src/hdfe.
 * ======================================================================== */

/* Structure for FE factor */
typedef struct {
    ST_int num_levels;
    ST_int *levels;      /* Level for each observation (0-based) */
    ST_double *counts;   /* Count per level */
    ST_double *weighted_counts;  /* Weighted count per level */
} FE_Factor_Simple;

/* Initialize factor from FE variable */
static ST_retcode init_factor(
    FE_Factor_Simple *f,
    const ST_int *fe_var,
    const ST_double *weights,
    ST_int N
) {
    ST_int i, level, max_level = 0;

    /* Find max level */
    for (i = 0; i < N; i++) {
        if (fe_var[i] > max_level) max_level = fe_var[i];
    }
    f->num_levels = max_level;

    /* Allocate arrays */
    f->levels = (ST_int *)malloc(N * sizeof(ST_int));
    f->counts = (ST_double *)calloc(max_level, sizeof(ST_double));
    f->weighted_counts = (ST_double *)calloc(max_level, sizeof(ST_double));

    if (!f->levels || !f->counts || !f->weighted_counts) {
        return CBINSCATTER_ERR_MEMORY;
    }

    /* Copy levels (convert to 0-based) and compute counts */
    for (i = 0; i < N; i++) {
        level = fe_var[i] - 1;  /* Convert to 0-based */
        f->levels[i] = level;
        f->counts[level] += 1.0;
        f->weighted_counts[level] += (weights != NULL) ? weights[i] : 1.0;
    }

    return CBINSCATTER_OK;
}

/* Free factor */
static void free_factor(FE_Factor_Simple *f) {
    free(f->levels);
    free(f->counts);
    free(f->weighted_counts);
}

static void free_factors(FE_Factor_Simple *factors, ST_int G) {
    if (!factors) return;
    for (ST_int g = 0; g < G; g++) free_factor(&factors[g]);
    free(factors);
}

/* The G factors of the absorbed variables; *max_levels is the largest
 * number of levels. The caller frees the result with free_factors(). */
static ST_retcode build_factors(
    FE_Factor_Simple **out,
    const ST_int *fe_vars,
    const ST_double *weights,
    ST_int N,
    ST_int G,
    ST_int *max_levels
) {
    FE_Factor_Simple *factors = (FE_Factor_Simple *)calloc(G > 0 ? G : 1, sizeof(FE_Factor_Simple));
    *out = NULL;
    *max_levels = 0;
    if (!factors) return CBINSCATTER_ERR_MEMORY;
    for (ST_int g = 0; g < G; g++) {
        ST_retcode rc = init_factor(&factors[g], &fe_vars[(size_t)g * N], weights, N);
        if (rc != CBINSCATTER_OK) {
            free_factors(factors, G);
            return rc;
        }
        if (factors[g].num_levels > *max_levels) *max_levels = factors[g].num_levels;
    }
    *out = factors;
    return CBINSCATTER_OK;
}

/* Per-block level sums of demean_one_factor() may use at most this much memory */
#define CBS_BLOCKED_DEMEAN_MAX_BYTES ((size_t)256 << 20)

/*
 * Demean ncols columns within the levels of one factor: a single exact
 * projection, in parallel over blocks of rows. Each level's sums are formed
 * per block and added in block order, so results do not depend on the thread
 * count. Returns CBINSCATTER_OK, CBINSCATTER_ERR_MEMORY, or -1 when the
 * per-block sums would need too much memory (the caller falls back).
 */
static ST_retcode demean_one_factor(
    ST_double **cols,
    ST_int ncols,
    const FE_Factor_Simple *f,
    const ST_double *weights,
    ST_int N
) {
    const ST_int num_blocks = cbs_num_blocks(N);
    const ST_int L = f->num_levels;
    const size_t per_block = (size_t)ncols * (size_t)L;
    size_t cells;
    if (ncols < 1 || L < 1) return CBINSCATTER_OK;
    if (ctools_safe_mul_size((size_t)num_blocks, per_block, &cells) != 0 ||
        cells > CBS_BLOCKED_DEMEAN_MAX_BYTES / sizeof(ST_double)) {
        return -1;
    }
    ST_double *block_sums = (ST_double *)malloc(cells * sizeof(ST_double));
    ST_double *means = (ST_double *)malloc(per_block * sizeof(ST_double));
    if (!block_sums || !means) {
        free(block_sums); free(means);
        return CBINSCATTER_ERR_MEMORY;
    }
    ST_int blk, lev;

    /* Level sums per block */
    #pragma omp parallel for schedule(static) if (num_blocks > 1)
    for (blk = 0; blk < num_blocks; blk++) {
        ST_int i0, i1;
        cbs_block_range(blk, N, &i0, &i1);
        for (ST_int c = 0; c < ncols; c++) {
            ST_double *sc = block_sums + (size_t)blk * per_block + (size_t)c * L;
            const ST_double *col = cols[c];
            memset(sc, 0, (size_t)L * sizeof(ST_double));
            if (weights) {
                for (ST_int i = i0; i < i1; i++) sc[f->levels[i]] += weights[i] * col[i];
            } else {
                for (ST_int i = i0; i < i1; i++) sc[f->levels[i]] += col[i];
            }
        }
    }

    /* Level means: block sums added in block order */
    #pragma omp parallel for schedule(static) if ((size_t)L * ncols * num_blocks > 1000000)
    for (lev = 0; lev < L; lev++) {
        for (ST_int c = 0; c < ncols; c++) {
            ST_double sum = 0.0;
            for (ST_int bb = 0; bb < num_blocks; bb++)
                sum += block_sums[(size_t)bb * per_block + (size_t)c * L + lev];
            if (f->weighted_counts[lev] > 0) sum /= f->weighted_counts[lev];
            means[(size_t)c * L + lev] = sum;
        }
    }

    /* Subtract the means */
    #pragma omp parallel for schedule(static) if (num_blocks > 1)
    for (blk = 0; blk < num_blocks; blk++) {
        ST_int i0, i1;
        cbs_block_range(blk, N, &i0, &i1);
        for (ST_int c = 0; c < ncols; c++) {
            const ST_double *mc = means + (size_t)c * L;
            ST_double *col = cols[c];
            for (ST_int i = i0; i < i1; i++) col[i] -= mc[f->levels[i]];
        }
    }

    free(block_sums);
    free(means);
    return CBINSCATTER_OK;
}

/* The shared solver's per-thread buffers (three columns each) may use at
 * most this much memory */
#define CBS_CG_BUFFER_MAX_BYTES ((size_t)1 << 30)

/*
 * Partial G >= 2 absorbed factors out of ncols columns with the shared
 * conjugate-gradient solver (symmetric Kaczmarz preconditioning) that
 * creghdfe uses. Columns run in parallel and each is solved serially, so
 * results do not depend on the thread count. fe_vars is G x N column-major
 * with levels 1..L. Returns CBINSCATTER_OK, CBINSCATTER_ERR_HDFE when a
 * column does not converge within maxiter iterations, or
 * CBINSCATTER_ERR_MEMORY.
 */
static ST_retcode partial_out_shared(
    ST_double **cols,
    ST_int ncols,
    const ST_int *fe_vars,
    ST_int N,
    ST_int G,
    const ST_double *weights,
    ST_int maxiter,
    ST_double tolerance,
    ST_int verbose
) {
    HDFE_State S;
    ST_retcode rc = CBINSCATTER_OK;
    ST_int **levels = NULL, *level_max = NULL;

    memset(&S, 0, sizeof(S));
    if (ncols < 1 || G < 1) return CBINSCATTER_OK;
    /* Without iterations the solver would leave the columns as they are */
    if (maxiter < 1) return CBINSCATTER_ERR_HDFE;
    levels = (ST_int **)malloc((size_t)G * sizeof(ST_int *));
    level_max = (ST_int *)calloc((size_t)G, sizeof(ST_int));
    if (!levels || !level_max) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto done;
    }
    for (ST_int g = 0; g < G; g++) {
        const ST_int *col = &fe_vars[(size_t)g * N];
        ST_int mx = 0;
        for (ST_int i = 0; i < N; i++) if (col[i] > mx) mx = col[i];
        levels[g] = (ST_int *)col;
        level_max[g] = mx;
    }
    ST_retcode build_rc = hdfe_build_factors(&S, G, levels, level_max, N, NULL, N, weights, 0);
    if (build_rc != 0) {
        rc = (build_rc == 920) ? CBINSCATTER_ERR_MEMORY : CBINSCATTER_ERR_INVALID;
        goto done;
    }
    S.has_weights = (weights != NULL);
    S.weights = (ST_double *)weights;
    S.maxiter = maxiter;
    S.tolerance = tolerance;
    S.verbose = verbose;
    S.factors_initialized = 1;
    {
        /* One thread per column, within the memory cap */
        int nt = resid_threads();
        size_t per_thread = 3 * (size_t)N * sizeof(ST_double);
        size_t cap = per_thread ? CBS_CG_BUFFER_MAX_BYTES / per_thread : (size_t)nt;
        if (cap < 1) cap = 1;
        if ((size_t)nt > cap) nt = (int)cap;
        S.num_threads = nt;
    }
    if (ctools_hdfe_alloc_buffers(&S, 0, ncols) != 0) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto done;
    }

    int failed = 0;
    ST_int max_iters = 0, c;
    #pragma omp parallel for num_threads(S.num_threads) schedule(dynamic) \
        reduction(| : failed) reduction(max : max_iters) if (S.num_threads > 1)
    for (c = 0; c < ncols; c++) {
        ST_int iters = cg_solve_column_threaded(&S, cols[c], resid_thread_id());
        if (iters < 0) {
            failed |= 1;
            iters = -iters;
        }
        if (iters > max_iters) max_iters = iters;
    }
    if (verbose) {
        char msg[256];
        if (failed) snprintf(msg, sizeof(msg), "cbinscatter HDFE: max iterations reached\n");
        else snprintf(msg, sizeof(msg), "cbinscatter HDFE converged in %d iterations\n", max_iters);
        SF_display(msg);
    }
    if (failed) rc = CBINSCATTER_ERR_HDFE;

done:
    S.weights = NULL;  /* owned by the caller */
    ctools_hdfe_state_cleanup(&S);
    free(levels);
    free(level_max);
    return rc;
}

/*
 * Partial the G absorbed factors out of ncols columns: one factor by a single
 * block-parallel projection (or by the shared solver when its per-block sums
 * would be too large), several by the shared solver.
 */
static ST_retcode partial_out_factors(
    ST_double **cols,
    ST_int ncols,
    const ST_int *fe_vars,
    ST_int N,
    ST_int G,
    const ST_double *weights,
    ST_int maxiter,
    ST_double tolerance,
    ST_int verbose
) {
    if (ncols < 1 || G < 1) return CBINSCATTER_OK;
    if (G == 1) {
        FE_Factor_Simple *factors = NULL;
        ST_int max_levels = 0;
        ST_retcode rc = build_factors(&factors, fe_vars, weights, N, 1, &max_levels);
        if (rc != CBINSCATTER_OK) return rc;
        rc = demean_one_factor(cols, ncols, &factors[0], weights, N);
        free_factors(factors, 1);
        if (rc != -1) return rc;
    }
    return partial_out_shared(cols, ncols, fe_vars, N, G, weights, maxiter, tolerance, verbose);
}

ST_retcode hdfe_residualize(
    ST_double *y,
    ST_double *x,
    const ST_int *fe_vars,
    ST_int N,
    ST_int G,
    const ST_double *weights,
    ST_int weight_type,
    ST_int maxiter,
    ST_double tolerance,
    ST_int verbose,
    ST_int *dropped
) {
    ST_double *cols[2] = {y, x};

    (void)weight_type;  /* Currently unused */
    *dropped = 0;

    return partial_out_factors(cols, (y == x) ? 1 : 2, fe_vars, N, G, weights, maxiter,
                               tolerance, verbose);
}

/* ========================================================================
 * Y-only OLS Residualization (binsreg method)
 * ======================================================================== */

/* ========================================================================
 * Y-only HDFE Residualization (binsreg method)
 * ======================================================================== */

/* ========================================================================
 * Combined Residualization
 * ======================================================================== */

ST_retcode combined_residualize(
    ST_double *y,
    ST_double *x,
    ST_double *controls,
    const ST_int *fe_vars,
    ST_int N,
    ST_int K_ctrl,
    ST_int G,
    const ST_double *weights,
    ST_int weight_type,
    ST_int maxiter,
    ST_double tolerance,
    ST_int verbose,
    ST_int *dropped
) {
    ST_retcode rc;

    *dropped = 0;

    /* Strategy: first partial the absorbed effects out of y, x and every
     * control together, then OLS on the residuals (FWL) */
    ST_double **cols = (ST_double **)malloc(((size_t)K_ctrl + 2) * sizeof(ST_double *));
    if (!cols) return CBINSCATTER_ERR_MEMORY;
    cols[0] = y;
    cols[1] = x;
    for (ST_int k = 0; k < K_ctrl; k++) cols[2 + k] = &controls[(size_t)k * N];
    rc = partial_out_factors(cols, K_ctrl + 2, fe_vars, N, G, weights, maxiter, tolerance, verbose);
    free(cols);
    if (rc != CBINSCATTER_OK) return rc;

    /* OLS residualize partialled y,x on partialled controls */
    return ols_residualize(y, x, controls, N, K_ctrl, weights, weight_type);
}

/* ========================================================================
 * Batch HDFE Residualization
 *
 * Residualizes multiple variables with the same absorbed factors.
 * ======================================================================== */

ST_retcode hdfe_residualize_batch(
    ST_double **vars,         /* K_vars pointers to N-length arrays to residualize */
    ST_int K_vars,
    const ST_int *fe_vars,
    ST_int N,
    ST_int G,
    const ST_double *weights,
    ST_int weight_type,
    ST_int maxiter,
    ST_double tolerance,
    ST_int verbose,
    ST_int *dropped
) {
    (void)weight_type;
    *dropped = 0;

    if (K_vars == 0) return CBINSCATTER_OK;
    return partial_out_factors(vars, K_vars, fe_vars, N, G, weights, maxiter, tolerance, verbose);
}

/* ========================================================================
 * Combined Y-only Residualization (binsreg method)
 * ======================================================================== */

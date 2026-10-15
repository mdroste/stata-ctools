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
 * Simplified implementation using Gauss-Seidel/Kaczmarz method
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

/* Project vector onto one FE (demean within groups); means is a scratch
 * buffer of at least f->num_levels values */
static void project_one_fe(
    ST_double *y,
    const FE_Factor_Simple *f,
    const ST_double *weights,
    ST_int N,
    ST_double *means
) {
    ST_int i, level;
    ST_double w;

    /* Zero out means */
    memset(means, 0, f->num_levels * sizeof(ST_double));

    /* Accumulate weighted sums */
    for (i = 0; i < N; i++) {
        level = f->levels[i];
        w = (weights != NULL) ? weights[i] : 1.0;
        means[level] += w * y[i];
    }

    /* Divide by weighted counts to get means */
    for (level = 0; level < f->num_levels; level++) {
        if (f->weighted_counts[level] > 0) {
            means[level] /= f->weighted_counts[level];
        }
    }

    /* Subtract means */
    for (i = 0; i < N; i++) {
        level = f->levels[i];
        y[i] -= means[level];
    }
}

/* Kaczmarz sweep (one pass through all FEs) */
static void kaczmarz_sweep(
    ST_double *y,
    const FE_Factor_Simple *factors,
    ST_int G,
    const ST_double *weights,
    ST_int N,
    ST_int forward,
    ST_double *means
) {
    ST_int g;

    if (forward) {
        for (g = 0; g < G; g++) {
            project_one_fe(y, &factors[g], weights, N, means);
        }
    } else {
        for (g = G - 1; g >= 0; g--) {
            project_one_fe(y, &factors[g], weights, N, means);
        }
    }
}

/* Compute norm of vector */
static ST_double vec_norm(const ST_double *y, ST_int N) {
    ST_double sum = 0.0;
    ST_int i;
    for (i = 0; i < N; i++) {
        sum += y[i] * y[i];
    }
    return sqrt(sum);
}

/* One iteration on one column: `sweeps` symmetric (forward, backward)
 * Kaczmarz sweeps, and the relative change of the column over them.
 * prev (N values) and means (max_levels values) are scratch. */
static ST_double demean_iteration(
    ST_double *col,
    const FE_Factor_Simple *factors,
    ST_int G,
    const ST_double *weights,
    ST_int N,
    ST_int sweeps,
    ST_double *prev,
    ST_double *means
) {
    memcpy(prev, col, N * sizeof(ST_double));
    ST_double norm_prev = vec_norm(col, N);

    for (ST_int s = 0; s < sweeps; s++) {
        kaczmarz_sweep(col, factors, G, weights, N, 1, means);  /* Forward */
        kaczmarz_sweep(col, factors, G, weights, N, 0, means);  /* Backward */
    }

    if (!(norm_prev > 0)) return 0.0;
    ST_double diff = 0.0;
    for (ST_int i = 0; i < N; i++) {
        ST_double d = col[i] - prev[i];
        diff += d * d;
    }
    return sqrt(diff) / norm_prev;
}

/* Iterative demeaning (Gauss-Seidel / alternating projections) of one
 * column that is swept twice per iteration: what the y = x call of the
 * pairwise routine below does to it. Returns CBINSCATTER_ERR_HDFE when the
 * iteration limit is reached. */
static ST_retcode demean_column_twice(
    ST_double *col,
    const FE_Factor_Simple *factors,
    ST_int G,
    const ST_double *weights,
    ST_int N,
    ST_int maxiter,
    ST_double tolerance,
    ST_double *prev,
    ST_double *means
) {
    for (ST_int iter = 0; iter < maxiter; iter++) {
        if (demean_iteration(col, factors, G, weights, N, 2, prev, means) < tolerance)
            return CBINSCATTER_OK;
    }
    return CBINSCATTER_ERR_HDFE;
}

/* y and x demeaned together: both are swept until both have converged.
 * scratch holds 2 * N + 2 * max_levels values. */
static ST_retcode demean_pair(
    ST_double *y,
    ST_double *x,
    const FE_Factor_Simple *factors,
    ST_int G,
    const ST_double *weights,
    ST_int N,
    ST_int max_levels,
    ST_int maxiter,
    ST_double tolerance,
    ST_int verbose,
    ST_double *scratch
) {
    ST_double *y_prev = scratch, *x_prev = scratch + N;
    ST_double *y_means = scratch + 2 * (size_t)N, *x_means = y_means + max_levels;
    ST_int iter;

    for (iter = 0; iter < maxiter; iter++) {
        ST_double change_y = 0.0, change_x = 0.0;

        /* The two columns are independent within an iteration */
        #pragma omp parallel sections if (N > 100000)
        {
            #pragma omp section
            change_y = demean_iteration(y, factors, G, weights, N, 1, y_prev, y_means);
            #pragma omp section
            change_x = demean_iteration(x, factors, G, weights, N, 1, x_prev, x_means);
        }

        if (verbose && (iter % 100 == 0 || iter < 10)) {
            char msg[256];
            snprintf(msg, sizeof(msg), "cbinscatter HDFE iter %d: change_y=%.2e change_x=%.2e\n",
                     iter, change_y, change_x);
            SF_display(msg);
        }

        if (change_y < tolerance && change_x < tolerance) {
            if (verbose) {
                char msg[256];
                snprintf(msg, sizeof(msg), "cbinscatter HDFE converged in %d iterations\n", iter + 1);
                SF_display(msg);
            }
            return CBINSCATTER_OK;
        }
    }

    if (verbose) SF_display("cbinscatter HDFE: max iterations reached\n");
    return CBINSCATTER_ERR_HDFE;
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
    FE_Factor_Simple *factors = NULL;
    ST_double *scratch = NULL;
    ST_int max_levels = 0;
    ST_retcode rc;

    (void)weight_type;  /* Currently unused */
    *dropped = 0;

    rc = build_factors(&factors, fe_vars, weights, N, G, &max_levels);
    if (rc != CBINSCATTER_OK) return rc;

    scratch = (ST_double *)ctools_safe_malloc2(2 * ((size_t)N + (size_t)max_levels) + 1, sizeof(ST_double));
    if (!scratch) {
        rc = CBINSCATTER_ERR_MEMORY;
    } else if (y == x) {
        rc = demean_column_twice(y, factors, G, weights, N, maxiter, tolerance,
                                 scratch, scratch + N);
        if (rc == CBINSCATTER_ERR_HDFE && verbose)
            SF_display("cbinscatter HDFE: max iterations reached\n");
    } else {
        rc = demean_pair(y, x, factors, G, weights, N, max_levels, maxiter, tolerance,
                         verbose, scratch);
    }

    free(scratch);
    free_factors(factors, G);
    return rc;
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
    FE_Factor_Simple *factors = NULL;
    ST_double *scratch = NULL;
    ST_int max_levels = 0;

    (void)weight_type;
    *dropped = 0;

    /* Strategy: First HDFE residualize everything, then OLS on residuals.
     * The factors are built once and shared by every column. */
    rc = build_factors(&factors, fe_vars, weights, N, G, &max_levels);
    if (rc != CBINSCATTER_OK) return rc;

    /* Step 1: HDFE residualize y and x */
    scratch = (ST_double *)ctools_safe_malloc2(2 * ((size_t)N + (size_t)max_levels) + 1, sizeof(ST_double));
    if (!scratch) {
        free_factors(factors, G);
        return CBINSCATTER_ERR_MEMORY;
    }
    rc = demean_pair(y, x, factors, G, weights, N, max_levels, maxiter, tolerance, verbose, scratch);
    free(scratch);
    if (rc != CBINSCATTER_OK) {
        free_factors(factors, G);
        return rc;
    }

    /* Step 2: HDFE residualize each control variable. The controls are
     * independent of one another, so they run one per thread, each with its
     * own scratch (N + max_levels values). */
    if (K_ctrl > 0) {
        int nt = resid_threads();
        if (nt > K_ctrl) nt = K_ctrl;
        const size_t per = (size_t)N + (size_t)max_levels + 1;
        ST_double *ctrl_scratch = (ST_double *)ctools_safe_malloc3((size_t)nt, per, sizeof(ST_double));
        while (!ctrl_scratch && nt > 1) {
            /* fewer threads need less scratch */
            nt = nt / 2;
            ctrl_scratch = (ST_double *)ctools_safe_malloc3((size_t)nt, per, sizeof(ST_double));
        }
        if (!ctrl_scratch) {
            free_factors(factors, G);
            return CBINSCATTER_ERR_MEMORY;
        }
        int failed = 0;
        ST_int k;
        #pragma omp parallel for num_threads(nt) schedule(dynamic) reduction(| : failed) if (nt > 1)
        for (k = 0; k < K_ctrl; k++) {
            ST_double *buf = ctrl_scratch + (size_t)resid_thread_id() * per;
            if (demean_column_twice(&controls[(size_t)k * N], factors, G, weights, N, maxiter,
                                    tolerance, buf, buf + N) != CBINSCATTER_OK)
                failed |= 1;
        }
        free(ctrl_scratch);
        if (failed) {
            free_factors(factors, G);
            return CBINSCATTER_ERR_HDFE;
        }
    }
    free_factors(factors, G);

    /* Step 3: OLS residualize partialled y,x on partialled controls */
    rc = ols_residualize(y, x, controls, N, K_ctrl, weights, weight_type);

    return rc;
}

/* ========================================================================
 * Batch HDFE Residualization
 *
 * Residualizes multiple variables simultaneously with shared FE structure.
 * Instead of K separate iteration loops (each doing 2*G*N accesses to
 * FE level arrays), this does a single iteration loop where each sweep
 * demeanes all K variables together, reducing FE level array accesses
 * from K * iters * 2 * G * N to iters * 2 * G * N.
 * ======================================================================== */

/* Demean one variable within the groups of one FE, adding the squared
 * weighted norm of the projection to *change_ss. means is a scratch buffer
 * of at least f->num_levels values. */
static void project_one_fe_tracked(
    ST_double *var,
    const FE_Factor_Simple *f,
    const ST_double *weights,
    ST_int N,
    ST_double *means,
    ST_double *change_ss
) {
    ST_int L = f->num_levels;

    memset(means, 0, (size_t)L * sizeof(ST_double));

    /* Accumulate weighted sums */
    for (ST_int i = 0; i < N; i++) {
        ST_double w = (weights != NULL) ? weights[i] : 1.0;
        means[f->levels[i]] += w * var[i];
    }

    /* Divide by weighted counts to get means */
    for (ST_int level = 0; level < L; level++) {
        if (f->weighted_counts[level] > 0) {
            ST_double inv_w = 1.0 / f->weighted_counts[level];
            means[level] *= inv_w;
            *change_ss += f->weighted_counts[level] * means[level] * means[level];
        }
    }

    /* Subtract means */
    for (ST_int i = 0; i < N; i++) {
        var[i] -= means[f->levels[i]];
    }
}

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
    FE_Factor_Simple *factors = NULL;
    ST_double *change_ss = NULL;
    ST_double *norms = NULL;
    ST_double *means_buf = NULL;
    ST_retcode rc = CBINSCATTER_OK;
    ST_int g, iter, v;
    ST_int max_levels = 0;

    (void)weight_type;
    *dropped = 0;

    if (K_vars == 0) return CBINSCATTER_OK;

    rc = build_factors(&factors, fe_vars, weights, N, G, &max_levels);
    if (rc != CBINSCATTER_OK) return rc;

    /* The variables are demeaned independently, one per thread, each thread
     * with its own buffer of group means. All of them keep iterating until
     * every one has converged, as when they were swept together. */
    int nt = resid_threads();
    if (nt > K_vars) nt = K_vars;
    const size_t per = (size_t)max_levels + 1;
    means_buf = (ST_double *)ctools_safe_malloc3((size_t)nt, per, sizeof(ST_double));
    /* Every column must converge, including bin indicators and controls when
     * the outcome is already orthogonal to the FEs. Keep the original scales
     * so columns fully explained by the FEs can converge to zero as well. */
    change_ss = (ST_double *)ctools_safe_calloc2((size_t)K_vars, sizeof(ST_double));
    norms = (ST_double *)ctools_safe_calloc2((size_t)K_vars, sizeof(ST_double));
    if (!means_buf || !change_ss || !norms) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }
    #pragma omp parallel for num_threads(nt) schedule(static) if (nt > 1 && (size_t)N * (size_t)K_vars > 200000)
    for (v = 0; v < K_vars; v++) {
        ST_double ss = 0.0;
        for (ST_int i = 0; i < N; i++) {
            ST_double w = weights ? weights[i] : 1.0;
            ss += w * vars[v][i] * vars[v][i];
        }
        norms[v] = sqrt(ss);
    }

    /* Iterative demeaning */
    for (iter = 0; iter < maxiter; iter++) {
        /* Forward then backward sweep over the FEs for each variable; the
         * squared norms of its projections are accumulated in that order */
        #pragma omp parallel for num_threads(nt) schedule(dynamic) private(g) if (nt > 1 && (size_t)N * (size_t)K_vars > 200000)
        for (v = 0; v < K_vars; v++) {
            ST_double *means = means_buf + (size_t)resid_thread_id() * per;
            ST_double ss = 0.0;
            for (g = 0; g < G; g++) {
                project_one_fe_tracked(vars[v], &factors[g], weights, N, means, &ss);
            }
            for (g = G - 1; g >= 0; g--) {
                project_one_fe_tracked(vars[v], &factors[g], weights, N, means, &ss);
            }
            change_ss[v] = ss;
        }

        /* For each column the update is the sum of 2*G projections. By
         * Cauchy-Schwarz its squared weighted norm is at most 2*G times the
         * sum of their squared norms. Accumulating those norms from the group
         * means avoids an N*K_vars copy and cancellation in norm differences. */
        int converged = 1;
        for (v = 0; v < K_vars; v++) {
            ST_double diff = 2.0 * G * change_ss[v];
            if (!isfinite(diff) || !isfinite(norms[v]) ||
                (norms[v] > 0 ? sqrt(diff) / norms[v] >= tolerance : diff > 0)) {
                converged = 0;
            }
        }

        if (converged) {
            if (verbose) {
                char msg[256];
                snprintf(msg, sizeof(msg), "cbinscatter HDFE (batch, %d vars) converged in %d iterations\n",
                         K_vars, iter + 1);
                SF_display(msg);
            }
            break;
        }
    }

    if (iter == maxiter) {
        if (verbose) SF_display("cbinscatter HDFE (batch): max iterations reached\n");
        rc = CBINSCATTER_ERR_HDFE;
    }

cleanup:
    free_factors(factors, G);
    free(change_ss);
    free(norms);
    free(means_buf);

    return rc;
}

/* ========================================================================
 * Combined Y-only Residualization (binsreg method)
 * ======================================================================== */

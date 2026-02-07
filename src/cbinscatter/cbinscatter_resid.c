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

/* Use SF_is_missing() from stplugin.h for missing value checks */

/* ========================================================================
 * Helper: Solve linear system with collinearity handling
 * Uses modified Cholesky to detect collinear columns, drops them,
 * solves reduced system. beta entries for collinear columns are set to 0.
 * Returns 0 on success, non-zero on failure.
 * ======================================================================== */

static ST_int solve_with_collinearity(
    const ST_double *XtX,
    const ST_double *rhs,
    ST_int K,
    ST_double *beta
) {
    ST_int i, j, k, a, b;
    ST_int *is_collinear = NULL;
    ST_int *keep = NULL;
    ST_double *L = NULL;
    ST_double *XtX_r = NULL, *rhs_r = NULL, *beta_r = NULL;
    ST_int rc = 0;
    ST_int num_collinear = 0;
    const ST_double tol = 1e-10;  /* Tolerance for collinearity detection */

    is_collinear = (ST_int *)calloc(K, sizeof(ST_int));
    L = (ST_double *)malloc(K * K * sizeof(ST_double));
    if (!is_collinear || !L) {
        free(is_collinear);
        free(L);
        return -1;
    }

    /* Inline modified Cholesky to detect collinear columns */
    memcpy(L, XtX, K * K * sizeof(ST_double));

    for (k = 0; k < K; k++) {
        ST_double sum = L[k * K + k];
        for (j = 0; j < k; j++) {
            sum -= L[k * K + j] * L[k * K + j];
        }

        if (sum < tol) {
            /* Column k is collinear */
            is_collinear[k] = 1;
            num_collinear++;
            L[k * K + k] = 0.0;
            for (i = k + 1; i < K; i++) {
                L[i * K + k] = 0.0;
            }
            continue;
        }

        ST_double diag = sqrt(sum);
        L[k * K + k] = diag;

        for (i = k + 1; i < K; i++) {
            sum = L[i * K + k];
            for (j = 0; j < k; j++) {
                sum -= L[i * K + j] * L[k * K + j];
            }
            L[i * K + k] = sum / diag;
        }
    }

    free(L);

    ST_int K_r = K - num_collinear;
    if (K_r == 0) {
        memset(beta, 0, K * sizeof(ST_double));
        free(is_collinear);
        return 0;
    }

    /* If no collinearity detected, solve directly */
    if (num_collinear == 0) {
        free(is_collinear);
        return ctools_solve_cholesky(XtX, rhs, K, beta);
    }

    /* Build reduced system excluding collinear columns */
    keep = (ST_int *)malloc(K_r * sizeof(ST_int));
    XtX_r = (ST_double *)malloc((size_t)K_r * K_r * sizeof(ST_double));
    rhs_r = (ST_double *)malloc(K_r * sizeof(ST_double));
    beta_r = (ST_double *)malloc(K_r * sizeof(ST_double));

    if (!keep || !XtX_r || !rhs_r || !beta_r) {
        rc = -1;
        goto done;
    }

    ST_int idx = 0;
    for (j = 0; j < K; j++) {
        if (!is_collinear[j]) keep[idx++] = j;
    }

    for (a = 0; a < K_r; a++) {
        rhs_r[a] = rhs[keep[a]];
        for (b = 0; b < K_r; b++) {
            XtX_r[a * K_r + b] = XtX[keep[a] * K + keep[b]];
        }
    }

    rc = ctools_solve_cholesky(XtX_r, rhs_r, K_r, beta_r);

    if (rc == 0) {
        memset(beta, 0, K * sizeof(ST_double));
        for (a = 0; a < K_r; a++) {
            beta[keep[a]] = beta_r[a];
        }
    }

done:
    free(is_collinear);
    free(keep);
    free(XtX_r);
    free(rhs_r);
    free(beta_r);
    return rc;
}

/* ========================================================================
 * Weighted X'X Computation
 * ======================================================================== */

ST_retcode compute_xtx_weighted(
    const ST_double *X,
    ST_int N,
    ST_int K,
    const ST_double *weights,
    ST_double *XtX
) {
    ST_int i, j, k;
    ST_double w, xj;

    /* Initialize to zero */
    memset(XtX, 0, K * K * sizeof(ST_double));

    /* Compute X'WX */
    for (i = 0; i < N; i++) {
        w = (weights != NULL) ? weights[i] : 1.0;
        for (j = 0; j < K; j++) {
            xj = X[j * N + i];  /* Column-major */
            ST_double wxj = w * xj;
            for (k = j; k < K; k++) {
                XtX[j * K + k] += wxj * X[k * N + i];
            }
        }
    }

    /* Mirror to lower triangle */
    for (j = 0; j < K; j++) {
        for (k = j + 1; k < K; k++) {
            XtX[k * K + j] = XtX[j * K + k];
        }
    }

    return CBINSCATTER_OK;
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
            X_full[(j + 1) * N + i] = controls[j * N + i];
        }
    }

    /* Compute X'X */
    rc = compute_xtx_weighted(X_full, N, K_full, weights, XtX);
    if (rc != CBINSCATTER_OK) goto cleanup;

    /* Compute X'y and X'x */
    for (i = 0; i < N; i++) {
        w = (weights != NULL) ? weights[i] : 1.0;
        ST_double wy = w * y[i];
        ST_double wx = w * x[i];
        for (j = 0; j < K_full; j++) {
            xj = X_full[j * N + i];
            Xty[j] += xj * wy;
            Xtx[j] += xj * wx;
        }
    }

    /* Solve for beta_y and beta_x, handling collinear controls */
    if (ctools_solve_cholesky(XtX, Xty, K_full, beta_y) != 0) {
        /* Singular matrix - detect and drop collinear columns */
        if (solve_with_collinearity(XtX, Xty, K_full, beta_y) != 0 ||
            solve_with_collinearity(XtX, Xtx, K_full, beta_x) != 0) {
            rc = CBINSCATTER_ERR_SINGULAR;
            goto cleanup;
        }
    } else if (ctools_solve_cholesky(XtX, Xtx, K_full, beta_x) != 0) {
        rc = CBINSCATTER_ERR_SINGULAR;
        goto cleanup;
    }

    /* Compute residuals: y - X*beta_y, x - X*beta_x */
    for (i = 0; i < N; i++) {
        ST_double fitted_y = 0.0, fitted_x = 0.0;
        for (j = 0; j < K_full; j++) {
            xj = X_full[j * N + i];
            fitted_y += xj * beta_y[j];
            fitted_x += xj * beta_x[j];
        }
        y[i] -= fitted_y;
        x[i] -= fitted_x;
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
    ST_double *means;    /* Mean per level (working buffer) */
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
    f->means = (ST_double *)calloc(max_level, sizeof(ST_double));

    if (!f->levels || !f->counts || !f->weighted_counts || !f->means) {
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
    free(f->means);
}

/* Project vector onto one FE (demean within groups) */
static void project_one_fe(
    ST_double *y,
    const FE_Factor_Simple *f,
    const ST_double *weights,
    ST_int N
) {
    ST_int i, level;
    ST_double w;

    /* Zero out means */
    memset(f->means, 0, f->num_levels * sizeof(ST_double));

    /* Accumulate weighted sums */
    for (i = 0; i < N; i++) {
        level = f->levels[i];
        w = (weights != NULL) ? weights[i] : 1.0;
        f->means[level] += w * y[i];
    }

    /* Divide by weighted counts to get means */
    for (level = 0; level < f->num_levels; level++) {
        if (f->weighted_counts[level] > 0) {
            f->means[level] /= f->weighted_counts[level];
        }
    }

    /* Subtract means */
    for (i = 0; i < N; i++) {
        level = f->levels[i];
        y[i] -= f->means[level];
    }
}

/* Kaczmarz sweep (one pass through all FEs) */
static void kaczmarz_sweep(
    ST_double *y,
    FE_Factor_Simple *factors,
    ST_int G,
    const ST_double *weights,
    ST_int N,
    ST_int forward
) {
    ST_int g;

    if (forward) {
        for (g = 0; g < G; g++) {
            project_one_fe(y, &factors[g], weights, N);
        }
    } else {
        for (g = G - 1; g >= 0; g--) {
            project_one_fe(y, &factors[g], weights, N);
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
    ST_double *y_prev = NULL, *x_prev = NULL;
    ST_retcode rc = CBINSCATTER_OK;
    ST_int g, iter;
    ST_double norm_y_prev, norm_x_prev;
    ST_double change_y, change_x;

    (void)weight_type;  /* Currently unused */
    *dropped = 0;

    /* Allocate factors */
    factors = (FE_Factor_Simple *)calloc(G, sizeof(FE_Factor_Simple));
    if (!factors) return CBINSCATTER_ERR_MEMORY;

    /* Initialize each factor */
    for (g = 0; g < G; g++) {
        rc = init_factor(&factors[g], &fe_vars[g * N], weights, N);
        if (rc != CBINSCATTER_OK) goto cleanup;
    }

    /* Allocate previous iteration buffers */
    y_prev = (ST_double *)malloc(N * sizeof(ST_double));
    x_prev = (ST_double *)malloc(N * sizeof(ST_double));
    if (!y_prev || !x_prev) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }

    /* Iterative demeaning (Gauss-Seidel / alternating projections) */
    for (iter = 0; iter < maxiter; iter++) {
        /* Save previous values */
        memcpy(y_prev, y, N * sizeof(ST_double));
        memcpy(x_prev, x, N * sizeof(ST_double));
        norm_y_prev = vec_norm(y, N);
        norm_x_prev = vec_norm(x, N);

        /* Forward and backward sweep for symmetric Kaczmarz */
        kaczmarz_sweep(y, factors, G, weights, N, 1);  /* Forward */
        kaczmarz_sweep(y, factors, G, weights, N, 0);  /* Backward */
        kaczmarz_sweep(x, factors, G, weights, N, 1);
        kaczmarz_sweep(x, factors, G, weights, N, 0);

        /* Check convergence - compute relative change */
        if (norm_y_prev > 0) {
            ST_double diff_y = 0.0;
            for (ST_int i = 0; i < N; i++) {
                ST_double d = y[i] - y_prev[i];
                diff_y += d * d;
            }
            change_y = sqrt(diff_y) / norm_y_prev;
        } else {
            change_y = 0.0;
        }

        if (norm_x_prev > 0) {
            ST_double diff_x = 0.0;
            for (ST_int i = 0; i < N; i++) {
                ST_double d = x[i] - x_prev[i];
                diff_x += d * d;
            }
            change_x = sqrt(diff_x) / norm_x_prev;
        } else {
            change_x = 0.0;
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
            break;
        }
    }

    if (iter == maxiter && verbose) {
        SF_display("cbinscatter HDFE: max iterations reached\n");
    }

cleanup:
    if (factors) {
        for (g = 0; g < G; g++) {
            free_factor(&factors[g]);
        }
        free(factors);
    }
    free(y_prev);
    free(x_prev);

    return rc;
}

/* ========================================================================
 * Y-only OLS Residualization (binsreg method)
 * ======================================================================== */

ST_retcode ols_residualize_y_only(
    ST_double *y,
    const ST_double *controls,
    ST_int N,
    ST_int K,
    const ST_double *weights,
    ST_int weight_type
) {
    ST_double *XtX = NULL;
    ST_double *Xty = NULL;
    ST_double *beta_y = NULL;
    ST_retcode rc = CBINSCATTER_OK;
    ST_int i, j;
    ST_double w, xj;

    (void)weight_type;  /* Currently unused, reserved for future */

    /* Allocate working arrays */
    XtX = (ST_double *)malloc((size_t)K * K * sizeof(ST_double));
    Xty = (ST_double *)calloc(K, sizeof(ST_double));
    beta_y = (ST_double *)malloc(K * sizeof(ST_double));

    if (!XtX || !Xty || !beta_y) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }

    /* Compute X'X */
    rc = compute_xtx_weighted(controls, N, K, weights, XtX);
    if (rc != CBINSCATTER_OK) goto cleanup;

    /* Compute X'y */
    for (i = 0; i < N; i++) {
        w = (weights != NULL) ? weights[i] : 1.0;
        ST_double wy = w * y[i];
        for (j = 0; j < K; j++) {
            xj = controls[j * N + i];
            Xty[j] += xj * wy;
        }
    }

    /* Solve for beta_y: (X'X) * beta_y = X'y, handling collinear controls */
    if (ctools_solve_cholesky(XtX, Xty, K, beta_y) != 0) {
        /* Singular matrix - detect and drop collinear columns */
        if (solve_with_collinearity(XtX, Xty, K, beta_y) != 0) {
            rc = CBINSCATTER_ERR_SINGULAR;
            goto cleanup;
        }
    }

    /* Compute residuals: y - X*beta_y */
    for (i = 0; i < N; i++) {
        ST_double fitted_y = 0.0;
        for (j = 0; j < K; j++) {
            xj = controls[j * N + i];
            fitted_y += xj * beta_y[j];
        }
        y[i] -= fitted_y;
    }

cleanup:
    free(XtX);
    free(Xty);
    free(beta_y);

    return rc;
}

/* ========================================================================
 * Y-only HDFE Residualization (binsreg method)
 * ======================================================================== */

ST_retcode hdfe_residualize_y_only(
    ST_double *y,
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
    ST_double *y_prev = NULL;
    ST_retcode rc = CBINSCATTER_OK;
    ST_int g, iter;
    ST_double norm_y_prev;
    ST_double change_y;

    (void)weight_type;  /* Currently unused */
    *dropped = 0;

    /* Allocate factors */
    factors = (FE_Factor_Simple *)calloc(G, sizeof(FE_Factor_Simple));
    if (!factors) return CBINSCATTER_ERR_MEMORY;

    /* Initialize each factor */
    for (g = 0; g < G; g++) {
        rc = init_factor(&factors[g], &fe_vars[g * N], weights, N);
        if (rc != CBINSCATTER_OK) goto cleanup;
    }

    /* Allocate previous iteration buffer */
    y_prev = (ST_double *)malloc(N * sizeof(ST_double));
    if (!y_prev) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }

    /* Iterative demeaning (Gauss-Seidel / alternating projections) */
    for (iter = 0; iter < maxiter; iter++) {
        /* Save previous values */
        memcpy(y_prev, y, N * sizeof(ST_double));
        norm_y_prev = vec_norm(y, N);

        /* Forward and backward sweep for symmetric Kaczmarz */
        kaczmarz_sweep(y, factors, G, weights, N, 1);  /* Forward */
        kaczmarz_sweep(y, factors, G, weights, N, 0);  /* Backward */

        /* Check convergence - compute relative change */
        if (norm_y_prev > 0) {
            ST_double diff_y = 0.0;
            for (ST_int i = 0; i < N; i++) {
                ST_double d = y[i] - y_prev[i];
                diff_y += d * d;
            }
            change_y = sqrt(diff_y) / norm_y_prev;
        } else {
            change_y = 0.0;
        }

        if (verbose && (iter % 100 == 0 || iter < 10)) {
            char msg[256];
            snprintf(msg, sizeof(msg), "cbinscatter HDFE (y-only) iter %d: change_y=%.2e\n",
                     iter, change_y);
            SF_display(msg);
        }

        if (change_y < tolerance) {
            if (verbose) {
                char msg[256];
                snprintf(msg, sizeof(msg), "cbinscatter HDFE (y-only) converged in %d iterations\n", iter + 1);
                SF_display(msg);
            }
            break;
        }
    }

    if (iter == maxiter && verbose) {
        SF_display("cbinscatter HDFE (y-only): max iterations reached\n");
    }

cleanup:
    if (factors) {
        for (g = 0; g < G; g++) {
            free_factor(&factors[g]);
        }
        free(factors);
    }
    free(y_prev);

    return rc;
}

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

    /* Strategy: First HDFE residualize everything, then OLS on residuals */

    /* Step 1: HDFE residualize y and x */
    rc = hdfe_residualize(y, x, fe_vars, N, G, weights, weight_type,
                          maxiter, tolerance, verbose, dropped);
    if (rc != CBINSCATTER_OK) return rc;

    /* Step 2: HDFE residualize each control variable */
    for (ST_int k = 0; k < K_ctrl; k++) {
        ST_double *ctrl_col = &controls[k * N];
        ST_int dummy_dropped;
        rc = hdfe_residualize(ctrl_col, ctrl_col, fe_vars, N, G, weights,
                              weight_type, maxiter, tolerance, 0, &dummy_dropped);
        if (rc != CBINSCATTER_OK) return rc;
    }

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

/* Batch version of project_one_fe: demean K variables within groups of one FE */
static void project_one_fe_batch(
    ST_double **vars,         /* K pointers to N-length arrays */
    ST_int K_vars,
    const FE_Factor_Simple *f,
    const ST_double *weights,
    ST_int N,
    ST_double *means_buf      /* num_levels * K_vars scratch buffer */
) {
    ST_int L = f->num_levels;

    memset(means_buf, 0, (size_t)L * K_vars * sizeof(ST_double));

    /* Accumulate weighted sums for all variables in one pass */
    for (ST_int i = 0; i < N; i++) {
        ST_int level = f->levels[i];
        ST_double w = (weights != NULL) ? weights[i] : 1.0;
        ST_double *row = &means_buf[(size_t)level * K_vars];
        for (ST_int v = 0; v < K_vars; v++) {
            row[v] += w * vars[v][i];
        }
    }

    /* Divide by weighted counts to get means */
    for (ST_int level = 0; level < L; level++) {
        if (f->weighted_counts[level] > 0) {
            ST_double inv_w = 1.0 / f->weighted_counts[level];
            ST_double *row = &means_buf[(size_t)level * K_vars];
            for (ST_int v = 0; v < K_vars; v++) {
                row[v] *= inv_w;
            }
        }
    }

    /* Subtract means for all variables in one pass */
    for (ST_int i = 0; i < N; i++) {
        ST_int level = f->levels[i];
        ST_double *row = &means_buf[(size_t)level * K_vars];
        for (ST_int v = 0; v < K_vars; v++) {
            vars[v][i] -= row[v];
        }
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
    ST_double *prev_buf = NULL;
    ST_double *means_buf = NULL;
    ST_retcode rc = CBINSCATTER_OK;
    ST_int g, iter;

    (void)weight_type;
    *dropped = 0;

    if (K_vars == 0) return CBINSCATTER_OK;

    /* Allocate factors */
    factors = (FE_Factor_Simple *)calloc(G, sizeof(FE_Factor_Simple));
    if (!factors) return CBINSCATTER_ERR_MEMORY;

    /* Initialize each factor and find max levels for means buffer */
    ST_int max_levels = 0;
    for (g = 0; g < G; g++) {
        rc = init_factor(&factors[g], &fe_vars[g * N], weights, N);
        if (rc != CBINSCATTER_OK) goto cleanup;
        if (factors[g].num_levels > max_levels)
            max_levels = factors[g].num_levels;
    }

    /* Allocate means buffer: max_levels * K_vars */
    means_buf = (ST_double *)malloc((size_t)max_levels * K_vars * sizeof(ST_double));
    /* Allocate buffer for convergence check: K_vars previous norms */
    prev_buf = (ST_double *)malloc(N * sizeof(ST_double));
    if (!means_buf || !prev_buf) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }

    /* Iterative demeaning */
    for (iter = 0; iter < maxiter; iter++) {
        /* Compute max norm before sweep (use first variable as representative) */
        ST_double norm_prev = 0.0;
        for (ST_int i = 0; i < N; i++) {
            norm_prev += vars[0][i] * vars[0][i];
        }
        norm_prev = sqrt(norm_prev);
        memcpy(prev_buf, vars[0], N * sizeof(ST_double));

        /* Forward sweep: demean all variables for each FE */
        for (g = 0; g < G; g++) {
            project_one_fe_batch(vars, K_vars, &factors[g], weights, N, means_buf);
        }
        /* Backward sweep */
        for (g = G - 1; g >= 0; g--) {
            project_one_fe_batch(vars, K_vars, &factors[g], weights, N, means_buf);
        }

        /* Check convergence using first variable */
        ST_double change = 0.0;
        if (norm_prev > 0) {
            ST_double diff = 0.0;
            for (ST_int i = 0; i < N; i++) {
                ST_double d = vars[0][i] - prev_buf[i];
                diff += d * d;
            }
            change = sqrt(diff) / norm_prev;
        }

        if (change < tolerance) {
            if (verbose) {
                char msg[256];
                snprintf(msg, sizeof(msg), "cbinscatter HDFE (batch, %d vars) converged in %d iterations\n",
                         K_vars, iter + 1);
                SF_display(msg);
            }
            break;
        }
    }

cleanup:
    if (factors) {
        for (g = 0; g < G; g++) {
            free_factor(&factors[g]);
        }
        free(factors);
    }
    free(prev_buf);
    free(means_buf);

    return rc;
}

/* ========================================================================
 * Combined Y-only Residualization (binsreg method)
 * ======================================================================== */

ST_retcode combined_residualize_y_only(
    ST_double *y,
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

    /* Strategy: First HDFE residualize y and controls, then OLS on residuals */

    /* Step 1: HDFE residualize y only */
    rc = hdfe_residualize_y_only(y, fe_vars, N, G, weights, weight_type,
                                  maxiter, tolerance, verbose, dropped);
    if (rc != CBINSCATTER_OK) return rc;

    /* Step 2: HDFE residualize each control variable */
    for (ST_int k = 0; k < K_ctrl; k++) {
        ST_double *ctrl_col = &controls[k * N];
        ST_int dummy_dropped;
        /* Use y_only version but pass ctrl_col for both - it only modifies first arg */
        rc = hdfe_residualize_y_only(ctrl_col, fe_vars, N, G, weights,
                                      weight_type, maxiter, tolerance, 0, &dummy_dropped);
        if (rc != CBINSCATTER_OK) return rc;
    }

    /* Step 3: OLS residualize partialled y on partialled controls (y only) */
    rc = ols_residualize_y_only(y, controls, N, K_ctrl, weights, weight_type);

    return rc;
}

/*
 * creghdfe_hdfe.c
 *
 * Global HDFE state (g_state) and its cleanup
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_hdfe.h"
#include "hdfe/hdfe_remap.h"
#include "../ctools_hdfe_utils.h"
#include "../ctools_config.h"
#include "../ctools_spi.h"  /* Error-checking SPI wrappers */

/* Define the global state pointer */
HDFE_State *g_state = NULL;

ST_retcode (*creghdfe_boot_retain_hook)(const struct creghdfe_boot_source_s *) = NULL;
void (*creghdfe_boot_release_hook)(void) = NULL;

/*
 * Clean up global state
 */
void cleanup_state(void)
{
    if (g_state == NULL) return;

    /* Free weights (ctools_hdfe_state_cleanup does not free weights) */
    if (g_state->weights != NULL) free(g_state->weights);
    g_state->weights = NULL;

    /* Free factors, thread buffers, etc. */
    ctools_hdfe_state_cleanup(g_state);

    free(g_state);
    g_state = NULL;
}

/* ============================================================================
 * Cleanup function for ctools_cleanup system
 * ============================================================================ */

void creghdfe_cleanup_state(void)
{
    /* Reuse existing cleanup function */
    cleanup_state();
    /* and the model retained for vce(bootstrap) */
    if (creghdfe_boot_release_hook) creghdfe_boot_release_hook();
}

/* ============================================================================
 * Heterogeneous slopes (see creghdfe_hdfe.h)
 * ============================================================================ */

#include <float.h>
#include "../ctools_ols.h" /* ctools_invsym */

/* The intercept factor of factor g's absorbed effect, or -1 (slopes only) */
static ST_int effect_intercept(const creghdfe_input *in, ST_int G, ST_int g)
{
    for (ST_int h = 0; h < G; h++)
        if (in->fe_absvar[h] == in->fe_absvar[g] && in->slope_idx[h] == 0)
            return h;
    return -1;
}

/* reghdfe's scale of a slope variable (compute_stdevs): the standard
 * deviation from uncentered moments, at least 1e-3 */
static ST_double slope_stdev(const ST_double *x, ST_int N)
{
    ST_double s = 0.0, ss = 0.0;
    for (ST_int i = 0; i < N; i++) {
        s += x[i];
        ss += x[i] * x[i];
    }
    ST_double v = N > 1 ? (ss - s * s / N) / (N - 1) : 0.0;
    ST_double sd = v > 0.0 ? sqrt(v) : 0.0;
    return sd > 1e-3 ? sd : 1e-3;
}

ST_retcode creghdfe_setup_slopes(HDFE_State *S, const creghdfe_input *in,
                                 const ST_double *const *slope_cols, const ST_int *mask,
                                 ST_int N_orig, const ST_double *weights)
{
    const ST_int N = S->N, G = S->G;
    ST_retcode rc = 0;
    if (!in->num_slopes) return 0;
    ST_double *w = NULL, *sw = NULL, *ss = NULL;
    /* The solver's weights: aweights and pweights sum to N */
    if (weights) {
        w = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(*w));
        if (!w) return 920;
        ST_double sum = 0.0;
        for (ST_int row = 0, i = 0; row < N_orig; row++)
            if (mask[row]) sum += (w[i++] = weights[row]);
        if (S->weight_type == 1 || S->weight_type == 3) {
            if (!(sum > 0.0)) { free(w); return 498; }
            for (ST_int i = 0; i < N; i++) w[i] *= N / sum;
        }
    }
    for (ST_int g = 0; g < G && !rc; g++) {
        if (!in->slope_idx[g]) continue;
        FE_Factor *f = &S->factors[g];
        const ST_int L = f->num_levels;
        f->slope = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(*f->slope));
        sw = (ST_double *)calloc((size_t)L, sizeof(*sw));
        ss = (ST_double *)calloc((size_t)L, sizeof(*ss));
        if (!f->slope || !sw || !ss || !slope_cols[g]) { rc = 920; break; }
        for (ST_int row = 0, i = 0; row < N_orig; row++)
            if (mask[row]) f->slope[i++] = slope_cols[g][row];
        f->has_intercept = 0;
        f->num_slopes = 1;
        if (effect_intercept(in, G, g) >= 0) {
            f->slope_center = (ST_double *)calloc((size_t)L, sizeof(*f->slope_center));
            if (!f->slope_center) { rc = 920; break; }
            for (ST_int i = 0; i < N; i++) {
                ST_int l = f->levels[i] - 1;
                ST_double wi = w ? w[i] : 1.0;
                sw[l] += wi;
                f->slope_center[l] += wi * f->slope[i];
            }
            for (ST_int l = 0; l < L; l++)
                f->slope_center[l] = sw[l] > 0.0 ? f->slope_center[l] / sw[l] : 0.0;
        }
        /* Sums of squares of the loading: in the solver's weights (for the
         * test below), unweighted (counts) and in the loaded weights
         * (weighted_counts, normalized with the other factors later) */
        if (f->weighted_counts) memset(f->weighted_counts, 0, (size_t)L * sizeof(ST_double));
        memset(f->counts, 0, (size_t)L * sizeof(ST_double));
        for (ST_int row = 0, i = 0; row < N_orig; row++) {
            if (!mask[row]) continue;
            ST_int l = f->levels[i] - 1;
            ST_double x = f->slope[i] - (f->slope_center ? f->slope_center[l] : 0.0);
            ss[l] += (w ? w[i] : 1.0) * x * x;
            f->counts[l] += x * x;
            if (f->weighted_counts) f->weighted_counts[l] += weights[row] * x * x;
            i++;
        }
        /* reghdfe inverts each level's cross product of the standardized
         * slope with invsym(), which drops it below 1e-19 */
        ST_double sd = slope_stdev(f->slope, N);
        for (ST_int l = 0; l < L; l++) {
            if (!(ss[l] / (sd * sd) >= 1e-19)) {
                f->counts[l] = 0.0;
                if (f->weighted_counts) f->weighted_counts[l] = 0.0;
            }
        }
        free(sw); free(ss);
        sw = ss = NULL;
    }
    free(sw); free(ss); free(w);
    return rc;
}

ST_retcode creghdfe_redundant(const HDFE_State *S, const creghdfe_input *in, const ST_int *nested,
                              ST_int dof_pairs, ST_int *redundant, ST_int *exact)
{
    const ST_int N = S->N, G = S->G;
    FE_Factor sub[10];
    ST_int map[10], sub_nested[10], sub_redundant[10], n = 0;
    if (G > 10) return 198;
    for (ST_int g = 0; g < G; g++) {
        redundant[g] = 0;
        if (exact) exact[g] = 0;
        if (in->slope_idx[g]) continue;
        sub[n] = S->factors[g];
        sub_nested[n] = nested ? nested[g] : 0;
        map[n++] = g;
    }
    if (ctools_hdfe_redundant_intercepts(sub, n, N, nested ? sub_nested : NULL, dof_pairs,
                                         sub_redundant))
        return 920;
    ST_int first = -1, second = -1;
    for (ST_int i = 0; i < n; i++) {
        ST_int g = map[i];
        redundant[g] = sub_redundant[i];
        if (nested && nested[g]) {
            if (exact) exact[g] = 1;
        } else if (first < 0) {
            first = g;
        } else if (second < 0) {
            second = g;
        }
    }
    if (exact && first >= 0) exact[first] = 1;
    if (exact && dof_pairs && second >= 0) exact[second] = 1;
    if (!in->dof_continuous) return 0;

    /* Slopes: levels where the standardized slope is constant (with an
     * intercept) or zero (without), with reghdfe's tolerance of 1e-6 */
    for (ST_int g = 0; g < G; g++) {
        if (!in->slope_idx[g]) continue;
        const FE_Factor *f = &S->factors[g];
        const ST_int L = f->num_levels;
        ST_double *lo = (ST_double *)malloc((size_t)L * sizeof(*lo));
        ST_double *hi = (ST_double *)malloc((size_t)L * sizeof(*hi));
        if (!lo || !hi) { free(lo); free(hi); return 920; }
        for (ST_int l = 0; l < L; l++) { lo[l] = DBL_MAX; hi[l] = -DBL_MAX; }
        ST_double sd = slope_stdev(f->slope, N);
        for (ST_int i = 0; i < N; i++) {
            ST_int l = f->levels[i] - 1;
            ST_double x = f->slope[i] / sd;
            if (x < lo[l]) lo[l] = x;
            if (x > hi[l]) hi[l] = x;
        }
        ST_int has_int = effect_intercept(in, G, g) >= 0;
        for (ST_int l = 0; l < L; l++)
            redundant[g] += has_int ? (hi[l] - lo[l] <= 1e-6) : (fabs(lo[l]) + fabs(hi[l]) <= 1e-6);
        free(lo); free(hi);
    }
    return 0;
}

/* One absorbed effect for the recovery: its intercept factor (or -1), its
 * k slope factors, and per level the (weighted) count, the means of the
 * standardized slopes (with an intercept) and the generalized inverse of
 * their cross product */
typedef struct {
    ST_int intercept, k, slopes[10], L;
    const ST_int *levels;
    ST_double sd[10];
    ST_double *wsum, *xbar, *inv;
} creghdfe_effect;

/* Slope j of effect f at row i, standardized and centered */
static inline ST_double effect_x(const HDFE_State *S, const creghdfe_effect *f, ST_int j, ST_int i,
                                 ST_int l)
{
    return S->factors[f->slopes[j]].slope[i] / f->sd[j] - f->xbar[(size_t)l * f->k + j];
}

ST_retcode creghdfe_fe_recover_slopes(const HDFE_State *S, const creghdfe_input *in,
                                      ST_double *error, ST_double **alpha)
{
    const ST_int N = S->N, G = S->G;
    const ST_double *wt = S->weights;
    ST_retcode rc = 0;
    creghdfe_effect eff[10];
    ST_int E = 0;
    ST_double *proj = NULL, **update = NULL, *rsum = NULL, *xr = NULL, *coef = NULL;
    memset(eff, 0, sizeof(eff));
    if (G > 10) return 198;

    /* The effects in the order of their factors */
    for (ST_int g = 0; g < G; g++) {
        ST_int e;
        for (e = 0; e < E; e++) {
            ST_int h = eff[e].intercept >= 0 ? eff[e].intercept : eff[e].slopes[0];
            if (in->fe_absvar[h] == in->fe_absvar[g]) break;
        }
        if (e == E) {
            eff[E].intercept = -1;
            eff[E].levels = S->factors[g].levels;
            eff[E].L = S->factors[g].num_levels;
            E++;
        }
        if (in->slope_idx[g]) eff[e].slopes[eff[e].k++] = g;
        else eff[e].intercept = g;
    }

    ST_int Lmax = 1, kmax = 1;
    for (ST_int e = 0; e < E; e++) {
        if (eff[e].L > Lmax) Lmax = eff[e].L;
        if (eff[e].k > kmax) kmax = eff[e].k;
    }
    proj = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(*proj));
    update = (ST_double **)calloc((size_t)G, sizeof(*update));
    rsum = (ST_double *)malloc((size_t)Lmax * sizeof(*rsum));
    xr = (ST_double *)ctools_safe_malloc2((size_t)Lmax, (size_t)kmax * sizeof(*xr));
    coef = (ST_double *)ctools_safe_malloc2((size_t)Lmax, (size_t)kmax * sizeof(*coef));
    if (!proj || !update || !rsum || !xr || !coef) { rc = 920; goto cleanup; }
    for (ST_int g = 0; g < G; g++) {
        update[g] = (ST_double *)calloc((size_t)S->factors[g].num_levels, sizeof(**update));
        if (!update[g]) { rc = 920; goto cleanup; }
    }

    /* Per level: weights, slope means and the inverse cross products */
    for (ST_int e = 0; e < E && !rc; e++) {
        creghdfe_effect *f = &eff[e];
        const ST_int L = f->L, k = f->k;
        f->wsum = (ST_double *)calloc((size_t)L, sizeof(*f->wsum));
        if (!f->wsum) { rc = 920; break; }
        for (ST_int i = 0; i < N; i++) f->wsum[f->levels[i] - 1] += wt ? wt[i] : 1.0;
        if (!k) continue;
        f->xbar = (ST_double *)calloc((size_t)L * k, sizeof(*f->xbar));
        f->inv = (ST_double *)calloc((size_t)L * k * k, sizeof(*f->inv));
        ST_double *cross = (ST_double *)calloc((size_t)L * k * k, sizeof(*cross));
        if (!f->xbar || !f->inv || !cross) { free(cross); rc = 920; break; }
        for (ST_int j = 0; j < k; j++) f->sd[j] = slope_stdev(S->factors[f->slopes[j]].slope, N);
        if (f->intercept >= 0) {
            for (ST_int i = 0; i < N; i++) {
                ST_int l = f->levels[i] - 1;
                for (ST_int j = 0; j < k; j++)
                    f->xbar[(size_t)l * k + j] += (wt ? wt[i] : 1.0) *
                        S->factors[f->slopes[j]].slope[i] / f->sd[j];
            }
            for (ST_int l = 0; l < L; l++)
                for (ST_int j = 0; j < k; j++)
                    f->xbar[(size_t)l * k + j] /= f->wsum[l];
        }
        for (ST_int i = 0; i < N; i++) {
            ST_int l = f->levels[i] - 1;
            ST_double wi = wt ? wt[i] : 1.0, d[10];
            for (ST_int j = 0; j < k; j++) d[j] = effect_x(S, f, j, i, l);
            for (ST_int a = 0; a < k; a++)
                for (ST_int b = 0; b < k; b++)
                    cross[((size_t)l * k + a) * k + b] += wi * d[a] * d[b];
        }
        for (ST_int l = 0; l < L && !rc; l++) {
            ST_int rank = 0;
            if (ctools_invsym(cross + (size_t)l * k * k, k, f->inv + (size_t)l * k * k, &rank))
                rc = 920;
        }
        free(cross);
    }
    if (rc) goto cleanup;

    ST_double scale = 1.0, sum_w = 0.0;
    for (ST_int i = 0; i < N; i++) {
        if (fabs(error[i]) > scale) scale = fabs(error[i]);
        sum_w += wt ? wt[i] : 1.0;
    }
    ST_int converged = 0;
    for (ST_int iteration = 0; iteration < S->maxiter; iteration++) {
        /* One forward sweep over the effects, as hdfe_fe_recover(), each
         * effect fitted on all its coefficients within each level */
        memcpy(proj, error, (size_t)N * sizeof(*proj));
        for (ST_int e = 0; e < E; e++) {
            creghdfe_effect *f = &eff[e];
            const ST_int L = f->L, k = f->k;
            memset(rsum, 0, (size_t)L * sizeof(*rsum));
            if (k) memset(xr, 0, (size_t)L * k * sizeof(*xr));
            for (ST_int i = 0; i < N; i++) {
                ST_int l = f->levels[i] - 1;
                ST_double wr = (wt ? wt[i] : 1.0) * proj[i];
                rsum[l] += wr;
                for (ST_int j = 0; j < k; j++) xr[(size_t)l * k + j] += wr * effect_x(S, f, j, i, l);
            }
            for (ST_int l = 0; l < L; l++) {
                ST_double mean = f->intercept >= 0 && f->wsum[l] > 0 ? rsum[l] / f->wsum[l] : 0.0;
                ST_double a = mean;
                for (ST_int j = 0; j < k; j++) {
                    ST_double b = 0.0;
                    for (ST_int m = 0; m < k; m++)
                        b += f->inv[((size_t)l * k + j) * k + m] * xr[(size_t)l * k + m];
                    coef[(size_t)l * k + j] = b;
                    a -= (f->xbar[(size_t)l * k + j]) * b;
                    update[f->slopes[j]][l] = b / f->sd[j];
                }
                if (f->intercept >= 0) update[f->intercept][l] = a;
                rsum[l] = mean;
            }
            for (ST_int i = 0; i < N; i++) {
                ST_int l = f->levels[i] - 1;
                ST_double fitted = rsum[l];
                for (ST_int j = 0; j < k; j++) fitted += coef[(size_t)l * k + j] * effect_x(S, f, j, i, l);
                proj[i] -= fitted;
            }
        }
        ST_double num = 0.0, den = 0.0, change = 0.0;
        for (ST_int i = 0; i < N; i++) {
            ST_double weight = wt ? wt[i] : 1.0;
            ST_double step = error[i] - proj[i];
            proj[i] = step;
            num += weight * error[i] * step;
            den += weight * step * step;
            change += weight * fabs(step) / (fabs(error[i]) + scale);
        }
        if (change / sum_w <= S->tolerance) {
            converged = 1;
            break;
        }
        ST_double t = den > 0.0 ? num / den : 1.0;
        if (!isfinite(t)) t = 1.0;
        for (ST_int i = 0; i < N; i++) error[i] -= t * proj[i];
        for (ST_int g = 0; g < G; g++)
            for (ST_int l = 0; l < S->factors[g].num_levels; l++)
                alpha[g][l] += t * update[g][l];
    }
    if (!converged) rc = 430;
cleanup:
    for (ST_int e = 0; e < E; e++) {
        free(eff[e].wsum); free(eff[e].xbar); free(eff[e].inv);
    }
    for (ST_int g = 0; update && g < G; g++) free(update[g]);
    free(update); free(proj); free(rsum); free(xr); free(coef);
    return rc;
}

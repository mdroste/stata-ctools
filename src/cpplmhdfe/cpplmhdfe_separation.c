/*
 * cpplmhdfe_separation.c
 *
 * Separation detection for PPML estimation.
 * Implements FE-level screening and mu-based post-iteration screening.
 * Part of the ctools Stata plugin suite
 */

#include "cpplmhdfe_separation.h"
#include <string.h>

/*
 * Detect separation via FE-level screening.
 * For each FE level, if ALL observations in that level have y = 0,
 * mark them as separated (their FE coefficient would go to -infinity).
 */
ST_int ppml_detect_separation_fe(
    const ST_double *y,
    const FE_Factor *factors,
    ST_int G,
    ST_int N,
    ST_int *sep_mask)
{
    ST_int g, i, lev;
    ST_int num_separated = 0;

    memset(sep_mask, 0, N * sizeof(ST_int));

    for (g = 0; g < G; g++) {
        ST_int num_levels = factors[g].num_levels;

        /* Compute sum of y per level */
        ST_double *y_sum = (ST_double *)calloc(num_levels, sizeof(ST_double));
        if (!y_sum) continue;

        for (i = 0; i < N; i++) {
            lev = factors[g].levels[i] - 1;  /* 1-based to 0-based */
            if (lev >= 0 && lev < num_levels) {
                y_sum[lev] += y[i];
            }
        }

        /* Mark observations in all-zero levels */
        for (i = 0; i < N; i++) {
            if (sep_mask[i]) continue;  /* Already marked */
            lev = factors[g].levels[i] - 1;
            if (lev >= 0 && lev < num_levels && y_sum[lev] == 0.0) {
                sep_mask[i] = 1;
                num_separated++;
            }
        }

        free(y_sum);
    }

    return num_separated;
}

/*
 * Detect separation via mu-based screening (post-iteration).
 * For observations where y=0 and mu < sep_tol, mark as separated.
 */
ST_int ppml_detect_separation_mu(
    const ST_double *y,
    const ST_double *mu,
    ST_double sep_tol,
    ST_int N,
    ST_int *sep_mask)
{
    ST_int i;
    ST_int num_new = 0;

    for (i = 0; i < N; i++) {
        if (!sep_mask[i] && y[i] == 0.0 && mu[i] < sep_tol) {
            sep_mask[i] = 1;
            num_new++;
        }
    }

    return num_new;
}

ST_int ppml_select_fe_sample(const ST_double *y, ST_int *const *levels,
    const ST_int *num_levels, ST_int G, ST_int N, ST_int *mask,
    ST_int *num_singletons, ST_int *num_separated)
{
    ST_int **counts = calloc((size_t)G, sizeof(*counts));
    unsigned char **positive = calloc((size_t)G, sizeof(*positive));
    ST_int rc = -1;
    *num_singletons = *num_separated = 0;
    if (!counts || !positive) goto cleanup;
    for (ST_int g = 0; g < G; g++) {
        if (num_levels[g] < 1) goto cleanup;
        counts[g] = calloc((size_t)num_levels[g], sizeof(**counts));
        positive[g] = calloc((size_t)num_levels[g], sizeof(**positive));
        if (!counts[g] || !positive[g]) goto cleanup;
        for (ST_int i = 0; i < N; i++)
            if (levels[g][i] < 1 || levels[g][i] > num_levels[g]) goto cleanup;
    }
    for (;;) {
        for (ST_int g = 0; g < G; g++) {
            memset(counts[g], 0, (size_t)num_levels[g] * sizeof(**counts));
            memset(positive[g], 0, (size_t)num_levels[g] * sizeof(**positive));
            for (ST_int i = 0; i < N; i++) {
                if (!mask[i]) continue;
                ST_int level = levels[g][i] - 1;
                counts[g][level]++;
                if (y[i] > 0) positive[g][level] = 1;
            }
        }
        ST_int removed = 0;
        for (ST_int i = 0; i < N; i++) {
            if (!mask[i]) continue;
            ST_int singleton = 0, separated = 0;
            for (ST_int g = 0; g < G; g++) {
                ST_int level = levels[g][i] - 1;
                singleton |= counts[g][level] == 1;
                separated |= !positive[g][level];
            }
            if (singleton || separated) {
                mask[i] = 0;
                removed++;
                if (singleton) (*num_singletons)++;
                else (*num_separated)++;
            }
        }
        if (!removed) break;
    }
    rc = 0;
cleanup:
    for (ST_int g = 0; g < G; g++) {
        if (counts) free(counts[g]);
        if (positive) free(positive[g]);
    }
    free(counts); free(positive);
    return rc;
}

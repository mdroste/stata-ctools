/*
 * cpplmhdfe_separation.c
 *
 * FE-based sample selection for PPML estimation: iteratively drops
 * singletons and FE groups whose outcome is all zero until stable.
 * Part of the ctools Stata plugin suite
 */

#include "cpplmhdfe_separation.h"
#include <string.h>

ST_int ppml_select_fe_sample(const ST_double *y, ST_int *const *levels,
    const ST_int *num_levels, ST_int G, ST_int N, ST_int *mask,
    ST_int *num_singletons, ST_int *num_separated)
{
    return ppml_select_fe_sample_ex(y, levels, num_levels, G, N, mask,
        num_singletons, num_separated, NULL, 0, 1);
}

ST_int ppml_select_fe_sample_ex(const ST_double *y, ST_int *const *levels,
    const ST_int *num_levels, ST_int G, ST_int N, ST_int *mask,
    ST_int *num_singletons, ST_int *num_separated,
    const ST_int *slope_only, ST_int keep_singletons, ST_int check_separation)
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
                singleton |= !keep_singletons && counts[g][level] == 1;
                separated |= check_separation && !(slope_only && slope_only[g]) && !positive[g][level];
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

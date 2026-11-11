/*
 * hdfe_fe_recover.c
 *
 * Fixed-effect coefficients from the fixed-effect component (see
 * hdfe_fe_recover.h). Moved unchanged from creghdfe's store_fixed_effects so
 * that bootstrap replicates can run it without the SPI.
 * Part of the ctools Stata plugin suite
 */

#include "hdfe_fe_recover.h"
#include "../ctools_cancel.h"

ST_retcode hdfe_fe_recover(const HDFE_State *S, ST_double *error, ST_double **alpha)
{
    ST_int N = S->N, G = S->G;
    ST_retcode rc = 0;
    ST_double *proj = malloc((size_t)N * sizeof(*proj));
    ST_double **update = calloc((size_t)G, sizeof(*update));
    ST_double scale = 1.0;
    if (!proj || !update)
    {
        rc = 920;
        goto cleanup;
    }
    for (ST_int g = 0; g < G; g++)
    {
        update[g] = calloc((size_t)S->factors[g].num_levels, sizeof(**update));
        if (!update[g])
        {
            rc = 920;
            goto cleanup;
        }
    }
    for (ST_int i = 0; i < N; i++)
    {
        if (fabs(error[i]) > scale)
            scale = fabs(error[i]);
    }
    ST_double sum_w = 0.0;
    for (ST_int i = 0; i < N; i++)
        sum_w += S->weights ? S->weights[i] : 1.0;
    ST_int converged = 0;
    for (ST_int iteration = 0; iteration < S->maxiter; iteration++)
    {
        if (ctools_cancel_requested()) { rc = 1; goto cleanup; }
        /* One forward sweep on the current error: proj = error - the error
         * left after projecting out each effect in turn */
        memcpy(proj, error, (size_t)N * sizeof(*proj));
        for (ST_int g = 0; g < G; g++)
        {
            FE_Factor *f = &S->factors[g];
            memset(update[g], 0, (size_t)f->num_levels * sizeof(**update));
            for (ST_int i = 0; i < N; i++)
            {
                ST_double weight = S->weights ? S->weights[i] : 1.0;
                update[g][f->levels[i] - 1] += weight * proj[i];
            }
            for (ST_int level = 0; level < f->num_levels; level++)
            {
                ST_double count = S->weights ? f->weighted_counts[level] : f->counts[level];
                update[g][level] = count > 0 ? update[g][level] / count : 0;
            }
            for (ST_int i = 0; i < N; i++)
                proj[i] -= update[g][f->levels[i] - 1];
        }
        ST_double num = 0.0, den = 0.0, change = 0.0;
        for (ST_int i = 0; i < N; i++)
        {
            ST_double weight = S->weights ? S->weights[i] : 1.0;
            ST_double step = error[i] - proj[i];
            proj[i] = step;
            num += weight * error[i] * step;
            den += weight * step * step;
            change += weight * fabs(step) / (fabs(error[i]) + scale);
        }
        if (change / sum_w <= S->tolerance)
        {
            converged = 1;
            break;
        }
        /* Steepest-descent step length along the sweep */
        ST_double t = den > 0.0 ? num / den : 1.0;
        if (!isfinite(t))
            t = 1.0;
        for (ST_int i = 0; i < N; i++)
            error[i] -= t * proj[i];
        for (ST_int g = 0; g < G; g++)
        {
            for (ST_int level = 0; level < S->factors[g].num_levels; level++)
                alpha[g][level] += t * update[g][level];
        }
    }
    if (!converged)
        rc = 430;
cleanup:
    for (ST_int g = 0; update && g < G; g++)
        free(update[g]);
    free(update);
    free(proj);
    return rc;
}

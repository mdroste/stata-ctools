/*
 * hdfe_sample.c
 *
 * Estimation-sample construction shared by the HDFE estimators
 * Part of the ctools Stata plugin suite
 */

#include "hdfe_sample.h"
#include "hdfe_remap.h"

ST_retcode hdfe_load_factors(HDFE_LoadedFactor *factors, ST_int G, const ST_double *const *values,
                             ST_int N)
{
    int failed = 0;
    ST_int g;

    for (g = 0; g < G; g++)
    {
        factors[g].levels = (ST_int *)malloc((size_t)(N > 0 ? N : 1) * sizeof(ST_int));
        factors[g].num_obs = N;
        factors[g].counts = NULL;
        factors[g].num_levels = 0;
        if (!factors[g].levels)
            failed = 1;
    }
    if (failed)
        return 920;

    /* Fused remap + counting: O(N) for integer codes, a sort otherwise */
#pragma omp parallel for schedule(dynamic) reduction(| : failed) if (G > 1)
    for (g = 0; g < G; g++)
    {
        ST_int *counts_g = NULL;
        if (remap_and_count(values[g], N, factors[g].levels, &factors[g].num_levels, &counts_g,
                            NULL, NULL) != 0 ||
            factors[g].num_levels == 0)
        {
            failed |= 1;
        }
        factors[g].counts = counts_g;
    }
    return failed ? 920 : 0;
}

ST_retcode hdfe_build_factors(HDFE_State *S, ST_int G, ST_int *const *levels,
                              const ST_int *level_max, ST_int N_in, const ST_int *keep, ST_int N,
                              const ST_double *weights, ST_int weighted)
{
    ST_retcode rc = 0;
    ST_int *remap = NULL;

    if (!S || S->factors || G < 0 || N_in < 0 || N < 0 || N > N_in || (G > 0 && (!levels || !level_max)))
    {
        return 198;
    }

    S->G = G;
    S->N = N;
    S->factors = (FE_Factor *)calloc((size_t)(G > 0 ? G : 1), sizeof(FE_Factor));
    if (!S->factors)
    {
        return 920;
    }

    for (ST_int g = 0; g < G && !rc; g++)
    {
        FE_Factor *f = &S->factors[g];
        const ST_int *src = levels[g];
        ST_int next_level = 0, idx = 0;

        if (level_max[g] < 0)
        {
            rc = 198;
            break;
        }
        remap = (ST_int *)calloc((size_t)level_max[g] + 1, sizeof(ST_int));
        f->levels = (ST_int *)malloc((size_t)(N > 0 ? N : 1) * sizeof(ST_int));
        if (!remap || !f->levels)
        {
            rc = 920;
            break;
        }

        /* Renumber the surviving levels in order of first appearance */
        for (ST_int i = 0; i < N_in; i++)
        {
            if (keep && !keep[i])
                continue;
            ST_int old_level = src[i];
            if (old_level < 1 || old_level > level_max[g] || idx >= N)
            {
                rc = 198;
                break;
            }
            if (remap[old_level] == 0)
                remap[old_level] = ++next_level;
            f->levels[idx++] = remap[old_level];
        }
        free(remap);
        remap = NULL;
        if (!rc && idx != N)
            rc = 198;
        if (rc)
            break;

        f->num_levels = next_level;
        f->max_level = next_level;
        f->has_intercept = 1;
        f->counts = (ST_double *)calloc((size_t)(next_level > 0 ? next_level : 1), sizeof(ST_double));
        if (weights || weighted)
        {
            f->weighted_counts =
                (ST_double *)calloc((size_t)(next_level > 0 ? next_level : 1), sizeof(ST_double));
        }
        if (!f->counts || ((weights || weighted) && !f->weighted_counts))
        {
            rc = 920;
            break;
        }

        /* Counts, and weights accumulated in row order */
        idx = 0;
        for (ST_int i = 0; i < N_in; i++)
        {
            if (keep && !keep[i])
                continue;
            ST_int level = f->levels[idx++] - 1;
            f->counts[level] += 1.0;
            if (weights)
                f->weighted_counts[level] += weights[i];
        }
    }

    if (rc)
    {
        free(remap);
        ctools_hdfe_state_cleanup(S);
    }
    return rc;
}

void hdfe_compact_columns(const ST_double *src, ST_int N_in, ST_int K, const ST_int *keep,
                          ST_double *dst, ST_int N)
{
    if (!keep)
    {
        memcpy(dst, src, (size_t)N * (size_t)K * sizeof(ST_double));
        return;
    }
    for (ST_int k = 0; k < K; k++)
    {
        const ST_double *src_col = src + (size_t)k * N_in;
        ST_double *dst_col = dst + (size_t)k * N;
        ST_int idx = 0;
        for (ST_int i = 0; i < N_in && idx < N; i++)
        {
            if (keep[i])
                dst_col[idx++] = src_col[i];
        }
    }
}

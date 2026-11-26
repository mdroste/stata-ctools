/*
 * hdfe_sample.c
 *
 * Estimation-sample construction shared by the HDFE estimators
 * Part of the ctools Stata plugin suite
 */

#include "hdfe_sample.h"
#include "hdfe_remap.h"
#include <stdint.h>

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

    /* Each remap runs in parallel over its rows; with too few rows for that,
     * the factors run in parallel instead (the remaps then run serially).
     * The remaps are not called from an inactive parallel region: their own
     * regions would then be nested, which costs a new thread team each. */
    if (G > 1 && ctools_hdfe_row_blocks(N) == 1)
    {
#pragma omp parallel for schedule(dynamic) reduction(| : failed)
        for (g = 0; g < G; g++)
        {
            if (hdfe_remap_levels(values[g], N, factors[g].levels, &factors[g].num_levels) != 0 ||
                factors[g].num_levels == 0)
            {
                failed |= 1;
            }
        }
    }
    else
    {
        for (g = 0; g < G; g++)
        {
            if (hdfe_remap_levels(values[g], N, factors[g].levels, &factors[g].num_levels) != 0 ||
                factors[g].num_levels == 0)
            {
                failed |= 1;
            }
        }
    }
    return failed ? 920 : 0;
}

/* Largest per-block bitset workspace of build_factors_parallel */
#define HDFE_BUILD_BITSET_BUDGET ((size_t)256 * 1024 * 1024)

/*
 * hdfe_build_factors() in B row blocks, with the same result: each factor's
 * surviving levels are renumbered in order of first appearance among the
 * kept rows, as the serial loop does it, without atomic operations:
 *   1. each block marks the old levels of its kept rows in a bitset of its
 *      own, and counts its kept rows
 *   2. a level stays marked only in the first block that has it; each
 *      block's marks are counted, which gives every block the number its
 *      first new level takes
 *   3. each block numbers its marked levels in row order (a level is
 *      numbered at its first kept row in that block, and unmarked)
 *   4. the kept rows read their new levels
 * Counts, and weighted counts summed in row order, take one serial pass per
 * factor, the factors in parallel. Returns 1 (S untouched) when the bitsets
 * would exceed their budget, so the caller runs the serial loop.
 */
static ST_retcode build_factors_parallel(HDFE_State *S, ST_int G, ST_int *const *levels,
                                         const ST_int *level_max, ST_int N_in, const ST_int *keep,
                                         ST_int N, const ST_double *weights, ST_int weighted, int B)
{
    ST_retcode rc = 0;
    const size_t stride = (size_t)B + 1;
    size_t *woff = (size_t *)calloc((size_t)G + 1, sizeof(size_t));
    if (!woff)
        return 920;
    for (ST_int g = 0; g < G; g++)
    {
        if (level_max[g] < 0)
        {
            free(woff);
            return 198;
        }
        woff[g + 1] = woff[g] + ((size_t)level_max[g] + 64) / 64;
    }
    const size_t words = woff[G];
    if (words > HDFE_BUILD_BITSET_BUDGET / 8 / (size_t)B)
    {
        free(woff);
        return 1;
    }

    uint64_t *bits = (uint64_t *)calloc((size_t)B * words, sizeof(uint64_t));
    ST_int **newid = (ST_int **)calloc((size_t)G, sizeof(ST_int *));
    ST_int *kept = (ST_int *)calloc(stride, sizeof(ST_int));
    ST_int *nfirst = (ST_int *)calloc((size_t)G * stride, sizeof(ST_int));
    int bad = 0, b;
    ST_int g;
    long long w;
    if (!bits || !newid || !kept || !nfirst)
    {
        rc = 920;
        goto done;
    }
    for (g = 0; g < G; g++)
    {
        newid[g] = (ST_int *)malloc(((size_t)level_max[g] + 1) * sizeof(ST_int));
        S->factors[g].levels = (ST_int *)malloc((size_t)(N > 0 ? N : 1) * sizeof(ST_int));
        if (!newid[g] || !S->factors[g].levels)
        {
            rc = 920;
            goto done;
        }
    }

    /* 1. Levels present in each block */
#pragma omp parallel for schedule(static) reduction(| : bad)
    for (b = 0; b < B; b++)
    {
        uint64_t *bb = bits + (size_t)b * words;
        ST_int i0, i1, n = 0;
        ctools_hdfe_block_range(N_in, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
        {
            if (keep && !keep[i])
                continue;
            n++;
            for (ST_int h = 0; h < G; h++)
            {
                ST_int old_level = levels[h][i];
                if (old_level < 1 || old_level > level_max[h])
                {
                    bad = 1;
                    continue;
                }
                bb[woff[h] + ((size_t)old_level >> 6)] |= (uint64_t)1 << (old_level & 63);
            }
        }
        kept[b + 1] = n;
    }
    for (b = 0; b < B; b++)
        kept[b + 1] += kept[b];
    if (bad || kept[B] != N)
    {
        rc = 198;
        goto done;
    }

    /* 2. Keep each level in its first block only, and count */
#pragma omp parallel for schedule(static)
    for (w = 0; w < (long long)words; w++)
    {
        uint64_t seen = 0;
        for (int t = 0; t < B; t++)
        {
            uint64_t v = bits[(size_t)t * words + (size_t)w];
            bits[(size_t)t * words + (size_t)w] = v & ~seen;
            seen |= v;
        }
    }
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        const uint64_t *bb = bits + (size_t)b * words;
        for (ST_int h = 0; h < G; h++)
        {
            ST_int n = 0;
            for (size_t k = woff[h]; k < woff[h + 1]; k++)
                n += __builtin_popcountll(bb[k]);
            nfirst[(size_t)h * stride + b + 1] = n;
        }
    }
    for (g = 0; g < G; g++)
        for (b = 0; b < B; b++)
            nfirst[(size_t)g * stride + b + 1] += nfirst[(size_t)g * stride + b];

    /* 3. New levels, at their first kept rows */
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        uint64_t *bb = bits + (size_t)b * words;
        ST_int i0, i1;
        ctools_hdfe_block_range(N_in, B, b, &i0, &i1);
        for (ST_int h = 0; h < G; h++)
        {
            const ST_int *src = levels[h];
            ST_int *id = newid[h];
            uint64_t *bh = bb + woff[h];
            ST_int next = nfirst[(size_t)h * stride + b];
            for (ST_int i = i0; i < i1; i++)
            {
                if (keep && !keep[i])
                    continue;
                const ST_int old_level = src[i];
                const uint64_t bit = (uint64_t)1 << (old_level & 63);
                uint64_t *word = &bh[(size_t)old_level >> 6];
                if (*word & bit)
                {
                    *word &= ~bit;
                    id[old_level] = ++next;
                }
            }
        }
    }

    /* 4. Levels of the kept rows */
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        ST_int i0, i1, idx = kept[b];
        ctools_hdfe_block_range(N_in, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
        {
            if (keep && !keep[i])
                continue;
            for (ST_int h = 0; h < G; h++)
                S->factors[h].levels[idx] = newid[h][levels[h][i]];
            idx++;
        }
    }

    for (g = 0; g < G; g++)
    {
        FE_Factor *f = &S->factors[g];
        const ST_int L = nfirst[(size_t)g * stride + B];
        f->num_levels = L;
        f->max_level = L;
        f->has_intercept = 1;
        f->counts = (ST_double *)calloc((size_t)(L > 0 ? L : 1), sizeof(ST_double));
        if (weights || weighted)
            f->weighted_counts = (ST_double *)calloc((size_t)(L > 0 ? L : 1), sizeof(ST_double));
        if (!f->counts || ((weights || weighted) && !f->weighted_counts))
        {
            rc = 920;
            goto done;
        }
    }

    /* Counts, and weights accumulated in row order */
#pragma omp parallel for schedule(dynamic) if (G > 1)
    for (g = 0; g < G; g++)
    {
        FE_Factor *f = &S->factors[g];
        if (!weights)
        {
            for (ST_int i = 0; i < N; i++)
                f->counts[f->levels[i] - 1] += 1.0;
            continue;
        }
        ST_int idx = 0;
        for (ST_int i = 0; i < N_in; i++)
        {
            if (keep && !keep[i])
                continue;
            ST_int level = f->levels[idx++] - 1;
            f->counts[level] += 1.0;
            f->weighted_counts[level] += weights[i];
        }
    }

done:
    if (newid)
        for (g = 0; g < G; g++)
            free(newid[g]);
    free(newid);
    free(bits);
    free(kept);
    free(nfirst);
    free(woff);
    return rc;
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

    const int B = ctools_hdfe_row_blocks(N_in);
    if (B > 1 && G > 0)
    {
        rc = build_factors_parallel(S, G, levels, level_max, N_in, keep, N, weights, weighted, B);
        if (rc != 1)
        {
            if (rc)
                ctools_hdfe_state_cleanup(S);
            return rc;
        }
        rc = 0;
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

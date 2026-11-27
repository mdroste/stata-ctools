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

/* ========================================================================
 * Load, singleton removal and factors of the sample, sharing their work
 * ======================================================================== */

ST_retcode hdfe_load_factors_ordered(HDFE_LoadedFactor *factors, HDFE_FactorOrder *order, ST_int G,
                                     const ST_double *const *values, ST_int N)
{
    int failed = 0;
    ST_int g;

    for (g = 0; g < G; g++)
    {
        factors[g].levels = (ST_int *)malloc((size_t)(N > 0 ? N : 1) * sizeof(ST_int));
        factors[g].num_obs = N;
        factors[g].counts = NULL;
        factors[g].num_levels = 0;
        order[g].ordered = 0;
        order[g].rows.first_row = NULL;
        order[g].rows.single_rows = NULL;
        order[g].rows.num_single = -1;
        order[g].rows.local = 0;
        if (!factors[g].levels)
            failed = 1;
    }
    if (failed)
        return 920;

    /* As in hdfe_load_factors(): the remaps are parallel over their rows, or
     * over the factors when the rows are few */
    if (G > 1 && ctools_hdfe_row_blocks(N) == 1)
    {
#pragma omp parallel for schedule(dynamic) reduction(| : failed)
        for (g = 0; g < G; g++)
        {
            int rc = hdfe_remap_levels_ordered(values[g], N, factors[g].levels,
                                               &factors[g].num_levels, &order[g].rows);
            if (rc < 0 || factors[g].num_levels == 0)
                failed |= 1;
            else
                order[g].ordered = rc == 0;
        }
    }
    else
    {
        for (g = 0; g < G; g++)
        {
            int rc = hdfe_remap_levels_ordered(values[g], N, factors[g].levels,
                                               &factors[g].num_levels, &order[g].rows);
            if (rc < 0 || factors[g].num_levels == 0)
                failed |= 1;
            else
                order[g].ordered = rc == 0;
        }
    }
    return failed ? 920 : 0;
}

void hdfe_factor_order_free(HDFE_FactorOrder *order, ST_int G)
{
    if (!order)
        return;
    for (ST_int g = 0; g < G; g++)
    {
        free(order[g].rows.first_row);
        free(order[g].rows.single_rows);
        order[g].rows.first_row = NULL;
        order[g].rows.single_rows = NULL;
        order[g].rows.num_single = -1;
        order[g].rows.local = 0;
    }
}

void hdfe_singletons_free(HDFE_Singletons *sing)
{
    if (!sing)
        return;
    free(sing->keep);
    free(sing->dropped);
    if (sing->affected)
        for (ST_int g = 0; g < sing->G; g++)
        {
            free(sing->affected[g].level);
            free(sing->affected[g].left);
        }
    free(sing->affected);
    if (sing->rows_left)
        for (ST_int g = 0; g < sing->G; g++)
            free(sing->rows_left[g]);
    free(sing->rows_left);
    memset(sing, 0, sizeof(*sing));
}

static inline int level_bit(const uint64_t *bits, ST_int level)
{
    return (int)((bits[(size_t)level >> 6] >> (level & 63)) & 1);
}

static inline void level_bit_set(uint64_t *bits, ST_int level)
{
    bits[(size_t)level >> 6] |= (uint64_t)1 << (level & 63);
}

/* Entries of the ascending array a[0..n) that are below x */
static ST_int ordered_below(const ST_int *a, ST_int n, ST_int x)
{
    ST_int lo = 0, hi = n;
    while (lo < hi)
    {
        ST_int mid = lo + (hi - lo) / 2;
        if (a[mid] < x)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}

/* The levels of one factor that singleton removal follows: those of the
 * dropped rows. A counted level has the rows it has left and the XOR of
 * their numbers (left with one row, it names that row); a level with rows
 * to spare in the first rows of the data is not counted, and only has a
 * lower bound of the rows left until that falls below two */
typedef struct
{
    ST_int *slot;      /* slot[level] = index + 1 in the arrays below, 0 = not followed */
    uint64_t *counted; /* The levels of the counting pass under way */
    ST_int n, cap;
    ST_int *level, *count, *xr;
    ST_int *state; /* 0 lower bound, 1 counted, 2 to be counted in this round */
    ST_int *place; /* Position in the list of the counting pass */
} sing_track;

#define SING_BOUND 0
#define SING_EXACT 1
#define SING_QUEUED 2

static int sing_track_add(sing_track *t, ST_int level, ST_int count, ST_int xr, ST_int state)
{
    if (t->n == t->cap)
    {
        const ST_int cap = t->cap ? 2 * t->cap : 256;
        ST_int **arrays[5] = {&t->level, &t->count, &t->xr, &t->state, &t->place};
        for (int a = 0; a < 5; a++)
        {
            ST_int *grown = (ST_int *)realloc(*arrays[a], (size_t)cap * sizeof(ST_int));
            if (!grown)
                return -1;
            *arrays[a] = grown;
        }
        t->cap = cap;
    }
    t->level[t->n] = level;
    t->count[t->n] = count;
    t->xr[t->n] = xr;
    t->state[t->n] = state;
    t->place[t->n] = 0;
    t->slot[level] = ++t->n;
    return 0;
}

/* Rows and XOR of row numbers of the m followed levels list[], over all N
 * rows, in B blocks */
static int sing_count_levels(const ST_int *lev, ST_int N, int B, sing_track *t, const ST_int *list,
                             ST_int m)
{
    ST_int *priv = (ST_int *)calloc((size_t)B * (size_t)m * 2, sizeof(ST_int));
    int b;
    if (!priv)
        return -1;
    for (ST_int k = 0; k < m; k++)
    {
        const ST_int level = t->level[list[k]];
        t->place[list[k]] = k;
        t->counted[(size_t)level >> 6] |= (uint64_t)1 << (level & 63);
    }
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        ST_int *p = priv + (size_t)b * (size_t)m * 2;
        const uint64_t *counted = t->counted;
        const ST_int *slot = t->slot, *place = t->place;
        ST_int i0, i1;
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
        {
            const ST_int l = lev[i];
            if (!((counted[(size_t)l >> 6] >> (l & 63)) & 1))
                continue;
            const size_t s = (size_t)place[slot[l] - 1];
            p[2 * s]++;
            p[2 * s + 1] ^= i;
        }
    }
    for (ST_int k = 0; k < m; k++)
    {
        const ST_int level = t->level[list[k]];
        ST_int n = 0, x = 0;
        for (b = 0; b < B; b++)
        {
            n += priv[((size_t)b * (size_t)m + (size_t)k) * 2];
            x ^= priv[((size_t)b * (size_t)m + (size_t)k) * 2 + 1];
        }
        t->count[list[k]] = n;
        t->xr[list[k]] = x;
        t->counted[(size_t)level >> 6] = 0;
    }
    free(priv);
    return 0;
}

static int sing_compare_rows(const void *a, const void *b)
{
    const ST_int x = *(const ST_int *)a, y = *(const ST_int *)b;
    return (x > y) - (x < y);
}

/* Appends a row to the dropped rows unless it is there already or, with
 * fweights, its weight is not 1. Returns -1 on allocation failure. */
static int sing_drop(ST_int row, uint64_t *gone, const ST_double *fweights, ST_int **rows,
                     size_t *n_rows, size_t *cap_rows)
{
    if (((gone[(size_t)row >> 6] >> (row & 63)) & 1) || (fweights && fweights[row] != 1.0))
        return 0;
    if (*n_rows == *cap_rows)
    {
        size_t cap = *cap_rows ? 2 * *cap_rows : 1024;
        ST_int *grown = (ST_int *)realloc(*rows, cap * sizeof(ST_int));
        if (!grown)
            return -1;
        *rows = grown;
        *cap_rows = cap;
    }
    gone[(size_t)row >> 6] |= (uint64_t)1 << (row & 63);
    (*rows)[(*n_rows)++] = row;
    return 0;
}

/* Most private count entries (levels x blocks) of one counting pass, and the
 * first rows looked at before counting */
#define HDFE_SINGLETON_PASS_ENTRIES ((size_t)1 << 22)
/* The removal lists what it drops while the rows alone in a level are at
 * most one row in 128, or this many */
#ifndef HDFE_LISTED_MIN_SEEDS /* tests set a small one */
#define HDFE_LISTED_MIN_SEEDS 4096
#endif
/* Dropped rows are listed (and the levels shifted over them) while they are
 * at most one row in 64, or this many */
#ifndef HDFE_LISTED_MIN_DROPPED /* tests set a small one */
#define HDFE_LISTED_MIN_DROPPED 8192
#endif
#ifndef HDFE_SINGLETON_PROBE_ROWS /* tests set a small one */
#define HDFE_SINGLETON_PROBE_ROWS 32768
#endif
/* Levels from which the new numbers are worked out in parallel ranges */
#ifndef HDFE_ORDER_RANGE_LEVELS
#define HDFE_ORDER_RANGE_LEVELS 262144
#endif

/*
 * Singleton removal from the rows alone in their level: only the levels of
 * dropped rows are followed. Returns 0 when done (out filled in), 1 when the
 * removals are too many or spread too widely for this to pay, or an allocation
 * failed: the caller then starts again with singletons_seeded().
 */
static int singletons_listed(HDFE_Singletons *out, ST_int *const *levels, const ST_int *num_levels,
                             const HDFE_FactorOrder *order, ST_int G, ST_int N,
                             const ST_double *fweights, int B)
{
    int result = 1;
    ST_int g, passes = 0;
    size_t seeds = 0, most_seeds = (size_t)N / 128, most_dropped = (size_t)N / 64;
    if (most_seeds < HDFE_LISTED_MIN_SEEDS)
        most_seeds = HDFE_LISTED_MIN_SEEDS;
    if (most_dropped < HDFE_LISTED_MIN_DROPPED)
        most_dropped = HDFE_LISTED_MIN_DROPPED;
    for (g = 0; g < G; g++)
        seeds += (size_t)order[g].rows.num_single;
    if (seeds > most_seeds)
        return 1;

    sing_track *track = (sing_track *)calloc((size_t)G, sizeof(sing_track));
    uint64_t *gone = (uint64_t *)calloc((size_t)N / 64 + 1, sizeof(uint64_t)); /* dropped rows */
    ST_int *rows = NULL;  /* the dropped rows, in the order they are found */
    ST_int *queue = NULL; /* followed levels to count in this round */
    size_t n_rows = 0, cap_rows = 0, n_queue = 0, cap_queue = 0;
    const ST_int probe = N < HDFE_SINGLETON_PROBE_ROWS ? N : HDFE_SINGLETON_PROBE_ROWS;
    if (!track || !gone)
        goto done;
    for (g = 0; g < G; g++)
    {
        track[g].slot = (ST_int *)calloc((size_t)num_levels[g] + 1, sizeof(ST_int));
        track[g].counted = (uint64_t *)calloc((size_t)num_levels[g] / 64 + 1, sizeof(uint64_t));
        if (!track[g].slot || !track[g].counted)
            goto done;
    }

    /* The rows alone in a level: their level has that one row */
    for (g = 0; g < G; g++)
    {
        const ST_int *single = order[g].rows.single_rows;
        for (ST_int k = 0; k < order[g].rows.num_single; k++)
            if (sing_track_add(&track[g], levels[g][single[k]], 1, single[k], SING_EXACT))
                goto done;
    }
    for (g = 0; g < G; g++)
    {
        const ST_int *single = order[g].rows.single_rows;
        for (ST_int k = 0; k < order[g].rows.num_single; k++)
            if (sing_drop(single[k], gone, fweights, &rows, &n_rows, &cap_rows))
                goto done;
    }

    /* Rounds: rows[f0..f1) were dropped in the last one. Their levels lose
     * them; a level left with one row drops that row in the next round */
    for (size_t f0 = 0, f1 = n_rows; f0 < f1; f0 = f1, f1 = n_rows)
    {
        if (n_rows > most_dropped)
            goto done;
        for (g = 0; g < G; g++)
        {
            sing_track *t = &track[g];
            const ST_int *lev = levels[g];
            const ST_int t0 = t->n;
            n_queue = 0;
#define SING_QUEUE(index) \
    do \
    { \
        if (n_queue == cap_queue) \
        { \
            size_t cap_ = cap_queue ? 2 * cap_queue : 256; \
            ST_int *grown_ = (ST_int *)realloc(queue, cap_ * sizeof(ST_int)); \
            if (!grown_) \
                goto done; \
            queue = grown_; \
            cap_queue = cap_; \
        } \
        t->state[index] = SING_QUEUED; \
        queue[n_queue++] = (index); \
    } while (0)

            for (size_t k = f0; k < f1; k++)
                if (!t->slot[lev[rows[k]]] && sing_track_add(t, lev[rows[k]], 0, 0, SING_BOUND))
                    goto done;
            if (t->n > t0)
            {
                /* Levels met for the first time: the rows they have among
                 * the first rows of the data, less all their dropped rows,
                 * are a lower bound of the rows they have left. With two or
                 * more they need no count yet */
                for (ST_int i = 0; i < probe; i++)
                {
                    const ST_int s = t->slot[lev[i]] - 1;
                    if (s >= t0)
                        t->count[s]++;
                }
                for (size_t k = 0; k < f1; k++)
                {
                    const ST_int s = t->slot[lev[rows[k]]] - 1;
                    if (s >= t0)
                        t->count[s]--;
                }
                for (ST_int s = t0; s < t->n; s++)
                    if (t->count[s] < 2)
                        SING_QUEUE(s);
            }
            /* The levels followed before lose the rows of this round */
            for (size_t k = f0; k < f1; k++)
            {
                const ST_int s = t->slot[lev[rows[k]]] - 1;
                if (s >= t0)
                    continue;
                t->count[s]--;
                t->xr[s] ^= rows[k];
                if (t->state[s] == SING_BOUND && t->count[s] < 2)
                    SING_QUEUE(s);
            }
#undef SING_QUEUE
            if (n_queue)
            {
                /* Count the queued levels over all rows, less the rows
                 * dropped up to this round */
                if (++passes > 2 * G + 2 || n_queue * (size_t)B > HDFE_SINGLETON_PASS_ENTRIES ||
                    sing_count_levels(lev, N, B, t, queue, (ST_int)n_queue))
                    goto done;
                for (size_t k = 0; k < f1; k++)
                {
                    const ST_int s = t->slot[lev[rows[k]]] - 1;
                    if (t->state[s] == SING_QUEUED)
                    {
                        t->count[s]--;
                        t->xr[s] ^= rows[k];
                    }
                }
                for (size_t k = 0; k < n_queue; k++)
                    t->state[queue[k]] = SING_EXACT;
            }
            /* A row dropped here by an earlier factor of this round is
             * taken out of its levels in the next one */
            for (size_t k = f0; k < f1; k++)
            {
                const ST_int s = t->slot[lev[rows[k]]] - 1;
                if (t->state[s] == SING_EXACT && t->count[s] == 1 &&
                    sing_drop(t->xr[s], gone, fweights, &rows, &n_rows, &cap_rows))
                    goto done;
            }
        }
    }

    out->affected = (HDFE_AffectedLevels *)calloc((size_t)(G > 0 ? G : 1), sizeof(HDFE_AffectedLevels));
    if (!out->affected)
        goto done;
    out->G = G;
    out->num_dropped = (ST_int)n_rows;
    if (n_rows)
    {
        ST_int i;
        out->keep = (ST_int *)malloc((size_t)N * sizeof(ST_int));
        if (!out->keep)
        {
            hdfe_singletons_free(out);
            goto done;
        }
#pragma omp parallel for schedule(static)
        for (i = 0; i < N; i++)
            out->keep[i] = 1;
        for (size_t k = 0; k < n_rows; k++)
            out->keep[rows[k]] = 0;
        qsort(rows, n_rows, sizeof(ST_int), sing_compare_rows);
        out->dropped = rows;
        rows = NULL;
        for (g = 0; g < G; g++)
        {
            out->affected[g].n = track[g].n;
            out->affected[g].level = track[g].level;
            out->affected[g].left = track[g].count;
            track[g].level = track[g].count = NULL;
        }
    }
    out->listed = 1;
    result = 0;

done:
    if (track)
        for (g = 0; g < G; g++)
        {
            free(track[g].slot);
            free(track[g].counted);
            free(track[g].level);
            free(track[g].count);
            free(track[g].xr);
            free(track[g].state);
            free(track[g].place);
        }
    free(track);
    free(gone);
    free(rows);
    free(queue);
    return result;
}

/* ------------------------------------------------------------------------
 * The same from any number of rows alone in their level
 * ------------------------------------------------------------------------ */

/* Every level of a factor is counted in row blocks, each with a table of its
 * own, while the tables hold this many entries at most */
#ifndef HDFE_SEEDED_TABLE_ENTRIES /* tests set a small one */
#define HDFE_SEEDED_TABLE_ENTRIES ((size_t)1 << 22)
#endif

/*
 * The counts of a factor in the removal below: cx[2 * l] and cx[2 * l + 1]
 * are the rows level l has left and the XOR of their numbers, so that a
 * level left with one row names it. Three ways to fill them from the rows
 * still in the sample (keep), each over all the threads; cx is zero on entry.
 */

/* Few levels: B row blocks with a table each. Returns 1 when the tables
 * cannot be allocated. */
static int seed_count_tables(const ST_int *lev, const ST_int *keep, ST_int N, ST_int L, int B,
                             ST_int *cx)
{
    const size_t entries = 2 * ((size_t)L + 1);
    /* tables at least four cache lines apart */
    const size_t stride = entries < 64 ? 64 : entries;
    ST_int *priv = (ST_int *)calloc((size_t)B * stride, sizeof(ST_int));
    int b;
    long long l;
    if (!priv)
        return 1;
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        ST_int *p = priv + (size_t)b * stride;
        ST_int i0, i1;
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
        {
            const ST_int in = keep[i] != 0;
            ST_int *slot = &p[2 * (size_t)lev[i]];
            slot[0] += in;
            slot[1] ^= i & -in;
        }
    }
#pragma omp parallel for schedule(static) if (L >= 65536)
    for (l = 0; l <= (long long)L; l++)
    {
        ST_int n = 0, x = 0;
        for (int t = 0; t < B; t++)
        {
            n += priv[(size_t)t * stride + 2 * (size_t)l];
            x ^= priv[(size_t)t * stride + 2 * (size_t)l + 1];
        }
        cx[2 * (size_t)l] = n;
        cx[2 * (size_t)l + 1] = x;
    }
    free(priv);
    return 0;
}

/* Rows in level order, the levels numbered in order of first appearance:
 * the levels first met in a block are that block's alone, which counts them
 * in place, and its rows of earlier levels are set aside and counted at the
 * end. Returns 1, with cx zero again, when a block has too many of the
 * latter or memory runs out. */
static int seed_count_local(const ST_int *lev, const ST_int *keep, ST_int N, ST_int L,
                            const ST_int *first_row, int B, ST_int *cx)
{
    ST_int **aside = (ST_int **)calloc((size_t)B, sizeof(ST_int *));
    /* each block's number of rows set aside, a cache line apart */
    ST_int *n_aside = (ST_int *)calloc((size_t)B * 16, sizeof(ST_int));
    int b, failed = 0;
    if (!aside || !n_aside)
    {
        free(aside);
        free(n_aside);
        return 1;
    }
#pragma omp parallel for schedule(static) reduction(| : failed)
    for (b = 0; b < B; b++)
    {
        ST_int i0, i1, n = 0;
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        /* the levels first met before this block */
        const ST_int earlier = ordered_below(first_row, L, i0);
        const ST_int cap = (i1 - i0) / 8 + 16;
        ST_int *list = (ST_int *)malloc((size_t)cap * sizeof(ST_int));
        aside[b] = list;
        if (!list)
        {
            failed = 1;
            continue;
        }
        for (ST_int i = i0; i < i1; i++)
        {
            const ST_int l = lev[i], in = keep[i] != 0;
            if (l > earlier)
            {
                ST_int *slot = &cx[2 * (size_t)l];
                slot[0] += in;
                slot[1] ^= i & -in;
            }
            else if (in)
            {
                if (n == cap)
                {
                    failed = 1;
                    break;
                }
                list[n++] = i;
            }
        }
        n_aside[(size_t)b * 16] = n;
    }
    if (failed)
    {
        long long k;
        const long long entries = 2 * ((long long)L + 1);
#pragma omp parallel for schedule(static)
        for (k = 0; k < entries; k++)
            cx[k] = 0;
    }
    else
    {
        for (b = 0; b < B; b++)
            for (ST_int k = 0; k < n_aside[(size_t)b * 16]; k++)
            {
                const ST_int i = aside[b][k];
                ST_int *slot = &cx[2 * (size_t)lev[i]];
                slot[0]++;
                slot[1] ^= i;
            }
    }
    for (b = 0; b < B; b++)
        free(aside[b]);
    free(aside);
    free(n_aside);
    return failed;
}

/* A row of a level, if it is still in the sample */
static inline void seed_count_row(ST_int *cx, ST_int level, ST_int in, ST_int row)
{
    if (in)
    {
        ST_int *slot = &cx[2 * (size_t)level];
        slot[0]++;
        slot[1] ^= row;
    }
}

/* Many levels in rows of no order: in T ranges of levels with a thread each
 * (CTOOLS_HDFE_ROWS_IN_RANGE). With the levels numbered in order of first
 * appearance (first_row given), a range has no row before the first row of
 * its lowest level. */
static void seed_count_ranges(const ST_int *lev, const ST_int *keep, ST_int N, ST_int L,
                              const ST_int *first_row, int T, ST_int *cx)
{
    int t;
#pragma omp parallel for schedule(static)
    for (t = 0; t < T; t++)
    {
        uint32_t lo, width;
        ctools_hdfe_level_range(L, T, t, &lo, &width);
        const ST_int start = (first_row && width > 0) ? first_row[lo - 1] : 0;
        CTOOLS_HDFE_ROWS_IN_RANGE(lev, N, lo, width, start, seed_count_row(cx, lev[row], keep[row], row));
    }
}

/*
 * Singleton removal from the rows alone in their level, which the load
 * lists, however many they are. Those rows go first. A factor in which every
 * row dropped so far was alone in its level has lost nothing else by it and
 * is left uncounted; once a row goes for another factor's sake, every level
 * of the factor is counted on the rows left, once, and from then on its
 * levels lose the rows of each round as they come. A level left with one row
 * drops that row in the next round; the rows a round drops stay in the
 * counts (flag 2) until that next round takes them out. The sample is the
 * fixed point, the same in whatever order the rows go. Returns 0 (out filled
 * in), or 1 when memory runs out: the caller then counts every level of
 * every factor.
 */
static int singletons_seeded(HDFE_Singletons *out, ST_int *const *levels, const ST_int *num_levels,
                             const HDFE_FactorOrder *order, ST_int G, ST_int N,
                             const ST_double *fweights, int B)
{
    int result = 1;
    ST_int g, i, seeds = 0;
    ST_int **cx = (ST_int **)calloc((size_t)G, sizeof(ST_int *));
    /* droppable rows alone in a level of each factor */
    ST_int *alone = (ST_int *)calloc((size_t)G, sizeof(ST_int));
    ST_int *keep = (ST_int *)malloc((size_t)N * sizeof(ST_int));
    ST_int *rows = NULL; /* the rows dropped after the first round, as they are found */
    size_t n_rows = 0, cap_rows = 0;
    if (!cx || !alone || !keep)
        goto done;

#pragma omp parallel for schedule(static)
    for (i = 0; i < N; i++)
        keep[i] = 1;

    /* The rows alone in a level (with fweights, those of weight 1). A
     * factor's are distinct rows; one alone in two factors is dropped once */
    for (g = 0; g < G; g++)
    {
        const ST_int *single = order[g].rows.single_rows;
        const long long n = order[g].rows.num_single;
        long long k;
        ST_int mine = 0, fresh = 0;
#pragma omp parallel for schedule(static) reduction(+ : mine, fresh) if (n >= 65536)
        for (k = 0; k < n; k++)
        {
            const ST_int row = single[k];
            if (fweights && fweights[row] != 1.0)
                continue;
            mine++;
            if (keep[row])
            {
                keep[row] = 0;
                fresh++;
            }
        }
        alone[g] = mine;
        seeds += fresh;
    }

/* Drops a row at the end of the round unless it is dropped already or, with
 * fweights, its weight is not 1 */
#define SEED_DROP(row) \
    do \
    { \
        const ST_int row_ = (row); \
        if (keep[row_] == 1 && (!fweights || fweights[row_] == 1.0)) \
        { \
            if (n_rows == cap_rows) \
            { \
                const size_t cap_ = cap_rows ? 2 * cap_rows : 1024; \
                ST_int *grown_ = (ST_int *)realloc(rows, cap_ * sizeof(ST_int)); \
                if (!grown_) \
                    goto done; \
                rows = grown_; \
                cap_rows = cap_; \
            } \
            keep[row_] = 2; \
            rows[n_rows++] = row_; \
        } \
    } while (0)

    /* Rounds: the first drops the rows above, each later one rows[f0..f1) */
    size_t f0 = 0, f1 = 0;
    for (int first = 1; first || f0 < f1; first = 0)
    {
        for (g = 0; g < G; g++)
        {
            const ST_int *lev = levels[g];
            const ST_int L = num_levels[g];
            if (cx[g])
            {
                /* counted: its levels lose this round's rows */
                ST_int *c = cx[g];
                for (size_t k = f0; k < f1; k++)
                {
                    const ST_int row = rows[k];
                    ST_int *slot = &c[2 * (size_t)lev[row]];
                    slot[0]--;
                    slot[1] ^= row;
                    if (slot[0] == 1)
                        SEED_DROP(slot[1]);
                }
                continue;
            }
            /* Not counted. In the first round the rows alone in a level of
             * this factor take nothing from its other levels; after it,
             * every row dropped comes from a counted factor */
            if (first && seeds == alone[g])
                continue;
            ST_int *c = (ST_int *)calloc(2 * ((size_t)L + 1), sizeof(ST_int));
            if (!c)
                goto done;
            cx[g] = c;
            /* on the rows left: this round's are out (flag 0) */
            const int tables = ((size_t)L + 1) * (size_t)B <= HDFE_SEEDED_TABLE_ENTRIES;
            const ST_int *first_row = order[g].ordered ? order[g].rows.first_row : NULL;
            if (!(tables && seed_count_tables(lev, keep, N, L, B, c) == 0) &&
                !(first_row && order[g].rows.local &&
                  seed_count_local(lev, keep, N, L, first_row, B, c) == 0))
                seed_count_ranges(lev, keep, N, L, first_row, B, c);
            for (ST_int l = 1; l <= L; l++)
                if (c[2 * (size_t)l] == 1)
                    SEED_DROP(c[2 * (size_t)l + 1]);
        }
        for (size_t k = f1; k < n_rows; k++)
            keep[rows[k]] = 0;
        f0 = f1;
        f1 = n_rows;
    }
#undef SEED_DROP

    out->G = G;
    out->num_dropped = seeds + (ST_int)n_rows;
    if (out->num_dropped == 0)
    {
        out->affected = (HDFE_AffectedLevels *)calloc((size_t)G, sizeof(HDFE_AffectedLevels));
        if (!out->affected)
        {
            memset(out, 0, sizeof(*out));
            goto done;
        }
        out->listed = 1;
    }
    else
    {
        out->counted = 1;
        out->rows_left = cx;
        cx = NULL;
        out->keep = keep;
        keep = NULL;
    }
    result = 0;

done:
    if (cx)
        for (g = 0; g < G; g++)
            free(cx[g]);
    free(cx);
    free(alone);
    free(keep);
    free(rows);
    return result;
}

ST_retcode hdfe_drop_singletons(HDFE_Singletons *out, ST_int *const *levels,
                                const ST_int *num_levels, const HDFE_FactorOrder *order, ST_int G,
                                ST_int N, const ST_double *fweights, ST_int verbose)
{
    if (!out || G < 0 || N < 0 || (G > 0 && (!levels || !num_levels)))
        return 198;
    memset(out, 0, sizeof(*out));
    if (G == 0 || N == 0)
    {
        out->listed = 1;
        out->affected = (HDFE_AffectedLevels *)calloc(1, sizeof(HDFE_AffectedLevels));
        return out->affected ? 0 : 920;
    }

    const int B = ctools_hdfe_row_blocks(N);
    int known = order != NULL && B > 1;
    for (ST_int g = 0; g < G && known; g++)
        known = order[g].rows.num_single >= 0 &&
                (order[g].rows.num_single == 0 || order[g].rows.single_rows);
    int done = known && singletons_listed(out, levels, num_levels, order, G, N, fweights, B) == 0;
    if (known && !done)
    {
        /* Too many rows alone in a level, or removals that spread, for the
         * list: the same start, with whole factors counted where need be */
        memset(out, 0, sizeof(*out));
        done = singletons_seeded(out, levels, num_levels, order, G, N, fweights, B) == 0;
    }
    if (done)
    {
        if (verbose && out->num_dropped > 0)
        {
            char buf[256];
            snprintf(buf, sizeof(buf), "ctools: Dropped %d singletons\n", (int)out->num_dropped);
            SF_display(buf);
        }
        return 0;
    }

    /* Every level of every factor counted */
    memset(out, 0, sizeof(*out));
    ST_int *keep = (ST_int *)malloc((size_t)N * sizeof(ST_int));
    if (!keep)
        return 920;
    ST_int dropped = ctools_remove_singletons((ST_int **)levels, G, N, keep, fweights, verbose);
    if (dropped < 0)
    {
        free(keep);
        return 920;
    }
    out->num_dropped = dropped;
    if (dropped == 0)
    {
        free(keep);
        out->affected = (HDFE_AffectedLevels *)calloc((size_t)G, sizeof(HDFE_AffectedLevels));
        if (!out->affected)
            return 920;
        out->G = G;
        out->listed = 1;
    }
    else
    {
        out->keep = keep;
    }
    return 0;
}

ST_retcode hdfe_singletons_from_keep(HDFE_Singletons *out, ST_int *const *levels,
                                     const ST_int *num_levels, ST_int G, ST_int N, ST_int *keep,
                                     ST_int num_dropped)
{
    if (!out || G < 0 || N < 0 || num_dropped < 0 || num_dropped > N ||
        (G > 0 && (!levels || !num_levels)) || (num_dropped > 0 && !keep))
        return 198;
    memset(out, 0, sizeof(*out));
    out->num_dropped = num_dropped;
    out->keep = keep;

    const int B = ctools_hdfe_row_blocks(N);
    size_t most_dropped = (size_t)N / 64;
    if (most_dropped < HDFE_LISTED_MIN_DROPPED)
        most_dropped = HDFE_LISTED_MIN_DROPPED;
    /* Many rows dropped, or few rows in all: the caller renumbers */
    if (num_dropped > 0 && (B == 1 || (size_t)num_dropped > most_dropped))
        return 0;

    ST_retcode rc = 920;
    int listed = 1, b;
    ST_int g;
    ST_int *rows = NULL, *start = NULL, *queue = NULL;
    sing_track *track = NULL;
    HDFE_AffectedLevels *affected =
        (HDFE_AffectedLevels *)calloc((size_t)(G > 0 ? G : 1), sizeof(HDFE_AffectedLevels));
    if (!affected)
        goto done;
    if (num_dropped > 0)
    {
        /* The dropped rows, ascending: each block's number, then the rows */
        rows = (ST_int *)malloc((size_t)num_dropped * sizeof(ST_int));
        start = (ST_int *)calloc((size_t)B + 1, sizeof(ST_int));
        track = (sing_track *)calloc((size_t)(G > 0 ? G : 1), sizeof(sing_track));
        queue = (ST_int *)malloc((size_t)num_dropped * sizeof(ST_int));
        if (!rows || !start || !track || !queue)
            goto done;
#pragma omp parallel for schedule(static)
        for (b = 0; b < B; b++)
        {
            ST_int i0, i1, n = 0;
            ctools_hdfe_block_range(N, B, b, &i0, &i1);
            for (ST_int i = i0; i < i1; i++)
                n += !keep[i];
            start[b + 1] = n;
        }
        for (b = 0; b < B; b++)
            start[b + 1] += start[b];
        if (start[B] != num_dropped)
        {
            rc = 198;
            goto done;
        }
#pragma omp parallel for schedule(static)
        for (b = 0; b < B; b++)
        {
            ST_int i0, i1, at = start[b];
            ctools_hdfe_block_range(N, B, b, &i0, &i1);
            for (ST_int i = i0; i < i1 && at < start[b + 1]; i++)
                if (!keep[i])
                    rows[at++] = i;
        }

        /* The levels of the dropped rows, and whether each has a row left: a
         * kept row among the first rows of the data settles it; the levels
         * without one are counted over all rows, less their dropped rows */
        const ST_int probe = N < HDFE_SINGLETON_PROBE_ROWS ? N : HDFE_SINGLETON_PROBE_ROWS;
        for (g = 0; g < G && listed; g++)
        {
            sing_track *t = &track[g];
            const ST_int *lev = levels[g];
            ST_int n_queue = 0;
            t->slot = (ST_int *)calloc((size_t)num_levels[g] + 1, sizeof(ST_int));
            t->counted = (uint64_t *)calloc((size_t)num_levels[g] / 64 + 1, sizeof(uint64_t));
            if (!t->slot || !t->counted)
                goto done;
            for (ST_int k = 0; k < num_dropped; k++)
            {
                const ST_int level = lev[rows[k]];
                if (level < 1 || level > num_levels[g])
                {
                    rc = 198;
                    goto done;
                }
                if (!t->slot[level] && sing_track_add(t, level, 0, 0, SING_BOUND))
                    goto done;
            }
            for (ST_int i = 0; i < probe; i++)
            {
                const ST_int s = t->slot[lev[i]] - 1;
                if (s >= 0 && keep[i])
                    t->count[s]++;
            }
            for (ST_int s = 0; s < t->n; s++)
                if (t->count[s] == 0)
                {
                    t->state[s] = SING_QUEUED;
                    queue[n_queue++] = s;
                }
            if (n_queue)
            {
                if ((size_t)n_queue * (size_t)B > HDFE_SINGLETON_PASS_ENTRIES)
                {
                    listed = 0; /* too many levels to follow: the caller renumbers */
                    break;
                }
                if (sing_count_levels(lev, N, B, t, queue, n_queue))
                    goto done;
                for (ST_int k = 0; k < num_dropped; k++)
                {
                    const ST_int s = t->slot[lev[rows[k]]] - 1;
                    if (t->state[s] == SING_QUEUED)
                        t->count[s]--;
                }
            }
        }
        if (listed)
        {
            for (g = 0; g < G; g++)
            {
                affected[g].n = track[g].n;
                affected[g].level = track[g].level;
                affected[g].left = track[g].count;
                track[g].level = track[g].count = NULL;
            }
            out->dropped = rows;
            rows = NULL;
        }
    }
    if (listed)
    {
        out->affected = affected;
        affected = NULL;
        out->G = G;
        out->listed = 1;
    }
    rc = 0;

done:
    if (track)
        for (g = 0; g < G; g++)
        {
            free(track[g].slot);
            free(track[g].counted);
            free(track[g].level);
            free(track[g].count);
            free(track[g].xr);
            free(track[g].state);
            free(track[g].place);
        }
    free(track);
    free(affected);
    free(rows);
    free(start);
    free(queue);
    if (rc)
        memset(out, 0, sizeof(*out));
    return rc;
}

void hdfe_singletons_counted(HDFE_Singletons *out, ST_int **rows_left, ST_int G)
{
    if (!rows_left)
        return;
    if (!out || out->listed || out->num_dropped == 0 || out->rows_left)
    {
        for (ST_int g = 0; g < G; g++)
            free(rows_left[g]);
        free(rows_left);
        return;
    }
    out->rows_left = rows_left;
    out->counted = 1;
    out->G = G;
}

/* A factor's rows per level are counted in row blocks, each with a table of
 * its own, while the tables hold this many entries at most */
#ifndef HDFE_COUNT_BLOCK_ENTRIES /* tests set a small one */
#define HDFE_COUNT_BLOCK_ENTRIES ((size_t)1 << 21)
#endif
/* Slots of a block's table of levels first met in earlier blocks (it gives
 * up when half are taken) */
#ifndef HDFE_COUNT_FAR_SLOTS /* tests set a small one; a power of two */
#define HDFE_COUNT_FAR_SLOTS 8192
#endif

/* Whether a factor's rows are counted already (the build takes them from
 * singleton removal when that counted every level): its first level has rows
 * in any sample, and the counts are zero until they are filled */
static int rows_counted(const FE_Factor *f)
{
    return f->num_levels > 0 && f->counts && f->counts[0] != 0.0;
}

int hdfe_count_rows_parallel(const HDFE_State *S, ST_int g, ST_int local)
{
    if (rows_counted(&S->factors[g]))
        return 1; /* nothing left for hdfe_count_rows() to do */
    const int B = ctools_hdfe_row_blocks(S->N);
    if (B <= 1)
        return 0;
    return local || (size_t)S->factors[g].num_levels * (size_t)B <= HDFE_COUNT_BLOCK_ENTRIES;
}

/* One serial pass, run by run: rows sorted by the factor then cost no chain
 * of updates of the same count (whole numbers add up exactly in any
 * grouping) */
static void count_rows_serial(FE_Factor *f, ST_int N)
{
    if (N <= 0)
        return;
    const ST_int *lev = f->levels;
    ST_double *counts = f->counts;
    ST_int level = lev[0], run = 1;
    for (ST_int i = 1; i < N; i++)
    {
        if (lev[i] == level)
        {
            run++;
            continue;
        }
        counts[level - 1] += (ST_double)run;
        level = lev[i];
        run = 1;
    }
    counts[level - 1] += (ST_double)run;
}

/* Few levels: a table per row block, summed */
static int count_rows_tables(FE_Factor *f, ST_int N, int B)
{
    const size_t L = (size_t)f->num_levels;
    /* tables at least two cache lines apart */
    const size_t stride = L < 32 ? 32 : L;
    ST_int *priv = (ST_int *)calloc((size_t)B * stride, sizeof(ST_int));
    const ST_int *lev = f->levels;
    int b;
    long long l;
    if (!priv)
        return 1;
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        ST_int *p = priv + (size_t)b * stride;
        ST_int i0, i1;
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
            p[lev[i] - 1]++;
    }
#pragma omp parallel for schedule(static) if (L >= 65536)
    for (l = 0; l < (long long)L; l++)
    {
        ST_int n = 0;
        for (int t = 0; t < B; t++)
            n += priv[(size_t)t * stride + (size_t)l];
        f->counts[l] = (ST_double)n;
    }
    free(priv);
    return 0;
}

/*
 * Many levels, in rows that mostly bring new levels block by block (rows
 * sorted by the factor, levels numbered in order of first appearance): the
 * levels above the largest of the earlier blocks belong to a block alone,
 * which counts them in place; the others, first met in an earlier block,
 * go to a small table of the block and are added at the end. Returns 1,
 * with the counts zero again, when a block has too many of the latter.
 */
static int count_rows_own_levels(FE_Factor *f, ST_int N, int B)
{
    const ST_int *lev = f->levels;
    ST_double *counts = f->counts;
    const size_t slots = HDFE_COUNT_FAR_SLOTS;
    ST_int *top = (ST_int *)calloc((size_t)B + 1, sizeof(ST_int));
    ST_int *earlier = (ST_int *)calloc((size_t)B * slots * 2, sizeof(ST_int)); /* level, rows */
    int b, failed = 0;
    if (!top || !earlier)
    {
        free(top);
        free(earlier);
        return 1;
    }
#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        ST_int i0, i1, most = 0;
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
            if (lev[i] > most)
                most = lev[i];
        top[b + 1] = most;
    }
    for (b = 0; b < B; b++)
        if (top[b + 1] < top[b])
            top[b + 1] = top[b];

#pragma omp parallel for schedule(static)
    for (b = 0; b < B; b++)
    {
        ST_int *table = earlier + (size_t)b * slots * 2;
        const ST_int floor_level = top[b];
        ST_int i0, i1, used = 0;
        ctools_hdfe_block_range(N, B, b, &i0, &i1);
        ST_int i = i0;
        while (i < i1)
        {
            const ST_int level = lev[i];
            ST_int run = 1;
            while (i + run < i1 && lev[i + run] == level)
                run++;
            i += run;
            if (level > floor_level)
            {
                counts[level - 1] += (ST_double)run;
                continue;
            }
            size_t h = ((size_t)(uint32_t)level * 2654435761u) & (slots - 1);
            while (table[2 * h] && table[2 * h] != level)
                h = (h + 1) & (slots - 1);
            if (!table[2 * h])
            {
                if ((size_t)++used > slots / 2 || CTOOLS_HDFE_LOAD(&failed))
                {
                    CTOOLS_HDFE_STORE(&failed, 1);
                    break;
                }
                table[2 * h] = level;
            }
            table[2 * h + 1] += run;
        }
    }
    if (failed)
    {
        long long l;
        const long long L = f->num_levels;
#pragma omp parallel for schedule(static)
        for (l = 0; l < L; l++)
            counts[l] = 0.0;
    }
    else
    {
        for (size_t k = 0; k < (size_t)B * slots; k++)
            if (earlier[2 * k])
                counts[earlier[2 * k] - 1] += (ST_double)earlier[2 * k + 1];
    }
    free(top);
    free(earlier);
    return failed;
}

/* Many levels in rows of no order, in T ranges of levels with a thread each
 * (whole numbers add up exactly in any grouping) */
static void count_rows_ranges(FE_Factor *f, ST_int N, int T)
{
    const ST_int *lev = f->levels;
    ST_double *counts = f->counts;
    const ST_int L = f->num_levels;
    int t;
#pragma omp parallel for schedule(static)
    for (t = 0; t < T; t++)
    {
        uint32_t lo, width;
        ctools_hdfe_level_range(L, T, t, &lo, &width);
        CTOOLS_HDFE_ROWS_IN_RANGE(lev, N, lo, width, 0, counts[lev[row] - 1] += 1.0);
    }
}

void hdfe_count_rows(HDFE_State *S, ST_int g, ST_int local)
{
    FE_Factor *f = &S->factors[g];
    const ST_int N = S->N;
    const int B = ctools_hdfe_row_blocks(N);
    if (rows_counted(f))
        return;
    if (B > 1 && (size_t)f->num_levels * (size_t)B <= HDFE_COUNT_BLOCK_ENTRIES)
    {
        if (count_rows_tables(f, N, B) == 0)
            return;
    }
    else if (B > 1 && local)
    {
        if (count_rows_own_levels(f, N, B) == 0)
            return;
    }
    count_rows_serial(f, N);
}

void hdfe_count_rows_serial(HDFE_State *S, ST_int g)
{
    if (!rows_counted(&S->factors[g]))
        count_rows_serial(&S->factors[g], S->N);
}

void hdfe_count_weights(HDFE_State *S, ST_int g, ST_int N_in, const ST_int *keep,
                        const ST_double *weights)
{
    FE_Factor *f = &S->factors[g];
    if (!weights || !f->weighted_counts)
        return;
    const ST_int *lev = f->levels;
    ST_double *sums = f->weighted_counts;
    if (!keep)
    {
        for (ST_int i = 0; i < N_in; i++)
            sums[lev[i] - 1] += weights[i];
        return;
    }
    ST_int idx = 0;
    for (ST_int i = 0; i < N_in; i++)
    {
        if (keep[i])
            sums[lev[idx++] - 1] += weights[i];
    }
}

void hdfe_count_factors(HDFE_State *S, const HDFE_FactorOrder *order, ST_int N_in,
                        const ST_int *keep, const ST_double *weights)
{
    const ST_int G = S->G;
    const int B = ctools_hdfe_row_blocks(S->N);
    ST_int g, job, jobs = weights ? 2 * G : G, serial_jobs = weights ? G : 0;
    /* Which factors' rows are counted in row blocks: decided here, outside
     * the parallel loop below (inside it there are no row blocks) */
    unsigned char *blocks = (unsigned char *)calloc((size_t)(G > 0 ? G : 1), 1);
    for (g = 0; g < G; g++)
    {
        if (blocks)
            blocks[g] = (unsigned char)hdfe_count_rows_parallel(S, g, order ? order[g].rows.local : 0);
        serial_jobs += !(blocks && blocks[g]);
    }
    if (B > 1 && serial_jobs == 1 && !weights)
    {
        /* One factor with many levels in rows of no order, and nothing to
         * run beside it: all the threads count it, by ranges of levels */
        for (g = 0; g < G; g++)
            if (!(blocks && blocks[g]) && !rows_counted(&S->factors[g]))
                count_rows_ranges(&S->factors[g], S->N, B);
        jobs = 0;
    }
    /* The serial passes side by side: rows of the factors with many
     * levels in no order, and the weights of every factor */
#pragma omp parallel for schedule(dynamic) if (B > 1 && serial_jobs > 1)
    for (job = 0; job < jobs; job++)
    {
        if (job >= G)
            hdfe_count_weights(S, job - G, N_in, keep, weights);
        else if (!(blocks && blocks[job]))
            hdfe_count_rows_serial(S, job);
    }
    for (g = 0; g < G; g++)
        if (blocks && blocks[g])
            hdfe_count_rows(S, g, order ? order[g].rows.local : 0);
    free(blocks);
}

/* A level whose first row was dropped but which has rows left: its first
 * kept row puts it in the new order */
typedef struct
{
    ST_int key;   /* first kept row */
    ST_int level; /* loaded number */
} ordered_moved;

static int ordered_compare_moved(const void *a, const void *b)
{
    const ST_int x = ((const ordered_moved *)a)->key, y = ((const ordered_moved *)b)->key;
    return (x > y) - (x < y);
}

/* The kept rows of lev[from..to) to lev[at...], through newid when given.
 * dropped[d...] are the dropped rows from `from` on. Rows are written at or
 * before where they are read. Returns the place after the last one. */
static ST_int ordered_shift_rows(ST_int *lev, const ST_int *src, ST_int src_first, ST_int from,
                                 ST_int to, ST_int at, const ST_int *dropped, ST_int n_dropped,
                                 ST_int d, const ST_int *newid)
{
    ST_int i = from;
    while (i < to)
    {
        const ST_int stop = (d < n_dropped && dropped[d] < to) ? dropped[d] : to;
        const ST_int *in = src + (i - src_first);
        const ST_int len = stop - i;
        if (newid)
        {
            for (ST_int r = 0; r < len; r++)
                lev[at + r] = newid[in[r]];
        }
        else if (len > 0 && lev + at != in)
        {
            memmove(lev + at, in, (size_t)len * sizeof(ST_int));
        }
        at += len;
        i = stop + 1; /* past the dropped row, if any */
        d++;
    }
    return at;
}

/*
 * The levels of one factor on the kept rows, renumbered in order of first
 * appearance among them, from levels already in order of first appearance
 * over all rows (lev, 1..L, first rows first_row), in place: lev then holds
 * the N kept rows. Only a level that lost every row (it goes) or its first
 * row (it moves to where its first kept row puts it) changes the order: the
 * others keep theirs, so the new numbers are the old ones shifted, found by
 * merging the moved levels into the rest by first row. Returns 0, 198 or
 * 920.
 */
static ST_retcode ordered_shift_levels(ST_int *lev, ST_int L, const ST_int *first_row,
                                       const HDFE_Singletons *sing,
                                       const HDFE_AffectedLevels *aff, ST_int N_in, int B,
                                       ST_int *num_levels_out)
{
    ST_retcode rc = 0;
    const ST_int *keep = sing->keep, *dropped = sing->dropped;
    const ST_int n_dropped = sing->num_dropped;
    ST_int nv = 0, nm = 0, k;
    ST_int *mark = NULL, *gone_ids = NULL, *moved_ids = NULL, *moved_keys = NULL, *newid = NULL;
    ST_int *first_kept = NULL, *heads = NULL, *head_at = NULL;
    uint64_t *bits = NULL;
    ordered_moved *moved = NULL;
    int b;

    for (k = 0; k < aff->n; k++)
    {
        if (aff->left[k] == 0)
            nv++;
        else if (!keep[first_row[aff->level[k] - 1]])
            nm++;
    }

    if (nv + nm > 0)
    {
        /* mark[level]: -1 gone, j + 1 moved (entry j of moved), 0 otherwise */
        mark = (ST_int *)calloc((size_t)L + 2, sizeof(ST_int));
        gone_ids = (ST_int *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(ST_int));
        moved_ids = (ST_int *)malloc((size_t)(nm > 0 ? nm : 1) * sizeof(ST_int));
        moved_keys = (ST_int *)malloc((size_t)(nm > 0 ? nm : 1) * sizeof(ST_int));
        moved = (ordered_moved *)malloc((size_t)(nm > 0 ? nm : 1) * sizeof(ordered_moved));
        newid = (ST_int *)malloc(((size_t)L + 1) * sizeof(ST_int));
        if (!mark || !gone_ids || !moved_ids || !moved_keys || !moved || !newid)
        {
            rc = 920;
            goto done;
        }
        ST_int iv = 0, im = 0;
        for (k = 0; k < aff->n; k++)
        {
            const ST_int level = aff->level[k];
            if (aff->left[k] == 0)
            {
                mark[level] = -1;
                gone_ids[iv++] = level;
            }
            else if (!keep[first_row[level - 1]])
            {
                moved[im].level = level;
                moved[im].key = -1;
                mark[level] = ++im;
            }
        }

        if (nm > 0)
        {
            /* First kept row of the moved levels: the first block that has
             * one, and its first */
            bits = (uint64_t *)calloc((size_t)L / 64 + 1, sizeof(uint64_t));
            first_kept = (ST_int *)malloc((size_t)B * (size_t)nm * sizeof(ST_int));
            if (!bits || !first_kept)
            {
                rc = 920;
                goto done;
            }
            for (k = 0; k < nm; k++)
                bits[(size_t)moved[k].level >> 6] |= (uint64_t)1 << (moved[k].level & 63);
#pragma omp parallel for schedule(static) if (B > 1)
            for (b = 0; b < B; b++)
            {
                ST_int *fk = first_kept + (size_t)b * (size_t)nm;
                ST_int i0, i1;
                ctools_hdfe_block_range(N_in, B, b, &i0, &i1);
                for (ST_int j = 0; j < nm; j++)
                    fk[j] = -1;
                for (ST_int i = i0; i < i1; i++)
                {
                    const ST_int l = lev[i];
                    if (!((bits[(size_t)l >> 6] >> (l & 63)) & 1))
                        continue;
                    const ST_int j = mark[l] - 1;
                    if (fk[j] < 0 && keep[i])
                        fk[j] = i;
                }
            }
            for (k = 0; k < nm; k++)
            {
                for (b = 0; b < B && moved[k].key < 0; b++)
                    moved[k].key = first_kept[(size_t)b * (size_t)nm + (size_t)k];
                if (moved[k].key < 0)
                {
                    rc = 198;
                    goto done;
                }
            }
            for (k = 0; k < nm; k++)
                moved_ids[k] = moved[k].level;
            qsort(moved_ids, (size_t)nm, sizeof(ST_int), sing_compare_rows);
            qsort(moved, (size_t)nm, sizeof(ordered_moved), ordered_compare_moved);
            for (k = 0; k < nm; k++)
                moved_keys[k] = moved[k].key;
        }
        qsort(gone_ids, (size_t)nv, sizeof(ST_int), sing_compare_rows);

        /* New numbers, the old ones in C ranges: before level a come the
         * kept unmoved levels below it and the moved levels whose first
         * kept row is before a's first row; each range ends with the moved
         * levels up to the next range's first row */
        const int C = (B > 1 && L >= HDFE_ORDER_RANGE_LEVELS) ? B : 1;
        int c;
#pragma omp parallel for schedule(static) if (C > 1)
        for (c = 0; c < C; c++)
        {
            const ST_int a = 1 + (ST_int)(((long long)L * c) / C);
            const ST_int end = 1 + (ST_int)(((long long)L * (c + 1)) / C);
            ST_int j = a <= L ? ordered_below(moved_keys, nm, first_row[a - 1]) : nm;
            ST_int next = (a - 1) - ordered_below(gone_ids, nv, a) - ordered_below(moved_ids, nm, a) + j;
            for (ST_int o = a; o < end; o++)
            {
                if (mark[o] != 0)
                {
                    if (mark[o] < 0)
                        newid[o] = 0;
                    continue;
                }
                while (j < nm && moved_keys[j] < first_row[o - 1])
                    newid[moved[j++].level] = ++next;
                newid[o] = ++next;
            }
            const ST_int j_end = end <= L ? ordered_below(moved_keys, nm, first_row[end - 1]) : nm;
            while (j < j_end)
                newid[moved[j++].level] = ++next;
        }
    }

    /* The kept rows' levels, shifted down over the dropped rows, in Bs row
     * blocks, each with more rows than there are dropped rows. A block's
     * first rows land in the block before it, where that block may not have
     * read its last rows yet: every block sets as many of its first rows
     * aside as rows are dropped before its end, shifts the others (they stay
     * inside the block, behind the row being read), and the rows set aside
     * are placed once all blocks are done. */
    int Bs = B;
    while (Bs > 1 && N_in / Bs <= n_dropped)
        Bs--;
    head_at = (ST_int *)calloc((size_t)Bs + 1, sizeof(ST_int));
    if (!head_at)
    {
        rc = 920;
        goto done;
    }
    for (b = 1; b < Bs; b++)
    {
        ST_int i0, i1;
        ctools_hdfe_block_range(N_in, Bs, b, &i0, &i1);
        head_at[b + 1] = head_at[b] + ordered_below(dropped, n_dropped, i1);
    }
    heads = (ST_int *)malloc((size_t)(head_at[Bs] > 0 ? head_at[Bs] : 1) * sizeof(ST_int));
    if (!heads)
    {
        rc = 920;
        goto done;
    }
#pragma omp parallel for schedule(static) if (Bs > 1)
    for (b = 0; b < Bs; b++)
    {
        ST_int i0, i1;
        ctools_hdfe_block_range(N_in, Bs, b, &i0, &i1);
        const ST_int head = head_at[b + 1] - head_at[b];
        if (ordered_below(dropped, n_dropped, i1) == 0 && !newid)
            continue; /* no row above moves */
        if (head > 0)
            memcpy(heads + head_at[b], lev + i0, (size_t)head * sizeof(ST_int));
        const ST_int from = i0 + head, d = ordered_below(dropped, n_dropped, from);
        ordered_shift_rows(lev, lev + from, from, from, i1, from - d, dropped, n_dropped, d, newid);
    }
#pragma omp parallel for schedule(static) if (Bs > 1)
    for (b = 1; b < Bs; b++)
    {
        ST_int i0, i1;
        ctools_hdfe_block_range(N_in, Bs, b, &i0, &i1);
        const ST_int head = head_at[b + 1] - head_at[b];
        if (head == 0)
            continue;
        const ST_int d = ordered_below(dropped, n_dropped, i0);
        ordered_shift_rows(lev, heads + head_at[b], i0, i0, i0 + head, i0 - d, dropped, n_dropped, d,
                           newid);
    }
    *num_levels_out = L - nv;

done:
    free(mark);
    free(gone_ids);
    free(moved_ids);
    free(moved_keys);
    free(moved);
    free(newid);
    free(first_kept);
    free(bits);
    free(heads);
    free(head_at);
    return rc;
}

/* One factor in build_factors_ranked() */
typedef struct
{
    const ST_int *lev;       /* loaded levels, 1..L in order of first appearance */
    const ST_int *first_row; /* first row of each */
    const ST_int *left;      /* rows left per level (entry 2 * l), or NULL */
    ST_int L;
    size_t words, stride;      /* of a set of marks (a bit per level); of a block's */
    uint64_t *gone, *moved;    /* levels without a row left; with rows but not the first */
    uint64_t *marks;           /* per row block: its moved levels, then those it has first */
    ST_int *unmoved;           /* [c]: levels keeping their first row in the ranges before c */
    ST_int *found;             /* [b]: moved levels first met in the blocks before b */
    ST_int *key, *moved_level; /* the moved levels by first kept row, and those rows */
    ST_int *newid;             /* new number of each loaded level with a row left */
    ST_int nm, L_new;
    int C; /* ranges of levels that hold whole words of the marks */
} rank_factor;

/* Levels [*from, *stop) of range c of f */
static inline void rank_range(const rank_factor *f, int c, ST_int *from, ST_int *stop)
{
    const size_t w0 = f->words * (size_t)c / (size_t)f->C;
    const size_t w1 = f->words * ((size_t)c + 1) / (size_t)f->C;
    *from = w0 == 0 ? 1 : (ST_int)(64 * w0);
    *stop = 64 * w1 > (size_t)f->L ? f->L + 1 : (ST_int)(64 * w1);
}

static void rank_factor_free(rank_factor *f)
{
    free(f->gone);
    free(f->moved);
    free(f->marks);
    free(f->unmoved);
    free(f->found);
    free(f->key);
    free(f->moved_level);
    free(f->newid);
}

/*
 * hdfe_build_factors() for levels loaded in order of first appearance, when
 * the rows dropped are many and counted (sing->counted), with the same
 * result. The new numbers follow the first kept rows. A level whose first
 * row is kept comes where it came among such levels. One that lost its
 * first row but has rows left comes where its first kept row puts it, and
 * only those are looked for among the rows: each block of rows marks the
 * ones it has, a level stays marked in the first block that has it, and that
 * block notes the level's first kept row there. Taken block after block
 * those rows are ascending, and the numbers follow from merging them into
 * the first rows of the other levels, in ranges of levels. The kept rows
 * then read their levels through the new numbers, and the rows per level are
 * taken from sing where it counted them (the other counts are left zero for
 * the caller). Three passes over the rows, for all factors at once. Returns
 * 0, 198, 920, or 1 with S untouched when the general build must do it (the
 * marks would exceed their budget, or sing does not add up).
 */
static ST_retcode build_factors_ranked(HDFE_State *S, ST_int G, const HDFE_LoadedFactor *factors,
                                       const HDFE_FactorOrder *order, const HDFE_Singletons *sing,
                                       ST_int N_in, ST_int N, const ST_double *weights,
                                       ST_int weighted, int B)
{
    ST_retcode rc = 0;
    const ST_int *keep = sing->keep;
    ST_int g, any_moved = 0;
    int b, c;
    rank_factor *rf = (rank_factor *)calloc((size_t)G, sizeof(rank_factor));
    ST_int *kept = (ST_int *)calloc((size_t)B + 1, sizeof(ST_int));
    ST_int **out = (ST_int **)calloc((size_t)G, sizeof(ST_int *));
    FE_Factor *built = NULL;
    if (!rf || !kept || !out)
    {
        rc = 920;
        goto done;
    }

    /* The levels with a row left: first row kept, or moved */
    for (g = 0; g < G && !rc; g++)
    {
        rank_factor *f = &rf[g];
        long long gone_levels = 0, rows_left = 0;
        f->lev = factors[g].levels;
        f->first_row = order[g].rows.first_row;
        f->left = sing->rows_left ? sing->rows_left[g] : NULL;
        f->L = factors[g].num_levels;
        f->words = (size_t)f->L / 64 + 1;
        /* a block's marks at least two cache lines from the next block's */
        f->stride = f->words < 16 ? 16 : f->words;
        f->C = (B > 1 && f->L >= HDFE_ORDER_RANGE_LEVELS) ? B : 1;
        if ((size_t)f->C > f->words)
            f->C = (int)f->words;
        f->gone = (uint64_t *)calloc(f->words, sizeof(uint64_t));
        f->moved = (uint64_t *)calloc(f->words, sizeof(uint64_t));
        f->unmoved = (ST_int *)calloc((size_t)f->C + 1, sizeof(ST_int));
        f->newid = (ST_int *)malloc(((size_t)f->L + 1) * sizeof(ST_int));
        ST_int *moved_in = (ST_int *)calloc((size_t)f->C + 1, sizeof(ST_int));
        if (!f->gone || !f->moved || !f->unmoved || !f->newid || !moved_in)
        {
            free(moved_in);
            rc = 920;
            break;
        }
#pragma omp parallel for schedule(static) reduction(+ : gone_levels, rows_left) if (f->C > 1)
        for (c = 0; c < f->C; c++)
        {
            ST_int a, end, u = 0, m = 0;
            rank_range(f, c, &a, &end);
            for (ST_int o = a; o < end; o++)
            {
                const int first_kept = keep[f->first_row[o - 1]] != 0;
                /* without counts every dropped row was alone in its level */
                const ST_int left = f->left ? f->left[2 * (size_t)o] : first_kept;
                rows_left += left;
                if (left == 0)
                {
                    level_bit_set(f->gone, o);
                    gone_levels++;
                }
                else if (first_kept)
                    u++;
                else
                {
                    level_bit_set(f->moved, o);
                    m++;
                }
            }
            f->unmoved[c + 1] = u;
            moved_in[c + 1] = m;
        }
        for (c = 0; c < f->C; c++)
        {
            f->unmoved[c + 1] += f->unmoved[c];
            f->nm += moved_in[c + 1];
        }
        free(moved_in);
        f->L_new = f->unmoved[f->C] + f->nm;
        /* What sing says must add up: the kept rows, or one level gone for
         * every row dropped */
        if (f->left ? rows_left != (long long)N : gone_levels != (long long)(N_in - N))
            rc = 1;
        if (f->nm > 0 && !rc)
        {
            if (f->words > HDFE_BUILD_BITSET_BUDGET / 8 / (size_t)B)
            {
                rc = 1;
                break;
            }
            any_moved = 1;
            f->marks = (uint64_t *)calloc((size_t)B * f->stride, sizeof(uint64_t));
            f->found = (ST_int *)calloc((size_t)B + 1, sizeof(ST_int));
            f->key = (ST_int *)malloc((size_t)f->nm * sizeof(ST_int));
            f->moved_level = (ST_int *)malloc((size_t)f->nm * sizeof(ST_int));
            if (!f->marks || !f->found || !f->key || !f->moved_level)
                rc = 920;
        }
    }
    if (rc)
        goto done;

    /* First pass: the kept rows of each block, and the moved levels among
     * them */
#pragma omp parallel for schedule(static) if (B > 1)
    for (b = 0; b < B; b++)
    {
        ST_int i0, i1, n = 0;
        ctools_hdfe_block_range(N_in, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
            n += keep[i] != 0;
        kept[b + 1] = n;
        for (ST_int h = 0; h < G; h++)
        {
            const rank_factor *f = &rf[h];
            if (f->nm == 0)
                continue;
            uint64_t *m = f->marks + (size_t)b * f->stride;
            const uint64_t *moved = f->moved;
            const ST_int *lev = f->lev;
            for (ST_int i = i0; i < i1; i++)
            {
                const ST_int l = lev[i];
                if (level_bit(moved, l) && keep[i] && !level_bit(m, l))
                    level_bit_set(m, l);
            }
        }
    }
    for (b = 0; b < B; b++)
        kept[b + 1] += kept[b];
    if (kept[B] != N)
    {
        rc = 198;
        goto done;
    }

    if (any_moved)
    {
        /* Each moved level in its first block only */
        for (g = 0; g < G; g++)
        {
            rank_factor *f = &rf[g];
            long long w;
            int bad = 0;
            if (f->nm == 0)
                continue;
#pragma omp parallel for schedule(static) reduction(| : bad) if (B > 1 && f->words >= 4096)
            for (w = 0; w < (long long)f->words; w++)
            {
                uint64_t seen = 0;
                for (int t = 0; t < B; t++)
                {
                    uint64_t *slot = &f->marks[(size_t)t * f->stride + (size_t)w];
                    const uint64_t v = *slot;
                    *slot = v & ~seen;
                    seen |= v;
                }
                bad |= seen != f->moved[w];
            }
#pragma omp parallel for schedule(static) if (B > 1 && f->words >= 4096)
            for (b = 0; b < B; b++)
            {
                const uint64_t *m = f->marks + (size_t)b * f->stride;
                ST_int n = 0;
                for (size_t k = 0; k < f->words; k++)
                    n += __builtin_popcountll(m[k]);
                f->found[b + 1] = n;
            }
            for (b = 0; b < B; b++)
                f->found[b + 1] += f->found[b];
            if (bad || f->found[B] != f->nm)
            {
                /* a level counted with rows left has no kept row */
                rc = 1;
                goto done;
            }
        }
        /* Second pass: their first kept rows, ascending block after block */
#pragma omp parallel for schedule(static) if (B > 1)
        for (b = 0; b < B; b++)
        {
            ST_int i0, i1;
            ctools_hdfe_block_range(N_in, B, b, &i0, &i1);
            for (ST_int h = 0; h < G; h++)
            {
                rank_factor *f = &rf[h];
                if (f->nm == 0)
                    continue;
                uint64_t *m = f->marks + (size_t)b * f->stride;
                const ST_int *lev = f->lev;
                const ST_int stop = f->found[b + 1];
                ST_int j = f->found[b];
                for (ST_int i = i0; i < i1 && j < stop; i++)
                {
                    const ST_int l = lev[i];
                    uint64_t *word = &m[(size_t)l >> 6];
                    const uint64_t bit = (uint64_t)1 << (l & 63);
                    if ((*word & bit) && keep[i])
                    {
                        *word &= ~bit;
                        f->key[j] = i;
                        f->moved_level[j++] = l;
                    }
                }
            }
        }
    }

    /* New numbers: before an unmoved level come the unmoved levels below it
     * and the moved levels whose first kept row is before its first row;
     * each range ends with the moved levels up to the next range's first row */
    for (g = 0; g < G; g++)
    {
        rank_factor *f = &rf[g];
#pragma omp parallel for schedule(static) if (f->C > 1)
        for (c = 0; c < f->C; c++)
        {
            ST_int a, end;
            rank_range(f, c, &a, &end);
            ST_int j = a <= f->L ? ordered_below(f->key, f->nm, f->first_row[a - 1]) : f->nm;
            ST_int next = f->unmoved[c] + j;
            for (ST_int o = a; o < end; o++)
            {
                if (level_bit(f->gone, o) || level_bit(f->moved, o))
                    continue;
                while (j < f->nm && f->key[j] < f->first_row[o - 1])
                    f->newid[f->moved_level[j++]] = ++next;
                f->newid[o] = ++next;
            }
            const ST_int j_end =
                end <= f->L ? ordered_below(f->key, f->nm, f->first_row[end - 1]) : f->nm;
            while (j < j_end)
                f->newid[f->moved_level[j++]] = ++next;
        }
    }

    built = (FE_Factor *)calloc((size_t)G, sizeof(FE_Factor));
    if (!built)
    {
        rc = 920;
        goto done;
    }
    for (g = 0; g < G; g++)
    {
        const size_t L = (size_t)(rf[g].L_new > 0 ? rf[g].L_new : 1);
        out[g] = (ST_int *)malloc((size_t)N * sizeof(ST_int));
        built[g].counts = (ST_double *)calloc(L, sizeof(ST_double));
        if (weights || weighted)
            built[g].weighted_counts = (ST_double *)calloc(L, sizeof(ST_double));
        if (!out[g] || !built[g].counts || ((weights || weighted) && !built[g].weighted_counts))
            rc = 920;
    }
    if (rc)
        goto done;

    /* Third pass: the kept rows' levels */
#pragma omp parallel for schedule(static) if (B > 1)
    for (b = 0; b < B; b++)
    {
        ST_int i0, i1, idx = kept[b];
        ctools_hdfe_block_range(N_in, B, b, &i0, &i1);
        for (ST_int i = i0; i < i1; i++)
        {
            if (!keep[i])
                continue;
            for (ST_int h = 0; h < G; h++)
                out[h][idx] = rf[h].newid[rf[h].lev[i]];
            idx++;
        }
    }

    for (g = 0; g < G; g++)
    {
        FE_Factor *f = &built[g];
        f->levels = out[g];
        out[g] = NULL;
        f->num_levels = rf[g].L_new;
        f->max_level = rf[g].L_new;
        f->has_intercept = 1;
        if (rf[g].left)
        {
            /* every level was counted by the removal */
            const ST_int *left = rf[g].left, *id = rf[g].newid;
            const uint64_t *gone = rf[g].gone;
            const long long L = rf[g].L;
            long long l;
#pragma omp parallel for schedule(static) if (B > 1 && L >= HDFE_ORDER_RANGE_LEVELS)
            for (l = 1; l <= L; l++)
                if (!level_bit(gone, (ST_int)l))
                    f->counts[id[l] - 1] = (ST_double)left[2 * (size_t)l];
        }
    }
    S->G = G;
    S->N = N;
    S->factors = built;
    built = NULL;

done:
    if (built)
        for (g = 0; g < G; g++)
        {
            free(built[g].counts);
            free(built[g].weighted_counts);
        }
    free(built);
    if (rf)
        for (g = 0; g < G; g++)
            rank_factor_free(&rf[g]);
    if (out)
        for (g = 0; g < G; g++)
            free(out[g]);
    free(rf);
    free(out);
    free(kept);
    return rc;
}

ST_retcode hdfe_build_factors_ordered(HDFE_State *S, ST_int G, HDFE_LoadedFactor *factors,
                                      const HDFE_FactorOrder *order, const HDFE_Singletons *sing,
                                      ST_int N_in, const ST_double *weights, ST_int weighted,
                                      ST_int *counts_pending)
{
    ST_retcode rc = 0;
    ST_int g;
    const ST_int n_dropped = sing ? sing->num_dropped : 0;
    const ST_int *keep = sing ? sing->keep : NULL;
    const ST_int N = N_in - n_dropped;

    if (counts_pending)
        *counts_pending = 0;
    if (!S || S->factors || G < 0 || N_in < 0 || N < 0 || (G > 0 && (!factors || !order)))
        return 198;

    /* The loaded levels serve when they are in order of first appearance
     * and the dropped rows, if any, are listed with what the load recorded */
    int direct = G > 0 && N > 0 &&
                 (n_dropped == 0 || (sing->listed && sing->dropped && sing->affected));
    for (g = 0; g < G && direct; g++)
        direct = order[g].ordered && factors[g].levels && (n_dropped == 0 || order[g].rows.first_row);
    if (!direct)
    {
        /* Many rows dropped, and counted, from levels in order of first
         * appearance */
        int ranked = G > 0 && N > 0 && n_dropped > 0 && keep && sing->counted;
        for (g = 0; g < G && ranked; g++)
            ranked = order[g].ordered && factors[g].levels && order[g].rows.first_row;
        if (ranked)
        {
            rc = build_factors_ranked(S, G, factors, order, sing, N_in, N, weights, weighted,
                                      ctools_hdfe_row_blocks(N_in));
            if (rc == 0)
            {
                if (counts_pending)
                    *counts_pending = 1;
                else
                    hdfe_count_factors(S, order, N_in, keep, weights);
                return 0;
            }
            if (rc != 1)
                return rc;
            rc = 0;
        }
        ST_int *level_ptrs[64], level_max[64];
        if (G > 64)
            return 198;
        for (g = 0; g < G; g++)
        {
            level_ptrs[g] = factors[g].levels;
            level_max[g] = factors[g].num_levels;
        }
        return hdfe_build_factors(S, G, level_ptrs, level_max, N_in, keep, N, weights, weighted);
    }

    S->G = G;
    S->N = N;
    S->factors = (FE_Factor *)calloc((size_t)G, sizeof(FE_Factor));
    if (!S->factors)
        return 920;

    const int B = ctools_hdfe_row_blocks(N_in);
    for (g = 0; g < G && !rc; g++)
    {
        FE_Factor *f = &S->factors[g];
        ST_int L = factors[g].num_levels;
        if (n_dropped > 0)
        {
            rc = ordered_shift_levels(factors[g].levels, L, order[g].rows.first_row, sing,
                                      &sing->affected[g], N_in, B, &L);
            if (rc)
                break;
        }
        /* the loaded array is the factor's from here on */
        f->levels = factors[g].levels;
        factors[g].levels = NULL;
        f->num_levels = L;
        f->max_level = L;
        f->has_intercept = 1;
        f->counts = (ST_double *)calloc((size_t)(L > 0 ? L : 1), sizeof(ST_double));
        if (weights || weighted)
            f->weighted_counts = (ST_double *)calloc((size_t)(L > 0 ? L : 1), sizeof(ST_double));
        if (!f->counts || ((weights || weighted) && !f->weighted_counts))
            rc = 920;
    }
    if (rc)
    {
        ctools_hdfe_state_cleanup(S);
        return rc;
    }

    if (counts_pending)
    {
        *counts_pending = 1;
        return 0;
    }
    hdfe_count_factors(S, order, N_in, keep, weights);
    return 0;
}

/*
    cjoinby_plan.c
    Match counts and output layout for cjoinby

    The master keys are probed against the using index in the untouched
    master, so a join that cannot be made (mismatched key types, too many
    output observations) fails with nothing to undo.

    The output keeps the master's order: each master observation is followed
    by its pairings, in the using data's order. Key groups that no master
    observation has (unmatched(using|both)) follow the master observations in
    key order. joinby is a sortpreserve program: when the master is sorted,
    its result is sorted the same way, with the using-only observations
    placed by their key values and missing values of the other sort
    variables. cj_plan_place reproduces that order in the layout itself, so
    no sort of the joined data is needed afterwards.
*/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ctools_config.h"
#include "cjoinby_internal.h"

/* Rows a worker should have before another worker is worth waking */
#define CJ_PLAN_ROWS_PER_THREAD 32768
/* Team sizes the per-worker totals are laid out for */
#define CJ_PLAN_MAX_TEAM 256

void cj_plan_init(cj_plan *p)
{
    if (p)
        memset(p, 0, sizeof(*p));
}

void cj_plan_free(cj_plan *p)
{
    if (!p)
        return;
    free(p->gid);
    free(p->off);
    free(p->row);
    cj_plan_init(p);
}

static int cj_plan_threads(size_t n)
{
    int threads = ctools_get_max_threads();
    if (threads > CJ_PLAN_MAX_TEAM)
        threads = CJ_PLAN_MAX_TEAM;
    if (n / CJ_PLAN_ROWS_PER_THREAD < (size_t)threads)
        threads = (int)(n / CJ_PLAN_ROWS_PER_THREAD);
    return threads < 1 ? 1 : threads;
}

/* ----------------------------------------------------------------------------
 * Placement of using-only groups in a sorted master
 * -------------------------------------------------------------------------- */

typedef struct
{
    const cj_index *ix;
    const stata_data *mkeys;
    const cj_sortspec *sort;
} cj_place_ctx;

/* qsort has no context argument; plan calls are made on Stata's thread only */
static const cj_place_ctx *g_place;

/* Using-only groups a and b by the master's sort variables: a key compares by
 * its value; any other sort variable is missing in both, hence equal. */
static int cj_place_group_cmp(const void *pa, const void *pb)
{
    const uint32_t a = *(const uint32_t *)pa, b = *(const uint32_t *)pb;
    const cj_index *ix = g_place->ix;
    const cj_sortspec *sort = g_place->sort;
    for (int s = 0; s < sort->nsort; s++)
    {
        const int col = sort->col[s];
        if (col >= ix->nkeys)
            continue;
        const cj_keycol *kc = &ix->keytab[col];
        int c;
        if (kc->is_str)
        {
            c = strcmp(kc->str + (size_t)a * kc->slot, kc->str + (size_t)b * kc->slot);
        }
        else
        {
            double u = kc->num[a], v = kc->num[b];
            c = (u > v) - (u < v);
        }
        if (c)
            return c;
    }
    return (a > b) - (a < b);
}

/* Master row i against using-only group g by the master's sort variables. A
 * using-only observation holds its key values and, in every other master
 * variable, a missing value: . for numbers, "" for strings. */
static int cj_place_master_cmp(const cj_place_ctx *ctx, size_t i, uint32_t g)
{
    const cj_index *ix = ctx->ix;
    const cj_sortspec *sort = ctx->sort;
    for (int s = 0; s < sort->nsort; s++)
    {
        const int col = sort->col[s];
        const stata_variable *mv = &ctx->mkeys->vars[col];
        int c;
        if (mv->type == STATA_TYPE_STRING)
        {
            const char *other =
                col < ix->nkeys ? ix->keytab[col].str + (size_t)g * ix->keytab[col].slot : "";
            c = strcmp(mv->data.str[i], other);
        }
        else
        {
            double u = mv->data.dbl[i];
            double v = col < ix->nkeys ? ix->keytab[col].num[g] : SV_missval;
            c = (u > v) - (u < v);
        }
        if (c)
            return c;
    }
    return 0;
}

/* Merge the nt using-only groups in tail[] into the master sequence: a stable
 * sort on the sort variables of the master rows followed by the using-only
 * rows. Replaces p->gid by the merged drivers and sets p->row. */
static stata_retcode cj_plan_place(cj_plan *p, const cj_index *ix, const stata_data *mkeys,
                                   const cj_sortspec *sort, uint32_t *tail, size_t nt)
{
    const size_t nm = p->nm, ndrv = nm + nt;
    cj_place_ctx ctx = {ix, mkeys, sort};

    g_place = &ctx;
    int ordered = 1;
    for (size_t t = 1; t < nt && ordered; t++)
        ordered = cj_place_group_cmp(&tail[t - 1], &tail[t]) < 0;
    if (!ordered)
        qsort(tail, nt, sizeof(uint32_t), cj_place_group_cmp);
    g_place = NULL;

    /* pos[t]: master rows that sort at or before the group (ties go first) */
    uint32_t *pos = (uint32_t *)ctools_safe_malloc2(nt ? nt : 1, sizeof(uint32_t));
    uint32_t *gid = (uint32_t *)ctools_safe_malloc2(ndrv ? ndrv : 1, sizeof(uint32_t));
    uint32_t *row = (uint32_t *)ctools_safe_malloc2(ndrv ? ndrv : 1, sizeof(uint32_t));
    if (!pos || !gid || !row)
    {
        free(pos);
        free(gid);
        free(row);
        return STATA_ERR_MEMORY;
    }
    int threads = cj_plan_threads(nt * 16);
    (void)threads;
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
    for (size_t t = 0; t < nt; t++)
    {
        size_t lo = 0, hi = nm;
        while (lo < hi)
        {
            size_t mid = lo + (hi - lo) / 2;
            if (cj_place_master_cmp(&ctx, mid, tail[t]) <= 0)
                lo = mid + 1;
            else
                hi = mid;
        }
        pos[t] = (uint32_t)lo;
    }

    size_t d = 0, t = 0;
    for (size_t i = 0; i <= nm; i++)
    {
        while (t < nt && pos[t] == i)
        {
            gid[d] = tail[t++];
            row[d++] = CJ_NONE;
        }
        if (i < nm)
        {
            gid[d] = p->gid[i];
            row[d++] = (uint32_t)i;
        }
    }
    free(pos);
    free(p->gid);
    p->gid = gid;
    p->row = row;
    p->placed = 1;
    return STATA_OK;
}

/* ----------------------------------------------------------------------------
 * Layout
 * -------------------------------------------------------------------------- */

/* Output rows of every driver as prefix sums, and the first driver that is
 * not a master row staying where it is. Returns 0, or 920 past the limit. */
static int cj_plan_offsets(cj_plan *p, const cj_index *ix, int keep_master)
{
    const size_t ndrv = p->ndrv;
    const uint32_t *gid = p->gid, *row = p->row;
    const uint32_t *gstart = ix->gstart;
    const uint32_t lone = keep_master ? 1 : 0;
    uint32_t *off = p->off;
    uint64_t total[CJ_PLAN_MAX_TEAM + 1];
    size_t first[CJ_PLAN_MAX_TEAM];
    int threads = cj_plan_threads(ndrv), team_used = 1;
    (void)threads;

    memset(total, 0, sizeof(total));
    for (int t = 0; t < CJ_PLAN_MAX_TEAM; t++)
        first[t] = ndrv;
    int overflow = 0;
#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        size_t team = 1, tid = 0;
#ifdef _OPENMP
        team = (size_t)omp_get_num_threads();
        tid = (size_t)omp_get_thread_num();
#endif
        const size_t base = ndrv / team, extra = ndrv % team;
        const size_t start = tid * base + (tid < extra ? tid : extra);
        const size_t end = start + base + (tid < extra);
        uint64_t sum = 0;
        size_t moved = ndrv;
        /* counts first (off[d + 1] holds driver d's rows for now) */
        for (size_t d = start; d < end; d++)
        {
            const uint32_t g = gid[d];
            const uint32_t count = g != CJ_NONE ? gstart[g + 1] - gstart[g] : lone;
            const uint32_t source = row ? row[d] : (d < p->nm ? (uint32_t)d : CJ_NONE);
            off[d + 1] = count;
            sum += count;
            if (moved == ndrv && (count != 1 || source != d))
                moved = d;
        }
        total[tid + 1] = sum;
        first[tid] = moved;
#pragma omp barrier
#pragma omp single
        {
            team_used = (int)team;
            for (size_t t = 0; t < team; t++)
                total[t + 1] += total[t];
            if (total[team] > (uint64_t)INT_MAX)
                overflow = 1;
        }
        if (!overflow)
        {
            uint64_t running = total[tid];
            for (size_t d = start; d < end; d++)
            {
                running += off[d + 1];
                off[d + 1] = (uint32_t)running;
            }
        }
    }
    if (overflow)
        return 920;
    off[0] = 0;
    p->nout = (size_t)total[team_used];

    /* An output row stays in place only while every driver before it did */
    size_t stage_from = ndrv;
    for (int t = 0; t < team_used; t++)
        if (first[t] < stage_from)
            stage_from = first[t];
    if (stage_from > p->nm)
        stage_from = p->nm;
    p->stage_from = stage_from;
    return 0;
}

int cj_plan_build(cj_plan *p, const cj_index *ix, const stata_data *mkeys, size_t nm, int unmatched,
                  const cj_sortspec *sort)
{
    if (!p || !ix || (nm && !mkeys))
        return 198;
    cj_plan_init(p);
    p->nm = nm;
    const size_t ng = ix->ngroups;
    const int keep_master = (unmatched & CJ_KEEP_MASTER) != 0;
    const int keep_using = (unmatched & CJ_KEEP_USING) != 0;
    int rc = 0;

    unsigned char *matched = (unsigned char *)calloc(ng ? ng : 1, 1);
    p->gid = (uint32_t *)ctools_safe_malloc2(nm ? nm : 1, sizeof(uint32_t));
    if (!matched || !p->gid)
    {
        free(matched);
        cj_plan_free(p);
        return 920;
    }
    if (nm)
    {
        stata_retcode prc = cj_index_probe(ix, mkeys, p->gid, matched);
        if (prc)
        {
            free(matched);
            cj_plan_free(p);
            return prc == STATA_ERR_INVALID_INPUT ? 106 : ctools_stata_rc(prc);
        }
    }

    /* Key groups without a master observation, in key order */
    size_t nt = 0;
    uint32_t *tail = NULL;
    if (keep_using && ng)
    {
        for (size_t g = 0; g < ng; g++)
            nt += !matched[g];
        tail = (uint32_t *)ctools_safe_malloc2(nt ? nt : 1, sizeof(uint32_t));
        if (!tail)
        {
            free(matched);
            cj_plan_free(p);
            return 920;
        }
        size_t t = 0;
        for (size_t g = 0; g < ng; g++)
            if (!matched[g])
                tail[t++] = (uint32_t)g;
    }
    free(matched);

    p->ndrv = nm + nt;
    if (nt && sort && sort->nsort > 0)
    {
        stata_retcode prc = cj_plan_place(p, ix, mkeys, sort, tail, nt);
        if (prc)
            rc = ctools_stata_rc(prc);
    }
    else if (nt)
    {
        uint32_t *grown = (uint32_t *)realloc(p->gid, p->ndrv * sizeof(uint32_t));
        if (!grown)
            rc = 920;
        else
        {
            p->gid = grown;
            memcpy(grown + nm, tail, nt * sizeof(uint32_t));
        }
    }
    if (!rc)
    {
        p->off = (uint32_t *)ctools_safe_malloc2(p->ndrv + 1, sizeof(uint32_t));
        if (!p->off)
            rc = 920;
    }
    if (!rc)
        rc = cj_plan_offsets(p, ix, keep_master);
    if (!rc)
    {
        /* using-only rows, then the master rows left alone, then the pairs */
        uint64_t n_using = 0, lone = 0;
        for (size_t t = 0; t < nt; t++)
            n_using += cj_group_count(ix, tail[t]);
        int threads = cj_plan_threads(p->ndrv);
        (void)threads;
        const uint32_t *gid = p->gid;
        const size_t ndrv = p->ndrv;
#pragma omp parallel for schedule(static) reduction(+ : lone) num_threads(threads) if (threads > 1)
        for (size_t d = 0; d < ndrv; d++)
            lone += gid[d] == CJ_NONE;
        p->n_using = (size_t)n_using;
        p->n_master = keep_master ? (size_t)lone : 0;
        p->n_both = p->nout - p->n_using - p->n_master;
    }
    free(tail);
    if (rc)
        cj_plan_free(p);
    return rc;
}

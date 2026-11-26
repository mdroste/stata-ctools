/*
    cjoinby_write.c
    Fused, row-tiled writer for cjoinby

    Everything the join writes goes out in one pass over the output rows.
    Each worker owns a run of rows and moves it in tiles, writing every
    variable of a tile before the next, so each Stata observation is touched
    while it is in cache instead of once per variable. Per tile:

    1. resolve each row's master and using observation from the plan;
    2. gather the using observations (keys and payload) into tile slots;
    3. where master observations move, gather them too, give using-only rows
       their key values, apply update/replace to the shared variables in the
       tile, and store the master variables;
    4. store the new variables, _merge and the source observation numbers.

    Output rows before the plan's stage_from are master observations that stay
    where they are: their master variables are not rewritten.
*/

#include <stdlib.h>
#include <string.h>
#include <stdatomic.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ctools_transport_internal.h"
#include "cjoinby_internal.h"

#define CJ_TILE_BYTES ((size_t)64 * 1024)
#define CJ_TILE_MAX_ROWS 256
#define CJ_TILE_MIN_ROWS 16
/* SPI calls a worker should have before another worker is worth waking */
#define CJ_CELLS_PER_THREAD 65536
/* Using rows are gathered at random: fetch a few rows ahead */
#define CJ_PREFETCH_ROWS 8
#define CJ_PREFETCH_BYTES 1024

/* Keep the first nonzero code */
static void cj_record_error(atomic_int *error, int rc)
{
    int expected = 0;
    if (rc)
        atomic_compare_exchange_strong(error, &expected, rc);
}

/* Store column col of table t, held in cnt consecutive tile slots starting
 * at buf, into plugin variable var from observation obs0 on. */
static inline int cj_store_column(const cj_table *t, size_t col, ST_int var, const char *buf,
                                  size_t cnt, ST_int obs0)
{
    const size_t stride = t->stride;
    const char *slot = buf + t->offset[col];
    int rc = 0;
    if (t->width[col] == 0)
    {
        ST_IIID store_fn = IO_VSTORE_FN;
        for (size_t i = 0; i < cnt && !rc; i++, slot += stride)
            rc = store_fn(var, obs0 + (ST_int)i, *(const double *)(const void *)slot);
    }
    else
    {
        ST_IIIS store_fn = (_stata_)->sstore;
        for (size_t i = 0; i < cnt && !rc; i++, slot += stride)
            rc = store_fn(var, obs0 + (ST_int)i, (char *)slot);
    }
    return rc;
}

/* Copy a using value into a master slot of the same generic type. A string
 * destination holds dst_width bytes and a terminator. */
static inline void cj_copy_value(char *dst, int dst_width, const char *src, int src_width)
{
    if (dst_width == 0)
    {
        memcpy(dst, src, sizeof(double));
    }
    else
    {
        size_t len = strnlen(src, (size_t)src_width);
        if (len > (size_t)dst_width)
            len = (size_t)dst_width;
        memcpy(dst, src, len);
        dst[len] = '\0';
    }
}

/* update: a missing master value takes the using value, whatever it is;
 * replace: so does any master value when the using value is not missing. */
static inline int cj_takes_using(const char *mslot, const char *uslot, int is_str, int replace)
{
    if (is_str)
        return mslot[0] == '\0' || (replace && uslot[0] != '\0');
    double m, u;
    memcpy(&m, mslot, sizeof(double));
    memcpy(&u, uslot, sizeof(double));
    return m >= SV_missval || (replace && u < SV_missval);
}

static stata_retcode cj_write_check(const cj_writer *w)
{
    const cj_table *ut = w->utab, *mt = w->mtab;
    if (!w->ix || !w->plan || !ut || w->limit > w->plan->nout)
        return STATA_ERR_INVALID_INPUT;
    if ((size_t)SF_nobs() < w->limit)
        return STATA_ERR_INVALID_INPUT;
    if (ut->nvars != (size_t)(w->nkeys + w->npay))
        return STATA_ERR_INVALID_INPUT;
    for (int k = 0; k < w->nkeys; k++)
    {
        int var = w->key_mcol[k] + 1;
        if (var < 1 || var > SF_nvars())
            return STATA_ERR_INVALID_INPUT;
        if (mt && (size_t)w->key_mcol[k] >= mt->nvars)
            return STATA_ERR_INVALID_INPUT;
        if (SF_var_is_strl(var) || (SF_var_is_string(var) != 0) != (ut->width[k] > 0))
            return STATA_ERR_UNSUPPORTED_TYPE;
    }
    for (int c = 0; c < w->npay; c++)
    {
        const int ucol = w->nkeys + c;
        int var = w->pay_mcol[c] >= 0 ? w->pay_mcol[c] + 1 : w->pay_var[c];
        if (var < 1 || var > SF_nvars())
            return STATA_ERR_INVALID_INPUT;
        if (mt && w->pay_mcol[c] >= 0 && (size_t)w->pay_mcol[c] >= mt->nvars)
            return STATA_ERR_INVALID_INPUT;
        if (SF_var_is_strl(var) || (SF_var_is_string(var) != 0) != (ut->width[ucol] > 0))
            return STATA_ERR_UNSUPPORTED_TYPE;
    }
    const int extra[3] = {w->merge_var, w->mobs_var, w->uobs_var};
    for (int e = 0; e < 3; e++)
    {
        if (extra[e] == 0)
            continue;
        if (extra[e] < 1 || extra[e] > SF_nvars())
            return STATA_ERR_INVALID_INPUT;
        if (SF_var_is_string(extra[e]))
            return STATA_ERR_UNSUPPORTED_TYPE;
    }
    if (mt)
        for (size_t j = 0; j < mt->nvars; j++)
            if (mt->var_idx[j] != (int)j + 1)
                return STATA_ERR_INVALID_INPUT;
    return STATA_OK;
}

stata_retcode cj_write_output(const cj_writer *w)
{
    stata_retcode rc = cj_write_check(w);
    if (rc)
        return rc;

    const cj_index *ix = w->ix;
    const cj_plan *p = w->plan;
    const cj_table *ut = w->utab, *mt = w->mtab;
    const size_t limit = w->limit;
    const size_t stage_from = w->stage_from;
    const size_t ndrv = p->ndrv;
    const uint32_t *off = p->off, *gid = p->gid, *gstart = ix->gstart;
    const perm_idx_t *order = ix->keys.sort_order;
    const int nkeys = w->nkeys, npay = w->npay;
    const ST_int merge_var = (ST_int)w->merge_var;
    const ST_int mobs_var = (ST_int)w->mobs_var, uobs_var = (ST_int)w->uobs_var;
    if (limit == 0)
        return STATA_OK;

    int n_new = 0, n_shared = 0;
    for (int c = 0; c < npay; c++)
    {
        if (w->pay_mcol[c] < 0)
            n_new++;
        else
            n_shared++;
    }
    const int patch_shared = w->update && n_shared > 0;
    const int has_tail = p->n_using > 0;
    /* the using rows are needed in every tile, or only where using-only rows are */
    const int gather_always = n_new > 0 || patch_shared;

    const size_t mstride = mt ? mt->stride : 0, ustride = ut->stride;
    size_t tile = CJ_TILE_MAX_ROWS;
    {
        size_t widest = mstride > ustride ? mstride : ustride;
        if (widest && CJ_TILE_BYTES / widest < tile)
            tile = CJ_TILE_BYTES / widest;
        if (tile < CJ_TILE_MIN_ROWS)
            tile = CJ_TILE_MIN_ROWS;
    }

    int threads = 1;
#ifdef _OPENMP
    {
        size_t row_cells = (mt ? mt->nvars : 0) + (size_t)n_new + (merge_var > 0) + (mobs_var > 0) +
                           (uobs_var > 0) + 1;
        size_t cells;
        threads = ctools_get_max_threads();
        if (threads < 1)
            threads = 1;
        if (ctools_safe_mul_size(limit, row_cells, &cells) == 0 &&
            cells / CJ_CELLS_PER_THREAD < (size_t)threads)
            threads = (int)(cells / CJ_CELLS_PER_THREAD);
        if (threads < 1)
            threads = 1;
    }
#endif
    (void)threads;

    /* What a row without a source observation receives */
    char *mblank = mt ? (char *)malloc(mstride) : NULL;
    char *ublank = (char *)malloc(ustride);
    if ((mt && !mblank) || !ublank)
    {
        free(mblank);
        free(ublank);
        return STATA_ERR_MEMORY;
    }
    if (mt)
        cj_table_blank(mt, mblank);
    cj_table_blank(ut, ublank);

    const size_t uprefetch = ustride < CJ_PREFETCH_BYTES ? ustride : CJ_PREFETCH_BYTES;
    atomic_int error = 0;

#pragma omp parallel num_threads(threads) if (threads > 1)
    {
        size_t team = 1, tid = 0;
#ifdef _OPENMP
        team = (size_t)omp_get_num_threads();
        tid = (size_t)omp_get_thread_num();
#endif
        const size_t base = limit / team, extra = limit % team;
        const size_t start = tid * base + (tid < extra ? tid : extra);
        const size_t end = start + base + (tid < extra);
        char *mbuf = mt ? (char *)ctools_safe_malloc2(tile, mstride) : NULL;
        char *ubuf = (char *)ctools_safe_malloc2(tile, ustride);
        uint32_t mrow[CJ_TILE_MAX_ROWS], urow[CJ_TILE_MAX_ROWS];
        if ((mt && !mbuf) || !ubuf)
            cj_record_error(&error, STATA_ERR_MEMORY);

        /* The driver of the worker's first row: the last one starting at or
         * before it (drivers without rows are skipped as the rows advance) */
        size_t d = 0;
        if (start < end)
        {
            size_t lo = 0, hi = ndrv;
            while (lo < hi)
            {
                size_t mid = lo + (hi - lo) / 2;
                if (off[mid + 1] <= start)
                    lo = mid + 1;
                else
                    hi = mid;
            }
            d = lo;
        }

        size_t cnt = 0;
        for (size_t t0 = start; t0 < end && !atomic_load(&error); t0 += cnt)
        {
            cnt = end - t0 < tile ? end - t0 : tile;
            /* a tile lies wholly before or wholly after stage_from */
            if (t0 < stage_from && t0 + cnt > stage_from)
                cnt = stage_from - t0;
            const ST_int obs0 = (ST_int)(t0 + 1);
            const int staged = mt != NULL && t0 >= stage_from;
            int io_rc = 0, any_tail = 0;

            /* 1. sources */
            for (size_t i = 0; i < cnt; i++)
            {
                const size_t j = t0 + i;
                while (j >= off[d + 1])
                    d++;
                const uint32_t g = gid[d];
                mrow[i] = cj_plan_master_row(p, d);
                urow[i] = g != CJ_NONE ? (uint32_t)order[gstart[g] + (j - off[d])] : CJ_NONE;
                any_tail |= mrow[i] == CJ_NONE;
            }

            /* 2. using observations */
            const int have_u = gather_always || (has_tail && any_tail);
            if (have_u)
            {
                for (size_t i = 0; i < cnt; i++)
                {
#if defined(__GNUC__) || defined(__clang__)
                    if (ut->rowmajor && i + CJ_PREFETCH_ROWS < cnt &&
                        urow[i + CJ_PREFETCH_ROWS] != CJ_NONE)
                    {
                        const char *ahead =
                            ut->rs.rows + (size_t)urow[i + CJ_PREFETCH_ROWS] * ustride;
                        for (size_t b = 0; b < uprefetch; b += CACHE_LINE_SIZE)
                            __builtin_prefetch(ahead + b, 0, 0);
                    }
#endif
                    if (urow[i] == CJ_NONE)
                        memcpy(ubuf + i * ustride, ublank, ustride);
                    else
                        cj_table_gather(ut, urow[i], ubuf + i * ustride);
                }
            }

            /* 3. master variables */
            if (staged)
            {
                for (size_t i = 0; i < cnt; i++)
                {
                    if (mrow[i] == CJ_NONE)
                        memcpy(mbuf + i * mstride, mblank, mstride);
                    else
                        cj_table_gather(mt, (size_t)mrow[i] - stage_from, mbuf + i * mstride);
                }
                /* using-only rows hold the key values of their group */
                for (size_t i = 0; any_tail && i < cnt; i++)
                {
                    if (mrow[i] != CJ_NONE)
                        continue;
                    for (int k = 0; k < nkeys; k++)
                    {
                        const int mc = w->key_mcol[k];
                        cj_copy_value(mbuf + i * mstride + mt->offset[mc], mt->width[mc],
                                      ubuf + i * ustride + ut->offset[k], ut->width[k]);
                    }
                }
                /* update and replace, decided on the values in the tile */
                for (int c = 0; patch_shared && c < npay; c++)
                {
                    const int mc = w->pay_mcol[c];
                    if (mc < 0)
                        continue;
                    const int ucol = nkeys + c;
                    const int is_str = mt->width[mc] > 0;
                    char *mslot = mbuf + mt->offset[mc];
                    const char *uslot = ubuf + ut->offset[ucol];
                    for (size_t i = 0; i < cnt; i++, mslot += mstride, uslot += ustride)
                    {
                        if (urow[i] == CJ_NONE)
                            continue;
                        if (cj_takes_using(mslot, uslot, is_str, w->replace))
                            cj_copy_value(mslot, mt->width[mc], uslot, ut->width[ucol]);
                    }
                }
                for (size_t j = 0; j < mt->nvars && !io_rc; j++)
                    io_rc = cj_store_column(mt, j, (ST_int)mt->var_idx[j], mbuf, cnt, obs0);
            }
            else if (any_tail)
            {
                /* No master observation moves: a using-only row is a new,
                 * empty observation that needs its keys and, under update,
                 * the using values of the shared variables. */
                for (size_t i = 0; i < cnt && !io_rc; i++)
                {
                    if (mrow[i] != CJ_NONE)
                        continue;
                    for (int k = 0; k < nkeys && !io_rc; k++)
                        io_rc = cj_store_column(ut, (size_t)k, (ST_int)(w->key_mcol[k] + 1),
                                                ubuf + i * ustride, 1, obs0 + (ST_int)i);
                    for (int c = 0; patch_shared && c < npay && !io_rc; c++)
                        if (w->pay_mcol[c] >= 0)
                            io_rc = cj_store_column(ut, (size_t)(nkeys + c),
                                                    (ST_int)(w->pay_mcol[c] + 1),
                                                    ubuf + i * ustride, 1, obs0 + (ST_int)i);
                }
            }

            /* 4. new variables, _merge, source observation numbers */
            for (int c = 0; c < npay && !io_rc; c++)
                if (w->pay_mcol[c] < 0)
                    io_rc = cj_store_column(ut, (size_t)(nkeys + c), (ST_int)w->pay_var[c], ubuf,
                                            cnt, obs0);
            if (merge_var > 0)
            {
                ST_IIID store_fn = IO_VSTORE_FN;
                for (size_t i = 0; i < cnt && !io_rc; i++)
                    io_rc = store_fn(merge_var, obs0 + (ST_int)i,
                                     mrow[i] == CJ_NONE   ? (double)CJ_MERGE_USING
                                     : urow[i] == CJ_NONE ? (double)CJ_MERGE_MASTER
                                                          : (double)CJ_MERGE_BOTH);
            }
            if (mobs_var > 0)
            {
                ST_IIID store_fn = IO_VSTORE_FN;
                for (size_t i = 0; i < cnt && !io_rc; i++)
                    io_rc = store_fn(mobs_var, obs0 + (ST_int)i,
                                     mrow[i] == CJ_NONE ? SV_missval : (double)mrow[i] + 1.0);
            }
            if (uobs_var > 0)
            {
                ST_IIID store_fn = IO_VSTORE_FN;
                for (size_t i = 0; i < cnt && !io_rc; i++)
                    io_rc = store_fn(uobs_var, obs0 + (ST_int)i,
                                     urow[i] == CJ_NONE ? SV_missval : (double)urow[i] + 1.0);
            }
            if (io_rc)
                cj_record_error(&error, ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE));
        }
        free(mbuf);
        free(ubuf);
    }
    free(mblank);
    free(ublank);
    ctools_memory_barrier();
    return (stata_retcode)atomic_load(&error);
}

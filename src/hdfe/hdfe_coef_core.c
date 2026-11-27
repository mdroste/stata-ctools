/*
 * hdfe_coef_core.c
 *
 * Shared pieces of the coefficient-space projection (hdfe_coef.h): rows
 * grouped by an effect, the cross-tabulation of two effects, conjugate
 * gradients for several right-hand sides, and a dense Cholesky factorization
 * Part of the ctools Stata plugin suite
 */

#include "hdfe_coef_internal.h"

/* A pivot at most this fraction of its row's original diagonal is taken to
 * be rounding and its row is fixed at zero. The pivots of the directions the
 * callers know to be singular (which they fix by name) come out between
 * 1e-16 and 4e-11 of the diagonal in either sign (benchmark_hdfe_coef.c, up
 * to 5,000 levels and a million rows), so the threshold catches most but
 * not all of the rest; one that survives only shifts the solution along a
 * direction the residuals do not see. Fixing a row whose small pivot is
 * genuine would be the costlier mistake, hence no larger. */
#define HDFE_COEF_PIVOT_TOL 1e-11

/* The recurrence for the residual sum of squares subtracts each step's
 * decrease from f0, so it resolves nothing below about 1e-15 * f0. */
#define HDFE_COEF_POTENTIAL_FLOOR 1e-14

/* Rounding of a solve relative to the column it starts from. The level
 * means of B, the right-hand side and the residual recurrence are all
 * rounded at the size of the raw column, and S is singular: the part of
 * that rounding outside its range is never reduced, and an iteration that
 * reaches it starts to follow it instead of converging. Measured at about
 * one unit of 2.2e-16 of the column's root sum of squares
 * (benchmark_hdfe_coef.c, three effects, tolerance 1e-13); a level mean
 * over n rows can be off by sqrt(n / 3) units, which the margin covers up
 * to levels of a few million rows. */
#define HDFE_COEF_ROUNDING (1024.0 * 2.220446049250313e-16)

/* Per-thread stamp arrays of the cell merge may use this much memory. */
#define HDFE_COEF_STAMP_BYTES ((size_t)256 << 20)

/* ========================================================================
 * Rows grouped by one effect
 * ======================================================================== */

/* The checks of hdfe_coef_group_rows(): 198 when a level (or a level of
 * second) is out of range; otherwise *unsorted tells whether some level is
 * below the one before it and *repeats whether, among rows in ascending
 * order, a level of second fails to increase within a level. */
static ST_retcode hcc_group_check(ST_int N, const ST_int *levels, ST_int nlev,
                                  const ST_int *second, ST_int nsecond, int *unsorted_out,
                                  int *repeats_out)
{
    int bad = 0, unsorted = 0, repeats = 0;
#pragma omp parallel for schedule(static)                                                          \
    reduction(| : bad, unsorted, repeats) if (N > HDFE_COEF_PARALLEL_MIN)
    for (ST_int i = 0; i < N; i++)
    {
        const ST_int l = levels[i];
        if (l < 1 || l > nlev)
            bad |= 1;
        if (second && (second[i] < 1 || second[i] > nsecond))
            bad |= 1;
        if (i > 0)
        {
            if (l < levels[i - 1])
                unsorted |= 1;
            else if (second && l == levels[i - 1] && second[i] <= second[i - 1])
                repeats |= 1;
        }
    }
    *unsorted_out = unsorted;
    *repeats_out = repeats;
    return bad ? 198 : 0;
}

/* rowptr (nlev + 1) of rows that arrive grouped: the first row of every
 * level; a level without rows starts where the next one does. */
static void hcc_group_sorted(ST_int N, const ST_int *levels, ST_int nlev, ST_int *rowptr)
{
#pragma omp parallel for schedule(static) if (nlev > HDFE_COEF_PARALLEL_MIN)
    for (ST_int l = 0; l < nlev; l++)
        rowptr[l] = -1;
#pragma omp parallel for schedule(static) if (N > HDFE_COEF_PARALLEL_MIN)
    for (ST_int i = 0; i < N; i++)
    {
        if (i == 0 || levels[i] != levels[i - 1])
            rowptr[levels[i] - 1] = i;
    }
    rowptr[nlev] = N;
    for (ST_int l = nlev - 1; l >= 0; l--)
    {
        if (rowptr[l] < 0)
            rowptr[l] = rowptr[l + 1];
    }
}

/* Stable counting sort of rows that do not arrive grouped: the position of
 * every row among the rows ordered by level (*dest_out, N, the caller's to
 * free) and rowptr (nlev + 1, zeroed by the caller). Chunks of rows count
 * and place in parallel; the positions are those of the serial sort whatever
 * the number of chunks. Returns 0 or 920. */
static ST_retcode hcc_group_positions(ST_int N, const ST_int *levels, ST_int nlev, ST_int *rowptr,
                                      ST_int **dest_out)
{
    ST_int *dest = (ST_int *)malloc((size_t)N * sizeof(ST_int));
    *dest_out = NULL;
    if (!dest)
        return 920;
    int nc = hdfe_coef_threads();
    if (nc > 64)
        nc = 64;
    if (N <= HDFE_COEF_PARALLEL_MIN)
        nc = 1;
    while (nc > 1 && (size_t)nc * (size_t)nlev > 2 * (size_t)N + ((size_t)1 << 20))
        nc /= 2;
    ST_int *hist = nc > 1 ? (ST_int *)calloc((size_t)nc * (size_t)nlev, sizeof(ST_int)) : NULL;
    if (hist)
    {
#pragma omp parallel for schedule(static, 1)
        for (int c = 0; c < nc; c++)
        {
            ST_int *h = hist + (size_t)c * (size_t)nlev;
            const ST_int i0 = (ST_int)((long long)N * c / nc);
            const ST_int i1 = (ST_int)((long long)N * (c + 1) / nc);
            for (ST_int i = i0; i < i1; i++)
                h[levels[i] - 1]++;
        }
#pragma omp parallel for schedule(static)
        for (ST_int l = 0; l < nlev; l++)
        {
            ST_int total = 0;
            for (int c = 0; c < nc; c++)
                total += hist[(size_t)c * (size_t)nlev + (size_t)l];
            rowptr[l + 1] = total;
        }
        for (ST_int l = 0; l < nlev; l++)
            rowptr[l + 1] += rowptr[l];
#pragma omp parallel for schedule(static)
        for (ST_int l = 0; l < nlev; l++)
        {
            ST_int run = rowptr[l];
            for (int c = 0; c < nc; c++)
            {
                ST_int *cell = &hist[(size_t)c * (size_t)nlev + (size_t)l];
                const ST_int count = *cell;
                *cell = run;
                run += count;
            }
        }
#pragma omp parallel for schedule(static, 1)
        for (int c = 0; c < nc; c++)
        {
            ST_int *h = hist + (size_t)c * (size_t)nlev;
            const ST_int i0 = (ST_int)((long long)N * c / nc);
            const ST_int i1 = (ST_int)((long long)N * (c + 1) / nc);
            for (ST_int i = i0; i < i1; i++)
                dest[i] = h[levels[i] - 1]++;
        }
        free(hist);
    }
    else
    {
        for (ST_int i = 0; i < N; i++)
            rowptr[levels[i]]++;
        for (ST_int l = 0; l < nlev; l++)
            rowptr[l + 1] += rowptr[l];
        /* rowptr[l] walks through level l's positions, then is restored */
        for (ST_int i = 0; i < N; i++)
            dest[i] = rowptr[levels[i] - 1]++;
        for (ST_int l = nlev - 1; l > 0; l--)
            rowptr[l] = rowptr[l - 1];
        rowptr[0] = 0;
    }
    *dest_out = dest;
    return 0;
}

ST_retcode hdfe_coef_group_rows(ST_int N, const ST_int *levels, ST_int nlev, const ST_int *second,
                                ST_int nsecond, ST_int **rowptr_out, ST_int **dest_out,
                                int *sorted_out, int *increasing_out)
{
    *rowptr_out = NULL;
    *dest_out = NULL;
    *sorted_out = 0;
    if (increasing_out)
        *increasing_out = 0;
    if (N < 1 || nlev < 1 || !levels)
        return 198;

    int unsorted = 0, repeats = 0;
    const ST_retcode check =
        hcc_group_check(N, levels, nlev, second, nsecond, &unsorted, &repeats);
    if (check)
        return check;

    ST_int *rowptr = (ST_int *)calloc((size_t)nlev + 1, sizeof(ST_int));
    if (!rowptr)
        return 920;

    if (!unsorted)
    {
        hcc_group_sorted(N, levels, nlev, rowptr);
        *rowptr_out = rowptr;
        *sorted_out = 1;
        if (increasing_out)
            *increasing_out = second && !repeats;
        return 0;
    }

    const ST_retcode placed = hcc_group_positions(N, levels, nlev, rowptr, dest_out);
    if (placed)
    {
        free(rowptr);
        return placed;
    }
    *rowptr_out = rowptr;
    return 0;
}

/* ========================================================================
 * Rows grouped by one effect, their fields moved along
 * ======================================================================== */

/* Ranges of levels that the rows are dealt into, at most, and the rows aimed
 * for in one: a range is put in order within a window of that many rows,
 * which should stay in the second-level cache, while every thread that deals
 * keeps a line of each destination open per range. */
#define HDFE_COEF_MOVE_PARTS 512
#define HDFE_COEF_MOVE_PART_ROWS 16384

/* Fewer rows than this are placed directly through their positions: the
 * destinations then stay in the cache as they are, and dealing first only
 * adds a pass (measured with 16 threads: the two even at 2 million rows,
 * dealing ahead from 4 million; with a level for every ten rows it is ahead
 * at 2 million already). */
#define HDFE_COEF_MOVE_MIN_ROWS ((ST_int)1 << 21)

/* The fields of row i to position at */
HDFE_COEF_KERNEL void hcc_move_row(const hdfe_coef_fields *F, ST_int i, ST_int at)
{
    for (int a = 0; a < F->nint; a++)
        F->idst[a][at] = F->isrc[a][i];
    if (F->want_order)
        F->order[at] = i;
    if (F->dsrc)
        F->ddst[at] = F->dsrc[i];
}

/*
 * The first write to every page of a block that malloc() just returned, a
 * stretch of the block per thread. The system maps a page at its first
 * write, and where the threads that deal reach new pages of several hundred
 * destinations in no order, those faults hold one another up: dealing ten
 * million rows into two fresh arrays took 9.1 ms with 16 threads against 4.7
 * ms into arrays written before, while this pass over them costs 1 ms
 * (Windows 11). The values written are overwritten.
 */
static void hcc_first_touch(void *block, size_t bytes)
{
    char *RESTRICT p = (char *)block;
    const long long pages = (long long)(bytes >> 12);
    if (!p || !bytes)
        return;
#pragma omp parallel for schedule(static) if (bytes > ((size_t)4 << 20))
    for (long long g = 0; g < pages; g++)
        p[(size_t)g << 12] = 0;
    p[bytes - 1] = 0;
}

/* x[0..n) rearranged: the value at j goes to at[j], a permutation of 0..n-1.
 * scratch holds n values. */
HDFE_COEF_KERNEL void hcc_move_ints(ST_int *RESTRICT x, const ST_int *RESTRICT at, ST_int n,
                                    ST_int *RESTRICT scratch)
{
    memcpy(scratch, x, (size_t)n * sizeof(ST_int));
    for (ST_int j = 0; j < n; j++)
        x[at[j]] = scratch[j];
}

HDFE_COEF_KERNEL void hcc_move_doubles(ST_double *RESTRICT x, const ST_int *RESTRICT at, ST_int n,
                                       ST_double *RESTRICT scratch)
{
    memcpy(scratch, x, (size_t)n * sizeof(ST_double));
    for (ST_int j = 0; j < n; j++)
        x[at[j]] = scratch[j];
}

static void hcc_move_free(hdfe_coef_fields *F)
{
    for (int a = 0; a < HDFE_COEF_MAX_EFFECTS; a++)
    {
        free(F->idst[a]);
        F->idst[a] = NULL;
    }
    free(F->ddst);
    free(F->order);
    F->ddst = NULL;
    F->order = NULL;
}

ST_retcode hdfe_coef_group_move(ST_int N, const ST_int *levels, ST_int nlev, const ST_int *second,
                                ST_int nsecond, hdfe_coef_fields *F, ST_int **rowptr_out,
                                int *sorted_out, int *increasing_out)
{
    *rowptr_out = NULL;
    *sorted_out = 0;
    if (increasing_out)
        *increasing_out = 0;
    for (int a = 0; a < HDFE_COEF_MAX_EFFECTS; a++)
        F->idst[a] = NULL;
    F->ddst = NULL;
    F->order = NULL;
    if (N < 1 || nlev < 1 || !levels || F->nint < 0 || F->nint > HDFE_COEF_MAX_EFFECTS)
        return 198;

    int unsorted = 0, repeats = 0;
    const ST_retcode check =
        hcc_group_check(N, levels, nlev, second, nsecond, &unsorted, &repeats);
    if (check)
        return check;

    ST_int *rowptr = (ST_int *)calloc((size_t)nlev + 1, sizeof(ST_int));
    if (!rowptr)
        return 920;
    if (!unsorted)
    {
        hcc_group_sorted(N, levels, nlev, rowptr);
        *rowptr_out = rowptr;
        *sorted_out = 1;
        if (increasing_out)
            *increasing_out = second && !repeats;
        return 0;
    }

    const int nint = F->nint;
    int ok = 1;
    for (int a = 0; a < nint; a++)
    {
        F->idst[a] = (ST_int *)malloc((size_t)N * sizeof(ST_int));
        ok &= F->idst[a] != NULL;
    }
    if (F->dsrc)
    {
        F->ddst = (ST_double *)malloc((size_t)N * sizeof(ST_double));
        ok &= F->ddst != NULL;
    }
    if (F->want_order)
    {
        F->order = (ST_int *)malloc((size_t)N * sizeof(ST_int));
        ok &= F->order != NULL;
    }
    if (!ok)
    {
        free(rowptr);
        hcc_move_free(F);
        return 920;
    }
    /* a copy that the loops below can keep in registers */
    const hdfe_coef_fields L = *F;

    if (N < HDFE_COEF_MOVE_MIN_ROWS)
    {
        ST_int *dest = NULL;
        const ST_retcode placed = hcc_group_positions(N, levels, nlev, rowptr, &dest);
        if (placed)
        {
            free(rowptr);
            hcc_move_free(F);
            return placed;
        }
#pragma omp parallel for schedule(static) if (N > HDFE_COEF_PARALLEL_MIN)
        for (ST_int i = 0; i < N; i++)
            hcc_move_row(&L, i, dest[i]);
        free(dest);
        *rowptr_out = rowptr;
        return 0;
    }

    int nc = hdfe_coef_threads();
    if (nc > 64)
        nc = 64;

    /* The ranges: the levels l with the same (l - 1) >> shift */
    ST_int want = N / HDFE_COEF_MOVE_PART_ROWS;
    if (want > HDFE_COEF_MOVE_PARTS)
        want = HDFE_COEF_MOVE_PARTS;
    if (want < 1)
        want = 1;
    int shift = 0;
    while (((nlev - 1) >> shift) + 1 > want)
        shift++;
    const ST_int P = ((nlev - 1) >> shift) + 1;
    const ST_int span = (ST_int)1 << shift;
    /* per chunk of rows, its rows of every range; the chunks' a few lines apart */
    const size_t stride = (size_t)P + 32;
    ST_int *cnt = (ST_int *)calloc((size_t)nc * stride, sizeof(ST_int));
    ST_int *tstart = (ST_int *)malloc(((size_t)P + 1) * sizeof(ST_int));
    /* the level of every dealt row, unless a range is one level */
    ST_int *tlev = shift > 0 ? (ST_int *)malloc((size_t)N * sizeof(ST_int)) : NULL;
    if (!cnt || !tstart || (shift > 0 && !tlev))
    {
        /* No memory to deal with: a counting sort that places every row
         * directly, in row order. */
        free(cnt);
        free(tstart);
        free(tlev);
        for (ST_int i = 0; i < N; i++)
            rowptr[levels[i]]++;
        for (ST_int l = 0; l < nlev; l++)
            rowptr[l + 1] += rowptr[l];
        /* rowptr[l] walks through level l's positions, then is restored */
        for (ST_int i = 0; i < N; i++)
            hcc_move_row(&L, i, rowptr[levels[i] - 1]++);
        for (ST_int l = nlev - 1; l > 0; l--)
            rowptr[l] = rowptr[l - 1];
        rowptr[0] = 0;
        *rowptr_out = rowptr;
        return 0;
    }

    for (int a = 0; a < nint; a++)
        hcc_first_touch(L.idst[a], (size_t)N * sizeof(ST_int));
    hcc_first_touch(L.order, (size_t)N * sizeof(ST_int));
    hcc_first_touch(L.ddst, (size_t)N * sizeof(ST_double));
    hcc_first_touch(tlev, (size_t)N * sizeof(ST_int));
#pragma omp parallel for schedule(static, 1)
    for (int c = 0; c < nc; c++)
    {
        ST_int *RESTRICT h = cnt + (size_t)c * stride;
        const ST_int i0 = (ST_int)((long long)N * c / nc);
        const ST_int i1 = (ST_int)((long long)N * (c + 1) / nc);
        for (ST_int i = i0; i < i1; i++)
            h[(levels[i] - 1) >> shift]++;
    }
    /* A range's rows lie together, chunk after chunk: in row order, and
     * where they will stay, since the ranges follow one another by level */
    {
        ST_int run = 0;
        for (ST_int p = 0; p < P; p++)
        {
            tstart[p] = run;
            for (int c = 0; c < nc; c++)
            {
                ST_int *cell = &cnt[(size_t)c * stride + (size_t)p];
                const ST_int count = *cell;
                *cell = run;
                run += count;
            }
        }
        tstart[P] = N;
    }
#pragma omp parallel for schedule(static, 1)
    for (int c = 0; c < nc; c++)
    {
        ST_int *RESTRICT cur = cnt + (size_t)c * stride;
        const ST_int i0 = (ST_int)((long long)N * c / nc);
        const ST_int i1 = (ST_int)((long long)N * (c + 1) / nc);
        if (tlev)
        {
            for (ST_int i = i0; i < i1; i++)
            {
                const ST_int l = levels[i];
                const ST_int at = cur[(l - 1) >> shift]++;
                tlev[at] = l;
                hcc_move_row(&L, i, at);
            }
        }
        else
        {
            for (ST_int i = i0; i < i1; i++)
                hcc_move_row(&L, i, cur[levels[i] - 1]++);
        }
    }
    /* Every range by itself: the rows of each of its levels counted, then
     * every field rearranged within the range's window, a row going behind
     * the rows of its level dealt before it, which are the rows before it.
     * A range's window and its levels' entries of rowptr are its own. */
    int failed = 0;
    if (!tlev)
    {
        for (ST_int l = 0; l < nlev; l++)
            rowptr[l] = tstart[l];
    }
    else
    {
#pragma omp parallel reduction(| : failed)
        {
            ST_int *RESTRICT cur = (ST_int *)malloc((size_t)span * sizeof(ST_int));
            ST_int *RESTRICT at = NULL;
            ST_double *RESTRICT scratch = NULL;
            ST_int room = 0;
            if (!cur)
                failed |= 1;
#pragma omp for schedule(dynamic, 1)
            for (ST_int p = 0; p < P; p++)
            {
                if (failed)
                    continue;
                const ST_int l0 = (ST_int)((long long)p << shift);
                const ST_int nl = nlev - l0 < span ? nlev - l0 : span;
                const ST_int r0 = tstart[p], n = tstart[p + 1] - r0;
                if (nl == 1)
                {
                    rowptr[l0] = r0;
                    continue;
                }
                if (n > room)
                {
                    /* a window larger than any before (levels of very
                     * unequal size make some far larger than aimed for) */
                    free(at);
                    free(scratch);
                    at = (ST_int *)malloc((size_t)n * sizeof(ST_int));
                    scratch = (ST_double *)malloc((size_t)n * sizeof(ST_double));
                    room = n;
                    if (!at || !scratch)
                    {
                        room = 0;
                        failed |= 1;
                        continue;
                    }
                }
                const ST_int *RESTRICT lv = tlev + r0;
                memset(cur, 0, (size_t)nl * sizeof(ST_int));
                for (ST_int j = 0; j < n; j++)
                    cur[lv[j] - 1 - l0]++;
                ST_int pos = 0;
                for (ST_int j = 0; j < nl; j++)
                {
                    const ST_int count = cur[j];
                    rowptr[l0 + j] = r0 + pos;
                    cur[j] = pos;
                    pos += count;
                }
                for (ST_int j = 0; j < n; j++)
                    at[j] = cur[lv[j] - 1 - l0]++;
                for (int a = 0; a < nint; a++)
                    hcc_move_ints(L.idst[a] + r0, at, n, (ST_int *)scratch);
                if (L.order)
                    hcc_move_ints(L.order + r0, at, n, (ST_int *)scratch);
                if (L.ddst)
                    hcc_move_doubles(L.ddst + r0, at, n, scratch);
            }
            free(cur);
            free(at);
            free(scratch);
        }
    }
    rowptr[nlev] = N;
    free(cnt);
    free(tstart);
    free(tlev);
    if (failed)
    {
        free(rowptr);
        hcc_move_free(F);
        return 920;
    }
    *rowptr_out = rowptr;
    return 0;
}

/* ========================================================================
 * Cross-tabulation of two effects
 * ======================================================================== */

void hdfe_crosstab_free(hdfe_crosstab *C)
{
    if (!C)
        return;
    free(C->rowptr);
    free(C->col_own);
    free(C->val_own);
    memset(C, 0, sizeof(*C));
}

/* Threads that may each hold `arrays` stamp arrays over ncol levels */
static int hcc_stamp_threads(ST_int ncol, int arrays)
{
    int t = hdfe_coef_threads();
    const size_t each = (size_t)ncol * sizeof(ST_int) * (size_t)arrays;
    while (t > 1 && each * (size_t)t > HDFE_COEF_STAMP_BYTES)
        t--;
    return t;
}

ST_retcode hdfe_crosstab_build(hdfe_crosstab *C, ST_int N, const ST_int *row_levels, ST_int nrow,
                               const ST_int *col_levels, ST_int ncol, const ST_double *w)
{
    memset(C, 0, sizeof(*C));
    C->nrow = nrow;
    C->ncol = ncol;
    if (!col_levels || ncol < 1)
        return 198;

    ST_int *len = NULL, *newptr = NULL, *mcol = NULL;
    ST_double *mval = NULL;
    int sorted = 0, increasing = 0;
    /* The rows in the order of the row effect */
    hdfe_coef_fields moved;
    memset(&moved, 0, sizeof(moved));
    moved.nint = 1;
    moved.isrc[0] = col_levels;
    moved.dsrc = w;
    ST_retcode rc = hdfe_coef_group_move(N, row_levels, nrow, col_levels, ncol, &moved,
                                         &C->rowptr, &sorted, &increasing);
    if (rc)
        return rc;
    C->sorted = sorted;
    const int par = N > HDFE_COEF_PARALLEL_MIN;

    const ST_int *src_col = col_levels;
    const ST_double *src_val = w;
    if (!sorted)
    {
        C->col_own = moved.idst[0];
        C->val_own = moved.ddst;
        src_col = C->col_own;
        src_val = C->val_own;
    }
    C->col = src_col;
    C->val = src_val;
    C->nnz = N;
    if (sorted && increasing)
        return 0;

    /* Distinct cells of every row */
    const ST_int *rowptr = C->rowptr;
    len = (ST_int *)malloc((size_t)nrow * sizeof(ST_int));
    if (!len)
    {
        rc = 920;
        goto fail;
    }
    long long total = 0;
    int failed = 0;
    {
        const int tuse = hcc_stamp_threads(ncol, 1);
#pragma omp parallel num_threads(tuse) reduction(+ : total) reduction(| : failed) if (par)
        {
            ST_int *stamp = (ST_int *)malloc((size_t)ncol * sizeof(ST_int));
            if (!stamp)
                failed |= 1;
            else
                memset(stamp, 0xFF, (size_t)ncol * sizeof(ST_int));
#pragma omp for schedule(dynamic, 1024)
            for (ST_int r = 0; r < nrow; r++)
            {
                if (!stamp)
                    continue;
                ST_int distinct = 0;
                for (ST_int e = rowptr[r]; e < rowptr[r + 1]; e++)
                {
                    const ST_int j = src_col[e] - 1;
                    if (stamp[j] != r)
                    {
                        stamp[j] = r;
                        distinct++;
                    }
                }
                len[r] = distinct;
                total += distinct;
            }
            free(stamp);
        }
    }
    if (failed)
    {
        rc = 920;
        goto fail;
    }
    if (total == (long long)N)
    {
        free(len);
        return 0;
    }

    /* Merge the rows of each cell, cells in order of first appearance and
     * each cell's weights added in row order */
    newptr = (ST_int *)malloc(((size_t)nrow + 1) * sizeof(ST_int));
    mcol = (ST_int *)malloc((size_t)total * sizeof(ST_int));
    mval = (ST_double *)malloc((size_t)total * sizeof(ST_double));
    if (!newptr || !mcol || !mval)
    {
        rc = 920;
        goto fail;
    }
    newptr[0] = 0;
    for (ST_int r = 0; r < nrow; r++)
        newptr[r + 1] = newptr[r] + len[r];
    {
        const int tuse = hcc_stamp_threads(ncol, 2);
#pragma omp parallel num_threads(tuse) reduction(| : failed) if (par)
        {
            ST_int *stamp = (ST_int *)malloc((size_t)ncol * sizeof(ST_int));
            ST_int *slot = (ST_int *)malloc((size_t)ncol * sizeof(ST_int));
            if (!stamp || !slot)
                failed |= 1;
            else
                memset(stamp, 0xFF, (size_t)ncol * sizeof(ST_int));
#pragma omp for schedule(dynamic, 1024)
            for (ST_int r = 0; r < nrow; r++)
            {
                if (!stamp || !slot)
                    continue;
                ST_int out = newptr[r];
                for (ST_int e = rowptr[r]; e < rowptr[r + 1]; e++)
                {
                    const ST_int j = src_col[e] - 1;
                    const ST_double v = src_val ? src_val[e] : 1.0;
                    if (stamp[j] == r)
                        mval[slot[j]] += v;
                    else
                    {
                        stamp[j] = r;
                        slot[j] = out;
                        mcol[out] = j + 1;
                        mval[out] = v;
                        out++;
                    }
                }
            }
            free(stamp);
            free(slot);
        }
    }
    if (failed)
    {
        rc = 920;
        goto fail;
    }
    free(len);
    free(C->rowptr);
    free(C->col_own);
    free(C->val_own);
    C->rowptr = newptr;
    C->col_own = mcol;
    C->val_own = mval;
    C->col = mcol;
    C->val = mval;
    C->nnz = (ST_int)total;
    return 0;

fail:
    free(len);
    free(newptr);
    free(mcol);
    free(mval);
    hdfe_crosstab_free(C);
    return rc;
}

ST_int *hdfe_coef_chunks(const ST_int *rowptr, const ST_int *rows, ST_int n, int nchunks)
{
    ST_int *bound = (ST_int *)malloc(((size_t)nchunks + 1) * sizeof(ST_int));
    if (!bound)
        return NULL;
    long long total = 0;
    for (ST_int a = 0; a < n; a++)
    {
        const ST_int r = rows ? rows[a] : a;
        total += rowptr[r + 1] - rowptr[r];
    }
    bound[0] = 0;
    long long seen = 0;
    int c = 1;
    for (ST_int a = 0; a < n && c < nchunks; a++)
    {
        const ST_int r = rows ? rows[a] : a;
        seen += rowptr[r + 1] - rowptr[r];
        while (c < nchunks && seen * nchunks >= total * c)
            bound[c++] = a + 1;
    }
    while (c <= nchunks)
        bound[c++] = n;
    return bound;
}

/* ========================================================================
 * A cross-tabulation kept while the weights change
 * ======================================================================== */

/* Chunks of levels of the row effect that the weights are summed in */
#define HDFE_COEF_VALUE_CHUNKS 256
#define HDFE_COEF_VALUE_CHUNK_ROWS 8192

void hdfe_crossplan_free(hdfe_crossplan *P)
{
    if (!P)
        return;
    free(P->grouprow);
    free(P->cellptr);
    free(P->order);
    free(P->cell);
    free(P->col);
    free(P->val);
    free(P->chunk);
    memset(P, 0, sizeof(*P));
}

/* As hdfe_crosstab_build(), cell for cell, without the weights: the cells of
 * a level in the order their first rows appear. */
ST_retcode hdfe_crossplan_build(hdfe_crossplan *P, ST_int N, const ST_int *row_levels,
                                ST_int nrow, const ST_int *col_levels, ST_int ncol)
{
    memset(P, 0, sizeof(*P));
    hdfe_crosstab *C = &P->C;
    C->nrow = nrow;
    C->ncol = ncol;
    P->N = N;
    if (!col_levels || ncol < 1)
        return 198;

    ST_int *len = NULL;
    int sorted = 0, increasing = 0;
    /* the rows in the order of the row effect, and where each came from */
    hdfe_coef_fields moved;
    memset(&moved, 0, sizeof(moved));
    moved.nint = 1;
    moved.isrc[0] = col_levels;
    moved.want_order = 1;
    ST_retcode rc = hdfe_coef_group_move(N, row_levels, nrow, col_levels, ncol, &moved,
                                         &P->grouprow, &sorted, &increasing);
    if (rc)
        return rc;
    C->sorted = sorted;
    const int par = N > HDFE_COEF_PARALLEL_MIN;

    const ST_int *src_col = col_levels;
    if (!sorted)
    {
        P->order = moved.order;
        P->col = moved.idst[0];
        src_col = P->col;
    }
    C->rowptr = P->grouprow;
    C->col = src_col;
    C->nnz = N;

    {
        int nchunk = N / HDFE_COEF_VALUE_CHUNK_ROWS;
        if (nchunk > HDFE_COEF_VALUE_CHUNKS)
            nchunk = HDFE_COEF_VALUE_CHUNKS;
        if (nchunk < 1)
            nchunk = 1;
        P->chunk = hdfe_coef_chunks(P->grouprow, NULL, nrow, nchunk);
        P->nchunk = nchunk;
        if (!P->chunk)
        {
            rc = 920;
            goto fail;
        }
    }

    const ST_int *rowptr = P->grouprow;
    long long total = N;
    int failed = 0;
    if (!(sorted && increasing))
    {
        /* Distinct cells of every row */
        len = (ST_int *)malloc((size_t)nrow * sizeof(ST_int));
        if (!len)
        {
            rc = 920;
            goto fail;
        }
        total = 0;
        const int tuse = hcc_stamp_threads(ncol, 1);
#pragma omp parallel num_threads(tuse) reduction(+ : total) reduction(| : failed) if (par)
        {
            ST_int *stamp = (ST_int *)malloc((size_t)ncol * sizeof(ST_int));
            if (!stamp)
                failed |= 1;
            else
                memset(stamp, 0xFF, (size_t)ncol * sizeof(ST_int));
#pragma omp for schedule(dynamic, 1024)
            for (ST_int r = 0; r < nrow; r++)
            {
                if (!stamp)
                    continue;
                ST_int distinct = 0;
                for (ST_int e = rowptr[r]; e < rowptr[r + 1]; e++)
                {
                    const ST_int j = src_col[e] - 1;
                    if (stamp[j] != r)
                    {
                        stamp[j] = r;
                        distinct++;
                    }
                }
                len[r] = distinct;
                total += distinct;
            }
            free(stamp);
        }
        if (failed)
        {
            rc = 920;
            goto fail;
        }
    }
    if (total == (long long)N)
    {
        /* one row per cell: the weights in grouped order, or as they are */
        free(len);
        len = NULL;
        if (!sorted)
        {
            P->val = (ST_double *)malloc((size_t)N * sizeof(ST_double));
            if (!P->val)
            {
                rc = 920;
                goto fail;
            }
        }
        return 0;
    }

    /* The cell of every row, cells in order of first appearance */
    {
        ST_int *newptr = (ST_int *)malloc(((size_t)nrow + 1) * sizeof(ST_int));
        ST_int *mcol = (ST_int *)malloc((size_t)total * sizeof(ST_int));
        P->cell = (ST_int *)malloc((size_t)N * sizeof(ST_int));
        P->val = (ST_double *)malloc((size_t)total * sizeof(ST_double));
        if (!newptr || !mcol || !P->cell || !P->val)
        {
            free(newptr);
            free(mcol);
            rc = 920;
            goto fail;
        }
        newptr[0] = 0;
        for (ST_int r = 0; r < nrow; r++)
            newptr[r + 1] = newptr[r] + len[r];
        ST_int *cell = P->cell;
        const int tuse = hcc_stamp_threads(ncol, 2);
#pragma omp parallel num_threads(tuse) reduction(| : failed) if (par)
        {
            ST_int *stamp = (ST_int *)malloc((size_t)ncol * sizeof(ST_int));
            ST_int *slot = (ST_int *)malloc((size_t)ncol * sizeof(ST_int));
            if (!stamp || !slot)
                failed |= 1;
            else
                memset(stamp, 0xFF, (size_t)ncol * sizeof(ST_int));
#pragma omp for schedule(dynamic, 1024)
            for (ST_int r = 0; r < nrow; r++)
            {
                if (!stamp || !slot)
                    continue;
                ST_int out = newptr[r];
                for (ST_int e = rowptr[r]; e < rowptr[r + 1]; e++)
                {
                    const ST_int j = src_col[e] - 1;
                    if (stamp[j] == r)
                        cell[e] = slot[j];
                    else
                    {
                        stamp[j] = r;
                        slot[j] = out;
                        mcol[out] = j + 1;
                        cell[e] = out;
                        out++;
                    }
                }
            }
            free(stamp);
            free(slot);
        }
        if (failed)
        {
            free(newptr);
            free(mcol);
            rc = 920;
            goto fail;
        }
        free(len);
        free(P->col);
        P->col = mcol;
        P->cellptr = newptr;
        C->rowptr = newptr;
        C->col = mcol;
        C->nnz = (ST_int)total;
    }
    return 0;

fail:
    free(len);
    hdfe_crossplan_free(P);
    return rc;
}

void hdfe_crossplan_values(hdfe_crossplan *P, const ST_double *w)
{
    hdfe_crosstab *C = &P->C;
    if (!P->val)
    {
        C->val = w;
        return;
    }
    const ST_int N = P->N;
    const ST_int *RESTRICT order = P->order, *RESTRICT cell = P->cell;
    ST_double *RESTRICT val = P->val;
    const int par = N > HDFE_COEF_PARALLEL_MIN;
    if (!cell)
    {
#pragma omp parallel for schedule(static) if (par)
        for (ST_int e = 0; e < N; e++)
            val[e] = w[order[e]];
    }
    else
    {
        /* A level's cells belong to that level alone, and its rows come in
         * their original order: every cell's weights are added in the order
         * hdfe_crosstab_build() adds them, whichever thread has the level. */
        const ST_int *RESTRICT grouprow = P->grouprow, *RESTRICT cellptr = P->cellptr;
        const ST_int *chunk = P->chunk;
#pragma omp parallel for schedule(dynamic, 1) if (par && P->nchunk > 1)
        for (int c = 0; c < P->nchunk; c++)
        {
            const ST_int r0 = chunk[c], r1 = chunk[c + 1];
            for (ST_int x = cellptr[r0]; x < cellptr[r1]; x++)
                val[x] = 0.0;
            if (order)
            {
                for (ST_int e = grouprow[r0]; e < grouprow[r1]; e++)
                    val[cell[e]] += w[order[e]];
            }
            else
            {
                for (ST_int e = grouprow[r0]; e < grouprow[r1]; e++)
                    val[cell[e]] += w[e];
            }
        }
    }
    C->val = val;
}

/* ========================================================================
 * Preconditioned conjugate gradients for several right-hand sides
 * ======================================================================== */

/* s[k] = sum over j of a[j,k] * b[j,k] (* m[j] when m != NULL): fixed blocks
 * summed in order, so the result does not depend on the thread count */
static void hcc_dots(const ST_double *a, const ST_double *b, const ST_double *m, ST_int n, int K,
                     ST_double *part, ST_double *s, int par)
{
    const ST_int nblk = (n + HDFE_COEF_PCG_BLOCK - 1) / HDFE_COEF_PCG_BLOCK;
#pragma omp parallel for schedule(static) if (par)
    for (ST_int blk = 0; blk < nblk; blk++)
    {
        const ST_int j0 = blk * HDFE_COEF_PCG_BLOCK;
        const ST_int j1 = n - j0 > HDFE_COEF_PCG_BLOCK ? j0 + HDFE_COEF_PCG_BLOCK : n;
        ST_double acc[HDFE_COEF_KMAX];
        for (int k = 0; k < K; k++)
            acc[k] = 0.0;
        for (ST_int j = j0; j < j1; j++)
        {
            const ST_double *aj = a + (size_t)j * K, *bj = b + (size_t)j * K;
            const ST_double mj = m ? m[j] : 1.0;
            for (int k = 0; k < K; k++)
                acc[k] += aj[k] * bj[k] * mj;
        }
        for (int k = 0; k < K; k++)
            part[(size_t)blk * K + k] = acc[k];
    }
    for (int k = 0; k < K; k++)
    {
        ST_double total = 0.0;
        for (ST_int blk = 0; blk < nblk; blk++)
            total += part[(size_t)blk * K + k];
        s[k] = total;
    }
}

/* Whether a column whose solve ended at rounding level should be solved once
 * more from its residuals: the rounding of this round, a fraction of the
 * column's sum of squares scale, has to exceed what the tolerance allows,
 * and the next round has to start from a clearly smaller column (its sum of
 * squares will be pot, the residual sum of squares of this one). */
static int hcc_refine(ST_double tol2, ST_double f0, ST_double pot, ST_double scale)
{
    const ST_double least = HDFE_COEF_POTENTIAL_FLOOR * f0;
    const ST_double left = pot > least ? pot : least;
    return HDFE_COEF_ROUNDING * HDFE_COEF_ROUNDING * scale > tol2 * left && 4.0 * pot < scale;
}

ST_retcode hdfe_coef_pcg(const hdfe_coef_pcg_problem *P, const ST_double *g, ST_double *x,
                         const ST_double *f0, const ST_double *scale, ST_double *work,
                         unsigned char *state, ST_int *iterations)
{
    const ST_int n = P->n;
    const int K = P->K;
    const size_t len = (size_t)n * (size_t)K;
    ST_double *r = work, *p = work + len, *q = work + 2 * len, *z = work + 3 * len;
    ST_double *part = work + 4 * len;
    const ST_double *minv = P->minv;
    const ST_double tol2 = P->tolerance * P->tolerance;
    const int par = len > HDFE_COEF_PARALLEL_MIN;
    ST_double rz[HDFE_COEF_KMAX], pq[HDFE_COEF_KMAX], alpha[HDFE_COEF_KMAX], beta[HDFE_COEF_KMAX];
    ST_double pot[HDFE_COEF_KMAX], next[HDFE_COEF_KMAX], noise[HDFE_COEF_KMAX];
    int done[HDFE_COEF_KMAX];
    ST_int stopped[HDFE_COEF_KMAX];

    *iterations = 0;
    for (int k = 0; k < K; k++)
    {
        pot[k] = f0[k];
        noise[k] = HDFE_COEF_ROUNDING * HDFE_COEF_ROUNDING * scale[k];
        done[k] = state[k] != HDFE_COEF_COLUMN_ACTIVE;
        stopped[k] = 0;
    }
    if (P->warm)
    {
        /* residual and residual sum of squares at the starting point:
         * E(x) = f0 - 2 x'g + x'S x */
        ST_double xq[HDFE_COEF_KMAX], xg[HDFE_COEF_KMAX];
        P->apply(P->apply_ctx, x, q);
        hcc_dots(x, q, NULL, n, K, part, xq, par);
        hcc_dots(x, g, NULL, n, K, part, xg, par);
        for (int k = 0; k < K; k++)
            pot[k] += xq[k] - 2.0 * xg[k];
#pragma omp parallel for schedule(static) if (par)
        for (ST_int j = 0; j < n; j++)
        {
            for (int k = 0; k < K; k++)
                r[(size_t)j * K + k] = g[(size_t)j * K + k] - q[(size_t)j * K + k];
        }
    }
    else
    {
        memset(x, 0, len * sizeof(ST_double));
        memcpy(r, g, len * sizeof(ST_double));
    }

    /* p = z = M^-1 r, rz = r'z */
    if (P->precond)
    {
        P->precond(P->precond_ctx, r, z);
        hcc_dots(r, z, NULL, n, K, part, rz, par);
        memcpy(p, z, len * sizeof(ST_double));
    }
    else
    {
        hcc_dots(r, r, minv, n, K, part, rz, par);
#pragma omp parallel for schedule(static) if (par)
        for (ST_int j = 0; j < n; j++)
        {
            const ST_double mj = minv ? minv[j] : 1.0;
            for (int k = 0; k < K; k++)
                p[(size_t)j * K + k] = mj * r[(size_t)j * K + k];
        }
    }
    int active = 0;
    for (int k = 0; k < K; k++)
    {
        if (done[k])
            continue;
        if (!(rz[k] > noise[k]))
        {
            /* nothing above this round's rounding to solve for */
            done[k] = 1;
            state[k] = hcc_refine(tol2, f0[k], pot[k], scale[k]) ? HDFE_COEF_COLUMN_ROUNDING
                                                                 : HDFE_COEF_COLUMN_CONVERGED;
        }
        else
            active++;
    }
    if (P->warm || active < K)
    {
        /* columns that take no part keep their x and get no direction */
#pragma omp parallel for schedule(static) if (par)
        for (ST_int j = 0; j < n; j++)
        {
            for (int k = 0; k < K; k++)
            {
                if (done[k])
                    p[(size_t)j * K + k] = 0.0;
            }
        }
    }

    ST_int it = 0;
    int closing = 0;
    while (active > 0 && (it < P->maxiter || (closing && it - P->maxiter < P->extend)))
    {
        if (ctools_cancel_requested())
            return 1;
        it++;
        closing = 1;
        P->apply(P->apply_ctx, p, q);
        hcc_dots(p, q, NULL, n, K, part, pq, par);
        for (int k = 0; k < K; k++)
            alpha[k] = (!done[k] && pq[k] > 0.0) ? rz[k] / pq[k] : 0.0;

        /* x += alpha p, r -= alpha q */
#pragma omp parallel for schedule(static) if (par)
        for (ST_int j = 0; j < n; j++)
        {
            for (int k = 0; k < K; k++)
            {
                const size_t ix = (size_t)j * K + k;
                x[ix] += alpha[k] * p[ix];
                r[ix] -= alpha[k] * q[ix];
            }
        }
        if (P->precond)
        {
            P->precond(P->precond_ctx, r, z);
            hcc_dots(r, z, NULL, n, K, part, next, par);
        }
        else
            hcc_dots(r, r, minv, n, K, part, next, par);

        for (int k = 0; k < K; k++)
        {
            beta[k] = 0.0;
            if (done[k])
                continue;
            /* the step reduced the squared distance to the projection by
             * alpha * rz; stop when that is at most tol^2 of what is left */
            const ST_double decrease = alpha[k] * rz[k];
            pot[k] -= decrease;
            stopped[k] = it;
            const ST_double least = HDFE_COEF_POTENTIAL_FLOOR * f0[k];
            const ST_double left = pot[k] > least ? pot[k] : least;
            if (pq[k] > 0.0 && decrease <= tol2 * left)
            {
                done[k] = 1;
                state[k] = HDFE_COEF_COLUMN_CONVERGED;
            }
            else if (!(pq[k] > 0.0) || !(decrease > noise[k]) || !(next[k] > noise[k]))
            {
                /* The residual is down to the rounding of the recurrence.
                 * Steps from here on divide noise by noise and can undo the
                 * digits already gained, so this is as far as the round
                 * goes. */
                done[k] = 1;
                state[k] = hcc_refine(tol2, f0[k], pot[k], scale[k]) ? HDFE_COEF_COLUMN_ROUNDING
                                                                     : HDFE_COEF_COLUMN_CONVERGED;
            }
            else
            {
                beta[k] = next[k] / rz[k];
                if (decrease > HDFE_COEF_NEAR * tol2 * left)
                    closing = 0;
            }
            if (done[k])
                active--;
            rz[k] = next[k];
        }

        /* p = z + beta p; a finished column's direction is zeroed */
#pragma omp parallel for schedule(static) if (par)
        for (ST_int j = 0; j < n; j++)
        {
            const ST_double mj = minv ? minv[j] : 1.0;
            for (int k = 0; k < K; k++)
            {
                const size_t ix = (size_t)j * K + k;
                const ST_double zz = P->precond ? z[ix] : mj * r[ix];
                p[ix] = done[k] ? 0.0 : zz + beta[k] * p[ix];
            }
        }
    }
    for (int k = 0; k < K; k++)
    {
        if (stopped[k] > *iterations)
            *iterations = stopped[k];
    }
    return active > 0 ? 430 : 0;
}

/* ========================================================================
 * Dense Cholesky of a positive semidefinite matrix
 * ======================================================================== */

ST_retcode hdfe_coef_chol_factor(ST_double *A, ST_int n, const unsigned char *pin,
                                 unsigned char *dead, ST_double *weakest)
{
    ST_double *d0 = (ST_double *)malloc((size_t)(n > 0 ? n : 1) * sizeof(ST_double));
    if (!d0)
        return 920;
    for (ST_int j = 0; j < n; j++)
        d0[j] = A[(size_t)j * n + j];
    ST_double low = 1.0;

    /* Column by column: the pivot on one thread, then the column below it
     * shared out. Every entry is one dot product of two row prefixes. */
#pragma omp parallel if (n >= 256)
    {
        for (ST_int j = 0; j < n; j++)
        {
            ST_double *RESTRICT Aj = A + (size_t)j * n;
#pragma omp single
            {
                ST_double acc = 0.0;
#pragma omp simd reduction(+ : acc)
                for (ST_int k = 0; k < j; k++)
                    acc += Aj[k] * Aj[k];
                const ST_double d = Aj[j] - acc;
                if ((pin && pin[j]) || !(d0[j] > 0.0) || !(d > HDFE_COEF_PIVOT_TOL * d0[j]))
                {
                    dead[j] = 1;
                    Aj[j] = 1.0;
                }
                else
                {
                    dead[j] = 0;
                    Aj[j] = sqrt(d);
                    if (d < low * d0[j])
                        low = d / d0[j];
                }
            }
            const int skip = dead[j];
            const ST_double inv = 1.0 / Aj[j];
#pragma omp for schedule(static)
            for (ST_int i = j + 1; i < n; i++)
            {
                ST_double *RESTRICT Ai = A + (size_t)i * n;
                if (skip)
                {
                    Ai[j] = 0.0;
                    continue;
                }
                ST_double s = 0.0;
#pragma omp simd reduction(+ : s)
                for (ST_int k = 0; k < j; k++)
                    s += Ai[k] * Aj[k];
                Ai[j] = (Ai[j] - s) * inv;
            }
        }
    }
    free(d0);
    *weakest = low;
    return 0;
}

void hdfe_coef_direct_state(ST_int n, int K, const ST_double *g, const ST_double *b,
                            const ST_double *f0, const ST_double *scale, ST_double tolerance,
                            ST_double weakest, unsigned char *state)
{
    const ST_double tol2 = tolerance * tolerance;
    const ST_double unit = HDFE_COEF_ROUNDING * HDFE_COEF_ROUNDING;
    for (int k = 0; k < K; k++)
    {
        /* b'g = b'S b, the part of f0 the kept effects explain */
        ST_double explained = 0.0;
        for (ST_int j = 0; j < n; j++)
            explained += b[(size_t)j * K + k] * g[(size_t)j * K + k];
        if (!(explained > 0.0))
            explained = 0.0;
        const ST_double pot = f0[k] > explained ? f0[k] - explained : 0.0;
        const ST_double least = HDFE_COEF_POTENTIAL_FLOOR * f0[k];
        const ST_double left = pot > least ? pot : least;
        /* the solve's own error is relative to what it explains and grows
         * with the weakest pivot; a second round explains almost nothing */
        const int solve = unit * explained > tol2 * left * weakest * weakest;
        state[k] = (solve || hcc_refine(tol2, f0[k], pot, scale[k])) ? HDFE_COEF_COLUMN_ROUNDING
                                                                     : HDFE_COEF_COLUMN_CONVERGED;
    }
}

int hdfe_coef_refine_span(const unsigned char *state, int K, int *first, int *count)
{
    int lo = K, hi = -1;
    for (int k = 0; k < K; k++)
    {
        if (state[k] == HDFE_COEF_COLUMN_ROUNDING)
        {
            if (k < lo)
                lo = k;
            hi = k;
        }
    }
    if (hi < 0)
        return 0;
    *first += lo;
    *count = hi - lo + 1;
    return 1;
}

void hdfe_coef_chol_solve(const ST_double *L, ST_int n, const unsigned char *dead, ST_double *x,
                          int K)
{
    for (ST_int i = 0; i < n; i++)
    {
        const ST_double *Li = L + (size_t)i * n;
        ST_double *xi = x + (size_t)i * K;
        if (dead[i])
        {
            for (int k = 0; k < K; k++)
                xi[k] = 0.0;
            continue;
        }
        for (ST_int j = 0; j < i; j++)
        {
            const ST_double l = Li[j];
            const ST_double *xj = x + (size_t)j * K;
            for (int k = 0; k < K; k++)
                xi[k] -= l * xj[k];
        }
        const ST_double inv = 1.0 / Li[i];
        for (int k = 0; k < K; k++)
            xi[k] *= inv;
    }
    for (ST_int i = n - 1; i >= 0; i--)
    {
        ST_double *xi = x + (size_t)i * K;
        if (dead[i])
            continue;
        for (ST_int j = i + 1; j < n; j++)
        {
            const ST_double l = L[(size_t)j * n + i];
            const ST_double *xj = x + (size_t)j * K;
            for (int k = 0; k < K; k++)
                xi[k] -= l * xj[k];
        }
        const ST_double inv = 1.0 / L[(size_t)i * n + i];
        for (int k = 0; k < K; k++)
            xi[k] *= inv;
    }
}

/* ========================================================================
 * Passes over the data rows
 * ======================================================================== */

/* t[level of B] += w x over the rows i0..i1-1 */
HDFE_COEF_KERNEL void hcc_rows_sum(const hdfe_coef_rows *D, const ST_double *RESTRICT X, ST_int i0,
                                   ST_int i1, ST_double *RESTRICT t, const int K)
{
    const ST_int *RESTRICT lvB = D->lvB;
    const ST_double *RESTRICT w = D->w;
    const size_t N = (size_t)D->N;
    if (w)
    {
        for (ST_int i = i0; i < i1; i++)
        {
            ST_double *tf = t + (size_t)(lvB[i] - 1) * K;
            const ST_double wi = w[i];
            for (int k = 0; k < K; k++)
                tf[k] += wi * X[(size_t)k * N + i];
        }
    }
    else
    {
        for (ST_int i = i0; i < i1; i++)
        {
            ST_double *tf = t + (size_t)(lvB[i] - 1) * K;
            for (int k = 0; k < K; k++)
                tf[k] += X[(size_t)k * N + i];
        }
    }
}

/* With v = x - a0[level of B]: g[level of each kept effect] += w v,
 * ss[k] += w v^2 and ss[K + k] += w x^2, over the rows i0..i1-1 */
HDFE_COEF_KERNEL void hcc_rows_rhs(const hdfe_coef_rows *D, const ST_double *RESTRICT X, ST_int i0,
                                   ST_int i1, const ST_double *RESTRICT a0, ST_double *RESTRICT g,
                                   ST_double *RESTRICT ss, const int nk, const int K)
{
    const ST_int *RESTRICT lvB = D->lvB;
    const ST_double *RESTRICT w = D->w;
    const size_t N = (size_t)D->N;
    ST_double wv[HDFE_COEF_KMAX];
    for (ST_int i = i0; i < i1; i++)
    {
        const ST_double *af = a0 + (size_t)(lvB[i] - 1) * K;
        const ST_double wi = w ? w[i] : 1.0;
        for (int k = 0; k < K; k++)
        {
            const ST_double x = X[(size_t)k * N + i];
            const ST_double v = x - af[k];
            wv[k] = wi * v;
            ss[k] += wv[k] * v;
            ss[K + k] += wi * x * x;
        }
        for (int h = 0; h < nk; h++)
        {
            ST_double *gj = g + (size_t)(D->off[h] + D->lvk[h][i] - 1) * K;
            for (int k = 0; k < K; k++)
                gj[k] += wv[k];
        }
    }
}

/* Values of the K columns in one stretch of rows of the fused pass below: a
 * stretch stays in the second-level cache between its two readings. */
#define HDFE_COEF_FUSED_VALUES 16384

/*
 * Rows grouped by B: hcc_rows_sum(), the division by the level counts and
 * hcc_rows_rhs() over the rows i0..i1-1 (i0 and i1 are level boundaries),
 * a stretch of whole levels at a time. A level's sum is complete when its
 * rows end, so the pass that subtracts the level means can follow at once
 * and finds the rows in the cache: the columns come from memory once
 * instead of twice. Every sum adds the same values in the same order as the
 * two separate passes.
 *
 * a0 is set here for every level from the one after row i0-1's through row
 * i1-1's (through the last one when i1 is N); levels without rows get zero.
 */
HDFE_COEF_KERNEL void hcc_rows_rhs_grouped(const hdfe_coef_rows *D, const ST_double *RESTRICT X,
                                           ST_int i0, ST_int i1, ST_double *RESTRICT a0,
                                           ST_double *RESTRICT g, ST_double *RESTRICT ss,
                                           const int nk, const int K)
{
    const ST_int *RESTRICT lvB = D->lvB;
    const ST_double *RESTRICT invB = D->invB;
    const ST_int step = HDFE_COEF_FUSED_VALUES / K > 256 ? HDFE_COEF_FUSED_VALUES / K : 256;
    /* levels below lo are done (0-based) */
    ST_int lo = i0 > 0 ? lvB[i0 - 1] : 0;
    ST_int s = i0;
    while (s < i1)
    {
        ST_int e = i1 - s > step ? s + step : i1;
        while (e < i1 && lvB[e] == lvB[e - 1])
            e++;
        const ST_int hi = lvB[e - 1];
        memset(a0 + (size_t)lo * K, 0, (size_t)(hi - lo) * K * sizeof(ST_double));
        hcc_rows_sum(D, X, s, e, a0, K);
        for (ST_int f = lo; f < hi; f++)
        {
            const ST_double ib = invB[f];
            for (int k = 0; k < K; k++)
                a0[(size_t)f * K + k] *= ib;
        }
        hcc_rows_rhs(D, X, s, e, a0, g, ss, nk, K);
        lo = hi;
        s = e;
    }
    if (i1 == D->N && lo < D->nB)
        memset(a0 + (size_t)lo * K, 0, (size_t)(D->nB - lo) * K * sizeof(ST_double));
}

/* x -= a[level of B] + the coefficients of the row's kept levels */
HDFE_COEF_KERNEL void hcc_rows_residual(const hdfe_coef_rows *D, ST_double *RESTRICT X, ST_int i0,
                                        ST_int i1, const ST_double *RESTRICT a,
                                        const ST_double *RESTRICT beta, const int nk, const int K)
{
    const ST_int *RESTRICT lvB = D->lvB;
    const size_t N = (size_t)D->N;
    ST_double s[HDFE_COEF_KMAX];
    for (ST_int i = i0; i < i1; i++)
    {
        const ST_double *af = a + (size_t)(lvB[i] - 1) * K;
        for (int k = 0; k < K; k++)
            s[k] = af[k];
        for (int h = 0; h < nk; h++)
        {
            const ST_double *bj = beta + (size_t)(D->off[h] + D->lvk[h][i] - 1) * K;
            for (int k = 0; k < K; k++)
                s[k] += bj[k];
        }
        for (int k = 0; k < K; k++)
            X[(size_t)k * N + i] -= s[k];
    }
}

void hdfe_coef_rowplan_free(hdfe_coef_rowplan *R)
{
    free(R->tbuf);
    free(R->gbuf);
    free(R->ssbuf);
    memset(R, 0, sizeof(*R));
}

/* buf holds at least bytes: kept when it already does (contents are scratch) */
static int hcc_grow(ST_double **buf, size_t *have, size_t bytes)
{
    if (*buf && *have >= bytes)
        return 1;
    free(*buf);
    *buf = (ST_double *)malloc(bytes ? bytes : 1);
    *have = *buf ? bytes : 0;
    return *buf != NULL;
}

ST_retcode hdfe_coef_rowplan_init(hdfe_coef_rowplan *R, const hdfe_coef_rows *D, int K)
{
    const ST_int N = D->N;
    const size_t per_k = (size_t)D->m * (size_t)K * sizeof(ST_double);
    const size_t per_b = (size_t)D->nB * (size_t)K * sizeof(ST_double);
    int nrc = N / HDFE_COEF_ROW_CHUNK_ROWS;
    if (nrc > HDFE_COEF_ROW_CHUNKS)
        nrc = HDFE_COEF_ROW_CHUNKS;
    if (nrc < 1)
        nrc = 1;
    while (nrc > 1 && (per_k * (size_t)nrc > HDFE_COEF_BUFFER_BYTES ||
                       (!D->sorted && per_b * (size_t)nrc > HDFE_COEF_BUFFER_BYTES)))
        nrc /= 2;
    if (nrc < 4 && N >= 4 * HDFE_COEF_ROW_CHUNK_ROWS && K > 1)
    {
        /* The buffers leave too few chunks to be worth it: one column per
         * thread, in column-major scratch. */
        R->nrc = 0;
        R->private_b = 0;
        if (!hcc_grow(&R->gbuf, &R->gbytes,
                      ((size_t)D->nB + (size_t)D->m) * (size_t)K * sizeof(ST_double)))
        {
            hdfe_coef_rowplan_free(R);
            return 920;
        }
        return 0;
    }
    R->nrc = nrc;
    R->private_b = !D->sorted && nrc > 1;
    R->bound[0] = 0;
    for (int c = 1; c < nrc; c++)
    {
        ST_int i = (ST_int)((long long)N * c / nrc);
        if (i < R->bound[c - 1])
            i = R->bound[c - 1];
        /* grouped rows: a level of B stays within one chunk */
        if (D->sorted)
        {
            while (i > 0 && i < N && D->lvB[i] == D->lvB[i - 1])
                i++;
        }
        R->bound[c] = i;
    }
    R->bound[nrc] = N;
    int ok = hcc_grow(&R->gbuf, &R->gbytes, (size_t)nrc * per_k);
    ok &= hcc_grow(&R->ssbuf, &R->ssbytes, (size_t)nrc * 2 * (size_t)K * sizeof(ST_double));
    if (R->private_b)
        ok &= hcc_grow(&R->tbuf, &R->tbytes, (size_t)nrc * per_b);
    if (!ok)
    {
        hdfe_coef_rowplan_free(R);
        return 920;
    }
    return 0;
}

void hdfe_coef_rows_rhs(const hdfe_coef_rows *D, const hdfe_coef_rowplan *R, const ST_double *X,
                        int K, ST_double *a0, ST_double *g, ST_double *f0, ST_double *scale)
{
    const ST_int N = D->N, nB = D->nB, m = D->m;
    const int nk = D->nk;
    const ST_double *invB = D->invB;
    const int par = (size_t)N * (size_t)K > HDFE_COEF_PARALLEL_MIN;

    if (R->nrc == 0)
    {
        ST_double *ca = R->gbuf, *cg = R->gbuf + (size_t)nB * (size_t)K;
#pragma omp parallel for schedule(dynamic, 1) if (par)
        for (int k = 0; k < K; k++)
        {
            const ST_double *x = X + (size_t)k * (size_t)N;
            ST_double *ak = ca + (size_t)k * nB, *gk = cg + (size_t)k * m;
            ST_double ss[2] = {0.0, 0.0};
            memset(gk, 0, (size_t)m * sizeof(ST_double));
            if (D->sorted)
                hcc_rows_rhs_grouped(D, x, 0, N, ak, gk, ss, nk, 1);
            else
            {
                memset(ak, 0, (size_t)nB * sizeof(ST_double));
                hcc_rows_sum(D, x, 0, N, ak, 1);
                for (ST_int f = 0; f < nB; f++)
                    ak[f] *= invB[f];
                hcc_rows_rhs(D, x, 0, N, ak, gk, ss, nk, 1);
            }
            f0[k] = ss[0];
            scale[k] = ss[1];
        }
#pragma omp parallel for schedule(static) if (par)
        for (ST_int f = 0; f < nB; f++)
        {
            for (int k = 0; k < K; k++)
                a0[(size_t)f * K + k] = ca[(size_t)k * nB + f];
        }
#pragma omp parallel for schedule(static) if (par)
        for (ST_int j = 0; j < m; j++)
        {
            for (int k = 0; k < K; k++)
                g[(size_t)j * K + k] = cg[(size_t)k * m + j];
        }
        return;
    }

    const int nrc = R->nrc;
    const size_t stride_b = (size_t)nB * (size_t)K, stride_k = (size_t)m * (size_t)K;
    /* level means of B */
    if (R->private_b)
    {
#pragma omp parallel for schedule(dynamic, 1) if (par)
        for (int c = 0; c < nrc; c++)
        {
            ST_double *t = R->tbuf + (size_t)c * stride_b;
            memset(t, 0, stride_b * sizeof(ST_double));
            HDFE_COEF_K_DISPATCH(K, hcc_rows_sum(D, X, R->bound[c], R->bound[c + 1], t, KK))
        }
#pragma omp parallel for schedule(static) if (par)
        for (ST_int f = 0; f < nB; f++)
        {
            for (int k = 0; k < K; k++)
            {
                ST_double total = 0.0;
                for (int c = 0; c < nrc; c++)
                    total += R->tbuf[(size_t)c * stride_b + (size_t)f * K + k];
                a0[(size_t)f * K + k] = total * invB[f];
            }
        }
    }
    else if (!D->sorted)
    {
        /* one chunk */
        memset(a0, 0, stride_b * sizeof(ST_double));
#pragma omp parallel for schedule(dynamic, 1) if (par && nrc > 1)
        for (int c = 0; c < nrc; c++)
        {
            HDFE_COEF_K_DISPATCH(K, hcc_rows_sum(D, X, R->bound[c], R->bound[c + 1], a0, KK))
        }
#pragma omp parallel for schedule(static) if (par)
        for (ST_int f = 0; f < nB; f++)
        {
            for (int k = 0; k < K; k++)
                a0[(size_t)f * K + k] *= invB[f];
        }
    }
    /* Level sums of the kept effects of the B-demeaned columns. Rows grouped
     * by B (a level of B then lies within one chunk) take the level means of
     * B in the same pass. */
    const int fused = D->sorted && !R->private_b;
#pragma omp parallel for schedule(dynamic, 1) if (par && nrc > 1)
    for (int c = 0; c < nrc; c++)
    {
        ST_double *gc = R->gbuf + (size_t)c * stride_k, *ssc = R->ssbuf + (size_t)c * 2 * K;
        memset(gc, 0, stride_k * sizeof(ST_double));
        for (int k = 0; k < 2 * K; k++)
            ssc[k] = 0.0;
        if (fused && nk == 1)
        {
            HDFE_COEF_K_DISPATCH(
                K, hcc_rows_rhs_grouped(D, X, R->bound[c], R->bound[c + 1], a0, gc, ssc, 1, KK))
        }
        else if (fused)
        {
            HDFE_COEF_K_DISPATCH(
                K, hcc_rows_rhs_grouped(D, X, R->bound[c], R->bound[c + 1], a0, gc, ssc, nk, KK))
        }
        else if (nk == 1)
        {
            HDFE_COEF_K_DISPATCH(
                K, hcc_rows_rhs(D, X, R->bound[c], R->bound[c + 1], a0, gc, ssc, 1, KK))
        }
        else
        {
            HDFE_COEF_K_DISPATCH(
                K, hcc_rows_rhs(D, X, R->bound[c], R->bound[c + 1], a0, gc, ssc, nk, KK))
        }
    }
#pragma omp parallel for schedule(static) if (par && stride_k > 4096)
    for (ST_int j = 0; j < m; j++)
    {
        for (int k = 0; k < K; k++)
        {
            ST_double total = 0.0;
            for (int c = 0; c < nrc; c++)
                total += R->gbuf[(size_t)c * stride_k + (size_t)j * K + k];
            g[(size_t)j * K + k] = total;
        }
    }
    for (int k = 0; k < K; k++)
    {
        ST_double centred = 0.0, raw = 0.0;
        for (int c = 0; c < nrc; c++)
        {
            centred += R->ssbuf[(size_t)c * 2 * K + k];
            raw += R->ssbuf[(size_t)c * 2 * K + K + k];
        }
        f0[k] = centred;
        scale[k] = raw;
    }
}

void hdfe_coef_rows_residual(const hdfe_coef_rows *D, ST_double *X, int K, const ST_double *a,
                             const ST_double *beta)
{
    const ST_int N = D->N;
    const int nk = D->nk;
    const ST_int nblk = (N + 65535) / 65536;
#pragma omp parallel for schedule(dynamic, 2) if ((size_t)N * (size_t)K > HDFE_COEF_PARALLEL_MIN)
    for (ST_int blk = 0; blk < nblk; blk++)
    {
        const ST_int i0 = blk * 65536, i1 = N - i0 > 65536 ? i0 + 65536 : N;
        if (nk == 1)
        {
            HDFE_COEF_K_DISPATCH(K, hcc_rows_residual(D, X, i0, i1, a, beta, 1, KK))
        }
        else
        {
            HDFE_COEF_K_DISPATCH(K, hcc_rows_residual(D, X, i0, i1, a, beta, nk, KK))
        }
    }
}

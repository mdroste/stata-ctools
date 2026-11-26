/*
 * hdfe_coef_lap.c
 *
 * Approximate Cholesky factorization of a graph Laplacian, the
 * preconditioner of the coefficient-space projection for two effects whose
 * graph is poorly connected (hdfe_coef.h)
 * Part of the ctools Stata plugin suite
 *
 * Kyng and Sachdeva, "Approximate Gaussian elimination for Laplacians"
 * (FOCS 2016), in the edge-by-edge form of Spielman's Laplacians.jl
 * (approxChol). Vertices are eliminated in order of (approximate) degree.
 * Eliminating a vertex with neighbors j1..jm exactly would connect every
 * pair of them; here neighbor jt is instead joined to one later neighbor,
 * drawn in proportion to its weight, with the weight that keeps the
 * expected Laplacian equal to the exact Schur complement. The factor
 * therefore has about as many entries as the graph has edges. A vertex with
 * one or two neighbors is eliminated exactly, so trees and chains, the
 * graphs on which conjugate gradients is slowest, are solved outright.
 *
 * The draws come from a generator with a fixed seed held in the call, so
 * the factor is the same on every run and thread count.
 */

#include "hdfe_coef_internal.h"

struct hdfe_lap
{
    ST_int n;
    ST_int ncol;     /* eliminated vertices, in order */
    ST_int *col;     /* ncol: the vertex eliminated at each step */
    size_t *colptr;  /* ncol + 1: its entries in rowval/fval */
    ST_int *rowval;  /* neighbor of each entry */
    ST_double *fval; /* fraction of the vertex's remaining weight on it */
    ST_double *d;    /* n: pivot of each eliminated vertex, 0 otherwise */
    size_t nnz, cap;
};

void hdfe_lap_free(hdfe_lap *L)
{
    if (!L)
        return;
    free(L->col);
    free(L->colptr);
    free(L->rowval);
    free(L->fval);
    free(L->d);
    free(L);
}

/* xorshift64: uniform on [0, 1) */
static inline ST_double hcl_uniform(uint64_t *state)
{
    uint64_t x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *state = x;
    return (ST_double)(x >> 11) * (1.0 / 9007199254740992.0);
}

static int hcl_by_row(const void *a, const void *b)
{
    const hdfe_lap_leaf *x = (const hdfe_lap_leaf *)a, *y = (const hdfe_lap_leaf *)b;
    if (x->row != y->row)
        return (x->row > y->row) - (x->row < y->row);
    return (x->entry > y->entry) - (x->entry < y->entry);
}

static int hcl_by_weight(const void *a, const void *b)
{
    const hdfe_lap_leaf *x = (const hdfe_lap_leaf *)a, *y = (const hdfe_lap_leaf *)b;
    if (x->val != y->val)
        return (x->val > y->val) - (x->val < y->val);
    return (x->row > y->row) - (x->row < y->row);
}

static void hcl_sort(hdfe_lap_leaf *c, ST_int len, int by_weight)
{
    int (*cmp)(const void *, const void *) = by_weight ? hcl_by_weight : hcl_by_row;
    if (len > 24)
    {
        qsort(c, (size_t)len, sizeof(*c), cmp);
        return;
    }
    for (ST_int i = 1; i < len; i++)
    {
        const hdfe_lap_leaf x = c[i];
        ST_int j = i - 1;
        while (j >= 0 && cmp(&c[j], &x) > 0)
        {
            c[j + 1] = c[j];
            j--;
        }
        c[j + 1] = x;
    }
}

/* First leaf after `from` whose cumulative weight reaches a uniform draw on
 * the weight that follows `from` */
static inline ST_int hcl_draw_later(const ST_double *cum, ST_int from, ST_int len, uint64_t *state)
{
    const ST_double r = hcl_uniform(state) * (cum[len - 1] - cum[from]) + cum[from];
    ST_int lo = from + 1, hi = len - 1;
    while (lo < hi)
    {
        const ST_int mid = lo + (hi - lo) / 2;
        if (cum[mid] >= r)
            hi = mid;
        else
            lo = mid + 1;
    }
    return lo;
}

ST_int hdfe_lap_sample_star(hdfe_lap_leaf *cs, ST_double *cum, ST_int len, uint64_t *state,
                            ST_int *eu, ST_int *ev, ST_double *ew)
{
    hcl_sort(cs, len, 1);
    ST_double total = 0.0;
    for (ST_int x = 0; x < len; x++)
    {
        total += cs[x].val;
        cum[x] = total;
    }
    ST_double wdeg = total, scale = 1.0;
    ST_int ne = 0;
    for (ST_int at = 0; at < len - 1; at++)
    {
        const ST_double f = cs[at].val * scale / wdeg;
        const ST_int other = hcl_draw_later(cum, at, len, state);
        eu[ne] = cs[at].row;
        ev[ne] = cs[other].row;
        ew[ne] = f * (1.0 - f) * wdeg;
        ne++;
        scale *= 1.0 - f;
        wdeg *= (1.0 - f) * (1.0 - f);
    }
    return ne;
}

/* Bucket queue on the vertex degrees (Laplacians.jl's ApproxCholPQ): degrees
 * above n share buckets, which keeps the order approximate and the queue
 * linear in size. */
typedef struct
{
    ST_int *prev, *next, *key, *lists;
    ST_int nlists, minlist, n;
} hcl_queue;

static inline ST_int hcl_bucket(ST_int key, ST_int n)
{
    return key <= n ? key : n + key / n;
}

static void hcl_move(hcl_queue *q, ST_int i, ST_int newkey, ST_int oldlist, ST_int newlist)
{
    const ST_int prev = q->prev[i], next = q->next[i];
    if (next >= 0)
        q->prev[next] = prev;
    if (prev >= 0)
        q->next[prev] = next;
    else
        q->lists[oldlist] = next;
    const ST_int head = q->lists[newlist];
    if (head >= 0)
        q->prev[head] = i;
    q->lists[newlist] = i;
    q->prev[i] = -1;
    q->next[i] = head;
    q->key[i] = newkey;
}

static void hcl_increase(hcl_queue *q, ST_int i)
{
    const ST_int oldlist = hcl_bucket(q->key[i], q->n), newlist = hcl_bucket(q->key[i] + 1, q->n);
    if (newlist != oldlist)
        hcl_move(q, i, q->key[i] + 1, oldlist, newlist);
    else
        q->key[i]++;
}

static void hcl_decrease(hcl_queue *q, ST_int i)
{
    const ST_int oldlist = hcl_bucket(q->key[i], q->n), newlist = hcl_bucket(q->key[i] - 1, q->n);
    if (newlist != oldlist)
    {
        hcl_move(q, i, q->key[i] - 1, oldlist, newlist);
        if (newlist < q->minlist)
            q->minlist = newlist;
    }
    else
        q->key[i]--;
}

static ST_int hcl_pop(hcl_queue *q)
{
    while (q->minlist < q->nlists && q->lists[q->minlist] < 0)
        q->minlist++;
    if (q->minlist >= q->nlists)
        return -1;
    const ST_int i = q->lists[q->minlist];
    const ST_int next = q->next[i];
    q->lists[q->minlist] = next;
    if (next >= 0)
        q->prev[next] = -1;
    return i;
}

static int hcl_push(hdfe_lap *L, ST_int row, ST_double f)
{
    if (L->nnz == L->cap)
    {
        const size_t cap = L->cap ? 2 * L->cap : 1024;
        ST_int *rowval = (ST_int *)malloc(cap * sizeof(ST_int));
        ST_double *fval = (ST_double *)malloc(cap * sizeof(ST_double));
        if (!rowval || !fval)
        {
            free(rowval);
            free(fval);
            return -1;
        }
        if (L->nnz)
        {
            memcpy(rowval, L->rowval, L->nnz * sizeof(ST_int));
            memcpy(fval, L->fval, L->nnz * sizeof(ST_double));
        }
        free(L->rowval);
        free(L->fval);
        L->rowval = rowval;
        L->fval = fval;
        L->cap = cap;
    }
    L->rowval[L->nnz] = row;
    L->fval[L->nnz] = f;
    L->nnz++;
    return 0;
}

/* Parallel edges added up, self loops and nonpositive weights dropped, in
 * place. Returns the new number of edges, or -1 when memory runs out. */
static ST_int hcl_merge_edges(ST_int n, ST_int m, ST_int *eu, ST_int *ev, ST_double *ew)
{
    const size_t mm = (size_t)(m > 0 ? m : 1);
    ST_int *ptr = (ST_int *)calloc((size_t)n + 2, sizeof(ST_int));
    ST_int *tv = (ST_int *)malloc(mm * sizeof(ST_int));
    ST_double *tw = (ST_double *)malloc(mm * sizeof(ST_double));
    ST_int *stamp = (ST_int *)malloc(((size_t)n + 1) * sizeof(ST_int));
    ST_int *slot = (ST_int *)malloc(((size_t)n + 1) * sizeof(ST_int));
    ST_int out = -1;
    if (!ptr || !tv || !tw || !stamp || !slot)
        goto done;
    /* edges grouped by their smaller endpoint, in input order */
    for (ST_int e = 0; e < m; e++)
        ptr[(eu[e] < ev[e] ? eu[e] : ev[e]) + 2]++;
    for (ST_int v = 0; v < n; v++)
        ptr[v + 2] += ptr[v + 1];
    for (ST_int e = 0; e < m; e++)
    {
        const ST_int a = eu[e] < ev[e] ? eu[e] : ev[e], b = eu[e] < ev[e] ? ev[e] : eu[e];
        const ST_int at = ptr[a + 1]++;
        tv[at] = b;
        tw[at] = ew[e];
    }
    for (ST_int v = 0; v < n; v++)
        stamp[v] = -1;
    out = 0;
    for (ST_int a = 0; a < n; a++)
    {
        for (ST_int at = ptr[a]; at < ptr[a + 1]; at++)
        {
            const ST_int b = tv[at];
            if (b == a || !(tw[at] > 0.0))
                continue;
            if (stamp[b] == a)
                ew[slot[b]] += tw[at];
            else
            {
                stamp[b] = a;
                slot[b] = out;
                eu[out] = a;
                ev[out] = b;
                ew[out] = tw[at];
                out++;
            }
        }
    }
done:
    free(ptr);
    free(tv);
    free(tw);
    free(stamp);
    free(slot);
    return out;
}

hdfe_lap *hdfe_lap_factor(ST_int n, ST_int m, ST_int *eu, ST_int *ev, ST_double *ew)
{
    if (n < 1 || m < 0)
        return NULL;
    m = hcl_merge_edges(n, m, eu, ev, ew);
    if (m < 0)
        return NULL;

    hdfe_lap *L = (hdfe_lap *)calloc(1, sizeof(*L));
    if (!L)
        return NULL;
    L->n = n;
    /* Every edge is two entries, one in the list of each endpoint, each
     * holding the other endpoint, the weight and its twin. An entry whose
     * weight is zero is dead. */
    const size_t entries = (size_t)2 * (size_t)(m > 0 ? m : 1);
    ST_int *row = (ST_int *)malloc(entries * sizeof(ST_int));
    ST_int *next = (ST_int *)malloc(entries * sizeof(ST_int));
    ST_int *twin = (ST_int *)malloc(entries * sizeof(ST_int));
    ST_double *val = (ST_double *)malloc(entries * sizeof(ST_double));
    ST_int *head = (ST_int *)malloc((size_t)n * sizeof(ST_int));
    hcl_queue q;
    memset(&q, 0, sizeof(q));
    q.n = n;
    q.nlists = n + 2 + (ST_int)(((long long)2 * m) / n);
    q.prev = (ST_int *)malloc((size_t)n * sizeof(ST_int));
    q.next = (ST_int *)malloc((size_t)n * sizeof(ST_int));
    q.key = (ST_int *)calloc((size_t)n, sizeof(ST_int));
    q.lists = (ST_int *)malloc((size_t)q.nlists * sizeof(ST_int));
    ST_int cap = 64;
    hdfe_lap_leaf *cs = (hdfe_lap_leaf *)malloc((size_t)cap * sizeof(*cs));
    ST_double *cum = (ST_double *)malloc((size_t)cap * sizeof(ST_double));
    L->col = (ST_int *)malloc((size_t)n * sizeof(ST_int));
    L->colptr = (size_t *)malloc(((size_t)n + 1) * sizeof(size_t));
    L->d = (ST_double *)calloc((size_t)n, sizeof(ST_double));
    int ok = row && next && twin && val && head && q.prev && q.next && q.key && q.lists && cs &&
             cum && L->col && L->colptr && L->d;
    if (!ok)
        goto fail;

    for (ST_int v = 0; v < n; v++)
        head[v] = -1;
    for (ST_int e = 0; e < m; e++)
    {
        const ST_int u = eu[e], v = ev[e], a = 2 * e, b = 2 * e + 1;
        row[a] = v;
        val[a] = ew[e];
        next[a] = head[u];
        head[u] = a;
        twin[a] = b;
        row[b] = u;
        val[b] = ew[e];
        next[b] = head[v];
        head[v] = b;
        twin[b] = a;
        q.key[u]++;
        q.key[v]++;
    }
    for (ST_int l = 0; l < q.nlists; l++)
        q.lists[l] = -1;
    for (ST_int v = 0; v < n; v++)
    {
        const ST_int l = hcl_bucket(q.key[v], n);
        const ST_int first = q.lists[l];
        q.prev[v] = -1;
        q.next[v] = first;
        if (first >= 0)
            q.prev[first] = v;
        q.lists[l] = v;
    }

    uint64_t state = UINT64_C(0x9E3779B97F4A7C15);
    ST_int ncol = 0;
    for (ST_int step = 0; step < n; step++)
    {
        const ST_int i = hcl_pop(&q);
        if (i < 0)
            break;
        /* live entries of vertex i */
        ST_int len = 0;
        for (ST_int e = head[i]; e >= 0; e = next[e])
        {
            if (!(val[e] > 0.0))
                continue;
            if (len == cap)
            {
                const ST_int grown = 2 * cap;
                hdfe_lap_leaf *cs2 = (hdfe_lap_leaf *)malloc((size_t)grown * sizeof(*cs2));
                ST_double *cum2 = (ST_double *)malloc((size_t)grown * sizeof(ST_double));
                if (!cs2 || !cum2)
                {
                    free(cs2);
                    free(cum2);
                    goto fail;
                }
                memcpy(cs2, cs, (size_t)len * sizeof(*cs));
                free(cs);
                free(cum);
                cs = cs2;
                cum = cum2;
                cap = grown;
            }
            cs[len].entry = e;
            cs[len].row = row[e];
            cs[len].val = val[e];
            len++;
        }
        if (len == 0)
            continue; /* the last vertex of its component: pivot 0, held at zero */

        /* parallel edges merge; the twin of each absorbed entry dies */
        hcl_sort(cs, len, 0);
        ST_int last = -1, current = -1;
        for (ST_int x = 0; x < len; x++)
        {
            if (cs[x].row != current)
            {
                current = cs[x].row;
                cs[++last] = cs[x];
            }
            else
            {
                cs[last].val += cs[x].val;
                val[twin[cs[x].entry]] = 0.0;
                hcl_decrease(&q, current);
            }
        }
        len = last + 1;
        hcl_sort(cs, len, 1);

        L->col[ncol] = i;
        L->colptr[ncol] = L->nnz;
        ncol++;
        ST_double total = 0.0;
        for (ST_int x = 0; x < len; x++)
        {
            total += cs[x].val;
            cum[x] = total;
        }
        ST_double wdeg = total, scale = 1.0;
        for (ST_int at = 0; at < len - 1; at++)
        {
            const ST_int e = cs[at].entry, j = cs[at].row, back = twin[e];
            const ST_double f = cs[at].val * scale / wdeg;
            const ST_int k = cs[hcl_draw_later(cum, at, len, &state)].row;
            hcl_increase(&q, k);
            const ST_double fill = f * (1.0 - f) * wdeg;
            /* j's entry for i now names k; i's entry moves to k's list and names j */
            row[back] = k;
            val[back] = fill;
            twin[back] = e;
            next[e] = head[k];
            head[k] = e;
            twin[e] = back;
            val[e] = fill;
            row[e] = j;
            scale *= 1.0 - f;
            wdeg *= (1.0 - f) * (1.0 - f);
            if (hcl_push(L, j, f))
                goto fail;
        }
        {
            const ST_int e = cs[len - 1].entry, j = cs[len - 1].row;
            hcl_decrease(&q, j);
            val[twin[e]] = 0.0;
            if (hcl_push(L, j, 1.0))
                goto fail;
            L->d[i] = cs[len - 1].val * scale;
        }
    }
    L->ncol = ncol;
    L->colptr[ncol] = L->nnz;
    free(row);
    free(next);
    free(twin);
    free(val);
    free(head);
    free(q.prev);
    free(q.next);
    free(q.key);
    free(q.lists);
    free(cs);
    free(cum);
    return L;

fail:
    free(row);
    free(next);
    free(twin);
    free(val);
    free(head);
    free(q.prev);
    free(q.next);
    free(q.key);
    free(q.lists);
    free(cs);
    free(cum);
    hdfe_lap_free(L);
    return NULL;
}

HDFE_COEF_KERNEL void hcl_solve_kernel(const hdfe_lap *L, ST_double *y, const int K)
{
    const ST_int ncol = L->ncol;
    const ST_int *RESTRICT rowval = L->rowval;
    const ST_double *RESTRICT fval = L->fval;
    ST_double yi[HDFE_COEF_KMAX];
    /* forward */
    for (ST_int c = 0; c < ncol; c++)
    {
        const size_t j0 = L->colptr[c], j1 = L->colptr[c + 1] - 1;
        ST_double *pi = y + (size_t)L->col[c] * K;
        for (int k = 0; k < K; k++)
            yi[k] = pi[k];
        for (size_t jj = j0; jj < j1; jj++)
        {
            ST_double *pj = y + (size_t)rowval[jj] * K;
            const ST_double f = fval[jj];
            for (int k = 0; k < K; k++)
            {
                pj[k] += f * yi[k];
                yi[k] *= 1.0 - f;
            }
        }
        ST_double *pl = y + (size_t)rowval[j1] * K;
        for (int k = 0; k < K; k++)
        {
            pl[k] += yi[k];
            pi[k] = yi[k];
        }
    }
    /* pivots; the vertex left over in each component stays at zero */
    for (ST_int v = 0; v < L->n; v++)
    {
        ST_double *pv = y + (size_t)v * K;
        const ST_double d = L->d[v];
        if (d > 0.0)
        {
            const ST_double inv = 1.0 / d;
            for (int k = 0; k < K; k++)
                pv[k] *= inv;
        }
        else
        {
            for (int k = 0; k < K; k++)
                pv[k] = 0.0;
        }
    }
    /* backward */
    for (ST_int c = ncol - 1; c >= 0; c--)
    {
        const size_t j0 = L->colptr[c], j1 = L->colptr[c + 1] - 1;
        ST_double *pi = y + (size_t)L->col[c] * K;
        const ST_double *pl = y + (size_t)rowval[j1] * K;
        for (int k = 0; k < K; k++)
            yi[k] = pi[k] + pl[k];
        for (size_t jj = j1; jj > j0; jj--)
        {
            const ST_double *pj = y + (size_t)rowval[jj - 1] * K;
            const ST_double f = fval[jj - 1];
            for (int k = 0; k < K; k++)
                yi[k] = (1.0 - f) * yi[k] + f * pj[k];
        }
        for (int k = 0; k < K; k++)
            pi[k] = yi[k];
    }
}

void hdfe_lap_solve(const hdfe_lap *L, ST_double *y, int K)
{
    HDFE_COEF_K_DISPATCH(K, hcl_solve_kernel(L, y, KK))
}

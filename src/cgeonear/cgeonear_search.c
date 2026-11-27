/* cgeonear's search; see cgeonear_search.h.
 *
 * For each base location, a depth-first search of the neighbors' k-d tree,
 * nearer child first. A branch is skipped when its box is farther (in a
 * straight line) than the bound the candidates found so far allow; a point
 * that passes the same test gets its exact distance. The bound is derived
 * from the k-th nearest exact distance found so far (and from within() in
 * long form), with margins for rounding, so the result is the one a
 * comparison with every neighbor gives. */
#include <math.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "cgeonear_search.h"
#include "ctools_cancel.h"
#include "ctools_sort_pairs.h"

#define CGEO_PI 3.14159265358979323846
#define CGEO_MI2KM 1.609344
#define CGEO_CHUNK 256

typedef struct
{
    double key; /* ranking distance; NaN when Vincenty failed */
    uint32_t loc;
} cgeo_cand;

/* a ranks before b: nearer, then lower identifier; failed distances last */
static inline int cand_less(cgeo_cand a, cgeo_cand b)
{
    int an = isnan(a.key), bn = isnan(b.key);
    if (an != bn)
        return bn;
    if (!an && a.key != b.key)
        return a.key < b.key;
    return a.loc < b.loc;
}

static int cand_cmp(const void *pa, const void *pb)
{
    cgeo_cand a = *(const cgeo_cand *)pa, b = *(const cgeo_cand *)pb;
    return cand_less(a, b) ? -1 : (cand_less(b, a) ? 1 : 0);
}

static void cand_sort(cgeo_cand *v, size_t n)
{
    if (n < 24)
    {
        for (size_t i = 1; i < n; i++)
        {
            cgeo_cand x = v[i];
            size_t j = i;
            for (; j > 0 && cand_less(x, v[j - 1]); j--)
                v[j] = v[j - 1];
            v[j] = x;
        }
        return;
    }
    qsort(v, n, sizeof(*v), cand_cmp);
}

/* Max-heap on rank: h[0] is the worst of the n kept */
static void heap_up(cgeo_cand *h, size_t i)
{
    cgeo_cand x = h[i];
    while (i)
    {
        size_t p = (i - 1) / 2;
        if (!cand_less(h[p], x))
            break;
        h[i] = h[p];
        i = p;
    }
    h[i] = x;
}

static void heap_down(cgeo_cand *h, size_t n)
{
    cgeo_cand x = h[0];
    size_t i = 0;
    for (;;)
    {
        size_t c = 2 * i + 1;
        if (c >= n)
            break;
        if (c + 1 < n && cand_less(h[c], h[c + 1]))
            c++;
        if (!cand_less(x, h[c]))
            break;
        h[i] = h[c];
        i = c;
    }
    h[i] = x;
}

/* The squared chord (tree units) beyond which no location can be at a
 * ranking distance of `key` or less. The chord is never longer than the arc
 * on the unit sphere, nor than the geodesic between earth-centred positions,
 * which Vincenty's result matches to well under a millimetre except near
 * antipodes (where nothing is pruned). The margins cover rounding. */
static inline double chord_bound2(const cgeo_cfg *c, double key)
{
    if (isnan(key))
        return INFINITY;
    double b;
    if (c->sphere)
    {
        double angle = c->wide ? key : key / c->r;
        b = angle * (1 + 1e-12) + 1e-12;
    }
    else
    {
        if (key > 15000)
            return INFINITY;
        b = key * (1 + 1e-9) + 1e-6;
    }
    return b * b;
}

/* Ranking distance from base b to neighbor p */
static inline double exact_key(const cgeo_cfg *c, const cgeo_point *b, const cgeo_point *p)
{
    if (c->sphere)
    {
        double theta = cgeo_haversine_rad(b->lat, b->a1, b->lon, p->lat, p->a1, p->lon);
        return c->wide ? theta : theta * c->r;
    }
    double ba, smd;
    int passes;
    int st =
        cgeo_vincenty(&c->e, b->a1, b->a2, p->a1, p->a2, p->lon - b->lon, 24, &ba, &smd, &passes);
    if (st == CGEO_VINCENTY_CONVERGED)
        return ba * (smd / 1000);
    return st == CGEO_VINCENTY_COINCIDENT ? 0 : NAN;
}

double cgeo_report_value(const cgeo_cfg *c, double key)
{
    if (isnan(key))
        return c->missing;
    /* geonear: cd * R in wide form on the sphere, then / 1.609344 */
    double v = c->sphere && c->wide ? key * c->r : key;
    if (c->miles)
        v = v / CGEO_MI2KM;
    return v;
}

static void prep_point(const cgeo_cfg *c, double latdeg, double londeg, uint32_t loc, cgeo_point *p)
{
    /* geonear: lat * _pi / 180, which Stata evaluates as lat * (_pi / 180) */
    const double d2r = CGEO_PI / 180;
    double lat = latdeg * d2r, lon = londeg * d2r;
    double cl = cos(lat), sl = sin(lat), co = cos(lon), so = sin(lon);
    p->lat = lat;
    p->lon = lon;
    p->loc = loc;
    p->pad = 0;
    if (c->sphere)
    {
        p->a1 = cgeo_cos(lat);
        p->a2 = 0;
        p->c[0] = cl * co;
        p->c[1] = cl * so;
        p->c[2] = sl;
    }
    else
    {
        cgeo_reduced(&c->e, lat, &p->a1, &p->a2);
        double e2 = c->e.f * (2 - c->e.f);
        double nu = c->e.a / sqrt(1 - e2 * sl * sl) / 1000;
        p->c[0] = nu * cl * co;
        p->c[1] = nu * cl * so;
        p->c[2] = nu * (1 - e2) * sl;
    }
}

double cgeo_rank_distance(const cgeo_cfg *c, double lat1, double lon1, double lat2, double lon2)
{
    cgeo_point b, p;
    prep_point(c, lat1, lon1, 0, &b);
    prep_point(c, lat2, lon2, 0, &p);
    return exact_key(c, &b, &p);
}

void cgeo_index_init(cgeo_index *x)
{
    cgeo_tree_init(&x->tree);
}
void cgeo_index_free(cgeo_index *x)
{
    cgeo_tree_free(&x->tree);
}

int cgeo_index_build(cgeo_index *x, const cgeo_cfg *c, const double *lat, const double *lon,
                     size_t n, int threads)
{
    cgeo_index_free(x);
    cgeo_point *pts = malloc((n ? n : 1) * sizeof(*pts));
    if (!pts)
        return 920;
#pragma omp parallel for schedule(static) num_threads(threads) if (n > 20000)
    for (size_t i = 0; i < n; i++)
        prep_point(c, lat[i], lon[i], (uint32_t)i, &pts[i]);
    /* the tree takes the points over */
    return cgeo_tree_build(&x->tree, pts, n, threads);
}

/* A location's position as prep_point places it, for the Morton order only */
static void position(const cgeo_cfg *c, double latdeg, double londeg, double out[3])
{
    const double d2r = CGEO_PI / 180;
    double lat = latdeg * d2r, lon = londeg * d2r;
    double cl = cos(lat), sl = sin(lat), co = cos(lon), so = sin(lon);
    double r = 1, rz = 1;
    if (!c->sphere)
    {
        double e2 = c->e.f * (2 - c->e.f);
        r = c->e.a / sqrt(1 - e2 * sl * sl) / 1000;
        rz = r * (1 - e2);
    }
    out[0] = r * cl * co;
    out[1] = r * cl * so;
    out[2] = rz * sl;
}

void cgeo_result_init(cgeo_result *r)
{
    memset(r, 0, sizeof(*r));
}

void cgeo_result_free(cgeo_result *r)
{
    free(r->loc);
    free(r->val);
    free(r->row_start);
    free(r->nbor);
    free(r->dist);
    cgeo_result_init(r);
}

/* Growing candidate list (long form) */
typedef struct
{
    cgeo_cand *v;
    size_t n, cap;
} cand_vec;

static int vec_push(cand_vec *v, cgeo_cand x)
{
    if (v->n == v->cap)
    {
        size_t cap = v->cap ? 2 * v->cap : 64;
        cgeo_cand *nv = realloc(v->v, cap * sizeof(*nv));
        if (!nv)
            return 1;
        v->v = nv;
        v->cap = cap;
    }
    v->v[v->n++] = x;
    return 0;
}

typedef struct
{
    size_t node, lo, hi;
    int depth;
    double d2;
} stack_entry;

/* Search for base location q. heap holds up to c->k candidates; in long
 * form every candidate that ranked among the nearest when it was found, or
 * lies within(), is also appended to vec. Returns the number of distances
 * computed, or SIZE_MAX when memory ran out. */
static size_t search_one(const cgeo_tree *t, const cgeo_cfg *c, const cgeo_point *q, uint32_t self,
                         cgeo_cand *heap, size_t *hn_out, cand_vec *vec, double within_bound2)
{
    const int longform = !c->wide;
    const size_t k = c->k;
    size_t hn = 0, evals = 0;
    double bound2 = k ? INFINITY : within_bound2;
    stack_entry st[160];
    int sp = 0;
    *hn_out = 0;
    if (!t->n)
        return 0;
    st[sp++] = (stack_entry){0, 0, t->n, 0, cgeo_box_dist2(t, 0, q->c)};
    while (sp)
    {
        stack_entry e = st[--sp];
        if (e.d2 > bound2)
            continue;
        if (e.depth == t->depth)
        {
            for (size_t i = e.lo; i < e.hi; i++)
            {
                const cgeo_point *p = &t->pts[i];
                double dx = p->c[0] - q->c[0], dy = p->c[1] - q->c[1], dz = p->c[2] - q->c[2];
                double c2 = dx * dx + dy * dy + dz * dz;
                if (c2 > bound2 || p->loc == self)
                    continue;
                cgeo_cand x = {exact_key(c, q, p), p->loc};
                evals++;
                int keep = 0;
                if (k)
                {
                    if (hn < k)
                    {
                        heap[hn++] = x;
                        heap_up(heap, hn - 1);
                        keep = 1;
                    }
                    else if (cand_less(x, heap[0]))
                    {
                        heap[0] = x;
                        heap_down(heap, hn);
                        keep = 1;
                    }
                    if (keep && hn == k)
                    {
                        double b2 = chord_bound2(c, heap[0].key);
                        bound2 = longform && within_bound2 > b2 ? within_bound2 : b2;
                    }
                }
                if (longform)
                {
                    if (!keep && x.key <= c->within)
                        keep = 1;
                    if (keep && vec_push(vec, x))
                        return SIZE_MAX;
                }
            }
            continue;
        }
        size_t mid = e.lo + (e.hi - e.lo) / 2;
        size_t n1 = 2 * e.node + 1, n2 = n1 + 1;
        double d1 = cgeo_box_dist2(t, n1, q->c), d2 = cgeo_box_dist2(t, n2, q->c);
        stack_entry a = {n1, e.lo, mid, e.depth + 1, d1}, b = {n2, mid, e.hi, e.depth + 1, d2};
        if (d1 <= d2)
        {
            if (d2 <= bound2)
                st[sp++] = b;
            st[sp++] = a;
        }
        else
        {
            if (d1 <= bound2)
                st[sp++] = a;
            st[sp++] = b;
        }
    }
    *hn_out = hn;
    return evals;
}

int cgeo_search_all(const cgeo_cfg *c, const cgeo_index *x, const double *lat, const double *lon,
                    size_t nb, const uint32_t *self, int threads, cgeo_result *out)
{
    cgeo_result_free(out);
    if (threads < 1)
        threads = 1;
    const size_t k = c->k;
    const size_t nchunks = (nb + CGEO_CHUNK - 1) / CGEO_CHUNK;
    ctools_order_pair_t *order = malloc((nb ? nb : 1) * sizeof(*order));
    size_t *counts = NULL;
    cand_vec *chunk_vec = NULL;
    int rc = 0;
    if (!order)
    {
        rc = 920;
        goto done;
    }
    /* Visit the bases along a Morton curve, so that consecutive searches
     * share the tree's nodes in cache. A base's point is made when it is
     * searched, which keeps 64 bytes per base out of memory. */
    double lo = c->sphere ? -1.0 : -c->e.a / 1000 * 1.01, hi = -lo;
#pragma omp parallel for schedule(static) num_threads(threads) if (nb > 20000)
    for (size_t i = 0; i < nb; i++)
    {
        double pos[3];
        position(c, lat[i], lon[i], pos);
        order[i].order_key = (size_t)(cgeo_morton(pos, lo, hi) >> 31);
        order[i].orig_idx = i;
    }
    if (nb)
        ctools_sort_order_pairs(order, nb);
    if (c->wide)
    {
        size_t cells = nb * k > 0 ? nb * k : 1;
        out->loc = malloc(cells * sizeof(*out->loc));
        out->val = malloc(cells * sizeof(*out->val));
        if (!out->loc || !out->val)
        {
            rc = 920;
            goto done;
        }
    }
    else
    {
        counts = calloc(nb ? nb : 1, sizeof(*counts));
        chunk_vec = calloc(nchunks ? nchunks : 1, sizeof(*chunk_vec));
        if (!counts || !chunk_vec)
        {
            rc = 920;
            goto done;
        }
    }
    double within_bound2 = c->wide ? -1 : chord_bound2(c, c->within);
    size_t total_evals = 0;
    atomic_int oom = 0;
    ctools_cancel_team team = {0};
#pragma omp parallel num_threads(threads) reduction(+ : total_evals)
    {
        cgeo_cand *heap = malloc((k ? k : 1) * sizeof(*heap));
        cand_vec vec = {0};
        if (!heap)
            oom = 1;
#pragma omp for schedule(dynamic, 1) nowait
        for (size_t ch = 0; ch < nchunks; ch++)
        {
            if (oom || ctools_cancel_requested())
                continue;
            size_t begin = ch * CGEO_CHUNK, end = begin + CGEO_CHUNK < nb ? begin + CGEO_CHUNK : nb;
            cand_vec kept = {0};
            for (size_t pos = begin; pos < end && !oom; pos++)
            {
                size_t b = order[pos].orig_idx;
                size_t hn = 0;
                cgeo_point q;
                prep_point(c, lat[b], lon[b], (uint32_t)b, &q);
                vec.n = 0;
                size_t ev = search_one(&x->tree, c, &q, self ? self[b] : CGEO_NONE, heap, &hn, &vec,
                                       within_bound2);
                if (ev == SIZE_MAX)
                {
                    oom = 1;
                    break;
                }
                total_evals += ev;
                if (c->wide)
                {
                    cand_sort(heap, hn);
                    for (size_t j = 0; j < k; j++)
                    {
                        out->loc[b * k + j] = j < hn ? heap[j].loc : CGEO_NONE;
                        out->val[b * k + j] =
                            j < hn ? cgeo_report_value(c, heap[j].key) : c->missing;
                    }
                    continue;
                }
                /* geonear keeps the first nearcount() and every one within(),
                 * which in rank order is the longer of two prefixes, then at
                 * most limit() */
                cand_sort(vec.v, vec.n);
                size_t m = k < vec.n ? k : vec.n, nw = 0;
                while (nw < vec.n && vec.v[nw].key <= c->within)
                    nw++;
                if (nw > m)
                    m = nw;
                if (m > c->limit)
                    m = c->limit;
                counts[b] = m;
                for (size_t j = 0; j < m && !oom; j++)
                    if (vec_push(&kept, vec.v[j]))
                        oom = 1;
            }
            if (!c->wide)
                chunk_vec[ch] = kept;
            else
                free(kept.v);
        }
        free(heap);
        free(vec.v);
        ctools_cancel_wait(&team);
    }
    out->evals = total_evals;
    if (ctools_cancel_requested())
    {
        rc = 1;
        goto done;
    }
    if (oom)
    {
        rc = 920;
        goto done;
    }
    if (c->wide)
        goto done;
    /* Long form: lay the pairs out in base order (ascending identifier) */
    out->row_start = malloc((nb + 1) * sizeof(*out->row_start));
    if (!out->row_start)
    {
        rc = 920;
        goto done;
    }
    size_t total = 0;
    for (size_t b = 0; b < nb; b++)
    {
        out->row_start[b] = total;
        total += counts[b];
    }
    out->row_start[nb] = total;
    out->nbor = malloc((total ? total : 1) * sizeof(*out->nbor));
    out->dist = malloc((total ? total : 1) * sizeof(*out->dist));
    if (!out->nbor || !out->dist)
    {
        rc = 920;
        goto done;
    }
#pragma omp parallel for schedule(dynamic, 4) num_threads(threads) if (nchunks > 8)
    for (size_t ch = 0; ch < nchunks; ch++)
    {
        size_t begin = ch * CGEO_CHUNK, end = begin + CGEO_CHUNK < nb ? begin + CGEO_CHUNK : nb;
        size_t src = 0;
        for (size_t pos = begin; pos < end; pos++)
        {
            size_t b = order[pos].orig_idx, dst = out->row_start[b];
            for (size_t j = 0; j < counts[b]; j++, src++)
            {
                out->nbor[dst + j] = chunk_vec[ch].v[src].loc;
                out->dist[dst + j] = cgeo_report_value(c, chunk_vec[ch].v[src].key);
            }
        }
    }
    out->nout = total;
done:
    if (chunk_vec)
        for (size_t ch = 0; ch < nchunks; ch++)
            free(chunk_vec[ch].v);
    free(chunk_vec);
    free(counts);
    free(order);
    if (rc)
    {
        size_t evals = out->evals;
        cgeo_result_free(out);
        out->evals = evals;
    }
    return rc;
}

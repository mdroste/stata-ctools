/* Balanced k-d tree over neighbor locations; see cgeonear_tree.h. */
#include <stdlib.h>
#include <string.h>
#include "cgeonear_tree.h"

void cgeo_tree_init(cgeo_tree *t)
{
    memset(t, 0, sizeof(*t));
}

void cgeo_tree_free(cgeo_tree *t)
{
    free(t->pts);
    free(t->box);
    cgeo_tree_init(t);
}

static void point_swap(cgeo_point *a, cgeo_point *b)
{
    cgeo_point t = *a;
    *a = *b;
    *b = t;
}

/* Place the element of rank k (by coordinate dim) at a[k], smaller or equal
 * ones before it and larger or equal ones after it. Hoare partitioning
 * splits runs of equal keys evenly, so co-located points cost no more. */
static void select_rank(cgeo_point *a, size_t n, size_t k, int dim)
{
    size_t lo = 0, hi = n - 1;
    while (hi > lo)
    {
        if (hi - lo < 16)
        {
            for (size_t i = lo + 1; i <= hi; i++)
                for (size_t j = i; j > lo && a[j - 1].c[dim] > a[j].c[dim]; j--)
                    point_swap(&a[j - 1], &a[j]);
            return;
        }
        size_t m = lo + (hi - lo) / 2;
        double x = a[lo].c[dim], y = a[m].c[dim], z = a[hi].c[dim];
        double pivot = x < y ? (y < z ? y : (x < z ? z : x)) : (x < z ? x : (y < z ? z : y));
        size_t i = lo, j = hi;
        for (;;)
        {
            while (a[i].c[dim] < pivot)
                i++;
            while (a[j].c[dim] > pivot)
                j--;
            if (i >= j)
                break;
            point_swap(&a[i], &a[j]);
            i++;
            j--;
        }
        /* a[lo..j] <= pivot <= a[j+1..hi] */
        if (k <= j)
            hi = j;
        else
            lo = j + 1;
    }
}

static void bounding_box(const cgeo_point *a, size_t n, double *b)
{
    for (int k = 0; k < 3; k++)
    {
        b[k] = a[0].c[k];
        b[k + 3] = a[0].c[k];
    }
    for (size_t i = 1; i < n; i++)
        for (int k = 0; k < 3; k++)
        {
            double v = a[i].c[k];
            if (v < b[k])
                b[k] = v;
            if (v > b[k + 3])
                b[k + 3] = v;
        }
}

static void build_node(cgeo_tree *t, size_t node, size_t lo, size_t hi, int depth)
{
    double *b = t->box + 6 * node;
    bounding_box(t->pts + lo, hi - lo, b);
    if (depth == t->depth)
        return;
    int dim = 0;
    double ext = b[3] - b[0];
    for (int k = 1; k < 3; k++)
        if (b[k + 3] - b[k] > ext)
        {
            ext = b[k + 3] - b[k];
            dim = k;
        }
    size_t mid = lo + (hi - lo) / 2;
    select_rank(t->pts + lo, hi - lo, mid - lo, dim);
    if (hi - lo > 65536)
    {
#pragma omp task default(shared) firstprivate(node, lo, mid, depth)
        build_node(t, 2 * node + 1, lo, mid, depth + 1);
        build_node(t, 2 * node + 2, mid, hi, depth + 1);
#pragma omp taskwait
    }
    else
    {
        build_node(t, 2 * node + 1, lo, mid, depth + 1);
        build_node(t, 2 * node + 2, mid, hi, depth + 1);
    }
}

int cgeo_tree_build(cgeo_tree *t, cgeo_point *pts, size_t n, int threads)
{
    cgeo_tree_free(t);
    t->pts = pts;
    if (!n)
        return 0;
    int depth = 0;
    while (((n + ((size_t)1 << depth) - 1) >> depth) > CGEO_TREE_LEAF)
        depth++;
    size_t nnodes = ((size_t)2 << depth) - 1;
    t->box = malloc(nnodes * 6 * sizeof(*t->box));
    if (!t->pts || !t->box)
    {
        cgeo_tree_free(t);
        return 920;
    }
    t->n = n;
    t->nnodes = nnodes;
    t->depth = depth;
    if (threads < 1)
        threads = 1;
#pragma omp parallel num_threads(threads) if (n > 65536 && threads > 1)
    {
#pragma omp single
        build_node(t, 0, 0, n, 0);
    }
    return 0;
}

static uint64_t spread21(uint64_t v)
{
    v &= 0x1fffff;
    v = (v | v << 32) & 0x1f00000000ffffULL;
    v = (v | v << 16) & 0x1f0000ff0000ffULL;
    v = (v | v << 8) & 0x100f00f00f00f00fULL;
    v = (v | v << 4) & 0x10c30c30c30c30c3ULL;
    v = (v | v << 2) & 0x1249249249249249ULL;
    return v;
}

uint64_t cgeo_morton(const double c[3], double lo, double hi)
{
    double scale = hi > lo ? 2097151.0 / (hi - lo) : 0;
    uint64_t key = 0;
    for (int k = 0; k < 3; k++)
    {
        double v = (c[k] - lo) * scale;
        uint64_t q = v <= 0 ? 0 : (v >= 2097151.0 ? 2097151 : (uint64_t)v);
        key |= spread21(q) << k;
    }
    return key;
}

/* A balanced k-d tree over neighbor locations in three dimensions.
 *
 * Points are positions on the unit sphere (great-circle distances) or on the
 * reference ellipsoid in earth-centred coordinates, in km (Vincenty). The
 * straight-line (chord) distance between two positions never exceeds the
 * length of the shortest path between them along the surface, so a chord
 * bound derived from a surface distance keeps every location that can be
 * that near: searches prune on chords and rank on the exact distances.
 *
 * Layout: the tree is implicit. Node 0 holds points [0, n); a node holding
 * [lo, hi) below the leaf depth has children 2i+1 with [lo, mid) and 2i+2
 * with [mid, hi), mid = lo + (hi - lo) / 2. Every leaf is at depth `depth`.
 * The points are stored in tree order; box[6 * node] holds a node's bounding
 * box (minimum x, y, z, then maximum x, y, z). */
#ifndef CGEONEAR_TREE_H
#define CGEONEAR_TREE_H

#include <stddef.h>
#include <stdint.h>

/* A neighbor location in tree order: its position, what the exact distance
 * needs (latitude and longitude in radians; cos(lat) on the sphere, the sine
 * and cosine of the reduced latitude on the ellipsoid) and its location
 * number (the rank of its identifier, zero based). 64 bytes. */
typedef struct
{
    double c[3];
    double lat, lon;
    double a1, a2;
    uint32_t loc;
    uint32_t pad;
} cgeo_point;

typedef struct
{
    cgeo_point *pts; /* owned; n points in tree order */
    double *box;     /* owned; 6 doubles per node */
    size_t n;
    size_t nnodes;
    int depth;
} cgeo_tree;

#define CGEO_TREE_LEAF 16

void cgeo_tree_init(cgeo_tree *t);
void cgeo_tree_free(cgeo_tree *t);

/* Build from n points, which the tree takes over and reorders (pts must come
 * from malloc; the tree frees it, also when the build fails). Returns 0, or
 * 920 when memory runs out. Uses up to `threads` OpenMP threads. */
int cgeo_tree_build(cgeo_tree *t, cgeo_point *pts, size_t n, int threads);

/* Squared distance from q to the box of node `node` (0 inside it) */
static inline double cgeo_box_dist2(const cgeo_tree *t, size_t node, const double q[3])
{
    const double *b = t->box + 6 * node;
    double d2 = 0;
    for (int k = 0; k < 3; k++)
    {
        double lo = b[k] - q[k], hi = q[k] - b[k + 3];
        double d = lo > 0 ? lo : (hi > 0 ? hi : 0);
        d2 += d * d;
    }
    return d2;
}

/* Morton (Z-order) key of a position, for visiting queries in an order that
 * keeps the tree's nodes in cache. lo/hi bound the coordinates. */
uint64_t cgeo_morton(const double c[3], double lo, double hi);

#endif

/* cgeonear's search: the nearest neighbors of every base location.
 *
 * No Stata Plugin Interface calls, so that the native tests can compare it
 * with a comparison of every base with every neighbor. Locations are given
 * in degrees, numbered by the rank of their identifiers; ties in distance
 * go to the lower number (geonear's documented order: ascending nborid). */
#ifndef CGEONEAR_SEARCH_H
#define CGEONEAR_SEARCH_H

#include <stddef.h>
#include <stdint.h>
#include "cgeo_math.h"
#include "cgeonear_tree.h"

#define CGEO_NONE UINT32_MAX

typedef struct
{
    int wide;     /* wide form (else long) */
    int sphere;   /* great-circle distances (else the ellipsoid) */
    int miles;    /* report miles */
    size_t k;     /* nearcount() */
    size_t limit; /* long form: limit(), SIZE_MAX when not given */
    double r;     /* sphere radius, km */
    cgeo_ellipsoid e;
    double within;  /* long form: km (0 by default, which keeps distance 0) */
    double missing; /* the value written for a distance that failed */
} cgeo_cfg;

/* The neighbor locations in a k-d tree, with what the distances need */
typedef struct
{
    cgeo_tree tree;
} cgeo_index;

void cgeo_index_init(cgeo_index *x);
void cgeo_index_free(cgeo_index *x);
/* Returns 0 or 920 (memory) */
int cgeo_index_build(cgeo_index *x, const cgeo_cfg *c, const double *lat, const double *lon,
                     size_t n, int threads);

/* Results. Wide form: k entries per base location b at b * k: the neighbor
 * (CGEO_NONE when there are fewer) and its distance in the reported units.
 * Long form: base b's pairs are rows row_start[b] .. row_start[b + 1] - 1,
 * nearest first. evals counts the distances computed. Owns every array. */
typedef struct
{
    uint32_t *loc;
    double *val;
    size_t *row_start;
    uint32_t *nbor;
    double *dist;
    size_t nout;
    size_t evals;
} cgeo_result;

void cgeo_result_init(cgeo_result *r);
void cgeo_result_free(cgeo_result *r);

/* Search every base location. self[b], when self is not NULL, is the
 * neighbor that ignoreself leaves out for base b (or CGEO_NONE). Returns 0,
 * 920 (memory), or 1 when Break was pressed (ctools_cancel_requested()). */
int cgeo_search_all(const cgeo_cfg *c, const cgeo_index *x, const double *lat, const double *lon,
                    size_t nb, const uint32_t *self, int threads, cgeo_result *out);

/* The ranking distance of base (lat1, lon1) to neighbor (lat2, lon2), all in
 * degrees, as the search computes it: geonear's haversine angle in wide form
 * on the sphere, else its distance in km; NaN when Vincenty failed. */
double cgeo_rank_distance(const cgeo_cfg *c, double lat1, double lon1, double lat2, double lon2);
/* The value geonear stores for a ranking distance */
double cgeo_report_value(const cgeo_cfg *c, double key);

#endif

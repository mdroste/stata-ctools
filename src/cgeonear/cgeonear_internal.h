/* cgeonear internals shared across the module's files. */
#ifndef CGEONEAR_INTERNAL_H
#define CGEONEAR_INTERNAL_H

#include <stddef.h>
#include <stdint.h>
#include "stplugin.h"
#include "cgeonear_search.h" /* CGEO_NONE */

/* Locations as geonear_prep_locations makes them: one per identifier that
 * has at least one observation with nonmissing identifier, latitude and
 * longitude, in ascending order of the identifier (Stata's sort order), so
 * that a location's number is the rank of its identifier. Coordinates are
 * in degrees as stored. Owns every array. */
typedef struct
{
    size_t n;
    int id_is_str;
    double *id_num; /* n identifiers (numeric) */
    char **id_str;  /* n identifiers (string), pointers into id_buf */
    char *id_buf;
    double *lat, *lon;
} cgeo_locs;

void cgeo_locs_init(cgeo_locs *l);
void cgeo_locs_free(cgeo_locs *l);

/* Why a location table could not be made; the ado prints geonear's message */
typedef enum
{
    CGEO_FAIL_NONE = 0,
    CGEO_FAIL_NOTCONST, /* r(459) lat or lon not constant within id group */
    CGEO_FAIL_NOOBS,    /* r(459) no observation for ... locations */
    CGEO_FAIL_STRING,   /* r(7) lat/lon variables not numeric */
    CGEO_FAIL_LAT,      /* r(198) latitude must be between -90 and 90 */
    CGEO_FAIL_LON       /* r(198) longitude must be between -180 and 180 */
} cgeo_fail;

/* Read plugin variables 1-3 (identifier, latitude, longitude) of every
 * observation and make the locations. When obs_loc is not NULL it receives
 * nobs entries: the location of each observation's identifier, or
 * CGEO_NONE (missing identifier, or one without coordinates). *in_order is
 * 1 when the observations are already in identifier order (Stata's sort
 * order). On failure *fail says which of geonear's checks failed; returns a
 * Stata return code. */
int cgeo_locs_read(cgeo_locs *l, uint32_t **obs_loc, size_t *nobs, int *in_order, cgeo_fail *fail);

/* Same identifier (both numeric or both string) */
int cgeo_same_id(const cgeo_locs *a, size_t i, const cgeo_locs *b, size_t j);

/* Compare identifiers in Stata's sort order: < 0, 0 or > 0 */
int cgeo_cmp_id(const cgeo_locs *a, size_t i, const cgeo_locs *b, size_t j);

/* The distance form, "cgeonear dist vars|scalar ..." (cgeonear_dist.c) */
ST_retcode cgeonear_dist_main(const char *args);

#endif

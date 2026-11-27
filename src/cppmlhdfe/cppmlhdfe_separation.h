/*
 * cppmlhdfe_separation.h
 *
 * Separation detection for PPML estimation
 * Part of the ctools Stata plugin suite
 */

#ifndef CPPMLHDFE_SEPARATION_H
#define CPPMLHDFE_SEPARATION_H

#include "cppmlhdfe_types.h"

/* Select on original rows until singleton and all-zero FE removals stabilize.
 * levels are 1-based, mask is initialized by the caller; -1 means allocation
 * or invalid-level failure. No compact/original index spaces are mixed. */
ST_int ppml_select_fe_sample(const ST_double *y, ST_int *const *levels,
    const ST_int *num_levels, ST_int G, ST_int N, ST_int *mask,
    ST_int *num_singletons, ST_int *num_separated);

ST_int ppml_select_fe_sample_ex(const ST_double *y, ST_int *const *levels,
    const ST_int *num_levels, ST_int G, ST_int N, ST_int *mask,
    ST_int *num_singletons, ST_int *num_separated,
    const ST_int *slope_only, ST_int keep_singletons, ST_int check_separation);

/* The same, handing out its counts when rows_left is not NULL: on success
 * (*rows_left)[g][2 * l] is the number of rows level l of factor g has left
 * (the odd entries are the selection's own), G arrays of 2 * (num_levels[g]
 * + 1) entries and the array of them, all the caller's to free; NULL after a
 * failure. This is the layout HDFE_Singletons.rows_left takes
 * (hdfe_singletons_counted). */
ST_int ppml_select_fe_sample_counted(const ST_double *y, ST_int *const *levels,
    const ST_int *num_levels, ST_int G, ST_int N, ST_int *mask,
    ST_int *num_singletons, ST_int *num_separated,
    const ST_int *slope_only, ST_int keep_singletons, ST_int check_separation,
    ST_int ***rows_left);

/* Whether ppml_select_fe_sample_ex() would keep every row of a sample whose
 * mask is all ones, found in one pass: no factor has a row alone in a level
 * (num_single[g], the load's count of such rows, is 0 for every factor;
 * negative means not known; not looked at with keep_singletons), and with
 * check_separation every level of every factor that is not slope-only has a
 * positive outcome. Returns 1 when so; 0 when a row would be dropped or this
 * cannot tell (the selection then decides). */
int ppml_fe_sample_is_whole(const ST_double *y, ST_int *const *levels,
    const ST_int *num_levels, const ST_int *num_single, ST_int G, ST_int N,
    const ST_int *slope_only, ST_int keep_singletons, ST_int check_separation);

#endif /* CPPMLHDFE_SEPARATION_H */

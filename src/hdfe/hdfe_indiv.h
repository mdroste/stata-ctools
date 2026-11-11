/*
 * hdfe_indiv.h
 *
 * Individual fixed effects in a group setting, as reghdfe's group(),
 * individual() and aggregation() fit them: the regression has one row per
 * group (a patent), and an absorbed individual effect (its inventors)
 * enters each group's row as the mean (or the sum) of the effects of the
 * individuals in the group.
 * Part of the ctools Stata plugin suite
 *
 * Rows are the regression rows (one per group). The edges list, for each
 * row, the individuals of its group; an individual may be in many groups.
 * reghdfe solves the projection on the absorbed effects with LSMR, since
 * the individual effect has no per-level closed form; so does
 * hdfe_indiv_partial_out(), with the standard factors of an HDFE_State
 * (intercepts only) and the individual effect in one design.
 */

#ifndef HDFE_INDIV_H
#define HDFE_INDIV_H

#include "../ctools_hdfe_utils.h"

typedef struct {
    int mean;      /* aggregation(mean) 1, aggregation(sum) 0 */
    ST_int N_rows; /* rows the edges refer to */
    ST_int NI;     /* individuals (0-based indices below NI) */
    ST_int *ptr;   /* N_rows + 1 offsets into ind */
    ST_int *ind;   /* the individual of each edge, grouped by row */
} HDFE_Indiv;

/* The edges from three plugin variables (1-based plugin variable indices)
 * read over the whole in-range: the group number (1..#groups), the
 * individual number (1..#individuals) and the individual-level sample
 * marker. row_group[r] is the group number of regression row r (r <
 * N_rows); an edge whose group has no regression row is left out. Returns
 * 0, 920 (memory), 198 (a group number on two rows, or out of range) or the
 * loader's return code. */
ST_retcode hdfe_indiv_load(HDFE_Indiv *I, const ST_double *row_group, ST_int N_rows, int gvar,
                           int ivar, int tvar, int mean);

/* reghdfe's singletons with individual effects, to their fixed point: a row
 * alone in a level of a standard factor (levels[g][r], 1-based, g < G), and
 * a group with an individual that is in no other remaining group (the
 * 1-core of the group-individual graph with every group counted twice, as
 * reghdfe's Individual_Factor::drop_singletons computes it). mask (N_rows)
 * is set to 1 for the rows kept. Returns the rows dropped, or -1 (memory). */
ST_int hdfe_indiv_singletons(ST_int *const *levels, ST_int G, const HDFE_Indiv *I, ST_int *mask);

/* Keep the rows with mask[r] != 0 (in order) and number the individuals that
 * remain in them 0..NI-1. Returns 0 or 920. */
ST_retcode hdfe_indiv_compact(HDFE_Indiv *I, const ST_int *mask);

/* Replace each of the K columns of data (N x K, column-major, N = I->N_rows
 * = S->N) by its residual from the weighted least-squares projection on the
 * intercepts of S's factors and the individual effect (S->weights, or
 * unweighted). LSMR with diagonal preconditioning, stopped by its tests
 * with atol = btol = tol, as reghdfe's lsmr(). iterations: the most any
 * column took. Returns 0, 1 (cancelled), 920 or 430 (no convergence in maxiter). */
ST_retcode hdfe_indiv_partial_out(const HDFE_State *S, const HDFE_Indiv *I, ST_double *data, ST_int N,
                                  ST_int K, ST_double tol, ST_int maxiter, ST_int *iterations);

void hdfe_indiv_free(HDFE_Indiv *I);

#endif

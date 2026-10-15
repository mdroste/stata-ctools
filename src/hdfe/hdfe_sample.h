/*
 * hdfe_sample.h
 *
 * Estimation-sample construction shared by the HDFE estimators: the absorbed
 * factors of an HDFE_State on the retained rows, and row compaction.
 * Part of the ctools Stata plugin suite
 */

#ifndef CTOOLS_HDFE_SAMPLE_H
#define CTOOLS_HDFE_SAMPLE_H

#include "hdfe_types.h"

/* One absorbed factor on the loaded rows, before the estimation sample is
 * selected. */
typedef struct
{
    ST_int *levels;    /* Level of each loaded row (1-based) */
    ST_int *counts;    /* Rows per level */
    ST_int num_levels; /* Number of distinct levels */
    ST_int num_obs;    /* Number of loaded rows */
} HDFE_LoadedFactor;

/*
 * Remap G raw fixed-effect columns to 1-based levels with counts
 * (remap_and_count), one factor per thread. values[g] holds the N loaded
 * values of factor g, without missing values. factors must be zero-filled on
 * entry; its levels and counts are allocated here and stay owned by the
 * caller, also after a failure.
 *
 * Returns 0, or 920 when an allocation fails or a factor has no levels.
 */
ST_retcode hdfe_load_factors(HDFE_LoadedFactor *factors, ST_int G, const ST_double *const *values,
                             ST_int N);

/*
 * Build S->factors for the estimation sample.
 *
 *   G          number of absorbed factors
 *   levels     levels[g][i], i < N_in: 1-based level code of factor g on the
 *              loaded rows (as produced by remap_and_count)
 *   level_max  level_max[g]: largest code that occurs in levels[g]
 *   N_in       number of loaded rows
 *   keep       N_in flags, nonzero for a row of the estimation sample;
 *              NULL keeps every row
 *   N          number of kept rows
 *   weights    N_in observation weights aligned with levels, or NULL
 *   weighted   nonzero allocates weighted_counts even when weights is NULL
 *              (zero-filled, for callers that set them later)
 *
 * The surviving levels of each factor are renumbered 1..num_levels in order
 * of first appearance among the kept rows. counts holds the kept rows per
 * level; weighted_counts accumulates weights in row order. Every factor has
 * has_intercept = 1 and max_level = num_levels; inv_counts and the solver
 * buffers are left to ctools_hdfe_alloc_buffers().
 *
 * Sets S->G, S->N and S->factors, which must be unset on entry. Returns 0, 198
 * for an invalid argument or level code, or 920 on allocation failure; after
 * a failure S holds no factors.
 */
ST_retcode hdfe_build_factors(HDFE_State *S, ST_int G, ST_int *const *levels,
                              const ST_int *level_max, ST_int N_in, const ST_int *keep, ST_int N,
                              const ST_double *weights, ST_int weighted);

/*
 * Copy the kept rows of K columns. src is N_in x K and dst is N x K, both
 * column-major and not overlapping; keep is as in hdfe_build_factors() and
 * must flag exactly N rows. Row order is preserved.
 */
void hdfe_compact_columns(const ST_double *src, ST_int N_in, ST_int K, const ST_int *keep,
                          ST_double *dst, ST_int N);

#endif /* CTOOLS_HDFE_SAMPLE_H */

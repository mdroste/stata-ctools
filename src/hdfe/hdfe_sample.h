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
    ST_int *counts;    /* Unused, NULL (hdfe_build_factors counts the kept rows) */
    ST_int num_levels; /* Number of distinct levels */
    ST_int num_obs;    /* Number of loaded rows */
} HDFE_LoadedFactor;

/*
 * Remap G raw fixed-effect columns to 1-based levels (hdfe_remap_levels),
 * each in parallel over its rows (one factor per thread when the rows are
 * few). The numbering is arbitrary but depends on the values alone; the
 * consumers (singleton removal, hdfe_build_factors) do not depend on it.
 * values[g] holds the N loaded values of factor g, without missing values.
 * factors must be zero-filled on entry; its levels are allocated here and
 * stay owned by the caller, also after a failure.
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
 * Outside parallel regions large samples are processed in parallel row
 * blocks with the same result (weighted counts are still summed in row
 * order, one factor per thread).
 *
 * Sets S->G, S->N and S->factors, which must be unset on entry. Returns 0, 198
 * for an invalid argument or level code, or 920 on allocation failure; after
 * a failure S holds no factors.
 */
ST_retcode hdfe_build_factors(HDFE_State *S, ST_int G, ST_int *const *levels,
                              const ST_int *level_max, ST_int N_in, const ST_int *keep, ST_int N,
                              const ST_double *weights, ST_int weighted);

/* ========================================================================
 * The same three steps (load, singleton removal, factors of the sample),
 * sharing what each learns: the load numbers the levels in order of first
 * appearance and notes the rows alone in their level, so that singleton
 * removal starts from those rows, and a sample that loses no row (or few)
 * keeps the loaded levels instead of renumbering them.
 * ======================================================================== */

#include "hdfe_remap.h"

/* What the load knows about a factor besides its levels */
typedef struct
{
    ST_int ordered;       /* Levels are numbered in order of first appearance */
    HDFE_RemapOrder rows; /* First rows and single rows, when recorded */
} HDFE_FactorOrder;

/*
 * hdfe_load_factors() with the levels numbered in order of first appearance
 * where that is possible (order[g].ordered), and order[g].rows filled in when
 * the rows were loaded in parallel blocks. order holds G entries; free them
 * with hdfe_factor_order_free().
 */
ST_retcode hdfe_load_factors_ordered(HDFE_LoadedFactor *factors, HDFE_FactorOrder *order, ST_int G,
                                     const ST_double *const *values, ST_int N);
void hdfe_factor_order_free(HDFE_FactorOrder *order, ST_int G);

/* The levels of one factor that lost a row to singleton removal */
typedef struct
{
    ST_int n;      /* Number of such levels */
    ST_int *level; /* Their (loaded) numbers */
    ST_int *left;  /* Rows each has left in the sample */
} HDFE_AffectedLevels;

typedef struct
{
    ST_int num_dropped;            /* Rows dropped */
    ST_int *keep;                  /* N flags, nonzero for a kept row; NULL when none is dropped */
    ST_int listed;                 /* dropped and affected are filled in */
    ST_int *dropped;               /* The rows dropped, ascending */
    HDFE_AffectedLevels *affected; /* Per factor */
    ST_int G;
} HDFE_Singletons;

/*
 * Singleton removal to the fixed point, the sample of
 * ctools_remove_singletons() (same fweights rule). levels[g] are 1-based
 * with num_levels[g] the largest; order is hdfe_load_factors_ordered()'s, or
 * NULL. When the load recorded the rows alone in their level and they are
 * few, only the levels those rows touch are counted (one parallel pass per
 * factor and round); otherwise, or when the removals spread, every level is
 * counted as before. out is zero-filled on entry; release it with
 * hdfe_singletons_free(). Returns 0, or 920.
 */
ST_retcode hdfe_drop_singletons(HDFE_Singletons *out, ST_int *const *levels,
                                const ST_int *num_levels, const HDFE_FactorOrder *order, ST_int G,
                                ST_int N, const ST_double *fweights, ST_int verbose);
void hdfe_singletons_free(HDFE_Singletons *sing);

/*
 * The same record for keep flags made elsewhere (a sample selection of the
 * caller's own): num_dropped of the N flags are zero. When those rows are
 * few they are listed with the levels they leave (out->listed), so that
 * hdfe_build_factors_ordered() can shift the levels; when they are many, or
 * touch too many levels, only num_dropped and keep are set and the build
 * renumbers. out->keep is the caller's array, borrowed: set it to NULL
 * before hdfe_singletons_free() unless out is to free it. Returns 0, 198 or
 * 920 (out then zero).
 */
ST_retcode hdfe_singletons_from_keep(HDFE_Singletons *out, ST_int *const *levels,
                                     const ST_int *num_levels, ST_int G, ST_int N, ST_int *keep,
                                     ST_int num_dropped);

/*
 * hdfe_build_factors() from the loaded factors, with the same result. When
 * every factor is ordered and no row was dropped, the loaded level arrays
 * become the factors' (factors[g].levels is then NULL on return); when the
 * dropped rows are listed, the levels are shifted in one pass instead of
 * renumbered from scratch. sing may be NULL (no row dropped).
 *
 * With counts_pending non-NULL the per-level counts may be left for the
 * caller: *counts_pending is then nonzero, the count arrays are allocated
 * and zero, and the caller fills them before they are read, either with
 * hdfe_count_factors() or job by job (below), at the same time as other
 * work if it likes. keep and weights of those calls are sing->keep and the
 * weights given here.
 */
ST_retcode hdfe_build_factors_ordered(HDFE_State *S, ST_int G, HDFE_LoadedFactor *factors,
                                      const HDFE_FactorOrder *order, const HDFE_Singletons *sing,
                                      ST_int N_in, const ST_double *weights, ST_int weighted,
                                      ST_int *counts_pending);

/* Every factor's rows per level, and weights per level summed in row order.
 * order is the load's (it says which factors arrive in level order), or NULL */
void hdfe_count_factors(HDFE_State *S, const HDFE_FactorOrder *order, ST_int N_in,
                        const ST_int *keep, const ST_double *weights);

/*
 * The same job by job, for a caller that runs the serial passes on other
 * threads while it works. hdfe_count_rows() fills the rows per level of
 * factor g from the main thread, in parallel row blocks where it can; local
 * is order[g].rows.local, or 0. hdfe_count_rows_parallel() says whether it
 * can (few levels, or rows in level order); when not, it is one serial pass
 * and hdfe_count_rows_serial() does the same on any thread.
 * hdfe_count_weights() sums the weights of factor g's levels in row order,
 * serially, on any thread; it does nothing without weights. Each count
 * array is filled by one of these once.
 */
int hdfe_count_rows_parallel(const HDFE_State *S, ST_int g, ST_int local);
void hdfe_count_rows(HDFE_State *S, ST_int g, ST_int local);
void hdfe_count_rows_serial(HDFE_State *S, ST_int g);
void hdfe_count_weights(HDFE_State *S, ST_int g, ST_int N_in, const ST_int *keep,
                        const ST_double *weights);

/*
 * Copy the kept rows of K columns. src is N_in x K and dst is N x K, both
 * column-major and not overlapping; keep is as in hdfe_build_factors() and
 * must flag exactly N rows. Row order is preserved.
 */
void hdfe_compact_columns(const ST_double *src, ST_int N_in, ST_int K, const ST_int *keep,
                          ST_double *dst, ST_int N);

#endif /* CTOOLS_HDFE_SAMPLE_H */

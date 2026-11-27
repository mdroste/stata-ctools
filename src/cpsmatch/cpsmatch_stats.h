#ifndef CPSMATCH_STATS_H
#define CPSMATCH_STATS_H

#include <stddef.h>
#include <stdint.h>
#include "stplugin.h"

/* Inputs are borrowed loaded columns in filtered-row order. treatment is the
 * stored byte _treated column, not the original (possibly nonbinary) input.
 * NULL outcome requests support counts only. NULL weight/support describe the
 * no-treated early return: zero weights and support == 1 for every loaded row.
 * No allocations or host callbacks occur while computing the statistics. */
typedef struct
{
    double y1, y0, n1, n0, var1, var0, weight_squared;
    double y1_u, y0_u, n1_u, n0_u, var1_u, var0_u;
    double n0_on, n0_off, n1_on, n1_off;
} cpsmatch_statistics;

void cpsmatch_compute_statistics(size_t n, const double *treatment, const double *outcome,
                                 const double *weight, const double *support,
                                 cpsmatch_statistics *result);

/* Publish scalar sufficient statistics on the calling thread; return SPI errors. */
ST_retcode cpsmatch_save_statistics(const cpsmatch_statistics *result);

/* Exact sums of doubles: a fixed-point integer over 32-bit digits, so the
 * sum is independent of the order of its terms and rounds once to double.
 * Stata's summarize gives r(sum) as that rounded exact sum (its quad
 * precision accumulator), r(mean) as r(sum)/r(N) and r(Var) as the rounded
 * exact sum of fl(fl(x - r(mean))^2) over r(N) - 1. Terms must be finite. */
#define CPSMATCH_XDIGITS 68
typedef struct
{
    int64_t d[CPSMATCH_XDIGITS];
    uint32_t pending;
} cpsmatch_xsum;

void cpsmatch_xsum_add(cpsmatch_xsum *a, double x);
void cpsmatch_xsum_merge(cpsmatch_xsum *into, const cpsmatch_xsum *from);
double cpsmatch_xsum_value(const cpsmatch_xsum *a);

/* What _cpsmatch_effects reads with summarize and generate for one outcome,
 * computed on the matching outputs: names follow its locals and r()
 * results (see build/cpsmatch.ado). The weight sums do not depend on the
 * outcome. Missing results are SV_missval. */
typedef struct
{
    double u1u, n1u, var1u, u0u, n0u, var0u; /* summarize v if _treated == 1 / 0 */
    double m1t, n1, m0t;                     /* v and _v, treated on support */
    double m1u, m0u, n0;                     /* _v and v, controls on support */
    double var1, var0;                       /* Var of v, treated / controls on support */
    double sumvt, sumvu, sumve;              /* ai(): sums of the float terms */
} cpsmatch_effects;

typedef struct
{
    double w2c, w2t, w1t, w1c; /* sums of _weight^2 and (1 + _weight)^2, controls / treated */
} cpsmatch_effect_weights;

/* Accumulators the effect passes need per thread */
#define CPSMATCH_EFFECT_SUMS 14

/* arm: the stored _treated; support, weight: final _support and _weight;
 * y the outcome, matched _<outcome>, self _self_<outcome> (NULL without
 * ai). ai: the ai() neighbors (0 without); work holds threads *
 * CPSMATCH_EFFECT_SUMS accumulators. */
void cpsmatch_compute_effects(size_t n, const double *arm, const double *support,
                              const double *weight, const double *y, const double *matched,
                              const double *self, int ai, int ate, int samplevar, int altvariance,
                              int neighbors, int threads, cpsmatch_xsum *work,
                              cpsmatch_effects *out);
void cpsmatch_compute_effect_weights(size_t n, const double *arm, const double *support,
                                     const double *weight, int threads, cpsmatch_xsum *work,
                                     cpsmatch_effect_weights *out);

#endif

#ifndef CPSMATCH_STATS_H
#define CPSMATCH_STATS_H

#include <stddef.h>
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

#endif

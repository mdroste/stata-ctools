/*
 * ccollapse_impl.h
 *
 * Summary statistics by group for Stata: the plugin side of ccollapse, which
 * follows Stata's collapse.
 * Part of the ctools suite
 *
 * The ado converts in two calls, so that nothing in the dataset changes
 * until every statistic has been computed:
 *
 *   compute nby=<b> nsrc=<s> wtype=<w> tlen=<n> [marker] [cw] [verbose]
 *       The plugin varlist holds the by() variables, the distinct source
 *       variables, the weight variable when wtype is not none, and last,
 *       with marker, a sample marker in place of an if condition. An in
 *       range is passed on the plugin call itself.
 *
 *       The targets are in the calling ado's local ccollapse_targets, of tlen
 *       bytes (Stata cuts a plugin's argument string off after a few thousand
 *       bytes): a comma-separated list of <source>:<stat>, where <source> is
 *       the one-based position among the source variables and <stat> one of
 *       mean sum rawsum count percent sd semean sebinomial sepoisson min
 *       max first last firstnm lastnm iqr p1 ... p99.
 *
 *       Reads the selected observations, leaves out those with a zero or
 *       missing weight (and, with cw, those with a missing source value),
 *       groups them and computes one value per group and target. Keeps the
 *       group keys and the values; writes nothing to the dataset. Sets
 *       scalars _ccollapse_ngroups, _ccollapse_nobs, _ccollapse_threads,
 *       _ccollapse_method (how the groups were found: 1 one group, 2 rows
 *       already in order, 3 dense integer keys, 4 hash, 5 sort) and
 *       _ccollapse_time_*. Returns 2000 when no observation is left, 402 for
 *       a negative weight (iweights excepted) and 401 for a noninteger
 *       frequency weight, without a message: the ado issues Stata's own.
 *
 *   store
 *       The plugin varlist holds the by() variables and one variable per
 *       target, in the order of the compute call. Writes the keys and values
 *       of group g to observation g and frees them.
 *
 * Any other call, and any failure, frees what compute kept.
 */

#ifndef CCOLLAPSE_IMPL_H
#define CCOLLAPSE_IMPL_H

#include "stplugin.h"

/*
 * Main entry point for the ccollapse command.
 * Called from the ctools_plugin.c dispatcher.
 *
 * @param args  "compute ..." or "store" (see above)
 * @return      0 on success, Stata error code on failure
 */
ST_retcode ccollapse_main(const char *args);

/* Frees the groups and values compute kept. Idempotent. */
void ccollapse_cleanup_state(void);

#endif /* CCOLLAPSE_IMPL_H */

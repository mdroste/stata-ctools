/*
 * cbinscatter_impl.h
 *
 * Main header for cbinscatter module
 * Part of the ctools Stata plugin suite
 */

#ifndef CBINSCATTER_IMPL_H
#define CBINSCATTER_IMPL_H

#include "stplugin.h"

/*
 * Main entry point for cbinscatter plugin command
 *
 * Called from ctools_plugin.c dispatcher.
 * Expects args to specify subcommand: "compute_bins"
 *
 * Variables expected in plugin call (in order):
 *   y, x, [controls...], [absorb_vars...], [by_var], [weight_var]
 *
 * Configuration read from Stata scalars:
 *   __cbinscatter_nquantiles, __cbinscatter_linetype, etc.
 *
 * Results stored to Stata:
 *   - Matrix __cbinscatter_bins filled with bin data
 *   - Matrix __cbinscatter_coefs filled with fit coefficients
 *   - Scalars __cbinscatter_N, __cbinscatter_N_dropped, etc.
 */
ST_retcode cbinscatter_main(const char *args);

#endif /* CBINSCATTER_IMPL_H */

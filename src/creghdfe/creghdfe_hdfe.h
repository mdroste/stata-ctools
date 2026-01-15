/*
 * creghdfe_hdfe.h
 *
 * HDFE initialization, singleton detection, factor management
 * Part of the ctools Stata plugin suite
 */

#ifndef CREGHDFE_HDFE_H
#define CREGHDFE_HDFE_H

#include "creghdfe_types.h"

/* Clean up global state */
void cleanup_state(void);

/* ========================================================================
 * Plugin command handlers
 * ======================================================================== */

/* HDFE initialization command */
ST_retcode do_hdfe_init(int argc, char *argv[]);

/* Mobility groups computation command */
ST_retcode do_mobility_groups(int argc, char *argv[]);

#endif /* CREGHDFE_HDFE_H */

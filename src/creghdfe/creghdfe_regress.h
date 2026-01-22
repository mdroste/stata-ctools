/*
 * creghdfe_regress.h
 *
 * Full regression command for creghdfe
 * Part of the ctools Stata plugin suite
 */

#ifndef CREGHDFE_REGRESS_H
#define CREGHDFE_REGRESS_H

#include "creghdfe_types.h"

/* Full regression command: HDFE init + partial out + OLS in one shot */
ST_retcode do_full_regression(int argc, char *argv[]);

#endif /* CREGHDFE_REGRESS_H */

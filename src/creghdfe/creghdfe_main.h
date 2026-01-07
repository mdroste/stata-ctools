/*
 * creghdfe_main.h
 *
 * Main entry point and full_regression command
 * Part of the ctools Stata plugin suite
 */

#ifndef CREGHDFE_MAIN_H
#define CREGHDFE_MAIN_H

#include "creghdfe_types.h"

/* Full regression command: HDFE init + partial out + OLS in one shot */
ST_retcode do_full_regression(int argc, char *argv[]);

#endif /* CREGHDFE_MAIN_H */

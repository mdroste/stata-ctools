/*
    cbsample_impl.h
    High-performance bootstrap sampling for Stata

    Author: ctools project
    License: MIT
*/

#ifndef CBSAMPLE_IMPL_H
#define CBSAMPLE_IMPL_H

#include "stplugin.h"

/*
    Main entry point for cbsample command.

    Performs bootstrap sampling with replacement.
    Writes frequency weights to a new variable.

    Arguments format: "freq_var_idx [n=N] [cluster=ncl cl_idx1 ...] [strata=nst st_idx1 ...] [verbose]"

    @param args   Command arguments string
    @return       Stata return code (0 = success)
*/
ST_retcode cbsample_main(const char *args);

#endif /* CBSAMPLE_IMPL_H */

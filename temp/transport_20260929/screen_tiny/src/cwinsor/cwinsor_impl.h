/*
    cwinsor_impl.h
    Header file for cwinsor - High-performance winsorization for Stata

    cwinsor is a C-accelerated replacement for winsor2, providing fast parallel
    winsorization of numeric variables. Key features:

    - Parallel processing across variables using OpenMP
    - O(n) average-case percentile computation via quickselect
    - By-group winsorization support
    - Options for trimming (replace with missing) vs winsorizing
    - Support for generating new variables or replacing in place

    Usage from Stata (via plugin call):
        plugin call ctools_plugin ..., "cwinsor <args>"

    Arguments format:
        "nvars ngroups var_indices... group_sizes... [options]"

    Options:
        p(#)        - Lower percentile cutoff (default 1)
        q(#)        - Upper percentile cutoff (default 99), or use cuts(# #)
        cuts(# #)   - Explicit lower and upper percentiles
        trim        - Trim (replace with missing) instead of winsorize
        verbose     - Display timing breakdown

    Author: ctools project
    License: MIT
*/

#ifndef CWINSOR_IMPL_H
#define CWINSOR_IMPL_H

#include "stplugin.h"

/*
    Main entry point for cwinsor command.

    @param args  Command arguments string from Stata
    @return      ST_retcode: 0 on success, Stata error code on failure
                 198 = syntax error
                 920 = memory allocation error
*/
ST_retcode cwinsor_main(const char *args);

#endif /* CWINSOR_IMPL_H */

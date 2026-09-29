/*
    crangestat_impl.h
    Header file for crangestat - High-performance range statistics for Stata

    crangestat is a C-accelerated replacement for rangestat, computing statistics
    for each observation using all observations where a numeric key variable
    is within the low and high bounds defined for the current observation.

    Key features:
    - Parallel processing using OpenMP
    - Efficient sorted-window computation
    - Multiple statistics in a single pass
    - By-group support for panel data

    Usage from Stata (via plugin call):
        plugin call ctools_plugin ..., "crangestat <args>"

    Statistics supported:
        count     - number of non-missing observations
        mean      - arithmetic mean
        sum       - sum
        min       - minimum
        max       - maximum
        sd        - standard deviation
        variance  - variance
        median    - median (p50)
        iqr       - interquartile range (p75 - p25)
        first     - first value (by key order)
        last      - last value (by key order)
        firstnm   - first non-missing value
        lastnm    - last non-missing value

    Author: ctools project
    License: MIT
*/

#ifndef CRANGESTAT_IMPL_H
#define CRANGESTAT_IMPL_H

#include "stplugin.h"

/*
    Main entry point for crangestat command.

    @param args  Command arguments string from Stata
    @return      ST_retcode: 0 on success, Stata error code on failure
                 198 = syntax error
                 920 = memory allocation error
*/
ST_retcode crangestat_main(const char *args);

#endif /* CRANGESTAT_IMPL_H */

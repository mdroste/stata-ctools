/*
    csample_impl.h
    High-performance random sampling for Stata

    Author: ctools project
    License: MIT
*/

#ifndef CSAMPLE_IMPL_H
#define CSAMPLE_IMPL_H

#include "stplugin.h"

/*
    Main entry point for csample command.

    Performs random sampling without replacement.
    Marks non-selected observations for dropping.

    Arguments format: "percent|count=N [by=nby by_idx1 by_idx2 ...] [verbose]"

    @param args   Command arguments string
    @return       Stata return code (0 = success)
*/
ST_retcode csample_main(const char *args);

#endif /* CSAMPLE_IMPL_H */

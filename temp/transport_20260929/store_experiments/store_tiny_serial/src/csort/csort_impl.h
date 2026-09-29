/*
    csort_impl.h
    csort command implementation

    High-performance radix sort for Stata datasets.
*/

#ifndef CSORT_IMPL_H
#define CSORT_IMPL_H

#include "stplugin.h"

/*
    Main entry point for csort command.

    Performs high-performance radix sort on the current dataset.

    @param args  Space-separated list of 1-based variable indices to sort by

    @return 0 on success, Stata error code on failure
*/
ST_retcode csort_main(const char *args);

#endif /* CSORT_IMPL_H */

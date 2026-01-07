/*
    cmerge_impl.h
    C-Accelerated Merge for Stata

    High-performance drop-in replacement for Stata's merge command.
    Uses parallel data loading, radix sort (if needed), and parallel merge.
*/

#ifndef CMERGE_IMPL_H
#define CMERGE_IMPL_H

#include "stplugin.h"

/*
    Main entry point for cmerge command.

    Performs a merge operation between master (in memory) and using datasets.
    Supports 1:1, m:1, 1:m, and m:m merge types.

    Arguments format:
        "merge_type keyvar_indices using_nobs using_nvars [options]"

    Where:
        merge_type: 0=1:1, 1=m:1, 2=1:m, 3=m:m
        keyvar_indices: space-separated 1-based variable indices for key variables
        using_nobs: number of observations in using dataset
        using_nvars: number of variables in using dataset
        options: verbose, assert(), keep(), generate(), etc.

    The .ado file handles:
        - Loading using dataset into a temporary frame
        - Passing using data via Stata's data interface
        - Creating the _merge variable
        - Setting sort order after merge

    Returns:
        0 on success, Stata error code on failure
*/
ST_retcode cmerge_main(const char *args);

#endif /* CMERGE_IMPL_H */

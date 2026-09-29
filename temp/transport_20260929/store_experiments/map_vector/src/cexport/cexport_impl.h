/*
    cexport_impl.h
    cexport command interface

    High-performance delimited text file export for Stata datasets.
    Part of the ctools suite.
*/

#ifndef CEXPORT_IMPL_H
#define CEXPORT_IMPL_H

#include "stplugin.h"

/*
    Main entry point for cexport command.

    Arguments (passed as space-separated string):
        filename    - Output file path
        delimiter   - Field delimiter character
        Options:
            noheader   - Don't write variable names as header row
            quote      - Quote all string fields
            quoteif    - Quote strings only if they contain delimiter/newline
            replace    - Overwrite existing file
            verbose    - Print progress information

    Returns:
        0 on success, Stata error code on failure
*/
ST_retcode cexport_main(const char *args);

/*
 * Cleanup function for cexport persistent state.
 * Currently a no-op since cexport uses stack-allocated context,
 * but provided for interface consistency.
 * Safe to call multiple times (idempotent).
 */
void cexport_cleanup_state(void);

#endif /* CEXPORT_IMPL_H */

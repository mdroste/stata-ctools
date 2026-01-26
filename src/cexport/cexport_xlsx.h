/*
    cexport_xlsx.h
    Excel (.xlsx) export interface

    High-performance Excel export for Stata datasets.
    Part of the ctools suite.
*/

#ifndef CEXPORT_XLSX_H
#define CEXPORT_XLSX_H

#include "stplugin.h"
#include <stdbool.h>
#include <stdint.h>

/*
    Main entry point for cexport excel command.

    Arguments (passed as space-separated string):
        filename    - Output file path (.xlsx)
        Options:
            sheet=name     - Worksheet name (default: "Sheet1")
            firstrow       - Write variable names as first row
            replace        - Overwrite existing file
            nolabel        - Export numeric values instead of value labels
            dateformat=fmt - Date format string
            verbose        - Print progress information

    Returns:
        0 on success, Stata error code on failure
*/
ST_retcode cexport_xlsx_main(const char *args);

#endif /* CEXPORT_XLSX_H */

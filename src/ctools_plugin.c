/*
    ctools_plugin.c
    Main plugin entry point for the ctools Stata package

    This is the single entry point for all ctools commands. It dispatches
    to the appropriate command handler based on the first argument.

    Supported commands:
    - csort: High-performance radix sort
    - creghdfe: C-accelerated reghdfe (high-dimensional fixed effects regression)
    - civreghdfe: C-accelerated ivreghdfe (IV with high-dimensional fixed effects)
    - cimport: High-performance CSV import (replaces import delimited)
    - cexport: High-performance CSV export (replaces export delimited)
    - cmerge: C-accelerated merge (replaces merge)
    - cqreg: C-accelerated quantile regression (replaces qreg)
    - cbinscatter: C-accelerated binned scatter plots

    Usage from Stata:
        plugin call ctools_plugin ..., "csort <args>"
        plugin call ctools_plugin ..., "creghdfe <args>"
        plugin call ctools_plugin ..., "civreghdfe <args>"
        plugin call ctools_plugin ..., "cimport <args>"
        plugin call ctools_plugin ..., "cexport <args>"
        plugin call ctools_plugin ..., "cmerge <args>"
        plugin call ctools_plugin ..., "cqreg <args>"
        plugin call ctools_plugin ..., "cbinscatter <args>"
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "stplugin.h"
#include "csort_impl.h"
#include "creghdfe_impl.h"
#include "civreghdfe_impl.h"
#include "cimport_impl.h"
#include "cexport_impl.h"
#include "cmerge_impl.h"
#include "cqreg_impl.h"
#include "cbinscatter_impl.h"

/*
    Main plugin entry point.

    Arguments (passed from Stata):
        argv[0]: Command string in format "command args..."
                 First word is the command name, rest are command-specific args

    Returns:
        0 on success, Stata error code on failure
*/
STDLL stata_call(int argc, char *argv[])
{
    char *cmd_str;
    char *cmd_name;
    char *cmd_args;
    char *space_pos;
    ST_retcode rc;


    /* Check for arguments */
    if (argc < 1 || argv[0] == NULL || strlen(argv[0]) == 0) {
        SF_error("ctools: no command specified\n");
        return 198;
    }



    /* Make a copy of the command string to parse */
    cmd_str = strdup(argv[0]);
    if (cmd_str == NULL) {
        SF_error("ctools: memory allocation failed\n");
        return 920;
    }


    /* Find the first space to separate command from args */
    space_pos = strchr(cmd_str, ' ');
    if (space_pos != NULL) {
        *space_pos = '\0';
        cmd_name = cmd_str;
        cmd_args = space_pos + 1;
        /* Skip leading whitespace in args */
        while (*cmd_args == ' ' || *cmd_args == '\t') cmd_args++;
    } else {
        cmd_name = cmd_str;
        cmd_args = "";
    }


    /* Dispatch to appropriate command handler */
    if (strcmp(cmd_name, "csort") == 0) {
        rc = csort_main(cmd_args);
    }
    else if (strcmp(cmd_name, "creghdfe") == 0) {
        rc = creghdfe_main(cmd_args);
    }
    else if (strcmp(cmd_name, "cimport") == 0) {
        rc = cimport_main(cmd_args);
    }
    else if (strcmp(cmd_name, "cexport") == 0) {
        rc = cexport_main(cmd_args);
    }
    else if (strcmp(cmd_name, "cmerge") == 0) {
        rc = cmerge_main(cmd_args);
    }
    else if (strcmp(cmd_name, "cqreg") == 0) {
        rc = cqreg_main(cmd_args);
    }
    else if (strcmp(cmd_name, "cbinscatter") == 0) {
        rc = cbinscatter_main(cmd_args);
    }
    else if (strcmp(cmd_name, "civreghdfe") == 0) {
        rc = civreghdfe_main(cmd_args);
    }
    else {
        char msg[256];
        snprintf(msg, sizeof(msg), "ctools: unknown command '%s'\n", cmd_name);
        SF_error(msg);
        free(cmd_str);
        return 198;
    }

    free(cmd_str);
    return rc;
}

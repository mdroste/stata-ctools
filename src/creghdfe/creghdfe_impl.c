/*
 * creghdfe_impl.c
 *
 * Main entry point and dispatcher for creghdfe command
 * Part of the ctools Stata plugin suite
 *
 * This file provides the main entry point that dispatches to:
 *   - creghdfe_hdfe.c: HDFE initialization, singleton detection
 *   - creghdfe_main.c: Full regression (init + partial out + OLS + VCE)
 */

#include "creghdfe_impl.h"
#include "creghdfe_hdfe.h"
#include "creghdfe_main.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*
 * Main entry point for creghdfe command.
 * Parses the subcommand and dispatches to appropriate handler.
 */
ST_retcode creghdfe_main(const char *args)
{
    char *args_copy;
    char *subcommand;
    char *space_pos;
    ST_retcode rc;

    /* Handle empty args */
    if (args == NULL || strlen(args) == 0) {
        SF_error("creghdfe: no subcommand specified\n");
        return 198;
    }

    /* Make a copy to parse */
    args_copy = strdup(args);
    if (args_copy == NULL) {
        SF_error("creghdfe: memory allocation failed\n");
        return 920;
    }

    /* Extract subcommand (first word) */
    space_pos = strchr(args_copy, ' ');
    if (space_pos != NULL) {
        *space_pos = '\0';
    }
    subcommand = args_copy;

    /* Dispatch to appropriate handler */
    if (strcmp(subcommand, "hdfe_init") == 0) {
        rc = do_hdfe_init(0, NULL);
    }
    else if (strcmp(subcommand, "mobility_groups") == 0) {
        rc = do_mobility_groups(0, NULL);
    }
    else if (strcmp(subcommand, "full_regression") == 0) {
        rc = do_full_regression(0, NULL);
    }
    else {
        char msg[256];
        snprintf(msg, sizeof(msg), "creghdfe: unknown subcommand '%s'\n", subcommand);
        SF_error(msg);
        free(args_copy);
        return 198;
    }

    free(args_copy);
    return rc;
}

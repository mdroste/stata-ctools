/*
 * cppmlhdfe_impl.c
 *
 * Main entry point and dispatcher for cppmlhdfe command
 * Part of the ctools Stata plugin suite
 */

#include "cppmlhdfe_impl.h"
#include "cppmlhdfe_irls.h"
#include "cppmlhdfe_boot.h"
#include "../ctools_parse.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*
 * Main entry point for cppmlhdfe command.
 * Parses the subcommand and dispatches to appropriate handler.
 */
ST_retcode cppmlhdfe_main(const char *args)
{
    char *args_copy;
    char *subcommand;
    char *space_pos;
    ST_retcode rc;

    /* Handle empty args */
    if (args == NULL || strlen(args) == 0) {
        SF_error("cppmlhdfe: no subcommand specified\n");
        return 198;
    }

    /* Make a copy to parse */
    args_copy = strdup(args);
    if (args_copy == NULL) {
        SF_error("cppmlhdfe: memory allocation failed\n");
        return 920;
    }

    /* Extract subcommand (first word) */
    space_pos = strchr(args_copy, ' ');
    if (space_pos != NULL) {
        *space_pos = '\0';
    }
    subcommand = args_copy;

    /* vce(bootstrap) is linked here, not into the fit (cppmlhdfe_boot.h) */
    cppmlhdfe_boot_install();

    /* Dispatch to appropriate handler */
    if (strcmp(subcommand, "full_regression") == 0) {
        /* retain keeps the model for vce(bootstrap) (ctools_runtime.c holds
         * the command in phase 1 until boot_finish) */
        cppmlhdfe_boot_retain_requested = ctools_parse_bool_option(args, "retain");
        rc = do_ppml_regression(0, NULL);
        cppmlhdfe_boot_retain_requested = 0;
    }
    else if (strcmp(subcommand, "boot_chunk") == 0) {
        rc = cppmlhdfe_boot_chunk(args);
    }
    else if (strcmp(subcommand, "boot_finish") == 0) {
        rc = cppmlhdfe_boot_finish(args);
    }
    else if (strcmp(subcommand, "wild_setup") == 0) {
        rc = cppmlhdfe_wild_setup(args);
    }
    else if (strcmp(subcommand, "wild_chunk") == 0) {
        rc = cppmlhdfe_wild_chunk(args);
    }
    else if (strcmp(subcommand, "wild_finish") == 0) {
        rc = cppmlhdfe_wild_finish(args);
    }
    else if (strcmp(subcommand, "boot_release") == 0) {
        /* The phase reset at the start of this call freed the model */
        cppmlhdfe_boot_release();
        rc = 0;
    }
    else {
        char msg[256];
        snprintf(msg, sizeof(msg), "cppmlhdfe: unknown subcommand '%s'\n", subcommand);
        SF_error(msg);
        free(args_copy);
        return 198;
    }

    free(args_copy);
    return rc;
}

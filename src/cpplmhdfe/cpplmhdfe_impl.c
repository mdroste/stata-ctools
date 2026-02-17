/*
 * cpplmhdfe_impl.c
 *
 * Main entry point and dispatcher for cpplmhdfe command
 * Part of the ctools Stata plugin suite
 */

#include "cpplmhdfe_impl.h"
#include "cpplmhdfe_irls.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*
 * Main entry point for cpplmhdfe command.
 * Parses the subcommand and dispatches to appropriate handler.
 */
ST_retcode cpplmhdfe_main(const char *args)
{
    char *args_copy;
    char *subcommand;
    char *space_pos;
    ST_retcode rc;

    /* Handle empty args */
    if (args == NULL || strlen(args) == 0) {
        SF_error("cpplmhdfe: no subcommand specified\n");
        return 198;
    }

    /* Make a copy to parse */
    args_copy = strdup(args);
    if (args_copy == NULL) {
        SF_error("cpplmhdfe: memory allocation failed\n");
        return 920;
    }

    /* Extract subcommand (first word) */
    space_pos = strchr(args_copy, ' ');
    if (space_pos != NULL) {
        *space_pos = '\0';
    }
    subcommand = args_copy;

    /* Dispatch to appropriate handler */
    if (strcmp(subcommand, "full_regression") == 0) {
        rc = do_ppml_regression(0, NULL);
    }
    else {
        char msg[256];
        snprintf(msg, sizeof(msg), "cpplmhdfe: unknown subcommand '%s'\n", subcommand);
        SF_error(msg);
        free(args_copy);
        return 198;
    }

    free(args_copy);
    return rc;
}

/*
 * Cleanup function (idempotent).
 * cpplmhdfe doesn't maintain persistent global state,
 * but we provide this for the cleanup infrastructure.
 */
void cpplmhdfe_cleanup_state(void)
{
    /* No persistent global state to clean up */
}

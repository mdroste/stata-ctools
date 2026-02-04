/*
 * ctools_cleanup.c
 * Centralized cleanup implementation for persistent global state
 */

#include <stdlib.h>
#include "ctools_cleanup.h"

/* Track the current command for multi-phase operation support */
static ctools_command_t g_current_command = CTOOLS_CMD_NONE;
static ctools_command_t g_previous_command = CTOOLS_CMD_NONE;

void ctools_set_current_command(ctools_command_t cmd)
{
    g_previous_command = g_current_command;
    g_current_command = cmd;
}

ctools_command_t ctools_get_current_command(void)
{
    return g_current_command;
}

/*
 * Clean up stale state from other commands.
 *
 * Logic:
 * - If the current command is the same as the previous command, preserve state
 *   (supports multi-phase operations like cmerge load_using -> execute)
 * - If switching to a different command, clean up the previous command's state
 * - Special case: CTOOLS_CMD_NONE cleans everything (used at plugin exit)
 */
void ctools_cleanup_stale_state(ctools_command_t current_cmd)
{
    /* If same command as before, preserve state for multi-phase operations */
    if (current_cmd != CTOOLS_CMD_NONE && current_cmd == g_previous_command) {
        return;
    }

    /* Clean up state from other commands */

    /* cmerge cache - only clean if NOT running cmerge */
    if (current_cmd != CTOOLS_CMD_CMERGE) {
        cmerge_cleanup_cache();
    }

    /* cimport cache - only clean if NOT running cimport */
    if (current_cmd != CTOOLS_CMD_CIMPORT) {
        cimport_cleanup_cache();
    }

    /* creghdfe and civreghdfe clean themselves up, but call for safety */
    if (current_cmd != CTOOLS_CMD_CREGHDFE) {
        creghdfe_cleanup_state();
    }
    if (current_cmd != CTOOLS_CMD_CIVREGHDFE) {
        civreghdfe_cleanup_state();
    }

    /* cexport usually has no persistent state, but call for completeness */
    if (current_cmd != CTOOLS_CMD_CEXPORT) {
        cexport_cleanup_state();
    }
}

void ctools_cleanup_all(void)
{
    /* Force cleanup of everything */
    cmerge_cleanup_cache();
    cimport_cleanup_cache();
    creghdfe_cleanup_state();
    civreghdfe_cleanup_state();
    cexport_cleanup_state();

    /* Reset command tracking */
    g_current_command = CTOOLS_CMD_NONE;
    g_previous_command = CTOOLS_CMD_NONE;
}

/*
 * ctools_cleanup.h
 * Centralized cleanup for persistent global state
 *
 * This module provides cleanup functions for all ctools commands that maintain
 * persistent state between plugin calls. Called automatically at plugin exit
 * to prevent memory leaks from interrupted operations.
 *
 * Design principles:
 * - All cleanup functions are idempotent (safe to call multiple times)
 * - Cleanup functions handle NULL/uninitialized state gracefully
 * - State required for multi-phase operations (cmerge, cimport) is preserved
 *   when calling the same command, but cleaned when switching commands
 */

#ifndef CTOOLS_CLEANUP_H
#define CTOOLS_CLEANUP_H

#include "stplugin.h"

/*
 * Command identifiers for tracking which command is running.
 * Used to preserve state for multi-phase operations within the same command.
 */
typedef enum {
    CTOOLS_CMD_NONE = 0,
    CTOOLS_CMD_CSORT,
    CTOOLS_CMD_CMERGE,
    CTOOLS_CMD_CIMPORT,
    CTOOLS_CMD_CEXPORT,
    CTOOLS_CMD_CREGHDFE,
    CTOOLS_CMD_CIVREGHDFE,
    CTOOLS_CMD_CQREG,
    CTOOLS_CMD_CBINSCATTER,
    CTOOLS_CMD_CENCODE,
    CTOOLS_CMD_CDECODE,
    CTOOLS_CMD_CWINSOR,
    CTOOLS_CMD_CDESTRING,
    CTOOLS_CMD_CSAMPLE,
    CTOOLS_CMD_CBSAMPLE,
    CTOOLS_CMD_CRANGESTAT,
    CTOOLS_CMD_CPSMATCH,
    CTOOLS_CMD_OTHER
} ctools_command_t;

/*
 * Set the current command being executed.
 * Call at the start of each command dispatch.
 */
void ctools_set_current_command(ctools_command_t cmd);

/*
 * Get the current command being executed.
 */
ctools_command_t ctools_get_current_command(void);

/*
 * Clean up stale state from other commands.
 * Called at the start of each plugin invocation to free state from
 * previously interrupted commands, while preserving state for the
 * current command's multi-phase operations.
 *
 * @param current_cmd  The command about to execute (its state is preserved)
 */
void ctools_cleanup_stale_state(ctools_command_t current_cmd);

/*
 * Force cleanup of all persistent state.
 * Called on explicit user request or when Stata is exiting.
 * Cleans up ALL state regardless of which command is running.
 */
void ctools_cleanup_all(void);

/*
 * Individual command cleanup functions.
 * These are called by ctools_cleanup_stale_state() and ctools_cleanup_all().
 * Each function is idempotent and handles uninitialized state gracefully.
 */

/* cmerge: Clears the using dataset cache (g_using_cache) */
void cmerge_cleanup_cache(void);

/* cimport: Clears the parsed CSV context cache (g_cimport_ctx) */
void cimport_cleanup_cache(void);

/* creghdfe: Clears the HDFE state (g_state) - usually already cleaned */
void creghdfe_cleanup_state(void);

/* civreghdfe: Clears the IV-HDFE state - usually already cleaned */
void civreghdfe_cleanup_state(void);

/* cexport: Clears any export context - usually stack-allocated */
void cexport_cleanup_state(void);

#endif /* CTOOLS_CLEANUP_H */

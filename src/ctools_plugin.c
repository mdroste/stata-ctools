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
    - cencode: C-accelerated string encoding (replaces encode)
    - cwinsor: C-accelerated winsorization (replaces winsor2)
    - cdestring: C-accelerated string to numeric conversion (replaces destring)
    - cdecode: C-accelerated numeric to string decoding (replaces decode)
    - crangestat: C-accelerated range statistics (replaces rangestat)

    Usage from Stata:
        plugin call ctools_plugin ..., "csort <args>"
        plugin call ctools_plugin ..., "creghdfe <args>"
        plugin call ctools_plugin ..., "civreghdfe <args>"
        plugin call ctools_plugin ..., "cimport <args>"
        plugin call ctools_plugin ..., "cexport <args>"
        plugin call ctools_plugin ..., "cmerge <args>"
        plugin call ctools_plugin ..., "cqreg <args>"
        plugin call ctools_plugin ..., "cbinscatter <args>"
        plugin call ctools_plugin ..., "cencode <args>"
        plugin call ctools_plugin ..., "cwinsor <args>"
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_config.h"
#include "ctools_threads.h"
#include "ctools_cleanup.h"

/*
    Initialize OpenMP safely to avoid conflicts with other OpenMP runtimes
    that may already be loaded (e.g., by Stata's Java VM).

    The KMP_DUPLICATE_LIB_OK environment variable tells the LLVM/Intel OpenMP
    runtime to tolerate multiple copies of the runtime library being loaded.
    Without this, the runtime aborts when it detects another copy.
*/
static int ctools_omp_initialized = 0;

static void ctools_init_openmp(void)
{
    if (ctools_omp_initialized) return;
    ctools_omp_initialized = 1;

#ifdef _OPENMP
    /* Set KMP_DUPLICATE_LIB_OK to prevent abort when another OpenMP runtime
       is already loaded (common with Stata's Java VM) */
    #if defined(__APPLE__) || defined(__linux__)
    setenv("KMP_DUPLICATE_LIB_OK", "TRUE", 0);  /* 0 = don't overwrite if set */
    #elif defined(_WIN32)
    _putenv_s("KMP_DUPLICATE_LIB_OK", "TRUE");
    #endif
#endif
}

/*
    Parse threads(N) argument from command string.
    If found, sets the global thread limit and removes the argument from the string.

    Format: threads(N) where N is a positive integer

    @param cmd_args   [in/out] Command arguments string (modified in place)
    @return           Thread count if found (>0), 0 if not found, -1 on error
*/
static int parse_threads_arg(char *cmd_args)
{
    char *threads_start;
    char *paren_open;
    char *paren_close;
    int num_threads;
    char *rest;

    if (cmd_args == NULL || *cmd_args == '\0') {
        return 0;  /* No args */
    }

    /* Find "threads(" in the argument string */
    threads_start = strstr(cmd_args, "threads(");
    if (threads_start == NULL) {
        return 0;  /* Not found */
    }

    paren_open = threads_start + 7;  /* Point to '(' */
    if (*paren_open != '(') {
        return 0;  /* Malformed */
    }

    paren_close = strchr(paren_open, ')');
    if (paren_close == NULL) {
        return -1;  /* Missing closing paren */
    }

    /* Parse the number - temporarily null-terminate at closing paren */
    *paren_close = '\0';
    int parse_ok = ctools_safe_atoi(paren_open + 1, &num_threads);
    *paren_close = ')';  /* Restore for memmove below */
    if (!parse_ok) {
        return -1;  /* Invalid number format */
    }
    if (num_threads <= 0) {
        return -1;  /* Invalid number */
    }

    /* Remove the threads() argument from the string:
       Copy everything after "threads(N)" over the threads(...) portion */
    rest = paren_close + 1;
    while (*rest == ' ' || *rest == '\t') rest++;  /* Skip whitespace after ) */

    /* Move the rest of the string over threads(...) */
    memmove(threads_start, rest, strlen(rest) + 1);

    /* Trim trailing whitespace if threads() was at the end */
    size_t len = strlen(cmd_args);
    while (len > 0 && (cmd_args[len-1] == ' ' || cmd_args[len-1] == '\t')) {
        cmd_args[--len] = '\0';
    }

    return num_threads;
}
#include "csort_impl.h"
#include "creghdfe_impl.h"
#include "civreghdfe_impl.h"
#include "cimport_impl.h"
#include "cexport_impl.h"
#include "cexport_xlsx.h"
#include "cmerge_impl.h"
#include "cqreg_impl.h"
#include "cbinscatter_impl.h"
#include "cencode_impl.h"
#include "cwinsor_impl.h"
#include "cdestring_impl.h"
#include "cdecode_impl.h"
#include "csample_impl.h"
#include "cbsample_impl.h"
#include "crangestat_impl.h"
#include "cpsmatch_impl.h"

/*
    Map command name to command ID for cleanup system.
    Returns the command type for tracking multi-phase operations.
*/
static ctools_command_t get_command_type(const char *cmd_name)
{
    if (strcmp(cmd_name, "csort") == 0) return CTOOLS_CMD_CSORT;
    if (strcmp(cmd_name, "cmerge") == 0) return CTOOLS_CMD_CMERGE;
    if (strcmp(cmd_name, "cimport") == 0) return CTOOLS_CMD_CIMPORT;
    if (strcmp(cmd_name, "cexport") == 0) return CTOOLS_CMD_CEXPORT;
    if (strcmp(cmd_name, "cexport_xlsx") == 0) return CTOOLS_CMD_CEXPORT;
    if (strcmp(cmd_name, "creghdfe") == 0) return CTOOLS_CMD_CREGHDFE;
    if (strcmp(cmd_name, "civreghdfe") == 0) return CTOOLS_CMD_CIVREGHDFE;
    if (strcmp(cmd_name, "cqreg") == 0) return CTOOLS_CMD_CQREG;
    if (strcmp(cmd_name, "cbinscatter") == 0) return CTOOLS_CMD_CBINSCATTER;
    if (strcmp(cmd_name, "cencode") == 0) return CTOOLS_CMD_CENCODE;
    if (strcmp(cmd_name, "cdecode") == 0) return CTOOLS_CMD_CDECODE;
    if (strcmp(cmd_name, "cwinsor") == 0) return CTOOLS_CMD_CWINSOR;
    if (strcmp(cmd_name, "cdestring") == 0) return CTOOLS_CMD_CDESTRING;
    if (strcmp(cmd_name, "csample") == 0) return CTOOLS_CMD_CSAMPLE;
    if (strcmp(cmd_name, "cbsample") == 0) return CTOOLS_CMD_CBSAMPLE;
    if (strcmp(cmd_name, "crangestat") == 0) return CTOOLS_CMD_CRANGESTAT;
    if (strcmp(cmd_name, "cpsmatch") == 0) return CTOOLS_CMD_CPSMATCH;
    return CTOOLS_CMD_OTHER;
}

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
    int thread_count;

    /* Initialize OpenMP safely before any parallel code runs */
    ctools_init_openmp();

    /* Reset thread limit to default at start of each call */
    ctools_reset_max_threads();

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

    /* Parse threads(N) argument and set thread limit if specified */
    thread_count = parse_threads_arg(cmd_args);
    if (thread_count < 0) {
        SF_error("ctools: invalid threads() argument\n");
        free(cmd_str);
        return 198;
    }
    if (thread_count > 0) {
        ctools_set_max_threads(thread_count);
    }

    /* Cleanup stale state from other commands before dispatch.
     * This preserves state for multi-phase operations (like cmerge load_using
     * followed by cmerge execute) while freeing state from interrupted
     * operations when switching to a different command. */
    ctools_command_t current_cmd = get_command_type(cmd_name);
    ctools_cleanup_stale_state(current_cmd);
    ctools_set_current_command(current_cmd);

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
    else if (strcmp(cmd_name, "cexport_xlsx") == 0) {
        rc = cexport_xlsx_main(cmd_args);
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
    else if (strcmp(cmd_name, "cencode") == 0) {
        rc = cencode_main(cmd_args);
    }
    else if (strcmp(cmd_name, "cwinsor") == 0) {
        rc = cwinsor_main(cmd_args);
    }
    else if (strcmp(cmd_name, "cdestring") == 0) {
        rc = cdestring_main(cmd_args);
    }
    else if (strcmp(cmd_name, "cdecode") == 0) {
        rc = cdecode_main(cmd_args);
    }
    else if (strcmp(cmd_name, "csample") == 0) {
        rc = csample_main(cmd_args);
    }
    else if (strcmp(cmd_name, "cbsample") == 0) {
        rc = cbsample_main(cmd_args);
    }
    else if (strcmp(cmd_name, "crangestat") == 0) {
        rc = crangestat_main(cmd_args);
    }
    else if (strcmp(cmd_name, "cpsmatch") == 0) {
        rc = cpsmatch_main(cmd_args);
    }
    else {
        char msg[256];
        snprintf(msg, sizeof(msg), "ctools: unknown command '%s'\n", cmd_name);
        SF_error(msg);
        free(cmd_str);
        return 198;
    }

    free(cmd_str);

    /* Clean up thread pool at end of each plugin call */
    ctools_destroy_global_pool();

    return rc;
}

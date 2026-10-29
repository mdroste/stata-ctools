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

#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE  /* dladdr */
#endif
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_config.h"
#include "ctools_threads.h"
#include "ctools_runtime.h"
#include "ctools_version.h"

/*
    Initialize OpenMP safely to avoid conflicts with other OpenMP runtimes
    that may already be loaded (e.g., by Stata's Java VM).

    The KMP_DUPLICATE_LIB_OK environment variable tells the LLVM/Intel OpenMP
    runtime to tolerate multiple copies of the runtime library being loaded.
    Without this, the runtime aborts when it detects another copy.
*/
static int ctools_omp_initialized = 0;

/*
    Keep this image mapped for the life of the process. The OpenMP runtime and
    the persistent thread pool register pthread_atfork handlers and park worker
    threads inside the plugin; if Stata unloads it (e.g. clear all), the next
    fork (shell, !) or pool wakeup would execute unmapped code.
*/
static void ctools_pin_image(void)
{
#ifdef _WIN32
    HMODULE self;
    GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_PIN,
                       (LPCSTR)(void *)&ctools_pin_image, &self);
#else
    Dl_info info;
    if (dladdr((void *)&ctools_pin_image, &info) && info.dli_fname)
        (void)dlopen(info.dli_fname, RTLD_NOW | RTLD_NOLOAD | RTLD_NODELETE);
#endif
}

static void ctools_init_openmp(void)
{
    if (ctools_omp_initialized) return;
    ctools_omp_initialized = 1;
    ctools_pin_image();

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

    /* Only the leading option token belongs to the dispatcher. */
    threads_start = cmd_args;
    while (*threads_start == ' ' || *threads_start == '\t') threads_start++;
    if (strncmp(threads_start, "threads(", 8) != 0) return 0;

    paren_open = threads_start + 7;  /* Point to '(' */
    if (*paren_open != '(') {
        return 0;  /* Malformed */
    }

    paren_close = strchr(paren_open, ')');
    if (paren_close == NULL || (paren_close[1] != '\0' &&
                               paren_close[1] != ' ' && paren_close[1] != '\t')) {
        return -1;  /* Missing closing paren */
    }

    /* Parse the number - temporarily null-terminate at closing paren */
    *paren_close = '\0';
    int parse_ok = ctools_safe_atoi(paren_open + 1, &num_threads);
    *paren_close = ')';  /* Restore for memmove below */
    if (!parse_ok) {
        return -1;  /* Invalid number format */
    }
    if (num_threads <= 0 || num_threads > CTOOLS_MAX_THREADS) {
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

/*
    Session thread setting from ctools, threads(#): global CTOOLS_THREADS.
    Encoded like ctools_parse_thread_value (>0 set, 0 unset, -1 invalid).
*/
static int session_threads(void)
{
    char value[64] = "";
    if (SF_macro_use("CTOOLS_THREADS", value, (int)sizeof(value)) != 0) return 0;
    return ctools_parse_thread_value(value, 0);
}

static ST_retcode save_thread_local(char *name, int value)
{
    char text[16] = "";
    if (value < 0) snprintf(text, sizeof(text), "invalid");
    else if (value > 0) snprintf(text, sizeof(text), "%d", value);
    return SF_macro_save(name, text);
}

/*
    Thread settings for ctools.ado: "threads_info" reports the policy in
    locals; "threads_save N" saves N to the per-user settings file, or
    removes the saved value when N is 0. No worker threads are started.
*/
static ST_retcode thread_settings_call(const char *args)
{
    const ctools_thread_policy *p;
    char text[32];
    int n = 0, session, effective;
    ST_retcode rc = 0;

    if (strncmp(args, "save ", 5) == 0) {
        if (!ctools_safe_atoi(args + 5, &n) || n < 0 || n > CTOOLS_MAX_THREADS) {
            SF_error("ctools: threads() must be an integer from 1 through 256, or default\n");
            return 198;
        }
        if (ctools_save_thread_setting(n) != 0) {
            char msg[4400];
            const char *path = ctools_thread_policy_get()->settings_path;
            snprintf(msg, sizeof(msg), "ctools: could not update thread setting in %s (%s)\n",
                     path[0] ? path : "home directory", strerror(errno));
            SF_error(msg);
            return 603;
        }
    } else if (strcmp(args, "info") != 0) {
        SF_error("ctools: unknown thread settings request\n");
        return 198;
    }

    p = ctools_thread_policy_get();
    session = session_threads();
    effective = session > 0 ? ctools_clamp_thread_request(session) : p->threads;
    snprintf(text, sizeof(text), "%d", effective);
    rc = rc ? rc : SF_macro_save("_ctools_threads_effective", text);
    snprintf(text, sizeof(text), "%d", p->threads);
    rc = rc ? rc : SF_macro_save("_ctools_threads_default", text);
    rc = rc ? rc : SF_macro_save("_ctools_threads_source", (char *)p->source);
    snprintf(text, sizeof(text), "%d", p->detected);
    rc = rc ? rc : SF_macro_save("_ctools_threads_detected", text);
    snprintf(text, sizeof(text), "%d", p->available);
    rc = rc ? rc : SF_macro_save("_ctools_threads_available", text);
    snprintf(text, sizeof(text), "%d", p->limit);
    rc = rc ? rc : SF_macro_save("_ctools_threads_limit", text);
    snprintf(text, sizeof(text), "%d", CTOOLS_DEFAULT_THREAD_CAP);
    rc = rc ? rc : SF_macro_save("_ctools_threads_cap", text);
    snprintf(text, sizeof(text), "%d", CTOOLS_OPENMP_ENABLED);
    rc = rc ? rc : SF_macro_save("_ctools_threads_openmp", text);
    rc = rc ? rc : save_thread_local("_ctools_threads_session", session);
    rc = rc ? rc : save_thread_local("_ctools_threads_saved", p->saved);
    rc = rc ? rc : save_thread_local("_ctools_threads_env", p->env_ctools);
    rc = rc ? rc : save_thread_local("_ctools_threads_omp", p->env_omp);
    rc = rc ? rc : save_thread_local("_ctools_threads_omp_limit", p->omp_limit);
    rc = rc ? rc : SF_macro_save("_ctools_threads_file", (char *)p->settings_path);
    return rc;
}
#include "csort/csort_impl.h"
#include "creghdfe/creghdfe_impl.h"
#include "civreghdfe/civreghdfe_impl.h"
#include "cimport/cimport_impl.h"
#include "io/cio.h"
#include "io/cexport_xls.h"
#include "cexport/cexport_impl.h"
#include "cexport/cexport_xlsx.h"
#include "cmerge/cmerge_impl.h"
#include "cqreg/cqreg_impl.h"
#include "cbinscatter/cbinscatter_impl.h"
#include "cencode/cencode_impl.h"
#include "cwinsor/cwinsor_impl.h"
#include "cdestring/cdestring_impl.h"
#include "cdecode/cdecode_impl.h"
#include "csample/csample_impl.h"
#include "cbsample/cbsample_impl.h"
#include "crangestat/crangestat_impl.h"
#include "cpsmatch/cpsmatch_impl.h"
#include "cppmlhdfe/cppmlhdfe_impl.h"
#include "cipolate/cipolate_impl.h"
#include "csplit/csplit_impl.h"
#include "crangejoin/crangejoin_impl.h"
#include "cgroup/cgroup_impl.h"

static const ctools_command_descriptor commands[] = {
#define CTOOLS_COMMAND(name, owner, handler, cleanup) {name, owner, handler, cleanup},
#include "ctools_commands.def"
#undef CTOOLS_COMMAND
};

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

    /* Identity queries do not start threads or alter command state. */
    if (argc > 0 && argv[0] && strcmp(argv[0], "version") == 0) {
        ST_retcode version_rc = SF_macro_save("_ctools_plugin_version", CTOOLS_VERSION);
        if (version_rc) return version_rc;
        return SF_macro_save("_ctools_plugin_revision", CTOOLS_BUILD_REVISION);
    }

    /* Read-only monotonic clock query: preserve phase caches and caller timers. */
    if (argc > 0 && argv[0] && strcmp(argv[0], "clock") == 0) {
        char seconds[32];
        snprintf(seconds, sizeof(seconds), "%.17g", ctools_timer_seconds());
        return SF_macro_save("_ctools_clock_value", seconds);
    }

    /* Initialize OpenMP safely before any parallel code runs */
    ctools_init_openmp();

    /* Capture the true system thread count on first call, before any
     * external plugin (reghdfe, gtools, etc.) can pollute it */
    ctools_init_thread_count();

    /* Reset thread limit to default at start of each call */
    ctools_reset_max_threads();

    if (argc > 0 && argv[0] && strncmp(argv[0], "threads_", 8) == 0) {
        return thread_settings_call(argv[0] + 8);
    }

    /* Check for arguments */
    if (argc < 1 || argv[0] == NULL || strlen(argv[0]) == 0) {
        SF_error("ctools: no command specified\n");
        ctools_command_reset();
        return 198;
    }

    /* Make a copy of the command string to parse */
    cmd_str = strdup(argv[0]);
    if (cmd_str == NULL) {
        SF_error("ctools: memory allocation failed\n");
        ctools_command_reset();
        return 920;
    }

    /* Find the first space to separate command from args */
    space_pos = strpbrk(cmd_str, " \t\r\n");
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
        SF_error("ctools: threads() must be an integer from 1 through 256\n");
        free(cmd_str);
        ctools_command_reset();
        return 198;
    }
    if (thread_count == 0) {
        /* No threads() option: the session setting applies, if any. */
        thread_count = session_threads();
        if (thread_count < 0) {
            SF_error("ctools: global CTOOLS_THREADS must be an integer from 1 through 256; "
                     "clear it with ctools, threads(default)\n");
            free(cmd_str);
            ctools_command_reset();
            return 198;
        }
    }
    if (thread_count > 0) {
        /* More threads than usable CPUs only adds contention. */
        ctools_set_max_threads(ctools_clamp_thread_request(thread_count));
    }

    const ctools_command_descriptor *command = NULL;
    for (size_t i = 0; i < sizeof(commands) / sizeof(commands[0]); i++) {
        if (!strcmp(cmd_name, commands[i].name)) { command = &commands[i]; break; }
    }
    rc = ctools_command_begin(command, cmd_args);
    if (!rc) rc = command->handler(cmd_args);
    else SF_error("ctools: unknown command or invalid command phase\n");

    free(cmd_str);

    /* Clean up thread pool at end of each plugin call */
    ctools_destroy_global_pool();
    ctools_command_finish(rc);

    return rc;
}

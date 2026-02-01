/*
    csort_impl.c
    csort command implementation

    High-performance parallel sort for Stata datasets.
    Orchestrates the sorting process:
    1. Load data from Stata to C
    2. Sort the data using IPS4O (default) or other algorithms
    3. Transfer sorted data back to Stata
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_timer.h"
#include "csort_impl.h"
#include "csort_stream.h"

/* Debug flag - set to 1 to enable debug output */
#define CSORT_DEBUG 0

/*
    Parse streaming mode option from argument string.
    Looks for "stream" or "stream(#)" to enable streaming mode.

    Streaming mode:
    - Only loads key (sort) variables into C memory
    - Streams permutation application to non-key variables
    - Reduces memory usage for wide datasets
    - Best when: many non-key columns, limited memory

    The optional number specifies how many variables to load at a time
    (1-16, default 1). Higher values use more memory but may be faster.

    Returns 0 if disabled, 1-16 for the batch size.
*/
static int parse_stream_option(const char *args)
{
    const char *p;
    int batch_size = 1;  /* Default: process 1 variable at a time */

    /* Look for "stream" in the arguments */
    p = strstr(args, "stream");
    if (p == NULL) {
        return 0;
    }

    /* Check if it's "stream=0" (explicitly disabled) */
    if (p[6] == '=' && p[7] == '0') {
        return 0;
    }

    /* Check for "stream(#)" format */
    if (p[6] == '(') {
        char *endptr;
        long val = strtol(p + 7, &endptr, 10);

        /* Validate: must be followed by ')' and be a positive integer */
        if (endptr != p + 7 && *endptr == ')') {
            if (val < 1) {
                batch_size = 1;
            } else if (val > 16) {
                batch_size = 16;  /* Cap at 16 */
            } else {
                batch_size = (int)val;
            }
        }
        /* If parsing fails, use default of 1 */
    }

    return batch_size;
}

/*
    Parse algorithm option from argument string.
    Looks for "alg=X" where X is:
      0 or "lsd"      -> SORT_ALG_LSD
      1 or "msd"      -> SORT_ALG_MSD
      2 or "timsort"  -> SORT_ALG_TIMSORT
      3 or "sample"   -> SORT_ALG_SAMPLE
      4 or "counting" -> SORT_ALG_COUNTING
      5 or "merge"    -> SORT_ALG_MERGE
      6 or "ips4o"    -> SORT_ALG_IPS4O
      7 or "auto"     -> SORT_ALG_AUTO (default)
*/
static sort_algorithm_t parse_algorithm(const char *args)
{
    const char *p;

    /* Look for "alg=" in the arguments */
    p = strstr(args, "alg=");
    if (p == NULL) {
        return SORT_ALG_AUTO;  /* Default: auto-select best algorithm */
    }

    p += 4;  /* Skip "alg=" */

    /* Check for numeric or string algorithm specifier */
    if (*p == '0' || strncmp(p, "lsd", 3) == 0) {
        return SORT_ALG_LSD;
    } else if (*p == '1' || strncmp(p, "msd", 3) == 0) {
        return SORT_ALG_MSD;
    } else if (*p == '2' || strncmp(p, "timsort", 7) == 0) {
        return SORT_ALG_TIMSORT;
    } else if (*p == '3' || strncmp(p, "sample", 6) == 0) {
        return SORT_ALG_SAMPLE;
    } else if (*p == '4' || strncmp(p, "counting", 8) == 0) {
        return SORT_ALG_COUNTING;
    } else if (*p == '5' || strncmp(p, "merge", 5) == 0) {
        return SORT_ALG_MERGE;
    } else if (*p == '6' || strncmp(p, "ips4o", 5) == 0) {
        return SORT_ALG_IPS4O;
    } else if (*p == '7' || strncmp(p, "auto", 4) == 0) {
        return SORT_ALG_AUTO;
    }

    return SORT_ALG_AUTO;  /* Default for unrecognized */
}

/* Parse the sort variable indices from the argument string.
   Only parse numbers that appear before options (alg=, stream, threads) to avoid
   picking up numbers from option values. */
static int parse_sort_vars(const char *args, int **sort_vars, size_t *nsort)
{
    const char *p;
    const char *opts_start;
    char *endptr;
    size_t count = 0;
    size_t capacity = 16;
    int *vars;
    long val;

    vars = (int *)malloc(capacity * sizeof(int));
    if (vars == NULL) {
        return -1;
    }

    /* Find where options start - stop at any of: alg=, stream, threads */
    opts_start = strstr(args, "alg=");
    const char *stream_start = strstr(args, "stream");
    const char *threads_start = strstr(args, "threads");

    /* Use the earliest option as the stopping point */
    if (stream_start != NULL && (opts_start == NULL || stream_start < opts_start)) {
        opts_start = stream_start;
    }
    if (threads_start != NULL && (opts_start == NULL || threads_start < opts_start)) {
        opts_start = threads_start;
    }

    p = args;
    while (*p != '\0') {
        /* Stop if we've reached options */
        if (opts_start != NULL && p >= opts_start) {
            break;
        }

        /* Skip whitespace */
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0') break;

        /* Stop if we've reached options after skipping spaces */
        if (opts_start != NULL && p >= opts_start) {
            break;
        }

        /* Parse integer */
        val = strtol(p, &endptr, 10);
        if (endptr == p) {
            /* No number found, skip this character */
            p++;
            continue;
        }

        /* Add to array */
        if (count >= capacity) {
            capacity *= 2;
            int *new_vars = (int *)realloc(vars, capacity * sizeof(int));
            if (new_vars == NULL) {
                free(vars);
                return -1;
            }
            vars = new_vars;
        }
        vars[count++] = (int)val;
        p = endptr;
    }

    *sort_vars = vars;
    *nsort = count;
    return 0;
}

/*
    Main entry point for csort command.
*/
ST_retcode csort_main(const char *args)
{
    stata_data data;
    stata_timer timer;
    stata_retcode rc;
    int *sort_vars = NULL;
    size_t nsort = 0;
    size_t obs1, nvars;
    double t_start, t_end;
    char msg[256];
    sort_algorithm_t algorithm;

    /* Initialize timer */
    timer.load_time = 0.0;
    timer.sort_time = 0.0;
    timer.store_time = 0.0;
    timer.total_time = 0.0;

    /* Start total timer */
    t_start = ctools_timer_seconds();

    /* Check arguments */
    if (args == NULL || strlen(args) == 0) {
        SF_error("csort: no sort variables specified\n");
        return 198;
    }

    /* Parse sort variable indices from arguments */
    if (parse_sort_vars(args, &sort_vars, &nsort) != 0) {
        SF_error("csort: memory allocation failed\n");
        return 920;
    }

    if (nsort == 0) {
        SF_error("csort: no valid sort variables specified\n");
        free(sort_vars);
        return 198;
    }

    /* Parse algorithm option */
    algorithm = parse_algorithm(args);

    /* Get observation range and variable count from Stata */
    obs1 = SF_in1();
    nvars = SF_nvars();

    if (nvars == 0) {
        SF_error("csort: no variables in dataset\n");
        free(sort_vars);
        return 2000;
    }

    /* Check for streaming mode option (returns 0 if disabled, 1-16 for batch size) */
    int stream_batch_size = parse_stream_option(args);

    /* ================================================================
       MEMORY-EFFICIENT MODE: For large datasets with many columns
       Only loads key variables, streams permutation to non-keys
       ================================================================ */
    if (stream_batch_size > 0) {
        csort_stream_timings stream_timings = {0};

        /* Build array of all variable indices (1-based) */
        int *all_var_indices = (int *)malloc(nvars * sizeof(int));
        if (!all_var_indices) {
            SF_error("csort: memory allocation failed\n");
            free(sort_vars);
            return 920;
        }
        for (size_t j = 0; j < nvars; j++) {
            all_var_indices[j] = (int)(j + 1);
        }

        /* Call streaming sort with batch size */
        rc = csort_stream_sort(sort_vars, nsort, all_var_indices, nvars,
                                algorithm, 0, stream_batch_size, &stream_timings);

        free(all_var_indices);
        free(sort_vars);

        if (rc != STATA_OK) {
            snprintf(msg, sizeof(msg), "csort: streaming sort failed (error %d)\n", rc);
            SF_error(msg);
            return 920;
        }

        /* Calculate total time */
        t_end = ctools_timer_seconds();
        timer.total_time = t_end - t_start;

        /* Store timing results in Stata scalars for streaming mode */
        SF_scal_save("_csort_stream", 1.0);  /* Flag indicating streaming mode was used */
        SF_scal_save("_csort_time_load", stream_timings.load_keys_time);
        SF_scal_save("_csort_time_sort", stream_timings.sort_time);
        SF_scal_save("_csort_time_permute", stream_timings.permute_keys_time);
        SF_scal_save("_csort_time_store", stream_timings.store_keys_time);
        SF_scal_save("_csort_time_stream", stream_timings.stream_nonkeys_time);
        SF_scal_save("_csort_time_cleanup", 0.0);  /* Minimal cleanup in streaming mode */
        SF_scal_save("_csort_time_total", timer.total_time);

        /* Store thread diagnostics */
        CTOOLS_SAVE_THREAD_INFO("_csort");

        return 0;
    }

    /* ================================================================
       STANDARD MODE: Load all data, sort, store
       ================================================================ */

    /* ================================================================
       PHASE 1: Load data from Stata to C
       ================================================================ */
#if CSORT_DEBUG
    {
        char dbg[256];
        snprintf(dbg, sizeof(dbg), "csort: Starting data load, nvars=%zu algorithm=%d\n",
                 nvars, algorithm);
        SF_display(dbg);
    }
#endif
    timer.load_time = ctools_timer_seconds();

    rc = ctools_data_load(&data, nvars);

    timer.load_time = ctools_timer_seconds() - timer.load_time;

    if (rc != STATA_OK) {
        snprintf(msg, sizeof(msg), "csort: failed to load data (error %d)\n", rc);
        SF_error(msg);
        stata_data_free(&data);
        free(sort_vars);
        return 920;
    }

    /* ================================================================
       PHASE 2: Sort the data
       ================================================================ */
#if CSORT_DEBUG
    {
        char dbg[512];
        snprintf(dbg, sizeof(dbg), "csort debug: nobs=%zu nvars=%zu nsort=%zu\n",
                 data.nobs, data.nvars, nsort);
        SF_display(dbg);
        snprintf(dbg, sizeof(dbg), "csort debug: sort_vars = [");
        SF_display(dbg);
        for (size_t i = 0; i < nsort; i++) {
            snprintf(dbg, sizeof(dbg), "%d ", sort_vars[i]);
            SF_display(dbg);
        }
        SF_display("]\n");

        /* Show first 5 observations of first string var (make) */
        if (data.vars[0].type == STATA_TYPE_STRING) {
            SF_display("csort debug: First 5 make values (before sort):\n");
            for (size_t i = 0; i < 5 && i < data.nobs; i++) {
                snprintf(dbg, sizeof(dbg), "  [%zu] %s\n", i, data.vars[0].data.str[i]);
                SF_display(dbg);
            }
        }

        /* Show initial sort_order */
        SF_display("csort debug: Initial sort_order (first 10): ");
        for (size_t i = 0; i < 10 && i < data.nobs; i++) {
            snprintf(dbg, sizeof(dbg), "%zu ", data.sort_order[i]);
            SF_display(dbg);
        }
        SF_display("\n");
    }
#endif

    timer.sort_time = ctools_timer_seconds();
    double t_permute = 0.0;  /* Permutation time */

    /* Call unified sort dispatcher (computes sort_order only) */
    rc = ctools_sort_dispatch(&data, sort_vars, nsort, algorithm);

    /* Record sort computation time */
    timer.sort_time = ctools_timer_seconds() - timer.sort_time;

    /* Apply permutation separately (timed) */
    if (rc == STATA_OK) {
        t_permute = ctools_timer_seconds();
        rc = ctools_apply_permutation(&data);
        t_permute = ctools_timer_seconds() - t_permute;
    }

    if (rc != STATA_OK) {
        snprintf(msg, sizeof(msg), "csort: sort failed (error %d)\n", rc);
        SF_error(msg);
        stata_data_free(&data);
        free(sort_vars);
        return 920;
    }

#if CSORT_DEBUG
    {
        char dbg[512];
        /* Show sort_order after sort */
        SF_display("csort debug: Final sort_order (first 10): ");
        for (size_t i = 0; i < 10 && i < data.nobs; i++) {
            snprintf(dbg, sizeof(dbg), "%zu ", data.sort_order[i]);
            SF_display(dbg);
        }
        SF_display("\n");

        /* Show first 5 make values after sort (data is now permuted) */
        if (data.vars[0].type == STATA_TYPE_STRING) {
            SF_display("csort debug: First 5 make values (after sort):\n");
            for (size_t i = 0; i < 5 && i < data.nobs; i++) {
                snprintf(dbg, sizeof(dbg), "  [%zu] %s\n", i, data.vars[0].data.str[i]);
                SF_display(dbg);
            }
        }
    }
#endif

    /* ================================================================
       PHASE 3: Transfer sorted data back to Stata
       ================================================================ */
    timer.store_time = ctools_timer_seconds();

    rc = ctools_data_store(&data, obs1);

    timer.store_time = ctools_timer_seconds() - timer.store_time;

    if (rc != STATA_OK) {
        snprintf(msg, sizeof(msg), "csort: failed to store data (error %d)\n", rc);
        SF_error(msg);
        stata_data_free(&data);
        free(sort_vars);
        return 920;
    }

    /* Clean up C memory (this can be slow for large datasets!) */
    double t_cleanup = ctools_timer_seconds();
    stata_data_free(&data);
    free(sort_vars);
    t_cleanup = ctools_timer_seconds() - t_cleanup;

    /* Calculate total time (AFTER cleanup, so it's fully accounted) */
    t_end = ctools_timer_seconds();
    timer.total_time = t_end - t_start;

    /* Store timing results in Stata scalars (standard mode) */
    SF_scal_save("_csort_stream", 0.0);  /* Flag indicating standard mode was used */
    SF_scal_save("_csort_time_load", timer.load_time);
    SF_scal_save("_csort_time_sort", timer.sort_time);
    SF_scal_save("_csort_time_permute", t_permute);
    SF_scal_save("_csort_time_store", timer.store_time);
    SF_scal_save("_csort_time_stream", 0.0);  /* No streaming in standard mode */
    SF_scal_save("_csort_time_cleanup", t_cleanup);
    SF_scal_save("_csort_time_total", timer.total_time);

    /* Store thread diagnostics */
    CTOOLS_SAVE_THREAD_INFO("_csort");

    return 0;
}

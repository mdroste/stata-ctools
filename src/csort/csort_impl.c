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
#include "ctools_timer.h"
#include "csort_impl.h"
#include "csort_stream.h"

/* Debug flag - set to 1 to enable debug output */
#define CSORT_DEBUG 0

/*
    Parse memory-efficient mode option from argument string.
    Looks for "memeff" to enable memory-efficient streaming mode.

    Memory-efficient mode:
    - Only loads key (sort) variables into C memory
    - Streams permutation application to non-key variables in blocks
    - Dramatically reduces memory usage for wide datasets
    - Best when: many non-key columns, limited memory

    Returns 1 if memory-efficient mode enabled, 0 otherwise.
*/
static int parse_memeff_option(const char *args)
{
    const char *p;

    /* Look for "memeff" in the arguments */
    p = strstr(args, "memeff");
    if (p == NULL) {
        return 0;
    }

    /* Check if it's "memeff=0" (explicitly disabled) */
    if (p[6] == '=' && p[7] == '0') {
        return 0;
    }

    return 1;
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
      6 or "ips4o"    -> SORT_ALG_IPS4O (default)
*/
static sort_algorithm_t parse_algorithm(const char *args)
{
    const char *p;

    /* Look for "alg=" in the arguments */
    p = strstr(args, "alg=");
    if (p == NULL) {
        return SORT_ALG_IPS4O;  /* Default */
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
    }

    return SORT_ALG_IPS4O;  /* Default for unrecognized */
}

/* Parse the sort variable indices from the argument string.
   Only parse numbers that appear before "alg=" to avoid picking up
   numbers from algorithm names like "ips4o". */
static int parse_sort_vars(const char *args, int **sort_vars, size_t *nsort)
{
    const char *p;
    const char *alg_start;
    char *endptr;
    size_t count = 0;
    size_t capacity = 16;
    int *vars;
    long val;

    vars = (int *)malloc(capacity * sizeof(int));
    if (vars == NULL) {
        return -1;
    }

    /* Find where options start (at "alg=") to stop parsing numbers there */
    alg_start = strstr(args, "alg=");

    p = args;
    while (*p != '\0') {
        /* Stop if we've reached the algorithm option */
        if (alg_start != NULL && p >= alg_start) {
            break;
        }

        /* Skip whitespace */
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0') break;

        /* Stop if we've reached the algorithm option after skipping spaces */
        if (alg_start != NULL && p >= alg_start) {
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

    /* Check for memory-efficient mode option */
    int use_memeff = parse_memeff_option(args);

    /* ================================================================
       MEMORY-EFFICIENT MODE: For large datasets with many columns
       Only loads key variables, streams permutation to non-keys
       ================================================================ */
    if (use_memeff) {
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

        /* Call streaming sort */
        rc = csort_stream_sort(sort_vars, nsort, all_var_indices, nvars,
                                algorithm, 0, &stream_timings);

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

        /* Store timing results in Stata scalars for memeff mode */
        SF_scal_save("_csort_memeff", 1.0);  /* Flag indicating memeff mode was used */
        SF_scal_save("_csort_time_load", stream_timings.load_keys_time);
        SF_scal_save("_csort_time_sort", stream_timings.sort_time);
        SF_scal_save("_csort_time_permute", stream_timings.permute_keys_time);
        SF_scal_save("_csort_time_store", stream_timings.store_keys_time);
        SF_scal_save("_csort_time_stream", stream_timings.stream_nonkeys_time);
        SF_scal_save("_csort_time_cleanup", 0.0);  /* Minimal cleanup in memeff mode */
        SF_scal_save("_csort_time_total", timer.total_time);

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

    /* Call the selected sort algorithm using _order_only variants
       to track sort vs permutation timing separately for all algorithms. */
    switch (algorithm) {
        case SORT_ALG_MSD:
            rc = ctools_sort_radix_msd_order_only(&data, sort_vars, nsort);
            break;
        case SORT_ALG_TIMSORT:
            rc = ctools_sort_timsort_order_only(&data, sort_vars, nsort);
            break;
        case SORT_ALG_SAMPLE:
            rc = ctools_sort_sample_order_only(&data, sort_vars, nsort);
            break;
        case SORT_ALG_COUNTING:
            rc = ctools_sort_counting_order_only(&data, sort_vars, nsort);
            /* If counting sort not suitable, fall back to LSD radix */
            if (rc == STATA_ERR_UNSUPPORTED_TYPE) {
                rc = ctools_sort_radix_lsd_order_only(&data, sort_vars, nsort);
            }
            break;
        case SORT_ALG_MERGE:
            rc = ctools_sort_merge_order_only(&data, sort_vars, nsort);
            break;
        case SORT_ALG_LSD:
            rc = ctools_sort_radix_lsd_order_only(&data, sort_vars, nsort);
            break;
        case SORT_ALG_IPS4O:
        default:
            rc = ctools_sort_ips4o_order_only(&data, sort_vars, nsort);
            break;
    }

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
    SF_scal_save("_csort_memeff", 0.0);  /* Flag indicating standard mode was used */
    SF_scal_save("_csort_time_load", timer.load_time);
    SF_scal_save("_csort_time_sort", timer.sort_time);
    SF_scal_save("_csort_time_permute", t_permute);
    SF_scal_save("_csort_time_store", timer.store_time);
    SF_scal_save("_csort_time_stream", 0.0);  /* No streaming in standard mode */
    SF_scal_save("_csort_time_cleanup", t_cleanup);
    SF_scal_save("_csort_time_total", timer.total_time);

    return 0;
}

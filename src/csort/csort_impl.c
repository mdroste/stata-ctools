/*
    csort_impl.c
    csort command implementation

    High-performance radix sort for Stata datasets.
    Orchestrates the sorting process:
    1. Load data from Stata to C
    2. Sort the data using LSD radix sort
    3. Transfer sorted data back to Stata
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_timer.h"
#include "csort_impl.h"

/* Parse the sort variable indices from the argument string */
static int parse_sort_vars(const char *args, int **sort_vars, size_t *nsort)
{
    const char *p;
    char *endptr;
    size_t count = 0;
    size_t capacity = 16;
    int *vars;
    long val;

    vars = (int *)malloc(capacity * sizeof(int));
    if (vars == NULL) {
        return -1;
    }

    p = args;
    while (*p != '\0') {
        /* Skip whitespace */
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0') break;

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
            vars = (int *)realloc(vars, capacity * sizeof(int));
            if (vars == NULL) {
                return -1;
            }
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

    /* Get observation range and variable count from Stata */
    obs1 = SF_in1();
    nvars = SF_nvars();

    if (nvars == 0) {
        SF_error("csort: no variables in dataset\n");
        free(sort_vars);
        return 2000;
    }

    /* ================================================================
       PHASE 1: Load data from Stata to C
       ================================================================ */
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
    timer.sort_time = ctools_timer_seconds();

    rc = ctools_sort_radix_lsd(&data, sort_vars, nsort);

    timer.sort_time = ctools_timer_seconds() - timer.sort_time;

    if (rc != STATA_OK) {
        snprintf(msg, sizeof(msg), "csort: sort failed (error %d)\n", rc);
        SF_error(msg);
        stata_data_free(&data);
        free(sort_vars);
        return 920;
    }

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

    /* Calculate total time */
    t_end = ctools_timer_seconds();
    timer.total_time = t_end - t_start;

    /* Store timing results in Stata scalars */
    SF_scal_save("_csort_time_load", timer.load_time);
    SF_scal_save("_csort_time_sort", timer.sort_time);
    SF_scal_save("_csort_time_store", timer.store_time);
    SF_scal_save("_csort_time_total", timer.total_time);

    /* Clean up */
    stata_data_free(&data);
    free(sort_vars);

    return 0;
}

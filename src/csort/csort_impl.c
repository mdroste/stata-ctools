/*
    csort_impl.c
    csort command implementation

    High-performance parallel sort for Stata datasets.
    Orchestrates the sorting process:
    1. Load data from Stata to C: whole observations, row by row, when every
       variable can be staged that way; otherwise column by column (or keys
       only in streaming mode)
    2. Sort the keys with the selected algorithm
    3. Transfer the data back to Stata in sorted order
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "ctools_rowstage.h"
#include "csort_impl.h"
#include "csort_stream.h"

/*
    Parse streaming mode option from argument string.
    Looks for "stream" or "stream(#)" to enable streaming mode,
    and "nostream" to explicitly disable auto-streaming.

    Streaming mode:
    - Only loads key (sort) variables into C memory
    - Streams permutation application to non-key variables
    - Reduces memory usage for wide datasets
    - Best when: many non-key columns, limited memory

    The optional number specifies how many variables to load at a time
    (1-16, default 1). Higher values use more memory but may be faster.

    Returns:
      -1       if explicitly disabled via "nostream"
       0       if not specified (auto-detect eligible)
       1-16    for the batch size (explicitly enabled)
*/
static int parse_stream_option(const char *args)
{
    const char *p;
    int batch_size = 1;  /* Default: process 1 variable at a time */

    /* Check for "nostream" first (must check before "stream" since
       "nostream" contains "stream" as a substring) */
    if (strstr(args, "nostream") != NULL) {
        return -1;  /* Explicitly disabled */
    }

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

    /* Find where options start - stop at any of: alg=, stream, threads, strw: */
    opts_start = strstr(args, "alg=");
    const char *stream_start = strstr(args, "stream");
    const char *threads_start = strstr(args, "threads");
    const char *strw_start = strstr(args, "strw:");

    /* Use the earliest option as the stopping point */
    if (stream_start != NULL && (opts_start == NULL || stream_start < opts_start)) {
        opts_start = stream_start;
    }
    if (threads_start != NULL && (opts_start == NULL || threads_start < opts_start)) {
        opts_start = threads_start;
    }
    if (strw_start != NULL && (opts_start == NULL || strw_start < opts_start)) {
        opts_start = strw_start;
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
    Parse string widths from args.  Format: "strw:17,0,0,8,0,..."
    Returns malloc'd array of nvars ints (0 = numeric), or NULL if not found.
    Caller must free.
*/
static int *parse_string_widths(const char *args, size_t nvars)
{
    const char *p = strstr(args, "strw:");
    if (!p) return NULL;
    p += 5;  /* skip "strw:" */

    int *widths = (int *)calloc(nvars, sizeof(int));
    if (!widths) return NULL;

    for (size_t j = 0; j < nvars && *p != '\0' && *p != ' '; j++) {
        char *endptr;
        long val = strtol(p, &endptr, 10);
        widths[j] = (int)val;
        p = endptr;
        if (*p == ',') p++;
    }

    return widths;
}

/* An identity order leaves every Stata observation untouched. Random orders
   usually differ at the first entry, so this scan adds no full pass unless
   rows are already in order (or have a long unchanged prefix). */
static int csort_order_is_identity(const perm_idx_t *order, size_t nobs)
{
    for (size_t i = 0; i < nobs; i++) {
        if (order[i] != (perm_idx_t)i) return 0;
    }
    return 1;
}

/* sort_vars holds one-based columns. Compare adjacent pairs whose right row
   lies in [start, end); start must be at least one. Every numeric comparison
   uses the engines' sortable keys, including signed zero and extended missings.
   strcmp supplies unsigned byte ordering for strings. */
static int csort_key_range_is_ordered(const stata_data *data, const int *sort_vars,
                                       size_t nsort, size_t start, size_t end)
{
    for (size_t i = start; i < end; i++) {
        for (size_t k = 0; k < nsort; k++) {
            const stata_variable *var = &data->vars[sort_vars[k] - 1];
            int cmp;
            if (var->type == STATA_TYPE_DOUBLE) {
                uint64_t left = ctools_double_to_sortable(var->data.dbl[i - 1]);
                uint64_t right = ctools_double_to_sortable(var->data.dbl[i]);
                cmp = (left > right) - (left < right);
            } else {
                const char *left = var->data.str[i - 1];
                const char *right = var->data.str[i];
                cmp = strcmp(left ? left : "", right ? right : "");
            }
            if (cmp > 0) return 0;
            if (cmp < 0) break;
        }
    }
    return 1;
}

static int csort_keys_are_ordered(const stata_data *data, const int *sort_vars, size_t nsort)
{
    /* Reject ordinary shuffled data before waking workers. Large ordered or
       nearly ordered inputs share the remaining pairs across the configured
       team, with at least the sort policy's minimum observations per worker. */
    const size_t prefix_end = data->nobs < 1025 ? data->nobs : 1025;
    if (!csort_key_range_is_ordered(data, sort_vars, nsort, 1, prefix_end)) return 0;
    if (prefix_end == data->nobs) return 1;

    /* Appended observations often contain the only inversion. Screen the last
       64 pairs before a full scan, excluding pairs the prefix already checked.
       The middle scan below excludes this tail, so every ordered pair is still
       checked exactly once even when the prefix and tail meet. */
    const size_t tail_start = data->nobs - prefix_end > 64 ? data->nobs - 64 : prefix_end;
    if (!csort_key_range_is_ordered(data, sort_vars, nsort, tail_start, data->nobs)) return 0;
    if (tail_start == prefix_end) return 1;

#ifdef _OPENMP
    const size_t remaining = tail_start - prefix_end;
    int threads = ctools_get_max_threads();
    if (remaining / MIN_OBS_PER_THREAD < (size_t)threads)
        threads = (int)(remaining / MIN_OBS_PER_THREAD);
    if (threads > 1) {
        int ordered = 1;
        #pragma omp parallel num_threads(threads) reduction(&:ordered)
        {
            const size_t team = (size_t)omp_get_num_threads();
            const size_t tid = (size_t)omp_get_thread_num();
            const size_t base = remaining / team, extra = remaining % team;
            const size_t start = prefix_end + tid * base + (tid < extra ? tid : extra);
            const size_t end = start + base + (tid < extra);
            /* Partition pair indices, rather than row ranges: the first pair
               of each worker still compares against its preceding row. */
            ordered &= csort_key_range_is_ordered(data, sort_vars, nsort, start, end);
        }
        return ordered;
    }
#endif
    return csort_key_range_is_ordered(data, sort_vars, nsort, prefix_end, tail_start);
}

/* Both standard-mode loaders initialize sort_order to identity. Only auto
   bypasses the requested engine when the input is already ordered. */
static stata_retcode csort_compute_order(stata_data *data, int *sort_vars, size_t nsort,
                                         sort_algorithm_t algorithm, int *identity)
{
    if (algorithm == SORT_ALG_AUTO && csort_keys_are_ordered(data, sort_vars, nsort)) {
        *identity = 1;
        return STATA_OK;
    }
    stata_retcode rc = ctools_sort_dispatch(data, sort_vars, nsort, algorithm);
    *identity = rc == STATA_OK && csort_order_is_identity(data->sort_order, data->nobs);
    return rc;
}

/*
    Sort through row-major staging: read every observation once, sort copies
    of the key columns, then write every observation once in sorted order.
    sort_vars holds one-based plugin variable indices and widths is indexed by
    plugin variable index minus one.

    *modified is set once Stata may have been written. A failure before that
    (a variable that cannot be staged, or no memory for the staged rows) leaves
    the dataset untouched, so the caller can use the column transport instead.
*/
static stata_retcode csort_rowstaged(const int *widths, const int *sort_vars, size_t nsort,
                                     sort_algorithm_t algorithm, size_t obs1, size_t nobs,
                                     size_t nvars, stata_timer *timer, int *modified)
{
    ctools_rowstage staged;
    stata_data keys;
    stata_retcode rc = STATA_OK;
    int identity = 0;
    int *all_vars = (int *)ctools_safe_malloc2(nvars, sizeof(int));
    int *key_vars = (int *)ctools_safe_malloc2(nsort, sizeof(int));

    *modified = 0;
    ctools_rowstage_init(&staged);
    stata_data_init(&keys);
    if (!all_vars || !key_vars) rc = STATA_ERR_MEMORY;

    if (rc == STATA_OK) {
        for (size_t j = 0; j < nvars; j++) all_vars[j] = (int)j + 1;
        timer->load_time = ctools_timer_seconds();
        rc = ctools_rowstage_load(&staged, all_vars, nvars, widths, obs1, nobs);
    }

    /* Keys are sorted as contiguous columns, exactly as the loader lays them out. */
    if (rc == STATA_OK) {
        keys.vars = (stata_variable *)ctools_safe_cacheline_alloc2(nsort, sizeof(stata_variable));
        keys.sort_order = (perm_idx_t *)ctools_safe_cacheline_alloc2(nobs, sizeof(perm_idx_t));
        if (!keys.vars || !keys.sort_order) {
            rc = STATA_ERR_MEMORY;
        } else {
            memset(keys.vars, 0, nsort * sizeof(stata_variable));
            keys.nvars = nsort;
            keys.nobs = nobs;
            for (size_t i = 0; i < nobs; i++) keys.sort_order[i] = (perm_idx_t)i;
            for (size_t k = 0; k < nsort && rc == STATA_OK; k++) {
                key_vars[k] = (int)k + 1;
                rc = ctools_rowstage_column(&staged, (size_t)sort_vars[k] - 1, &keys.vars[k]);
            }
        }
    }
    timer->load_time = ctools_timer_seconds() - timer->load_time;

    if (rc == STATA_OK) {
        timer->sort_time = ctools_timer_seconds();
        rc = csort_compute_order(&keys, key_vars, nsort, algorithm, &identity);
        timer->sort_time = ctools_timer_seconds() - timer->sort_time;
    }

    if (rc == STATA_OK && !identity) {
        timer->store_time = ctools_timer_seconds();
        *modified = 1;
        rc = ctools_rowstage_store(&staged, keys.sort_order, nobs, obs1);
        timer->store_time = ctools_timer_seconds() - timer->store_time;
    }

    stata_data_free(&keys);
    ctools_rowstage_free(&staged);
    free(key_vars);
    free(all_vars);
    return rc;
}

/*
    Main entry point for csort command.
*/
ST_retcode csort_main(const char *args)
{
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
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
    size_t obs2 = SF_in2();
    nvars = SF_nvars();

    /* Validate every column before empty/ordered-input shortcuts. Loading
       short strL text is supported, but sorting cannot publish strL rows;
       whether the permutation is identity must not change that contract. */
    for (size_t v = 1; v <= nvars; v++) {
        if (SF_var_is_string((ST_int)v) && SF_var_is_strl((ST_int)v)) {
            SF_error("csort: strL I/O is unsupported; recast to str2045 when lossless\n");
            free(sort_vars);
            return ctools_stata_rc(STATA_ERR_UNSUPPORTED_TYPE);
        }
    }

    /* Handle empty dataset case - nothing to sort, return success */
    if (obs2 < obs1) {
        free(sort_vars);
        return 0;  /* Success - nothing to do */
    }

    if (nvars == 0) {
        SF_error("csort: no variables in dataset\n");
        free(sort_vars);
        return 2000;
    }

    /* Both paths need valid keys, each listed once: sorting by (x, x) is
       sorting by x, and the streaming path counts every key as a separate
       variable. Keep the first occurrence of each. */
    size_t nunique = 0;
    for (size_t k = 0; k < nsort; k++) {
        if (sort_vars[k] < 1 || (size_t)sort_vars[k] > nvars) {
            SF_error("csort: invalid sort variable index\n");
            free(sort_vars);
            return 198;
        }
        size_t prev = 0;
        while (prev < nunique && sort_vars[prev] != sort_vars[k]) prev++;
        if (prev == nunique) sort_vars[nunique++] = sort_vars[k];
    }
    nsort = nunique;

    /* Check for streaming mode option (returns -1 if nostream, 0 if auto, 1-16 for batch size) */
    int stream_batch_size = parse_stream_option(args);

    /* Auto-detect streaming mode: enable when K > 30 and N > 10M */
    if (stream_batch_size == 0) {
        size_t nobs = obs2 - obs1 + 1;
        if (nvars > 30 && nobs > 10000000) {
            int nthreads = ctools_get_max_threads();
            stream_batch_size = nthreads > 16 ? 16 : nthreads;
            if (stream_batch_size < 1) stream_batch_size = 1;
        }
    } else if (stream_batch_size == -1) {
        /* Explicit nostream: force standard mode */
        stream_batch_size = 0;
    }

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

        /* Parse string widths for streaming mode */
        int *strwidths = parse_string_widths(args, nvars);

        /* Call streaming sort with batch size */
        rc = csort_stream_sort(sort_vars, nsort, all_var_indices, nvars,
                                algorithm, 0, stream_batch_size, strwidths,
                                &stream_timings);

        free(strwidths);
        free(all_var_indices);
        free(sort_vars);

        if (rc != STATA_OK) {
            snprintf(msg, sizeof(msg), "csort: streaming sort failed (error %d)\n", rc);
            SF_error(msg);
            return ctools_stata_rc(rc);
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

    /* Whole observations move fastest row by row. Datasets that cannot be
       staged that way (strings without a usable width) and failures
       before any write continue with the column transport below. */
    {
        int modified = 0;
        int *strwidths = parse_string_widths(args, nvars);
        rc = csort_rowstaged(strwidths, sort_vars, nsort, algorithm, obs1, obs2 - obs1 + 1,
                             nvars, &timer, &modified);
        free(strwidths);
        if (rc == STATA_OK) {
            free(sort_vars);
            timer.total_time = ctools_timer_seconds() - t_start;
            SF_scal_save("_csort_stream", 0.0);
            SF_scal_save("_csort_time_load", timer.load_time);
            SF_scal_save("_csort_time_sort", timer.sort_time);
            SF_scal_save("_csort_time_permute", 0.0);
            SF_scal_save("_csort_time_store", timer.store_time);
            SF_scal_save("_csort_time_stream", 0.0);
            SF_scal_save("_csort_time_cleanup", 0.0);
            SF_scal_save("_csort_time_total", timer.total_time);
            CTOOLS_SAVE_THREAD_INFO("_csort");
            return 0;
        }
        if (modified || (rc != STATA_ERR_UNSUPPORTED_TYPE && rc != STATA_ERR_MEMORY)) {
            snprintf(msg, sizeof(msg), "csort: sort failed (error %d)\n", rc);
            SF_error(msg);
            free(sort_vars);
            return ctools_stata_rc(rc);
        }
        timer.load_time = timer.sort_time = timer.store_time = 0.0;
    }

    /* ================================================================
       PHASE 1: Load data from Stata to C
       ================================================================ */
    timer.load_time = ctools_timer_seconds();

    /* ctools_data_load auto-detects string widths from the __ctools_strw
       Stata local (set by _ctools_strw.ado) for flat buffer optimization. */
    rc = ctools_data_load(&filtered, NULL, 0, 0, 0, CTOOLS_LOAD_SKIP_IF);

    timer.load_time = ctools_timer_seconds() - timer.load_time;

    if (rc != STATA_OK) {
        snprintf(msg, sizeof(msg), "csort: failed to load data (error %d)\n", rc);
        SF_error(msg);
        ctools_filtered_data_free(&filtered);
        free(sort_vars);
        return ctools_stata_rc(rc);
    }

    /* ================================================================
       PHASE 2: Sort the data
       ================================================================ */

    timer.sort_time = ctools_timer_seconds();
    double t_permute = 0.0;  /* Permutation time */

    /* Call unified sort dispatcher (computes sort_order only) */
    int identity = 0;
    rc = csort_compute_order(&filtered.data, sort_vars, nsort, algorithm, &identity);

    /* Record sort computation time */
    timer.sort_time = ctools_timer_seconds() - timer.sort_time;

    /* Apply the order while storing, avoiding full-column permutation buffers. */

    if (rc != STATA_OK) {
        snprintf(msg, sizeof(msg), "csort: sort failed (error %d)\n", rc);
        SF_error(msg);
        ctools_filtered_data_free(&filtered);
        free(sort_vars);
        return ctools_stata_rc(rc);
    }


    /* ================================================================
       PHASE 3: Transfer sorted data back to Stata
       ================================================================ */
    if (!identity) {
        timer.store_time = ctools_timer_seconds();
        rc = ctools_data_store_sorted(&filtered.data, obs1);
        timer.store_time = ctools_timer_seconds() - timer.store_time;
    }

    if (rc != STATA_OK) {
        snprintf(msg, sizeof(msg), "csort: failed to store data (error %d)\n", rc);
        SF_error(msg);
        ctools_filtered_data_free(&filtered);
        free(sort_vars);
        return ctools_stata_rc(rc);
    }

    /* Clean up C memory (this can be slow for large datasets!) */
    double t_cleanup = ctools_timer_seconds();
    ctools_filtered_data_free(&filtered);
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

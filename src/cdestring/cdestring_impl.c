/*
 * cdestring_impl.c
 *
 * High-performance string-to-numeric conversion for Stata
 * Part of the ctools suite
 *
 * Algorithm (3-phase parallel pipeline):
 * 1. Bulk load (parallel across variables): Use ctools_data_load_selective()
 *    to load all source string variables from Stata into C memory.
 * 2. Parse (parallel across observations): OpenMP parallel loop converts
 *    strings to doubles in pure C memory with no SPI calls.
 * 3. Bulk store (parallel across variables): Write numeric results back to
 *    Stata via SF_vstore, one thread per variable.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdbool.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_timer.h"
#include "ctools_config.h"
#include "ctools_arena.h"
#include "ctools_threads.h"
#include "cdestring_impl.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CDESTRING_STR_BUF_SIZE  2048
#define CDESTRING_MAX_IGNORE    256
#define CDESTRING_MAX_VARS      1000

/* ============================================================================
 * Option Parsing
 * ============================================================================ */

typedef struct {
    int *src_indices;       /* Source string variable indices (1-based) */
    int *dst_indices;       /* Destination numeric variable indices (1-based) */
    int nvars;              /* Number of variables to process */
    char ignore_chars[CDESTRING_MAX_IGNORE]; /* Characters to ignore */
    int ignore_len;         /* Length of ignore_chars string */
    int force;              /* 1 = convert non-numeric to missing */
    int percent;            /* 1 = strip % and divide by 100 */
    int dpcomma;            /* 1 = use comma as decimal separator */
    int verbose;            /* 1 = print timing info */
} cdestring_options;

/*
 * Check if a character should be ignored during parsing.
 * Uses a lookup table for O(1) checking.
 */
static int ignore_table[256];

static void build_ignore_table(const char *ignore_chars, int len)
{
    memset(ignore_table, 0, sizeof(ignore_table));
    for (int i = 0; i < len; i++) {
        ignore_table[(unsigned char)ignore_chars[i]] = 1;
    }
}

/*
 * Strip ignored characters from a string in-place.
 * Returns the new length.
 */
static int strip_ignored_chars(char *str, int len)
{
    int write_pos = 0;
    for (int i = 0; i < len; i++) {
        unsigned char c = (unsigned char)str[i];
        if (!ignore_table[c]) {
            str[write_pos++] = str[i];
        }
    }
    str[write_pos] = '\0';
    return write_pos;
}

/*
 * Strip percent sign from string and return 1 if found.
 */
static int strip_percent(char *str, int *len)
{
    int found = 0;
    int write_pos = 0;
    for (int i = 0; i < *len; i++) {
        if (str[i] == '%') {
            found = 1;
        } else {
            str[write_pos++] = str[i];
        }
    }
    str[write_pos] = '\0';
    *len = write_pos;
    return found;
}

/*
 * Parse variable indices from command string.
 * Format: "src1 dst1 src2 dst2 ... nvars=N [options]"
 *
 * Returns:
 *   0 on success
 *  -1 on syntax error (invalid arguments)
 *  -2 on memory allocation failure
 */
static int parse_options(const char *args, cdestring_options *opts)
{
    memset(opts, 0, sizeof(cdestring_options));

    if (args == NULL || *args == '\0') {
        return -1;
    }

    /* Make a working copy */
    char *args_copy = strdup(args);
    if (!args_copy) return -2;  /* Memory allocation failure */

    /* First, find nvars= to know how many variables */
    const char *nvars_ptr = strstr(args, "nvars=");
    if (nvars_ptr == NULL) {
        free(args_copy);
        return -1;
    }

    /* Extract just the number part (until whitespace or end) */
    const char *num_start = nvars_ptr + 6;
    char nvars_buf[32];
    int i = 0;
    while (num_start[i] && num_start[i] != ' ' && num_start[i] != '\t' && i < 31) {
        nvars_buf[i] = num_start[i];
        i++;
    }
    nvars_buf[i] = '\0';

    if (!ctools_safe_atoi(nvars_buf, &opts->nvars)) {
        free(args_copy);
        return -1;  /* Invalid nvars value */
    }
    if (opts->nvars <= 0 || opts->nvars > CDESTRING_MAX_VARS) {
        free(args_copy);
        return -1;
    }

    /* Allocate index arrays */
    opts->src_indices = malloc(opts->nvars * sizeof(int));
    opts->dst_indices = malloc(opts->nvars * sizeof(int));
    if (!opts->src_indices || !opts->dst_indices) {
        free(opts->src_indices);
        free(opts->dst_indices);
        free(args_copy);
        return -2;  /* Memory allocation failure */
    }

    /* Parse variable indices (pairs of src, dst) */
    char *saveptr = NULL;
    char *token = strtok_r(args_copy, " \t", &saveptr);
    int idx = 0;

    while (token != NULL && idx < opts->nvars * 2) {
        /* Stop when we hit an option keyword */
        if (strncmp(token, "nvars=", 6) == 0 ||
            strncmp(token, "ignore=", 7) == 0 ||
            strcmp(token, "force") == 0 ||
            strcmp(token, "percent") == 0 ||
            strcmp(token, "dpcomma") == 0 ||
            strcmp(token, "verbose") == 0) {
            break;
        }

        int val;
        if (!ctools_safe_atoi(token, &val)) {
            /* Invalid number - skip this token */
            token = strtok_r(NULL, " \t", &saveptr);
            continue;
        }
        if (val > 0) {
            if (idx % 2 == 0) {
                opts->src_indices[idx / 2] = val;
            } else {
                opts->dst_indices[idx / 2] = val;
            }
            idx++;
        }
        token = strtok_r(NULL, " \t", &saveptr);
    }

    free(args_copy);

    if (idx < opts->nvars * 2) {
        /* Not enough indices provided */
        free(opts->src_indices);
        free(opts->dst_indices);
        return -1;
    }

    /* Parse other options from original args */

    /* ignore= option */
    const char *ignore_ptr = strstr(args, "ignore=");
    if (ignore_ptr != NULL) {
        ignore_ptr += 7;
        int i = 0;
        while (*ignore_ptr && *ignore_ptr != ' ' && *ignore_ptr != '\t' &&
               i < CDESTRING_MAX_IGNORE - 1) {
            /* Handle escape sequences */
            if (*ignore_ptr == '\\' && *(ignore_ptr + 1)) {
                ignore_ptr++;
                switch (*ignore_ptr) {
                    case 'n': opts->ignore_chars[i++] = '\n'; break;
                    case 't': opts->ignore_chars[i++] = '\t'; break;
                    case 'r': opts->ignore_chars[i++] = '\r'; break;
                    case 's': opts->ignore_chars[i++] = ' '; break;
                    case '\\': opts->ignore_chars[i++] = '\\'; break;
                    default: opts->ignore_chars[i++] = *ignore_ptr; break;
                }
            } else {
                opts->ignore_chars[i++] = *ignore_ptr;
            }
            ignore_ptr++;
        }
        opts->ignore_chars[i] = '\0';
        opts->ignore_len = i;
    }

    /* Boolean options */
    opts->force = (strstr(args, "force") != NULL) ? 1 : 0;
    opts->percent = (strstr(args, "percent") != NULL) ? 1 : 0;
    opts->dpcomma = (strstr(args, "dpcomma") != NULL) ? 1 : 0;
    opts->verbose = (strstr(args, "verbose") != NULL) ? 1 : 0;

    return 0;
}

static void free_options(cdestring_options *opts)
{
    free(opts->src_indices);
    free(opts->dst_indices);
    opts->src_indices = NULL;
    opts->dst_indices = NULL;
}

/* ============================================================================
 * Public Entry Point
 * ============================================================================ */

ST_retcode cdestring_main(const char *args)
{
    double t_start, t_parse, t_load, t_convert, t_store, t_total;
    cdestring_options opts;
    ST_retcode rc = 0;
    int total_converted = 0;
    int total_failed = 0;
    int had_nonnumeric = 0;
    int force_flag;

    t_start = ctools_timer_seconds();

    /* Parse arguments */
    int parse_rc = parse_options(args, &opts);
    if (parse_rc != 0) {
        if (parse_rc == -2) {
            SF_error("cdestring: memory allocation failed\n");
            return 920;
        }
        SF_error("cdestring: invalid arguments\n");
        return 198;
    }

    /* Save force flag before we potentially free opts */
    force_flag = opts.force;

    /* Build ignore character lookup table */
    if (opts.ignore_len > 0) {
        build_ignore_table(opts.ignore_chars, opts.ignore_len);
    }

    t_parse = ctools_timer_seconds();

    /* Get observation range */
    ST_int obs1 = SF_in1();
    ST_int obs2 = SF_in2();
    size_t nobs = (size_t)(obs2 - obs1 + 1);

    if (nobs == 0) {
        free_options(&opts);
        SF_scal_save("_cdestring_n_converted", 0);
        SF_scal_save("_cdestring_n_failed", 0);
        return 0;
    }

    int nvars = opts.nvars;

    /* Determine decimal and group separators */
    char dec_sep = opts.dpcomma ? ',' : '.';
    char grp_sep = opts.dpcomma ? '.' : '\0';

    /* ====================================================================
     * Phase 1: Bulk load source string variables (parallel across vars)
     * ==================================================================== */

    stata_data data;
    stata_retcode load_rc = ctools_data_load_selective(&data, opts.src_indices,
                                                        (size_t)nvars, (size_t)obs1, (size_t)obs2);
    if (load_rc != STATA_OK) {
        free_options(&opts);
        SF_error("cdestring: failed to load string data\n");
        return 920;
    }

    t_load = ctools_timer_seconds();

    /* ====================================================================
     * Phase 2: Parse strings to doubles (parallel across observations)
     * ==================================================================== */

    /* Build if-condition mask (sequential — SPI calls) */
    unsigned char *if_mask = malloc(nobs);
    if (!if_mask) {
        stata_data_free(&data);
        free_options(&opts);
        SF_error("cdestring: memory allocation failed\n");
        return 920;
    }
    for (size_t i = 0; i < nobs; i++) {
        if_mask[i] = SF_ifobs((ST_int)(obs1 + (ST_int)i)) ? 1 : 0;
    }

    /* Allocate result arrays: one double array per variable */
    double **results = malloc((size_t)nvars * sizeof(double *));
    if (!results) {
        stata_data_free(&data);
        free_options(&opts);
        SF_error("cdestring: memory allocation failed\n");
        return 920;
    }

    for (int v = 0; v < nvars; v++) {
        results[v] = malloc(nobs * sizeof(double));
        if (!results[v]) {
            for (int j = 0; j < v; j++) free(results[j]);
            free(results);
            stata_data_free(&data);
            free_options(&opts);
            SF_error("cdestring: memory allocation failed\n");
            return 920;
        }
    }

    /* Per-variable counters */
    int *var_converted = calloc(nvars, sizeof(int));
    int *var_failed = calloc(nvars, sizeof(int));
    if (!var_converted || !var_failed) {
        for (int v = 0; v < nvars; v++) free(results[v]);
        free(results);
        free(var_converted);
        free(var_failed);
        stata_data_free(&data);
        free_options(&opts);
        SF_error("cdestring: memory allocation failed\n");
        return 920;
    }

    /* Process each variable: parallel over observations */
    for (int v = 0; v < nvars; v++) {
        char **str_data = data.vars[v].data.str;
        double *res = results[v];
        int local_converted = 0;
        int local_failed = 0;
        int local_nonnumeric = 0;

        #pragma omp parallel for reduction(+:local_converted,local_failed,local_nonnumeric) schedule(static)
        for (size_t i = 0; i < nobs; i++) {
            /* Skip observations that don't meet the if/in condition */
            if (!if_mask[i]) {
                res[i] = SV_missval;
                continue;
            }

            const char *src = str_data[i];

            /* Empty or NULL string -> missing */
            if (src == NULL || src[0] == '\0') {
                res[i] = SV_missval;
                local_converted++;
                continue;
            }

            /* Work on a thread-local copy for in-place modifications */
            char work_buf[CDESTRING_STR_BUF_SIZE + 1];
            int len = (int)strlen(src);
            if (len > CDESTRING_STR_BUF_SIZE) len = CDESTRING_STR_BUF_SIZE;
            memcpy(work_buf, src, len);
            work_buf[len] = '\0';

            /* Strip ignored characters */
            if (opts.ignore_len > 0) {
                len = strip_ignored_chars(work_buf, len);
            }

            /* Strip percent sign */
            int has_percent = 0;
            if (opts.percent) {
                has_percent = strip_percent(work_buf, &len);
            }

            /* Empty after stripping -> missing */
            if (len == 0) {
                res[i] = SV_missval;
                local_converted++;
                continue;
            }

            /* Parse the number */
            double val;
            int parsed = ctools_parse_double_with_separators(
                work_buf, len, &val, SV_missval, dec_sep, grp_sep);

            if (parsed) {
                if (opts.percent && has_percent && val != SV_missval) {
                    val /= 100.0;
                }
                res[i] = val;
                local_converted++;
            } else {
                res[i] = SV_missval;
                if (opts.force) {
                    local_converted++;
                } else {
                    local_failed++;
                    local_nonnumeric = 1;
                }
            }
        }

        var_converted[v] = local_converted;
        var_failed[v] = local_failed;
        if (local_nonnumeric) had_nonnumeric = 1;
    }

    /* Free loaded string data and if-mask — no longer needed */
    stata_data_free(&data);
    free(if_mask);

    t_convert = ctools_timer_seconds();

    /* ====================================================================
     * Phase 3: Store results back to Stata (parallel across variables)
     * ==================================================================== */

    /* Each thread writes one variable's results via SF_vstore.
     * Different threads access different variable columns — SPI-safe. */
    #pragma omp parallel for schedule(static)
    for (int v = 0; v < nvars; v++) {
        int dst_idx = opts.dst_indices[v];
        double *res = results[v];
        for (size_t i = 0; i < nobs; i++) {
            SF_vstore(dst_idx, (ST_int)(obs1 + (ST_int)i), res[i]);
        }
    }

    t_store = ctools_timer_seconds();

    /* Aggregate results */
    for (int v = 0; v < nvars; v++) {
        total_converted += var_converted[v];
        total_failed += var_failed[v];
    }

    /* Free result arrays */
    for (int v = 0; v < nvars; v++) free(results[v]);
    free(results);
    free(var_converted);
    free(var_failed);

    t_total = t_store - t_start;

    /* Store results */
    SF_scal_save("_cdestring_n_converted", (double)total_converted);
    SF_scal_save("_cdestring_n_failed", (double)total_failed);
    SF_scal_save("_cdestring_time_parse", t_parse - t_start);
    SF_scal_save("_cdestring_time_load", t_load - t_parse);
    SF_scal_save("_cdestring_time_convert", t_convert - t_load);
    SF_scal_save("_cdestring_time_store", t_store - t_convert);
    SF_scal_save("_cdestring_time_total", t_total);

#ifdef _OPENMP
    SF_scal_save("_cdestring_openmp_enabled", 1.0);
    SF_scal_save("_cdestring_threads_max", (double)omp_get_max_threads());
#else
    SF_scal_save("_cdestring_openmp_enabled", 0.0);
    SF_scal_save("_cdestring_threads_max", 1.0);
#endif

    free_options(&opts);

    /* If non-numeric values found and force not specified, let ado handle message */
    if (had_nonnumeric && !force_flag) {
        rc = 0;
    }

    return rc;
}

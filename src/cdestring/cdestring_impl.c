/*
 * cdestring_impl.c
 *
 * High-performance string-to-numeric conversion for Stata
 * Part of the ctools suite
 *
 * Algorithm (3-phase parallel pipeline):
 * 1. Bulk load (parallel across variables): Use ctools_data_load()
 *    to load all source string variables from Stata into C memory.
 * 2. Parse (parallel across observations): OpenMP parallel loop converts
 *    strings to doubles in pure C memory with no SPI calls.
 * 3. Bulk store (parallel across variables): Write numeric results back to
 *    Stata via ctools_store_filtered, one variable at a time.
 */

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_runtime.h"
#include "ctools_config.h"
#include "ctools_parse.h"
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
static unsigned char ignore_table[256];

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
 * Parse the ignore= option value, handling escape sequences.
 * This is kept as a local helper because escape sequences like \s, \n, \t
 * are specific to cdestring and not handled by ctools_parse_string_option.
 */
static int parse_ignore_option(const char *args, char *ignore_chars, int max_len)
{
    const char *ignore_ptr = strstr(args, "ignore=");
    if (!ignore_ptr) return 0;

    ignore_ptr += 7;
    int i = 0;
    while (*ignore_ptr && *ignore_ptr != ' ' && *ignore_ptr != '\t' &&
           i < max_len - 1) {
        if (*ignore_ptr == '\\' && *(ignore_ptr + 1)) {
            ignore_ptr++;
            switch (*ignore_ptr) {
                case 'n': ignore_chars[i++] = '\n'; break;
                case 't': ignore_chars[i++] = '\t'; break;
                case 'r': ignore_chars[i++] = '\r'; break;
                case 's': ignore_chars[i++] = ' '; break;
                case '\\': ignore_chars[i++] = '\\'; break;
                default: ignore_chars[i++] = *ignore_ptr; break;
            }
        } else {
            ignore_chars[i++] = *ignore_ptr;
        }
        ignore_ptr++;
    }
    ignore_chars[i] = '\0';
    return i;
}

/*
 * Parse all options from args string.
 * Returns 0 on success, -1 on syntax error, -2 on memory failure.
 */
static int parse_options(const char *args, cdestring_options *opts)
{
    memset(opts, 0, sizeof(cdestring_options));

    if (args == NULL || *args == '\0') {
        return -1;
    }

    /* Parse nvars= */
    if (!ctools_parse_int_option(args, "nvars", &opts->nvars)) {
        return -1;
    }
    if (opts->nvars <= 0 || opts->nvars > CDESTRING_MAX_VARS) {
        return -1;
    }

    /* Allocate index arrays */
    opts->src_indices = malloc(opts->nvars * sizeof(int));
    opts->dst_indices = malloc(opts->nvars * sizeof(int));
    if (!opts->src_indices || !opts->dst_indices) {
        free(opts->src_indices);
        free(opts->dst_indices);
        return -2;
    }

    /* Parse variable index pairs from start of args (stack-allocated) */
    const char *cursor = args;
    int indices[CDESTRING_MAX_VARS * 2];

    if (ctools_parse_int_array(indices, (size_t)(opts->nvars * 2), &cursor) != 0) {
        free(opts->src_indices);
        free(opts->dst_indices);
        return -1;
    }

    for (int i = 0; i < opts->nvars; i++) {
        opts->src_indices[i] = indices[i * 2];
        opts->dst_indices[i] = indices[i * 2 + 1];
    }

    /* Parse ignore= with escape sequence handling */
    opts->ignore_len = parse_ignore_option(args, opts->ignore_chars,
                                            CDESTRING_MAX_IGNORE);

    /* Boolean options */
    opts->force = ctools_parse_bool_option(args, "force");
    opts->percent = ctools_parse_bool_option(args, "percent");
    opts->dpcomma = ctools_parse_bool_option(args, "dpcomma");
    opts->verbose = ctools_parse_bool_option(args, "verbose");

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
    int total_converted = 0;
    int total_failed = 0;

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



    /* Build ignore character lookup table */
    if (opts.ignore_len > 0) {
        build_ignore_table(opts.ignore_chars, opts.ignore_len);
    }

    t_parse = ctools_timer_seconds();

    int nvars = opts.nvars;

    /* Determine decimal and group separators */
    char dec_sep = opts.dpcomma ? ',' : '.';
    char grp_sep = opts.dpcomma ? '.' : '\0';

    /* ====================================================================
     * Phase 1: Bulk load source string variables with if/in filtering
     * Uses ctools_data_load() to load only filtered observations.
     * ==================================================================== */

    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    stata_retcode load_rc = ctools_data_load(&filtered, opts.src_indices,
                                                       (size_t)nvars, 0, 0, 0);
    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        free_options(&opts);
        SF_error("cdestring: failed to load string data\n");
        return 920;
    }

    size_t nobs = filtered.data.nobs;
    perm_idx_t *obs_map = filtered.obs_map;

    if (nobs == 0) {
        ctools_filtered_data_free(&filtered);
        free_options(&opts);
        SF_scal_save("_cdestring_n_converted", 0);
        SF_scal_save("_cdestring_n_failed", 0);
        return 0;
    }

    t_load = ctools_timer_seconds();

    /* ====================================================================
     * Phase 2: Parse strings to doubles (parallel across observations)
     * Note: No if_mask needed since data is pre-filtered.
     * ==================================================================== */

    /* Allocate result arrays: single contiguous block for all variables */
    double *results_block = malloc((size_t)nvars * nobs * sizeof(double));
    double *results_ptrs[CDESTRING_MAX_VARS];
    double **results = results_ptrs;
    if (!results_block) {
        ctools_filtered_data_free(&filtered);
        free_options(&opts);
        SF_error("cdestring: memory allocation failed\n");
        return 920;
    }

    for (int v = 0; v < nvars; v++) {
        results[v] = results_block + (size_t)v * nobs;
    }

    /* Per-variable counters (stack-allocated, nvars <= CDESTRING_MAX_VARS = 1000) */
    int var_converted[CDESTRING_MAX_VARS];
    int var_failed[CDESTRING_MAX_VARS];
    memset(var_converted, 0, nvars * sizeof(int));
    memset(var_failed, 0, nvars * sizeof(int));

    /* Process each variable: parallel over observations.
     * Note: All observations are valid (passed if/in filter). */
    for (int v = 0; v < nvars; v++) {
        char **str_data = filtered.data.vars[v].data.str;
        double *res = results[v];
        int local_converted = 0;
        int local_failed = 0;

        #pragma omp parallel for reduction(+:local_converted,local_failed) schedule(static)
        for (size_t i = 0; i < nobs; i++) {
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
                }
            }
        }

        var_converted[v] = local_converted;
        var_failed[v] = local_failed;
    }

    t_convert = ctools_timer_seconds();

    /* ====================================================================
     * Phase 3: Store results back to Stata via ctools_store_filtered
     * Uses obs_map to write to correct Stata observations.
     * ==================================================================== */

    #pragma omp parallel for schedule(static)
    for (int v = 0; v < nvars; v++) {
        ctools_store_filtered(results[v], nobs, opts.dst_indices[v], obs_map);
    }

    /* Free loaded string data — no longer needed after store */
    ctools_filtered_data_free(&filtered);

    t_store = ctools_timer_seconds();

    /* Aggregate results */
    for (int v = 0; v < nvars; v++) {
        total_converted += var_converted[v];
        total_failed += var_failed[v];
    }

    /* Free contiguous results block (results pointer array is stack-allocated) */
    free(results_block);
    /* var_converted and var_failed are stack-allocated */

    t_total = t_store - t_start;

    /* Store results */
    SF_scal_save("_cdestring_n_converted", (double)total_converted);
    SF_scal_save("_cdestring_n_failed", (double)total_failed);
    SF_scal_save("_cdestring_time_parse", t_parse - t_start);
    SF_scal_save("_cdestring_time_load", t_load - t_parse);
    SF_scal_save("_cdestring_time_convert", t_convert - t_load);
    SF_scal_save("_cdestring_time_store", t_store - t_convert);
    SF_scal_save("_cdestring_time_total", t_total);

    CTOOLS_SAVE_THREAD_INFO("_cdestring");

    free_options(&opts);

    return 0;
}

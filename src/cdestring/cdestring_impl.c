/*
 * cdestring_impl.c
 *
 * High-performance string-to-numeric conversion for Stata
 * Part of the ctools suite
 *
 * Algorithm:
 * 1. Parse command arguments (variable indices, options)
 * 2. For each variable (parallelized across variables):
 *    a. Read all string values from Stata into C memory
 *    b. For each observation (sequential):
 *       - Strip ignore characters if specified
 *       - Strip % and divide by 100 if percent option
 *       - Parse string to double using fast parser
 *       - Store result in C array
 *    c. Write numeric values back to Stata (sequential)
 *
 * Parallelization strategy:
 * - Processing is parallelized across variables (not observations)
 * - This avoids concurrent access to Stata's SPI for the same variable
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

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
    opts->nvars = atoi(nvars_ptr + 6);
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

        int val = atoi(token);
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
 * Main Implementation
 * ============================================================================ */

/*
 * Process a single variable: convert string to numeric.
 * This function runs entirely sequentially for Stata SPI safety.
 *
 * Returns:
 *   0 = success
 *   1 = had non-numeric values (only returned if force=0)
 *  -1 = error
 */
static int process_variable(int src_idx, int dst_idx,
                            ST_int obs1, size_t nobs,
                            cdestring_options *opts,
                            int *nconverted, int *nfailed)
{
    int local_converted = 0;
    int local_failed = 0;
    int had_nonnumeric = 0;

    /* Determine decimal and group separators */
    char dec_sep = opts->dpcomma ? ',' : '.';
    char grp_sep = opts->dpcomma ? '.' : '\0';

    /* Check if source is a string variable */
    if (!SF_var_is_string(src_idx)) {
        return -1;
    }

    int is_strl = SF_var_is_strl(src_idx);

    /* Allocate work buffer */
    char *buf = malloc(CDESTRING_STR_BUF_SIZE + 1);
    char *work_buf = malloc(CDESTRING_STR_BUF_SIZE + 1);
    char *large_buf = NULL;
    size_t large_buf_size = 0;

    if (!buf || !work_buf) {
        free(buf);
        free(work_buf);
        return -1;
    }

    /* Process each observation sequentially */
    for (size_t i = 0; i < nobs; i++) {
        ST_int obs = obs1 + (ST_int)i;

        /* Skip observations that don't meet the if/in condition */
        if (!SF_ifobs(obs)) {
            SF_vstore(dst_idx, obs, SV_missval);
            continue;
        }

        int slen = SF_sdatalen(src_idx, obs);

        /* Empty string -> missing */
        if (slen <= 0) {
            SF_vstore(dst_idx, obs, SV_missval);
            local_converted++;
            continue;
        }

        /* Get the string */
        char *str_ptr;
        if ((size_t)slen <= CDESTRING_STR_BUF_SIZE) {
            if (is_strl) {
                SF_strldata(src_idx, obs, buf, slen + 1);
            } else {
                SF_sdata(src_idx, obs, buf);
            }
            buf[slen] = '\0';
            str_ptr = buf;
        } else {
            /* Need larger buffer */
            if ((size_t)slen > large_buf_size) {
                free(large_buf);
                large_buf_size = (size_t)slen + 256;
                large_buf = malloc(large_buf_size + 1);
                if (!large_buf) {
                    SF_vstore(dst_idx, obs, SV_missval);
                    local_failed++;
                    had_nonnumeric = 1;
                    continue;
                }
            }
            if (is_strl) {
                SF_strldata(src_idx, obs, large_buf, slen + 1);
            } else {
                SF_sdata(src_idx, obs, large_buf);
            }
            large_buf[slen] = '\0';
            str_ptr = large_buf;
        }

        /* Copy to work buffer for modification */
        int len = slen;
        if ((size_t)len < CDESTRING_STR_BUF_SIZE) {
            memcpy(work_buf, str_ptr, len + 1);
            str_ptr = work_buf;
        }

        /* Strip ignored characters */
        if (opts->ignore_len > 0) {
            len = strip_ignored_chars(str_ptr, len);
        }

        /* Strip percent sign */
        int has_percent = 0;
        if (opts->percent) {
            has_percent = strip_percent(str_ptr, &len);
        }

        /* Parse the number */
        double result;
        int parsed = ctools_parse_double_with_separators(
            str_ptr, len, &result, SV_missval, dec_sep, grp_sep);

        if (parsed) {
            /* Apply percent transformation */
            if (opts->percent && has_percent && result != SV_missval) {
                result /= 100.0;
            }
            SF_vstore(dst_idx, obs, result);
            local_converted++;
        } else {
            /* Parse failed */
            SF_vstore(dst_idx, obs, SV_missval);
            if (opts->force) {
                local_converted++;
            } else {
                local_failed++;
                had_nonnumeric = 1;
            }
        }
    }

    free(large_buf);
    free(buf);
    free(work_buf);

    *nconverted = local_converted;
    *nfailed = local_failed;

    return had_nonnumeric ? 1 : 0;
}

/* ============================================================================
 * Public Entry Point
 * ============================================================================ */

ST_retcode cdestring_main(const char *args)
{
    double t_start, t_parse, t_convert, t_total;
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

    /* Allocate per-variable result arrays for OpenMP */
    int *var_converted = calloc(opts.nvars, sizeof(int));
    int *var_failed = calloc(opts.nvars, sizeof(int));
    int *var_result = calloc(opts.nvars, sizeof(int));

    if (!var_converted || !var_failed || !var_result) {
        free(var_converted);
        free(var_failed);
        free(var_result);
        free_options(&opts);
        SF_error("cdestring: memory allocation failed\n");
        return 920;
    }

    /* Process variables in parallel */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) if(opts.nvars > 1)
    #endif
    for (int v = 0; v < opts.nvars; v++) {
        int src_idx = opts.src_indices[v];
        int dst_idx = opts.dst_indices[v];

        var_result[v] = process_variable(src_idx, dst_idx, obs1, nobs,
                                          &opts, &var_converted[v], &var_failed[v]);
    }

    /* Aggregate results */
    for (int v = 0; v < opts.nvars; v++) {
        if (var_result[v] < 0) {
            rc = 920;
        }
        if (var_result[v] == 1) {
            had_nonnumeric = 1;
        }
        total_converted += var_converted[v];
        total_failed += var_failed[v];
    }

    free(var_converted);
    free(var_failed);
    free(var_result);

    if (rc != 0) {
        free_options(&opts);
        return rc;
    }

    t_convert = ctools_timer_seconds();
    t_total = t_convert - t_start;

    /* Store results */
    SF_scal_save("_cdestring_n_converted", (double)total_converted);
    SF_scal_save("_cdestring_n_failed", (double)total_failed);
    SF_scal_save("_cdestring_time_parse", t_parse - t_start);
    SF_scal_save("_cdestring_time_convert", t_convert - t_parse);
    SF_scal_save("_cdestring_time_total", t_total);

#ifdef _OPENMP
    SF_scal_save("_cdestring_openmp_enabled", 1.0);
    SF_scal_save("_cdestring_threads_max", (double)omp_get_max_threads());
#else
    SF_scal_save("_cdestring_openmp_enabled", 0.0);
    SF_scal_save("_cdestring_threads_max", 1.0);
#endif

    free_options(&opts);

    /* If non-numeric values found and force not specified, return error */
    if (had_nonnumeric && !force_flag) {
        /* Note: We've already stored missing values, the ado file will
           report which variables had issues */
        rc = 0;  /* Don't error, let ado file handle the message */
    }

    return rc;
}

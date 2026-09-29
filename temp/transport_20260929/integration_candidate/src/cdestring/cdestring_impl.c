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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

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
#define CDESTRING_MAX_FLOAT     0x1.fffffep+126  /* maxfloat(); larger is missing in a float */

/* ============================================================================
 * Option Parsing
 * ============================================================================ */

typedef struct {
    int *src_indices;       /* Source string variable indices (1-based) */
    int *dst_indices;       /* Destination numeric variable indices (1-based) */
    int nvars;              /* Number of variables to process */
    char ignore_chars[CDESTRING_MAX_IGNORE]; /* Characters to ignore */
    int ignore_len;         /* Length of ignore_chars string */
    int strip;              /* 1 = remove ignore_table bytes (ignore() and, with percent, "%") */
    int force;              /* 1 = keep conversions of variables with nonnumeric values */
    int percent;            /* 1 = strip %; divide the variable by 100 if any value had % */
    int dpcomma;            /* 1 = first comma is the decimal point; any "." is nonnumeric */
    int float_out;          /* 1 = outputs are float (range check, percent rounding) */
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

/* ASCII whitespace removed by ustrtrim() */
static inline int is_trim_space(unsigned char c)
{
    return c == ' ' || (c >= '\t' && c <= '\r');
}

/*
 * real() also reads a d or D exponent ("1d3" is 1000, "1.5D-2" is .015):
 * rewrite it as e when t[0..len) is exactly [+-]digits[.digits](d|D)[+-]digits.
 */
static int d_exponent_to_e(char *t, int len)
{
    int i = 0, digits = 0, exp_digits = 0;
    if (i < len && (t[i] == '+' || t[i] == '-')) i++;
    while (i < len && t[i] >= '0' && t[i] <= '9') { i++; digits++; }
    if (i < len && t[i] == '.') {
        i++;
        while (i < len && t[i] >= '0' && t[i] <= '9') { i++; digits++; }
    }
    if (!digits || i >= len || (t[i] != 'd' && t[i] != 'D')) return 0;
    int d = i++;
    if (i < len && (t[i] == '+' || t[i] == '-')) i++;
    while (i < len && t[i] >= '0' && t[i] <= '9') { i++; exp_digits++; }
    if (!exp_digits || i != len) return 0;
    t[d] = 'e';
    return 1;
}

/*
 * Convert one string as native destring does: remove ignore() bytes (and "%"
 * with percent), trim, keep "", "." and ".a"-".z" as missing values, apply
 * dpcomma, then parse like real(). Returns 0 for a nonnumeric value, which is
 * stored as missing; without force the ado then leaves the variable untouched.
 */
static int convert_value(const char *src, const cdestring_options *opts,
                         double missval, double *out)
{
    char buf[CDESTRING_STR_BUF_SIZE + 1];
    const char *s = src;
    int len = (int)strlen(src);

    if (opts->strip || opts->dpcomma) {
        int n = 0;
        for (int i = 0; i < len && n < CDESTRING_STR_BUF_SIZE; i++) {
            if (!ignore_table[(unsigned char)src[i]]) buf[n++] = src[i];
        }
        s = buf;
        len = n;
    }
    while (len > 0 && is_trim_space((unsigned char)s[0])) { s++; len--; }
    while (len > 0 && is_trim_space((unsigned char)s[len - 1])) len--;

    *out = missval;
    if (len == 0 || (len == 1 && s[0] == '.')) return 1;
    if (len == 2 && s[0] == '.' && s[1] >= 'a' && s[1] <= 'z') {
        /* Extended missing .a-.z: the bit pattern of . plus (letter << 40) */
        uint64_t bits;
        memcpy(&bits, &missval, sizeof(bits));
        bits += (uint64_t)(s[1] - 'a' + 1) << 40;
        memcpy(out, &bits, sizeof(bits));
        return 1;
    }

    if (opts->dpcomma) {
        /* Any "." is nonnumeric; the first "," is the decimal point unless
         * the value is "," or a "," precedes a lowercase letter. */
        int comma = -1, before_letter = 0;
        for (int i = 0; i < len; i++) {
            if (s[i] == '.') return 0;
            if (s[i] != ',') continue;
            if (comma < 0) comma = i;
            if (i + 1 < len && s[i + 1] >= 'a' && s[i + 1] <= 'z') before_letter = 1;
        }
        if (comma >= 0 && len > 1 && !before_letter) buf[(s - buf) + comma] = '.';
    }

    /* The parser returns infinite results, numbers at or above maxdouble and
     * "NA"/"NaN" as missing, which destring counts as nonnumeric. Finite
     * values below -maxdouble stay numbers, as real() keeps them. */
    double val;
    if (!ctools_parse_double_fast(s, len, &val, missval)) {
        char dexp[CDESTRING_STR_BUF_SIZE + 1];
        if (len > CDESTRING_STR_BUF_SIZE) return 0;
        memcpy(dexp, s, (size_t)len);
        if (!d_exponent_to_e(dexp, len) ||
            !ctools_parse_double_fast(dexp, len, &val, missval))
            return 0;
    }
    if (!(val < missval))
        return 0;
    if (opts->float_out && fabs(val) > CDESTRING_MAX_FLOAT) return 0;
    *out = val;
    return 1;
}

/*
 * Save the per-variable nonnumeric counts to local cdestring_failed.
 */
static ST_retcode save_failed_counts(const int *failed, int nvars)
{
    char buf[CDESTRING_MAX_VARS * 12 + 1];
    size_t used = 0;
    buf[0] = '\0';
    for (int v = 0; v < nvars; v++) {
        used += (size_t)snprintf(buf + used, sizeof(buf) - used, "%s%d",
                                 v ? " " : "", failed[v]);
    }
    return SF_macro_save("_cdestring_failed", buf);
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
    if (ctools_parse_int_option(args, "nvars", &opts->nvars) != CTOOLS_PARSE_OK) {
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
    opts->float_out = ctools_parse_bool_option(args, "float");
    opts->verbose = ctools_parse_bool_option(args, "verbose");
    opts->strip = (opts->ignore_len > 0 || opts->percent);

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

    /* Build ignore character lookup table; percent implies ignore("%") */
    build_ignore_table(opts.ignore_chars, opts.ignore_len);
    if (opts.percent) ignore_table['%'] = 1;

    t_parse = ctools_timer_seconds();

    int nvars = opts.nvars;
    const double missval = SV_missval;

    /* Per-variable counters (stack-allocated, nvars <= CDESTRING_MAX_VARS = 1000) */
    int var_converted[CDESTRING_MAX_VARS];
    int var_failed[CDESTRING_MAX_VARS];
    int var_percent[CDESTRING_MAX_VARS];
    memset(var_converted, 0, nvars * sizeof(int));
    memset(var_failed, 0, nvars * sizeof(int));
    memset(var_percent, 0, nvars * sizeof(int));

    /* ====================================================================
     * Phase 1: Bulk load source string variables with if/in filtering
     * Uses ctools_data_load() to load only filtered observations.
     * ==================================================================== */

    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    stata_retcode load_rc = ctools_data_load(&filtered, opts.src_indices,
                                              (size_t)nvars, 0, 0, CTOOLS_LOAD_CHECK_IF | CTOOLS_LOAD_NO_SORT_ORDER);
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
        return save_failed_counts(var_failed, nvars);
    }

    t_load = ctools_timer_seconds();

    /* ====================================================================
     * Phase 2: Parse strings to doubles (parallel across observations)
     * Note: No if_mask needed since data is pre-filtered.
     * ==================================================================== */

    /* Allocate result arrays: single contiguous block for all variables */
    double *results_block = ctools_safe_malloc3((size_t)nvars, nobs, sizeof(double));
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

    /* Transposed loop: single OMP region with outer loop over observations,
     * inner loop over variables. This gives one fork/join barrier instead of
     * nvars barriers, and keeps per-observation data cache-hot across variables. */
    #pragma omp parallel
    {
        /* Thread-local per-variable counters */
        int tl_converted[CDESTRING_MAX_VARS];
        int tl_failed[CDESTRING_MAX_VARS];
        int tl_percent[CDESTRING_MAX_VARS];
        memset(tl_converted, 0, nvars * sizeof(int));
        memset(tl_failed, 0, nvars * sizeof(int));
        memset(tl_percent, 0, nvars * sizeof(int));

        #pragma omp for schedule(static)
        for (size_t i = 0; i < nobs; i++) {
            for (int v = 0; v < nvars; v++) {
                const char *src = filtered.data.vars[v].data.str[i];
                if (src == NULL) src = "";

                if (opts.percent && strchr(src, '%')) tl_percent[v]++;
                if (convert_value(src, &opts, missval, &results[v][i])) {
                    tl_converted[v]++;
                } else {
                    tl_failed[v]++;
                }
            }
        }

        /* Aggregate thread-local counters */
        #pragma omp critical
        {
            for (int v = 0; v < nvars; v++) {
                var_converted[v] += tl_converted[v];
                var_failed[v] += tl_failed[v];
                var_percent[v] += tl_percent[v];
            }
        }
    }

    /* percent: like native, divide the whole variable by 100 when any value
     * had "%" (after rounding to float under float); ./100 and .a/100 are . */
    int any_percent = 0;
    for (int v = 0; v < nvars; v++) any_percent |= (var_percent[v] > 0);
    if (any_percent) {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < nobs; i++) {
            for (int v = 0; v < nvars; v++) {
                if (!var_percent[v]) continue;
                double x = results[v][i];
                if (x >= missval) results[v][i] = missval;
                else results[v][i] = (opts.float_out ? (double)(float)x : x) / 100.0;
            }
        }
    }

    t_convert = ctools_timer_seconds();

    /* Without force, the ado leaves variables with nonnumeric values untouched */
    ST_retcode mac_rc = save_failed_counts(var_failed, nvars);
    if (mac_rc) {
        ctools_filtered_data_free(&filtered);
        free(results_block);
        free_options(&opts);
        return mac_rc;
    }

    /* ====================================================================
     * Phase 3: Store results back to Stata via ctools_store_filtered
     * Uses obs_map to write to correct Stata observations.
     * Variables the ado will leave untouched are not stored.
     * ==================================================================== */

    int store_rc = STATA_OK;
    if (nvars == 1) {
        if (opts.force || var_failed[0] == 0)
            store_rc = ctools_store_filtered_rowpar(results[0], nobs, opts.dst_indices[0], obs_map);
    } else {
        #pragma omp parallel for schedule(static) reduction(max:store_rc)
        for (int v = 0; v < nvars; v++) {
            if (!opts.force && var_failed[v] > 0) continue;
            int rc = ctools_store_filtered(results[v], nobs, opts.dst_indices[v], obs_map);
            if (rc > store_rc) store_rc = rc;
        }
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
    if (store_rc) { free_options(&opts); return 459; }

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

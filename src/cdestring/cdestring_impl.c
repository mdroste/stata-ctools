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
 *
 * The ado runs this as two calls. "parse" does steps 1 and 2, keeps the
 * numbers, and reports for each variable its nonnumeric count, missing count,
 * the storage type compress would give it and the ignore() characters it
 * contained. The ado then creates each destination with its final type and
 * "store" does step 3, so Stata never counts, scans or recasts the result.
 * A call without either word converts and stores at once into destinations
 * that already exist.
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
#define CDESTRING_MAX_FLOAT     0x1.fffffep+126  /* maxfloat(); larger is missing in a float */

/* ============================================================================
 * Option Parsing
 * ============================================================================ */

typedef struct {
    int *src_indices;       /* Source string variable indices (1-based) */
    int *dst_indices;       /* Destination numeric variable indices (1-based) */
    int nvars;              /* Number of variables to process */
    char *ignore_chars;      /* Complete decoded ignore() list */
    int ignore_len;         /* Length of ignore_chars string */
    int strip;              /* 1 = remove ignore_table bytes (ignore() and, with percent, "%") */
    int asbytes;            /* 1 = ignore() lists bytes, not UTF-8 characters */
    int illegal;            /* 1 = also remove bytes outside well-formed UTF-8 */
    int force;              /* 1 = keep conversions of variables with nonnumeric values */
    int percent;            /* 1 = strip %; divide the variable by 100 if any value had % */
    int dpcomma;            /* 1 = first comma is the decimal point; any "." is nonnumeric */
    int float_out;          /* 1 = outputs are float (range check, percent rounding) */
    int verbose;            /* 1 = print timing info */
    int selection;          /* CTOOLS_LOAD_SKIP_IF without an if condition,
                               CTOOLS_LOAD_MARKER with a sample marker last */
} cdestring_options;

/*
 * Check if a character should be ignored during parsing.
 * Uses a lookup table for O(1) checking.
 */
static unsigned char ignore_table[256];
/* Where each ignore_table byte first occurs in the ignore() list (-1: not
 * listed), and where U+FFFD first occurs when the list holds characters.
 * convert_value() reports the characters a value contained by these offsets. */
static int ignore_first[256];
static int ignore_fffd;

static void build_ignore_table(const char *ignore_chars, int len, int asbytes)
{
    memset(ignore_table, 0, sizeof(ignore_table));
    for (int c = 0; c < 256; c++) ignore_first[c] = -1;
    ignore_fffd = -1;
    for (int i = 0; i < len; i++) {
        unsigned char c = (unsigned char)ignore_chars[i];
        /* As characters, a byte of 128 or more is part of a multibyte
         * sequence, matched whole in convert_value(). */
        if (asbytes || c < 128) {
            ignore_table[c] = 1;
            if (ignore_first[c] < 0) ignore_first[c] = i;
        }
        if (!asbytes && ignore_fffd < 0 && i + 3 <= len &&
            !memcmp(ignore_chars + i, "\xef\xbf\xbd", 3)) ignore_fffd = i;
    }
}

/* Length of the well-formed UTF-8 sequence at s (Unicode Table 3-7), or 0
 * when the byte at s does not start one. Stata's Unicode functions read each
 * byte of an ill-formed sequence as U+FFFD, which ignore(, illegal) removes. */
static int utf8_sequence_len(const unsigned char *s, int avail)
{
    unsigned c = s[0];
    if (c < 0x80) return 1;
    if (c < 0xc2 || c > 0xf4) return 0;
    int width = c >= 0xf0 ? 4 : c >= 0xe0 ? 3 : 2;
    if (avail < width) return 0;
    unsigned lo = 0x80, hi = 0xbf;
    if (c == 0xe0) lo = 0xa0;
    else if (c == 0xed) hi = 0x9f;
    else if (c == 0xf0) lo = 0x90;
    else if (c == 0xf4) hi = 0x8f;
    if (s[1] < lo || s[1] > hi) return 0;
    for (int k = 2; k < width; k++)
        if ((s[k] & 0xc0) != 0x80) return 0;
    return width;
}

/* Unicode White_Space, as accepted by ustrtrim(), encoded as whole UTF-8
 * sequences. Invalid/truncated sequences are left for numeric rejection. */
static int trim_space_len(const char *text, int len)
{
    const unsigned char *s = (const unsigned char *)text;
    if (!len) return 0;
    if (s[0] == ' ' || (s[0] >= 9 && s[0] <= 13)) return 1;
    if (len >= 2 && s[0] == 0xc2 && (s[1] == 0x85 || s[1] == 0xa0)) return 2;
    if (len < 3) return 0;
    if ((s[0] == 0xe1 && s[1] == 0x9a && s[2] == 0x80) ||
        (s[0] == 0xe2 && s[1] == 0x80 &&
         ((s[2] >= 0x80 && s[2] <= 0x8a) || s[2] == 0xa8 || s[2] == 0xa9 || s[2] == 0xaf)) ||
        (s[0] == 0xe2 && s[1] == 0x81 && s[2] == 0x9f) ||
        (s[0] == 0xe3 && s[1] == 0x80 && s[2] == 0x80)) return 3;
    return 0;
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

static int hex_digit(char c)
{
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

/*
 * real() also reads Stata's hexadecimal format (%21x), [+-]1.hhhX[+-]hhh:
 * one to 13 mantissa digits and one to three exponent digits, either case.
 * The parts are the IEEE bits: the exponent -3ff gives a subnormal (0.hhh
 * times 2^-3fe), and a mantissa written 0.000... is only accepted there, as
 * zero. Exponents past +-3ff are not numbers; 1.0X+3ff and above are Stata's
 * missing encodings, retained under force but counted as nonnumeric.
 * destring first cuts leading
 * zeros to one, so "00.0X-3ff" is zero too.
 */
static int parse_stata_hex(const char *s, int len, double *out)
{
    int i = 0, negative = 0;
    if (len > 0 && (s[0] == '+' || s[0] == '-')) negative = (s[i++] == '-');
    else while (len - i > 1 && s[i] == '0' && s[i + 1] == '0') i++;
    if (i >= len || (s[i] != '0' && s[i] != '1')) return 0;
    int lead = s[i++] - '0';
    if (i >= len || s[i++] != '.') return 0;
    uint64_t mantissa = 0;
    int digits = 0, d;
    while (i < len && (d = hex_digit(s[i])) >= 0) {
        if (++digits > 13) return 0;
        mantissa = (mantissa << 4) | (uint64_t)d;
        i++;
    }
    if (!digits || i >= len || (s[i] != 'X' && s[i] != 'x')) return 0;
    mantissa <<= 4 * (13 - digits);
    i++;
    if (i >= len || (s[i] != '+' && s[i] != '-')) return 0;
    int exp_negative = s[i++] == '-';
    int exponent = 0, exp_digits = 0;
    while (i < len && (d = hex_digit(s[i])) >= 0) {
        if (++exp_digits > 3) return 0;
        exponent = exponent * 16 + d;
        i++;
    }
    if (!exp_digits || i != len) return 0;
    if (exp_negative) exponent = -exponent;
    uint64_t bits;
    if (!lead) {
        if (mantissa || exponent != -0x3ff) return 0;
        bits = 0;
    } else {
        if (exponent < -0x3ff || exponent > 0x3ff) return 0;
        bits = ((uint64_t)(exponent + 0x3ff) << 52) | mantissa;
    }
    if (negative) bits |= (uint64_t)1 << 63;
    memcpy(out, &bits, sizeof(bits));
    return 1;
}

/*
 * Convert one string as native destring does: remove ignore() bytes (and "%"
 * with percent), trim, keep "", "." and ".a"-".z" as missing values, apply
 * dpcomma, then parse like real(). Returns 0 for a nonnumeric value, which is
 * stored as missing; without force the ado then leaves the variable untouched.
 * Returns -1 for allocation failure, which must never be treated as bad text.
 *
 * seen, when not NULL, has one byte per byte of the ignore() list: the byte at
 * the offset where a character first occurs in the list is set when the value
 * contains that character. The byte of U+FFFD is also set for a byte outside
 * well-formed UTF-8, which Stata's Unicode functions read as U+FFFD.
 */
static int convert_value(const char *src, const cdestring_options *opts,
                         double missval, double *out, unsigned char *seen)
{
    char buf[CDESTRING_STR_BUF_SIZE + 1];
    const char *s = src;
    int len = (int)strlen(src);

    if (opts->strip || opts->dpcomma) {
        int n = 0;
        for (int i = 0; i < len && n < CDESTRING_STR_BUF_SIZE; i++) {
            unsigned char c = (unsigned char)src[i];
            if (ignore_table[c]) {
                if (seen && ignore_first[c] >= 0) seen[ignore_first[c]] = 1;
                continue;
            }
            if (opts->asbytes) {
                buf[n++] = src[i];
                continue;
            }
            /* Match complete UTF-8 codepoints, never a union of their bytes. */
            int width = utf8_sequence_len((const unsigned char *)src + i, len - i);
            if (!width) {
                if (seen && ignore_fffd >= 0) seen[ignore_fffd] = 1;
                if (!opts->illegal) buf[n++] = src[i];
                continue;
            }
            int ignored = 0;
            if (width > 1) {
                for (int j = 0; j + width <= opts->ignore_len; j++) {
                    if (!memcmp(src + i, opts->ignore_chars + j, (size_t)width)) {
                        ignored = 1;
                        if (seen) seen[j] = 1;
                        break;
                    }
                }
            }
            if (!ignored) {
                /* Preserve a valid sequence whole. Revisiting its continuation
                 * bytes would classify them as illegal standalone bytes. */
                if (width > CDESTRING_STR_BUF_SIZE - n) {
                    *out = missval;
                    return 0;
                }
                memcpy(buf + n, src + i, (size_t)width);
                n += width;
            }
            i += width - 1;
        }
        s = buf;
        len = n;
    }
    int space;
    while ((space = trim_space_len(s, len))) { s += space; len -= space; }
    while (len > 0) {
        int at = len - 1;
        while (at > 0 && len - at < 3 && ((unsigned char)s[at] & 0xc0) == 0x80) at--;
        if (trim_space_len(s + at, len - at) != len - at) break;
        len = at;
    }

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
    ctools_number_status parsed = ctools_parse_double_status(s, len, &val, missval);
    if (parsed == CTOOLS_NUMBER_MEMORY) return -1;
    if (parsed == CTOOLS_NUMBER_INVALID) {
        char dexp[CDESTRING_STR_BUF_SIZE + 1];
        if (len > CDESTRING_STR_BUF_SIZE) return 0;
        memcpy(dexp, s, (size_t)len);
        parsed = d_exponent_to_e(dexp, len)
                     ? ctools_parse_double_status(dexp, len, &val, missval)
                     : CTOOLS_NUMBER_INVALID;
        if (parsed == CTOOLS_NUMBER_MEMORY) return -1;
        if (parsed == CTOOLS_NUMBER_INVALID) {
            if (!parse_stata_hex(s, len, &val)) return 0;
            if (val >= missval) {
                /* Native destring counts these as nonnumeric, but force
                 * retains real()'s missing encoding before storage/compress. */
                *out = val;
                return 0;
            }
        }
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
    size_t size = (size_t)nvars * 12 + 1;
    char *buf = (char *)malloc(size);
    if (!buf) return 920;
    size_t used = 0;
    buf[0] = '\0';
    for (int v = 0; v < nvars; v++) {
        used += (size_t)snprintf(buf + used, size - used, "%s%d",
                                 v ? " " : "", failed[v]);
    }
    ST_retcode rc = SF_macro_save("_cdestring_failed", buf);
    free(buf);
    return rc;
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
static int parse_options(const char *args, cdestring_options *opts, int staged)
{
    memset(opts, 0, sizeof(cdestring_options));

    if (args == NULL || *args == '\0') {
        return -1;
    }

    /* Parse nvars= */
    if (ctools_parse_int_option(args, "nvars", &opts->nvars) != CTOOLS_PARSE_OK) {
        return -1;
    }
    if (opts->nvars <= 0 || opts->nvars > SF_nvars()) {
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

    if (staged) {
        /* "parse": the plugin varlist holds the sources only. The
         * destinations do not exist until the ado knows their types. */
        for (int i = 0; i < opts->nvars; i++) {
            opts->src_indices[i] = i + 1;
            opts->dst_indices[i] = 0;
        }
    } else if (ctools_parse_bool_option(args, "seqpairs")) {
        /* The plugin varlist holds the sources, then the destinations: no
         * index list, which on wide data exceeds Stata's plugin argument
         * length. */
        if (2 * (size_t)opts->nvars > (size_t)SF_nvars()) {
            free(opts->src_indices);
            free(opts->dst_indices);
            return -1;
        }
        for (int i = 0; i < opts->nvars; i++) {
            opts->src_indices[i] = i + 1;
            opts->dst_indices[i] = opts->nvars + i + 1;
        }
    } else {
        /* Variable index pairs from the start of args */
        const char *cursor = args;
        int *indices = malloc((size_t)opts->nvars * 2 * sizeof(int));
        if (!indices) {
            free(opts->src_indices);
            free(opts->dst_indices);
            return -2;
        }
        if (ctools_parse_int_array(indices, (size_t)(opts->nvars * 2), &cursor) != 0) {
            free(indices);
            free(opts->src_indices);
            free(opts->dst_indices);
            return -1;
        }
        for (int i = 0; i < opts->nvars; i++) {
            opts->src_indices[i] = indices[i * 2];
            opts->dst_indices[i] = indices[i * 2 + 1];
        }
        free(indices);
    }

    /* The ado passes raw bytes through a local macro, never command syntax.
     * Keep ignore= for older callers and native fixtures. */
    int raw_len = 0;
    int raw_rc = ctools_parse_int_option(args, "ignorelen", &raw_len);
    if (raw_rc < 0 || raw_len < 0 || raw_len == INT_MAX) {
        free(opts->src_indices);
        free(opts->dst_indices);
        return -1;
    }
    const char *ignore = strstr(args, "ignore=");
    size_t capacity = raw_rc == CTOOLS_PARSE_OK ? (size_t)raw_len + 1
                       : ignore ? strcspn(ignore + 7, " \t") + 1 : 1;
    if (capacity > INT_MAX || !(opts->ignore_chars = calloc(capacity, 1))) {
        free(opts->src_indices);
        free(opts->dst_indices);
        return -2;
    }
    if (raw_rc == CTOOLS_PARSE_OK) {
        if (SF_macro_use("___cds_ignore", opts->ignore_chars, (ST_int)capacity) ||
            strlen(opts->ignore_chars) != (size_t)raw_len ||
            SF_macro_save("_cdestring_ignore_read", "1")) {
            free(opts->ignore_chars);
            free(opts->src_indices);
            free(opts->dst_indices);
            return -1;
        }
        opts->ignore_len = raw_len;
    } else {
        opts->ignore_len = parse_ignore_option(args, opts->ignore_chars, (int)capacity);
    }

    /* Boolean options */
    opts->force = ctools_parse_bool_option(args, "force");
    opts->percent = ctools_parse_bool_option(args, "percent");
    opts->dpcomma = ctools_parse_bool_option(args, "dpcomma");
    opts->float_out = ctools_parse_bool_option(args, "float");
    opts->verbose = ctools_parse_bool_option(args, "verbose");
    opts->selection = ctools_parse_bool_option(args, "skipif") ? CTOOLS_LOAD_SKIP_IF
                    : ctools_parse_bool_option(args, "touse") ? CTOOLS_LOAD_MARKER : CTOOLS_LOAD_CHECK_IF;
    opts->asbytes = ctools_parse_bool_option(args, "asbytes");
    opts->illegal = ctools_parse_bool_option(args, "illegal");
    opts->strip = (opts->ignore_len > 0 || opts->percent || opts->illegal);

    return 0;
}

static void free_options(cdestring_options *opts)
{
    free(opts->src_indices);
    free(opts->dst_indices);
    free(opts->ignore_chars);
    opts->ignore_chars = NULL;
    opts->src_indices = NULL;
    opts->dst_indices = NULL;
}

/* ============================================================================
 * Public Entry Point
 * ============================================================================ */

/* Each column fits Stata's observation limit; their combined count need not. */
static uint64_t cdestring_total(const int *counts, int nvars)
{
    uint64_t total = 0;
    for (int v = 0; v < nvars; v++) total += (uint64_t)counts[v];
    return total;
}

/* ============================================================================
 * Staged conversion: "parse" keeps the numbers until "store" writes them
 * ============================================================================ */

typedef struct {
    int active;
    double *block;           /* nvars columns of nobs converted values */
    perm_idx_t *obs_map;     /* Stata observation of each loaded row (a copy) */
    unsigned char *publish;  /* 1: the ado creates a destination, store fills it */
    size_t nobs;
    int nvars;
    int n_publish;
    ST_int data_nobs;        /* SF_nobs() at parse; store requires the same */
    double elapsed;          /* seconds parse took */
} cdestring_stage;

static cdestring_stage g_stage;

void cdestring_cleanup_state(void)
{
    free(g_stage.block);
    free(g_stage.obs_map);
    free(g_stage.publish);
    memset(&g_stage, 0, sizeof(g_stage));
}

typedef struct {
    double min, max;         /* over the nonmissing values */
    uint64_t missing;
    int fraction;            /* some nonmissing value is not an integer */
} cdestring_stat;

static void cdestring_stat_init(cdestring_stat *st, int nvars)
{
    for (int v = 0; v < nvars; v++) {
        st[v].min = HUGE_VAL;
        st[v].max = -HUGE_VAL;
        st[v].missing = 0;
        st[v].fraction = 0;
    }
}

/*
 * The values each published variable will hold. Under float they are rounded
 * to float here, as storing into a float variable would round them, so that
 * the type chosen from them is exact. Returns -1 on allocation failure.
 */
static int cdestring_describe(double **results, size_t nobs, int nvars,
                              const unsigned char *publish, int float_out,
                              double missval, cdestring_stat *stats)
{
    int alloc_failed = 0;
    cdestring_stat_init(stats, nvars);
    #pragma omp parallel
    {
        cdestring_stat *tl = (cdestring_stat *)malloc((size_t)nvars * sizeof(cdestring_stat));
        if (tl) {
            cdestring_stat_init(tl, nvars);
        } else {
            #pragma omp atomic write
            alloc_failed = 1;
        }

        #pragma omp for schedule(static)
        for (size_t i = 0; i < nobs; i++) {
            if (!tl) continue;
            for (int v = 0; v < nvars; v++) {
                if (!publish[v]) continue;
                double x = results[v][i];
                if (!(x < missval)) {
                    tl[v].missing++;
                    continue;
                }
                if (float_out) {
                    x = (double)(float)x;
                    results[v][i] = x;
                }
                /* Every double of magnitude 2^52 or more is an integer */
                if (fabs(x) < 4503599627370496.0 && x != (double)(int64_t)x)
                    tl[v].fraction = 1;
                if (x < tl[v].min) tl[v].min = x;
                if (x > tl[v].max) tl[v].max = x;
            }
        }

        if (tl) {
            #pragma omp critical
            {
                for (int v = 0; v < nvars; v++) {
                    stats[v].missing += tl[v].missing;
                    stats[v].fraction |= tl[v].fraction;
                    if (tl[v].min < stats[v].min) stats[v].min = tl[v].min;
                    if (tl[v].max > stats[v].max) stats[v].max = tl[v].max;
                }
            }
        }
        free(tl);
    }
    return alloc_failed ? -1 : 0;
}

/*
 * The storage type compress gives a double (under float: a float) with these
 * values: byte, int or long when every value is an integer the type holds
 * (a float goes no further than int), byte when every value is missing.
 */
static const char *cdestring_storage_type(const cdestring_stat *st, size_t nobs, int float_out)
{
    if (st->missing == (uint64_t)nobs) return "byte";
    if (!st->fraction) {
        if (st->min >= -127.0 && st->max <= 100.0) return "byte";
        if (st->min >= -32767.0 && st->max <= 32740.0) return "int";
        if (!float_out && st->min >= -2147483647.0 && st->max <= 2147483620.0) return "long";
    }
    return float_out ? "float" : "double";
}

/*
 * Save what the ado needs to create and describe each destination: locals
 * cdestring_types and cdestring_missing (one word per variable; "." for a
 * variable left untouched) and cdestring_seen (per variable, the offsets in
 * the ignore() list of the characters its values contained, joined by
 * commas, or "." for none).
 */
static ST_retcode cdestring_report(const cdestring_stat *stats, const unsigned char *publish,
                                   const unsigned char *seen, size_t seen_width,
                                   int nvars, size_t nobs, int float_out)
{
    size_t seen_marks = 0;
    for (size_t k = 0; seen && k < (size_t)nvars * seen_width; k++) seen_marks += seen[k];
    size_t types_cap = (size_t)nvars * 8 + 1;
    size_t missing_cap = (size_t)nvars * 22 + 1;
    size_t seen_cap = seen_marks * 12 + (size_t)nvars * 2 + 1;
    char *types = (char *)malloc(types_cap);
    char *missing = (char *)malloc(missing_cap);
    char *marks = (char *)malloc(seen_cap);
    ST_retcode rc = 920;
    if (types && missing && marks) {
        size_t tu = 0, mu = 0, su = 0;
        types[0] = missing[0] = marks[0] = '\0';
        for (int v = 0; v < nvars; v++) {
            const char *sep = v ? " " : "";
            if (publish[v]) {
                tu += (size_t)snprintf(types + tu, types_cap - tu, "%s%s", sep,
                                       cdestring_storage_type(&stats[v], nobs, float_out));
                mu += (size_t)snprintf(missing + mu, missing_cap - mu, "%s%llu", sep,
                                       (unsigned long long)stats[v].missing);
            } else {
                tu += (size_t)snprintf(types + tu, types_cap - tu, "%s.", sep);
                mu += (size_t)snprintf(missing + mu, missing_cap - mu, "%s.", sep);
            }
            su += (size_t)snprintf(marks + su, seen_cap - su, "%s", sep);
            int any = 0;
            for (size_t k = 0; seen && k < seen_width; k++) {
                if (!seen[(size_t)v * seen_width + k]) continue;
                su += (size_t)snprintf(marks + su, seen_cap - su, "%s%zu", any ? "," : "", k);
                any = 1;
            }
            if (!any) su += (size_t)snprintf(marks + su, seen_cap - su, ".");
        }
        rc = SF_macro_save("_cdestring_types", types);
        if (!rc) rc = SF_macro_save("_cdestring_missing", missing);
        if (!rc) rc = SF_macro_save("_cdestring_seen", marks);
    }
    free(types);
    free(missing);
    free(marks);
    return rc;
}

/*
 * "store": write the staged numbers into the destinations the ado created,
 * one plugin variable per published variable in their order.
 */
static ST_retcode cdestring_store_staged(void)
{
    double t_start = ctools_timer_seconds();
    if (!g_stage.active || SF_nvars() != g_stage.n_publish) {
        SF_error("cdestring: store does not match the parsed variables\n");
        cdestring_cleanup_state();
        return 198;
    }
    if (SF_nobs() != g_stage.data_nobs) {
        SF_error("cdestring: the data changed between parse and store\n");
        cdestring_cleanup_state();
        return 459;
    }

    size_t nobs = g_stage.nobs;
    int nvars = g_stage.nvars;
    int store_rc = STATA_OK;
    int *dst = (int *)malloc((size_t)(nvars ? nvars : 1) * sizeof(int));
    if (!dst) {
        cdestring_cleanup_state();
        SF_error("cdestring: memory allocation failed\n");
        return 920;
    }
    int next = 0;
    for (int v = 0; v < nvars; v++) dst[v] = g_stage.publish[v] ? ++next : 0;

    if (!nobs) {
        /* nothing was selected: the destinations stay missing */
    } else if (g_stage.n_publish == 1) {
        for (int v = 0; v < nvars; v++) {
            if (dst[v])
                store_rc = ctools_store_filtered_rowpar(g_stage.block + (size_t)v * nobs,
                                                        nobs, dst[v], g_stage.obs_map);
        }
    } else {
        #pragma omp parallel for schedule(static) reduction(max:store_rc)
        for (int v = 0; v < nvars; v++) {
            if (!dst[v]) continue;
            int rc = ctools_store_filtered(g_stage.block + (size_t)v * nobs, nobs,
                                           dst[v], g_stage.obs_map);
            if (rc > store_rc) store_rc = rc;
        }
    }
    free(dst);

    double t_store = ctools_timer_seconds() - t_start;
    double t_total = g_stage.elapsed + t_store;
    cdestring_cleanup_state();
    if (store_rc) return ctools_stata_rc(store_rc);

    ST_retcode rc;
    if ((rc = SF_scal_save("_cdestring_time_store", t_store))) return rc;
    if ((rc = SF_scal_save("_cdestring_time_total", t_total))) return rc;
    if ((rc = SF_scal_save("_cdestring_threads_max", (double)ctools_get_max_threads()))) return rc;
    return SF_scal_save("_cdestring_openmp_enabled", (double)CTOOLS_OPENMP_ENABLED);
}

/*
 * Finish "parse": describe the converted values, report them to the ado and
 * keep them for "store". block passes to the stage, which also frees it when
 * anything here fails; the observation map is copied, so the caller still
 * frees loaded.
 */
static ST_retcode cdestring_stage_finish(double *block, double **results,
                                         const ctools_filtered_data *loaded, size_t nobs,
                                         const cdestring_options *opts, const int *var_failed,
                                         const unsigned char *seen, size_t seen_width,
                                         double missval, uint64_t total_converted,
                                         uint64_t total_failed, double t_start,
                                         double t_parse, double t_load)
{
    int nvars = opts->nvars;
    cdestring_cleanup_state();
    g_stage.block = block;
    if (nobs) {
        g_stage.obs_map = (perm_idx_t *)malloc(nobs * sizeof(perm_idx_t));
        if (g_stage.obs_map) memcpy(g_stage.obs_map, loaded->obs_map, nobs * sizeof(perm_idx_t));
    }
    g_stage.publish = (unsigned char *)malloc((size_t)nvars);
    cdestring_stat *stats = (cdestring_stat *)malloc((size_t)nvars * sizeof(cdestring_stat));
    ST_retcode rc = 920;
    if (g_stage.publish && stats && (!nobs || g_stage.obs_map)) {
        g_stage.n_publish = 0;
        for (int v = 0; v < nvars; v++) {
            g_stage.publish[v] = (unsigned char)(opts->force || var_failed[v] == 0);
            g_stage.n_publish += g_stage.publish[v];
        }
        if (cdestring_describe(results, nobs, nvars, g_stage.publish, opts->float_out,
                               missval, stats) == 0)
            rc = cdestring_report(stats, g_stage.publish, seen, seen_width, nvars, nobs,
                                  opts->float_out);
    }
    free(stats);
    if (rc == 920) SF_error("cdestring: memory allocation failed\n");

    double t_done = ctools_timer_seconds();
    if (!rc) rc = SF_scal_save("_cdestring_n_converted", (double)total_converted);
    if (!rc) rc = SF_scal_save("_cdestring_n_failed", (double)total_failed);
    if (!rc) rc = SF_scal_save("_cdestring_nobs", (double)nobs);
    if (!rc) rc = SF_scal_save("_cdestring_time_parse", t_parse - t_start);
    if (!rc) rc = SF_scal_save("_cdestring_time_load", t_load - t_parse);
    if (!rc) rc = SF_scal_save("_cdestring_time_convert", t_done - t_load);
    /* store replaces these two; they stand when no variable is stored */
    if (!rc) rc = SF_scal_save("_cdestring_time_store", 0.0);
    if (!rc) rc = SF_scal_save("_cdestring_time_total", t_done - t_start);
    if (!rc) rc = SF_scal_save("_cdestring_threads_max", (double)ctools_get_max_threads());
    if (!rc) rc = SF_scal_save("_cdestring_openmp_enabled", (double)CTOOLS_OPENMP_ENABLED);
    if (rc) {
        cdestring_cleanup_state();
        return rc;
    }
    g_stage.nobs = nobs;
    g_stage.nvars = nvars;
    g_stage.data_nobs = SF_nobs();
    g_stage.elapsed = t_done - t_start;
    g_stage.active = 1;
    return 0;
}

static ST_retcode cdestring_convert(const char *args, int staged);

ST_retcode cdestring_main(const char *args)
{
    char op[16] = {0};
    if (args) sscanf(args, "%15s", op);
    if (!strcmp(op, "store")) return cdestring_store_staged();
    cdestring_cleanup_state();
    if (!strcmp(op, "discard")) return 0;
    return cdestring_convert(args, !strcmp(op, "parse"));
}

static ST_retcode cdestring_convert(const char *args, int staged)
{
    double t_start, t_parse, t_load, t_convert, t_store, t_total;
    cdestring_options opts;

    t_start = ctools_timer_seconds();

    /* Parse arguments */
    int parse_rc = parse_options(args, &opts, staged);
    if (parse_rc != 0) {
        if (parse_rc == -2) {
            SF_error("cdestring: memory allocation failed\n");
            return 920;
        }
        SF_error("cdestring: invalid arguments\n");
        return 198;
    }

    /* Build ignore character lookup table; percent implies ignore("%") */
    build_ignore_table(opts.ignore_chars, opts.ignore_len, opts.asbytes);
    if (opts.percent) ignore_table['%'] = 1;

    t_parse = ctools_timer_seconds();

    int nvars = opts.nvars;
    const double missval = SV_missval;

    /* Per-variable counters: converted, failed and percent, nvars each */
    int *var_counts = (int *)calloc((size_t)nvars * 3, sizeof(int));
    if (!var_counts) {
        free_options(&opts);
        SF_error("cdestring: memory allocation failed\n");
        return 920;
    }
    int *var_converted = var_counts;
    int *var_failed = var_counts + nvars;
    int *var_percent = var_counts + 2 * (size_t)nvars;

    /* ====================================================================
     * Phase 1: Bulk load source string variables with if/in filtering
     * Uses ctools_data_load() to load only filtered observations.
     * ==================================================================== */

    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    stata_retcode load_rc = ctools_data_load(&filtered, opts.src_indices, (size_t)nvars, 0, 0,
                                              opts.selection | CTOOLS_LOAD_NO_SORT_ORDER);
    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        free_options(&opts);
        free(var_counts);
        SF_error("cdestring: failed to load string data\n");
        return ctools_stata_rc(load_rc);
    }

    size_t nobs = filtered.data.nobs;
    perm_idx_t *obs_map = filtered.obs_map;

    if (nobs == 0 && staged) {
        /* Nothing was selected: every destination is all missing. */
        ST_retcode empty_rc = save_failed_counts(var_failed, nvars);
        if (!empty_rc)
            empty_rc = cdestring_stage_finish(NULL, NULL, &filtered, 0, &opts, var_failed, NULL, 0,
                                              missval, 0, 0, t_start, t_parse,
                                              ctools_timer_seconds());
        ctools_filtered_data_free(&filtered);
        free_options(&opts);
        free(var_counts);
        return empty_rc;
    }
    if (nobs == 0) {
        ctools_filtered_data_free(&filtered);
        free_options(&opts);
        ST_retcode empty_rc = save_failed_counts(var_failed, nvars);
        free(var_counts);
        const char *names[] = {"_cdestring_n_converted", "_cdestring_n_failed",
            "_cdestring_time_parse", "_cdestring_time_load", "_cdestring_time_convert",
            "_cdestring_time_store", "_cdestring_time_total"};
        for (size_t j = 0; j < sizeof(names) / sizeof(names[0]) && !empty_rc; j++)
            empty_rc = SF_scal_save((char *)names[j], 0);
        return empty_rc;
    }

    t_load = ctools_timer_seconds();

    /* ====================================================================
     * Phase 2: Parse strings to doubles (parallel across observations)
     * Note: No if_mask needed since data is pre-filtered.
     * ==================================================================== */

    /* Allocate result arrays: single contiguous block for all variables */
    double *results_block = ctools_safe_malloc3((size_t)nvars, nobs, sizeof(double));
    double **results = (double **)malloc((size_t)nvars * sizeof(double *));
    if (!results_block || !results) {
        free(results_block);
        free(results);
        ctools_filtered_data_free(&filtered);
        free_options(&opts);
        free(var_counts);
        SF_error("cdestring: memory allocation failed\n");
        return 920;
    }

    for (int v = 0; v < nvars; v++) {
        results[v] = results_block + (size_t)v * nobs;
    }

    /* "parse" also notes which ignore() characters each variable contained */
    size_t seen_width = staged ? (size_t)opts.ignore_len : 0;
    unsigned char *var_seen = NULL;
    if (seen_width && !(var_seen = (unsigned char *)calloc((size_t)nvars, seen_width))) {
        free(results_block);
        free(results);
        ctools_filtered_data_free(&filtered);
        free_options(&opts);
        free(var_counts);
        SF_error("cdestring: memory allocation failed\n");
        return 920;
    }

    /* Transposed loop: single OMP region with outer loop over observations,
     * inner loop over variables. This gives one fork/join barrier instead of
     * nvars barriers, and keeps per-observation data cache-hot across variables. */
    int tl_alloc_failed = 0;
    #pragma omp parallel
    {
        /* Thread-local per-variable counters */
        int *tl_counts = (int *)calloc((size_t)nvars * 3, sizeof(int));
        unsigned char *tl_seen = seen_width ? (unsigned char *)calloc((size_t)nvars, seen_width) : NULL;
        if (seen_width && !tl_seen) {
            free(tl_counts);
            tl_counts = NULL;
        }
        int tl_ok = tl_counts != NULL;
        int *tl_converted = tl_counts;
        int *tl_failed = tl_ok ? tl_counts + nvars : NULL;
        int *tl_percent = tl_ok ? tl_counts + 2 * (size_t)nvars : NULL;
        if (!tl_ok) {
            #pragma omp atomic write
            tl_alloc_failed = 1;
        }

        #pragma omp for schedule(static)
        for (size_t i = 0; i < nobs; i++) {
            if (!tl_ok) continue;
            for (int v = 0; v < nvars; v++) {
                const char *src = filtered.data.vars[v].data.str[i];
                if (src == NULL) src = "";

                if (opts.percent && strchr(src, '%')) tl_percent[v]++;
                int converted = convert_value(src, &opts, missval, &results[v][i],
                                              tl_seen ? tl_seen + (size_t)v * seen_width : NULL);
                if (converted < 0) {
                    #pragma omp atomic write
                    tl_alloc_failed = 1;
                } else if (converted) {
                    tl_converted[v]++;
                } else {
                    tl_failed[v]++;
                }
            }
        }

        /* Aggregate thread-local counters */
        if (tl_ok) {
            #pragma omp critical
            {
                for (int v = 0; v < nvars; v++) {
                    var_converted[v] += tl_converted[v];
                    var_failed[v] += tl_failed[v];
                    var_percent[v] += tl_percent[v];
                }
                for (size_t k = 0; tl_seen && k < (size_t)nvars * seen_width; k++)
                    var_seen[k] |= tl_seen[k];
            }
        }
        free(tl_counts);
        free(tl_seen);
    }
    if (tl_alloc_failed) {
        ctools_filtered_data_free(&filtered);
        free(results_block);
        free(results);
        free(var_seen);
        free(var_counts);
        free_options(&opts);
        SF_error("cdestring: memory allocation failed\n");
        return 920;
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
        free(results);
        free(var_seen);
        free(var_counts);
        free_options(&opts);
        return mac_rc;
    }

    if (staged) {
        /* The strings are no longer needed; the numbers stay for "store". */
        ST_retcode stage_rc = cdestring_stage_finish(
            results_block, results, &filtered, nobs, &opts, var_failed, var_seen, seen_width,
            missval, cdestring_total(var_converted, nvars), cdestring_total(var_failed, nvars),
            t_start, t_parse, t_load);
        ctools_filtered_data_free(&filtered);
        free(results);
        free(var_seen);
        free(var_counts);
        free_options(&opts);
        return stage_rc;
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
    uint64_t total_converted = cdestring_total(var_converted, nvars);
    uint64_t total_failed = cdestring_total(var_failed, nvars);

    free(results_block);
    free(results);
    free(var_counts);
    if (store_rc) { free_options(&opts); return ctools_stata_rc(store_rc); }

    t_total = t_store - t_start;

    /* Store results */
    ST_retcode result_rc;
    if ((result_rc = SF_scal_save("_cdestring_n_converted", (double)total_converted))) { free_options(&opts); return result_rc; }
    if ((result_rc = SF_scal_save("_cdestring_n_failed", (double)total_failed))) { free_options(&opts); return result_rc; }
    if ((result_rc = SF_scal_save("_cdestring_time_parse", t_parse - t_start))) { free_options(&opts); return result_rc; }
    if ((result_rc = SF_scal_save("_cdestring_time_load", t_load - t_parse))) { free_options(&opts); return result_rc; }
    if ((result_rc = SF_scal_save("_cdestring_time_convert", t_convert - t_load))) { free_options(&opts); return result_rc; }
    if ((result_rc = SF_scal_save("_cdestring_time_store", t_store - t_convert))) { free_options(&opts); return result_rc; }
    if ((result_rc = SF_scal_save("_cdestring_time_total", t_total))) { free_options(&opts); return result_rc; }

    if ((result_rc = SF_scal_save("_cdestring_threads_max", (double)ctools_get_max_threads()))) { free_options(&opts); return result_rc; }
    if ((result_rc = SF_scal_save("_cdestring_openmp_enabled", (double)CTOOLS_OPENMP_ENABLED))) { free_options(&opts); return result_rc; }

    free_options(&opts);

    return 0;
}

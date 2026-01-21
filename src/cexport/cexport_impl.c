/*
    cexport_impl.c
    High-Performance Delimited Text Export for Stata

    Exports Stata datasets to CSV/delimited text files with parallel processing.
    Designed for maximum throughput on large datasets.

    Architecture:
    1. PHASE 1: Load data from Stata into C memory (parallel per-variable)
    2. PHASE 2: Format data into text buffers (parallel per-chunk)
    3. PHASE 3: Write buffers to file (sequential, buffered I/O)

    Performance Optimizations:
    - Parallel data loading via ctools_data_load
    - Chunked parallel formatting: each thread formats a range of rows
    - Custom fast double-to-string conversion (avoids sprintf overhead)
    - Large I/O buffers (64KB) to minimize syscalls
    - Memory-mapped output option for very large files
    - Pre-computed field widths for aligned output (optional)

    Thread Safety:
    - Each formatting thread works on disjoint row ranges
    - No shared mutable state during formatting phase
    - Single writer thread for I/O (sequential)
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <errno.h>

#include "ctools_threads.h"

#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_timer.h"
#include "cexport_impl.h"
#include "cexport_io.h"

/* ========================================================================
   Configuration Constants
   Use shared constants from ctools_config.h where applicable.
   ======================================================================== */

/* Shared constants (from ctools_config.h):
   - CTOOLS_IO_BUFFER_SIZE: I/O buffer size (64KB)
   - CTOOLS_EXPORT_CHUNK_SIZE: rows per parallel chunk (10K)
   - CTOOLS_IO_MAX_THREADS: max threads for I/O (16)
*/

/* Local constants specific to cexport */
#define CEXPORT_DBL_BUFFER_SIZE  32             /* Buffer for double formatting */
#define CEXPORT_INITIAL_ROW_BUF  4096           /* Initial row buffer size */

/* ========================================================================
   Variable Storage Types (from Stata)
   ======================================================================== */

typedef enum {
    VARTYPE_STRING = 0,
    VARTYPE_BYTE   = 1,
    VARTYPE_INT    = 2,
    VARTYPE_LONG   = 3,
    VARTYPE_FLOAT  = 4,
    VARTYPE_DOUBLE = 5
} vartype_t;

/* ========================================================================
   Export Context and Configuration
   ======================================================================== */

typedef struct {
    /* Output file */
    char *filename;
    FILE *fp;                    /* Legacy: only used for fallback single-threaded path */

    /* Formatting options */
    char delimiter;
    bool write_header;
    bool quote_strings;
    bool quote_if_needed;
    bool verbose;

    /* Line ending configuration (for binary mode with explicit endings) */
    char line_ending[3];         /* "\n" or "\r\n" + null terminator */
    int line_ending_len;         /* 1 for LF, 2 for CRLF */
    bool use_crlf;               /* True for Windows-style CRLF */

    /* I/O backend configuration */
    cexport_io_backend_t io_backend;  /* PWRITE or MMAP */
    cexport_io_flags_t io_flags;      /* I/O flags (DIRECT, NOFSYNC, PREFAULT) */
    bool use_parallel_io;        /* True to use parallel offset writes */

    /* Adaptive sizing (computed at runtime) */
    size_t actual_avg_row_size;  /* Sampled average row size */
    size_t adaptive_chunk_size;  /* Rows per chunk (adaptive) */

    /* Data from Stata */
    stata_data data;

    /* Observation range for if/in filtering */
    size_t obs1;    /* First observation (1-based Stata index) */
    size_t nobs_loaded;  /* Number of observations loaded */

    /* Variable names (for header) */
    char **varnames;
    size_t nvars;

    /* Variable storage types (for numeric precision) */
    vartype_t *vartypes;

    /* Optimization flags */
    bool all_numeric;    /* True if all variables are numeric (fast path) */

    /* Timing */
    double time_load;
    double time_format;
    double time_write;
    double time_total;
} cexport_context;

/* Global context for this export operation */
static cexport_context g_ctx;

/* Timing uses ctools_timer.h */

/* ========================================================================
   Fast Double to String Conversion

   Custom implementation that avoids sprintf overhead.
   Handles integers, decimals, scientific notation, and Stata missing values.

   OPTIMIZATION: Uses custom fast integer conversion to avoid snprintf,
   which provides 5-10x speedup for integer values (common in data).
   ======================================================================== */

/* Check if value is a Stata missing value */
static inline bool is_stata_missing(double val)
{
    /* Stata missing values: . .a .b ... .z are encoded as specific large values */
    return SF_is_missing(val);
}

/*
    NOTE: Integer-to-string conversion functions are now provided by ctools_types.h:
    - ctools_uctools_int64_to_str()
    - ctools_ctools_int64_to_str()
    - ctools_ctools_uint64_to_str_fast()
*/

/*
    Two-digit lookup table for fast digit pair output in format_decimal_fast.
    "00", "01", "02", ... "99"
    Kept local for specialized fractional digit formatting.
*/
static const char DIGIT_PAIRS[200] = {
    '0','0', '0','1', '0','2', '0','3', '0','4', '0','5', '0','6', '0','7', '0','8', '0','9',
    '1','0', '1','1', '1','2', '1','3', '1','4', '1','5', '1','6', '1','7', '1','8', '1','9',
    '2','0', '2','1', '2','2', '2','3', '2','4', '2','5', '2','6', '2','7', '2','8', '2','9',
    '3','0', '3','1', '3','2', '3','3', '3','4', '3','5', '3','6', '3','7', '3','8', '3','9',
    '4','0', '4','1', '4','2', '4','3', '4','4', '4','5', '4','6', '4','7', '4','8', '4','9',
    '5','0', '5','1', '5','2', '5','3', '5','4', '5','5', '5','6', '5','7', '5','8', '5','9',
    '6','0', '6','1', '6','2', '6','3', '6','4', '6','5', '6','6', '6','7', '6','8', '6','9',
    '7','0', '7','1', '7','2', '7','3', '7','4', '7','5', '7','6', '7','7', '7','8', '7','9',
    '8','0', '8','1', '8','2', '8','3', '8','4', '8','5', '8','6', '8','7', '8','8', '8','9',
    '9','0', '9','1', '9','2', '9','3', '9','4', '9','5', '9','6', '9','7', '9','8', '9','9'
};

/*
    Fast decimal formatting with full double precision (15 significant digits).
    Avoids snprintf entirely for common float values.

    Returns: number of characters written, or -1 if fallback needed

    NOTE: Currently unused - Stata requires 16 digits for doubles, this uses 15.
    Kept for potential future use with float variables.
*/
#if defined(__GNUC__) || defined(__clang__)
__attribute__((unused))
#endif
static int format_decimal_fast(double val, char *buf)
{
    /* Handle negative numbers */
    int pos = 0;
    if (val < 0) {
        buf[pos++] = '-';
        val = -val;
    }

    /* Split into integer and fractional parts */
    double int_part_d = floor(val);
    double frac_part = val - int_part_d;

    /* Integer part must fit in uint64 */
    if (int_part_d > 9007199254740992.0) {
        return -1;  /* Fallback to snprintf */
    }

    uint64_t int_part = (uint64_t)int_part_d;

    /* Format integer part */
    pos += ctools_uint64_to_str_fast(int_part, buf + pos);

    /* Handle fractional part - use 15 decimal places for full precision */
    if (frac_part > 5e-16) {
        buf[pos++] = '.';

        /* Multiply by 10^15 and round
         * 10^15 = 1,000,000,000,000,000 fits in uint64 */
        uint64_t frac_scaled = (uint64_t)(frac_part * 1000000000000000.0 + 0.5);

        /* Format 15 digits using digit pairs, then strip trailing zeros */
        char frac_buf[16];
        int frac_len = 15;

        /* Format pairs from least significant */
        int idx;
        idx = (frac_scaled % 100) * 2; frac_scaled /= 100;
        frac_buf[14] = DIGIT_PAIRS[idx + 1]; frac_buf[13] = DIGIT_PAIRS[idx];
        idx = (frac_scaled % 100) * 2; frac_scaled /= 100;
        frac_buf[12] = DIGIT_PAIRS[idx + 1]; frac_buf[11] = DIGIT_PAIRS[idx];
        idx = (frac_scaled % 100) * 2; frac_scaled /= 100;
        frac_buf[10] = DIGIT_PAIRS[idx + 1]; frac_buf[9] = DIGIT_PAIRS[idx];
        idx = (frac_scaled % 100) * 2; frac_scaled /= 100;
        frac_buf[8] = DIGIT_PAIRS[idx + 1]; frac_buf[7] = DIGIT_PAIRS[idx];
        idx = (frac_scaled % 100) * 2; frac_scaled /= 100;
        frac_buf[6] = DIGIT_PAIRS[idx + 1]; frac_buf[5] = DIGIT_PAIRS[idx];
        idx = (frac_scaled % 100) * 2; frac_scaled /= 100;
        frac_buf[4] = DIGIT_PAIRS[idx + 1]; frac_buf[3] = DIGIT_PAIRS[idx];
        idx = (frac_scaled % 100) * 2; frac_scaled /= 100;
        frac_buf[2] = DIGIT_PAIRS[idx + 1]; frac_buf[1] = DIGIT_PAIRS[idx];
        /* Last digit (15th) */
        frac_buf[0] = '0' + (frac_scaled % 10);

        /* Strip trailing zeros */
        while (frac_len > 1 && frac_buf[frac_len - 1] == '0') {
            frac_len--;
        }

        /* Copy fractional digits */
        for (int i = 0; i < frac_len; i++) {
            buf[pos++] = frac_buf[i];
        }
    }

    buf[pos] = '\0';
    return pos;
}

/*
    Powers of 10 for fast decimal formatting.
    POW10[i] = 10^i for i in [0, 16]
*/
static const uint64_t POW10[17] = {
    1ULL,
    10ULL,
    100ULL,
    1000ULL,
    10000ULL,
    100000ULL,
    1000000ULL,
    10000000ULL,
    100000000ULL,
    1000000000ULL,
    10000000000ULL,
    100000000000ULL,
    1000000000000ULL,
    10000000000000ULL,
    100000000000000ULL,
    1000000000000000ULL,
    10000000000000000ULL
};

/*
    Format unsigned integer with leading zeros to exactly 'width' digits.
    Used for fractional part formatting.
    Returns number of characters written.
*/
static inline int format_frac_digits(uint64_t val, char *buf, int width)
{
    /* Format right-to-left */
    char temp[20];
    int pos = width;
    temp[pos] = '\0';

    while (pos > 0) {
        pos--;
        temp[pos] = '0' + (val % 10);
        val /= 10;
    }

    /* Copy and strip trailing zeros */
    int len = width;
    while (len > 1 && temp[len - 1] == '0') {
        len--;
    }

    memcpy(buf, temp, len);
    return len;
}

/*
    Fast double to string - optimized version.

    Key optimizations:
    1. Fast path for integers (handles ~80% of typical data)
    2. Fast path for "simple decimals" with few decimal places
    3. Direct formatting without leading zero (Stata style)
    4. Single-pass trailing zero removal
    5. Fallback to snprintf only for complex cases (scientific notation)

    Returns: number of characters written (not including null terminator)
*/
static int double_to_str(double val, char *buf, int buf_size, bool missing_as_dot, vartype_t vtype)
{
    (void)buf_size;  /* Used for safety in fallback path */

    /* Handle Stata missing values */
    if (is_stata_missing(val)) {
        if (missing_as_dot) {
            buf[0] = '.';
            buf[1] = '\0';
            return 1;
        } else {
            buf[0] = '\0';
            return 0;
        }
    }

    /* For float storage type, truncate to single precision before formatting.
     * This ensures we match Stata's native export delimited output exactly.
     * Without this, the float-to-double conversion introduces spurious precision
     * (e.g., 0.1f becomes 0.10000000149011612 as double). */
    if (vtype == VARTYPE_FLOAT) {
        val = (double)(float)val;
    }

    /* Handle special floating point values */
    if (isnan(val)) {
        buf[0] = '.';
        buf[1] = '\0';
        return 1;
    }
    if (isinf(val)) {
        if (val > 0) {
            memcpy(buf, "inf", 4);
            return 3;
        } else {
            memcpy(buf, "-inf", 5);
            return 4;
        }
    }

    /* FAST PATH 1: Exact integers
     * Handles ~80% of typical Stata data (IDs, counts, years, etc.)
     */
    if (val == (double)(int64_t)val && val >= -9007199254740992.0 && val <= 9007199254740992.0) {
        return ctools_int64_to_str((int64_t)val, buf);
    }

    /* Determine precision based on storage type */
    int max_frac_digits = (vtype == VARTYPE_FLOAT) ? 8 : 16;

    /* Handle sign */
    int pos = 0;
    double abs_val = val;
    if (val < 0) {
        buf[pos++] = '-';
        abs_val = -val;
    }

    /* FAST PATH 2: Simple decimals (most common in real data)
     * Values that can be represented exactly with limited decimal places.
     * Handles values like 1.5, 0.25, 123.456, etc.
     *
     * Upper bound: 1e9 to avoid precision issues with large floats.
     * Values >= 1e9 may not convert back exactly when reconstructed.
     */
    if (abs_val < 1e9 && abs_val > 1e-10) {
        /* Extract integer and fractional parts */
        double int_part_d = floor(abs_val);
        double frac_part = abs_val - int_part_d;

        /* Check if fractional part is "simple" (few decimal places needed) */
        if (frac_part > 0) {
            /* Try progressively more decimal places until we get exact match */
            for (int decimals = 1; decimals <= max_frac_digits && decimals <= 15; decimals++) {
                double scale = (double)POW10[decimals];
                double scaled_frac = frac_part * scale;
                uint64_t frac_int = (uint64_t)(scaled_frac + 0.5);

                /* Check for rounding to next integer */
                if (frac_int >= POW10[decimals]) {
                    int_part_d += 1.0;
                    frac_int = 0;
                    frac_part = 0;
                    break;
                }

                /* Check if this representation is exact enough */
                double reconstructed = int_part_d + (double)frac_int / scale;
                double rel_error = fabs(reconstructed - abs_val) / abs_val;

                if (rel_error < 1e-15 || (vtype == VARTYPE_FLOAT && rel_error < 1e-7)) {
                    /* Found exact representation! Format directly. */
                    uint64_t int_part = (uint64_t)int_part_d;

                    if (int_part == 0) {
                        /* Value between -1 and 1: Stata omits leading zero */
                        buf[pos++] = '.';
                    } else {
                        /* Format integer part */
                        pos += ctools_uint64_to_str_fast(int_part, buf + pos);
                        buf[pos++] = '.';
                    }

                    /* Format fractional part with trailing zeros stripped */
                    if (frac_int == 0) {
                        /* No fractional part after rounding - remove decimal point */
                        pos--;
                    } else {
                        pos += format_frac_digits(frac_int, buf + pos, decimals);
                    }

                    buf[pos] = '\0';
                    return pos;
                }
            }
        } else {
            /* No fractional part - should have been caught by integer path,
             * but handle it anyway */
            uint64_t int_part = (uint64_t)int_part_d;
            pos += ctools_uint64_to_str_fast(int_part, buf + pos);
            buf[pos] = '\0';
            return pos;
        }
    }

    /* FALLBACK: Use snprintf for complex cases
     * - Very small numbers (scientific notation)
     * - Very large numbers
     * - Numbers requiring many decimal places
     *
     * Stata uses scientific notation for very small numbers (< 0.0001)
     * with format like "1.00000000000e-10" (11 decimal places in mantissa)
     */
    char temp[48];
    int len;

    /* Check if this will use scientific notation (very small numbers) */
    if (abs_val > 0 && abs_val < 0.0001) {
        /* Use Stata-compatible scientific notation format
         * Stata uses 11 decimal places for doubles, 7 for floats */
        if (vtype == VARTYPE_FLOAT) {
            len = snprintf(temp, sizeof(temp), "%.7e", val);
        } else {
            len = snprintf(temp, sizeof(temp), "%.11e", val);
        }
        /* Copy directly - no post-processing needed for scientific notation */
        memcpy(buf + pos, temp + (temp[0] == '-' ? 1 : 0), len - (temp[0] == '-' ? 1 : 0));
        pos += len - (temp[0] == '-' ? 1 : 0);
        buf[pos] = '\0';
        return pos;
    }

    /* For other cases, use %g format */
    if (vtype == VARTYPE_FLOAT) {
        len = snprintf(temp, sizeof(temp), "%.8g", val);
    } else {
        len = snprintf(temp, sizeof(temp), "%.16g", val);
    }

    /* Process snprintf output in single pass:
     * - Skip leading zero for values between -1 and 1
     * - Strip trailing zeros after decimal point
     */
    int src = 0;
    int has_sign = (temp[0] == '-');
    if (has_sign) {
        src++;
    }

    /* Skip leading zero: "0." -> "." */
    int skip_zero = (temp[src] == '0' && temp[src + 1] == '.');
    if (skip_zero) {
        src++;
    }

    /* Find decimal point and exponent positions */
    int dot_pos = -1;
    int exp_pos = -1;
    for (int i = src; i < len; i++) {
        if (temp[i] == '.') dot_pos = i;
        else if (temp[i] == 'e' || temp[i] == 'E') {
            exp_pos = i;
            break;
        }
    }

    /* Copy to output, stripping trailing zeros */
    pos = 0;
    if (has_sign) {
        buf[pos++] = '-';
    }

    if (dot_pos >= 0 && exp_pos < 0) {
        /* Has decimal, no exponent - can strip trailing zeros */
        int last_nonzero = len - 1;
        while (last_nonzero > dot_pos && temp[last_nonzero] == '0') {
            last_nonzero--;
        }
        /* If only decimal point left, remove it too */
        if (last_nonzero == dot_pos) {
            last_nonzero--;
        }
        for (int i = src; i <= last_nonzero; i++) {
            buf[pos++] = temp[i];
        }
    } else if (dot_pos >= 0 && exp_pos > 0) {
        /* Has decimal and exponent - copy as-is (Stata keeps trailing zeros in sci notation) */
        for (int i = src; i < len; i++) {
            buf[pos++] = temp[i];
        }
    } else {
        /* No decimal point - copy as-is */
        for (int i = src; i < len; i++) {
            buf[pos++] = temp[i];
        }
    }

    buf[pos] = '\0';
    return pos;
}

/* ========================================================================
   String Quoting and Escaping
   ======================================================================== */

/*
    Check if a string needs quoting.

    A string needs quoting if it contains:
    - The delimiter character
    - A double quote
    - A newline or carriage return
*/
static inline bool string_needs_quoting(const char *str, char delimiter)
{
    if (str == NULL) return false;

    for (const char *p = str; *p != '\0'; p++) {
        if (*p == delimiter || *p == '"' || *p == '\n' || *p == '\r') {
            return true;
        }
    }
    return false;
}

/*
    Write a quoted string to a buffer.

    Doubles any internal quotes (CSV escaping convention).
    Returns number of characters written.
*/
static int write_quoted_string(const char *str, char *buf, int buf_size)
{
    /* Handle edge cases: need at least 3 bytes for empty quoted string '""' + null */
    if (buf_size < 1) {
        return 0;
    }
    if (buf_size < 3) {
        buf[0] = '\0';
        return 0;
    }
    if (str == NULL) {
        buf[0] = '"';
        buf[1] = '"';
        buf[2] = '\0';
        return 2;
    }

    int pos = 0;
    buf[pos++] = '"';

    for (const char *p = str; *p != '\0' && pos < buf_size - 2; p++) {
        if (*p == '"') {
            /* Escape quote by doubling it */
            if (pos < buf_size - 3) {
                buf[pos++] = '"';
                buf[pos++] = '"';
            }
        } else {
            buf[pos++] = *p;
        }
    }

    buf[pos++] = '"';
    buf[pos] = '\0';
    return pos;
}

/* ========================================================================
   Row Formatting
   ======================================================================== */

/*
    Thread arguments for parallel row formatting.
*/
typedef struct {
    size_t start_row;        /* First row to format (0-based) */
    size_t end_row;          /* Last row to format (exclusive) */
    char *output_buffer;     /* Output buffer for formatted rows */
    size_t buffer_size;      /* Size of output buffer */
    size_t bytes_written;    /* Actual bytes written to buffer */
    int success;             /* 1 on success, 0 on failure */
} format_chunk_args_t;

/*
    Thread arguments for parallel offset writes.
*/
typedef struct {
    cexport_io_file *file;   /* Output file handle */
    const char *buffer;      /* Data to write */
    size_t len;              /* Number of bytes to write */
    size_t offset;           /* File offset to write at */
    int success;             /* 1 on success, 0 on failure */
} write_chunk_args_t;

/*
    Thread function: Write a chunk at a specific offset.
    Used for parallel I/O with non-overlapping ranges.
*/
static void *write_chunk_thread(void *arg)
{
    write_chunk_args_t *args = (write_chunk_args_t *)arg;

    if (args->len == 0) {
        args->success = 1;
        return NULL;
    }

    ssize_t written = cexport_io_write_at(args->file, args->buffer,
                                           args->len, args->offset);
    args->success = (written == (ssize_t)args->len);

    return args->success ? NULL : (void *)(intptr_t)-1;
}

/*
    Format a single row of ALL NUMERIC variables into the output buffer.

    FAST PATH: No type checking, no string handling.
    Used when g_ctx.all_numeric is true for maximum throughput.

    Returns: number of bytes written, or -1 on buffer overflow
*/
static int format_row_numeric(size_t row_idx, char *buf, size_t buf_size)
{
    size_t pos = 0;
    size_t nvars = g_ctx.data.nvars;
    char delimiter = g_ctx.delimiter;
    const char *line_end = g_ctx.line_ending;
    int line_end_len = g_ctx.line_ending_len;

    for (size_t j = 0; j < nvars; j++) {
        double val = g_ctx.data.vars[j].data.dbl[row_idx];
        vartype_t vtype = g_ctx.vartypes[j];

        /* Check buffer space (assume max 32 chars for number + delimiter/newline) */
        if (pos + 34 + line_end_len > buf_size) {
            return -1;
        }

        /* Write directly to output buffer */
        int field_len = double_to_str(val, buf + pos, 32, false, vtype);
        pos += field_len;

        /* Add delimiter or line ending */
        if (j < nvars - 1) {
            buf[pos++] = delimiter;
        } else {
            /* Copy line ending (1 or 2 bytes) */
            memcpy(buf + pos, line_end, line_end_len);
            pos += line_end_len;
        }
    }

    return (int)pos;
}

/*
    Format a single row into the output buffer.

    OPTIMIZED:
    - Small stack buffer for field formatting (256 bytes handles 99.9% of cases)
    - Inline numeric path (most common) to avoid function call overhead
    - Write directly to output buffer when possible

    Returns: number of bytes written, or -1 on buffer overflow
*/
static int format_row(size_t row_idx, char *buf, size_t buf_size)
{
    char field_buf[256];  /* Small buffer - most fields are <50 chars */
    size_t pos = 0;
    size_t nvars = g_ctx.data.nvars;
    char delimiter = g_ctx.delimiter;
    const char *line_end = g_ctx.line_ending;
    int line_end_len = g_ctx.line_ending_len;

    for (size_t j = 0; j < nvars; j++) {
        stata_variable *var = &g_ctx.data.vars[j];
        int field_len;

        if (var->type == STATA_TYPE_DOUBLE) {
            /* FAST PATH: Numeric variable - inline for speed */
            double val = var->data.dbl[row_idx];
            vartype_t vtype = g_ctx.vartypes[j];

            /* Check buffer space (assume max 32 chars for number + line ending) */
            if (pos + 34 + line_end_len > buf_size) {
                return -1;
            }

            /* Write directly to output buffer */
            field_len = double_to_str(val, buf + pos, sizeof(field_buf), false, vtype);
            pos += field_len;
        } else {
            /* String variable */
            const char *str = var->data.str[row_idx];
            if (str == NULL) str = "";

            bool need_quote = g_ctx.quote_strings ||
                             (g_ctx.quote_if_needed && string_needs_quoting(str, delimiter));

            if (need_quote) {
                field_len = write_quoted_string(str, field_buf, sizeof(field_buf));
                if (pos + field_len + 2 + line_end_len > buf_size) return -1;
                memcpy(buf + pos, field_buf, field_len);
            } else {
                field_len = (int)strlen(str);
                if (pos + field_len + 2 + line_end_len > buf_size) return -1;
                memcpy(buf + pos, str, field_len);
            }
            pos += field_len;
        }

        /* Add delimiter or line ending */
        if (j < nvars - 1) {
            buf[pos++] = delimiter;
        } else {
            /* Copy line ending (1 or 2 bytes) */
            memcpy(buf + pos, line_end, line_end_len);
            pos += line_end_len;
        }
    }

    return (int)pos;
}

/*
    Thread function: Format a chunk of rows.
    Uses fast path for all-numeric datasets.
*/
static void *format_chunk_thread(void *arg)
{
    format_chunk_args_t *args = (format_chunk_args_t *)arg;
    size_t pos = 0;
    bool all_numeric = g_ctx.all_numeric;

    for (size_t i = args->start_row; i < args->end_row; i++) {
        /* Check if this observation satisfies the if condition */
        ST_int stata_obs = (ST_int)(g_ctx.obs1 + i);
        if (!SF_ifobs(stata_obs)) {
            continue;  /* Skip observations that don't match if condition */
        }

        int row_len;
        if (all_numeric) {
            row_len = format_row_numeric(i, args->output_buffer + pos, args->buffer_size - pos);
        } else {
            row_len = format_row(i, args->output_buffer + pos, args->buffer_size - pos);
        }
        if (row_len < 0) {
            args->success = 0;
            args->bytes_written = pos;
            return NULL;
        }
        pos += row_len;
    }

    args->success = 1;
    args->bytes_written = pos;
    return NULL;
}

/* ========================================================================
   Header Writing
   ======================================================================== */

/*
    Format the header row into a buffer.
    Returns the number of bytes written, or -1 on error.
    The caller must provide a buffer of sufficient size.
*/
static int format_header(char *buf, size_t buf_size)
{
    size_t pos = 0;
    const char *line_end = g_ctx.line_ending;
    int line_end_len = g_ctx.line_ending_len;

    for (size_t j = 0; j < g_ctx.nvars; j++) {
        const char *name = g_ctx.varnames[j];
        size_t name_len = strlen(name);

        /* Check if name needs quoting */
        bool need_quote = g_ctx.quote_strings ||
                         (g_ctx.quote_if_needed && string_needs_quoting(name, g_ctx.delimiter));

        if (need_quote) {
            if (pos + name_len * 2 + 4 + line_end_len > buf_size) {
                return -1;
            }
            int quoted_len = write_quoted_string(name, buf + pos, buf_size - pos);
            pos += quoted_len;
        } else {
            if (pos + name_len + 2 + line_end_len > buf_size) {
                return -1;
            }
            memcpy(buf + pos, name, name_len);
            pos += name_len;
        }

        if (j < g_ctx.nvars - 1) {
            buf[pos++] = g_ctx.delimiter;
        } else {
            /* Copy line ending (1 or 2 bytes) */
            memcpy(buf + pos, line_end, line_end_len);
            pos += line_end_len;
        }
    }

    return (int)pos;
}


/* ========================================================================
   Main Export Function
   ======================================================================== */

/*
    Parse command arguments.

    Format: filename delimiter [options...]
    Options: noheader, quote, quoteif, verbose, crlf, mmap
*/
static int parse_args(const char *args)
{
    char *args_copy = strdup(args);
    if (args_copy == NULL) return -1;

    char *saveptr;
    char *token;
    int arg_idx = 0;

    /* Initialize defaults */
    g_ctx.filename = NULL;
    g_ctx.delimiter = ',';
    g_ctx.write_header = true;
    g_ctx.quote_strings = false;
    g_ctx.quote_if_needed = true;  /* Default: quote only if needed */
    g_ctx.verbose = false;

    /* Line ending defaults: LF (Unix-style) for consistency with Stata */
    g_ctx.use_crlf = false;
    g_ctx.line_ending[0] = '\n';
    g_ctx.line_ending[1] = '\0';
    g_ctx.line_ending[2] = '\0';
    g_ctx.line_ending_len = 1;

    /* I/O backend defaults: pwrite with parallel I/O enabled */
    g_ctx.io_backend = CEXPORT_IO_PWRITE;
    g_ctx.io_flags = CEXPORT_IO_FLAG_NONE;
    g_ctx.use_parallel_io = true;

    /* Adaptive sizing defaults (computed later) */
    g_ctx.actual_avg_row_size = 0;
    g_ctx.adaptive_chunk_size = CTOOLS_EXPORT_CHUNK_SIZE;

    token = strtok_r(args_copy, " \t", &saveptr);
    while (token != NULL) {
        if (arg_idx == 0) {
            /* First arg: filename */
            g_ctx.filename = strdup(token);
        } else if (arg_idx == 1) {
            /* Second arg: delimiter */
            if (strlen(token) == 1) {
                g_ctx.delimiter = token[0];
            } else if (strcmp(token, "tab") == 0) {
                g_ctx.delimiter = '\t';
            }
        } else {
            /* Options */
            if (strcmp(token, "noheader") == 0) {
                g_ctx.write_header = false;
            } else if (strcmp(token, "quote") == 0) {
                g_ctx.quote_strings = true;
            } else if (strcmp(token, "noquoteif") == 0) {
                g_ctx.quote_if_needed = false;
            } else if (strcmp(token, "verbose") == 0) {
                g_ctx.verbose = true;
            } else if (strcmp(token, "crlf") == 0) {
                /* Windows-style line endings */
                g_ctx.use_crlf = true;
                g_ctx.line_ending[0] = '\r';
                g_ctx.line_ending[1] = '\n';
                g_ctx.line_ending[2] = '\0';
                g_ctx.line_ending_len = 2;
            } else if (strcmp(token, "mmap") == 0) {
                /* Use memory-mapped I/O backend */
                g_ctx.io_backend = CEXPORT_IO_MMAP;
            } else if (strcmp(token, "noparallel") == 0) {
                /* Disable parallel I/O (for debugging/comparison) */
                g_ctx.use_parallel_io = false;
            } else if (strcmp(token, "nofsync") == 0) {
                /* Skip final fsync for faster but less durable writes */
                g_ctx.io_flags |= CEXPORT_IO_FLAG_NOFSYNC;
            } else if (strcmp(token, "direct") == 0) {
                /* Direct I/O bypasses OS cache (for very large files) */
                g_ctx.io_flags |= CEXPORT_IO_FLAG_DIRECT;
            } else if (strcmp(token, "prefault") == 0) {
                /* Pre-fault mmap pages to avoid page fault latency */
                g_ctx.io_flags |= CEXPORT_IO_FLAG_PREFAULT;
            }
        }

        arg_idx++;
        token = strtok_r(NULL, " \t", &saveptr);
    }

    free(args_copy);

    if (g_ctx.filename == NULL) {
        SF_error("cexport: no output filename specified\n");
        return -1;
    }

    return 0;
}

/*
    Load variable names from Stata macro.
    The .ado file sets _cexport_varnames macro with space-separated variable names.
*/
static int load_varnames(void)
{
    size_t nvars = SF_nvars();
    g_ctx.varnames = (char **)ctools_safe_calloc2(nvars, sizeof(char *));
    if (g_ctx.varnames == NULL) return -1;

    g_ctx.nvars = nvars;

    /* Try to get variable names from macro */
    char varnames_buf[32768];
    size_t names_parsed = 0;
    if (SF_macro_use("CEXPORT_VARNAMES", varnames_buf, sizeof(varnames_buf)) == 0 &&
        strlen(varnames_buf) > 0) {
        /* Parse space-separated names */
        char *p = varnames_buf;
        for (size_t j = 0; j < nvars; j++) {
            /* Skip whitespace */
            while (*p == ' ' || *p == '\t') p++;
            if (*p == '\0') break;

            /* Find end of name */
            char *start = p;
            while (*p != ' ' && *p != '\t' && *p != '\0') p++;

            size_t len = p - start;
            g_ctx.varnames[j] = (char *)malloc(len + 1);
            if (g_ctx.varnames[j] == NULL) return -1;
            memcpy(g_ctx.varnames[j], start, len);
            g_ctx.varnames[j][len] = '\0';
            names_parsed++;
        }

        /* Fill remaining entries with fallback names if macro was short */
        for (size_t j = names_parsed; j < nvars; j++) {
            char namebuf[32];
            snprintf(namebuf, sizeof(namebuf), "v%zu", j + 1);
            g_ctx.varnames[j] = strdup(namebuf);
            if (g_ctx.varnames[j] == NULL) return -1;
        }
    } else {
        /* Fallback: use generic names */
        for (size_t j = 0; j < nvars; j++) {
            char namebuf[32];
            snprintf(namebuf, sizeof(namebuf), "v%zu", j + 1);
            g_ctx.varnames[j] = strdup(namebuf);
            if (g_ctx.varnames[j] == NULL) return -1;
        }
    }

    return 0;
}

/*
    Load variable types from Stata macro.
    The .ado file sets CEXPORT_VARTYPES macro with space-separated type codes:
    0=string, 1=byte, 2=int, 3=long, 4=float, 5=double
*/
static int load_vartypes(void)
{
    size_t nvars = g_ctx.nvars;
    g_ctx.vartypes = (vartype_t *)ctools_safe_malloc2(nvars, sizeof(vartype_t));
    if (g_ctx.vartypes == NULL) return -1;

    /* Initialize to double (safest default - most precision) */
    for (size_t j = 0; j < nvars; j++) {
        g_ctx.vartypes[j] = VARTYPE_DOUBLE;
    }

    /* Try to get variable types from macro */
    char vartypes_buf[32768];
    if (SF_macro_use("CEXPORT_VARTYPES", vartypes_buf, sizeof(vartypes_buf)) == 0 &&
        strlen(vartypes_buf) > 0) {
        /* Parse space-separated type codes */
        char *p = vartypes_buf;
        for (size_t j = 0; j < nvars; j++) {
            /* Skip whitespace */
            while (*p == ' ' || *p == '\t') p++;
            if (*p == '\0') break;

            /* Parse type code */
            int type_code = (int)strtol(p, &p, 10);
            if (type_code >= 0 && type_code <= 5) {
                g_ctx.vartypes[j] = (vartype_t)type_code;
            }
        }
    }

    return 0;
}

/*
    Free all resources in context.
*/
static void cleanup_context(void)
{
    if (g_ctx.filename != NULL) {
        free(g_ctx.filename);
        g_ctx.filename = NULL;
    }

    if (g_ctx.fp != NULL) {
        fclose(g_ctx.fp);
        g_ctx.fp = NULL;
    }

    if (g_ctx.varnames != NULL) {
        for (size_t j = 0; j < g_ctx.nvars; j++) {
            if (g_ctx.varnames[j] != NULL) {
                free(g_ctx.varnames[j]);
            }
        }
        free(g_ctx.varnames);
        g_ctx.varnames = NULL;
    }

    if (g_ctx.vartypes != NULL) {
        free(g_ctx.vartypes);
        g_ctx.vartypes = NULL;
    }

    stata_data_free(&g_ctx.data);
}

/* ========================================================================
   Dynamic Buffer Sizing

   Sample first N rows to compute actual average row size,
   avoiding the conservative 2x memory allocation.
   ======================================================================== */

#define CEXPORT_SAMPLE_ROWS 100  /* Number of rows to sample for sizing */

/*
    Sample rows to compute actual average row size.
    Returns average bytes per row, or 0 on error.
*/
static size_t sample_row_sizes(size_t nobs)
{
    char sample_buf[4096];
    size_t sample_count = (nobs < CEXPORT_SAMPLE_ROWS) ? nobs : CEXPORT_SAMPLE_ROWS;
    size_t total_size = 0;
    size_t rows_sampled = 0;
    bool all_numeric = g_ctx.all_numeric;

    /* Sample evenly distributed rows */
    size_t step = (nobs > sample_count) ? (nobs / sample_count) : 1;

    for (size_t i = 0; i < nobs && rows_sampled < sample_count; i += step) {
        /* Check if this observation satisfies the if condition */
        ST_int stata_obs = (ST_int)(g_ctx.obs1 + i);
        if (!SF_ifobs(stata_obs)) {
            continue;
        }

        int row_len;
        if (all_numeric) {
            row_len = format_row_numeric(i, sample_buf, sizeof(sample_buf));
        } else {
            row_len = format_row(i, sample_buf, sizeof(sample_buf));
        }

        if (row_len > 0) {
            total_size += row_len;
            rows_sampled++;
        }
    }

    if (rows_sampled == 0) {
        return 0;
    }

    /* Return average with 20% safety margin (much better than 100% margin) */
    return (total_size / rows_sampled) * 12 / 10 + 64;
}

/*
    Compute adaptive chunk size based on data characteristics.
    All-numeric datasets can use larger chunks (less memory per row).
*/
static size_t compute_adaptive_chunk_size(size_t nobs, bool all_numeric)
{
    size_t base_chunk_size = CTOOLS_EXPORT_CHUNK_SIZE;

    if (all_numeric) {
        /* All-numeric: use 5x larger chunks (50K rows instead of 10K) */
        base_chunk_size = base_chunk_size * 5;
    }

    /* Cap at 10% of total rows or 100K rows, whichever is smaller */
    size_t max_chunk = nobs / 10;
    if (max_chunk < CTOOLS_EXPORT_CHUNK_SIZE) {
        max_chunk = CTOOLS_EXPORT_CHUNK_SIZE;
    }
    if (max_chunk > 100000) {
        max_chunk = 100000;
    }

    return (base_chunk_size < max_chunk) ? base_chunk_size : max_chunk;
}

/* ========================================================================
   Zero-Copy mmap Formatting

   Thread arguments and functions for formatting directly into mapped memory.
   ======================================================================== */

/*
    Thread arguments for zero-copy mmap formatting.
    Each thread formats directly into mapped_data + offset.
*/
typedef struct {
    char *dest;              /* Destination pointer in mapped memory */
    size_t dest_size;        /* Available space at destination */
    size_t start_row;        /* First row to format (0-based) */
    size_t end_row;          /* Last row to format (exclusive) */
    size_t bytes_written;    /* Actual bytes written */
    int success;             /* 1 on success, 0 on failure */
} format_mmap_args_t;

/*
    Thread function: Format rows directly into mapped memory (zero-copy).
    Eliminates the memcpy step in the mmap write path.
*/
static void *format_mmap_thread(void *arg)
{
    format_mmap_args_t *args = (format_mmap_args_t *)arg;
    size_t pos = 0;
    bool all_numeric = g_ctx.all_numeric;

    for (size_t i = args->start_row; i < args->end_row; i++) {
        /* Check if this observation satisfies the if condition */
        ST_int stata_obs = (ST_int)(g_ctx.obs1 + i);
        if (!SF_ifobs(stata_obs)) {
            continue;
        }

        int row_len;
        if (all_numeric) {
            row_len = format_row_numeric(i, args->dest + pos, args->dest_size - pos);
        } else {
            row_len = format_row(i, args->dest + pos, args->dest_size - pos);
        }

        if (row_len < 0) {
            args->success = 0;
            args->bytes_written = pos;
            return NULL;
        }
        pos += row_len;
    }

    args->success = 1;
    args->bytes_written = pos;
    return NULL;
}

/* ========================================================================
   Single-Threaded Buffered Writer

   Buffers multiple rows before calling fwrite to reduce syscall overhead.
   ======================================================================== */

#define CEXPORT_ROW_BUFFER_ROWS 1000  /* Buffer this many rows before fwrite */

/*
    Write rows with buffering for the single-threaded path.
    Buffers multiple rows and writes in batches to reduce syscalls.
*/
static int write_rows_buffered(FILE *fp, size_t nobs, size_t avg_row_size, size_t *rows_written_out)
{
    bool all_numeric = g_ctx.all_numeric;
    size_t rows_written = 0;

    /* Allocate buffer for multiple rows */
    size_t buffer_size = avg_row_size * CEXPORT_ROW_BUFFER_ROWS;
    if (buffer_size < 65536) buffer_size = 65536;
    if (buffer_size > 4 * 1024 * 1024) buffer_size = 4 * 1024 * 1024;

    char *buffer = (char *)malloc(buffer_size);
    if (buffer == NULL) {
        return -1;
    }

    size_t buf_pos = 0;

    for (size_t i = 0; i < nobs; i++) {
        /* Check if this observation satisfies the if condition */
        ST_int stata_obs = (ST_int)(g_ctx.obs1 + i);
        if (!SF_ifobs(stata_obs)) {
            continue;
        }

        /* Check if we need to flush buffer */
        if (buf_pos > buffer_size - avg_row_size * 2) {
            if (fwrite(buffer, 1, buf_pos, fp) != buf_pos) {
                free(buffer);
                return -1;
            }
            buf_pos = 0;
        }

        /* Format row into buffer */
        int row_len;
        if (all_numeric) {
            row_len = format_row_numeric(i, buffer + buf_pos, buffer_size - buf_pos);
        } else {
            row_len = format_row(i, buffer + buf_pos, buffer_size - buf_pos);
        }

        if (row_len < 0) {
            /* Buffer full - flush and retry */
            if (buf_pos > 0) {
                if (fwrite(buffer, 1, buf_pos, fp) != buf_pos) {
                    free(buffer);
                    return -1;
                }
                buf_pos = 0;
            }

            /* Retry formatting */
            if (all_numeric) {
                row_len = format_row_numeric(i, buffer + buf_pos, buffer_size - buf_pos);
            } else {
                row_len = format_row(i, buffer + buf_pos, buffer_size - buf_pos);
            }

            if (row_len < 0) {
                /* Row too large for buffer - shouldn't happen */
                free(buffer);
                return -1;
            }
        }

        buf_pos += row_len;
        rows_written++;
    }

    /* Flush remaining data */
    if (buf_pos > 0) {
        if (fwrite(buffer, 1, buf_pos, fp) != buf_pos) {
            free(buffer);
            return -1;
        }
    }

    free(buffer);
    *rows_written_out = rows_written;
    return 0;
}

/*
    Main entry point for cexport command.
*/
ST_retcode cexport_main(const char *args)
{
    double t_start, t_phase;
    char msg[512];
    ST_retcode rc;
    size_t nobs, nvars;

    /* Initialize context */
    memset(&g_ctx, 0, sizeof(g_ctx));
    stata_data_init(&g_ctx.data);

    t_start = ctools_timer_seconds();

    /* Parse arguments */
    if (parse_args(args) != 0) {
        cleanup_context();
        return 198;
    }

    /* Get dataset dimensions and store observation range for if/in filtering */
    g_ctx.obs1 = (size_t)SF_in1();
    nobs = SF_in2() - SF_in1() + 1;
    g_ctx.nobs_loaded = nobs;
    nvars = SF_nvars();

    if (nobs == 0) {
        SF_error("cexport: no observations to export\n");
        cleanup_context();
        return 2000;
    }

    if (nvars == 0) {
        SF_error("cexport: no variables to export\n");
        cleanup_context();
        return 2000;
    }

    (void)msg;  /* Used only for error messages now */

    /* ================================================================
       PHASE 1: Load data from Stata
       ================================================================ */
    t_phase = ctools_timer_seconds();

    /* Load variable names */
    if (load_varnames() != 0) {
        SF_error("cexport: failed to load variable names\n");
        cleanup_context();
        return 920;
    }

    /* Load variable types (for numeric precision) */
    if (load_vartypes() != 0) {
        SF_error("cexport: failed to load variable types\n");
        cleanup_context();
        return 920;
    }

    /* Load data using shared ctools infrastructure */
    rc = ctools_data_load(&g_ctx.data, nvars);
    if (rc != STATA_OK) {
        snprintf(msg, sizeof(msg), "cexport: failed to load data (error %d)\n", rc);
        SF_error(msg);
        cleanup_context();
        return 920;
    }

    g_ctx.time_load = ctools_timer_seconds() - t_phase;

    /* Detect if all variables are numeric for fast path */
    g_ctx.all_numeric = true;
    for (size_t j = 0; j < nvars; j++) {
        if (g_ctx.data.vars[j].type != STATA_TYPE_DOUBLE) {
            g_ctx.all_numeric = false;
            break;
        }
    }


    /* ================================================================
       PHASE 2: Format data into chunks (parallel)
       ================================================================ */
    t_phase = ctools_timer_seconds();

    /* Compute adaptive chunk size based on data characteristics */
    g_ctx.adaptive_chunk_size = compute_adaptive_chunk_size(nobs, g_ctx.all_numeric);
    size_t chunk_size = g_ctx.adaptive_chunk_size;

    /* Determine number of chunks and threads */
    size_t num_chunks = (nobs + chunk_size - 1) / chunk_size;
    size_t num_threads = num_chunks;
    if (num_threads > CTOOLS_IO_MAX_THREADS) num_threads = CTOOLS_IO_MAX_THREADS;
    if (num_threads < 1) num_threads = 1;

    /* Dynamic buffer sizing: sample actual row sizes instead of conservative estimate */
    size_t avg_row_size = sample_row_sizes(nobs);
    if (avg_row_size == 0) {
        /* Fallback to conservative estimate if sampling fails */
        if (nvars > SIZE_MAX / 22) {
            SF_error("cexport: variable count too large\n");
            cleanup_context();
            return 920;
        }
        avg_row_size = nvars * 22;
    }
    g_ctx.actual_avg_row_size = avg_row_size;

    /* Compute chunk buffer size with tighter margin (sampling provides real data) */
    size_t chunk_buffer_size;
    if (avg_row_size > SIZE_MAX / chunk_size) {
        SF_error("cexport: buffer size overflow\n");
        cleanup_context();
        return 920;
    }
    /* Only 1.2x margin needed since we sampled actual sizes */
    chunk_buffer_size = chunk_size * avg_row_size * 12 / 10 + 4096;

    /* Format header into a buffer */
    char *header_buf = NULL;
    size_t header_len = 0;
    if (g_ctx.write_header) {
        header_buf = (char *)malloc(65536);
        if (header_buf == NULL) {
            SF_error("cexport: memory allocation failed for header\n");
            cleanup_context();
            return 920;
        }
        int hlen = format_header(header_buf, 65536);
        if (hlen < 0) {
            SF_error("cexport: failed to format header\n");
            free(header_buf);
            cleanup_context();
            return 920;
        }
        header_len = (size_t)hlen;
    }

    /* For small datasets or single-threaded mode, use legacy fwrite path */
    if (!g_ctx.use_parallel_io || nobs < CTOOLS_EXPORT_CHUNK_SIZE || num_threads == 1) {
        /* Single-threaded: format and write directly using legacy FILE* */
        g_ctx.fp = fopen(g_ctx.filename, "wb");  /* Binary mode for explicit line endings */
        if (g_ctx.fp == NULL) {
            snprintf(msg, sizeof(msg), "cexport: cannot open file '%s' for writing: %s\n",
                    g_ctx.filename, strerror(errno));
            SF_error(msg);
            free(header_buf);
            cleanup_context();
            return 603;
        }
        setvbuf(g_ctx.fp, NULL, _IOFBF, CTOOLS_IO_BUFFER_SIZE);

        /* Write header */
        if (header_buf != NULL && header_len > 0) {
            if (fwrite(header_buf, 1, header_len, g_ctx.fp) != header_len) {
                SF_error("cexport: failed to write header\n");
                free(header_buf);
                fclose(g_ctx.fp);
                g_ctx.fp = NULL;
                cleanup_context();
                return 693;
            }
        }
        free(header_buf);
        header_buf = NULL;

        /* Single-threaded with buffering: format multiple rows before fwrite */
        size_t rows_written = 0;
        if (write_rows_buffered(g_ctx.fp, nobs, avg_row_size, &rows_written) != 0) {
            SF_error("cexport: write error during buffered write\n");
            fclose(g_ctx.fp);
            g_ctx.fp = NULL;
            cleanup_context();
            return 693;
        }
        nobs = rows_written;  /* Update count for reporting */

        g_ctx.time_format = ctools_timer_seconds() - t_phase;
        g_ctx.time_write = g_ctx.time_format;

        fclose(g_ctx.fp);
        g_ctx.fp = NULL;
    } else {
        /*
         * Multi-threaded parallel offset write pipeline:
         * 1. Format all chunks in parallel (existing approach)
         * 2. Compute prefix-sum of bytes_written to get file offsets
         * 3. Pre-size the file using cexport_io
         * 4. Write header at offset 0
         * 5. Parallel writes of each chunk at computed offsets
         */
        ctools_persistent_pool *pool = ctools_get_global_pool();
        if (pool == NULL) {
            /* Fall back to single-threaded if pool unavailable */
            if (g_ctx.verbose) {
                SF_display("  Thread pool unavailable, using single-threaded fallback\n");
            }

            g_ctx.fp = fopen(g_ctx.filename, "wb");
            if (g_ctx.fp == NULL) {
                snprintf(msg, sizeof(msg), "cexport: cannot open file '%s': %s\n",
                        g_ctx.filename, strerror(errno));
                SF_error(msg);
                free(header_buf);
                cleanup_context();
                return 603;
            }
            setvbuf(g_ctx.fp, NULL, _IOFBF, CTOOLS_IO_BUFFER_SIZE);

            if (header_buf != NULL && header_len > 0) {
                fwrite(header_buf, 1, header_len, g_ctx.fp);
            }
            free(header_buf);
            header_buf = NULL;

            /* Use buffered writer for fallback path too */
            size_t rows_written = 0;
            if (write_rows_buffered(g_ctx.fp, nobs, avg_row_size, &rows_written) != 0) {
                fclose(g_ctx.fp);
                g_ctx.fp = NULL;
                cleanup_context();
                return 920;
            }

            g_ctx.time_format = ctools_timer_seconds() - t_phase;
            g_ctx.time_write = g_ctx.time_format;

            fclose(g_ctx.fp);
            g_ctx.fp = NULL;
        } else if (g_ctx.io_backend == CEXPORT_IO_MMAP) {
            /* ============================================================
               ZERO-COPY MMAP PATH: Format directly into mapped memory
               Eliminates the buffer allocation and memcpy steps.
               ============================================================ */

            /* Estimate total file size (conservative to ensure we have space) */
            size_t estimated_total = header_len + (nobs * avg_row_size * 15 / 10);

            /* Open and pre-size file with mmap */
            cexport_io_file outfile;
            cexport_io_init(&outfile);

            if (cexport_io_open(&outfile, g_ctx.filename, CEXPORT_IO_MMAP, g_ctx.io_flags) != 0) {
                snprintf(msg, sizeof(msg), "cexport: cannot open file '%s': %s\n",
                        g_ctx.filename, outfile.error_message);
                SF_error(msg);
                free(header_buf);
                cleanup_context();
                return 603;
            }

            if (cexport_io_presize(&outfile, estimated_total) != 0) {
                snprintf(msg, sizeof(msg), "cexport: failed to pre-size file: %s\n",
                        outfile.error_message);
                SF_error(msg);
                cexport_io_close(&outfile, 0);
                free(header_buf);
                cleanup_context();
                return 693;
            }

            /* Get pointer to mapped memory */
            char *mapped_base = cexport_io_get_mapped_ptr(&outfile, 0);
            if (mapped_base == NULL) {
                SF_error("cexport: failed to get mapped memory pointer\n");
                cexport_io_close(&outfile, 0);
                free(header_buf);
                cleanup_context();
                return 693;
            }

            /* Write header directly to mapped memory */
            if (header_buf != NULL && header_len > 0) {
                memcpy(mapped_base, header_buf, header_len);
            }
            free(header_buf);
            header_buf = NULL;

            /* Compute estimated chunk offsets for parallel formatting */
            size_t *chunk_offsets = (size_t *)ctools_safe_malloc2(num_chunks + 1, sizeof(size_t));
            if (chunk_offsets == NULL) {
                SF_error("cexport: memory allocation failed\n");
                cexport_io_close(&outfile, 0);
                cleanup_context();
                return 920;
            }

            /* Estimate offsets based on sampled average */
            chunk_offsets[0] = header_len;
            for (size_t c = 0; c < num_chunks; c++) {
                size_t rows_in_chunk = chunk_size;
                if ((c + 1) * chunk_size > nobs) {
                    rows_in_chunk = nobs - c * chunk_size;
                }
                /* Use 1.5x estimate to ensure enough space per chunk */
                chunk_offsets[c + 1] = chunk_offsets[c] + (rows_in_chunk * avg_row_size * 15 / 10);
            }

            /* Allocate mmap format args */
            format_mmap_args_t *mmap_args = (format_mmap_args_t *)ctools_safe_malloc2(num_chunks, sizeof(format_mmap_args_t));
            if (mmap_args == NULL) {
                SF_error("cexport: memory allocation failed\n");
                free(chunk_offsets);
                cexport_io_close(&outfile, 0);
                cleanup_context();
                return 920;
            }

            /* Setup mmap format arguments */
            for (size_t c = 0; c < num_chunks; c++) {
                mmap_args[c].dest = mapped_base + chunk_offsets[c];
                mmap_args[c].dest_size = chunk_offsets[c + 1] - chunk_offsets[c];
                mmap_args[c].start_row = c * chunk_size;
                mmap_args[c].end_row = (c + 1) * chunk_size;
                if (mmap_args[c].end_row > nobs) mmap_args[c].end_row = nobs;
                mmap_args[c].bytes_written = 0;
                mmap_args[c].success = 0;
            }

            /* Submit all formatting tasks to write directly into mapped memory */
            for (size_t c = 0; c < num_chunks; c++) {
                ctools_persistent_pool_submit(pool, format_mmap_thread, &mmap_args[c]);
            }

            /* Wait for all formatting to complete */
            ctools_persistent_pool_wait(pool);

            g_ctx.time_format = ctools_timer_seconds() - t_phase;

            /* Check results and compute actual total size */
            size_t actual_total = header_len;
            int format_failed = 0;
            for (size_t c = 0; c < num_chunks; c++) {
                if (!mmap_args[c].success) {
                    format_failed = 1;
                    break;
                }
                actual_total = chunk_offsets[c] + mmap_args[c].bytes_written;
            }
            /* The last chunk's end position is the actual total */
            if (!format_failed && num_chunks > 0) {
                actual_total = chunk_offsets[num_chunks - 1] + mmap_args[num_chunks - 1].bytes_written;
            }

            if (format_failed) {
                SF_error("cexport: formatting to mmap failed\n");
                free(mmap_args);
                free(chunk_offsets);
                cexport_io_close(&outfile, 0);
                cleanup_context();
                return 920;
            }

            /* Since we formatted with gaps between chunks (overestimated offsets),
               we need to compact the data by shifting chunks together */
            size_t write_pos = header_len;
            for (size_t c = 0; c < num_chunks; c++) {
                if (mmap_args[c].bytes_written > 0 && chunk_offsets[c] != write_pos) {
                    /* Move this chunk's data to the correct position */
                    memmove(mapped_base + write_pos,
                            mapped_base + chunk_offsets[c],
                            mmap_args[c].bytes_written);
                }
                write_pos += mmap_args[c].bytes_written;
            }
            actual_total = write_pos;

            g_ctx.time_write = ctools_timer_seconds() - t_phase - g_ctx.time_format;

            /* Close and truncate to actual size */
            if (cexport_io_close(&outfile, actual_total) != 0) {
                SF_error("cexport: failed to close file\n");
                free(mmap_args);
                free(chunk_offsets);
                cleanup_context();
                return 693;
            }

            free(mmap_args);
            free(chunk_offsets);

        } else {
            /* ============================================================
               STEP 1: Allocate chunk buffers and format in parallel
               (PWRITE backend - format to buffers, then parallel write)
               ============================================================ */
            format_chunk_args_t *chunk_args = (format_chunk_args_t *)ctools_safe_malloc2(num_chunks, sizeof(format_chunk_args_t));
            char **chunk_buffers = (char **)ctools_safe_malloc2(num_chunks, sizeof(char *));

            if (chunk_args == NULL || chunk_buffers == NULL) {
                free(chunk_args);
                free(chunk_buffers);
                free(header_buf);
                SF_error("cexport: memory allocation failed\n");
                cleanup_context();
                return 920;
            }

            /* Initialize chunk_buffers to NULL for safe cleanup */
            for (size_t c = 0; c < num_chunks; c++) {
                chunk_buffers[c] = NULL;
            }

            /* Allocate buffers for all chunks */
            for (size_t c = 0; c < num_chunks; c++) {
                chunk_buffers[c] = (char *)malloc(chunk_buffer_size);
                if (chunk_buffers[c] == NULL) {
                    for (size_t u = 0; u < num_chunks; u++) {
                        if (chunk_buffers[u]) free(chunk_buffers[u]);
                    }
                    free(chunk_args);
                    free(chunk_buffers);
                    free(header_buf);
                    SF_error("cexport: memory allocation failed for chunk buffers\n");
                    cleanup_context();
                    return 920;
                }
            }

            /* Setup all chunk arguments using adaptive chunk size */
            for (size_t c = 0; c < num_chunks; c++) {
                chunk_args[c].start_row = c * chunk_size;
                chunk_args[c].end_row = (c + 1) * chunk_size;
                if (chunk_args[c].end_row > nobs) chunk_args[c].end_row = nobs;
                chunk_args[c].output_buffer = chunk_buffers[c];
                chunk_args[c].buffer_size = chunk_buffer_size;
                chunk_args[c].bytes_written = 0;
                chunk_args[c].success = 0;
            }

            /* Submit all formatting tasks */
            for (size_t c = 0; c < num_chunks; c++) {
                ctools_persistent_pool_submit(pool, format_chunk_thread, &chunk_args[c]);
            }

            /* Wait for all formatting to complete */
            ctools_persistent_pool_wait(pool);

            /* Check formatting results */
            for (size_t c = 0; c < num_chunks; c++) {
                if (!chunk_args[c].success) {
                    for (size_t u = 0; u < num_chunks; u++) {
                        if (chunk_buffers[u]) free(chunk_buffers[u]);
                    }
                    free(chunk_args);
                    free(chunk_buffers);
                    free(header_buf);
                    SF_error("cexport: formatting failed\n");
                    cleanup_context();
                    return 920;
                }
            }

            g_ctx.time_format = ctools_timer_seconds() - t_phase;

            /* ============================================================
               STEP 2: Compute prefix-sum for file offsets
               ============================================================ */
            double t_write = ctools_timer_seconds();

            size_t *offsets = (size_t *)ctools_safe_malloc2(num_chunks + 1, sizeof(size_t));
            if (offsets == NULL) {
                for (size_t u = 0; u < num_chunks; u++) {
                    if (chunk_buffers[u]) free(chunk_buffers[u]);
                }
                free(chunk_args);
                free(chunk_buffers);
                free(header_buf);
                SF_error("cexport: memory allocation failed for offsets\n");
                cleanup_context();
                return 920;
            }

            /* First chunk starts after header */
            offsets[0] = header_len;
            for (size_t c = 0; c < num_chunks; c++) {
                offsets[c + 1] = offsets[c] + chunk_args[c].bytes_written;
            }
            size_t total_file_size = offsets[num_chunks];

            /* ============================================================
               STEP 3: Open and pre-size file using cexport_io
               ============================================================ */
            cexport_io_file outfile;
            cexport_io_init(&outfile);

            if (cexport_io_open(&outfile, g_ctx.filename, g_ctx.io_backend, g_ctx.io_flags) != 0) {
                snprintf(msg, sizeof(msg), "cexport: cannot open file '%s': %s\n",
                        g_ctx.filename, outfile.error_message);
                SF_error(msg);
                free(offsets);
                for (size_t u = 0; u < num_chunks; u++) {
                    if (chunk_buffers[u]) free(chunk_buffers[u]);
                }
                free(chunk_args);
                free(chunk_buffers);
                free(header_buf);
                cleanup_context();
                return 603;
            }

            if (cexport_io_presize(&outfile, total_file_size) != 0) {
                snprintf(msg, sizeof(msg), "cexport: failed to pre-size file: %s\n",
                        outfile.error_message);
                SF_error(msg);
                cexport_io_close(&outfile, 0);
                free(offsets);
                for (size_t u = 0; u < num_chunks; u++) {
                    if (chunk_buffers[u]) free(chunk_buffers[u]);
                }
                free(chunk_args);
                free(chunk_buffers);
                free(header_buf);
                cleanup_context();
                return 693;
            }

            /* ============================================================
               STEP 4: Write header at offset 0
               ============================================================ */
            if (header_buf != NULL && header_len > 0) {
                ssize_t written = cexport_io_write_at(&outfile, header_buf, header_len, 0);
                if (written != (ssize_t)header_len) {
                    snprintf(msg, sizeof(msg), "cexport: failed to write header: %s\n",
                            outfile.error_message);
                    SF_error(msg);
                    cexport_io_close(&outfile, 0);
                    free(offsets);
                    for (size_t u = 0; u < num_chunks; u++) {
                        if (chunk_buffers[u]) free(chunk_buffers[u]);
                    }
                    free(chunk_args);
                    free(chunk_buffers);
                    free(header_buf);
                    cleanup_context();
                    return 693;
                }
            }
            free(header_buf);
            header_buf = NULL;

            /* ============================================================
               STEP 5: Parallel offset writes for all chunks
               ============================================================ */
            write_chunk_args_t *write_args = (write_chunk_args_t *)ctools_safe_malloc2(num_chunks, sizeof(write_chunk_args_t));
            if (write_args == NULL) {
                SF_error("cexport: memory allocation failed for write args\n");
                cexport_io_close(&outfile, 0);
                free(offsets);
                for (size_t u = 0; u < num_chunks; u++) {
                    if (chunk_buffers[u]) free(chunk_buffers[u]);
                }
                free(chunk_args);
                free(chunk_buffers);
                cleanup_context();
                return 920;
            }

            /* Setup write arguments for all chunks */
            for (size_t c = 0; c < num_chunks; c++) {
                write_args[c].file = &outfile;
                write_args[c].buffer = chunk_buffers[c];
                write_args[c].len = chunk_args[c].bytes_written;
                write_args[c].offset = offsets[c];
                write_args[c].success = 0;
            }

            /* Submit all write tasks in parallel */
            for (size_t c = 0; c < num_chunks; c++) {
                ctools_persistent_pool_submit(pool, write_chunk_thread, &write_args[c]);
            }

            /* Wait for all writes to complete */
            ctools_persistent_pool_wait(pool);

            /* Check write results */
            int write_failed = 0;
            for (size_t c = 0; c < num_chunks; c++) {
                if (!write_args[c].success) {
                    write_failed = 1;
                    break;
                }
            }

            /* Close file (truncate to actual size) */
            if (cexport_io_close(&outfile, total_file_size) != 0 || write_failed) {
                SF_error("cexport: write error during parallel I/O\n");
                free(write_args);
                free(offsets);
                for (size_t u = 0; u < num_chunks; u++) {
                    if (chunk_buffers[u]) free(chunk_buffers[u]);
                }
                free(chunk_args);
                free(chunk_buffers);
                cleanup_context();
                return 693;
            }

            g_ctx.time_write = ctools_timer_seconds() - t_write;

            /* Cleanup */
            free(write_args);
            free(offsets);
            for (size_t c = 0; c < num_chunks; c++) {
                if (chunk_buffers[c]) free(chunk_buffers[c]);
            }
            free(chunk_args);
            free(chunk_buffers);
        }
    }

    /* Calculate total time */
    g_ctx.time_total = ctools_timer_seconds() - t_start;

    /* Store timing in Stata scalars */
    SF_scal_save("_cexport_time_load", g_ctx.time_load);
    SF_scal_save("_cexport_time_format", g_ctx.time_format);
    SF_scal_save("_cexport_time_write", g_ctx.time_write);
    SF_scal_save("_cexport_time_total", g_ctx.time_total);
    CTOOLS_SAVE_THREAD_INFO("_cexport");

    /* Cleanup and return success */
    cleanup_context();
    return 0;
}

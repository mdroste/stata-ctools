/*
    ctools_types.c
    Implementation of common utility functions for Stata-C data structures
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ctools_types.h"
#include "ctools_config.h"

/*
    NOTE: Use ctools_aligned_free() from ctools_config.h for freeing aligned memory.
    This ensures consistent behavior across all platforms including Windows.
*/
#define aligned_free_internal ctools_aligned_free

void stata_data_init(stata_data *data)
{
    if (data == NULL) return;

    data->nobs = 0;
    data->nvars = 0;
    data->vars = NULL;
    data->sort_order = NULL;
}

/*
    Internal: Free a string arena (matches arena_create in ctools_data_io.c).
    This is duplicated here to avoid circular dependencies.
*/
typedef struct {
    char *base;
    size_t capacity;
    size_t used;
} string_arena_internal;

static void arena_free_internal(void *arena_ptr)
{
    if (arena_ptr != NULL) {
        string_arena_internal *arena = (string_arena_internal *)arena_ptr;
        free(arena->base);
        free(arena);
    }
}

void stata_data_free(stata_data *data)
{
    size_t i;

    if (data == NULL) return;

    /* Free each variable's data (allocated with aligned_alloc_cacheline) */
    if (data->vars != NULL) {
        for (i = 0; i < data->nvars; i++) {
            if (data->vars[i].type == STATA_TYPE_DOUBLE) {
                /* Numeric data uses aligned allocation */
                aligned_free_internal(data->vars[i].data.dbl);
            } else if (data->vars[i].type == STATA_TYPE_STRING) {
                if (data->vars[i].data.str != NULL) {
                    /* Check if strings were allocated via arena (fast path)
                       or individually via strdup (slow path) */
                    if (data->vars[i]._arena != NULL) {
                        /* Arena allocation: O(1) bulk free - just free the arena.
                           All strings are inside the arena's contiguous memory block. */
                        arena_free_internal(data->vars[i]._arena);
                        data->vars[i]._arena = NULL;
                    } else {
                        /* Legacy path: individual strdup allocations - O(n) frees.
                           This path is only used if arena allocation failed. */
                        for (size_t j = 0; j < data->vars[i].nobs; j++) {
                            free(data->vars[i].data.str[j]);
                        }
                    }
                    /* Free the pointer array (aligned allocation) */
                    aligned_free_internal(data->vars[i].data.str);
                }
            }
        }
        /* Free vars array (aligned) */
        aligned_free_internal(data->vars);
    }

    /* Free sort order array (aligned) */
    aligned_free_internal(data->sort_order);

    /* Reset the structure */
    stata_data_init(data);
}

/* ===========================================================================
   Type Conversion Utilities Implementation
   =========================================================================== */

/*
    Power of 10 lookup table for fast float parsing/formatting.
*/
const double ctools_pow10_table[23] = {
    1e0,  1e1,  1e2,  1e3,  1e4,  1e5,  1e6,  1e7,  1e8,  1e9,
    1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19,
    1e20, 1e21, 1e22
};

/*
    Two-digit lookup table for fast digit pair output.
    "00", "01", "02", ... "99"
*/
static const char CTOOLS_DIGIT_PAIRS[200] = {
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

int ctools_uint64_to_str(uint64_t val, char *buf)
{
    if (val == 0) {
        buf[0] = '0';
        buf[1] = '\0';
        return 1;
    }

    char tmp[24];
    int len = 0;

    while (val > 0) {
        tmp[len++] = '0' + (val % 10);
        val /= 10;
    }

    /* Reverse into output buffer */
    for (int i = 0; i < len; i++) {
        buf[i] = tmp[len - 1 - i];
    }
    buf[len] = '\0';
    return len;
}

int ctools_int64_to_str(int64_t val, char *buf)
{
    if (val < 0) {
        buf[0] = '-';
        return 1 + ctools_uint64_to_str((uint64_t)(-val), buf + 1);
    }
    return ctools_uint64_to_str((uint64_t)val, buf);
}

int ctools_uint64_to_str_fast(uint64_t val, char *buf)
{
    if (val == 0) {
        buf[0] = '0';
        buf[1] = '\0';
        return 1;
    }

    char tmp[24];
    int len = 0;

    /* Process two digits at a time using lookup table */
    while (val >= 100) {
        int idx = (val % 100) * 2;
        val /= 100;
        tmp[len++] = CTOOLS_DIGIT_PAIRS[idx + 1];
        tmp[len++] = CTOOLS_DIGIT_PAIRS[idx];
    }

    /* Handle remaining 1-2 digits */
    if (val >= 10) {
        int idx = (int)val * 2;
        tmp[len++] = CTOOLS_DIGIT_PAIRS[idx + 1];
        tmp[len++] = CTOOLS_DIGIT_PAIRS[idx];
    } else {
        tmp[len++] = '0' + (int)val;
    }

    /* Reverse into output buffer */
    for (int i = 0; i < len; i++) {
        buf[i] = tmp[len - 1 - i];
    }
    buf[len] = '\0';
    return len;
}

bool ctools_parse_double_fast(const char *str, int len, double *result, double missval)
{
    const char *p = str;
    const char *end = str + len;

    /* Skip leading whitespace */
    while (p < end && (*p == ' ' || *p == '\t')) p++;
    /* Skip trailing whitespace */
    while (end > p && (end[-1] == ' ' || end[-1] == '\t' || end[-1] == '\r' || end[-1] == '\n')) end--;

    len = (int)(end - p);

    /* Empty string -> missing */
    if (len == 0) {
        *result = missval;
        return true;
    }

    /* Single "." -> Stata missing */
    if (len == 1 && *p == '.') {
        *result = missval;
        return true;
    }

    /* "NA" -> missing */
    if (len == 2 && (p[0] == 'N' || p[0] == 'n') && (p[1] == 'A' || p[1] == 'a')) {
        *result = missval;
        return true;
    }

    /* "NaN" -> missing */
    if (len == 3 && (p[0] == 'N' || p[0] == 'n') && (p[1] == 'a' || p[1] == 'A') && (p[2] == 'N' || p[2] == 'n')) {
        *result = missval;
        return true;
    }

    /* Parse sign */
    bool negative = false;
    if (*p == '-') {
        negative = true;
        p++;
    } else if (*p == '+') {
        p++;
    }

    if (p >= end) return false;

    /* Parse mantissa */
    uint64_t mantissa = 0;
    int decimal_pos = -1;
    int digit_count = 0;
    int total_digits = 0;
    bool has_digits = false;

    while (p < end) {
        char c = *p;
        if (c >= '0' && c <= '9') {
            has_digits = true;
            total_digits++;
            if (digit_count < 18) {
                mantissa = mantissa * 10 + (c - '0');
                digit_count++;
            }
            p++;
        } else if (c == '.' && decimal_pos < 0) {
            decimal_pos = total_digits;
            p++;
        } else if (c == 'e' || c == 'E') {
            /* Scientific notation */
            p++;
            if (p >= end) return false;

            bool exp_negative = false;
            if (*p == '-') {
                exp_negative = true;
                p++;
            } else if (*p == '+') {
                p++;
            }

            int exponent = 0;
            while (p < end && *p >= '0' && *p <= '9') {
                exponent = exponent * 10 + (*p - '0');
                p++;
            }

            if (exp_negative) exponent = -exponent;

            double value = (double)mantissa;
            int decimal_shift = (decimal_pos >= 0) ? (total_digits - decimal_pos) : 0;
            int extra_digits = (total_digits > 18) ? (total_digits - 18) : 0;
            if (decimal_pos >= 0 && decimal_pos < total_digits) {
                extra_digits -= (total_digits - decimal_pos);
                if (extra_digits < 0) extra_digits = 0;
            }

            int final_exp = exponent - decimal_shift + extra_digits;

            if (final_exp > 0 && final_exp <= 22) {
                value *= ctools_pow10_table[final_exp];
            } else if (final_exp < 0 && final_exp >= -22) {
                value /= ctools_pow10_table[-final_exp];
            } else if (final_exp != 0) {
                value *= pow(10.0, final_exp);
            }

            *result = negative ? -value : value;
            return (p == end);
        } else {
            break;
        }
    }

    if (!has_digits) return false;
    if (p != end) return false;

    double value = (double)mantissa;

    /* Apply decimal shift */
    if (decimal_pos >= 0) {
        int decimal_digits = total_digits - decimal_pos;
        if (total_digits > 18 && decimal_pos < 18) {
            decimal_digits -= (total_digits - 18);
            if (decimal_digits < 0) decimal_digits = 0;
        }
        if (decimal_digits > 0 && decimal_digits <= 22) {
            value /= ctools_pow10_table[decimal_digits];
        } else if (decimal_digits > 22) {
            value /= pow(10.0, decimal_digits);
        }
    } else if (total_digits > 18) {
        int extra = total_digits - 18;
        if (extra <= 22) {
            value *= ctools_pow10_table[extra];
        } else {
            value *= pow(10.0, extra);
        }
    }

    *result = negative ? -value : value;
    return true;
}

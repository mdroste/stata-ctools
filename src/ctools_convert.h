#ifndef CTOOLS_CONVERT_H
#define CTOOLS_CONVERT_H
#include "ctools_data.h"

/* ---------------------------------------------------------------------------
   Type Conversion Utilities

   Shared functions for fast numeric/string conversion used across ctools.
   These avoid sprintf/sscanf overhead for performance-critical code paths.
   --------------------------------------------------------------------------- */

/*
    Convert double to sortable uint64 for radix sort.

    IEEE 754 doubles don't sort correctly as raw bit patterns because:
    - Negative numbers have sign bit set but are "less than" positive
    - Negative numbers sort in reverse order when compared as unsigned

    This function transforms the bit pattern so that:
    - Positive numbers: flip sign bit (0x8000... becomes 0x0000...)
    - Negative numbers: flip all bits (preserves ordering)

    Stata's missing values are finite doubles above every nonmissing value
    (. < .a < ... < .z), so they need no special case: the transform keeps
    them after all numbers and in Stata's order. -0 is mapped to +0 because
    Stata treats them as equal.

    @param d  Input double value
    @return   Sortable uint64 representation
*/
static inline uint64_t ctools_double_to_sortable(double d)
{
    uint64_t bits;
    if (d == 0)
        d = 0;
    memcpy(&bits, &d, sizeof(bits));
    return bits ^ ((bits >> 63) ? UINT64_MAX : UINT64_C(0x8000000000000000));
}

/*
    Fast unsigned integer to string conversion.
    Writes digits in reverse order then reverses.

    @param val  Value to convert
    @param buf  Output buffer (must hold at least 21 chars)
    @return     Number of characters written (not including null terminator)
*/
int ctools_uint64_to_str(uint64_t val, char *buf);

/*
    Fast signed integer to string conversion.

    @param val  Value to convert
    @param buf  Output buffer (must hold at least 22 chars for sign + digits)
    @return     Number of characters written (not including null terminator)
*/
int ctools_int64_to_str(int64_t val, char *buf);

/*
    Fast string to double parser.
    Handles integers, decimals, scientific notation, and common missing values.
    Much faster than strtod/atof for typical numeric data.

    Recognized missing value representations:
    - Empty string or whitespace only
    - "." (Stata missing)
    - "NA", "na", "NaN", "nan"

    @param str     Input string (not necessarily null-terminated)
    @param len     Length of input string
    @param result  Output: parsed double value (SV_missval for missing)
    @param missval The value to use for missing (typically SV_missval from stplugin.h)
    @return        true if successfully parsed, false if invalid format
*/
bool ctools_parse_double_fast(const char *str, int len, double *result, double missval);

/*
    Fast double parser with custom decimal and group separators.
    Handles European formats like "1.234,56" (decimal=',', group='.').

    @param str       Input string (need not be null-terminated)
    @param len       Length of input
    @param result    [out] Parsed double value
    @param missval   The value to use for missing
    @param dec_sep   Decimal separator character (default '.')
    @param grp_sep   Group/thousands separator character ('\0' = none)
    @return          true if successfully parsed, false if invalid format
*/
bool ctools_parse_double_with_separators(const char *str, int len, double *result, double missval,
                                         char dec_sep, char grp_sep);

/*
    Power of 10 lookup table for fast float parsing/formatting.
    Covers 10^0 through 10^22 (full double precision range without overflow).
*/
extern const double ctools_pow10_table[23];

/*
    Two-digit lookup table for fast digit pair output.
    "00", "01", "02", ... "99" stored as pairs of characters.
    Usage: buf[0] = CTOOLS_DIGIT_PAIRS[val*2]; buf[1] = CTOOLS_DIGIT_PAIRS[val*2+1];
*/
extern const char CTOOLS_DIGIT_PAIRS[200];

/* ---------------------------------------------------------------------------
   Safe String-to-Number Parsing

   These functions replace unsafe atoi()/atof() which return 0 on parse failure,
   making it impossible to distinguish "0" from invalid input like "abc".
   --------------------------------------------------------------------------- */

#include <stdlib.h>
#include <errno.h>
#include <limits.h>

/*
    Safe string to int conversion.
    Uses strtol internally with full error checking.

    @param str     Input string to parse
    @param result  Output: parsed integer value (unchanged on failure)
    @return        true if successfully parsed, false if:
                   - str is NULL or empty
                   - contains non-numeric characters
                   - value overflows int range
*/
static inline bool ctools_safe_atoi(const char *str, int *result)
{
    if (str == NULL || *str == '\0')
    {
        return false;
    }

    char *endptr;
    errno = 0;
    long val = strtol(str, &endptr, 10);

    /* Check for conversion errors */
    if (errno == ERANGE || val > INT_MAX || val < INT_MIN)
    {
        return false; /* Overflow */
    }
    if (endptr == str)
    {
        return false; /* No digits found */
    }
    /* Allow trailing whitespace but not other characters */
    while (*endptr == ' ' || *endptr == '\t')
    {
        endptr++;
    }
    if (*endptr != '\0')
    {
        return false; /* Trailing non-whitespace characters */
    }

    *result = (int)val;
    return true;
}

/*
    Safe string to size_t conversion.

    @param str     Input string to parse
    @param result  Output: parsed size_t value (unchanged on failure)
    @return        true if successfully parsed, false on error or negative value
*/
static inline bool ctools_safe_atozu(const char *str, size_t *result)
{
    if (str == NULL || *str == '\0')
    {
        return false;
    }

    /* Reject negative numbers for size_t */
    const char *p = str;
    while (*p == ' ' || *p == '\t')
        p++;
    if (*p == '-')
    {
        return false;
    }

    char *endptr;
    errno = 0;
    unsigned long long val = strtoull(str, &endptr, 10);

    if (errno == ERANGE || val > SIZE_MAX)
    {
        return false;
    }
    if (endptr == str)
    {
        return false;
    }
    while (*endptr == ' ' || *endptr == '\t')
    {
        endptr++;
    }
    if (*endptr != '\0')
    {
        return false;
    }

    *result = (size_t)val;
    return true;
}

#endif /* CTOOLS_CONVERT_H */

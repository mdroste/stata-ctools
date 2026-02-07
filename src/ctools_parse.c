/*
 * ctools_parse.c
 * Unified argument parsing utilities for ctools
 *
 * Common parsing functions extracted from cwinsor, cbsample, crangestat,
 * and other ctools commands.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <stdint.h>
#include "ctools_parse.h"
#include "ctools_types.h"  /* For ctools_safe_atof, ctools_safe_atoi, etc. */

/*
 * Check if a boolean option is present in the argument string.
 */
int ctools_parse_bool_option(const char *args, const char *name)
{
    if (!args || !name) return 0;

    const char *p = strstr(args, name);
    if (!p) return 0;

    /* Must be at start of string or preceded by whitespace */
    if (p > args && *(p - 1) != ' ' && *(p - 1) != '\t') return 0;

    /* Must be followed by end-of-string, whitespace, or '=' */
    char after = *(p + strlen(name));
    if (after != '\0' && after != ' ' && after != '\t' && after != '=') return 0;

    return 1;
}

/*
 * Parse a string option of the form "key=value".
 */
int ctools_parse_string_option(const char *args, const char *key,
                                char *buf, size_t bufsize)
{
    if (!args || !key || !buf || bufsize == 0) return 0;

    /* Build pattern "key=" */
    char pat[128];
    size_t keylen = strlen(key);
    if (keylen >= sizeof(pat) - 1) return 0;

    memcpy(pat, key, keylen);
    pat[keylen] = '=';
    pat[keylen + 1] = '\0';

    const char *p = strstr(args, pat);
    if (!p) return 0;

    /* Move past "key=" */
    p += keylen + 1;

    /* Copy value until whitespace or end-of-string */
    size_t i = 0;
    while (*p && *p != ' ' && *p != '\t' && i < bufsize - 1) {
        buf[i++] = *p++;
    }
    buf[i] = '\0';

    return (i > 0) ? 1 : 0;
}

/*
 * Parse an integer option of the form "key=value".
 */
int ctools_parse_int_option(const char *args, const char *key, int *value)
{
    if (!args || !key || !value) return 0;

    /* Build pattern "key=" */
    char pat[128];
    size_t keylen = strlen(key);
    if (keylen >= sizeof(pat) - 1) return 0;

    snprintf(pat, sizeof(pat), "%s=", key);

    const char *p = strstr(args, pat);
    if (!p) return 0;

    /* Move past "key=" */
    p += keylen + 1;

    /* Parse integer value (stops at first non-digit character).
     * We use strtol directly instead of ctools_safe_atoi because the value
     * may be followed by other options (e.g., "nvars=1 force"). */
    char *endptr;
    errno = 0;
    long result = strtol(p, &endptr, 10);
    if (endptr == p || errno == ERANGE || result > INT_MAX || result < INT_MIN) {
        return 0;
    }

    *value = (int)result;
    return 1;
}

/*
 * Parse a double option of the form "key=value", with default.
 */
double ctools_parse_double_option(const char *args, const char *key,
                                   double default_val)
{
    if (!args || !key) return default_val;

    /* Build pattern "key=" */
    char pat[128];
    snprintf(pat, sizeof(pat), "%s=", key);

    const char *p = strstr(args, pat);
    if (!p) return default_val;

    /* Move past "key=" */
    p += strlen(pat);

    /* Parse double value (may be followed by other options, so use strtod) */
    char *endptr;
    errno = 0;
    double result = strtod(p, &endptr);
    if (endptr == p || errno == ERANGE) {
        return default_val;
    }

    return result;
}

/*
 * Parse a size_t option of the form "key=value", with default.
 */
size_t ctools_parse_size_option(const char *args, const char *key,
                                 size_t default_val)
{
    if (!args || !key) return default_val;

    /* Build pattern "key=" */
    char pat[128];
    snprintf(pat, sizeof(pat), "%s=", key);

    const char *p = strstr(args, pat);
    if (!p) return default_val;

    /* Move past "key=" */
    p += strlen(pat);

    /* Parse size_t value (may be followed by other options, so use strtoul) */
    char *endptr;
    errno = 0;
    unsigned long val = strtoul(p, &endptr, 10);
    if (endptr == p || errno == ERANGE) {
        return default_val;
    }

    return (size_t)val;
}

/*
 * Parse a sequence of whitespace-separated integers from cursor position.
 */
int ctools_parse_int_array(int *arr, size_t count, const char **cursor)
{
    if (!arr || !cursor || !*cursor) return -1;

    const char *p = *cursor;

    for (size_t i = 0; i < count; i++) {
        /* Skip whitespace */
        while (*p == ' ' || *p == '\t') p++;

        if (*p == '\0') return -1;  /* Not enough integers */

        /* Parse integer */
        char *end;
        long val = strtol(p, &end, 10);
        if (end == p) return -1;  /* No digits found */

        arr[i] = (int)val;
        p = end;
    }

    *cursor = p;
    return 0;
}

/*
 * Skip leading whitespace and parse a single size_t value.
 */
int ctools_parse_next_size(const char **cursor, size_t *value)
{
    if (!cursor || !*cursor || !value) return -1;

    const char *p = *cursor;

    /* Skip whitespace */
    while (*p == ' ' || *p == '\t') p++;

    if (*p == '\0') return -1;

    /* Parse unsigned integer */
    char *end;
    unsigned long long val = strtoull(p, &end, 10);
    if (end == p) return -1;

    *value = (size_t)val;
    *cursor = end;
    return 0;
}

/*
 * Skip leading whitespace and parse a single int value.
 */
int ctools_parse_next_int(const char **cursor, int *value)
{
    if (!cursor || !*cursor || !value) return -1;

    const char *p = *cursor;

    /* Skip whitespace */
    while (*p == ' ' || *p == '\t') p++;

    if (*p == '\0') return -1;

    /* Parse integer */
    char *end;
    long val = strtol(p, &end, 10);
    if (end == p) return -1;

    *value = (int)val;
    *cursor = end;
    return 0;
}

/*
 * ctools_parse.h
 * Unified argument parsing utilities for ctools
 *
 * Provides common parsing functions used across ctools commands:
 * - Boolean option parsing (e.g., "verbose", "force")
 * - String option parsing (e.g., "name=value")
 * - Integer/double/size_t option parsing with defaults
 * - Integer array parsing from argument strings
 *
 * Usage:
 *   int verbose = ctools_parse_bool_option(args, "verbose");
 *   double pctl = ctools_parse_double_option(args, "p", 1.0);
 */

#ifndef CTOOLS_PARSE_H
#define CTOOLS_PARSE_H

#include <stddef.h>

/*
 * Check if a boolean option is present in the argument string.
 * Returns 1 if found as a standalone token, 0 otherwise.
 *
 * The option must be:
 * - At the start of args, OR preceded by whitespace
 * - Followed by end-of-string, whitespace, or '='
 *
 * @param args  The argument string to search
 * @param name  The option name to find
 * @return      1 if option is present, 0 otherwise
 *
 * Example:
 *   ctools_parse_bool_option("verbose trim", "verbose") -> 1
 *   ctools_parse_bool_option("myverbose", "verbose") -> 0
 */
int ctools_parse_bool_option(const char *args, const char *name);

/*
 * Parse a string option of the form "key=value".
 * Returns 1 if found and copied, 0 otherwise.
 *
 * Value extraction stops at first whitespace character.
 * The result is null-terminated and truncated if necessary.
 *
 * @param args     The argument string to search
 * @param key      The option key (without '=')
 * @param buf      Buffer to store the value
 * @param bufsize  Size of the buffer
 * @return         1 if option found and parsed, 0 otherwise
 *
 * Example:
 *   char name[64];
 *   ctools_parse_string_option("name=foo other", "name", name, sizeof(name));
 *   // name = "foo"
 */
int ctools_parse_string_option(const char *args, const char *key,
                                char *buf, size_t bufsize);

/*
 * Parse an integer option of the form "key=value".
 * Returns 1 if found and valid, 0 otherwise.
 *
 * @param args   The argument string to search
 * @param key    The option key (without '=')
 * @param value  Pointer to store the parsed value (unchanged if not found)
 * @return       1 if option found and parsed, 0 otherwise
 *
 * Example:
 *   int n = 10;  // default
 *   ctools_parse_int_option("count=5", "count", &n);
 *   // n = 5
 */
int ctools_parse_int_option(const char *args, const char *key, int *value);

/*
 * Parse a double option of the form "key=value", with default.
 * Returns the parsed value if found, or default_val otherwise.
 *
 * @param args        The argument string to search
 * @param key         The option key (without '=')
 * @param default_val Value to return if option not found or invalid
 * @return            Parsed value or default
 *
 * Example:
 *   double pctl = ctools_parse_double_option(args, "p", 1.0);
 */
double ctools_parse_double_option(const char *args, const char *key,
                                   double default_val);

/*
 * Parse a size_t option of the form "key=value", with default.
 * Returns the parsed value if found, or default_val otherwise.
 * Rejects negative values.
 *
 * @param args        The argument string to search
 * @param key         The option key (without '=')
 * @param default_val Value to return if option not found or invalid
 * @return            Parsed value or default
 */
size_t ctools_parse_size_option(const char *args, const char *key,
                                 size_t default_val);

/*
 * Parse a sequence of whitespace-separated integers from cursor position.
 * Advances cursor past the parsed integers.
 *
 * @param arr    Array to store parsed integers
 * @param count  Number of integers to parse
 * @param cursor Pointer to cursor position in string (updated on success)
 * @return       0 on success, -1 if fewer than count integers found
 *
 * Example:
 *   const char *p = "1 2 3 verbose";
 *   int arr[3];
 *   ctools_parse_int_array(arr, 3, &p);
 *   // arr = {1, 2, 3}, p points to " verbose"
 */
int ctools_parse_int_array(int *arr, size_t count, const char **cursor);

/*
 * Skip leading whitespace and parse a single size_t value.
 * Advances cursor past the parsed value.
 *
 * @param cursor Pointer to cursor position in string (updated on success)
 * @param value  Pointer to store the parsed value
 * @return       0 on success, -1 on failure
 */
int ctools_parse_next_size(const char **cursor, size_t *value);

/*
 * Skip leading whitespace and parse a single int value.
 * Advances cursor past the parsed value.
 *
 * @param cursor Pointer to cursor position in string (updated on success)
 * @param value  Pointer to store the parsed value
 * @return       0 on success, -1 on failure
 */
int ctools_parse_next_int(const char **cursor, int *value);

#endif /* CTOOLS_PARSE_H */

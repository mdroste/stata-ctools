/*
 * ctools_args.h - Argument parsing utilities for ctools commands
 *
 * Provides common utilities for parsing command-line arguments passed
 * from Stata to C plugin functions.
 *
 * Part of the ctools suite for Stata.
 */

#ifndef CTOOLS_ARGS_H
#define CTOOLS_ARGS_H

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/*
 * Find a named option in an argument string.
 * Returns pointer to the value after "name=" or NULL if not found.
 *
 * Example:
 *   const char *val = ctools_args_find_option(args, "alg=");
 *   if (val) printf("Algorithm: %s\n", val);
 */
static inline const char *ctools_args_find_option(const char *args, const char *name)
{
    const char *p = strstr(args, name);
    if (p == NULL) return NULL;
    return p + strlen(name);
}

/*
 * Check if a flag is present in the argument string.
 * Matches whole word only (space or end of string delimited).
 *
 * Example:
 *   if (ctools_args_has_flag(args, "verbose")) { ... }
 */
static inline bool ctools_args_has_flag(const char *args, const char *flag)
{
    const char *p = args;
    size_t len = strlen(flag);

    while ((p = strstr(p, flag)) != NULL) {
        /* Check if this is a whole word match */
        bool start_ok = (p == args || p[-1] == ' ' || p[-1] == '\t');
        bool end_ok = (p[len] == '\0' || p[len] == ' ' || p[len] == '\t');
        if (start_ok && end_ok) return true;
        p += len;
    }
    return false;
}

/*
 * Parse an integer array from space-separated tokens.
 * Stops parsing at the specified stop_at option (e.g., "alg=").
 * Returns number of integers parsed, or -1 on allocation failure.
 *
 * Example:
 *   int *vars;
 *   size_t nvars;
 *   if (ctools_args_parse_int_array(args, "alg=", &vars, &nvars) == 0) {
 *       // use vars[0..nvars-1]
 *       free(vars);
 *   }
 */
static inline int ctools_args_parse_int_array(const char *args, const char *stop_at,
                                               int **out_arr, size_t *out_count)
{
    const char *p = args;
    const char *stop_ptr = stop_at ? strstr(args, stop_at) : NULL;
    char *endptr;
    size_t count = 0;
    size_t capacity = 16;
    int *arr;
    long val;

    arr = (int *)malloc(capacity * sizeof(int));
    if (arr == NULL) return -1;

    while (*p != '\0') {
        /* Stop at the specified option */
        if (stop_ptr != NULL && p >= stop_ptr) break;

        /* Skip whitespace */
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0') break;
        if (stop_ptr != NULL && p >= stop_ptr) break;

        /* Parse integer */
        val = strtol(p, &endptr, 10);
        if (endptr == p) {
            p++;
            continue;
        }

        /* Grow array if needed */
        if (count >= capacity) {
            capacity *= 2;
            int *new_arr = (int *)realloc(arr, capacity * sizeof(int));
            if (new_arr == NULL) {
                free(arr);
                return -1;
            }
            arr = new_arr;
        }

        arr[count++] = (int)val;
        p = endptr;
    }

    *out_arr = arr;
    *out_count = count;
    return 0;
}

/*
 * Parse a single integer option value.
 * Returns the integer value or default_val if not found.
 *
 * Example:
 *   int nthreads = ctools_args_get_int(args, "threads=", 4);
 */
static inline int ctools_args_get_int(const char *args, const char *name, int default_val)
{
    const char *p = ctools_args_find_option(args, name);
    if (p == NULL) return default_val;
    return (int)strtol(p, NULL, 10);
}

/*
 * Copy arguments to a mutable buffer safely.
 * Returns 0 on success, -1 if args too long.
 *
 * Example:
 *   char buf[4096];
 *   if (ctools_args_copy(args, buf, sizeof(buf)) == 0) { ... }
 */
static inline int ctools_args_copy(const char *args, char *buf, size_t buf_size)
{
    size_t len = strlen(args);
    if (len >= buf_size) return -1;
    memcpy(buf, args, len + 1);
    return 0;
}

#endif /* CTOOLS_ARGS_H */

/* Checked, bounded parsing of the private ado/plugin argument protocol. */
#ifndef CTOOLS_PARSE_H
#define CTOOLS_PARSE_H
#include <stddef.h>
#include <stdint.h>

typedef enum { CTOOLS_PARSE_INVALID = -1, CTOOLS_PARSE_ABSENT = 0,
               CTOOLS_PARSE_OK = 1 } ctools_parse_status;
typedef struct { const char *next, *end; } ctools_arg_cursor;
typedef struct { const char *data; size_t length; } ctools_arg_token;
/* Tokens are whitespace delimited; double quotes protect spaces. No truncation.
 * The cursor never reads beyond length. Unbalanced quotes are invalid. */
void ctools_args_init(ctools_arg_cursor *cursor, const char *args, size_t length);
ctools_parse_status ctools_args_next(ctools_arg_cursor *cursor, ctools_arg_token *token);
/* Named values use key=value, optionally double quoted. Duplicate named values
 * are rejected; flags may repeat. Outputs change only on success. */
ctools_parse_status ctools_parse_string_option(const char *, const char *, char *, size_t);
ctools_parse_status ctools_parse_int_option(const char *, const char *, int *);
ctools_parse_status ctools_parse_u64_option(const char *, const char *, uint64_t *);
ctools_parse_status ctools_parse_double_checked(const char *, const char *, double *);
ctools_parse_status ctools_parse_seed_option(const char *, uint64_t *);
int ctools_parse_bool_option(const char *, const char *);
/* Compatibility convenience: default only when absent; invalid input is NaN. */
double ctools_parse_double_option(const char *, const char *, double);
/* Positional integers: 0 on success, -1 on absent/invalid, cursor unchanged on
 * failure. Integers must occupy complete tokens and fit in int. */
int ctools_parse_next_int(const char **, int *);
int ctools_parse_int_array(int *, size_t, const char **);
#endif

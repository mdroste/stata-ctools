#include "ctools_parse.h"
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void ctools_args_init(ctools_arg_cursor *c, const char *args, size_t n)
{
    c->next = args;
    c->end = args ? args + n : NULL;
}

ctools_parse_status ctools_args_next(ctools_arg_cursor *c, ctools_arg_token *t)
{
    if (!c || !t || !c->next || !c->end) return CTOOLS_PARSE_INVALID;
    const char *p = c->next;
    while (p < c->end && isspace((unsigned char)*p)) p++;
    if (p == c->end) return CTOOLS_PARSE_ABSENT;
    const char *start = p;
    int quoted = 0;
    while (p < c->end) {
        if (*p == '"') quoted = !quoted;
        else if (!quoted && isspace((unsigned char)*p)) break;
        p++;
    }
    if (quoted) return CTOOLS_PARSE_INVALID;
    t->data = start;
    t->length = (size_t)(p - start);
    c->next = p;
    return CTOOLS_PARSE_OK;
}

static ctools_parse_status option(const char *args, const char *key, ctools_arg_token *value)
{
    if (!args || !key || !*key || !value) return CTOOLS_PARSE_INVALID;
    size_t n = strlen(key);
    ctools_arg_cursor c;
    ctools_arg_token t, found = {0};
    ctools_args_init(&c, args, strlen(args));
    int rc;
    while ((rc = ctools_args_next(&c, &t)) == CTOOLS_PARSE_OK) {
        if (t.length < n || memcmp(t.data, key, n)) continue;
        if (t.length == n) return CTOOLS_PARSE_INVALID; /* value required */
        if (t.data[n] != '=') continue;
        if (found.data) return CTOOLS_PARSE_INVALID;
        found.data = t.data + n + 1;
        found.length = t.length - n - 1;
        if (found.length && found.data[0] == '"') {
            if (found.length < 2 || found.data[found.length - 1] != '"') return CTOOLS_PARSE_INVALID;
            found.data++; found.length -= 2;
        }
    }
    if (rc == CTOOLS_PARSE_INVALID) return rc;
    if (!found.data) return CTOOLS_PARSE_ABSENT;
    *value = found;
    return CTOOLS_PARSE_OK;
}

int ctools_parse_bool_option(const char *args, const char *name)
{
    if (!args || !name || !*name) return 0;
    size_t n = strlen(name);
    ctools_arg_cursor c;
    ctools_arg_token t;
    ctools_args_init(&c, args, strlen(args));
    int found = 0, rc;
    while ((rc = ctools_args_next(&c, &t)) == CTOOLS_PARSE_OK) {
        if (t.length >= n && !memcmp(t.data, name, n) &&
            (t.length == n || t.data[n] == '=')) found = 1;
    }
    return rc == CTOOLS_PARSE_INVALID ? 0 : found;
}

ctools_parse_status ctools_parse_string_option(const char *args, const char *key, char *buf, size_t cap)
{
    if (!buf || !cap) return CTOOLS_PARSE_INVALID;
    ctools_arg_token t;
    int rc = option(args, key, &t);
    if (rc != CTOOLS_PARSE_OK) return rc;
    if (t.length >= cap) return CTOOLS_PARSE_INVALID;
    memcpy(buf, t.data, t.length);
    buf[t.length] = 0;
    return CTOOLS_PARSE_OK;
}

static int token_int(ctools_arg_token t, int *out)
{
    char buf[64], *end;
    if (!out || !t.length || t.length >= sizeof(buf)) return -1;
    memcpy(buf, t.data, t.length); buf[t.length] = 0;
    if (isspace((unsigned char)buf[0])) return -1;
    errno = 0;
    long v = strtol(buf, &end, 10);
    if (end != buf + t.length || errno == ERANGE || v < INT_MIN || v > INT_MAX) return -1;
    *out = (int)v;
    return 0;
}

ctools_parse_status ctools_parse_int_option(const char *args, const char *key, int *out)
{
    ctools_arg_token t;
    int rc = option(args, key, &t);
    if (rc != CTOOLS_PARSE_OK) return rc;
    return token_int(t, out) ? CTOOLS_PARSE_INVALID : CTOOLS_PARSE_OK;
}

ctools_parse_status ctools_parse_u64_option(const char *args, const char *key, uint64_t *out)
{
    ctools_arg_token t;
    int rc = option(args, key, &t);
    if (rc != CTOOLS_PARSE_OK) return rc;
    if (!out || !t.length) return CTOOLS_PARSE_INVALID;
    uint64_t v = 0;
    for (size_t i = 0; i < t.length; i++) {
        unsigned digit = (unsigned char)t.data[i] - '0';
        if (digit > 9 || v > (UINT64_MAX - digit) / 10) return CTOOLS_PARSE_INVALID;
        v = v * 10 + digit;
    }
    *out = v;
    return CTOOLS_PARSE_OK;
}

ctools_parse_status ctools_parse_double_checked(const char *args, const char *key, double *out)
{
    ctools_arg_token t;
    int rc = option(args, key, &t);
    if (rc != CTOOLS_PARSE_OK) return rc;
    char buf[128], *end;
    if (!out || !t.length || t.length >= sizeof(buf)) return CTOOLS_PARSE_INVALID;
    memcpy(buf, t.data, t.length); buf[t.length] = 0;
    if (isspace((unsigned char)buf[0])) return CTOOLS_PARSE_INVALID;
    errno = 0;
    double v = strtod(buf, &end);
    if (end != buf + t.length || errno == ERANGE || !isfinite(v)) return CTOOLS_PARSE_INVALID;
    *out = v;
    return CTOOLS_PARSE_OK;
}

double ctools_parse_double_option(const char *args, const char *key, double fallback)
{
    double v = fallback;
    return ctools_parse_double_checked(args, key, &v) < 0 ? NAN : v;
}

ctools_parse_status ctools_parse_seed_option(const char *args, uint64_t *seed)
{
    uint64_t hi = 0, lo = 0;
    int a = ctools_parse_u64_option(args, "seedhi", &hi);
    int b = ctools_parse_u64_option(args, "seedlo", &lo);
    if (!seed || a < 0 || b < 0 || a != b || hi > 0x7fffffffULL || lo > 0x7fffffffULL)
        return CTOOLS_PARSE_INVALID;
    if (!a) return CTOOLS_PARSE_ABSENT;
    *seed = (hi << 31) | lo;
    return CTOOLS_PARSE_OK;
}

int ctools_parse_next_int(const char **cursor, int *value)
{
    if (!cursor || !*cursor) return -1;
    ctools_arg_cursor c;
    ctools_arg_token t;
    ctools_args_init(&c, *cursor, strlen(*cursor));
    if (ctools_args_next(&c, &t) != CTOOLS_PARSE_OK || token_int(t, value)) return -1;
    *cursor = c.next;
    return 0;
}

int ctools_parse_int_array(int *arr, size_t n, const char **cursor)
{
    if ((!arr && n) || !cursor || !*cursor) return -1;
    const char *p = *cursor;
    for (size_t i = 0; i < n; i++) if (ctools_parse_next_int(&p, arr + i)) return -1;
    *cursor = p;
    return 0;
}

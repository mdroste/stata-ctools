/*
    cfuzzyjoin_impl.c
    C-accelerated fuzzy joins: plugin calls and shared helpers

    cfuzzyjoin link   reclink's record linkage    (cfuzzyjoin_link.c)
    cfuzzyjoin match  matchit's similarity joins  (cfuzzyjoin_match.c)
    cfuzzyjoin group  strgroup's string groups    (cfuzzyjoin_group.c)

    The first word of the plugin argument names the call; the files above
    document each. Every call saves its counts and timings in scalars named
    _cfuzzyjoin_*, which the ado file reads.
*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfuzzyjoin_impl.h"
#include "cfuzzyjoin_internal.h"
#include "ctools_runtime.h"

void cfuzzyjoin_cleanup_state(void)
{
    cfj_match_cleanup();
    cfj_link_cleanup();
}

stata_retcode cfj_load(ctools_filtered_data *fd, int first, int count, int flags)
{
    if (count < 1 || count > 4096)
        return STATA_ERR_INVALID_INPUT;
    int idx[4096];
    for (int j = 0; j < count; j++)
        idx[j] = first + j;
    ctools_filtered_data_free(fd);
    return ctools_data_load_ex(fd, idx, (size_t)count, 0, 0, flags | CTOOLS_LOAD_NO_SORT_ORDER,
                               NULL);
}

void cfj_say(const char *fmt, ...)
{
    char buf[1024];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf) - 2, fmt, ap);
    va_end(ap);
    /* one line per call: no newline before the end */
    for (char *p = buf; *p; p++)
        if (*p == '\n' || *p == '\r')
            *p = ' ';
    if (!strncmp(buf, "{err}", 5))
    {
        /* SF_error prints SMCL as typed, and shows under quietly */
        char *text = buf + 5;
        size_t l = strlen(text);
        text[l] = '\n';
        text[l + 1] = '\0';
        SF_error(text);
        return;
    }
    size_t l = strlen(buf);
    buf[l] = '\n';
    buf[l + 1] = '\0';
    SF_display(buf);
}

/* Start of option `name` as a whole word followed by '(' or '=' */
static const char *find_opt(const char *args, const char *name, char *delim)
{
    size_t n = strlen(name);
    for (const char *p = args; p && *p; p++)
    {
        if ((p == args || p[-1] == ' ' || p[-1] == '\t') && !strncmp(p, name, n) &&
            (p[n] == '(' || p[n] == '='))
        {
            *delim = p[n];
            return p + n + 1;
        }
    }
    return NULL;
}

int cfj_opt_text(const char *args, const char *name, char *out, size_t size)
{
    char delim = 0;
    const char *p = find_opt(args, name, &delim);
    if (!p || !size)
        return 0;
    size_t n = 0;
    if (delim == '(')
    {
        int depth = 1;
        while (*p)
        {
            if (*p == '(')
                depth++;
            else if (*p == ')' && --depth == 0)
                break;
            if (n + 1 < size)
                out[n++] = *p;
            p++;
        }
    }
    else
    {
        while (*p && *p != ' ' && *p != '\t')
        {
            if (n + 1 < size)
                out[n++] = *p;
            p++;
        }
    }
    out[n] = '\0';
    return 1;
}

int cfj_opt_double(const char *args, const char *name, double *out)
{
    char text[64];
    char *end;
    if (!cfj_opt_text(args, name, text, sizeof(text)))
        return -1;
    *out = strtod(text, &end);
    return end == text ? -1 : 0;
}

int cfj_opt_int(const char *args, const char *name, int *out)
{
    char text[32];
    char *end;
    if (!cfj_opt_text(args, name, text, sizeof(text)))
        return -1;
    long v = strtol(text, &end, 10);
    if (end == text || *end)
        return -1;
    *out = (int)v;
    return 0;
}

int cfj_flag(const char *args, const char *word)
{
    size_t n = strlen(word);
    for (const char *p = args; p && *p; p++)
    {
        if ((p == args || p[-1] == ' ' || p[-1] == '\t') && !strncmp(p, word, n) &&
            (p[n] == '\0' || p[n] == ' ' || p[n] == '\t'))
            return 1;
    }
    return 0;
}

ST_retcode cfuzzyjoin_main(const char *args)
{
    if (!args)
        return 198;
    while (*args == ' ' || *args == '\t')
        args++;
    char op[32];
    size_t n = 0;
    while (args[n] && args[n] != ' ' && args[n] != '\t' && n + 1 < sizeof(op))
    {
        op[n] = args[n];
        n++;
    }
    op[n] = '\0';
    const char *rest = args + n;
    while (*rest == ' ' || *rest == '\t')
        rest++;

    if (!strcmp(op, "clear"))
    {
        cfuzzyjoin_cleanup_state();
        return 0;
    }
    if (!strcmp(op, "match_master"))
        return cfj_match_master(rest);
    if (!strcmp(op, "match_weights"))
        return cfj_match_weights(rest, 0);
    if (!strcmp(op, "match_mweights"))
        return cfj_match_weights(rest, 1);
    if (!strcmp(op, "match_using"))
        return cfj_match_using(rest);
    if (!strcmp(op, "match_store"))
        return cfj_match_store(rest);
    if (!strcmp(op, "match_columns"))
        return cfj_match_columns(rest);
    if (!strcmp(op, "link_master"))
        return cfj_link_master(rest);
    if (!strcmp(op, "link_exclude"))
        return cfj_link_exclude(rest);
    if (!strcmp(op, "link_using"))
        return cfj_link_using(rest);
    if (!strcmp(op, "link_layout"))
        return cfj_link_layout(rest);
    if (!strcmp(op, "link_store"))
        return cfj_link_store(rest);
    if (!strcmp(op, "group"))
        return cfj_group(rest);
    return 198;
}

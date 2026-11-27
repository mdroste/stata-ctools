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
#include "ctools_alloc.h"
#include "ctools_runtime.h"

void cfuzzyjoin_cleanup_state(void)
{
    cfj_match_cleanup();
    cfj_link_cleanup();
}

/* The longest string the shared loader reads; cfuzzyjoin compares strL
 * values on their first CFJ_STR_MAXLEN bytes */
#define CFJ_STR_MAXLEN 2045

/* Does strL variable idx hold a text value longer than CFJ_STR_MAXLEN bytes
 * (and no binary value, which stays an error of the loader) among the
 * observations in range? */
static int strl_has_long(int idx)
{
    const ST_int in1 = SF_in1(), in2 = SF_in2();
    int found = 0;
    for (ST_int obs = in1; obs >= 1 && obs <= in2; obs++)
    {
        if (SF_var_is_binary(idx, obs))
            return 0;
        if (SF_sdatalen(idx, obs) > CFJ_STR_MAXLEN)
            found = 1;
    }
    return found;
}

/* Fill column j of fd with strL variable idx, each value cut to its first
 * CFJ_STR_MAXLEN bytes (the column the loader read in its place is freed) */
static stata_retcode load_strl_cut(ctools_filtered_data *fd, size_t j, int idx)
{
    stata_variable *v = &fd->data.vars[j];
    const size_t n = fd->data.nobs;
    stata_variable_free(v);
    v->type = STATA_TYPE_STRING;
    v->nobs = n;
    v->buffer_storage = CTOOLS_BUFFER_HEAP;
    v->string_storage = CTOOLS_STRINGS_OWNED;
    v->str_maxlen = CFJ_STR_MAXLEN;
    v->data.str = (char **)calloc(n ? n : 1, sizeof(char *));
    if (!v->data.str)
        return STATA_ERR_MEMORY;
    char *buf = NULL;
    size_t cap = 0;
    stata_retcode rc = STATA_OK;
    for (size_t i = 0; i < n && !rc; i++)
    {
        const ST_int obs = fd->obs_map ? (ST_int)fd->obs_map[i] : SF_in1() + (ST_int)i;
        ST_int len = SF_sdatalen(idx, obs);
        if (len < 0)
        {
            rc = STATA_ERR_STATA_READ;
            break;
        }
        if ((size_t)len + 1 > cap)
        {
            char *more = (char *)realloc(buf, (size_t)len + 1);
            if (!more)
            {
                rc = STATA_ERR_MEMORY;
                break;
            }
            buf = more;
            cap = (size_t)len + 1;
        }
        /* SF_strldata returns the number of bytes copied */
        if (SF_strldata(idx, obs, buf, len + 1) != len)
        {
            rc = STATA_ERR_STATA_READ;
            break;
        }
        size_t keep = (size_t)len < CFJ_STR_MAXLEN ? (size_t)len : CFJ_STR_MAXLEN;
        char *s = (char *)malloc(keep + 1);
        if (!s)
        {
            rc = STATA_ERR_MEMORY;
            break;
        }
        memcpy(s, buf, keep);
        s[keep] = '\0';
        v->data.str[i] = s;
    }
    free(buf);
    return rc;
}

stata_retcode cfj_load(ctools_filtered_data *fd, int first, int count, int flags)
{
    if (count < 1 || count > 4096)
        return STATA_ERR_INVALID_INPUT;
    int idx[4096];
    for (int j = 0; j < count; j++)
        idx[j] = first + j;
    ctools_filtered_data_free(fd);

    /* strL variables with values the loader does not read (over 2,045
     * bytes): the loader reads a stand-in variable in their place, then
     * they are read here cut to 2,045 bytes */
    unsigned char cut[4096];
    int ncut = 0, standin = 0;
    for (int j = 0; j < count; j++)
    {
        cut[j] = SF_var_is_strl(idx[j]) && strl_has_long(idx[j]);
        ncut += cut[j];
    }
    if (ncut)
    {
        for (int v = 1; v <= SF_nvars() && !standin; v++)
            if (!SF_var_is_string(v))
                standin = v;
        if (!standin)
            ncut = 0; /* no numeric variable to stand in: the loader's error */
        for (int j = 0; j < count && ncut; j++)
            if (cut[j])
                idx[j] = standin;
    }

    stata_retcode rc = ctools_data_load_ex(fd, idx, (size_t)count, 0, 0,
                                           flags | CTOOLS_LOAD_NO_SORT_ORDER, NULL);
    if (rc)
        return rc;

    /* No observation: the loader leaves the columns out (or untyped); give
     * them their types, so that the calls can check them and see 0 rows */
    if (fd->data.nobs == 0)
    {
        if (!fd->data.vars)
        {
            fd->data.vars =
                (stata_variable *)ctools_safe_cacheline_alloc2((size_t)count, sizeof(stata_variable));
            if (!fd->data.vars)
                return STATA_ERR_MEMORY;
            memset(fd->data.vars, 0, (size_t)count * sizeof(stata_variable));
            fd->data.nvars = (size_t)count;
        }
        for (int j = 0; j < count && (size_t)j < fd->data.nvars; j++)
            fd->data.vars[j].type =
                SF_var_is_string(first + j) ? STATA_TYPE_STRING : STATA_TYPE_DOUBLE;
        return STATA_OK;
    }

    for (int j = 0; j < count && ncut && !rc; j++)
        if (cut[j])
            rc = load_strl_cut(fd, (size_t)j, first + j);
    if (rc)
        ctools_filtered_data_free(fd);
    return rc;
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

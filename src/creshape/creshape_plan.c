/*
 * creshape_plan.c - the roles of the variables and the variables of the result.
 *
 * The ado file hands over what only Stata knows (the variable names and the
 * request, in locals); everything reshape derives from them is worked out
 * here: which variables carry which j value, each stub's variables, the
 * variables that are constant within i, and the result's variable list.
 * Messages follow reshape's wording.
 */
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
#include "creshape_internal.h"
#include "ctools_config.h"
#include "ctools_parse.h"

/* ------------------------------------------------------------------------
 * Talking to the ado file
 * ------------------------------------------------------------------------ */

char *creshape_local(const char *name)
{
    char spi[160];
    snprintf(spi, sizeof(spi), "_%s", name);
    size_t cap = 4096;
    for (;;)
    {
        char *buf = (char *)malloc(cap);
        if (!buf)
            return NULL;
        buf[0] = buf[cap - 1] = '\0';
        /* an unset local reads as empty */
        if (SF_macro_use(spi, buf, (ST_int)(cap - 1)) != 0)
            buf[0] = '\0';
        /* a value that reaches the end of the buffer may have been cut */
        if (strlen(buf) + 2 < cap)
            return buf;
        free(buf);
        if (cap > (size_t)256 * 1024 * 1024)
            return NULL;
        cap *= 2;
    }
}

/*
 * Messages. Notes go through SF_display, which reads SMCL. Errors go through
 * SF_error, so that they show under -quietly-, and SF_error prints its text
 * as it is: the messages below are written with a few SMCL directives for
 * readability and are turned into plain text here ({bf:x} and {it:x} leave
 * x, {p a b c} starts a paragraph indented a and then b columns, {break}
 * ends a line). Both functions take one line per call: a string that starts
 * with a newline, or holds one before its end, makes Stata print blank lines
 * without end (StataNow 19.5).
 */
#define CRESHAPE_LINE_WIDTH 79

static void creshape_error_line(const char *text, size_t len)
{
    char line[512];
    if (len > sizeof(line) - 2)
        len = sizeof(line) - 2;
    if (!len)
    {
        SF_error(" \n");
        return;
    }
    memcpy(line, text, len);
    line[len] = '\n';
    line[len + 1] = '\0';
    SF_error(line);
}

/* One line of marked-up text: strip the directives and wrap a paragraph. */
static void creshape_error_plain(const char *text, size_t len)
{
    char *plain = (char *)malloc(len + 2);
    if (!plain)
        return;
    int first = -1, rest = 0;
    size_t n = 0;
    for (size_t i = 0; i < len;)
    {
        if (text[i] != '{')
        {
            plain[n++] = text[i++];
            continue;
        }
        const char *close = (const char *)memchr(text + i, '}', len - i);
        if (!strncmp(text + i, "{bf:", 4) || !strncmp(text + i, "{it:", 4))
        {
            /* keep the argument, drop its closing brace */
            size_t depth = 1, k = i + 4;
            for (; k < len && depth; k++)
            {
                if (text[k] == '{')
                    depth++;
                else if (text[k] == '}' && --depth == 0)
                    break;
                plain[n++] = text[k];
            }
            i = k < len ? k + 1 : len;
        }
        else if (!strncmp(text + i, "{p ", 3) && close)
        {
            first = atoi(text + i + 3);
            const char *second =
                (const char *)memchr(text + i + 3, ' ', (size_t)(close - text) - i - 3);
            rest = second ? atoi(second + 1) : first;
            i = (size_t)(close - text) + 1;
        }
        else if (!strncmp(text + i, "{p_end}", 7))
            i += 7;
        else if (!strncmp(text + i, "{break}", 7))
        {
            plain[n++] = '\n';
            i += 7;
        }
        else
            plain[n++] = text[i++];
    }
    plain[n] = '\0';
    if (first < 0)
    {
        creshape_error_line(plain, n);
        free(plain);
        return;
    }
    /* a paragraph: words wrapped at the line width */
    char line[CRESHAPE_LINE_WIDTH + 2];
    size_t at = 0;
    int indent = first;
    const char *p = plain;
    for (;;)
    {
        while (*p == ' ')
            p++;
        const size_t word = strcspn(p, " \n");
        if (word)
        {
            if (at && at + 1 + word > CRESHAPE_LINE_WIDTH)
            {
                creshape_error_line(line, at);
                at = 0;
            }
            if (!at)
            {
                at = (size_t)indent < sizeof(line) / 2 ? (size_t)indent : 0;
                memset(line, ' ', at);
                indent = rest;
            }
            else
                line[at++] = ' ';
            const size_t room = sizeof(line) - 1 - at;
            memcpy(line + at, p, word < room ? word : room);
            at += word < room ? word : room;
            p += word;
        }
        if (*p == '\n' || !*p)
        {
            if (at)
                creshape_error_line(line, at);
            at = 0;
            if (!*p)
                break;
            p++;
        }
    }
    free(plain);
}

void creshape_error_text(const char *text)
{
    while (*text)
    {
        const char *end = strchr(text, '\n');
        creshape_error_plain(text, end ? (size_t)(end - text) : strlen(text));
        if (!end)
            break;
        text = end + 1;
    }
}

void creshape_error(const char *fmt, ...)
{
    char text[1024];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(text, sizeof(text), fmt, ap);
    va_end(ap);
    creshape_error_plain(text, strlen(text));
}

void creshape_note(const char *fmt, ...)
{
    char text[1024];
    va_list ap;
    va_start(ap, fmt);
    size_t at = (size_t)snprintf(text, sizeof(text), "{txt}(");
    int n = vsnprintf(text + at, sizeof(text) - at - 4, fmt, ap);
    va_end(ap);
    at += n < 0 ? 0 : (size_t)n < sizeof(text) - at - 4 ? (size_t)n : sizeof(text) - at - 5;
    snprintf(text + at, sizeof(text) - at, ")\n");
    SF_display(text);
}

static int creshape_oom(void)
{
    creshape_error("creshape: out of memory");
    return 920;
}

/* ------------------------------------------------------------------------
 * Levels
 * ------------------------------------------------------------------------ */

int creshape_level_cmp(const creshape_state *st, const creshape_level *a, const creshape_level *b)
{
    if (st->string_j)
        return strcmp(a->text, b->text);
    return (a->num > b->num) - (a->num < b->num);
}

/* qsort has no context argument; the comparison reads the kind from here.
 * Plans run one at a time on the Stata thread. */
static const creshape_state *creshape_sorting;
static const creshape_level *creshape_sort_levels;

static int creshape_cmp_level_ptr(const void *a, const void *b)
{
    return creshape_level_cmp(creshape_sorting, (const creshape_level *)a,
                              (const creshape_level *)b);
}

static int creshape_cmp_level_at(const void *a, const void *b)
{
    size_t x = *(const size_t *)a, y = *(const size_t *)b;
    int c =
        creshape_level_cmp(creshape_sorting, &creshape_sort_levels[x], &creshape_sort_levels[y]);
    return c ? c : (x > y) - (x < y);
}

int creshape_levels_report(creshape_state *st)
{
    free(st->report_levels);
    st->nreport = 0;
    st->report_levels = (size_t *)malloc((st->nlevels ? st->nlevels : 1) * sizeof(size_t));
    if (!st->report_levels)
        return creshape_oom();
    for (size_t l = 0; l < st->nlevels; l++)
        st->report_levels[st->nreport++] = l;
    creshape_sorting = st;
    creshape_sort_levels = st->levels;
    qsort(st->report_levels, st->nreport, sizeof(size_t), creshape_cmp_level_at);
    return 0;
}

static int creshape_push_level(creshape_state *st, size_t *cap, double num, const char *text)
{
    if (st->nlevels == *cap)
    {
        size_t grown_cap = *cap ? *cap * 2 : 16;
        creshape_level *grown = (creshape_level *)realloc(st->levels, grown_cap * sizeof(*grown));
        if (!grown)
            return 920;
        st->levels = grown;
        *cap = grown_cap;
    }
    creshape_level *lv = &st->levels[st->nlevels];
    lv->num = num;
    lv->found = 0;
    lv->text = creshape_strdup(text);
    if (!lv->text)
        return 920;
    st->nlevels++;
    return 0;
}

/* Sort levels[] and drop repeated values. */
static void creshape_levels_unique(creshape_state *st)
{
    creshape_sorting = st;
    qsort(st->levels, st->nlevels, sizeof(creshape_level), creshape_cmp_level_ptr);
    size_t kept = 0;
    for (size_t l = 0; l < st->nlevels; l++)
    {
        if (kept && !creshape_level_cmp(st, &st->levels[kept - 1], &st->levels[l]))
        {
            free(st->levels[l].text);
            continue;
        }
        st->levels[kept++] = st->levels[l];
    }
    st->nlevels = kept;
}

/* ------------------------------------------------------------------------
 * Reading the request
 * ------------------------------------------------------------------------ */

static int creshape_cols_of(creshape_state *st, const char *local, int **cols_out, size_t *n_out)
{
    char *list = creshape_local(local);
    char **words = NULL;
    if (!list)
        return creshape_oom();
    long n = creshape_split(list, &words);
    if (n < 0)
    {
        free(list);
        return creshape_oom();
    }
    int *cols = (int *)malloc((n ? (size_t)n : 1) * sizeof(int));
    if (!cols)
    {
        free(words);
        free(list);
        return creshape_oom();
    }
    int rc = 0;
    for (long k = 0; k < n && !rc; k++)
    {
        cols[k] = creshape_names_find(st, words[k]);
        if (cols[k] < 0)
        {
            creshape_error("variable {bf:%s} not found", words[k]);
            rc = 111;
        }
    }
    free(words);
    free(list);
    if (rc)
    {
        free(cols);
        return rc;
    }
    *cols_out = cols;
    *n_out = (size_t)n;
    return 0;
}

int creshape_plan_read(creshape_state *st, const char *args)
{
    int flag = 0;
    st->string_j = ctools_parse_int_option(args, "string", &flag) == CTOOLS_PARSE_OK && flag;
    flag = 0;
    st->verbose = ctools_parse_int_option(args, "verbose", &flag) == CTOOLS_PARSE_OK && flag;
    flag = 0;
    st->report = ctools_parse_int_option(args, "report", &flag) == CTOOLS_PARSE_OK && flag;
    flag = 0;
    st->xi_listed = ctools_parse_int_option(args, "xilist", &flag) == CTOOLS_PARSE_OK && flag;
    st->maxvar = 0;
    if (ctools_parse_int_option(args, "maxvar", &st->maxvar) != CTOOLS_PARSE_OK || st->maxvar < 1)
        return 198;

    /* the variable names, in dataset order */
    st->nvars = (size_t)SF_nvars();
    st->nobs = (size_t)SF_nobs();
    st->old_nobs = st->nobs;
    st->old_nvars = st->nvars;
    st->names_buf = creshape_local("__creshape_vars");
    if (!st->names_buf)
        return creshape_oom();
    long n = creshape_split(st->names_buf, &st->names);
    if (n < 0)
        return creshape_oom();
    if ((size_t)n != st->nvars)
        return 198;
    int rc = creshape_names_index(st);
    if (rc)
        return creshape_oom();
    st->is_str = (unsigned char *)malloc(st->nvars ? st->nvars : 1);
    if (!st->is_str)
        return creshape_oom();
    for (size_t j = 0; j < st->nvars; j++)
        st->is_str[j] = SF_var_is_string((ST_int)j + 1) ? 1 : 0;

    if ((rc = creshape_cols_of(st, "__creshape_i", &st->icols, &st->ni)))
        return rc;
    if (!st->ni)
        return 198;

    char *j = creshape_local("__creshape_j");
    if (!j)
        return creshape_oom();
    snprintf(st->jname, sizeof(st->jname), "%s", j);
    free(j);
    st->jcol = creshape_names_find(st, st->jname);
    st->atwl = creshape_local("__creshape_atwl");
    if (!st->atwl)
        return creshape_oom();

    /* the stubs, each split at its @ */
    char *xij = creshape_local("__creshape_xij");
    char **words = NULL;
    if (!xij)
        return creshape_oom();
    n = creshape_split(xij, &words);
    if (n <= 0)
    {
        free(words);
        free(xij);
        return n < 0 ? creshape_oom() : 198;
    }
    st->stubs = (creshape_stub *)calloc((size_t)n, sizeof(creshape_stub));
    if (!st->stubs)
    {
        free(words);
        free(xij);
        return creshape_oom();
    }
    st->nstubs = (size_t)n;
    rc = 0;
    for (long k = 0; k < n && !rc; k++)
    {
        creshape_stub *sb = &st->stubs[k];
        const char *at = strchr(words[k], '@');
        sb->spec = creshape_strdup(words[k]);
        sb->suffix = creshape_strdup(at ? at + 1 : "");
        sb->prefix = creshape_strdup(words[k]);
        if (!sb->spec || !sb->suffix || !sb->prefix)
        {
            rc = creshape_oom();
            break;
        }
        if (at)
            sb->prefix[at - words[k]] = '\0';
        sb->long_name = creshape_join(sb->prefix, st->atwl, sb->suffix);
        sb->long_col = -1;
        if (!sb->long_name)
            rc = creshape_oom();
    }
    free(words);
    free(xij);
    if (rc)
        return rc;

    /* j values listed in j() */
    char *jv = creshape_local("__creshape_jv");
    if (!jv)
        return creshape_oom();
    words = NULL;
    n = creshape_split(jv, &words);
    if (n < 0)
    {
        free(jv);
        return creshape_oom();
    }
    size_t cap = 0;
    for (long k = 0; k < n && !rc; k++)
    {
        double num = 0;
        if (!st->string_j)
        {
            char *end = NULL;
            num = strtod(words[k], &end);
            if (end == words[k] || *end)
            {
                creshape_error("{bf:%s} is not a number; specify option {bf:string}", words[k]);
                rc = 198;
                break;
            }
        }
        if (creshape_push_level(st, &cap, num, words[k]))
            rc = creshape_oom();
    }
    st->user_levels = n > 0;
    free(words);
    free(jv);
    return rc;
}

/* ------------------------------------------------------------------------
 * The xi variables: every variable that is not i, j or xij, or those listed
 * (used[] marks the i, j and xij columns, and on return the listed ones too)
 * ------------------------------------------------------------------------ */

static int creshape_plan_xi(creshape_state *st, unsigned char *used)
{
    st->xicols = (int *)malloc((st->nvars ? st->nvars : 1) * sizeof(int));
    if (!st->xicols)
        return creshape_oom();
    st->nxi = 0;
    if (!st->xi_listed)
    {
        /* every other variable, in dataset order */
        for (size_t c = 0; c < st->nvars; c++)
            if (!used[c])
                st->xicols[st->nxi++] = (int)c;
        return 0;
    }
    /* reshape xi: those variables, in the order listed; the rest are dropped */
    int *listed = NULL;
    size_t nlisted = 0;
    int rc = creshape_cols_of(st, "__creshape_xi", &listed, &nlisted);
    if (rc)
        return rc;
    for (size_t k = 0; k < nlisted; k++)
        if (!used[listed[k]])
        {
            used[listed[k]] = 1;
            st->xicols[st->nxi++] = listed[k];
        }
    free(listed);
    return 0;
}

/* ------------------------------------------------------------------------
 * Long: wide data, the j values are in the variable names
 * ------------------------------------------------------------------------ */

/* The text a stub leaves of name between its two halves, or NULL. */
static const char *creshape_middle(const creshape_stub *sb, const char *name, size_t *len)
{
    size_t n = strlen(name), a = strlen(sb->prefix), c = strlen(sb->suffix);
    if (n <= a + c || memcmp(name, sb->prefix, a) || memcmp(name + n - c, sb->suffix, c))
        return NULL;
    *len = n - a - c;
    return name + a;
}

static int creshape_long_levels(creshape_state *st)
{
    size_t cap = 0;
    char *middle = (char *)malloc(CRESHAPE_NAME_BYTES);
    if (!middle)
        return creshape_oom();
    int rc = 0;
    for (size_t c = 0; c < st->nvars && !rc; c++)
    {
        /* a name two stubs match is read by the later one, as reshape does */
        const char *text = NULL;
        size_t len = 0;
        for (size_t k = 0; k < st->nstubs; k++)
        {
            size_t this_len;
            const char *m = creshape_middle(&st->stubs[k], st->names[c], &this_len);
            if (m)
            {
                text = m;
                len = this_len;
            }
        }
        if (!text)
            continue;
        memcpy(middle, text, len);
        middle[len] = '\0';
        if (st->string_j)
            rc = creshape_push_level(st, &cap, 0, middle);
        else
        {
            double value;
            char number[40];
            if (creshape_parse_real(middle, &value))
                continue;
            creshape_format_number(value, number, sizeof(number));
            rc = creshape_push_level(st, &cap, value, number);
        }
    }
    free(middle);
    if (rc)
        return creshape_oom();
    creshape_levels_unique(st);
    return 0;
}

const char creshape_picture[] = "\n"
                                "         long                                wide\n"
                                "        +---------------+                   +------------------+\n"
                                "        | i   j   a   b |                   | i   a1 a2  b1 b2 |\n"
                                "        |---------------| <-- creshape ---> |------------------|\n"
                                "        | 1   1   1   2 |                   | 1   1   3   2  4 |\n"
                                "        | 1   2   3   4 |                   | 2   5   7   6  8 |\n"
                                "        | 2   1   5   6 |                   +------------------+\n"
                                "        | 2   2   7   8 |\n"
                                "        +---------------+\n";

static void creshape_show_levels(const creshape_state *st)
{
    /* (j = 1 2 3), as reshape reports the values it found */
    size_t need = 16;
    for (size_t l = 0; l < st->nlevels; l++)
        need += strlen(st->levels[l].text) + 1;
    char *line = (char *)malloc(need);
    if (!line)
        return;
    size_t at = (size_t)snprintf(line, need, "{txt}(j =");
    for (size_t l = 0; l < st->nlevels; l++)
        at += (size_t)snprintf(line + at, need - at, " %s", st->levels[l].text);
    snprintf(line + at, need - at, ")\n");
    SF_display(line);
    free(line);
}

int creshape_plan_long(creshape_state *st)
{
    int rc;
    if (!st->user_levels)
    {
        if ((rc = creshape_long_levels(st)))
            return rc;
        if (!st->nlevels)
        {
            creshape_error("no xij variables found");
            creshape_error_text(
                "{p 4 4 2}You typed something like {bf:creshape long a b, i(i) j(j)}.{break}"
                "{bf:creshape} looked for existing variables named {bf:a}{it:#} and "
                "{bf:b}{it:#} but could not find any. Remember this picture:{p_end}");
            creshape_error_text(creshape_picture);
            return 111;
        }
        creshape_show_levels(st);
    }
    else
    {
        /* the observations of the result are in order of j */
        creshape_sorting = st;
        qsort(st->levels, st->nlevels, sizeof(creshape_level), creshape_cmp_level_ptr);
    }
    if ((rc = creshape_levels_report(st)))
        return rc;

    unsigned char *used = (unsigned char *)calloc(st->nvars ? st->nvars : 1, 1);
    if (!used)
        return creshape_oom();
    for (size_t k = 0; k < st->ni; k++)
        used[st->icols[k]] = 1;

    /* each stub's variables, one per j value; all strings or all numbers */
    for (size_t k = 0; k < st->nstubs; k++)
    {
        creshape_stub *sb = &st->stubs[k];
        sb->cols = (int *)malloc(st->nlevels * sizeof(int));
        sb->names = (char **)calloc(st->nlevels, sizeof(char *));
        if (!sb->cols || !sb->names)
        {
            free(used);
            return creshape_oom();
        }
        int kind = -1;
        for (size_t l = 0; l < st->nlevels; l++)
        {
            sb->names[l] = creshape_join(sb->prefix, st->levels[l].text, sb->suffix);
            if (!sb->names[l])
            {
                free(used);
                return creshape_oom();
            }
            int col = creshape_names_find(st, sb->names[l]);
            sb->cols[l] = col;
            if (col < 0)
            {
                if (!creshape_name_ok(sb->names[l]))
                {
                    creshape_error("variable {bf:%s} implied name too long", sb->names[l]);
                    free(used);
                    return 198;
                }
                continue;
            }
            if (kind < 0)
                kind = st->is_str[col];
            else if (kind != st->is_str[col])
            {
                creshape_error("variable {bf:%s} type mismatch with other {bf:%s} variables",
                               sb->names[l], sb->spec);
                free(used);
                return 198;
            }
            used[col] = 1;
            st->levels[l].found = 1;
        }
        sb->is_string = kind > 0;
    }
    for (size_t l = 0; l < st->nlevels; l++)
        for (size_t k = 0; k < st->nstubs; k++)
            if (st->stubs[k].cols[l] < 0)
                creshape_note("variable {bf:%s} not found", st->stubs[k].names[l]);

    rc = creshape_plan_xi(st, used);
    free(used);
    if (rc)
        return rc;

    /* the storage type compress would give j */
    if (st->string_j)
    {
        size_t width = 1;
        for (size_t l = 0; l < st->nlevels; l++)
            if (strlen(st->levels[l].text) > width)
                width = strlen(st->levels[l].text);
        if (width > CTOOLS_DTA_CODE_STR_MAX)
        {
            creshape_error("values of {bf:%s} are longer than 2,045 characters", st->jname);
            return 198;
        }
        st->j_type = (uint16_t)width;
    }
    else
    {
        double lo = 0, hi = 0;
        int whole = 1;
        for (size_t l = 0; l < st->nlevels; l++)
        {
            double v = st->levels[l].num;
            if (v != floor(v))
                whole = 0;
            if (!l || v < lo)
                lo = v;
            if (!l || v > hi)
                hi = v;
        }
        st->j_type = !whole                                      ? CTOOLS_DTA_CODE_DOUBLE
                     : lo >= -127 && hi <= 100                   ? CTOOLS_DTA_CODE_BYTE
                     : lo >= -32767 && hi <= 32740               ? CTOOLS_DTA_CODE_INT
                     : lo >= -2147483647.0 && hi <= 2147483620.0 ? CTOOLS_DTA_CODE_LONG
                                                                 : CTOOLS_DTA_CODE_DOUBLE;
    }

    size_t out;
    if (ctools_safe_mul_size(st->nobs, st->nlevels, &out) != 0 || out > CRESHAPE_MAX_OBS)
    {
        creshape_error("creshape: the long data would have more than 2,147,483,647 observations");
        return 920;
    }
    st->out_nobs = out;
    return 0;
}

/* ------------------------------------------------------------------------
 * Wide: long data, the j values are in the observations
 * ------------------------------------------------------------------------ */

int creshape_plan_wide(creshape_state *st)
{
    if (st->jcol < 0)
    {
        creshape_error("variable {bf:%s} not found", st->jname);
        return 111;
    }
    if (st->is_str[st->jcol] != (st->string_j ? 1 : 0))
    {
        if (st->string_j)
            creshape_error("variable {bf:%s} is numeric", st->jname);
        else
            creshape_error("variable {bf:%s} is string; specify option {bf:string}", st->jname);
        return 109;
    }
    unsigned char *used = (unsigned char *)calloc(st->nvars ? st->nvars : 1, 1);
    if (!used)
        return creshape_oom();
    for (size_t k = 0; k < st->ni; k++)
        used[st->icols[k]] = 1;
    used[st->jcol] = 1;
    for (size_t k = 0; k < st->nstubs; k++)
    {
        creshape_stub *sb = &st->stubs[k];
        sb->long_col = creshape_names_find(st, sb->long_name);
        if (sb->long_col < 0)
        {
            creshape_error("variable {bf:%s} not found", sb->long_name);
            free(used);
            return 111;
        }
        sb->is_string = st->is_str[sb->long_col];
        used[sb->long_col] = 1;
    }
    int rc = creshape_plan_xi(st, used);
    free(used);
    return rc;
}

/* ------------------------------------------------------------------------
 * The variables of the result
 * ------------------------------------------------------------------------ */

static uint32_t creshape_name_hash(const char *s)
{
    uint32_t h = 2166136261u;
    for (; *s; s++)
        h = (h ^ (unsigned char)*s) * 16777619u;
    return h;
}

static void creshape_out_copy(const creshape_state *st, creshape_out *o, int col)
{
    memset(o, 0, sizeof(*o));
    o->kind = CRESHAPE_OUT_COPY;
    o->col = col;
    o->is_string = st->is_str[col];
    snprintf(o->name, sizeof(o->name), "%s", st->names[col]);
}

int creshape_plan_outputs(creshape_state *st)
{
    const int wide = st->dir == CRESHAPE_WIDE;
    size_t middle;
    if (wide)
    {
        if (ctools_safe_mul_size(st->nlevels, st->nstubs, &middle) != 0)
            return creshape_oom();
        for (size_t k = 0; k < st->nstubs; k++)
        {
            creshape_stub *sb = &st->stubs[k];
            sb->names = (char **)calloc(st->nlevels ? st->nlevels : 1, sizeof(char *));
            if (!sb->names)
                return creshape_oom();
            for (size_t l = 0; l < st->nlevels; l++)
                if (!(sb->names[l] = creshape_join(sb->prefix, st->levels[l].text, sb->suffix)))
                    return creshape_oom();
        }
    }
    else
        middle = 1 + st->nstubs;
    const size_t total = st->ni + middle + st->nxi;
    if (total > (size_t)st->maxvar)
    {
        /* reshape's message and code for a wide form that does not fit */
        if (wide)
            creshape_error("variable {bf:%s} takes on too many values", st->jname);
        else
            creshape_error("no room to add more variables");
        char text[256];
        snprintf(text, sizeof(text),
                 "{p 4 4 2}The result would have %zu variables, and {bf:maxvar} is %d; see "
                 "{bf:help maxvar}.{p_end}",
                 total, st->maxvar);
        creshape_error_text(text);
        return wide ? 134 : 900;
    }
    st->out = (creshape_out *)calloc(total ? total : 1, sizeof(creshape_out));
    if (!st->out)
        return creshape_oom();
    size_t n = 0;
    for (size_t k = 0; k < st->ni; k++)
        creshape_out_copy(st, &st->out[n++], st->icols[k]);
    if (wide)
    {
        for (size_t l = 0; l < st->nlevels; l++)
            for (size_t k = 0; k < st->nstubs; k++)
            {
                creshape_out *o = &st->out[n++];
                o->kind = CRESHAPE_OUT_STUB;
                o->col = (int)k;
                o->level = (int)l;
                o->is_string = st->stubs[k].is_string;
                if (!creshape_name_ok(st->stubs[k].names[l]))
                {
                    creshape_error("{bf:%s} invalid name", st->stubs[k].names[l]);
                    return 198;
                }
                snprintf(o->name, sizeof(o->name), "%s", st->stubs[k].names[l]);
            }
    }
    else
    {
        creshape_out *o = &st->out[n++];
        o->kind = CRESHAPE_OUT_J;
        o->is_string = st->string_j;
        snprintf(o->name, sizeof(o->name), "%s", st->jname);
        for (size_t k = 0; k < st->nstubs; k++)
        {
            o = &st->out[n++];
            o->kind = CRESHAPE_OUT_STUB;
            o->col = (int)k;
            o->is_string = st->stubs[k].is_string;
            if (!creshape_name_ok(st->stubs[k].long_name))
            {
                creshape_error("{bf:%s} invalid name", st->stubs[k].long_name);
                return 198;
            }
            snprintf(o->name, sizeof(o->name), "%s", st->stubs[k].long_name);
        }
    }
    for (size_t k = 0; k < st->nxi; k++)
        creshape_out_copy(st, &st->out[n++], st->xicols[k]);
    st->out_nvars = n;

    /* no name twice: a new variable may not take the name of one that stays */
    size_t cap = 16;
    while (cap < 2 * n + 1)
        cap *= 2;
    uint32_t *slots = (uint32_t *)calloc(cap, sizeof(uint32_t));
    if (!slots)
        return creshape_oom();
    int rc = 0;
    for (size_t v = 0; v < n && !rc; v++)
    {
        size_t at = creshape_name_hash(st->out[v].name) & (cap - 1);
        while (slots[at])
        {
            if (!strcmp(st->out[slots[at] - 1].name, st->out[v].name))
            {
                creshape_error("variable {bf:%s} already defined", st->out[v].name);
                rc = 110;
                break;
            }
            at = (at + 1) & (cap - 1);
        }
        slots[at] = (uint32_t)v + 1;
    }
    free(slots);
    return rc;
}

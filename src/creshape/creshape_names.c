/*
 * creshape_names.c - variable names and j values as text.
 *
 * reshape ties the wide and the long form together through names: the wide
 * variables of stub inc are inc<j>, or the text before and after @ around
 * <j>. Nothing here calls Stata, so the rules can be tested alone.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "creshape_internal.h"

char *creshape_strdup(const char *text)
{
    size_t n = strlen(text) + 1;
    char *copy = (char *)malloc(n);
    if (copy)
        memcpy(copy, text, n);
    return copy;
}

char *creshape_join(const char *prefix, const char *middle, const char *suffix)
{
    size_t a = strlen(prefix), b = strlen(middle), c = strlen(suffix);
    char *out = (char *)malloc(a + b + c + 1);
    if (!out)
        return NULL;
    memcpy(out, prefix, a);
    memcpy(out + a, middle, b);
    memcpy(out + a + b, suffix, c + 1);
    return out;
}

long creshape_split(char *list, char ***words)
{
    size_t n = 0, cap = 16;
    char **w = (char **)malloc(cap * sizeof(*w));
    if (!w)
        return -1;
    char *p = list;
    for (;;)
    {
        while (*p == ' ')
            p++;
        if (!*p)
            break;
        if (n == cap)
        {
            char **grown = (char **)realloc(w, 2 * cap * sizeof(*w));
            if (!grown)
            {
                free(w);
                return -1;
            }
            w = grown;
            cap *= 2;
        }
        w[n++] = p;
        while (*p && *p != ' ')
            p++;
        if (*p)
            *p++ = '\0';
    }
    *words = w;
    return (long)n;
}

/*
 * reshape reads the text between a stub's two halves with real(). Of what
 * real() accepts, a variable name can hold only digits and an exponent
 * marked e, E, d or D (inc1e3 is j = 1000). Leading zeros are dropped by the
 * conversion, so inc08 is j = 8, and reshape then looks for inc8.
 */
int creshape_parse_real(const char *text, double *value)
{
    char buf[64];
    size_t n = 0, digits = 0;
    const char *p = text;
    while (*p >= '0' && *p <= '9')
    {
        if (n + 1 >= sizeof(buf))
            return -1;
        buf[n++] = *p++;
        digits++;
    }
    if (!digits)
        return -1;
    if (*p == 'e' || *p == 'E' || *p == 'd' || *p == 'D')
    {
        if (n + 1 >= sizeof(buf))
            return -1;
        buf[n++] = 'e';
        p++;
        digits = 0;
        while (*p >= '0' && *p <= '9')
        {
            if (n + 1 >= sizeof(buf))
                return -1;
            buf[n++] = *p++;
            digits++;
        }
        if (!digits)
            return -1;
    }
    if (*p)
        return -1;
    buf[n] = '\0';
    double v = strtod(buf, NULL);
    /* past Stata's largest double the result is missing */
    if (!(v < 0x1p1023))
        return -1;
    *value = v;
    return 0;
}

void creshape_format_number(double value, char *buf, size_t cap)
{
    if (value == 0)
        snprintf(buf, cap, "0");
    else if (value == floor(value) && fabs(value) < 1e16)
        snprintf(buf, cap, "%.0f", value);
    else if (fabs(value) >= 1e16)
        snprintf(buf, cap, "%.11e", value);
    else
    {
        /* fractions print without the leading zero: .5, -.5 */
        char tmp[40];
        snprintf(tmp, sizeof(tmp), "%.16g", value);
        const char *digits = tmp;
        const char *sign = "";
        if (*digits == '-')
        {
            sign = "-";
            digits++;
        }
        if (digits[0] == '0' && digits[1] == '.')
            digits++;
        snprintf(buf, cap, "%s%s", sign, digits);
    }
}

static const char *const creshape_reserved[] = {
    "_all", "_b", "byte", "_coef", "_cons", "double", "float", "if",    "in",   "int", "long",
    "_n",   "_N", "_pi",  "_pred", "_rc",   "_skip",  "strL",  "using", "with", NULL};

int creshape_name_ok(const char *text)
{
    size_t bytes = strlen(text), chars = 0;
    if (!bytes || bytes >= CRESHAPE_NAME_BYTES)
        return 0;
    for (size_t i = 0; i < bytes; i++)
    {
        unsigned char c = (unsigned char)text[i];
        if (c >= 0x80)
        {
            /* a character outside ASCII counts once, at its first byte */
            if ((c & 0xC0) != 0x80)
                chars++;
            continue;
        }
        chars++;
        if (c == '_' || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'))
            continue;
        if (i && c >= '0' && c <= '9')
            continue;
        return 0;
    }
    if (chars > 32)
        return 0;
    for (int i = 0; creshape_reserved[i]; i++)
        if (!strcmp(text, creshape_reserved[i]))
            return 0;
    if (!strncmp(text, "str", 3) && text[3] >= '1' && text[3] <= '9')
    {
        const char *p = text + 3;
        while (*p >= '0' && *p <= '9')
            p++;
        if (!*p && atoi(text + 3) <= 2045)
            return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------------
 * Name table
 * ------------------------------------------------------------------------ */

static uint32_t creshape_hash(const char *s)
{
    uint32_t h = 2166136261u;
    for (; *s; s++)
        h = (h ^ (unsigned char)*s) * 16777619u;
    return h;
}

int creshape_names_index(creshape_state *st)
{
    size_t cap = 16;
    while (cap < 2 * st->nvars + 1)
        cap *= 2;
    st->name_slots = (uint32_t *)calloc(cap, sizeof(uint32_t));
    if (!st->name_slots)
        return 920;
    st->name_cap = cap;
    for (size_t j = 0; j < st->nvars; j++)
    {
        size_t at = creshape_hash(st->names[j]) & (cap - 1);
        while (st->name_slots[at])
            at = (at + 1) & (cap - 1);
        st->name_slots[at] = (uint32_t)j + 1;
    }
    return 0;
}

int creshape_names_find(const creshape_state *st, const char *name)
{
    if (!st->name_slots)
        return -1;
    size_t at = creshape_hash(name) & (st->name_cap - 1);
    while (st->name_slots[at])
    {
        uint32_t col = st->name_slots[at] - 1;
        if (!strcmp(st->names[col], name))
            return (int)col;
        at = (at + 1) & (st->name_cap - 1);
    }
    return -1;
}

/*
    cegen_strings.c
    concat() and ends(): string results

    Both build one new string per observation. A first parallel pass sizes
    every result, a prefix sum places them in one block, and a second
    parallel pass writes them, so no string is allocated on its own.

    concat() formats a numeric source itself only where the text cannot
    differ from Stata's string(): integers of up to seven digits and
    missing values. A source holding anything else is reported back, and
    the ado passes it again as a string that Stata formatted.

    License: MIT
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ctools_alloc.h"
#include "ctools_config.h"
#include "cegen_internal.h"

/* Spacing of Stata's missing codes: .a = . + 2^1011, .b = . + 2*2^1011, ... */
#define CEGEN_MISSING_STEP 0x1p1011

/* Longest text of a number concat() formats: "-9999999" */
#define CEGEN_NUM_TEXT 8

/* 1 when string(v) is the plain digits of v, or a missing code */
static inline int simple_number(double v, double miss)
{
    if (!(v < miss))
        return v >= miss;
    return v == trunc(v) && v > -1e7 && v < 1e7;
}

/* Text of a simple number; returns its length (no terminator written) */
static size_t simple_number_text(double v, double miss, char *buf)
{
    if (!(v < miss))
    {
        const double code = (v - miss) / CEGEN_MISSING_STEP;
        buf[0] = '.';
        if (code >= 1.0 && code <= 26.0 && code == trunc(code))
        {
            buf[1] = (char)('a' + (int)code - 1);
            return 2;
        }
        return 1;
    }
    char digits[CEGEN_NUM_TEXT];
    size_t nd = 0, len = 0;
    long long value = (long long)v;
    unsigned long long magnitude =
        value < 0 ? (unsigned long long)(-value) : (unsigned long long)value;
    do
    {
        digits[nd++] = (char)('0' + (int)(magnitude % 10));
        magnitude /= 10;
    } while (magnitude);
    if (value < 0)
        buf[len++] = '-';
    while (nd)
        buf[len++] = digits[--nd];
    return len;
}

/* Trim blanks at both ends of text[0..len); returns the new start and length */
static inline const char *trim_blanks(const char *text, size_t *len)
{
    size_t n = *len;
    while (n && text[n - 1] == ' ')
        n--;
    while (n && text[0] == ' ')
    {
        text++;
        n--;
    }
    *len = n;
    return text;
}

/* Offsets of the n results in one block, from their sizes (each plus a
 * terminator); returns the block, or NULL when there is no memory */
static char *arena_layout(const uint32_t *size, size_t n, size_t **offset_out)
{
    size_t *offset = (size_t *)ctools_safe_malloc2(n + 1, sizeof(size_t));
    if (!offset)
        return NULL;
    size_t total = 0;
    for (size_t i = 0; i < n; i++)
    {
        offset[i] = total;
        total += (size_t)size[i] + 1;
    }
    offset[n] = total;
    char *arena = total <= (size_t)PTRDIFF_MAX ? (char *)malloc(total ? total : 1) : NULL;
    if (!arena)
    {
        free(offset);
        return NULL;
    }
    *offset_out = offset;
    return arena;
}

/* ---------------------------------------------------------------------------
   concat()
   --------------------------------------------------------------------------- */

static stata_retcode strings_concat(const cegen_job *job, const stata_data *data,
                                    cegen_result *result, cegen_report *report)
{
    const size_t n = data->nobs, k = job->nsrc;
    const double miss = SV_missval;
    const size_t plen = job->punct_len;
    const char *punct = job->punct;

    /* Numeric sources must hold only numbers this file can format */
    char *needfmt = (char *)calloc(k + 1, 1);
    if (!needfmt)
        return STATA_ERR_MEMORY;
    int any_needfmt = 0;
    for (size_t j = 0; j < k; j++)
    {
        needfmt[j] = '0';
        if (data->vars[j].type != STATA_TYPE_DOUBLE)
            continue;
        const double *x = data->vars[j].data.dbl;
        int complex_value = 0;
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n)) reduction(| : complex_value)
        for (size_t i = 0; i < n; i++)
            complex_value |= !simple_number(x[i], miss);
        if (complex_value)
        {
            needfmt[j] = '1';
            any_needfmt = 1;
        }
    }
    if (any_needfmt)
    {
        report->needfmt = needfmt;
        return STATA_OK;
    }
    free(needfmt);

    /* Pass 1: the length of every value before the last separator is cut,
     * which is also the width egen's variable grows to */
    uint32_t *size = (uint32_t *)ctools_safe_malloc2(n ? n : 1, sizeof(uint32_t));
    if (!size)
        return STATA_ERR_MEMORY;
    size_t widest = 0;
    int too_long = 0;
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n)) reduction(max : widest) reduction(| : too_long)
    for (size_t i = 0; i < n; i++)
    {
        size_t len = k * plen;
        char buf[CEGEN_NUM_TEXT];
        for (size_t j = 0; j < k; j++)
        {
            const stata_variable *v = &data->vars[j];
            if (v->type == STATA_TYPE_DOUBLE)
                len += simple_number_text(v->data.dbl[i], miss, buf);
            else if (v->data.str[i])
                len += strlen(v->data.str[i]);
        }
        if (len > UINT32_MAX / 2)
        {
            too_long = 1;
            len = 0;
        }
        size[i] = (uint32_t)len;
        if (len > widest)
            widest = len;
    }
    if (too_long)
    {
        free(size);
        return STATA_ERR_UNSUPPORTED_TYPE;
    }

    size_t *offset = NULL;
    char *arena = arena_layout(size, n, &offset);
    char **out = (char **)ctools_safe_malloc2(n ? n : 1, sizeof(char *));
    if (!arena || !out)
    {
        free(size);
        free(offset);
        free(arena);
        free(out);
        return STATA_ERR_MEMORY;
    }

    /* Pass 2: write, cut the trailing separator, trim blanks at both ends */
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n))
    for (size_t i = 0; i < n; i++)
    {
        char *dest = arena + offset[i], *p = dest;
        for (size_t j = 0; j < k; j++)
        {
            const stata_variable *v = &data->vars[j];
            if (v->type == STATA_TYPE_DOUBLE)
                p += simple_number_text(v->data.dbl[i], miss, p);
            else if (v->data.str[i])
            {
                const size_t len = strlen(v->data.str[i]);
                memcpy(p, v->data.str[i], len);
                p += len;
            }
            if (plen)
            {
                memcpy(p, punct, plen);
                p += plen;
            }
        }
        size_t len = (size_t)(p - dest) - plen;
        const char *start = trim_blanks(dest, &len);
        if (start != dest)
            memmove(dest, start, len);
        dest[len] = '\0';
        out[i] = dest;
    }
    free(size);
    free(offset);
    result->str = out;
    result->arena = arena;
    report->strlen_max = widest ? widest : 1;
    report->ngroups = n;
    return STATA_OK;
}

/* ---------------------------------------------------------------------------
   ends()
   --------------------------------------------------------------------------- */

/* The part of text that ends() keeps: its start, with the length in *len */
static const char *ends_part(const cegen_job *job, const char *text, size_t *len)
{
    const char *punct = job->punct;
    const size_t plen = job->punct_len, total = strlen(text);
    const char *hit = strstr(text, punct);
    if (job->last)
    {
        /* the last separator, found by scanning on after each one */
        const char *found = NULL;
        while (hit)
        {
            found = hit;
            hit = strstr(hit + plen, punct);
        }
        const char *start = found ? found + plen : text;
        *len = total - (size_t)(start - text);
        return start;
    }
    const char *start = text;
    size_t n;
    if (job->tail)
    {
        /* after the first separator; empty when there is none */
        start = hit ? hit + plen : text + total;
        n = total - (size_t)(start - text);
    }
    else
        n = hit ? (size_t)(hit - text) : total;
    if (job->trim)
        start = trim_blanks(start, &n);
    *len = n;
    return start;
}

static stata_retcode strings_ends(const cegen_job *job, const stata_data *data,
                                  cegen_result *result, cegen_report *report)
{
    const size_t n = data->nobs;
    if (job->nsrc != 1 || data->vars[0].type != STATA_TYPE_STRING || !job->punct_len)
        return STATA_ERR_UNSUPPORTED_TYPE;
    char *const *source = data->vars[0].data.str;

    uint32_t *size = (uint32_t *)ctools_safe_malloc2(n ? n : 1, sizeof(uint32_t));
    if (!size)
        return STATA_ERR_MEMORY;
    size_t widest = 0;
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n)) reduction(max : widest)
    for (size_t i = 0; i < n; i++)
    {
        size_t len = 0;
        if (source[i])
            (void)ends_part(job, source[i], &len);
        size[i] = (uint32_t)len;
        if (len > widest)
            widest = len;
    }

    size_t *offset = NULL;
    char *arena = arena_layout(size, n, &offset);
    char **out = (char **)ctools_safe_malloc2(n ? n : 1, sizeof(char *));
    if (!arena || !out)
    {
        free(size);
        free(offset);
        free(arena);
        free(out);
        return STATA_ERR_MEMORY;
    }
#pragma omp parallel for schedule(static) num_threads(cegen_threads_for(n))
    for (size_t i = 0; i < n; i++)
    {
        char *dest = arena + offset[i];
        size_t len = 0;
        if (source[i])
        {
            const char *start = ends_part(job, source[i], &len);
            memcpy(dest, start, len);
        }
        dest[len] = '\0';
        out[i] = dest;
    }
    free(size);
    free(offset);
    result->str = out;
    result->arena = arena;
    report->strlen_max = widest ? widest : 1;
    report->ngroups = n;
    return STATA_OK;
}

stata_retcode cegen_strings_run(const cegen_job *job, stata_data *data, cegen_result *result,
                                cegen_report *report)
{
    if (job->fn == CEGEN_FN_CONCAT)
        return strings_concat(job, data, result, report);
    if (job->fn == CEGEN_FN_ENDS)
        return strings_ends(job, data, result, report);
    return STATA_ERR_INVALID_INPUT;
}

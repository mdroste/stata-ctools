#include "ctools_text.h"
#include "ctools_alloc.h"
#include "ctools_arena.h"
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Append a str# value to its column block, growing it geometrically. The
 * entry holds the value's offset until the column is complete. */
static int ctools_text_pack(ctools_string_arena *pack, char **entry, const char *value, size_t bytes)
{
    size_t need = pack->used + bytes + 1;
    if (need > pack->capacity)
    {
        size_t grown = pack->capacity <= SIZE_MAX / 2 ? pack->capacity * 2 : need;
        if (grown < need)
            grown = need;
        if (grown < 4096)
            grown = 4096;
        char *base = realloc(pack->base, grown);
        if (!base)
            return 920;
        pack->base = base;
        pack->capacity = grown;
    }
    memcpy(pack->base + pack->used, value, bytes + 1);
    *entry = (char *)(uintptr_t)pack->used;
    pack->used = need;
    return 0;
}

ST_retcode ctools_load_text_data(ctools_filtered_data *fd, ST_int nvars, ctools_text_policy policy)
{
    size_t n = 0;
    if (policy.select_flags & (CTOOLS_LOAD_SKIP_IF | CTOOLS_LOAD_MARKER))
    {
        /* Stata answers an if one observation at a time; the shared
         * selection reads a sample marker in parallel instead. */
        stata_retcode select_rc = ctools_data_select(fd, 0, 0, policy.select_flags);
        if (select_rc != STATA_OK)
            return ctools_stata_rc(select_rc);
        n = fd->data.nobs;
    }
    else
    {
        size_t range = SF_in2() >= SF_in1() ? (size_t)(SF_in2() - SF_in1()) + 1 : 0;
        for (int i = SF_in1(); i <= SF_in2(); i++)
            if (SF_ifobs(i))
                n++;
        fd->data.nobs = n;
        fd->n_range = range;
        fd->was_filtered = n != range;
        fd->obs_map = ctools_safe_cacheline_alloc2(n ? n : 1, sizeof(perm_idx_t));
        if (!fd->obs_map)
            return 920;
        size_t at = 0;
        for (int i = SF_in1(); i <= SF_in2(); i++)
            if (SF_ifobs(i))
                fd->obs_map[at++] = (perm_idx_t)i;
    }
    fd->data.nvars = (size_t)nvars;
    fd->data.vars = ctools_safe_cacheline_alloc2((size_t)nvars, sizeof(stata_variable));
    if (fd->data.vars)
        memset(fd->data.vars, 0, (size_t)nvars * sizeof(stata_variable));
    if (!fd->data.vars)
        return 920;
    char cell[2046]; /* str# values are at most 2045 bytes */
    for (int j = 0; j < nvars; j++)
    {
        stata_variable *v = &fd->data.vars[j];
        v->nobs = n;
        v->type = SF_var_is_string(j + 1) ? STATA_TYPE_STRING : STATA_TYPE_DOUBLE;
        int islong = 0;
        ctools_string_arena *pack = NULL;
        if (v->type == STATA_TYPE_STRING)
        {
            v->data.str = ctools_safe_cacheline_alloc2(n ? n : 1, sizeof(char *));
            if (!v->data.str)
                return 920;
            memset(v->data.str, 0, (n ? n : 1) * sizeof(char *));
            islong = SF_var_is_strl(j + 1);
            if (!islong)
            {
                /* Pack str# values into one block sized by the bytes read; the
                 * arena wrapper frees it in bulk, and stata_variable_free never
                 * reads the entries of a fallback-free arena. */
                pack = calloc(1, sizeof(*pack));
                if (!pack)
                    return 920;
                pack->mode = CTOOLS_STRING_ARENA_NO_FALLBACK;
                v->_arena = pack;
            }
        }
        else
        {
            v->data.dbl = ctools_safe_cacheline_alloc2(n ? n : 1, sizeof(double));
            if (!v->data.dbl)
                return 920;
        }
        for (size_t i = 0; i < n; i++)
        {
            int obs = (int)fd->obs_map[i];
            if ((i & 4095) == 0 && SF_poll())
                return 1;
            if (v->type == STATA_TYPE_DOUBLE)
            {
                int rc = (_stata_)->safevdata(j + 1, obs, &v->data.dbl[i]);
                if (rc)
                    return rc;
                continue;
            }
            const char *value = cell;
            if (islong)
            {
                int length = SF_sdatalen(j + 1, obs);
                if (length < 0 || length == INT_MAX || SF_var_is_binary(j + 1, obs))
                    return 109;
                v->data.str[i] = malloc((size_t)length + 1);
                if (!v->data.str[i])
                    return 920;
                if (SF_strldata(j + 1, obs, v->data.str[i], length + 1) != length)
                    return 609;
                v->data.str[i][length] = 0;
                value = v->data.str[i];
            }
            else
            {
                int rc = SF_sdata(j + 1, obs, cell);
                if (rc)
                    return rc;
            }
            size_t bytes = strlen(value), units = 0;
            /* UTF-8 continuation bytes do not start characters; supplementary
             * characters occupy two units in Excel's UTF-16 cell limit. */
            for (size_t k = 0; k < bytes; k++)
            {
                unsigned char ch = (unsigned char)value[k];
                if ((ch & 0xc0) != 0x80)
                    units += ch >= 0xf0 ? 2 : 1;
            }
            if ((policy.max_bytes && bytes > policy.max_bytes) ||
                (policy.max_utf16_units && units > policy.max_utf16_units))
                return 109;
            if (bytes > v->str_maxlen)
                v->str_maxlen = bytes;
            if (pack && ctools_text_pack(pack, &v->data.str[i], cell, bytes))
                return 920;
        }
        if (pack)
        {
            /* Trim the block to the bytes used, then turn offsets into pointers. */
            if (pack->used && pack->used < pack->capacity)
            {
                char *base = realloc(pack->base, pack->used);
                if (base)
                {
                    pack->base = base;
                    pack->capacity = pack->used;
                }
            }
            for (size_t i = 0; i < n; i++)
                v->data.str[i] = pack->base + (uintptr_t)v->data.str[i];
        }
    }
    return 0;
}

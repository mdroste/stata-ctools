/*
    cjoinby_table.c
    Staged observations for cjoinby

    cjoinby moves whole observations: every output row is a master
    observation beside a using observation. Both sides are therefore staged
    row-major whenever they can be (ctools_rowstage), so that reading the
    data touches each Stata record once and writing gathers one staged row
    per output row. Data the row stage refuses (a str2045 variable, or string
    columns past its 2 GiB bound) are loaded as columns instead. The writer
    sees one layout for both: cj_table_gather copies a row into a tile slot.
*/

#include <stdlib.h>
#include <string.h>

#include "ctools_transport_internal.h"
#include "ctools_release.h"
#include "cjoinby_internal.h"

void cj_table_init(cj_table *t)
{
    if (!t)
        return;
    memset(t, 0, sizeof(*t));
    ctools_rowstage_init(&t->rs);
    ctools_filtered_data_init(&t->cols);
}

static void cj_table_drop_layout(cj_table *t)
{
    if (t->owns_layout)
    {
        free(t->var_idx);
        free(t->width);
        free(t->offset);
    }
    t->var_idx = NULL;
    t->width = NULL;
    t->offset = NULL;
    t->owns_layout = 0;
}

void cj_table_free(cj_table *t)
{
    if (!t)
        return;
    cj_table_drop_layout(t);
    ctools_rowstage_free(&t->rs);
    ctools_filtered_data_free(&t->cols);
    cj_table_init(t);
}

void cj_table_release(cj_table *t)
{
    if (!t)
        return;
    cj_table_drop_layout(t);
    ctools_rowstage_release(&t->rs);
    ctools_filtered_data_release(&t->cols);
    cj_table_init(t);
}

/* The slot layout of the listed variables, as the row stage lays its rows
 * out: numeric columns first, so that every double is 8-byte aligned. Used
 * when the observations are held as columns, or when there are none. */
static stata_retcode cj_table_plan_layout(cj_table *t, const int *var_idx, size_t nvars)
{
    int *widths = ctools_io_ctools_read_strw_from_stata((size_t)SF_nvars());
    stata_retcode rc = STATA_OK;

    t->var_idx = (int *)ctools_safe_malloc2(nvars ? nvars : 1, sizeof(int));
    t->width = (int *)ctools_safe_malloc2(nvars ? nvars : 1, sizeof(int));
    t->offset = (size_t *)ctools_safe_malloc2(nvars ? nvars : 1, sizeof(size_t));
    t->owns_layout = 1;
    if (!t->var_idx || !t->width || !t->offset)
        rc = STATA_ERR_MEMORY;

    size_t stride = 0;
    for (int strings = 0; strings <= 1 && !rc; strings++)
    {
        for (size_t j = 0; j < nvars && !rc; j++)
        {
            int idx = var_idx[j];
            if (idx < 1 || idx > SF_nvars())
            {
                rc = STATA_ERR_INVALID_INPUT;
                break;
            }
            int kind = ctools_io_variable_kind(idx);
            if ((kind != IO_NUMERIC) != strings)
                continue;
            size_t bytes = sizeof(double);
            t->width[j] = 0;
            if (strings)
            {
                int width = widths ? widths[idx - 1] : 0;
                if (kind != IO_FIXED_STRING || width < 1 || width > STATA_STR_MAXLEN)
                {
                    rc = STATA_ERR_UNSUPPORTED_TYPE;
                    break;
                }
                t->width[j] = width;
                bytes = (size_t)width + 1;
            }
            t->var_idx[j] = idx;
            t->offset[j] = stride;
            stride += bytes;
        }
    }
    free(widths);
    if (rc)
    {
        cj_table_drop_layout(t);
        return rc;
    }
    t->stride = (stride + 7) & ~(size_t)7;
    if (t->stride == 0)
        t->stride = 8;
    return STATA_OK;
}

stata_retcode cj_table_load(cj_table *t, const int *var_idx, size_t nvars, size_t obs1, size_t nobs)
{
    if (!t || !var_idx || nvars == 0)
        return STATA_ERR_INVALID_INPUT;
    cj_table_init(t);
    t->nvars = nvars;

    if (nobs == 0)
        return cj_table_plan_layout(t, var_idx, nvars);

    stata_retcode rc = ctools_rowstage_load(&t->rs, var_idx, nvars, NULL, obs1, nobs);
    if (rc == STATA_OK)
    {
        t->rowmajor = 1;
        t->nobs = nobs;
        t->stride = t->rs.stride;
        t->var_idx = t->rs.var_idx;
        t->width = t->rs.width;
        t->offset = t->rs.offset;
        return STATA_OK;
    }
    if (rc != STATA_ERR_UNSUPPORTED_TYPE)
        return rc;

    /* Strings the row stage does not take: columns, packed by their length */
    rc = cj_table_plan_layout(t, var_idx, nvars);
    if (rc)
        return rc;
    int *indices = (int *)ctools_safe_malloc2(nvars, sizeof(int));
    if (!indices)
    {
        cj_table_free(t);
        return STATA_ERR_MEMORY;
    }
    memcpy(indices, var_idx, nvars * sizeof(int));
    rc = ctools_data_load_ex(&t->cols, indices, nvars, obs1, obs1 + nobs - 1,
                             CTOOLS_LOAD_SKIP_IF | CTOOLS_LOAD_NO_SORT_ORDER, NULL);
    free(indices);
    if (rc == STATA_OK && (t->cols.data.nobs != nobs || t->cols.data.nvars != nvars))
        rc = STATA_ERR_STATA_READ;
    if (rc)
    {
        cj_table_free(t);
        return rc;
    }
    t->nobs = nobs;
    return STATA_OK;
}

stata_retcode cj_table_restore(const cj_table *t, size_t obs1)
{
    if (!t)
        return STATA_ERR_INVALID_INPUT;
    if (t->nobs == 0)
        return STATA_OK;
    if (t->rowmajor)
        return ctools_rowstage_store(&t->rs, NULL, t->nobs, obs1);
    /* the store reads the columns only */
    return ctools_data_store_selective((stata_data *)&t->cols.data, t->var_idx, t->nvars, obs1);
}

stata_retcode cj_table_column(const cj_table *t, size_t col, stata_variable *out)
{
    if (!t || !out || col >= t->nvars)
        return STATA_ERR_INVALID_INPUT;
    memset(out, 0, sizeof(*out));
    if (t->nobs == 0)
    {
        out->type = t->width[col] ? STATA_TYPE_STRING : STATA_TYPE_DOUBLE;
        out->buffer_storage = CTOOLS_BUFFER_BORROWED;
        out->string_storage = CTOOLS_STRINGS_BORROWED;
        return STATA_OK;
    }
    if (t->rowmajor)
        return ctools_rowstage_column(&t->rs, col, out);
    *out = t->cols.data.vars[col];
    out->_arena = NULL;
    out->buffer_storage = CTOOLS_BUFFER_BORROWED;
    out->string_storage = CTOOLS_STRINGS_BORROWED;
    return STATA_OK;
}

void cj_table_gather_columns(const cj_table *t, size_t row, char *slot)
{
    const stata_variable *vars = t->cols.data.vars;
    for (size_t j = 0; j < t->nvars; j++)
    {
        char *dst = slot + t->offset[j];
        if (t->width[j] == 0)
        {
            memcpy(dst, &vars[j].data.dbl[row], sizeof(double));
        }
        else
        {
            const char *s = vars[j].data.str[row];
            size_t len = s ? strnlen(s, (size_t)t->width[j]) : 0;
            if (len)
                memcpy(dst, s, len);
            dst[len] = '\0';
        }
    }
}

void cj_table_blank(const cj_table *t, char *slot)
{
    const double missing = SV_missval;
    memset(slot, 0, t->stride);
    for (size_t j = 0; j < t->nvars; j++)
        if (t->width[j] == 0)
            memcpy(slot + t->offset[j], &missing, sizeof(double));
}

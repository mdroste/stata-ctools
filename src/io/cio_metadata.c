#include "cio_internal.h"
#include "ctools_cancel.h"
#include <stdatomic.h>

const char *cio_storage(cio_dataset *dataset, cio_column *c, char *out, size_t len)
{
    if (c->binary || (c->is_string && c->width > 2045))
        return "strL";
    if (c->is_string)
    {
        snprintf(out, len, "str%d", c->width);
        return out;
    }
    if (c->force_double)
        return "double";
    if (!strcmp(dataset->kind, "sasxport5"))
        return "double";
    /* As native, non-integers are double, even when float holds them. */
    int integral = 1;
    double low = 0, high = 0;
    for (size_t i = 0; i < dataset->nobs; i++)
    {
        double x = c->numbers[i];
        if (SF_is_missing(x))
            continue;
        if (floor(x) != x)
            integral = 0;
        if (x < low)
            low = x;
        if (x > high)
            high = x;
    }
    if (integral && low >= -127 && high <= 100)
        return "byte";
    if (integral && low >= -32767 && high <= 32740)
        return "int";
    if (integral && low >= -2147483647.0 && high <= 2147483620.0)
        return "long";
    return "double";
}

/* Every column's storage type, side by side, the first time one is asked
 * for: the ado asks for each column in turn (an Excel column twice), and a
 * numeric column's type takes a pass over its values. Nothing changes a
 * scanned column, and a new scan makes new columns. */
static void cio_storage_all(cio_dataset *dataset)
{
    const int split = dataset->nvar > 1 && dataset->nobs >= 4096;
#pragma omp parallel for schedule(dynamic, 1) if (split)
    for (int j = 0; j < dataset->nvar; j++)
    {
        cio_column *c = &dataset->cols[j];
        char s[32];
        if (!c->storage[0])
            cio_copy_text(c->storage, sizeof(c->storage), cio_storage(dataset, c, s, sizeof(s)));
    }
}

int cio_column_info(cio_dataset *dataset, int index)
{
    if (index < 1 || index > dataset->nvar)
        return 198;
    cio_column *c = &dataset->cols[index - 1];
    int rc;
    if (!c->storage[0])
        cio_storage_all(dataset);
    if ((rc = cio_save_number("__cio_numerictext", c->format_numbers)))
        return rc;
    if ((rc = cio_save_number("__cio_external", c->binary || (c->is_string && c->width > 2045))))
        return rc;
    if ((rc = cio_save("__cio_excel_column", c->excel_column)))
        return rc;
    if ((rc = cio_save("__cio_name", c->name)))
        return rc;
    if ((rc = cio_save("__cio_type", c->storage)))
        return rc;
    if ((rc = cio_save("__cio_label", c->label)))
        return rc;
    char number[40] = "";
    if (c->label_numeric)
        snprintf(number, sizeof(number), "%.17g", c->label_number);
    if ((rc = cio_save("__cio_labelnumber", number)))
        return rc;
    if ((rc = cio_save("__cio_format", c->format)))
        return rc;
    return cio_save("__cio_labelset", c->labelset);
}

int cio_label_info(cio_dataset *dataset, int index)
{
    if (index < 1 || (size_t)index > dataset->nlabels)
        return 198;
    cio_label *l = &dataset->labels[index - 1];
    char s[40];
    int rc;
    if (SF_is_missing(l->value))
    {
        uint64_t bits, base;
        double v = SV_missval;
        memcpy(&bits, &l->value, 8);
        memcpy(&base, &v, 8);
        int n = (int)((bits - base) >> 40);
        snprintf(s, sizeof(s), n ? ".%c" : ".", n ? 'a' + n - 1 : 0);
    }
    else
        snprintf(s, sizeof(s), "%.17g", l->value);
    if ((rc = cio_save("__cio_labelvalue", s)))
        return rc;
    if ((rc = cio_save("__cio_labeltext", l->text)))
        return rc;
    return cio_save("__cio_labelset", l->set);
}

/* Whether cio_load() stores column c (strL and binary columns are stored
 * by the ado) */
static int cio_load_column(const cio_column *c)
{
    return !(c->binary || (c->is_string && c->width > 2045));
}

/* Rows per unit of the parallel store. A unit stores its rows of every
 * column in turn, so the rows' records in Stata's observation-major data
 * stay in cache while each column is written. */
#define CIO_STORE_ROWS 8192

int cio_load(cio_dataset *dataset)
{
    if (SF_nvars() != dataset->nvar || SF_nobs() < (int)dataset->nobs)
        return 198;
    for (int j = 0; j < dataset->nvar; j++)
    {
        const cio_column *c = &dataset->cols[j];
        if (cio_load_column(c) &&
            (!!SF_var_is_string(j + 1) != c->is_string || SF_var_is_strl(j + 1)))
            return 109;
    }
    const size_t nobs = dataset->nobs, units = (nobs + CIO_STORE_ROWS - 1) / CIO_STORE_ROWS;
    atomic_int failed = 0;
    ctools_cancel_team team = {0};
#pragma omp parallel if (units > 1)
    {
#pragma omp for schedule(dynamic, 1) nowait
        for (size_t u = 0; u < units; u++)
        {
            if (atomic_load(&failed) || ctools_cancel_requested())
                continue;
            size_t i0 = u * CIO_STORE_ROWS;
            size_t i1 = i0 + CIO_STORE_ROWS < nobs ? i0 + CIO_STORE_ROWS : nobs;
            for (int j = 0; j < dataset->nvar; j++)
            {
                const cio_column *c = &dataset->cols[j];
                if (!cio_load_column(c))
                    continue;
                int rc = 0;
                if (c->is_string)
                    for (size_t i = i0; i < i1 && !rc; i++)
                        rc = SF_sstore(j + 1, (int)i + 1, c->strings[i] ? c->strings[i] : "");
                else
                    for (size_t i = i0; i < i1 && !rc; i++)
                        rc = (_stata_)->safestore(j + 1, (int)i + 1, c->numbers[i]);
                if (rc)
                {
                    int none = 0;
                    atomic_compare_exchange_strong(&failed, &none, rc);
                    break;
                }
            }
        }
        ctools_cancel_wait(&team);
    }
    if (ctools_cancel_requested())
        return 1;
    return atomic_load(&failed);
}

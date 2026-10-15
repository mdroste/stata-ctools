#include "cio_internal.h"

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
    int integral = 1, exactfloat = 1;
    double low = 0, high = 0;
    for (size_t i = 0; i < dataset->nobs; i++)
    {
        double x = c->numbers[i];
        if (SF_is_missing(x))
            continue;
        if (floor(x) != x)
            integral = 0;
        if ((double)(float)x != x)
            exactfloat = 0;
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
    return !strcmp(dataset->kind, "dbase") && exactfloat ? "float" : "double";
}

int cio_column_info(cio_dataset *dataset, int index)
{
    if (index < 1 || index > dataset->nvar)
        return 198;
    cio_column *c = &dataset->cols[index - 1];
    char s[32];
    int rc;
    if ((rc = cio_save_number("__cio_numerictext", c->format_numbers)))
        return rc;
    if ((rc = cio_save_number("__cio_external", c->binary || (c->is_string && c->width > 2045))))
        return rc;
    if ((rc = cio_save("__cio_excel_column", c->excel_column)))
        return rc;
    if ((rc = cio_save("__cio_name", c->name)))
        return rc;
    if ((rc = cio_save("__cio_type", (char *)cio_storage(dataset, c, s, sizeof(s)))))
        return rc;
    if ((rc = cio_save("__cio_label", c->label)))
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

int cio_load(cio_dataset *dataset)
{
    if (SF_nvars() != dataset->nvar || SF_nobs() < (int)dataset->nobs)
        return 198;
    for (int j = 0; j < dataset->nvar; j++)
    {
        cio_column *c = &dataset->cols[j];
        if (c->binary || (c->is_string && c->width > 2045))
            continue;
        if (!!SF_var_is_string(j + 1) != c->is_string || SF_var_is_strl(j + 1))
            return 109;
        for (size_t i = 0; i < dataset->nobs; i++)
        {
            if ((i & 4095) == 0 && SF_poll())
                return 1;
            int rc = c->is_string ? SF_sstore(j + 1, (int)i + 1, c->strings[i] ? c->strings[i] : "")
                                  : (_stata_)->safestore(j + 1, (int)i + 1, c->numbers[i]);
            if (rc)
                return rc;
        }
    }
    return 0;
}

#include "cio_internal.h"

void *cio_realloc(void *p, size_t n, size_t width)
{
    if (width && n > SIZE_MAX / width)
        return NULL;
    return realloc(p, n * width);
}

void cio_dataset_free(cio_dataset *d)
{
    for (int j = 0; d->cols && j < d->nvar; j++)
    {
        cio_column *c = &d->cols[j];
        if (c->strings)
            for (size_t i = 0; i < d->nobs; i++)
                free(c->strings[i]);
        free(c->strings);
        free(c->numbers);
        free(c->lengths);
    }
    for (size_t i = 0; i < d->nlabels; i++)
    {
        free(d->labels[i].set);
        free(d->labels[i].text);
    }
    free(d->cols);
    free(d->labels);
    memset(d, 0, sizeof(*d));
}

void cio_copy_text(char *out, size_t cap, const char *s)
{
    snprintf(out, cap, "%s", s ? s : "");
}

int cio_local_text(const char *key, char *out, int size)
{
    char name[128];
    snprintf(name, sizeof(name), "_%s", key);
    out[0] = 0;
    return SF_macro_use(name, out, size - 1);
}

int cio_save(const char *key, const char *value)
{
    char name[128];
    if (!strncmp(key, "__cio_", 6))
        snprintf(name, sizeof(name), "CTOOLS_CIO_%s", key + 6);
    else
        cio_copy_text(name, sizeof(name), key);
    return SF_macro_save(name, (char *)(value ? value : ""));
}

int cio_save_number(const char *key, size_t v)
{
    char s[40];
    snprintf(s, sizeof(s), "%zu", v);
    return cio_save((char *)key, s);
}

int cio_local_number(const char *key, size_t v)
{
    char name[128], s[40];
    snprintf(name, sizeof(name), "_%s", key);
    snprintf(s, sizeof(s), "%zu", v);
    return SF_macro_save(name, s);
}

int cio_export_rows(cio_rows *rows, const char *selection)
{
    memset(rows, 0, sizeof(*rows));
    if (!strcmp(selection, "touse"))
    {
        /* Stata answers an if one observation at a time; the shared
           selection reads the sample marker in parallel instead. */
        rows->marker = 1;
        stata_retcode rc = ctools_data_select(&rows->marked, 0, 0, CTOOLS_LOAD_MARKER);
        if (rc != STATA_OK)
            return ctools_stata_rc(rc);
        rows->map = rows->marked.obs_map;
        rows->n = rows->marked.data.nobs;
        return 0;
    }
    rows->first = SF_in1();
    if (rows->first < 1 || SF_in2() < rows->first)
        return 0;
    rows->n = (size_t)(SF_in2() - rows->first) + 1;
    if (!strcmp(selection, "skipif"))
        return 0;
    rows->checked = malloc(rows->n * sizeof(perm_idx_t));
    if (!rows->checked)
        return 909;
    size_t n = 0;
    for (int i = rows->first; i <= SF_in2(); i++)
        if (SF_ifobs(i))
            rows->checked[n++] = (perm_idx_t)i;
    rows->map = rows->checked;
    rows->n = n;
    return 0;
}

void cio_rows_free(cio_rows *rows)
{
    if (rows->marker)
        ctools_filtered_data_free(&rows->marked);
    free(rows->checked);
    memset(rows, 0, sizeof(*rows));
}

double cio_tagged_missing(char tag)
{
    int n = tag >= 'A' && tag <= 'Z' ? tag - 'A' + 1 : tag >= 'a' && tag <= 'z' ? tag - 'a' + 1 : 0;
    uint64_t bits;
    double v = SV_missval;
    memcpy(&bits, &v, 8);
    bits += (uint64_t)n << 40;
    memcpy(&v, &bits, 8);
    return v;
}

double cio_number_value(readstat_value_t v)
{
    if (readstat_value_is_tagged_missing(v))
        return cio_tagged_missing(readstat_value_tag(v));
    if (readstat_value_is_system_missing(v))
        return SV_missval;
    switch (readstat_value_type(v))
    {
    case READSTAT_TYPE_INT8:
        return readstat_int8_value(v);
    case READSTAT_TYPE_INT16:
        return readstat_int16_value(v);
    case READSTAT_TYPE_INT32:
        return readstat_int32_value(v);
    case READSTAT_TYPE_FLOAT:
        return readstat_float_value(v);
    default:
        return readstat_double_value(v);
    }
}

int cio_reserve_rows(cio_dataset *d, size_t rows)
{
    if (rows <= d->capacity)
        return 0;
    if (rows > INT_MAX)
        return 920;
    size_t cap = d->capacity ? d->capacity : 256;
    while (cap < rows)
        cap = cap > INT_MAX / 2 ? (size_t)INT_MAX : cap * 2;
    for (int j = 0; j < d->nvar; j++)
    {
        cio_column *c = &d->cols[j];
        if (c->is_string)
        {
            char **p = cio_realloc(c->strings, cap, sizeof(char *));
            if (!p)
                return 909;
            if (c->is_string)
            {
                size_t *lengths = cio_realloc(c->lengths, cap, sizeof(size_t));
                if (!lengths)
                {
                    c->strings = p;
                    return 909;
                }
                c->lengths = lengths;
                memset(lengths + d->capacity, 0, (cap - d->capacity) * sizeof(size_t));
            }
            c->strings = p;
            memset(p + d->capacity, 0, (cap - d->capacity) * sizeof(char *));
        }
        if (!c->is_string || c->format_numbers)
        {
            double *p = cio_realloc(c->numbers, cap, sizeof(double));
            if (!p)
                return 909;
            c->numbers = p;
            for (size_t i = d->capacity; i < cap; i++)
                p[i] = SV_missval;
        }
    }
    d->capacity = cap;
    return 0;
}

int cio_civil_days(int y, unsigned m, unsigned d)
{
    y -= m <= 2;
    int era = (y >= 0 ? y : y - 399) / 400;
    unsigned yy = (unsigned)(y - era * 400);
    unsigned doy = (153 * (m > 2 ? m - 3 : m + 9) + 2) / 5 + d - 1;
    return era * 146097 + (int)(yy * 365 + yy / 4 - yy / 100 + doy) - 715815;
}

void cio_days_civil(int days, int *y, unsigned *m, unsigned *d)
{
    int z = days + 715815, era = (z >= 0 ? z : z - 146096) / 146097;
    unsigned doe = (unsigned)(z - era * 146097),
             yy = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    *y = (int)yy + era * 400;
    unsigned doy = doe - (365 * yy + yy / 4 - yy / 100), mp = (5 * doy + 2) / 153;
    *d = doy - (153 * mp + 2) / 5 + 1;
    *m = mp < 10 ? mp + 3 : mp - 9;
    *y += *m <= 2;
}

void cio_excel_date_text(double days, uint8_t code, char *out, size_t cap)
{
    if (!isfinite(days) || days < (double)INT_MIN || days > (double)INT_MAX)
    {
        cio_copy_text(out, cap, ".");
        return;
    }
    double ms = days * 86400 * 1000;
    long long whole = 0, seconds = 0;
    if (cio_excel_date_kind(code) == 2)
    {
        long long t = (long long)floor(ms);
        long long sec = t / 1000 - (t % 1000 < 0);
        whole = sec / 86400 - (sec % 86400 < 0);
        seconds = sec - whole * 86400;
    }
    else
        whole = (long long)floor(ms / 86400000);
    int year;
    unsigned month, day;
    cio_days_civil((int)whole, &year, &month, &day);
    static const char *months[] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
    static const char *small[] = {"jan", "feb", "mar", "apr", "may", "jun",
                                  "jul", "aug", "sep", "oct", "nov", "dec"};
    int total = (int)seconds;
    int h = total / 3600, m = total / 60 % 60, s = total % 60;
    char text[100];
    switch (code)
    {
    case 14:
        snprintf(text, sizeof(text), "%u/%u/%04d", month, day, year);
        snprintf(out, cap, "%10s", text);
        break;
    case 15:
        snprintf(out, cap, "%2u-%s-%02d", day, months[month - 1], (year % 100 + 100) % 100);
        break;
    case 16:
        snprintf(out, cap, "%2u-%s", day, months[month - 1]);
        break;
    case 17:
        snprintf(out, cap, "%s-%02d", months[month - 1], (year % 100 + 100) % 100);
        break;
    case 18:
        snprintf(out, cap, "%2d:%02d %s", h % 12 ? h % 12 : 12, m, h < 12 ? "AM" : "PM");
        break;
    case 19:
        snprintf(out, cap, "%2d:%02d:%02d %s", h % 12 ? h % 12 : 12, m, s, h < 12 ? "AM" : "PM");
        break;
    case 20:
        snprintf(out, cap, "%2d:%02d", h, m);
        break;
    case 21:
    case 46:
        snprintf(out, cap, "%02d:%02d:%02d", h, m, s);
        break;
    case 22:
        snprintf(text, sizeof(text), "%u/%u/%04d %d:%02d", month, day, year, h % 12 ? h % 12 : 12,
                 m);
        snprintf(out, cap, "%16s", text);
        break;
    case 45:
    case 47:
        snprintf(out, cap, "%02d:%02d", m, s);
        break;
    case CIO_EXCEL_DATE:
        snprintf(out, cap, "%02u%s%04d", day, small[month - 1], year);
        break;
    default:
        snprintf(out, cap, "%02u%s%04d %02d:%02d:%02d", day, small[month - 1], year, h, m, s);
        break;
    }
}

int cio_export_meta(int index, const char *field, char *out, int cap)
{
    char key[100];
    snprintf(key, sizeof(key), "__cio_%s_%d", field, index);
    return cio_local_text(key, out, cap);
}

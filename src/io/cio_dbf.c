#include "ctools_file.h"
#include "cio_internal.h"

unsigned cio_dbf_u16(const unsigned char *p)
{
    return (unsigned)p[0] | (unsigned)p[1] << 8;
}

uint32_t cio_dbf_u32(const unsigned char *p)
{
    return (uint32_t)p[0] | (uint32_t)p[1] << 8 | (uint32_t)p[2] << 16 | (uint32_t)p[3] << 24;
}

void cio_dbf_p16(unsigned char *p, unsigned v)
{
    p[0] = (unsigned char)v;
    p[1] = (unsigned char)(v >> 8);
}

void cio_dbf_p32(unsigned char *p, uint32_t v)
{
    for (int i = 0; i < 4; i++)
        p[i] = (unsigned char)(v >> (8 * i));
}

int cio_scan_dbf(cio_dataset *dataset, const char *path)
{
    FILE *fp = ctools_fopen(path, "rb");
    if (!fp)
        return 601;
    unsigned char h[32];
    int rc = 610;
    unsigned char *fields = NULL, *row = NULL;
    if (fread(h, 1, 32, fp) != 32 || (h[0] != 3 && h[0] != 4 && h[0] != 0x83 && h[0] != 0x8b))
        goto done;
    unsigned header = cio_dbf_u16(h + 8), record = cio_dbf_u16(h + 10);
    uint32_t n = cio_dbf_u32(h + 4);
    if (header < 33 || record < 1 || (header - 33) % 32 || n > INT_MAX)
        goto done;
    const char *datekeys[] = {"version", "year", "month", "day"};
    for (int j = 0; j < 4; j++)
    {
        char key[80];
        snprintf(key, sizeof(key), "__cio_dbf_%s", datekeys[j]);
        if ((rc = cio_save_number(key, h[j])))
            goto done;
    }
    rc = 610;
    dataset->nvar = (int)((header - 33) / 32);
    if (dataset->nvar > 255)
        goto done;
    dataset->cols = calloc((size_t)dataset->nvar, sizeof(cio_column));
    fields = malloc((size_t)dataset->nvar * 32);
    row = malloc(record);
    if ((!dataset->cols && dataset->nvar) || !fields || !row)
    {
        rc = 909;
        goto done;
    }
    if (fread(fields, 32, (size_t)dataset->nvar, fp) != (size_t)dataset->nvar || fgetc(fp) != 13)
        goto done;
    size_t offset = 1;
    for (int j = 0; j < dataset->nvar; j++)
    {
        unsigned char *f = fields + 32 * j;
        cio_column *c = &dataset->cols[j];
        memcpy(c->name, f, 11);
        c->name[11] = 0;
        cio_copy_text(c->label, sizeof(c->label), c->name);
        if (!strchr("CNFLD", f[11]) || !f[16] || offset + f[16] > record)
            goto done;
        c->is_string = f[11] == 'C';
        c->width = 1;
        if (c->is_string)
            cio_copy_text(c->format, sizeof(c->format), "%9s");
        else if (f[11] == 'D')
            cio_copy_text(c->format, sizeof(c->format), "%td");
        else
            snprintf(c->format, sizeof(c->format), "%%%u.%uf", f[16], f[17]);
        offset += f[16];
    }
    for (uint32_t i = 0; i < n; i++)
    {
        if ((i & 4095) == 0 && SF_poll())
        {
            rc = 1;
            goto done;
        }
        if (fread(row, 1, record, fp) != record)
            goto done;
        if (row[0] == '*')
            continue;
        if (row[0] != ' ')
            goto done;
        rc = cio_reserve_rows(dataset, dataset->nobs + 1);
        if (rc)
            goto done;
        rc = 610;
        size_t obs = dataset->nobs++;
        offset = 1;
        for (int j = 0; j < dataset->nvar; j++)
        {
            unsigned char *f = fields + 32 * j;
            cio_column *c = &dataset->cols[j];
            char value[256];
            memcpy(value, row + offset, f[16]);
            value[f[16]] = 0;
            offset += f[16];
            char *s = value;
            while (*s == ' ')
                s++;
            size_t len = strlen(s);
            while (len && s[len - 1] == ' ')
                s[--len] = 0;
            if (c->is_string)
            {
                c->strings[obs] = strdup(s);
                if (!c->strings[obs])
                {
                    rc = 909;
                    goto done;
                }
                if ((int)len > c->width)
                    c->width = (int)len;
            }
            else
            {
                double x = SV_missval;
                if (f[11] == 'L')
                {
                    if (*s == 'T' || *s == 't' || *s == 'Y' || *s == 'y')
                        x = 1;
                    else if (*s == 'F' || *s == 'f' || *s == 'N' || *s == 'n')
                        x = 0;
                }
                else if (f[11] == 'D' && len == 8)
                {
                    int y, m, d;
                    if (sscanf(s, "%4d%2d%2d", &y, &m, &d) == 3 && m >= 1 && m <= 12 && d >= 1 &&
                        d <= 31)
                        x = cio_civil_days(y, (unsigned)m, (unsigned)d);
                }
                else if (len)
                {
                    if (s[0] == '.' && s[1] >= 'a' && s[1] <= 'z' && !s[2])
                        x = cio_tagged_missing(s[1]);
                    else
                    {
                        char *end;
                        double t = strtod(s, &end);
                        if (end != s && !*end && isfinite(t))
                            x = t;
                    }
                }
                c->numbers[obs] = x;
            }
        }
    }
    for (int j = 0; j < dataset->nvar; j++)
        if (dataset->cols[j].is_string)
            snprintf(dataset->cols[j].format, sizeof(dataset->cols[j].format), "%%%ds",
                     dataset->cols[j].width > 9 ? dataset->cols[j].width : 9);
    rc = 0;
done:
    free(fields);
    free(row);
    fclose(fp);
    return rc;
}

void cio_dbf_number(double x, int width, char *out, size_t cap)
{
    if (x == 0)
    {
        cio_copy_text(out, cap, "0");
        return;
    }
    int digits = width - 2;
    char text[64];
    snprintf(text, sizeof(text), "%.*e", digits - 1, x);
    int exponent = atoi(strchr(text, 'e') + 1),
        mantissa = width - (exponent <= -100 || exponent >= 100 ? 8 : 7);
    if (exponent >= digits || (exponent < 0 && digits + exponent < mantissa))
    {
        snprintf(out, cap, "%.*e", mantissa, x);
        return;
    }
    snprintf(out, cap, "%.*f", exponent >= 0 ? digits - 1 - exponent : digits, x);
    size_t n = strlen(out);
    if (strchr(out, '.'))
    {
        while (out[n - 1] == '0')
            out[--n] = 0;
        if (out[n - 1] == '.')
            out[--n] = 0;
    }
    if (out[0] == '0' && out[1] == '.')
        memmove(out, out + 1, n);
    else if (out[0] == '-' && out[1] == '0' && out[2] == '.')
        memmove(out + 1, out + 2, n - 1);
}

int cio_export_dbf(const char *selection)
{
    char path[4096], replace[8], version[12], datafmt[8];
    cio_local_text("__cio_filename", path, sizeof(path));
    cio_local_text("__cio_replace", replace, sizeof(replace));
    cio_local_text("__cio_version", version, sizeof(version));
    cio_local_text("__cio_datafmt", datafmt, sizeof(datafmt));
    char count[40];
    cio_local_text("__cio_nvar", count, sizeof(count));
    int k = atoi(count), nplugin = cio_export_nvar(selection), rc = 0;
    size_t record = 1;
    if (k < 1)
        return 102;
    if (k > 255)
        return 103;
    unsigned char *fields = calloc((size_t)k, 32), *row = NULL;
    /* Field kinds: 0 general number, 1 date, 2 int/long, 3 byte. */
    int *kinds = calloc((size_t)k, sizeof(int));
    ctools_output output;
    memset(&output, 0, sizeof(output));
    FILE *fp = NULL;
    cio_rows rows;
    memset(&rows, 0, sizeof(rows));
    if (!fields || !kinds)
    {
        rc = 909;
        goto done;
    }
    for (int j = 0; j < k; j++)
    {
        char name[301], type[32], format[257];
        unsigned char *f = fields + 32 * j;
        cio_export_meta(j + 1, "name", name, sizeof(name));
        cio_export_meta(j + 1, "type", type, sizeof(type));
        cio_export_meta(j + 1, "format", format, sizeof(format));
        if (strlen(name) > 10)
        {
            rc = 693;
            goto done;
        }
        memcpy(f, name, strlen(name));
        int width, decimal = 0;
        if (SF_var_is_string(j + 1))
        {
            width = atoi(type + 3);
            f[11] = 'C';
            if (SF_var_is_strl(j + 1))
            {
                rc = 693;
                goto done;
            }
            if (width > 254)
            {
                rc = 108;
                goto done;
            }
        }
        else if (!strncmp(format, "%td", 3) || !strcmp(format, "%d"))
        {
            width = 8;
            f[11] = 'D';
            kinds[j] = 1;
        }
        else if (!strncmp(format, "%t", 2))
        {
            rc = 120;
            goto done;
        }
        else
        {
            f[11] = 'N';
            kinds[j] = !strcmp(type, "byte")                           ? 3
                       : !strcmp(type, "int") || !strcmp(type, "long") ? 2
                                                                       : 0;
            /* Native widths: byte 3, int 5, long 10, float/double 20 columns */
            width = kinds[j] == 3 ? 3 : !strcmp(type, "int") ? 5 : !strcmp(type, "long") ? 10 : 20;
            if (!strcmp(datafmt, "1"))
            {
                kinds[j] = 0;
                sscanf(format, "%%%d.%d", &width, &decimal);
                if (width < 1 || width > 254 || decimal > width)
                {
                    rc = 120;
                    goto done;
                }
            }
        }
        f[16] = (unsigned char)width;
        f[17] = (unsigned char)decimal;
        record += (size_t)width;
    }
    if ((rc = cio_export_rows(&rows, selection)))
        goto done;
    for (size_t r = 0; r < rows.n; r++)
    {
        int i = cio_row(&rows, r);
        /* Native's integer fields truncate the most negative values: byte
           -100..-127 to "-10", int -10000..-32767 to "-1000", long up from
           -1000000000 to "-100000000". Widen a field only when such values are
           exported (by one column: byte 4, int 6, long 11). */
        for (int j = 0; j < k; j++)
            if (kinds[j] >= 2)
            {
                unsigned width = fields[32 * j + 16];
                double bound = 1;
                for (unsigned d = 1; d < width; d++)
                    bound *= 10;
                double x;
                if ((rc = (_stata_)->safevdata(j + 1, i, &x)))
                    goto done;
                if (!SF_is_missing(x) && x <= -bound)
                {
                    char text[40];
                    int need = snprintf(text, sizeof(text), "%.0f", x);
                    if (need > (int)width && need < 255)
                    {
                        fields[32 * j + 16] = (unsigned char)need;
                        record += (size_t)need - width;
                    }
                }
            }
    }
    if (record > 4000)
    {
        rc = 902;
        goto done;
    }
    row = malloc(record);
    if (!row)
    {
        rc = 909;
        goto done;
    }
    rc = ctools_output_prepare(&output, path, !strcmp(replace, "1"));
    if (rc)
        goto done;
    fp = ctools_fopen(output.temporary, "wb");
    if (!fp)
    {
        rc = 603;
        goto done;
    }
    unsigned char header[32] = {0};
    header[0] = !strcmp(version, "IV") ? 4 : 3;
    time_t now = time(NULL);
    struct tm *date = localtime(&now);
    if (date)
    {
        header[1] = (unsigned char)date->tm_year;
        header[2] = (unsigned char)(date->tm_mon + 1);
        header[3] = (unsigned char)date->tm_mday;
    }
    char orig[8];
    cio_local_text("__cio_origdbfdate", orig, sizeof(orig));
    if (!strcmp(orig, "1"))
    {
        const char *keys[] = {"year", "month", "day"};
        for (int j = 0; j < 3; j++)
        {
            char key[80], value[40];
            snprintf(key, sizeof(key), "__cio_dbf_%s", keys[j]);
            cio_local_text(key, value, sizeof(value));
            if (*value)
                header[j + 1] = (unsigned char)atoi(value);
        }
    }
    cio_dbf_p32(header + 4, (uint32_t)rows.n);
    cio_dbf_p16(header + 8, 33 + (unsigned)k * 32);
    cio_dbf_p16(header + 10, (unsigned)record);
    if (fwrite(header, 1, 32, fp) != 32 || fwrite(fields, 32, (size_t)k, fp) != (size_t)k ||
        fputc(13, fp) == EOF)
    {
        rc = 603;
        goto done;
    }
    for (size_t r = 0; r < rows.n; r++)
    {
        int i = cio_row(&rows, r);
        if ((i & 4095) == 0 && SF_poll())
        {
            rc = 1;
            goto done;
        }
        memset(row, ' ', record);
        size_t offset = 1;
        for (int j = 0; j < k; j++)
        {
            char value[2046] = {0};
            unsigned char *f = fields + 32 * j;
            if (SF_var_is_string(j + 1))
            {
                if ((rc = SF_sdata(j + 1, i, value)))
                    goto done;
            }
            else
            {
                double x;
                if ((rc = (_stata_)->safevdata(j + 1, i, &x)))
                    goto done;
                if (SF_is_missing(x))
                {
                    if (kinds[j] != 1)
                    {
                        uint64_t bits, base;
                        double mv = SV_missval;
                        memcpy(&bits, &x, 8);
                        memcpy(&base, &mv, 8);
                        int code = (int)((bits - base) >> 40);
                        snprintf(value, sizeof(value), code ? ".%c" : ".",
                                 code ? 'a' + code - 1 : 0);
                    }
                }
                else if (kinds[j] == 1)
                {
                    int year;
                    unsigned month, day;
                    /* %td shows the floor of a fractional day. The 8-byte D
                       field holds only 0000-01-01 .. 9999-12-31, so other
                       values fail before the int conversion. */
                    double days = floor(x);
                    if (!(days >= -715875.0 && days <= 2936549.0))
                    {
                        rc = 108;
                        goto done;
                    }
                    cio_days_civil((int)days, &year, &month, &day);
                    snprintf(value, sizeof(value), "%04d%02u%02u", year, month, day);
                }
                else if (kinds[j] >= 2)
                {
                    snprintf(value, sizeof(value), "%.0f", x);
                }
                else if (!strcmp(datafmt, "1"))
                {
                    char index[40];
                    cio_export_meta(j + 1, "fmtvar", index, sizeof(index));
                    int var = atoi(index);
                    if (var < k + 1 || var > nplugin)
                    {
                        rc = 198;
                        goto done;
                    }
                    if ((rc = SF_sdata(var, i, value)))
                        goto done;
                }
                else
                    cio_dbf_number(x, f[16], value, sizeof(value));
            }
            size_t len = strlen(value);
            if (len > f[16])
            {
                rc = 108;
                goto done;
            }
            memcpy(row + offset, value, len);
            offset += f[16];
        }
        if (fwrite(row, 1, record, fp) != record)
        {
            rc = 603;
            goto done;
        }
    }
    if (fputc(26, fp) == EOF)
    {
        rc = 603;
        goto done;
    }
    if ((rc = cio_local_number("__cio_N", rows.n)))
        goto done;
    if (fclose(fp))
    {
        fp = NULL;
        rc = 603;
        goto done;
    }
    fp = NULL;
    rc = ctools_output_commit(&output);
done:
    if (fp)
        fclose(fp);
    ctools_output_cleanup(&output);
    free(fields);
    free(row);
    free(kinds);
    cio_rows_free(&rows);
    return rc;
}

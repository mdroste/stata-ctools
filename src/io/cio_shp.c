#include "ctools_file.h"
#include "cio_internal.h"

/* ESRI shapefile geometry. rec_header matches Stata's binary strL layout:
   big-endian content length followed by the record prefix before XY vertices.
 */
uint32_t cio_shp_be32(const unsigned char *p)
{
    return (uint32_t)p[3] | (uint32_t)p[2] << 8 | (uint32_t)p[1] << 16 | (uint32_t)p[0] << 24;
}
void cio_shp_putbe32(unsigned char *p, uint32_t v)
{
    for (int j = 0; j < 4; j++)
        p[j] = (unsigned char)(v >> (8 * (3 - j)));
}
double cio_shp_double(const unsigned char *p)
{
    uint64_t b = 0;
    for (int j = 0; j < 8; j++)
        b |= (uint64_t)p[j] << (8 * j);
    double d;
    memcpy(&d, &b, 8);
    return d;
}
void cio_shp_putdouble(unsigned char *p, double x)
{
    uint64_t b;
    memcpy(&b, &x, 8);
    for (int j = 0; j < 8; j++)
        p[j] = (unsigned char)(b >> (8 * j));
}
static unsigned char shp_file_header[100];
int cio_shp_dims(unsigned type, int *z, int *m, int *point, int *multi)
{
    *z = type == 11 || type == 13 || type == 15 || type == 18 || type == 31;
    *m = *z || type == 21 || type == 23 || type == 25 || type == 28;
    *point = type == 0 || type == 1 || type == 11 || type == 21;
    *multi = type == 8 || type == 18 || type == 28;
    return type == 0 || type == 1 || type == 3 || type == 5 || type == 8 || *z || *m;
}
int cio_shp_add_row(cio_dataset *dataset, int id, double x, double y, double z, double m, int order,
                    const unsigned char *blob, size_t bytes)
{
    int rc = cio_reserve_rows(dataset, dataset->nobs + 1);
    if (rc)
        return rc;
    size_t i = dataset->nobs++;
    dataset->cols[0].numbers[i] = id;
    dataset->cols[1].numbers[i] = x;
    dataset->cols[2].numbers[i] = y;
    if (dataset->shp_z)
        dataset->cols[dataset->shp_z].numbers[i] = z;
    if (dataset->shp_m)
        dataset->cols[dataset->shp_m].numbers[i] = m;
    if (dataset->shp_order)
        dataset->cols[dataset->shp_order].numbers[i] = order;
    cio_column *c = &dataset->cols[dataset->shp_blob];
    if (bytes)
    {
        c->strings[i] = malloc(bytes);
        if (!c->strings[i])
            return 909;
        memcpy(c->strings[i], blob, bytes);
        c->lengths[i] = bytes;
    }
    return 0;
}
int cio_scan_shp(cio_dataset *dataset, const char *path)
{
    FILE *fp = ctools_fopen(path, "rb");
    if (!fp)
        return 601;
    int rc = 610;
    unsigned char *body = NULL, *blob = NULL;
    if (fread(shp_file_header, 1, 100, fp) != 100)
    {
        rc = 3300;
        goto done;
    }
    if (cio_shp_be32(shp_file_header) != 9994 || cio_dbf_u32(shp_file_header + 28) != 1000)
        goto done;
    uint64_t expected_bytes = (uint64_t)cio_shp_be32(shp_file_header + 24) * 2;
    if (expected_bytes < 100 || ctools_fseek(fp, 0, SEEK_END) != 0) goto done;
    int64_t actual_bytes = ctools_ftell(fp);
    if (actual_bytes < 0 || (uint64_t)actual_bytes != expected_bytes ||
        ctools_fseek(fp, 100, SEEK_SET) != 0) {
        rc = 3300; /* native import shp reports a truncated record stream */
        goto done;
    }
    unsigned type = cio_dbf_u32(shp_file_header + 32);
    int hasz, hasm, ispoint, ismulti;
    if (!cio_shp_dims(type, &hasz, &hasm, &ispoint, &ismulti))
    {
        rc = 198;
        goto done;
    }
    dataset->nvar = type ? 4 + hasz + hasm + !ispoint : 0;
    dataset->cols = calloc((size_t)dataset->nvar, sizeof(cio_column));
    if (!dataset->cols)
    {
        rc = 909;
        goto done;
    }
    const char *names[7] = {"_ID", "_X", "_Y"};
    int next = 3;
    if (hasz)
    {
        names[next] = "_Z";
        dataset->shp_z = next++;
    }
    if (hasm)
    {
        names[next] = "_M";
        dataset->shp_m = next++;
    }
    dataset->shp_blob = next;
    names[next++] = "rec_header";
    if (!ispoint)
    {
        dataset->shp_order = next;
        names[next] = "shape_order";
    }
    for (int j = 0; j < dataset->nvar; j++)
    {
        cio_column *c = &dataset->cols[j];
        cio_copy_text(c->name, sizeof(c->name), names[j]);
        cio_copy_text(c->format, sizeof(c->format),
                      j == dataset->shp_blob              ? "%9s"
                      : j == 0 || j == dataset->shp_order ? "%12.0g"
                                                          : "%10.0g");
        c->force_double = j > 0 && j != dataset->shp_blob && j != dataset->shp_order;
    }
    if (type)
        dataset->cols[dataset->shp_blob].is_string = 1;
    if (type)
        dataset->cols[dataset->shp_blob].binary = 1;
    uint32_t count = 0;
    for (;;)
    {
        unsigned char rh[8];
        size_t got = fread(rh, 1, 8, fp);
        if (!got && feof(fp))
            break;
        if (got != 8)
            goto done;
        uint32_t id = cio_shp_be32(rh), words = cio_shp_be32(rh + 4);
        if (id > INT_MAX || words < 2 || words > INT_MAX / 2)
            goto done;
        size_t bytes = (size_t)words * 2;
        int64_t position = ctools_ftell(fp);
        if (position < 0 || (uint64_t)position > expected_bytes ||
            bytes > expected_bytes - (uint64_t)position) goto done;
        body = malloc(bytes);
        if (!body)
        {
            rc = 909;
            goto done;
        }
        if (fread(body, 1, bytes, fp) != bytes)
            goto done;
        uint32_t kind = cio_dbf_u32(body);
        if (kind != type && kind != 0)
            goto done;
        count++;
        if (kind == 0 || ispoint)
        {
            if (bytes < (kind ? 20 + (hasz ? 8 : 0) : 4))
                goto done;
            unsigned char b[8];
            memcpy(b, rh + 4, 4);
            memcpy(b + 4, body, 4);
            double z = SV_missval, m = SV_missval;
            if (kind && hasz)
                z = cio_shp_double(body + 20);
            if (kind && hasm && bytes >= 20 + (size_t)hasz * 8 + 8)
                m = cio_shp_double(body + 20 + (size_t)hasz * 8);
            rc = type
                     ? cio_shp_add_row(dataset, (int)id,
                                       kind ? cio_shp_double(body + 4) : SV_missval,
                                       kind ? cio_shp_double(body + 12) : SV_missval, z, m, 1, b, 8)
                     : 0;
            if (rc)
                goto done;
        }
        else
        {
            if (bytes < (ismulti ? 40 : 44))
                goto done;
            uint32_t parts = ismulti ? 1 : cio_dbf_u32(body + 36),
                     points = cio_dbf_u32(body + (ismulti ? 36 : 40));
            size_t prefix = ismulti ? 40 : 44 + (size_t)parts * 4 * (kind == 31 ? 2 : 1);
            if (!parts || prefix > bytes || (size_t)points > (bytes - prefix) / 16)
                goto done;
            size_t rangez = prefix + (size_t)points * 16,
                   rangem = rangez + (hasz ? 16 + (size_t)points * 8 : 0);
            if (hasz && (rangez > bytes || 16 + (size_t)points * 8 > bytes - rangez))
                goto done;
            /* Stata 16+ omits measures for MultipointM and Multipatch. */
            int presentm = hasm && kind != 28 && kind != 31 && rangem <= bytes &&
                           16 + (size_t)points * 8 <= bytes - rangem;
            size_t headerbytes = prefix + 4 + (hasz ? 16 : 0) + (presentm ? 16 : 0);
            blob = malloc(headerbytes);
            if (!blob)
            {
                rc = 909;
                goto done;
            }
            memcpy(blob, rh + 4, 4);
            memcpy(blob + 4, body, prefix);
            if (hasz)
                memcpy(blob + 4 + prefix, body + rangez, 16);
            if (presentm)
                memcpy(blob + 4 + prefix + (hasz ? 16 : 0), body + rangem, 16);
            int order = 1;
            for (uint32_t part = 0; part < parts; part++)
            {
                uint32_t from = ismulti ? 0 : cio_dbf_u32(body + 44 + (size_t)part * 4),
                         to = part + 1 == parts ? points
                                                : cio_dbf_u32(body + 44 + (size_t)(part + 1) * 4);
                if (from > to || to > points || (part == 0 && from != 0))
                {
                    rc = 610;
                    goto done;
                }
                rc = cio_shp_add_row(dataset, (int)id, SV_missval, SV_missval, SV_missval,
                                     SV_missval, order++, part == 0 ? blob : NULL,
                                     part == 0 ? headerbytes : 0);
                if (rc)
                    goto done;
                for (uint32_t point = from; point < to; point++)
                {
                    const unsigned char *p = body + prefix + (size_t)point * 16;
                    double z = hasz ? cio_shp_double(body + rangez + 16 + (size_t)point * 8)
                                    : SV_missval,
                           m = presentm ? cio_shp_double(body + rangem + 16 + (size_t)point * 8)
                                        : SV_missval;
                    rc = cio_shp_add_row(dataset, (int)id, cio_shp_double(p), cio_shp_double(p + 8),
                                         z, m, order++, NULL, 0);
                    if (rc)
                        goto done;
                }
            }
            free(blob);
            blob = NULL;
        }
        free(body);
        body = NULL;
        rc = 610;
    }
    rc = cio_save_number("__cio_shp_num_records", count);
    if (!rc)
        rc = cio_save_number("__cio_shp_type", type);
    const char *box[] = {"xmin", "ymin", "xmax", "ymax", "zmin", "zmax", "mmin", "mmax"};
    for (int j = 0; !rc && j < 8; j++)
    {
        char key[80], val[64];
        snprintf(key, sizeof(key), "__cio_shp_%s", box[j]);
        snprintf(val, sizeof(val), "%.18g", cio_shp_double(shp_file_header + 36 + 8 * j));
        rc = cio_save(key, val);
    }
done:
    free(body);
    free(blob);
    fclose(fp);
    return rc;
}
int cio_blob_info(cio_dataset *dataset, int col, int obs, size_t offset)
{
    if (col < 1 || col > dataset->nvar || obs < 1 || (size_t)obs > dataset->nobs)
        return 198;
    cio_column *c = &dataset->cols[col - 1];
    if (!c->is_string)
        return 109;
    size_t n =
        c->binary ? c->lengths[obs - 1] : (c->strings[obs - 1] ? strlen(c->strings[obs - 1]) : 0);
    if (offset > n)
        return 198;
    size_t len = n - offset;
    if (len > 16000)
        len = 16000;
    char hex[32001];
    const char *digits = "0123456789abcdef";
    unsigned char *data = (unsigned char *)c->strings[obs - 1];
    for (size_t i = 0; i < len; i++)
    {
        unsigned byte = data[offset + i];
        hex[2 * i] = digits[byte >> 4];
        hex[2 * i + 1] = digits[byte & 15];
    }
    hex[2 * len] = 0;
    int rc = cio_save("__cio_blobhex", hex);
    if (rc)
        return rc;
    return cio_save_number("__cio_blobnext", offset + len < n ? offset + len : 0);
}
int cio_shp_export(void)
{
    char filename[4096], replace[8], shx[8], text[80];
    unsigned char h[100] = {0};
    FILE *fp = NULL, *index = NULL;
    ctools_output output, idx;
    memset(&output, 0, sizeof(output));
    memset(&idx, 0, sizeof(idx));
    int rc = 0;
    cio_local_text("__cio_filename", filename, sizeof(filename));
    cio_local_text("__cio_replace", replace, sizeof(replace));
    cio_local_text("__cio_shx", shx, sizeof(shx));
    cio_shp_putbe32(h, 9994);
    cio_dbf_p32(h + 28, 1000);
    cio_local_text("__cio_shp_type", text, sizeof(text));
    cio_dbf_p32(h + 32, (uint32_t)atoi(text));
    const char *box[] = {"xmin", "ymin", "xmax", "ymax", "zmin", "zmax", "mmin", "mmax"};
    for (int j = 0; j < 8; j++)
    {
        char key[80];
        snprintf(key, sizeof(key), "__cio_shp_%s", box[j]);
        cio_local_text(key, text, sizeof(text));
        cio_shp_putdouble(h + 36 + 8 * j, strtod(text, NULL));
    }
    if ((rc = ctools_output_prepare(&output, filename, !strcmp(replace, "1"))))
        goto done;
    fp = ctools_fopen(output.temporary, "wb");
    if (!fp)
    {
        rc = 603;
        goto done;
    }
    if (fwrite(h, 1, 100, fp) != 100)
    {
        rc = 603;
        goto done;
    }
    if (!strcmp(shx, "1"))
    {
        char path[4096];
        cio_copy_text(path, sizeof(path), filename);
        size_t len = strlen(path);
        if (len < 4)
        {
            rc = 198;
            goto done;
        }
        memcpy(path + len - 4, ".shx", 5);
        if ((rc = ctools_output_prepare(&idx, path, !strcmp(replace, "1"))))
            goto done;
        index = ctools_fopen(idx.temporary, "wb");
        if (!index)
        {
            rc = 603;
            goto done;
        }
        if (fwrite(h, 1, 100, index) != 100)
        {
            rc = 603;
            goto done;
        }
    }
    unsigned records = 0;
    for (int i = SF_in1(); i <= SF_in2(); i++)
    {
        int length = SF_sdatalen(4, i);
        if (length <= 0)
            continue;
        if (length < 8)
        {
            rc = 198;
            goto done;
        }
        unsigned char *blob = malloc((size_t)length + 1);
        if (!blob)
        {
            rc = 909;
            goto done;
        }
        if (SF_strldata(4, i, (char *)blob, length + 1) != length)
        {
            free(blob);
            rc = 609;
            goto done;
        }
        unsigned kind = cio_dbf_u32(blob + 4), words = cio_shp_be32(blob);
        double id, x, y;
        int z, m, point, multi;
        if (!cio_shp_dims(kind, &z, &m, &point, &multi))
        {
            free(blob);
            rc = 198;
            goto done;
        }
        size_t prefix = 4;
        if (kind && !point)
        {
            if (length < (multi ? 44 : 48))
            {
                free(blob);
                rc = 198;
                goto done;
            }
            prefix = multi ? 40 : 44 + (size_t)cio_dbf_u32(blob + 40) * 4 * (kind == 31 ? 2 : 1);
        }
        if (prefix + 4 > (size_t)length)
        {
            free(blob);
            rc = 198;
            goto done;
        }
        if ((rc = (_stata_)->safevdata(1, i, &id)))
        {
            free(blob);
            goto done;
        }
        int64_t offset = ctools_ftell(fp);
        if (offset < 0 || offset / 2 > UINT32_MAX)
        {
            free(blob);
            rc = 603;
            goto done;
        }
        unsigned char record[8];
        cio_shp_putbe32(record, (uint32_t)id);
        cio_shp_putbe32(record + 4, words);
        if (fwrite(record, 1, 8, fp) != 8 || fwrite(blob + 4, 1, prefix, fp) != prefix)
        {
            free(blob);
            rc = 603;
            goto done;
        }
        records++;
        if (index)
        {
            cio_shp_putbe32(record, (uint32_t)(offset / 2));
            if (fwrite(record, 1, 8, index) != 8)
            {
                free(blob);
                rc = 603;
                goto done;
            }
        }
        int limit = i, from = point ? i : i + 1;
        if (kind && !point)
            for (int r = i + 1; r <= SF_in2(); r++)
            {
                double rid;
                if ((rc = (_stata_)->safevdata(1, r, &rid)))
                {
                    free(blob);
                    goto done;
                }
                if (rid != id)
                    break;
                limit = r;
            }
        for (int dimension = 0; kind && dimension < 1 + z + m; dimension++)
        {
            if (dimension && !point)
            {
                size_t range = prefix + 4 + (size_t)(dimension - 1) * 16;
                if (range + 16 > (size_t)length)
                {
                    free(blob);
                    rc = 3300;
                    goto done;
                }
                if (fwrite(blob + range, 1, 16, fp) != 16)
                {
                    free(blob);
                    rc = 603;
                    goto done;
                }
            }
            for (int r = from; r <= limit; r++)
            {
                if ((rc = (_stata_)->safevdata(2, r, &x)) || (rc = (_stata_)->safevdata(3, r, &y)))
                {
                    free(blob);
                    goto done;
                }
                if (!point && (SF_is_missing(x) || SF_is_missing(y)))
                    continue;
                unsigned char value[16];
                size_t bytes;
                if (!dimension)
                {
                    cio_shp_putdouble(value, x);
                    cio_shp_putdouble(value + 8, y);
                    bytes = 16;
                }
                else
                {
                    double v;
                    int var = dimension == 1 ? 5 : 6;
                    if ((rc = (_stata_)->safevdata(var, r, &v)))
                    {
                        free(blob);
                        goto done;
                    }
                    cio_shp_putdouble(value, SF_is_missing(v) ? -1e39 : v);
                    bytes = 8;
                }
                if (fwrite(value, 1, bytes, fp) != bytes)
                {
                    free(blob);
                    rc = 603;
                    goto done;
                }
            }
        }
        free(blob);
        int64_t end = ctools_ftell(fp);
        if (end < 0 || end - offset - 8 != (int64_t)words * 2)
        {
            rc = 198;
            goto done;
        }
    }
    int64_t bytes = ctools_ftell(fp);
    if (bytes < 0 || bytes / 2 > UINT32_MAX)
    {
        rc = 603;
        goto done;
    }
    cio_shp_putbe32(h + 24, (uint32_t)(bytes / 2));
    if (ctools_fseek(fp, 0, SEEK_SET) || fwrite(h, 1, 100, fp) != 100)
    {
        rc = 603;
        goto done;
    }
    if (fclose(fp))
    {
        fp = NULL;
        rc = 603;
        goto done;
    }
    fp = NULL;
    if (index)
    {
        cio_shp_putbe32(h + 24, 50 + 4 * records);
        if (ctools_fseek(index, 0, SEEK_SET) || fwrite(h, 1, 100, index) != 100)
        {
            rc = 603;
            goto done;
        }
        if (fclose(index))
        {
            index = NULL;
            rc = 603;
            goto done;
        }
        index = NULL;
    }
    ctools_output *outputs[] = {&output, &idx};
    rc = ctools_output_commit_group(outputs, idx.temporary ? 2 : 1);
done:
    if (fp)
        fclose(fp);
    if (index)
        fclose(index);
    ctools_output_cleanup(&output);
    ctools_output_cleanup(&idx);
    return rc;
}

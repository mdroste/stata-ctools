#include "ctools_file.h"
#include "cio_internal.h"

int cio_metadata_cb(readstat_metadata_t *m, void *ctx)
{
    cio_dataset *d = ctx;
    d->nvar = readstat_get_var_count(m);
    if (d->nvar < 0 || d->nvar > 32767)
    {
        d->error = 103;
        return READSTAT_HANDLER_ABORT;
    }
    d->cols = calloc((size_t)d->nvar, sizeof(cio_column));
    if (!d->cols && d->nvar)
    {
        d->error = 909;
        return READSTAT_HANDLER_ABORT;
    }
    cio_copy_text(d->datalabel, sizeof(d->datalabel), readstat_get_file_label(m));
    return READSTAT_HANDLER_OK;
}

void cio_import_format(cio_column *c, const char *kind)
{
    const char *f = c->source_format;
    if (c->is_string)
    {
        snprintf(c->format, sizeof(c->format), "%%%ds",
                 c->width > 2045 ? 9
                 : c->width > 9  ? c->width
                                 : 9);
        return;
    }
    cio_copy_text(c->format, sizeof(c->format), "%10.0g");
    if (!strcmp(kind, "spss"))
    {
        char code[40];
        int width = 0, decimals = 0;
        size_t k = 0;
        while (isalpha((unsigned char)f[k]) && k < sizeof(code) - 1)
        {
            code[k] = f[k];
            k++;
        }
        code[k] = 0;
        if (sscanf(f + k, "%d.%d", &width, &decimals) < 1)
            width = 9;
        const char *datefmt = !strcmp(code, "DATE")    ? "%tcDD-Mon-CCYY"
                              : !strcmp(code, "ADATE") ? "%tcNN/DD/CCYY"
                              : !strcmp(code, "EDATE") ? "%tcDD.NN.CCYY"
                              : !strcmp(code, "SDATE") ? "%tcCCYY/NN/DD"
                              : !strcmp(code, "JDATE") ? "%tcCCYYJJJ"
                              : !strcmp(code, "MOYR")  ? "%tcMon_CCYY"
                              : !strcmp(code, "QYR")   ? "%tcq_!Q_CCYY"
                              : !strcmp(code, "WKYR")  ? "%tcww_!W!K_CCYY"
                                                       : NULL;
        if (datefmt)
        {
            c->date_kind = 2;
            cio_copy_text(c->format, sizeof(c->format), datefmt);
        }
        else if (!strcmp(code, "DATETIME"))
        {
            c->date_kind = 2;
            cio_copy_text(c->format, sizeof(c->format), "%tcdd-Mon-CCYY_HH:MM:SS.ss");
        }
        else if (!strcmp(code, "TIME"))
        {
            /* Native converts SPSS TIME using the calendar epoch, but its
               invalid generated display format leaves the numeric default. */
            c->date_kind = 2;
        }
        else if (!strcmp(code, "DTIME"))
        {
            c->date_kind = 3;
            cio_copy_text(c->format, sizeof(c->format), "%tcjjj_HH:MM:SS");
        }
        else if (!strcmp(code, "COMMA"))
            snprintf(c->format, sizeof(c->format), "%%%d.%dfc", width, decimals);
        else if (!strcmp(code, "E"))
            snprintf(c->format, sizeof(c->format), "%%%d.%de", width, decimals);
        /* Native assigns Z's string display format to numeric cio_storage, which
           Stata's public metadata interface rejects. Retain the numeric default. */
        else if (!strcmp(code, "F") || !strcmp(code, "DOT") || !strcmp(code, "PCT"))
            snprintf(c->format, sizeof(c->format), "%%%d.%df", width, decimals);
    }
    else if (!strcmp(kind, "sas"))
    {
        /* SAS7BDAT stores format names separately from widths. Prefix matching
           misclassifies user formats such as DATEAMPM and DATE9. */
        const char *datefmt = !strcmp(f, "DATE")       ? "%td"
                              : !strcmp(f, "YYMMDD")   ? "%tdYYNNDD"
                              : !strcmp(f, "MMDDYY")   ? "%tdNN/DD/YY"
                              : !strcmp(f, "DDMMYY")   ? "%tdDD/NN/YY"
                              : !strcmp(f, "MONYY")    ? "%tdMonthCCYY"
                              : !strcmp(f, "YEAR")     ? "%tdCCYY"
                              : !strcmp(f, "WEEKDATE") ? "%td_DAYNAME,_Month_DD,_CCYY"
                                                       : NULL;
        if (datefmt)
        {
            c->date_kind = 1;
            cio_copy_text(c->format, sizeof(c->format), datefmt);
        }
        else if (!strcmp(f, "JULIAN"))
        {
            cio_copy_text(c->format, sizeof(c->format), "%8.0f");
        }
        else if (!strcmp(f, "DATETIME"))
        {
            c->date_kind = 2;
            cio_copy_text(c->format, sizeof(c->format), "%tc");
        }
        else if (!strcmp(f, "TIME") || !strcmp(f, "TOD") || !strcmp(f, "HHMM"))
        {
            c->date_kind = 3;
            cio_copy_text(c->format, sizeof(c->format),
                          !strcmp(f, "HHMM") ? "%tcHH:MM" : "%tcHH:MM:SS");
        }
    }
    else if (!strncmp(f, "DATETIME", 8) || !strncmp(f, "E8601DT", 7))
    {
        c->date_kind = 2;
        cio_copy_text(c->format, sizeof(c->format), "%tc");
    }
    else if (!strncmp(f, "DATE", 4) || !strncmp(f, "YYMMDD", 6) || !strncmp(f, "MMDDYY", 6) ||
             !strncmp(f, "DDMMYY", 6) || !strncmp(f, "JULIAN", 6))
    {
        c->date_kind = 1;
        cio_copy_text(c->format, sizeof(c->format), "%td");
    }
    else if (!strncmp(f, "TIME", 4))
    {
        c->date_kind = 3;
        cio_copy_text(c->format, sizeof(c->format), "%tcHH:MM:SS");
    }
}

int cio_variable_cb(int index, readstat_variable_t *v, const char *set, void *ctx)
{
    cio_dataset *d = ctx;
    if (index < 0 || index >= d->nvar)
    {
        d->error = 610;
        return READSTAT_HANDLER_ABORT;
    }
    cio_column *c = &d->cols[index];
    cio_copy_text(c->name, sizeof(c->name), readstat_variable_get_name(v));
    cio_copy_text(c->label, sizeof(c->label), readstat_variable_get_label(v));
    cio_copy_text(c->source_format, sizeof(c->source_format), readstat_variable_get_format(v));
    cio_copy_text(c->labelset, sizeof(c->labelset), set);
    c->is_string = readstat_variable_get_type_class(v) == READSTAT_TYPE_CLASS_STRING;
    c->width = 1;
    if ((!strcmp(d->kind, "sasxport5") || d->raw_storage) && c->is_string)
        c->width = (int)readstat_variable_get_storage_width(v);
    if (d->raw_storage && !c->is_string)
        c->force_double = 1;
    cio_import_format(c, d->kind);
    if (!strcmp(d->kind, "sasxport5"))
    {
        cio_copy_text(c->labelset, sizeof(c->labelset), c->source_format);
        if (c->date_kind || !strncmp(c->source_format, "BEST", 4))
            c->labelset[0] = 0;
        if (c->date_kind == 1)
            cio_copy_text(c->format, sizeof(c->format), "%d");
        else if (!c->is_string)
            cio_copy_text(c->format, sizeof(c->format), "%10.0g");
        c->date_kind = 0;
    }
    return READSTAT_HANDLER_OK;
}

int cio_value_cb(int row, readstat_variable_t *v, readstat_value_t value, void *ctx)
{
    cio_dataset *d = ctx;
    int j = readstat_variable_get_index(v);
    if (row < 0 || j < 0 || j >= d->nvar)
    {
        d->error = 610;
        return READSTAT_HANDLER_ABORT;
    }
    int rc = cio_reserve_rows(d, (size_t)row + 1);
    if (rc)
    {
        d->error = rc;
        return READSTAT_HANDLER_ABORT;
    }
    if ((size_t)row >= d->nobs)
        d->nobs = (size_t)row + 1;
    cio_column *c = &d->cols[j];
    if (c->is_string)
    {
        const char *s = readstat_string_value(value);
        if (!strcmp(d->kind, "spss") && readstat_value_is_defined_missing(value, v))
            s = "";
        if (!s)
            s = "";
        size_t len = strlen(s);
        if (len > INT_MAX)
        {
            d->error = 109;
            return READSTAT_HANDLER_ABORT;
        }
        char *p = malloc(len + 1);
        if (!p)
        {
            d->error = 909;
            return READSTAT_HANDLER_ABORT;
        }
        memcpy(p, s, len + 1);
        free(c->strings[row]);
        c->strings[row] = p;
        c->lengths[row] = len;
        if ((int)len > c->width)
            c->width = (int)len;
    }
    else
    {
        double x = cio_number_value(value);
        if (!strcmp(d->kind, "spss") && readstat_value_is_defined_missing(value, v))
            x = SV_missval;
        /* Native XPORT import collapses SAS special missings. */
        if (strcmp(d->kind, "sasxport5") && SF_is_missing(x))
            x = SV_missval;
        if (!SF_is_missing(x) && c->date_kind)
        {
            if (!strcmp(d->kind, "spss"))
                x = c->date_kind == 3 ? x * 1000 : x * 1000 - 11903760000000.0;
            else if (!strncmp(d->kind, "sasxport", 8))
            { /* Native transport values retain their original
                 units. */
            }
            else if (c->date_kind >= 2)
                x *= 1000;
        }
        c->numbers[row] = x;
    }
    return READSTAT_HANDLER_OK;
}

int cio_label_cb(const char *set, readstat_value_t value, const char *text, void *ctx)
{
    cio_dataset *d = ctx;
    if (readstat_value_type_class(value) == READSTAT_TYPE_CLASS_STRING)
        return READSTAT_HANDLER_OK;
    cio_label *p = cio_realloc(d->labels, d->nlabels + 1, sizeof(cio_label));
    if (!p)
    {
        d->error = 909;
        return READSTAT_HANDLER_ABORT;
    }
    d->labels = p;
    cio_label *l = &p[d->nlabels];
    memset(l, 0, sizeof(*l));
    l->set = strdup(set ? set : "");
    l->text = strdup(text ? text : "");
    l->value = cio_number_value(value);
    d->nlabels++;
    if (!l->set || !l->text)
    {
        d->error = 909;
        return READSTAT_HANDLER_ABORT;
    }
    return READSTAT_HANDLER_OK;
}

void cio_error_cb(const char *msg, void *ctx)
{
    (void)ctx;
    ctools_error("cio", "%s", msg);
}

int cio_progress_cb(double progress, void *ctx)
{
    (void)progress;
    (void)ctx;
    return SF_poll() ? READSTAT_HANDLER_ABORT : READSTAT_HANDLER_OK;
}

int cio_read_xpf(cio_dataset *dataset, const char *path)
{
    cio_dataset labels;
    memset(&labels, 0, sizeof(labels));
    cio_copy_text(labels.kind, sizeof(labels.kind), "labels");
    readstat_parser_t *p = readstat_parser_init();
    if (!p)
        return 909;
    readstat_set_metadata_handler(p, cio_metadata_cb);
    readstat_set_variable_handler(p, cio_variable_cb);
    readstat_set_value_handler(p, cio_value_cb);
    readstat_error_t err = readstat_parse_xport(p, path, &labels);
    readstat_parser_free(p);
    int rc = err == READSTAT_OK ? 0 : 610, fmt = -1, start = -1, text = -1;
    if (!rc)
    {
        for (int j = 0; j < labels.nvar; j++)
        {
            char name[301];
            cio_copy_text(name, sizeof(name), labels.cols[j].name);
            for (char *c = name; *c; c++)
                *c = (char)toupper((unsigned char)*c);
            if (!strcmp(name, "FMTNAME"))
                fmt = j;
            else if (!strcmp(name, "START"))
                start = j;
            else if (!strcmp(name, "LABEL"))
                text = j;
        }
        if (fmt < 0 || start < 0 || text < 0 || !labels.cols[fmt].is_string ||
            !labels.cols[text].is_string)
            rc = 610;
    }
    for (size_t i = 0; !rc && i < labels.nobs; i++)
    {
        char *set = labels.cols[fmt].strings[i], *label = labels.cols[text].strings[i];
        double val;
        if (labels.cols[start].is_string)
        {
            const char *str = labels.cols[start].strings[i];
            if (!str)
                continue;
            if (str[0] == '.')
                val = cio_tagged_missing(str[1]);
            else
            {
                char *end;
                val = strtod(str, &end);
                if (end == str || *end)
                    continue;
            }
        }
        else
            val = labels.cols[start].numbers[i];
        cio_label *all = cio_realloc(dataset->labels, dataset->nlabels + 1, sizeof(cio_label));
        if (!all)
        {
            rc = 909;
            break;
        }
        dataset->labels = all;
        cio_label *l = &all[dataset->nlabels++];
        memset(l, 0, sizeof(*l));
        l->set = strdup(set ? set : "");
        l->text = strdup(label ? label : "");
        l->value = val;
        if (!l->set || !l->text)
            rc = 909;
    }
    cio_dataset_free(&labels);
    return rc;
}

int cio_scan(cio_dataset *dataset, const char *kind)
{
    char path[4096], encoding[128], catalog[4096], xpf[4096], raw_storage[16];
    cio_local_text("__cio_filename", path, sizeof(path));
    cio_local_text("__cio_encoding", encoding, sizeof(encoding));
    cio_local_text("__cio_bcat", catalog, sizeof(catalog));
    cio_local_text("__cio_xpf", xpf, sizeof(xpf));
    cio_local_text("__cio_rawstorage", raw_storage, sizeof(raw_storage));
    cio_dataset_free(dataset);
    cio_copy_text(dataset->kind, sizeof(dataset->kind), kind);
    dataset->raw_storage = !strcmp(kind, "sasxport8") && atoi(raw_storage) != 0;
    if (!strcmp(kind, "dbase") || !strcmp(kind, "xls") || !strcmp(kind, "xlsx") ||
        !strcmp(kind, "shp"))
    {
        int rc = !strcmp(kind, "xlsx")  ? cio_scan_xlsx(dataset, path)
                 : !strcmp(kind, "xls") ? cio_scan_xls(dataset, path)
                 : !strcmp(kind, "shp") ? cio_scan_shp(dataset, path)
                                        : cio_scan_dbf(dataset, path);
        if (rc)
        {
            cio_dataset_free(dataset);
            return rc;
        }
        if ((rc = cio_save_number("__cio_nobs", dataset->nobs)))
            return rc;
        if ((rc = cio_save_number("__cio_nvar", (size_t)dataset->nvar)))
            return rc;
        if ((rc = cio_save_number("__cio_nlabels", 0)))
            return rc;
        return cio_save("__cio_datalabel", "");
    }
    readstat_parser_t *p = readstat_parser_init();
    if (!p)
        return 909;
    readstat_set_metadata_handler(p, cio_metadata_cb);
    readstat_set_variable_handler(p, cio_variable_cb);
    readstat_set_value_handler(p, cio_value_cb);
    readstat_set_value_label_handler(p, cio_label_cb);
    readstat_set_error_handler(p, cio_error_cb);
    readstat_set_progress_handler(p, cio_progress_cb);
    if (*encoding)
        readstat_set_file_character_encoding(p, encoding);
    readstat_set_handler_character_encoding(p, "UTF-8");
    readstat_error_t err;
    cio_xport_view view;
    memset(&view, 0, sizeof(view));
    if (!strcmp(kind, "sasxport5"))
    {
        char member[301];
        cio_local_text("__cio_member", member, sizeof(member));
        int rc = cio_xport_select(path, member, p, &view);
        if (rc)
        {
            readstat_parser_free(p);
            cio_dataset_free(dataset);
            return rc;
        }
    }
    if (!strcmp(kind, "sas"))
        err = readstat_parse_sas7bdat(p, path, dataset);
    else if (!strcmp(kind, "spss"))
        err = readstat_parse_sav(p, path, dataset);
    else
        err = readstat_parse_xport(p, path, dataset);
    if (err == READSTAT_OK && *catalog)
    {
        readstat_set_metadata_handler(p, NULL);
        readstat_set_variable_handler(p, NULL);
        readstat_set_value_handler(p, NULL);
        err = readstat_parse_sas7bcat(p, catalog, dataset);
    }
    readstat_parser_free(p);
    if (err != READSTAT_OK)
    {
        int rc = dataset->error ? dataset->error : !strcmp(kind, "sasxport5") ? 459 : 692;
        if (!dataset->error && err == READSTAT_ERROR_OPEN)
            rc = 601;
        if (!dataset->error && err == READSTAT_ERROR_MALLOC)
            rc = 909;
        ctools_error("cimport", "%s", readstat_error_message(err));
        cio_dataset_free(dataset);
        return rc;
    }
    if (!strcmp(kind, "sasxport5") && *xpf)
    {
        int rc = cio_read_xpf(dataset, xpf);
        if (rc)
        {
            cio_dataset_free(dataset);
            return rc;
        }
    }
    for (int j = 0; j < dataset->nvar; j++)
        if (dataset->cols[j].is_string)
            cio_import_format(&dataset->cols[j], kind);
    int rc = cio_save_number("__cio_nobs", dataset->nobs);
    if (!rc)
        rc = cio_save_number("__cio_nvar", (size_t)dataset->nvar);
    if (!rc)
        rc = cio_save_number("__cio_nlabels", dataset->nlabels);
    if (!rc)
        rc = cio_save("__cio_datalabel", dataset->datalabel);
    if (!rc && !strcmp(kind, "sasxport5"))
    {
        size_t width = 4;
        for (int j = 0; j < dataset->nvar; j++)
            width += dataset->cols[j].is_string ? (size_t)dataset->cols[j].width : 8;
        if (dataset->nobs && width > SIZE_MAX / dataset->nobs)
            rc = 920;
        else
            rc = cio_save_number("__cio_size", width * dataset->nobs);
    }
    return rc;
}

ssize_t cio_write_cb(const void *data, size_t len, void *ctx)
{
    FILE *f = ctx;
    size_t n = fwrite(data, 1, len, f);
    return n == len ? (ssize_t)n : -1;
}

int cio_xport_names(int nvar, char (*names)[301], char (*sets)[301], int rename)
{
    char (*original_sets)[301] = calloc((size_t)nvar, 301);
    if (!original_sets)
        return 909;
    int rc = 0;
    for (int field = 0; field < 2 && !rc; field++)
        for (int j = 0; j < nvar; j++)
        {
            char original[301], base[301];
            cio_export_meta(j + 1, field ? "labelset" : "name", original, sizeof(original));
            char *out = field ? sets[j] : names[j];
            if (field)
            {
                cio_copy_text(original_sets[j], 301, original);
                int previous = -1;
                for (int k = 0; k < j; k++)
                    if (!strcmp(original, original_sets[k]))
                    {
                        previous = k;
                        break;
                    }
                if (previous >= 0)
                {
                    cio_copy_text(out, 301, sets[previous]);
                    continue;
                }
            }
            if (!*original)
            {
                out[0] = 0;
                continue;
            }
            cio_copy_text(base, sizeof(base), original);
            for (char *p = base; *p; p++)
                *p = (char)toupper((unsigned char)*p);
            if (strlen(base) > 8 && !rename)
            {
                rc = 198;
                break;
            }
            base[8] = 0;
            cio_copy_text(out, 301, base);
            for (int suffix = 2;; suffix++)
            {
                int collision = 0;
                for (int k = 0; k < j; k++)
                    if (!strcmp(out, field ? sets[k] : names[k]))
                    {
                        collision = 1;
                        break;
                    }
                if (!collision)
                    break;
                if (!rename)
                {
                    rc = 198;
                    break;
                }
                char number[24];
                snprintf(number, sizeof(number), "%d", suffix);
                size_t keep = 8 - strlen(number);
                snprintf(out, 301, "%.*s%s", (int)keep, base, number);
            }
        }
    free(original_sets);
    return rc;
}

int cio_do_export(const char *kind)
{
    if (!strcmp(kind, "shp"))
        return cio_shp_export();
    if (!strcmp(kind, "dbase"))
        return cio_export_dbf();
    char path[4096], replace[10], label[1025];
    ctools_output companions[2] = {{0}};
    ctools_output output;
    memset(&output, 0, sizeof(output));
    cio_local_text("__cio_filename", path, sizeof(path));
    cio_local_text("__cio_replace", replace, sizeof(replace));
    int nvar = SF_nvars(), rc = 0;
    size_t nobs = 0;
    if (nvar < 1)
        return 102;
    for (int i = SF_in1(); i <= SF_in2(); i++)
        if (SF_ifobs(i))
            nobs++;
    readstat_writer_t *w = readstat_writer_init();
    if (!w)
        return 909;
    char metacount[40];
    cio_local_text("__cio_nmeta", metacount, sizeof(metacount));
    int nmeta = atoi(metacount);
    if (nmeta < nvar)
        nmeta = nvar;
    char (*names)[301] = calloc((size_t)nmeta, 301), (*sets)[301] = calloc((size_t)nmeta, 301);
    readstat_variable_t **vars = calloc((size_t)nvar, sizeof(*vars));
    int *dates = calloc((size_t)nvar, sizeof(int));
    char *s = NULL;
    size_t stringcap = 0, stringneed = 0; /* stringcap is always the real size of s */
    FILE *fp = NULL;
    if (!vars || !dates || !names || !sets)
    {
        rc = 909;
        goto done;
    }
    if (!strcmp(kind, "sasxport5"))
    {
        char rename[8];
        cio_local_text("__cio_rename", rename, sizeof(rename));
        if ((rc = cio_xport_names(nmeta, names, sets, !strcmp(rename, "1"))))
            goto done;
    }
    if (strcmp(kind, "sasxport5"))
        for (int j = 0; j < nmeta; j++)
            cio_export_meta(j + 1, "labelset", sets[j], 301);
    readstat_writer_set_error_handler(w, cio_error_cb);
    readstat_set_data_writer(w, cio_write_cb);
    cio_local_text("__cio_datalabel", label, sizeof(label));
    readstat_writer_set_file_label(w, label);
    char table[40];
    cio_local_text("__cio_table", table, sizeof(table));
    if (*table)
    {
        if (!strcmp(kind, "sasxport5"))
            table[8] = 0;
        readstat_writer_set_table_name(w, table);
    }
    if (!strncmp(kind, "sasxport", 8))
        readstat_writer_set_file_format_version(w, !strcmp(kind, "sasxport5") ? 5 : 8);
    else if (!strcmp(kind, "spss"))
        readstat_writer_set_compression(w, READSTAT_COMPRESS_ROWS);
    for (int j = 0; j < nvar; j++)
    {
        char name[301], type[40], format[257], set[301], count[40];
        cio_export_meta(j + 1, "name", name, sizeof(name));
        cio_export_meta(j + 1, "type", type, sizeof(type));
        cio_export_meta(j + 1, "format", format, sizeof(format));
        cio_export_meta(j + 1, "label", label, sizeof(label));
        if (!strcmp(kind, "sasxport5"))
            cio_copy_text(name, sizeof(name), names[j]);
        int isstr = SF_var_is_string(j + 1), width = isstr ? atoi(type + 3) : 0;
        if (SF_var_is_strl(j + 1))
        {
            width = 1;
            for (int i = SF_in1(); i <= SF_in2(); i++)
                if (SF_ifobs(i))
                {
                    if (SF_var_is_binary(j + 1, i))
                    {
                        rc = 109;
                        goto done;
                    }
                    int length = SF_sdatalen(j + 1, i);
                    if (length < 0)
                    {
                        rc = 109;
                        goto done;
                    }
                    size_t cap = (size_t)length + 1;
                    if (cap > stringcap)
                    {
                        char *tmp = realloc(s, cap);
                        if (!tmp)
                        {
                            rc = 909;
                            goto done;
                        }
                        s = tmp;
                        stringcap = cap;
                    }
                    if (SF_strldata(j + 1, i, s, (int)stringcap) != length)
                    {
                        rc = 609;
                        goto done;
                    }
                    s[length] = 0;
                    size_t logical = strlen(s) + (!strcmp(kind, "sasxport8") ? 1 : 0);
                    if (logical > INT_MAX)
                    {
                        rc = 109;
                        goto done;
                    }
                    if ((int)logical > width)
                        width = (int)logical;
                }
        }
        vars[j] = readstat_add_variable(
            w, name, isstr ? READSTAT_TYPE_STRING : READSTAT_TYPE_DOUBLE, (size_t)width);
        if (!vars[j])
        {
            rc = 909;
            goto done;
        }
        readstat_variable_set_label(vars[j], label);
        if (isstr && !strcmp(kind, "sasxport8"))
        {
            char stringformat[80];
            int display = 0;
            sscanf(format, "%%%d", &display);
            if (SF_var_is_strl(j + 1))
                display = 8;
            snprintf(stringformat, sizeof(stringformat), "$%d", display > 0 ? display : 8);
            readstat_variable_set_format(vars[j], stringformat);
        }
        if (!isstr)
        {
            char fmt[80];
            if (!strncmp(format, "%td", 3))
                dates[j] = 1;
            else if (!strncmp(format, "%tc", 3) || !strncmp(format, "%tC", 3))
                dates[j] = 2;
            if (!strcmp(kind, "spss"))
            {
                if (dates[j] == 1)
                    cio_copy_text(fmt, sizeof(fmt), "DATE11");
                else if (dates[j] == 2)
                    cio_copy_text(fmt, sizeof(fmt), "DATETIME23.2");
                else
                {
                    int fw = 0, dec = 0;
                    sscanf(format, "%%%d.%d", &fw, &dec);
                    if (fw < 8)
                        fw = 8;
                    snprintf(fmt, sizeof(fmt), "F%d.%d", fw + 3, dec);
                }
            }
            else
            {
                if (dates[j] == 1)
                    cio_copy_text(fmt, sizeof(fmt), "DATE");
                else if (dates[j] == 2)
                    cio_copy_text(fmt, sizeof(fmt), !strcmp(kind, "sasxport5") ? "" : "DATETIME");
                else
                    cio_copy_text(fmt, sizeof(fmt), "");
            }
            readstat_variable_set_format(vars[j], fmt);
        }
        cio_export_meta(j + 1, "labelset", set, sizeof(set));
        cio_export_meta(j + 1, "lcount", count, sizeof(count));
        if (!strcmp(kind, "sasxport5"))
            cio_copy_text(set, sizeof(set), sets[j]);
        if (*set && !strcmp(kind, "sasxport5"))
        {
            for (char *p = set; *p; p++)
                *p = (char)toupper((unsigned char)*p);
            readstat_variable_set_format(vars[j], set);
        }
        if (*set && !strcmp(kind, "spss"))
        {
            readstat_label_set_t *ls = NULL;
            for (long k = 0; k < w->label_sets_count; k++)
                if (!strcmp(w->label_sets[k]->name, set))
                {
                    ls = w->label_sets[k];
                    break;
                }
            int new_set = !ls;
            if (!ls)
                ls = readstat_add_label_set(w, READSTAT_TYPE_DOUBLE, set);
            if (!ls)
            {
                rc = 909;
                goto done;
            }
            readstat_variable_set_label_set(vars[j], ls);
            int n = atoi(count);
            for (int k = 1; new_set && k <= n; k++)
            {
                char key[100], val[64], text[32768];
                snprintf(key, sizeof(key), "__cio_lv_%d_%d", j + 1, k);
                cio_local_text(key, val, sizeof(val));
                snprintf(key, sizeof(key), "__cio_lt_%d_%d", j + 1, k);
                cio_local_text(key, text, sizeof(text));
                if (*val != '.')
                    readstat_label_double_value(ls, strtod(val, NULL), text);
            }
        }
        if (isstr && (size_t)width + 1 > stringneed)
            stringneed = (size_t)width + 1;
    }
    if (stringneed > stringcap)
    {
        char *tmp = realloc(s, stringneed);
        if (!tmp)
        {
            rc = 909;
            goto done;
        }
        s = tmp;
        stringcap = stringneed;
    }
    /* Use a sibling temporary and atomic publish; failed writes never damage
       an existing destination. Never unlink the destination before writing. */
    rc = ctools_output_prepare(&output, path, !strcmp(replace, "1"));
    if (rc)
        goto done;
    fp = ctools_fopen(output.temporary, "wb");
    if (!fp)
    {
        rc = 603;
        goto done;
    }
    readstat_error_t err;
    if (!strcmp(kind, "sas"))
        err = readstat_begin_writing_sas7bdat(w, fp, (long)nobs);
    else if (!strcmp(kind, "spss"))
        err = readstat_begin_writing_sav(w, fp, (long)nobs);
    else
        err = readstat_begin_writing_xport(w, fp, (long)nobs);
    if (err != READSTAT_OK)
        goto codec_error;
    for (int i = SF_in1(); i <= SF_in2(); i++)
    {
        if (!SF_ifobs(i))
            continue;
        if ((i & 4095) == 0 && SF_poll())
        {
            rc = 1;
            goto done;
        }
        if ((err = readstat_begin_row(w)) != READSTAT_OK)
            goto codec_error;
        for (int j = 0; j < nvar; j++)
        {
            if (SF_var_is_string(j + 1))
            {
                if (SF_var_is_strl(j + 1))
                {
                    int length = SF_sdatalen(j + 1, i);
                    if (length < 0 || (size_t)length >= stringcap)
                    {
                        rc = 109;
                        goto done;
                    }
                    if (SF_strldata(j + 1, i, s, (int)stringcap) != length)
                    {
                        rc = 609;
                        goto done;
                    }
                    s[length] = 0;
                }
                else if ((rc = SF_sdata(j + 1, i, s)))
                    goto done;
                err = readstat_insert_string_value(w, vars[j], s);
            }
            else
            {
                double x;
                if ((rc = (_stata_)->safevdata(j + 1, i, &x)))
                    goto done;
                if (SF_is_missing(x))
                {
                    uint64_t bits, base;
                    double mv = SV_missval;
                    memcpy(&bits, &x, 8);
                    memcpy(&base, &mv, 8);
                    int n = (int)((bits - base) >> 40);
                    if (n && strcmp(kind, "spss"))
                        err = readstat_insert_tagged_missing_value(w, vars[j], (char)('A' + n - 1));
                    else
                        err = readstat_insert_missing_value(w, vars[j]);
                }
                else
                {
                    if (!strcmp(kind, "spss") && dates[j])
                        x = ((dates[j] == 1 ? x * 86400000 : x) + 11903760000000.0) / 1000;
                    else if (!strcmp(kind, "sasxport8") && dates[j] == 1)
                        x *= 86400000;
                    else if (!strncmp(kind, "sasxport", 8))
                    { /* Transport preserves Stata numeric units. */
                    }
                    else if (dates[j] >= 2)
                        x /= 1000;
                    err = readstat_insert_double_value(w, vars[j], x);
                }
            }
            if (err != READSTAT_OK)
                goto codec_error;
        }
        if ((err = readstat_end_row(w)) != READSTAT_OK)
            goto codec_error;
    }
    if ((err = readstat_end_writing(w)) != READSTAT_OK)
        goto codec_error;
    if (fclose(fp))
    {
        fp = NULL;
        rc = 603;
        goto done;
    }
    fp = NULL;
    if (!strcmp(kind, "sasxport5"))
    {
        char xpf[4096];
        cio_local_text("__cio_xpfoutput", xpf, sizeof(xpf));
        if (*xpf && (rc = cio_export_xpf(nmeta, xpf, !strcmp(replace, "1"), sets, &companions[0])))
            goto done;
    }
    char sascode[4096];
    cio_local_text("__cio_sascode", sascode, sizeof(sascode));
    if (*sascode && (rc = cio_export_sascode(nmeta, sets, sascode, kind, !strcmp(replace, "1"), &companions[1])))
        goto done;
    ctools_output *outputs[3] = {&output};
    size_t count = 1;
    for (int i = 0; i < 2; i++) if (companions[i].temporary) outputs[count++] = &companions[i];
    rc = ctools_output_commit_group(outputs, count);
    goto done;
codec_error:
    ctools_error("cexport", "%s", readstat_error_message(err));
    rc = 610;
done:
    if (fp)
        fclose(fp);
    ctools_output_cleanup(&output);
    for (int i = 0; i < 2; i++) ctools_output_cleanup(&companions[i]);
    free(s);
    free(vars);
    free(dates);
    free(names);
    free(sets);
    readstat_writer_free(w);
    return rc;
}

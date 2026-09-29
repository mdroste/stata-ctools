/* Value-label companion used by Stata's version-5 XPORT command. */
static int export_xpf(int nvar, const char *path, int replace, char (*mapped)[301]) {
    int rc = 0, nsets = 0;
    long rows = 0;
    char **sets = calloc((size_t)nvar, sizeof(char *));
    int *columns = calloc((size_t)nvar, sizeof(int));
    readstat_writer_t *w = NULL;
    FILE *fp = NULL;
    cexport_output output;
    memset(&output, 0, sizeof(output));
    if (!sets || !columns) {
        rc = 909;
        goto done;
    }
    for (int j = 1; j <= nvar; j++) {
        char set[301], count[40];
        copy_text(set, sizeof(set), mapped[j - 1]);
        if (!*set)
            continue;
        int known = 0;
        for (int k = 0; k < nsets; k++)
            if (!strcmp(set, sets[k])) {
                known = 1;
                break;
            }
        if (known)
            continue;
        if (strlen(set) > 8) {
            rc = 198;
            goto done;
        }
        sets[nsets] = strdup(set);
        if (!sets[nsets]) {
            rc = 909;
            goto done;
        }
        columns[nsets++] = j;
        export_meta(j, "lcount", count, sizeof(count));
        rows += atoi(count);
    }
    if (!rows)
        goto done;
    w = readstat_writer_init();
    if (!w) {
        rc = 909;
        goto done;
    }
    readstat_set_data_writer(w, write_cb);
    readstat_writer_set_file_format_version(w, 5);
    readstat_writer_set_table_name(w, "FORMATS");
    readstat_variable_t *fmt = readstat_add_variable(w, "FMTNAME", READSTAT_TYPE_STRING, 8),
                        *start = readstat_add_variable(w, "START", READSTAT_TYPE_STRING, 12),
                        *label = readstat_add_variable(w, "LABEL", READSTAT_TYPE_STRING, 40);
    if (!fmt || !start || !label) {
        rc = 909;
        goto done;
    }
    if ((rc = cexport_output_prepare(&output, path, replace)))
        goto done;
    fp = fopen(output.temporary, "wb");
    if (!fp) {
        rc = 603;
        goto done;
    }
    readstat_error_t err = readstat_begin_writing_xport(w, fp, rows);
    if (err != READSTAT_OK) {
        rc = 610;
        goto done;
    }
    for (int j = 0; j < nsets; j++) {
        char set[301], count[40];
        copy_text(set, sizeof(set), sets[j]);
        for (char *p = set; *p; p++)
            *p = (char)toupper((unsigned char)*p);
        export_meta(columns[j], "lcount", count, sizeof(count));
        for (int k = 1; k <= atoi(count); k++) {
            char key[100], value[64], text[32768], num[32];
            snprintf(key, sizeof(key), "__cio_lv_%d_%d", columns[j], k);
            local_text(key, value, sizeof(value));
            snprintf(key, sizeof(key), "__cio_lt_%d_%d", columns[j], k);
            local_text(key, text, sizeof(text));
            if (value[0] == '.') {
                copy_text(num, sizeof(num), value);
                if (num[1])
                    num[1] = (char)toupper((unsigned char)num[1]);
            } else
                snprintf(num, sizeof(num), "%.11g", strtod(value, NULL));
            /* Version-5 format labels are limited to 40 bytes. */
            text[40] = 0;
            if ((err = readstat_begin_row(w)) != READSTAT_OK ||
                (err = readstat_insert_string_value(w, fmt, set)) != READSTAT_OK ||
                (err = readstat_insert_string_value(w, start, num)) != READSTAT_OK ||
                (err = readstat_insert_string_value(w, label, text)) != READSTAT_OK ||
                (err = readstat_end_row(w)) != READSTAT_OK) {
                rc = 610;
                goto done;
            }
        }
    }
    if ((err = readstat_end_writing(w)) != READSTAT_OK) {
        rc = 610;
        goto done;
    }
    if (fclose(fp)) {
        fp = NULL;
        rc = 603;
        goto done;
    }
    fp = NULL;
    rc = cexport_output_commit(&output);
done:
    if (fp)
        fclose(fp);
    cexport_output_cleanup(&output);
    readstat_writer_free(w);
    for (int j = 0; j < nsets; j++)
        free(sets[j]);
    free(sets);
    free(columns);
    return rc;
}

/* SAS PROC FORMAT text is serialized here, with the same naming and order as
   the transport dictionaries. The ado wrapper supplies label metadata only. */
static int export_sascode(int nvar, char (*mapped)[301], const char *path, const char *kind,
                          int replace) {
    cexport_output output;
    memset(&output, 0, sizeof(output));
    FILE *fp = NULL;
    int rc = 0;
    int any = 0;
    for (int j = 0; j < nvar; j++)
        if (*mapped[j])
            any = 1;
    if (!any)
        return 0;
    if ((rc = cexport_output_prepare(&output, path, replace)))
        goto done;
    fp = fopen(output.temporary, "wb");
    if (!fp) {
        rc = 603;
        goto done;
    }
    char cwd[4096], transport[4096];
    local_text("__cio_cwd", cwd, sizeof(cwd));
    local_text("__cio_transportpath", transport, sizeof(transport));
    int v5 = !strcmp(kind, "sasxport5");
    if (fprintf(fp,
                "%slibname datapath '%s' ;\nlibname xptfile xport '%s' ;\n\nproc copy in = xptfile "
                "out = datapath ;\n\nproc format library = work ;\n",
                v5 ? "\n" : "", cwd, transport) < 0) {
        rc = 603;
        goto done;
    }
    for (int j = 0; j < nvar; j++) {
        if (!*mapped[j])
            continue;
        int seen = 0;
        for (int k = 0; k < j; k++)
            if (!strcmp(mapped[j], mapped[k])) {
                seen = 1;
                break;
            }
        if (seen)
            continue;
        char set[301], count[40];
        copy_text(set, sizeof(set), mapped[j]);
        for (char *p = set; *p; p++)
            *p = (char)toupper((unsigned char)*p);
        export_meta(j + 1, "lcount", count, sizeof(count));
        int n = atoi(count);
        if (fprintf(fp, "    value %s\n", set) < 0) {
            rc = 603;
            goto done;
        }
        for (int k = 1; k <= n; k++) {
            char key[100], value[80], text[32768], num[80];
            snprintf(key, sizeof(key), "__cio_lv_%d_%d", j + 1, k);
            local_text(key, value, sizeof(value));
            snprintf(key, sizeof(key), "__cio_lt_%d_%d", j + 1, k);
            local_text(key, text, sizeof(text));
            if (value[0] == '.') {
                copy_text(num, sizeof(num), value);
                if (num[1])
                    num[1] = (char)toupper((unsigned char)num[1]);
            } else
                snprintf(num, sizeof(num), "%.11g", strtod(value, NULL));
            if (fprintf(fp, "%12s = '%s'%s\n", num, text, k == n ? " ;" : "") < 0) {
                rc = 603;
                goto done;
            }
        }
        if (fputc('\n', fp) == EOF) {
            rc = 603;
            goto done;
        }
    }
    if (fprintf(fp, "quit ;\n\n") < 0 || ferror(fp)) {
        rc = 603;
        goto done;
    }
    if (fclose(fp)) {
        fp = NULL;
        rc = 603;
        goto done;
    }
    fp = NULL;
    rc = cexport_output_commit(&output);
done:
    if (fp)
        fclose(fp);
    cexport_output_cleanup(&output);
    return rc;
}

#include "cio_internal.h"

/* Reuse the streaming ZIP/XML parser, then transfer metadata and values through
   the same bounded interface as XLS and statistical formats. */
#include "cimport/cimport_xlsx.h"
#include "ctools_platform.h"
/* Dense import costs about 26 bytes a cell: 10 in the parser's storage, 8 in
 * the staged column and 8 in Stata. A shape beyond half of physical memory
 * fails with r(909) before any of it is allocated; a lone value at XFD1048576
 * would otherwise ask for 17 billion cells (about 450 GB). */
size_t cio_xlsx_cell_limit(void)
{
    uint64_t memory = ctools_physical_memory_bytes();
    if (!memory)
        memory = UINT64_C(16) << 30;
    uint64_t cells = memory / 2 / 26;
    return cells > SIZE_MAX ? SIZE_MAX : (size_t)cells;
}
/* One package part may inflate to at most half of physical memory, or to
 * c(max_memory) when that is set lower (the ado passes it). */
uint64_t cio_xlsx_part_limit(void)
{
    uint64_t memory = ctools_physical_memory_bytes();
    uint64_t limit = memory ? memory / 2 : UINT64_C(8) << 30;
    char text[40] = "";
    cio_local_text("__cio_maxmemory", text, sizeof(text));
    double cap = atof(text);
    if (cap > 0 && cap < (double)limit)
        limit = (uint64_t)cap;
    return limit;
}
/* Every stored cell lives in the parser's column-major arrays. */
int cio_xlsx_stored(XLSXContext *x, int row, int col)
{
    return x->cm_active && row >= 1 && (size_t)row <= x->cm_num_rows && col >= 0 &&
           col < x->cm_num_cols;
}
int cio_xlsx_raw(XLSXContext *x, int row, int col, double *number, const char **text)
{
    *number = 0;
    *text = "";
    if (!cio_xlsx_stored(x, row, col))
        return XLSX_CELL_EMPTY;
    int type = x->cm_types[col][row - 1];
    *number = x->cm_numeric[col][row - 1];
    if (type == XLSX_CELL_STRING)
    {
        uintptr_t pointer = 0;
        memcpy(&pointer, number, sizeof(pointer));
        *text = (const char *)pointer;
    }
    if (type == XLSX_CELL_SHARED_STRING)
    {
        int index = (int)*number;
        if (index < 0 || (uint32_t)index >= x->num_shared_strings)
            return XLSX_CELL_EMPTY;
        *text = x->shared_strings[index];
    }
    if (!*text)
        *text = "";
    return type;
}
uint8_t cio_xlsx_native_format(XLSXContext *x, int row, int col)
{
    return cio_xlsx_stored(x, row, col) ? x->cm_formats[col][row - 1] : 0;
}
int cio_xlsx_numeric(int type)
{
    return type == XLSX_CELL_NUMBER || type == XLSX_CELL_DATE || type == XLSX_CELL_DATETIME ||
           type == XLSX_CELL_BOOLEAN;
}
int cio_scan_xlsx(cio_dataset *dataset, const char *path)
{
    XLSXContext x;
    xlsx_context_init(&x);
    int rc = 0;
    char sheet[1024], range[100], first[8], all[8], allfmt[80], maxvar[40];
    x.max_cells = cio_xlsx_cell_limit();
    x.max_part_bytes = cio_xlsx_part_limit();
    cio_local_text("__cio_maxvar", maxvar, sizeof(maxvar));
    cio_local_text("__cio_sheet", sheet, sizeof(sheet));
    cio_local_text("__cio_cellrange", range, sizeof(range));
    cio_local_text("__cio_firstrow", first, sizeof(first));
    cio_local_text("__cio_allstring", all, sizeof(all));
    cio_local_text("__cio_allstringfmt", allfmt, sizeof(allfmt));
    x.firstrow = !strcmp(first, "1");
    x.allstring = !strcmp(all, "1");
    if (*range && !xlsx_parse_cellrange(range, &x.cell_range))
    {
        rc = 198;
        goto done;
    }
    if ((rc = xlsx_open_file(&x, path)))
    {
        if (rc == 601)
            rc = 603;
        goto done;
    }
    if ((rc = xlsx_parse_workbook(&x)))
        goto done;
    if (*sheet)
    {
        x.selected_sheet = xlsx_select_sheet_by_name(&x, sheet);
        if (x.selected_sheet < 0)
        {
            snprintf(x.error_message, sizeof(x.error_message), "worksheet %s not found", sheet);
            rc = 601;
            goto done;
        }
    }
    else if (x.num_sheets > 0 && x.sheets[x.selected_sheet].kind != XLSX_SHEET_WORKSHEET)
        goto done; /* only chart, dialog or macro sheets: native imports no data */
    if ((rc = xlsx_parse_shared_strings(&x)) || (rc = xlsx_parse_styles(&x)) ||
        (rc = xlsx_parse_worksheet(&x)))
        goto done;
    /* Like native, a cellrange() must lie within the worksheet's data area. */
    const XLSXCellRange *cr = &x.cell_range;
    const char *outside = cr->start_col > x.area_max_col ? "column range in option cellrange();\n    The upper left column is right of the data area in the worksheet."
                          : cr->end_col > x.area_max_col ? "column range in option cellrange();\n    The lower right column is right of the data area in the worksheet."
                          : cr->start_row > x.area_max_row ? "row range in option cellrange();\n    The upper left row is after the data area in the worksheet."
                          : cr->end_row > x.area_max_row ? "row range in option cellrange();\n    The lower right row is after the data area in the worksheet."
                                                         : NULL;
    if (outside)
    {
        snprintf(x.error_message, sizeof(x.error_message), "invalid %s", outside);
        rc = 198;
        goto done;
    }
    int r0 = x.cell_range.start_row < 0 ? x.min_row : x.cell_range.start_row,
        r1 = x.cell_range.end_row < 0 ? x.max_row : x.cell_range.end_row;
    int c0 = x.cell_range.start_col < 0 ? 0 : x.cell_range.start_col,
        c1 = x.cell_range.end_col < 0 ? x.max_col : x.cell_range.end_col;
    /* A used range starts at the first column holding a cell element (with
       or without a value), not necessarily A. */
    if (x.cell_range.start_col < 0)
        c0 = x.min_col == INT32_MAX ? c1 + 1 : x.min_col;
    if (r0 < 1 || r1 < r0 || c1 < c0)
        goto done;
    /* Fail before staging: too many variables (as native, r(900)), or more
     * cells than memory can hold. */
    int limit = atoi(maxvar);
    if (limit > 0 && c1 - c0 + 1 > limit)
    {
        char shown[40];
        if (limit >= 1000)
            snprintf(shown, sizeof(shown), "%d,%03d", limit / 1000, limit % 1000);
        else
            snprintf(shown, sizeof(shown), "%d", limit);
        snprintf(x.error_message, sizeof(x.error_message),
                 "no room to add more variables\n    Up to %s variables are currently allowed, although you could reset the\n    maximum using set maxvar; see help memory.",
                 shown);
        rc = 900;
        goto done;
    }
    if ((size_t)(r1 - r0 + 1) > x.max_cells / (size_t)(c1 - c0 + 1))
    {
        snprintf(x.error_message, sizeof(x.error_message),
                 "worksheet too large: %d rows by %d columns exceed the memory available", r1 - r0 + 1,
                 c1 - c0 + 1);
        rc = 909;
        goto done;
    }
    dataset->nvar = c1 - c0 + 1;
    dataset->cols = calloc((size_t)dataset->nvar, sizeof(cio_column));
    if (!dataset->cols)
    {
        rc = 909;
        goto done;
    }
    int start = r0 + x.firstrow;
    for (int j = 0; j < dataset->nvar; j++)
    {
        cio_column *c = &dataset->cols[j];
        cio_xls_column_name(c0 + j, c->name);
        cio_copy_text(c->excel_column, sizeof(c->excel_column), c->name);
        c->width = 1;
        cio_copy_text(c->format, sizeof(c->format), "%10.0g");
        if (x.firstrow)
        {
            double value;
            const char *text;
            int type = cio_xlsx_raw(&x, r0, c0 + j, &value, &text);
            char number[100];
            if (type == XLSX_CELL_DATE || type == XLSX_CELL_DATETIME)
            {
                /* A date-styled heading labels as its date text, as native. */
                cio_excel_date_text(cio_excel_days(value, x.date1904),
                                    cio_xlsx_native_format(&x, r0, c0 + j), number, sizeof(number));
                text = number;
            }
            else if (cio_xlsx_numeric(type))
            {
                snprintf(number, sizeof(number), "%.16g", value);
                text = number;
                c->label_numeric = type == XLSX_CELL_NUMBER;
                c->label_number = value;
            }
            if (*text)
            {
                cio_copy_text(c->label, sizeof(c->label), text);
                if (!cio_xlsx_numeric(type))
                    cio_copy_text(c->name, sizeof(c->name), text);
            }
        }
        cio_excel_column numeric = {0};
        for (int r = start; r <= r1; r++)
        {
            double value;
            const char *text;
            int type = cio_xlsx_raw(&x, r, c0 + j, &value, &text);
            if (*text && !cio_xlsx_numeric(type))
                c->is_string = 1;
            if (cio_xlsx_numeric(type))
                cio_excel_column_add(&numeric, cio_xlsx_native_format(&x, r, c0 + j),
                                     type == XLSX_CELL_BOOLEAN);
        }
        c->date_kind = numeric.kind;
        cio_copy_text(c->format, sizeof(c->format), cio_excel_column_format(&numeric));
        if (x.allstring)
            c->is_string = 1;
        /* Numbers in a string column take allstring()'s format, else native's
           %18.0g text: the ado formats them from the raw values. */
        c->format_numbers = c->is_string && (*allfmt || numeric.count || numeric.booleans);
    }
    if (r1 >= start)
    {
        if ((rc = cio_reserve_rows(dataset, (size_t)(r1 - start + 1))))
            goto done;
        dataset->nobs = (size_t)(r1 - start + 1);
    }
    for (int j = 0; j < dataset->nvar; j++)
        for (int r = start; r <= r1; r++)
        {
            cio_column *c = &dataset->cols[j];
            size_t i = (size_t)(r - start);
            double value;
            const char *text;
            int type = cio_xlsx_raw(&x, r, c0 + j, &value, &text);
            int date = type == XLSX_CELL_DATE ? 1 : type == XLSX_CELL_DATETIME ? 2 : 0;
            if (c->is_string)
            {
                char number[100];
                if (cio_xlsx_numeric(type))
                {
                    if (c->format_numbers && !date)
                        c->numbers[i] = value;
                    if (date)
                    {
                        /* Text keeps each cell's own display format */
                        cio_excel_date_text(cio_excel_days(value, x.date1904),
                                            cio_xlsx_native_format(&x, r, c0 + j), number,
                                            sizeof(number));
                    }
                    else
                        snprintf(number, sizeof(number), "%.16g", value);
                    /* Native general-format text omits zero before a decimal point. */
                    if (!date && number[0] == '0' && number[1] == '.')
                        memmove(number, number + 1, strlen(number));
                    else if (!date && number[0] == '-' && number[1] == '0' && number[2] == '.')
                        memmove(number + 1, number + 2, strlen(number + 1));
                    text = number;
                }
                size_t n = strlen(text);
                if (n > INT_MAX)
                {
                    rc = 109;
                    goto done;
                }
                c->strings[i] = strdup(text);
                if (!c->strings[i])
                {
                    rc = 909;
                    goto done;
                }
                if ((int)n > c->width)
                    c->width = (int)n;
            }
            else if (cio_xlsx_numeric(type) && !(c->date_kind && type == XLSX_CELL_BOOLEAN))
            {
                /* One unit per column: General serials in a date column convert too;
                   TRUE/FALSE in a date column stay missing, as in native */
                c->numbers[i] =
                    c->date_kind ? cio_excel_stata_value(value, x.date1904, c->date_kind) : value;
            }
        }
    for (int j = 0; j < dataset->nvar; j++)
    {
        /* Room for the %18.0g text that replaces each number */
        if (dataset->cols[j].format_numbers && !*allfmt && dataset->cols[j].width < 18)
            dataset->cols[j].width = 18;
        if (dataset->cols[j].is_string)
            cio_import_format(&dataset->cols[j], "xlsx");
    }
done:
    if (rc && *x.error_message)
        ctools_error("cimport", "%s", x.error_message);
    xlsx_context_free(&x);
    return rc;
}

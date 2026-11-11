#include "cio_internal.h"

/* Workbook inspection leaves the current Stata dataset untouched. Sheets
   without cells (and chart, dialog or macro sheets) report no range, as native
   does for empty worksheets. */
int cio_excel_sheet_info(int index, const char *name, int row0, int col0, int row1, int col1)
{
    char key[100], range[80] = "", from[16], to[16];
    if (row0 >= 1 && row1 >= row0 && col0 >= 0 && col1 >= col0)
    {
        cio_xls_column_name(col0, from);
        cio_xls_column_name(col1, to);
        snprintf(range, sizeof(range), "%s%d:%s%d", from, row0, to, row1);
    }
    snprintf(key, sizeof(key), "__cio_worksheet_%d", index);
    int rc = cio_save(key, name);
    if (rc)
        return rc;
    snprintf(key, sizeof(key), "__cio_range_%d", index);
    return cio_save(key, range);
}
int cio_excel_workbook(const char *kind)
{
    char path[4096];
    cio_local_text("__cio_filename", path, sizeof(path));
    int rc = 0;
    if (!strcmp(kind, "xls"))
    {
        xls_error_t error;
        xlsWorkBook *book = xls_open_file(path, "UTF-8", &error);
        if (!book)
            return 603;
        rc = cio_save_number("__cio_nworksheets", book->sheets.count);
        for (unsigned j = 0; !rc && j < book->sheets.count; j++)
        {
            if (!cio_xls_is_worksheet(book, j))
            {
                rc = cio_excel_sheet_info((int)j + 1, cio_xls_sheet_name(book, j), 0, 0, -1, -1);
                continue;
            }
            xlsWorkSheet *sheet = xls_getWorkSheet(book, (int)j);
            if (!sheet)
            {
                rc = 909;
                break;
            }
            if (xls_parseWorkSheet(sheet) != LIBXLS_OK)
            {
                xls_close_WS(sheet);
                rc = 610;
                break;
            }
            /* The same used range as the import, formatted blanks included */
            int r0, c0, r1, c1;
            cio_xls_used_range(sheet, &r0, &c0, &r1, &c1);
            rc = cio_excel_sheet_info((int)j + 1, cio_xls_sheet_name(book, j), r1 < 0 ? 1 : r0 + 1,
                                      c0, r1 + 1, c1);
            xls_close_WS(sheet);
        }
        xls_close_WB(book);
        return rc;
    }
    XLSXContext first;
    xlsx_context_init(&first);
    first.max_part_bytes = cio_xlsx_part_limit();
    if ((rc = xlsx_open_file(&first, path)))
    {
        xlsx_context_free(&first);
        return rc == 601 ? 603 : rc;
    }
    if ((rc = xlsx_parse_workbook(&first)))
    {
        xlsx_context_free(&first);
        return rc;
    }
    rc = cio_save_number("__cio_nworksheets", (size_t)first.num_sheets);
    for (int j = 0; !rc && j < first.num_sheets; j++)
    {
        if (first.sheets[j].kind != XLSX_SHEET_WORKSHEET)
        {
            rc = cio_excel_sheet_info(j + 1, first.sheets[j].name, 0, 0, -1, -1);
            continue;
        }
        XLSXContext sheet;
        xlsx_context_init(&sheet);
        sheet.max_cells = cio_xlsx_cell_limit();
        sheet.max_part_bytes = first.max_part_bytes;
        if (!(rc = xlsx_open_file(&sheet, path)) && !(rc = xlsx_parse_workbook(&sheet)))
        {
            sheet.selected_sheet = j;
            if (!(rc = xlsx_parse_shared_strings(&sheet)) && !(rc = xlsx_parse_styles(&sheet)) &&
                !(rc = xlsx_parse_worksheet(&sheet)))
            {
                /* The range starts at the first column with a cell element. */
                int c0 = sheet.min_col == INT32_MAX ? sheet.max_col + 1 : sheet.min_col;
                rc = cio_excel_sheet_info(j + 1, first.sheets[j].name, sheet.min_row, c0,
                                          sheet.max_row, sheet.max_col);
            }
        }
        if (rc && *sheet.error_message)
            ctools_error("cimport", "%s", sheet.error_message);
        xlsx_context_free(&sheet);
    }
    xlsx_context_free(&first);
    return rc;
}

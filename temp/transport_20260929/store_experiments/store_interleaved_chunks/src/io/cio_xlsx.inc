/* Reuse the streaming ZIP/XML parser, then transfer metadata and values through
   the same bounded interface as XLS and statistical formats. */
#include "cimport/cimport_xlsx.h"
/* Every stored cell lives in the parser's column-major arrays. */
static int xlsx_stored(XLSXContext *x, int row, int col) {
  return x->cm_active && row >= 1 && (size_t)row <= x->cm_num_rows &&
         col >= 0 && col < x->cm_num_cols;
}
static int xlsx_raw(XLSXContext *x, int row, int col, double *number,
                    const char **text) {
  *number = 0;
  *text = "";
  if (!xlsx_stored(x, row, col))
    return XLSX_CELL_EMPTY;
  int type = x->cm_types[col][row - 1];
  *number = x->cm_numeric[col][row - 1];
  if (type == XLSX_CELL_STRING) {
    uintptr_t pointer = 0;
    memcpy(&pointer, number, sizeof(pointer));
    *text = (const char *)pointer;
  }
  if (type == XLSX_CELL_SHARED_STRING) {
    int index = (int)*number;
    if (index < 0 || (uint32_t)index >= x->num_shared_strings)
      return XLSX_CELL_EMPTY;
    *text = x->shared_strings[index];
  }
  if (!*text)
    *text = "";
  return type;
}
static uint8_t xlsx_native_format(XLSXContext *x, int row, int col) {
  return xlsx_stored(x, row, col) ? x->cm_formats[col][row - 1] : 0;
}
static int xlsx_numeric(int type) {
  return type == XLSX_CELL_NUMBER || type == XLSX_CELL_DATE ||
         type == XLSX_CELL_DATETIME || type == XLSX_CELL_BOOLEAN;
}
static int scan_xlsx(const char *path) {
  XLSXContext x;
  xlsx_context_init(&x);
  int rc = 0;
  char sheet[1024], range[100], first[8], all[8], allfmt[80];
  local_text("__cio_sheet", sheet, sizeof(sheet));
  local_text("__cio_cellrange", range, sizeof(range));
  local_text("__cio_firstrow", first, sizeof(first));
  local_text("__cio_allstring", all, sizeof(all));
  local_text("__cio_allstringfmt", allfmt, sizeof(allfmt));
  x.firstrow = !strcmp(first, "1");
  x.allstring = !strcmp(all, "1");
  if (*range && !xlsx_parse_cellrange(range, &x.cell_range)) {
    rc = 198;
    goto done;
  }
  if ((rc = xlsx_open_file(&x, path))) {
    if (rc == 601)
      rc = 603;
    goto done;
  }
  if ((rc = xlsx_parse_workbook(&x)))
    goto done;
  if (*sheet) {
    x.selected_sheet = xlsx_select_sheet_by_name(&x, sheet);
    if (x.selected_sheet < 0) {
      snprintf(x.error_message, sizeof(x.error_message),
               "worksheet %s not found", sheet);
      rc = 601;
      goto done;
    }
  }
  if ((rc = xlsx_parse_shared_strings(&x)) || (rc = xlsx_parse_styles(&x)) ||
      (rc = xlsx_parse_worksheet(&x)))
    goto done;
  int r0 = x.cell_range.start_row < 0 ? x.min_row : x.cell_range.start_row,
      r1 = x.cell_range.end_row < 0 ? x.max_row : x.cell_range.end_row;
  int c0 = x.cell_range.start_col < 0 ? 0 : x.cell_range.start_col,
      c1 = x.cell_range.end_col < 0 ? x.max_col : x.cell_range.end_col;
  /* A used range starts at the first occupied column, not necessarily A. */
  if (x.cell_range.start_col < 0) {
    c0 = c1 + 1;
    for (int c = 0; c <= c1; c++)
      for (int r = r0; r <= r1; r++) {
        double value;
        const char *text;
        int type = xlsx_raw(&x, r, c, &value, &text);
        if (xlsx_numeric(type) || *text) {
          if (c < c0)
            c0 = c;
          break;
        }
      }
  }
  if (r0 < 1 || r1 < r0 || c1 < c0)
    goto done;
  cache.nvar = c1 - c0 + 1;
  cache.cols = calloc((size_t)cache.nvar, sizeof(cio_column));
  if (!cache.cols) {
    rc = 909;
    goto done;
  }
  int start = r0 + x.firstrow;
  for (int j = 0; j < cache.nvar; j++) {
    cio_column *c = &cache.cols[j];
    xls_column_name(c0 + j, c->name);
    copy_text(c->excel_column, sizeof(c->excel_column), c->name);
    c->width = 1;
    copy_text(c->format, sizeof(c->format), "%10.0g");
    if (x.firstrow) {
      double value;
      const char *text;
      int type = xlsx_raw(&x, r0, c0 + j, &value, &text);
      char number[80];
      if (xlsx_numeric(type)) {
        snprintf(number, sizeof(number), "%.16g", value);
        text = number;
      }
      if (*text) {
        copy_text(c->label, sizeof(c->label), text);
        if (!xlsx_numeric(type))
          copy_text(c->name, sizeof(c->name), text);
      }
    }
    cio_excel_column numeric = {0};
    for (int r = start; r <= r1; r++) {
      double value;
      const char *text;
      int type = xlsx_raw(&x, r, c0 + j, &value, &text);
      if (*text && !xlsx_numeric(type))
        c->is_string = 1;
      if (xlsx_numeric(type))
        cio_excel_column_add(&numeric, xlsx_native_format(&x, r, c0 + j),
                             type == XLSX_CELL_BOOLEAN);
    }
    c->date_kind = numeric.kind;
    copy_text(c->format, sizeof(c->format), cio_excel_column_format(&numeric));
    if (x.allstring)
      c->is_string = 1;
    c->format_numbers = c->is_string && *allfmt;
  }
  if (r1 >= start) {
    if ((rc = reserve_rows(&cache, (size_t)(r1 - start + 1))))
      goto done;
    cache.nobs = (size_t)(r1 - start + 1);
  }
  for (int j = 0; j < cache.nvar; j++)
    for (int r = start; r <= r1; r++) {
      cio_column *c = &cache.cols[j];
      size_t i = (size_t)(r - start);
      double value;
      const char *text;
      int type = xlsx_raw(&x, r, c0 + j, &value, &text);
      int date = type == XLSX_CELL_DATE       ? 1
                 : type == XLSX_CELL_DATETIME ? 2
                                              : 0;
      if (c->is_string) {
        char number[100];
        if (xlsx_numeric(type)) {
          if (c->format_numbers && !date)
            c->numbers[i] = value;
          if (date) {
            /* Text keeps each cell's own display format */
            excel_date_text(cio_excel_days(value, x.date1904),
                            xlsx_native_format(&x, r, c0 + j), number,
                            sizeof(number));
          } else
            snprintf(number, sizeof(number), "%.16g", value);
          /* Native general-format text omits zero before a decimal point. */
          if (!date && number[0] == '0' && number[1] == '.')
            memmove(number, number + 1, strlen(number));
          else if (!date && number[0] == '-' && number[1] == '0' &&
                   number[2] == '.')
            memmove(number + 1, number + 2, strlen(number + 1));
          text = number;
        }
        size_t n = strlen(text);
        if (n > INT_MAX) {
          rc = 109;
          goto done;
        }
        c->strings[i] = strdup(text);
        if (!c->strings[i]) {
          rc = 909;
          goto done;
        }
        if ((int)n > c->width)
          c->width = (int)n;
      } else if (xlsx_numeric(type)) {
        /* One unit per column: General serials in a date column convert too */
        c->numbers[i] = c->date_kind ? cio_excel_stata_value(value, x.date1904,
                                                             c->date_kind)
                                     : value;
      }
    }
  for (int j = 0; j < cache.nvar; j++)
    if (cache.cols[j].is_string)
      import_format(&cache.cols[j], "xlsx");
done:
  if (rc && *x.error_message)
    ctools_error("cimport", "%s", x.error_message);
  xlsx_context_free(&x);
  return rc;
}

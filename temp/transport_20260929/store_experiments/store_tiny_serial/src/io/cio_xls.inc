/* Legacy BIFF Excel input. libxls handles OLE chains, SST continuation, cached
   formula results and both BIFF5/BIFF8. Shared transfer logic stays in cio.c.
 */
#include "vendor/libxls/include/xls.h"
/* Cell records present in the used range (NUMBER, RK, MULRK, FORMULA, text). */
static int xls_nonempty(const xlsCell *c) {
  return c &&
         (c->id == 0x203 || c->id == 0x27e || c->id == 0xbd || c->id == 0x6 ||
          c->id == 0x206 || (c->str && *c->str && c->id != 0x201));
}
/* libxls reports a non-numeric formula result only as text ("bool", "error",
   "" or the string), which cannot be told from a string result with the
   same text. Record each result type from the FORMULA record as it parses. */
typedef struct {
  XLSWORD row, col;
  XLSBYTE type; /* 0 text, 1 boolean, 2 error, 3 empty text */
} xls_formula_result;
static struct {
  xls_formula_result *items;
  size_t count, capacity;
  int failed;
} xls_results;
static void xls_formula_record(XLSWORD id, XLSWORD len, XLSBYTE *record) {
  if ((id != XLS_RECORD_FORMULA && id != XLS_RECORD_FORMULA_ALT) ||
      len < sizeof(FORMULA))
    return;
  const FORMULA *f = (const FORMULA *)record;
  if (f->res != 0xffff)
    return;
  if (xls_results.count == xls_results.capacity) {
    size_t capacity = xls_results.capacity ? 2 * xls_results.capacity : 256;
    xls_formula_result *items =
        cio_realloc(xls_results.items, capacity, sizeof(*items));
    if (!items) {
      xls_results.failed = 1;
      return;
    }
    xls_results.items = items;
    xls_results.capacity = capacity;
  }
  xls_formula_result *r = &xls_results.items[xls_results.count++];
  r->row = f->row;
  r->col = f->col;
  r->type = f->resid;
}
static int xls_result_order(const void *a, const void *b) {
  const xls_formula_result *x = a, *y = b;
  if (x->row != y->row)
    return x->row < y->row ? -1 : 1;
  return x->col < y->col ? -1 : x->col > y->col;
}
static void xls_results_reset(void) {
  free(xls_results.items);
  memset(&xls_results, 0, sizeof(xls_results));
}
enum { XLS_EMPTY, XLS_NUMBER, XLS_BOOLEAN, XLS_STRING };
/* Value class of a cell, as native import excel reads it: text formula
   results are strings, TRUE/FALSE are 1/0 (cell->d), and error or empty
   results are missing. */
static int xls_class(const xlsCell *c, const char **text) {
  *text = "";
  if (!c || c->id == 0x201)
    return XLS_EMPTY;
  if (c->id == 0x203 || c->id == 0x27e || c->id == 0xbd || c->id == 0x206)
    return XLS_NUMBER;
  if (c->id == 0x205) /* BOOLERR */
    return c->str && !strcmp(c->str, "bool") ? XLS_BOOLEAN : XLS_EMPTY;
  if (c->id == 0x6 || c->id == 0x406) {
    if (c->l != 0xffff)
      return XLS_NUMBER;
    xls_formula_result key = {c->row, c->col, 0}, *found = NULL;
    if (xls_results.count)
      found = bsearch(&key, xls_results.items, xls_results.count, sizeof(key),
                      xls_result_order);
    if (!found || found->type > 1)
      return XLS_EMPTY;
    if (found->type == 1)
      return XLS_BOOLEAN;
  }
  *text = c->str ? c->str : "";
  return **text ? XLS_STRING : XLS_EMPTY;
}
static int xls_ref(const char *s, int *row, int *col, const char **end) {
  int c = 0, r = 0;
  const char *p = s;
  while (isalpha((unsigned char)*p)) {
    c = c * 26 + toupper((unsigned char)*p) - 'A' + 1;
    if (c > 256)
      return 198;
    p++;
  }
  if (!c || !isdigit((unsigned char)*p))
    return 198;
  while (isdigit((unsigned char)*p)) {
    r = r * 10 + *p - '0';
    if (r > 65536)
      return 198;
    p++;
  }
  if (!r)
    return 198;
  *row = r - 1;
  *col = c - 1;
  *end = p;
  return 0;
}
static void xls_column_name(int index, char *name) {
  char rev[12];
  int n = 0;
  index++;
  do {
    rev[n++] = (char)('A' + (index - 1) % 26);
    index = (index - 1) / 26;
  } while (index);
  for (int i = 0; i < n; i++)
    name[i] = rev[n - i - 1];
  name[n] = 0;
}
static uint8_t xls_native_format(xlsWorkBook *w, const xlsCell *cell) {
  if (!cell || cell->xf >= w->xfs.count)
    return 0;
  int id = w->xfs.xf[cell->xf].format;
  uint8_t code = cio_excel_format_code(id, NULL);
  if (code)
    return code;
  for (unsigned j = 0; j < w->formats.count; j++)
    if (w->formats.format[j].index == id)
      return cio_excel_format_code(id, w->formats.format[j].value);
  return 0;
}
static int scan_xls(const char *path) {
  xls_error_t err;
  int rc = 0;
  char sheet[1024], range[100], first[8], allstring[8], allfmt[80];
  local_text("__cio_sheet", sheet, sizeof(sheet));
  local_text("__cio_cellrange", range, sizeof(range));
  local_text("__cio_firstrow", first, sizeof(first));
  local_text("__cio_allstring", allstring, sizeof(allstring));
  local_text("__cio_allstringfmt", allfmt, sizeof(allfmt));
  xlsWorkBook *w = xls_open_file(path, "UTF-8", &err);
  if (!w) {
    ctools_error("cimport", "%s", xls_getError(err));
    return 603;
  }
  xlsWorkSheet *ws = NULL;
  int which = 0;
  if (*sheet) {
    which = -1;
    for (unsigned j = 0; j < w->sheets.count; j++)
      if (!strcmp(sheet, w->sheets.sheet[j].name)) {
        which = (int)j;
        break;
      }
    if (which < 0) {
      rc = 601;
      goto done;
    }
  }
  if (!w->sheets.count) {
    rc = 610;
    goto done;
  }
  ws = xls_getWorkSheet(w, which);
  if (!ws) {
    rc = 909;
    goto done;
  }
  xls_results_reset();
  xls_set_formula_hander(xls_formula_record);
  err = xls_parseWorkSheet(ws);
  xls_set_formula_hander(NULL);
  if (xls_results.failed) {
    rc = 909;
    goto done;
  }
  if (err != LIBXLS_OK) {
    rc = 610;
    goto done;
  }
  if (xls_results.count)
    qsort(xls_results.items, xls_results.count, sizeof(xls_formula_result),
          xls_result_order);
  int r0 = 65536, c0 = 256, r1 = -1, c1 = -1;
  for (int r = 0; r <= ws->rows.lastrow; r++)
    for (int c = 0; c <= ws->rows.lastcol; c++)
      if (xls_nonempty(xls_cell(ws, (XLSWORD)r, (XLSWORD)c))) {
        if (r < r0)
          r0 = r;
        if (r > r1)
          r1 = r;
        if (c < c0)
          c0 = c;
        if (c > c1)
          c1 = c;
      }
  /* Explicit blank cells and ROW records also belong to the used range. */
  if (r1 >= 0) {
    r1 = ws->rows.lastrow;
    c1 = ws->rows.lastcol;
  }
  if (*range) {
    const char *end;
    end = range;
    if (*end != ':' && (rc = xls_ref(range, &r0, &c0, &end)))
      goto done;
    if (*end == ':') {
      end++;
      if (*end && (rc = xls_ref(end, &r1, &c1, &end)))
        goto done;
    }
    if (*end || r1 < r0 || c1 < c0) {
      rc = 198;
      goto done;
    }
  }
  if (r1 < r0 || c1 < c0)
    goto done;
  cache.nvar = c1 - c0 + 1;
  cache.cols = calloc((size_t)cache.nvar, sizeof(cio_column));
  if (!cache.cols) {
    rc = 909;
    goto done;
  }
  int start = r0 + (!strcmp(first, "1"));
  for (int j = 0; j < cache.nvar; j++) {
    cio_column *c = &cache.cols[j];
    xls_column_name(c0 + j, c->name);
    copy_text(c->excel_column, sizeof(c->excel_column), c->name);
    c->width = 1;
    copy_text(c->format, sizeof(c->format), "%10.0g");
    if (!strcmp(first, "1")) {
      xlsCell *h = xls_cell(ws, (XLSWORD)r0, (XLSWORD)(c0 + j));
      const char *text;
      int kind = xls_class(h, &text);
      if (kind == XLS_STRING) {
        copy_text(c->name, sizeof(c->name), text);
        copy_text(c->label, sizeof(c->label), text);
      } else if (kind == XLS_BOOLEAN) {
        copy_text(c->label, sizeof(c->label), h->d != 0 ? "1" : "0");
      } else if (kind == XLS_NUMBER && h->str && *h->str)
        copy_text(c->label, sizeof(c->label), h->str);
    }
    cio_excel_column numeric = {0};
    for (int r = start; r <= r1; r++) {
      xlsCell *cell = xls_cell(ws, (XLSWORD)r, (XLSWORD)(c0 + j));
      const char *text;
      int kind = xls_class(cell, &text);
      if (kind == XLS_STRING)
        c->is_string = 1;
      else if (kind != XLS_EMPTY)
        cio_excel_column_add(&numeric, xls_native_format(w, cell),
                             kind == XLS_BOOLEAN);
    }
    c->date_kind = numeric.kind;
    copy_text(c->format, sizeof(c->format), cio_excel_column_format(&numeric));
    if (!strcmp(allstring, "1"))
      c->is_string = 1;
    c->format_numbers = c->is_string && *allfmt;
  }
  if (r1 >= start) {
    rc = reserve_rows(&cache, (size_t)(r1 - start + 1));
    if (rc)
      goto done;
    cache.nobs = (size_t)(r1 - start + 1);
  }
  for (int r = start; r <= r1; r++)
    for (int j = 0; j < cache.nvar; j++) {
      cio_column *c = &cache.cols[j];
      xlsCell *cell = xls_cell(ws, (XLSWORD)r, (XLSWORD)(c0 + j));
      size_t i = (size_t)(r - start);
      const char *text;
      int kind = xls_class(cell, &text);
      if (c->is_string) {
        const char *s = text;
        char numeric[80];
        if (kind == XLS_NUMBER || kind == XLS_BOOLEAN) {
          /* Text keeps each cell's own display format */
          uint8_t code = xls_native_format(w, cell);
          int date = kind == XLS_NUMBER ? cio_excel_date_kind(code) : 0;
          if (c->format_numbers && !date)
            c->numbers[i] = cell->d;
          if (date)
            excel_date_text(cio_excel_days(cell->d, w->is1904), code, numeric,
                            sizeof(numeric));
          else
            snprintf(numeric, sizeof(numeric), "%.16g", cell->d);
          if (!date && numeric[0] == '0' && numeric[1] == '.')
            memmove(numeric, numeric + 1, strlen(numeric));
          else if (!date && numeric[0] == '-' && numeric[1] == '0' &&
                   numeric[2] == '.')
            memmove(numeric + 1, numeric + 2, strlen(numeric + 1));
          s = numeric;
        }
        size_t n = strlen(s);
        if (n > INT_MAX) {
          rc = 109;
          goto done;
        }
        c->strings[i] = strdup(s);
        if (!c->strings[i]) {
          rc = 909;
          goto done;
        }
        if ((int)n > c->width)
          c->width = (int)n;
      } else if (kind == XLS_NUMBER || kind == XLS_BOOLEAN) {
        /* One unit per column: General serials in a date column convert too */
        c->numbers[i] = c->date_kind ? cio_excel_stata_value(cell->d, w->is1904,
                                                             c->date_kind)
                                     : cell->d;
      }
    }
  for (int j = 0; j < cache.nvar; j++)
    if (cache.cols[j].is_string)
      import_format(&cache.cols[j], "xls");
done:
  xls_set_formula_hander(NULL);
  xls_results_reset();
  if (ws)
    xls_close_WS(ws);
  xls_close_WB(w);
  return rc;
}

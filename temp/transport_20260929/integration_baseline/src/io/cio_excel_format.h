/* Native Excel display-format mappings shared by the XLS and XLSX readers. */
#ifndef CTOOLS_CIO_EXCEL_FORMAT_H
#define CTOOLS_CIO_EXCEL_FORMAT_H
#include <ctype.h>
#include <stdint.h>
#include <string.h>
#define CIO_EXCEL_DATE 90
#define CIO_EXCEL_DATETIME 91
static inline uint8_t cio_excel_format_code(int id, const char *custom) {
  if ((id >= 14 && id <= 22) || (id >= 45 && id <= 47))
    return (uint8_t)id;
  /* Built-in number formats that native maps to a Stata display format. As
     for dates, the built-in id decides even if the workbook redefines it;
     custom number formats display as %10.0g. */
  switch (id) {
  case 2:
  case 3:
  case 4:
  case 7:
  case 8:
  case 9:
  case 10:
  case 11:
  case 37:
  case 38:
  case 39:
  case 40:
  case 48:
    return (uint8_t)id;
  }
  if (!custom || !*custom)
    return 0;
  char lower[256];
  size_t n = strlen(custom);
  if (n >= sizeof(lower))
    n = sizeof(lower) - 1;
  for (size_t j = 0; j < n; j++)
    lower[j] = (char)tolower((unsigned char)custom[j]);
  lower[n] = 0;
  if (!strcmp(lower, "h:mm:ss") || !strcmp(lower, "hh:mm:ss"))
    return 21;
  if (!strcmp(lower, "h:mm") || !strcmp(lower, "hh:mm"))
    return 20;
  if (!strcmp(lower, "h:mm am/pm") || !strcmp(lower, "hh:mm am/pm"))
    return 18;
  if (!strcmp(lower, "h:mm:ss am/pm") || !strcmp(lower, "hh:mm:ss am/pm"))
    return 19;
  if (!strcmp(lower, "m/d/yy h:mm") || !strcmp(lower, "m/d/yyyy h:mm"))
    return 22;
  uint8_t kind = 0;
  int quoted = 0;
  for (const char *p = lower; *p; p++) {
    if (*p == '"') {
      quoted = !quoted;
      continue;
    }
    if (quoted)
      continue;
    if (*p == '\\' || *p == '_' || *p == '*') {
      if (p[1])
        p++;
      continue;
    }
    if (*p == '[') {
      if (p[1] == 'h' || p[1] == 's')
        kind = CIO_EXCEL_DATETIME;
      while (p[1] && *p != ']')
        p++;
      continue;
    }
    if (*p == 'h' || *p == 's')
      kind = CIO_EXCEL_DATETIME;
    else if (!kind && (*p == 'd' || *p == 'm' || *p == 'y'))
      kind = CIO_EXCEL_DATE;
  }
  return kind;
}
static inline uint8_t cio_excel_date_kind(uint8_t code) {
  return code == CIO_EXCEL_DATE || (code >= 14 && code <= 17) ? 1
         : code == CIO_EXCEL_DATETIME || (code >= 18 && code <= 22) ||
                 (code >= 45 && code <= 47)
             ? 2
             : 0;
}
static inline const char *cio_excel_stata_format(uint8_t code) {
  switch (code) {
  case 2:
  case 7:
  case 8:
    return "%14.2f";
  case 3:
  case 37:
  case 38:
    return "%10.0gc";
  case 4:
  case 39:
  case 40:
    return "%14.2fc";
  case 9:
    return "%4.2f";
  case 10:
    return "%6.4f";
  case 11:
    return "%10.2e";
  case 48:
    return "%10.1e";
  case 14:
    return "%tdnn/dd/CCYY";
  case 15:
    return "%tddd-Mon-YY";
  case 16:
    return "%tddd-Mon";
  case 17:
    return "%tdMon-YY";
  case 18:
    return "%tchh:MM_AM";
  case 19:
    return "%tchh:MM:SS_AM";
  case 20:
    return "%tchH:MM";
  case 21:
  case 46:
    return "%tcHH:MM:SS";
  case 22:
    return "%tcnn/dd/ccYY_hh:MM";
  case 45:
  case 47:
    return "%tcMM:SS";
  case CIO_EXCEL_DATE:
    return "%td";
  case CIO_EXCEL_DATETIME:
    return "%tc";
  default:
    return "%10.0g";
  }
}
/* Excel serial to days since 01jan1960. In the 1900 system the fictitious
   29feb1900 (serial 60) becomes 28feb1900, as in native import excel. */
static inline double cio_excel_days(double serial, int date1904) {
  return serial - (date1904 ? 20454 : (serial < 60 ? 21915 : 21916));
}
/* Stata value of a numeric cell in a date (1) or datetime (2) column. Every
   numeric cell in such a column converts, including General-format serials;
   native rounds through seconds before milliseconds. */
static inline double cio_excel_stata_value(double serial, int date1904,
                                           int kind) {
  double milliseconds = cio_excel_days(serial, date1904) * 86400 * 1000;
  return kind == 2 ? milliseconds : milliseconds / 86400000;
}
/* One display decision per numeric column, as native import excel makes it:
   any date or time cell makes the column a date (%td) if all such cells are
   dates only, else a datetime (%tc). The cells' own Stata format is used only
   when every number shares it; otherwise plain %td, %tc or %10.0g. TRUE/FALSE
   cells take no part in that decision (their style is ignored, and in a date
   column they import as missing); a column of only TRUE/FALSE cells displays
   as %1.0f. */
typedef struct {
  int count, booleans, mixed; /* count: numbers other than TRUE/FALSE */
  uint8_t code, kind;
} cio_excel_column;
static inline void cio_excel_column_add(cio_excel_column *c, uint8_t code,
                                        int boolean) {
  if (boolean) {
    c->booleans++;
    return;
  }
  uint8_t kind = cio_excel_date_kind(code);
  if (!c->count)
    c->code = code;
  else if (code != c->code)
    c->mixed = 1;
  if (kind > c->kind)
    c->kind = kind;
  c->count++;
}
static inline const char *cio_excel_column_format(const cio_excel_column *c) {
  if (!c->count)
    return c->booleans ? "%1.0f" : "%10.0g";
  if (c->kind)
    return cio_excel_stata_format(!c->mixed      ? c->code
                                  : c->kind == 2 ? CIO_EXCEL_DATETIME
                                                 : CIO_EXCEL_DATE);
  return cio_excel_stata_format(c->mixed ? 0 : c->code);
}
#endif

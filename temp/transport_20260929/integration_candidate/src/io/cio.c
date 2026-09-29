/* Binary statistical formats. All parsing and serialization stays in C.
   The ado layer supplies dataset metadata and creates Stata variables because
   the plugin API cannot create variables or access their names and labels. */
#include "cio.h"
#include "cexport/cexport_parse.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "ctools_types.h"
#include "vendor/readstat/readstat.h"
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
  char name[301], label[1025], source_format[257], format[80], labelset[301],
      excel_column[16];
  int is_string, width, date_kind, binary, force_double, format_numbers;
  size_t *lengths;
  double *numbers;
  char **strings;
} cio_column;
typedef struct {
  char *set, *text;
  double value;
} cio_label;
typedef struct {
  cio_column *cols;
  int nvar;
  size_t nobs, capacity;
  cio_label *labels;
  size_t nlabels;
  char datalabel[1025], kind[24];
  int error;
  int shp_blob, shp_z, shp_m, shp_order;
  int raw_storage;
} cio_dataset;
static cio_dataset cache;
static void *cio_realloc(void *p, size_t n, size_t width) {
  if (width && n > SIZE_MAX / width)
    return NULL;
  return realloc(p, n * width);
}

static void dataset_free(cio_dataset *d) {
  for (int j = 0; d->cols && j < d->nvar; j++) {
    cio_column *c = &d->cols[j];
    if (c->strings)
      for (size_t i = 0; i < d->nobs; i++)
        free(c->strings[i]);
    free(c->strings);
    free(c->numbers);
    free(c->lengths);
  }
  for (size_t i = 0; i < d->nlabels; i++) {
    free(d->labels[i].set);
    free(d->labels[i].text);
  }
  free(d->cols);
  free(d->labels);
  memset(d, 0, sizeof(*d));
}
void cio_cleanup(void) { dataset_free(&cache); }
static void copy_text(char *out, size_t cap, const char *s) {
  snprintf(out, cap, "%s", s ? s : "");
}
static int local_text(const char *key, char *out, int size) {
  char name[128];
  snprintf(name, sizeof(name), "_%s", key);
  out[0] = 0;
  return SF_macro_use(name, out, size - 1);
}
/* Persistent metadata must use actual global names: underscore-prefixed SPI
   macro names belong to the calling ado level and disappear in helper programs.
 */
static int cio_save(const char *key, const char *value) {
  char name[128];
  if (!strncmp(key, "__cio_", 6))
    snprintf(name, sizeof(name), "CTOOLS_CIO_%s", key + 6);
  else
    copy_text(name, sizeof(name), key);
  return SF_macro_save(name, (char *)value);
}
static int save_number(const char *key, size_t v) {
  char s[40];
  snprintf(s, sizeof(s), "%zu", v);
  return cio_save((char *)key, s);
}
static double tagged_missing(char tag) {
  int n = tag >= 'A' && tag <= 'Z'   ? tag - 'A' + 1
          : tag >= 'a' && tag <= 'z' ? tag - 'a' + 1
                                     : 0;
  uint64_t bits;
  double v = SV_missval;
  memcpy(&bits, &v, 8);
  bits += (uint64_t)n << 40;
  memcpy(&v, &bits, 8);
  return v;
}
static double number_value(readstat_value_t v) {
  if (readstat_value_is_tagged_missing(v))
    return tagged_missing(readstat_value_tag(v));
  if (readstat_value_is_system_missing(v))
    return SV_missval;
  switch (readstat_value_type(v)) {
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
static int reserve_rows(cio_dataset *d, size_t rows) {
  if (rows <= d->capacity)
    return 0;
  if (rows > INT_MAX)
    return 920;
  size_t cap = d->capacity ? d->capacity : 256;
  while (cap < rows)
    cap = cap > INT_MAX / 2 ? (size_t)INT_MAX : cap * 2;
  for (int j = 0; j < d->nvar; j++) {
    cio_column *c = &d->cols[j];
    if (c->is_string) {
      char **p = cio_realloc(c->strings, cap, sizeof(char *));
      if (!p)
        return 909;
      if (c->is_string) {
        size_t *lengths = cio_realloc(c->lengths, cap, sizeof(size_t));
        if (!lengths) {
          c->strings = p;
          return 909;
        }
        c->lengths = lengths;
        memset(lengths + d->capacity, 0, (cap - d->capacity) * sizeof(size_t));
      }
      c->strings = p;
      memset(p + d->capacity, 0, (cap - d->capacity) * sizeof(char *));
    }
    if (!c->is_string || c->format_numbers) {
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
static int metadata_cb(readstat_metadata_t *m, void *ctx) {
  cio_dataset *d = ctx;
  d->nvar = readstat_get_var_count(m);
  if (d->nvar < 0 || d->nvar > 32767) {
    d->error = 103;
    return READSTAT_HANDLER_ABORT;
  }
  d->cols = calloc((size_t)d->nvar, sizeof(cio_column));
  if (!d->cols && d->nvar) {
    d->error = 909;
    return READSTAT_HANDLER_ABORT;
  }
  copy_text(d->datalabel, sizeof(d->datalabel), readstat_get_file_label(m));
  return READSTAT_HANDLER_OK;
}
/* SAS epoch is the same as Stata's. SPSS seconds start at 14 October 1582.
   Date kinds: 1 = day, 2 = datetime, 3 = time in milliseconds. */
static void import_format(cio_column *c, const char *kind) {
  const char *f = c->source_format;
  if (c->is_string) {
    snprintf(c->format, sizeof(c->format), "%%%ds",
             c->width > 2045 ? 9
             : c->width > 9  ? c->width
                             : 9);
    return;
  }
  copy_text(c->format, sizeof(c->format), "%10.0g");
  if (!strcmp(kind, "spss")) {
    char code[40];
    int width = 0, decimals = 0;
    size_t k = 0;
    while (isalpha((unsigned char)f[k]) && k < sizeof(code) - 1) {
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
    if (datefmt) {
      c->date_kind = 2;
      copy_text(c->format, sizeof(c->format), datefmt);
    } else if (!strcmp(code, "DATETIME")) {
      c->date_kind = 2;
      copy_text(c->format, sizeof(c->format), "%tcdd-Mon-CCYY_HH:MM:SS.ss");
    } else if (!strcmp(code, "TIME")) {
      /* Native converts SPSS TIME using the calendar epoch, but its
         invalid generated display format leaves the numeric default. */
      c->date_kind = 2;
    } else if (!strcmp(code, "DTIME")) {
      c->date_kind = 3;
      copy_text(c->format, sizeof(c->format), "%tcjjj_HH:MM:SS");
    } else if (!strcmp(code, "COMMA"))
      snprintf(c->format, sizeof(c->format), "%%%d.%dfc", width, decimals);
    else if (!strcmp(code, "E"))
      snprintf(c->format, sizeof(c->format), "%%%d.%de", width, decimals);
    /* Native assigns Z's string display format to numeric storage, which
       Stata's public metadata interface rejects. Retain the numeric default. */
    else if (!strcmp(code, "F") || !strcmp(code, "DOT") || !strcmp(code, "PCT"))
      snprintf(c->format, sizeof(c->format), "%%%d.%df", width, decimals);
  } else if (!strcmp(kind, "sas")) {
    /* SAS7BDAT stores format names separately from widths. Prefix matching
       misclassifies user formats such as DATEAMPM and DATE9. */
    const char *datefmt = !strcmp(f, "DATE")     ? "%td"
                          : !strcmp(f, "YYMMDD") ? "%tdYYNNDD"
                          : !strcmp(f, "MMDDYY") ? "%tdNN/DD/YY"
                          : !strcmp(f, "DDMMYY") ? "%tdDD/NN/YY"
                          : !strcmp(f, "MONYY")  ? "%tdMonthCCYY"
                          : !strcmp(f, "YEAR")   ? "%tdCCYY"
                          : !strcmp(f, "WEEKDATE")
                              ? "%td_DAYNAME,_Month_DD,_CCYY"
                              : NULL;
    if (datefmt) {
      c->date_kind = 1;
      copy_text(c->format, sizeof(c->format), datefmt);
    } else if (!strcmp(f, "JULIAN")) {
      copy_text(c->format, sizeof(c->format), "%8.0f");
    } else if (!strcmp(f, "DATETIME")) {
      c->date_kind = 2;
      copy_text(c->format, sizeof(c->format), "%tc");
    } else if (!strcmp(f, "TIME") || !strcmp(f, "TOD") || !strcmp(f, "HHMM")) {
      c->date_kind = 3;
      copy_text(c->format, sizeof(c->format),
                !strcmp(f, "HHMM") ? "%tcHH:MM" : "%tcHH:MM:SS");
    }
  } else if (!strncmp(f, "DATETIME", 8) || !strncmp(f, "E8601DT", 7)) {
    c->date_kind = 2;
    copy_text(c->format, sizeof(c->format), "%tc");
  } else if (!strncmp(f, "DATE", 4) || !strncmp(f, "YYMMDD", 6) ||
             !strncmp(f, "MMDDYY", 6) || !strncmp(f, "DDMMYY", 6) ||
             !strncmp(f, "JULIAN", 6)) {
    c->date_kind = 1;
    copy_text(c->format, sizeof(c->format), "%td");
  } else if (!strncmp(f, "TIME", 4)) {
    c->date_kind = 3;
    copy_text(c->format, sizeof(c->format), "%tcHH:MM:SS");
  }
}
static int variable_cb(int index, readstat_variable_t *v, const char *set,
                       void *ctx) {
  cio_dataset *d = ctx;
  if (index < 0 || index >= d->nvar) {
    d->error = 610;
    return READSTAT_HANDLER_ABORT;
  }
  cio_column *c = &d->cols[index];
  copy_text(c->name, sizeof(c->name), readstat_variable_get_name(v));
  copy_text(c->label, sizeof(c->label), readstat_variable_get_label(v));
  copy_text(c->source_format, sizeof(c->source_format),
            readstat_variable_get_format(v));
  copy_text(c->labelset, sizeof(c->labelset), set);
  c->is_string =
      readstat_variable_get_type_class(v) == READSTAT_TYPE_CLASS_STRING;
  c->width = 1;
  if ((!strcmp(d->kind, "sasxport5") || d->raw_storage) && c->is_string)
    c->width = (int)readstat_variable_get_storage_width(v);
  if (d->raw_storage && !c->is_string)
    c->force_double = 1;
  import_format(c, d->kind);
  if (!strcmp(d->kind, "sasxport5")) {
    copy_text(c->labelset, sizeof(c->labelset), c->source_format);
    if (c->date_kind || !strncmp(c->source_format, "BEST", 4))
      c->labelset[0] = 0;
    if (c->date_kind == 1)
      copy_text(c->format, sizeof(c->format), "%d");
    else if (!c->is_string)
      copy_text(c->format, sizeof(c->format), "%10.0g");
    c->date_kind = 0;
  }
  return READSTAT_HANDLER_OK;
}
static int value_cb(int row, readstat_variable_t *v, readstat_value_t value,
                    void *ctx) {
  cio_dataset *d = ctx;
  int j = readstat_variable_get_index(v);
  if (row < 0 || j < 0 || j >= d->nvar) {
    d->error = 610;
    return READSTAT_HANDLER_ABORT;
  }
  int rc = reserve_rows(d, (size_t)row + 1);
  if (rc) {
    d->error = rc;
    return READSTAT_HANDLER_ABORT;
  }
  if ((size_t)row >= d->nobs)
    d->nobs = (size_t)row + 1;
  cio_column *c = &d->cols[j];
  if (c->is_string) {
    const char *s = readstat_string_value(value);
    if (!strcmp(d->kind, "spss") && readstat_value_is_defined_missing(value, v))
      s = "";
    if (!s)
      s = "";
    size_t len = strlen(s);
    if (len > INT_MAX) {
      d->error = 109;
      return READSTAT_HANDLER_ABORT;
    }
    char *p = malloc(len + 1);
    if (!p) {
      d->error = 909;
      return READSTAT_HANDLER_ABORT;
    }
    memcpy(p, s, len + 1);
    free(c->strings[row]);
    c->strings[row] = p;
    c->lengths[row] = len;
    if ((int)len > c->width)
      c->width = (int)len;
  } else {
    double x = number_value(value);
    if (!strcmp(d->kind, "spss") && readstat_value_is_defined_missing(value, v))
      x = SV_missval;
    /* Native XPORT import collapses SAS special missings. */
    if (strcmp(d->kind, "sasxport5") && SF_is_missing(x))
      x = SV_missval;
    if (!SF_is_missing(x) && c->date_kind) {
      if (!strcmp(d->kind, "spss"))
        x = c->date_kind == 3 ? x * 1000 : x * 1000 - 11903760000000.0;
      else if (!strncmp(d->kind, "sasxport",
                        8)) { /* Native transport values retain their original
                                 units. */
      } else if (c->date_kind >= 2)
        x *= 1000;
    }
    c->numbers[row] = x;
  }
  return READSTAT_HANDLER_OK;
}
static int label_cb(const char *set, readstat_value_t value, const char *text,
                    void *ctx) {
  cio_dataset *d = ctx;
  if (readstat_value_type_class(value) == READSTAT_TYPE_CLASS_STRING)
    return READSTAT_HANDLER_OK;
  cio_label *p = cio_realloc(d->labels, d->nlabels + 1, sizeof(cio_label));
  if (!p) {
    d->error = 909;
    return READSTAT_HANDLER_ABORT;
  }
  d->labels = p;
  cio_label *l = &p[d->nlabels];
  memset(l, 0, sizeof(*l));
  l->set = strdup(set ? set : "");
  l->text = strdup(text ? text : "");
  l->value = number_value(value);
  d->nlabels++;
  if (!l->set || !l->text) {
    d->error = 909;
    return READSTAT_HANDLER_ABORT;
  }
  return READSTAT_HANDLER_OK;
}
static void error_cb(const char *msg, void *ctx) {
  (void)ctx;
  ctools_error("cio", "%s", msg);
}
static int progress_cb(double progress, void *ctx) {
  (void)progress;
  (void)ctx;
  return SF_poll() ? READSTAT_HANDLER_ABORT : READSTAT_HANDLER_OK;
}
/* dBase III/IV: fixed records, little-endian directory, ASCII field values. */
static unsigned dbf_u16(const unsigned char *p) {
  return (unsigned)p[0] | (unsigned)p[1] << 8;
}
static uint32_t dbf_u32(const unsigned char *p) {
  return (uint32_t)p[0] | (uint32_t)p[1] << 8 | (uint32_t)p[2] << 16 |
         (uint32_t)p[3] << 24;
}
static void dbf_p16(unsigned char *p, unsigned v) {
  p[0] = (unsigned char)v;
  p[1] = (unsigned char)(v >> 8);
}
static void dbf_p32(unsigned char *p, uint32_t v) {
  for (int i = 0; i < 4; i++)
    p[i] = (unsigned char)(v >> (8 * i));
}
/* Gregorian civil calendar, independent of locale, timezone, and time_t. */
static int civil_days(int y, unsigned m, unsigned d) {
  y -= m <= 2;
  int era = (y >= 0 ? y : y - 399) / 400;
  unsigned yy = (unsigned)(y - era * 400);
  unsigned doy = (153 * (m > 2 ? m - 3 : m + 9) + 2) / 5 + d - 1;
  return era * 146097 + (int)(yy * 365 + yy / 4 - yy / 100 + doy) - 715815;
}
static void days_civil(int days, int *y, unsigned *m, unsigned *d) {
  int z = days + 715815, era = (z >= 0 ? z : z - 146096) / 146097;
  unsigned doe = (unsigned)(z - era * 146097),
           yy = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
  *y = (int)yy + era * 400;
  unsigned doy = doe - (365 * yy + yy / 4 - yy / 100), mp = (5 * doy + 2) / 153;
  *d = doy - (153 * mp + 2) / 5 + 1;
  *m = mp < 10 ? mp + 3 : mp - 9;
  *y += *m <= 2;
}
static int scan_dbf(const char *path) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return 601;
  unsigned char h[32];
  int rc = 610;
  unsigned char *fields = NULL, *row = NULL;
  if (fread(h, 1, 32, fp) != 32 ||
      (h[0] != 3 && h[0] != 4 && h[0] != 0x83 && h[0] != 0x8b))
    goto done;
  unsigned header = dbf_u16(h + 8), record = dbf_u16(h + 10);
  uint32_t n = dbf_u32(h + 4);
  if (header < 33 || record < 1 || (header - 33) % 32 || n > INT_MAX)
    goto done;
  const char *datekeys[] = {"version", "year", "month", "day"};
  for (int j = 0; j < 4; j++) {
    char key[80];
    snprintf(key, sizeof(key), "__cio_dbf_%s", datekeys[j]);
    if ((rc = save_number(key, h[j])))
      goto done;
  }
  rc = 610;
  cache.nvar = (int)((header - 33) / 32);
  if (cache.nvar > 255)
    goto done;
  cache.cols = calloc((size_t)cache.nvar, sizeof(cio_column));
  fields = malloc((size_t)cache.nvar * 32);
  row = malloc(record);
  if ((!cache.cols && cache.nvar) || !fields || !row) {
    rc = 909;
    goto done;
  }
  if (fread(fields, 32, (size_t)cache.nvar, fp) != (size_t)cache.nvar ||
      fgetc(fp) != 13)
    goto done;
  size_t offset = 1;
  for (int j = 0; j < cache.nvar; j++) {
    unsigned char *f = fields + 32 * j;
    cio_column *c = &cache.cols[j];
    memcpy(c->name, f, 11);
    c->name[11] = 0;
    copy_text(c->label, sizeof(c->label), c->name);
    if (!strchr("CNFLD", f[11]) || !f[16] || offset + f[16] > record)
      goto done;
    c->is_string = f[11] == 'C';
    c->width = 1;
    if (c->is_string)
      copy_text(c->format, sizeof(c->format), "%9s");
    else if (f[11] == 'D')
      copy_text(c->format, sizeof(c->format), "%td");
    else
      snprintf(c->format, sizeof(c->format), "%%%u.%uf", f[16], f[17]);
    offset += f[16];
  }
  for (uint32_t i = 0; i < n; i++) {
    if ((i & 4095) == 0 && SF_poll()) {
      rc = 1;
      goto done;
    }
    if (fread(row, 1, record, fp) != record)
      goto done;
    if (row[0] == '*')
      continue;
    if (row[0] != ' ')
      goto done;
    rc = reserve_rows(&cache, cache.nobs + 1);
    if (rc)
      goto done;
    rc = 610;
    size_t obs = cache.nobs++;
    offset = 1;
    for (int j = 0; j < cache.nvar; j++) {
      unsigned char *f = fields + 32 * j;
      cio_column *c = &cache.cols[j];
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
      if (c->is_string) {
        c->strings[obs] = strdup(s);
        if (!c->strings[obs]) {
          rc = 909;
          goto done;
        }
        if ((int)len > c->width)
          c->width = (int)len;
      } else {
        double x = SV_missval;
        if (f[11] == 'L') {
          if (*s == 'T' || *s == 't' || *s == 'Y' || *s == 'y')
            x = 1;
          else if (*s == 'F' || *s == 'f' || *s == 'N' || *s == 'n')
            x = 0;
        } else if (f[11] == 'D' && len == 8) {
          int y, m, d;
          if (sscanf(s, "%4d%2d%2d", &y, &m, &d) == 3 && m >= 1 && m <= 12 &&
              d >= 1 && d <= 31)
            x = civil_days(y, (unsigned)m, (unsigned)d);
        } else if (len) {
          if (s[0] == '.' && s[1] >= 'a' && s[1] <= 'z' && !s[2])
            x = tagged_missing(s[1]);
          else {
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
  for (int j = 0; j < cache.nvar; j++)
    if (cache.cols[j].is_string)
      snprintf(cache.cols[j].format, sizeof(cache.cols[j].format), "%%%ds",
               cache.cols[j].width > 9 ? cache.cols[j].width : 9);
  rc = 0;
done:
  free(fields);
  free(row);
  fclose(fp);
  return rc;
}
static int read_xpf(const char *path) {
  cio_dataset labels;
  memset(&labels, 0, sizeof(labels));
  copy_text(labels.kind, sizeof(labels.kind), "labels");
  readstat_parser_t *p = readstat_parser_init();
  if (!p)
    return 909;
  readstat_set_metadata_handler(p, metadata_cb);
  readstat_set_variable_handler(p, variable_cb);
  readstat_set_value_handler(p, value_cb);
  readstat_error_t err = readstat_parse_xport(p, path, &labels);
  readstat_parser_free(p);
  int rc = err == READSTAT_OK ? 0 : 610, fmt = -1, start = -1, text = -1;
  if (!rc) {
    for (int j = 0; j < labels.nvar; j++) {
      char name[301];
      copy_text(name, sizeof(name), labels.cols[j].name);
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
  for (size_t i = 0; !rc && i < labels.nobs; i++) {
    char *set = labels.cols[fmt].strings[i],
         *label = labels.cols[text].strings[i];
    double val;
    if (labels.cols[start].is_string) {
      const char *str = labels.cols[start].strings[i];
      if (!str)
        continue;
      if (str[0] == '.')
        val = tagged_missing(str[1]);
      else {
        char *end;
        val = strtod(str, &end);
        if (end == str || *end)
          continue;
      }
    } else
      val = labels.cols[start].numbers[i];
    cio_label *all =
        cio_realloc(cache.labels, cache.nlabels + 1, sizeof(cio_label));
    if (!all) {
      rc = 909;
      break;
    }
    cache.labels = all;
    cio_label *l = &all[cache.nlabels++];
    memset(l, 0, sizeof(*l));
    l->set = strdup(set ? set : "");
    l->text = strdup(label ? label : "");
    l->value = val;
    if (!l->set || !l->text)
      rc = 909;
  }
  dataset_free(&labels);
  return rc;
}
/* Implementations share static helpers: workbook inspection follows readers. */
#include "cio_excel_format.h"
/* Stata's built-in and default Excel date text and fixed field widths. Native
   text is the Stata value that the numeric import stores (milliseconds
   days*86400*1000, or that over 86400000 for daily formats) shown through the
   column format, which truncates: 0.999999999999993 of a day is 23:59:59 and
   0.1 (0.09999999999999854 as computed) is 02:23:59. */
static void excel_date_text(double days, uint8_t code, char *out, size_t cap) {
  if (!isfinite(days) || days < (double)INT_MIN || days > (double)INT_MAX) {
    copy_text(out, cap, ".");
    return;
  }
  double ms = days * 86400 * 1000;
  long long whole = 0, seconds = 0;
  if (cio_excel_date_kind(code) == 2) {
    long long t = (long long)floor(ms);
    long long sec = t / 1000 - (t % 1000 < 0);
    whole = sec / 86400 - (sec % 86400 < 0);
    seconds = sec - whole * 86400;
  } else
    whole = (long long)floor(ms / 86400000);
  int year;
  unsigned month, day;
  days_civil((int)whole, &year, &month, &day);
  static const char *months[] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
  static const char *small[] = {"jan", "feb", "mar", "apr", "may", "jun",
                                "jul", "aug", "sep", "oct", "nov", "dec"};
  int total = (int)seconds;
  int h = total / 3600, m = total / 60 % 60, s = total % 60;
  char text[100];
  switch (code) {
  case 14:
    snprintf(text, sizeof(text), "%u/%u/%04d", month, day, year);
    snprintf(out, cap, "%10s", text);
    break;
  case 15:
    snprintf(out, cap, "%2u-%s-%02d", day, months[month - 1],
             (year % 100 + 100) % 100);
    break;
  case 16:
    snprintf(out, cap, "%2u-%s", day, months[month - 1]);
    break;
  case 17:
    snprintf(out, cap, "%s-%02d", months[month - 1], (year % 100 + 100) % 100);
    break;
  case 18:
    snprintf(out, cap, "%2d:%02d %s", h % 12 ? h % 12 : 12, m,
             h < 12 ? "AM" : "PM");
    break;
  case 19:
    snprintf(out, cap, "%2d:%02d:%02d %s", h % 12 ? h % 12 : 12, m, s,
             h < 12 ? "AM" : "PM");
    break;
  case 20:
    snprintf(out, cap, "%2d:%02d", h, m);
    break;
  case 21:
  case 46:
    snprintf(out, cap, "%02d:%02d:%02d", h, m, s);
    break;
  case 22:
    snprintf(text, sizeof(text), "%u/%u/%04d %d:%02d", month, day, year,
             h % 12 ? h % 12 : 12, m);
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
    snprintf(out, cap, "%02u%s%04d %02d:%02d:%02d", day, small[month - 1], year,
             h, m, s);
    break;
  }
}
// clang-format off
#include "cio_shp.inc"
#include "cio_xls.inc"
#include "cio_xlsx.inc"
#include "cio_excel.inc"
#include "cio_xport.inc"
// clang-format on
static int scan(const char *kind) {
  char path[4096], encoding[128], catalog[4096], xpf[4096], raw_storage[16];
  local_text("__cio_filename", path, sizeof(path));
  local_text("__cio_encoding", encoding, sizeof(encoding));
  local_text("__cio_bcat", catalog, sizeof(catalog));
  local_text("__cio_xpf", xpf, sizeof(xpf));
  local_text("__cio_rawstorage", raw_storage, sizeof(raw_storage));
  cio_cleanup();
  copy_text(cache.kind, sizeof(cache.kind), kind);
  cache.raw_storage = !strcmp(kind, "sasxport8") && atoi(raw_storage) != 0;
  if (!strcmp(kind, "dbase") || !strcmp(kind, "xls") || !strcmp(kind, "xlsx") ||
      !strcmp(kind, "shp")) {
    int rc = !strcmp(kind, "xlsx")  ? scan_xlsx(path)
             : !strcmp(kind, "xls") ? scan_xls(path)
             : !strcmp(kind, "shp") ? scan_shp(path)
                                    : scan_dbf(path);
    if (rc) {
      cio_cleanup();
      return rc;
    }
    if ((rc = save_number("__cio_nobs", cache.nobs)))
      return rc;
    if ((rc = save_number("__cio_nvar", (size_t)cache.nvar)))
      return rc;
    if ((rc = save_number("__cio_nlabels", 0)))
      return rc;
    return cio_save("__cio_datalabel", "");
  }
  readstat_parser_t *p = readstat_parser_init();
  if (!p)
    return 909;
  readstat_set_metadata_handler(p, metadata_cb);
  readstat_set_variable_handler(p, variable_cb);
  readstat_set_value_handler(p, value_cb);
  readstat_set_value_label_handler(p, label_cb);
  readstat_set_error_handler(p, error_cb);
  readstat_set_progress_handler(p, progress_cb);
  if (*encoding)
    readstat_set_file_character_encoding(p, encoding);
  readstat_set_handler_character_encoding(p, "UTF-8");
  readstat_error_t err;
  cio_xport_view view;
  memset(&view, 0, sizeof(view));
  if (!strcmp(kind, "sasxport5")) {
    char member[301];
    local_text("__cio_member", member, sizeof(member));
    int rc = xport_select(path, member, p, &view);
    if (rc) {
      readstat_parser_free(p);
      cio_cleanup();
      return rc;
    }
  }
  if (!strcmp(kind, "sas"))
    err = readstat_parse_sas7bdat(p, path, &cache);
  else if (!strcmp(kind, "spss"))
    err = readstat_parse_sav(p, path, &cache);
  else
    err = readstat_parse_xport(p, path, &cache);
  if (err == READSTAT_OK && *catalog) {
    readstat_set_metadata_handler(p, NULL);
    readstat_set_variable_handler(p, NULL);
    readstat_set_value_handler(p, NULL);
    err = readstat_parse_sas7bcat(p, catalog, &cache);
  }
  readstat_parser_free(p);
  if (err != READSTAT_OK) {
    int rc = cache.error ? cache.error : !strcmp(kind, "sasxport5") ? 459 : 692;
    if (!cache.error && err == READSTAT_ERROR_OPEN)
      rc = 601;
    if (!cache.error && err == READSTAT_ERROR_MALLOC)
      rc = 909;
    ctools_error("cimport", "%s", readstat_error_message(err));
    cio_cleanup();
    return rc;
  }
  if (!strcmp(kind, "sasxport5") && *xpf) {
    int rc = read_xpf(xpf);
    if (rc) {
      cio_cleanup();
      return rc;
    }
  }
  for (int j = 0; j < cache.nvar; j++)
    if (cache.cols[j].is_string)
      import_format(&cache.cols[j], kind);
  int rc = save_number("__cio_nobs", cache.nobs);
  if (!rc)
    rc = save_number("__cio_nvar", (size_t)cache.nvar);
  if (!rc)
    rc = save_number("__cio_nlabels", cache.nlabels);
  if (!rc)
    rc = cio_save("__cio_datalabel", cache.datalabel);
  if (!rc && !strcmp(kind, "sasxport5")) {
    size_t width = 4;
    for (int j = 0; j < cache.nvar; j++)
      width += cache.cols[j].is_string ? (size_t)cache.cols[j].width : 8;
    if (cache.nobs && width > SIZE_MAX / cache.nobs)
      rc = 920;
    else
      rc = save_number("__cio_size", width * cache.nobs);
  }
  return rc;
}
static const char *storage(cio_column *c, char *out, size_t len) {
  if (c->binary || (c->is_string && c->width > 2045))
    return "strL";
  if (c->is_string) {
    snprintf(out, len, "str%d", c->width);
    return out;
  }
  if (c->force_double)
    return "double";
  if (!strcmp(cache.kind, "sasxport5"))
    return "double";
  int integral = 1, exactfloat = 1;
  double low = 0, high = 0;
  for (size_t i = 0; i < cache.nobs; i++) {
    double x = c->numbers[i];
    if (SF_is_missing(x))
      continue;
    if (floor(x) != x)
      integral = 0;
    if ((double)(float)x != x)
      exactfloat = 0;
    if (x < low)
      low = x;
    if (x > high)
      high = x;
  }
  if (integral && low >= -127 && high <= 100)
    return "byte";
  if (integral && low >= -32767 && high <= 32740)
    return "int";
  if (integral && low >= -2147483647.0 && high <= 2147483620.0)
    return "long";
  return !strcmp(cache.kind, "dbase") && exactfloat ? "float" : "double";
}
static int column_info(int index) {
  if (index < 1 || index > cache.nvar)
    return 198;
  cio_column *c = &cache.cols[index - 1];
  char s[32];
  int rc;
  if ((rc = save_number("__cio_numerictext", c->format_numbers)))
    return rc;
  if ((rc = save_number("__cio_external",
                        c->binary || (c->is_string && c->width > 2045))))
    return rc;
  if ((rc = cio_save("__cio_excel_column", c->excel_column)))
    return rc;
  if ((rc = cio_save("__cio_name", c->name)))
    return rc;
  if ((rc = cio_save("__cio_type", (char *)storage(c, s, sizeof(s)))))
    return rc;
  if ((rc = cio_save("__cio_label", c->label)))
    return rc;
  if ((rc = cio_save("__cio_format", c->format)))
    return rc;
  return cio_save("__cio_labelset", c->labelset);
}
static int label_info(int index) {
  if (index < 1 || (size_t)index > cache.nlabels)
    return 198;
  cio_label *l = &cache.labels[index - 1];
  char s[40];
  int rc;
  if (SF_is_missing(l->value)) {
    uint64_t bits, base;
    double v = SV_missval;
    memcpy(&bits, &l->value, 8);
    memcpy(&base, &v, 8);
    int n = (int)((bits - base) >> 40);
    snprintf(s, sizeof(s), n ? ".%c" : ".", n ? 'a' + n - 1 : 0);
  } else
    snprintf(s, sizeof(s), "%.17g", l->value);
  if ((rc = cio_save("__cio_labelvalue", s)))
    return rc;
  if ((rc = cio_save("__cio_labeltext", l->text)))
    return rc;
  return cio_save("__cio_labelset", l->set);
}
static int load(void) {
  if (SF_nvars() != cache.nvar || SF_nobs() < (int)cache.nobs)
    return 198;
  for (int j = 0; j < cache.nvar; j++) {
    cio_column *c = &cache.cols[j];
    if (c->binary || (c->is_string && c->width > 2045))
      continue;
    if (!!SF_var_is_string(j + 1) != c->is_string || SF_var_is_strl(j + 1))
      return 109;
    for (size_t i = 0; i < cache.nobs; i++) {
      if ((i & 4095) == 0 && SF_poll())
        return 1;
      int rc =
          c->is_string
              ? SF_sstore(j + 1, (int)i + 1, c->strings[i] ? c->strings[i] : "")
              : (_stata_)->safestore(j + 1, (int)i + 1, c->numbers[i]);
      if (rc)
        return rc;
    }
  }
  return 0;
}

static ssize_t write_cb(const void *data, size_t len, void *ctx) {
  FILE *f = ctx;
  size_t n = fwrite(data, 1, len, f);
  return n == len ? (ssize_t)n : -1;
}
static int export_meta(int index, const char *field, char *out, int cap) {
  char key[100];
  snprintf(key, sizeof(key), "__cio_%s_%d", field, index);
  return local_text(key, out, cap);
}
/* XPORT5 names are uppercase, eight bytes, and unique within each namespace.
   Repeated references to one value-label definition keep the same name. */
static int xport_names(int nvar, char (*names)[301], char (*sets)[301],
                       int rename) {
  char (*original_sets)[301] = calloc((size_t)nvar, 301);
  if (!original_sets)
    return 909;
  int rc = 0;
  for (int field = 0; field < 2 && !rc; field++)
    for (int j = 0; j < nvar; j++) {
      char original[301], base[301];
      export_meta(j + 1, field ? "labelset" : "name", original,
                  sizeof(original));
      char *out = field ? sets[j] : names[j];
      if (field) {
        copy_text(original_sets[j], 301, original);
        int previous = -1;
        for (int k = 0; k < j; k++)
          if (!strcmp(original, original_sets[k])) {
            previous = k;
            break;
          }
        if (previous >= 0) {
          copy_text(out, 301, sets[previous]);
          continue;
        }
      }
      if (!*original) {
        out[0] = 0;
        continue;
      }
      copy_text(base, sizeof(base), original);
      for (char *p = base; *p; p++)
        *p = (char)toupper((unsigned char)*p);
      if (strlen(base) > 8 && !rename) {
        rc = 198;
        break;
      }
      base[8] = 0;
      copy_text(out, 301, base);
      for (int suffix = 2;; suffix++) {
        int collision = 0;
        for (int k = 0; k < j; k++)
          if (!strcmp(out, field ? sets[k] : names[k])) {
            collision = 1;
            break;
          }
        if (!collision)
          break;
        if (!rename) {
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
/* Native export dbase writes numbers with Stata's %20.0g: at most width-2
   significant digits and decimals, trailing zeros and a leading zero dropped,
   and exponent form when the integer part needs more digits or fixed notation
   shows fewer significant digits than the exponent form would. */
static void dbf_number(double x, int width, char *out, size_t cap) {
  if (x == 0) {
    copy_text(out, cap, "0");
    return;
  }
  int digits = width - 2;
  char text[64];
  snprintf(text, sizeof(text), "%.*e", digits - 1, x);
  int exponent = atoi(strchr(text, 'e') + 1),
      mantissa = width - (exponent <= -100 || exponent >= 100 ? 8 : 7);
  if (exponent >= digits || (exponent < 0 && digits + exponent < mantissa)) {
    snprintf(out, cap, "%.*e", mantissa, x);
    return;
  }
  snprintf(out, cap, "%.*f", exponent >= 0 ? digits - 1 - exponent : digits, x);
  size_t n = strlen(out);
  if (strchr(out, '.')) {
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
static int export_dbf(void) {
  char path[4096], replace[8], version[12], datafmt[8];
  local_text("__cio_filename", path, sizeof(path));
  local_text("__cio_replace", replace, sizeof(replace));
  local_text("__cio_version", version, sizeof(version));
  local_text("__cio_datafmt", datafmt, sizeof(datafmt));
  char count[40];
  local_text("__cio_nvar", count, sizeof(count));
  int k = atoi(count), rc = 0;
  uint32_t n = 0;
  size_t record = 1;
  if (k < 1)
    return 102;
  if (k > 255)
    return 103;
  unsigned char *fields = calloc((size_t)k, 32), *row = NULL;
  /* Field kinds: 0 general number, 1 date, 2 int/long, 3 byte. */
  int *kinds = calloc((size_t)k, sizeof(int));
  cexport_output output;
  memset(&output, 0, sizeof(output));
  FILE *fp = NULL;
  if (!fields || !kinds) {
    rc = 909;
    goto done;
  }
  for (int j = 0; j < k; j++) {
    char name[301], type[32], format[257];
    unsigned char *f = fields + 32 * j;
    export_meta(j + 1, "name", name, sizeof(name));
    export_meta(j + 1, "type", type, sizeof(type));
    export_meta(j + 1, "format", format, sizeof(format));
    if (strlen(name) > 10) {
      rc = 693;
      goto done;
    }
    memcpy(f, name, strlen(name));
    int width, decimal = 0;
    if (SF_var_is_string(j + 1)) {
      width = atoi(type + 3);
      f[11] = 'C';
      if (SF_var_is_strl(j + 1)) {
        rc = 693;
        goto done;
      }
      if (width > 254) {
        rc = 108;
        goto done;
      }
    } else if (!strncmp(format, "%td", 3) || !strcmp(format, "%d")) {
      width = 8;
      f[11] = 'D';
      kinds[j] = 1;
    } else if (!strncmp(format, "%t", 2)) {
      rc = 120;
      goto done;
    } else {
      f[11] = 'N';
      kinds[j] = !strcmp(type, "byte")                           ? 3
                 : !strcmp(type, "int") || !strcmp(type, "long") ? 2
                                                                 : 0;
      /* Native widths: byte 3, int 5, long 10, float/double 20 columns */
      width = kinds[j] == 3           ? 3
              : !strcmp(type, "int")  ? 5
              : !strcmp(type, "long") ? 10
                                      : 20;
      if (!strcmp(datafmt, "1")) {
        kinds[j] = 0;
        sscanf(format, "%%%d.%d", &width, &decimal);
        if (width < 1 || width > 254 || decimal > width) {
          rc = 120;
          goto done;
        }
      }
    }
    f[16] = (unsigned char)width;
    f[17] = (unsigned char)decimal;
    record += (size_t)width;
  }
  for (int i = SF_in1(); i <= SF_in2(); i++)
    if (SF_ifobs(i)) {
      n++;
      /* Native's integer fields truncate the most negative values: byte
         -100..-127 to "-10", int -10000..-32767 to "-1000", long up from
         -1000000000 to "-100000000". Widen a field only when such values are
         exported (by one column: byte 4, int 6, long 11). */
      for (int j = 0; j < k; j++)
        if (kinds[j] >= 2) {
          unsigned width = fields[32 * j + 16];
          double bound = 1;
          for (unsigned d = 1; d < width; d++)
            bound *= 10;
          double x;
          if ((rc = (_stata_)->safevdata(j + 1, i, &x)))
            goto done;
          if (!SF_is_missing(x) && x <= -bound) {
            char text[40];
            int need = snprintf(text, sizeof(text), "%.0f", x);
            if (need > (int)width && need < 255) {
              fields[32 * j + 16] = (unsigned char)need;
              record += (size_t)need - width;
            }
          }
        }
    }
  if (record > 4000) {
    rc = 902;
    goto done;
  }
  row = malloc(record);
  if (!row) {
    rc = 909;
    goto done;
  }
  rc = cexport_output_prepare(&output, path, !strcmp(replace, "1"));
  if (rc)
    goto done;
  fp = fopen(output.temporary, "wb");
  if (!fp) {
    rc = 603;
    goto done;
  }
  unsigned char header[32] = {0};
  header[0] = !strcmp(version, "IV") ? 4 : 3;
  time_t now = time(NULL);
  struct tm *date = localtime(&now);
  if (date) {
    header[1] = (unsigned char)date->tm_year;
    header[2] = (unsigned char)(date->tm_mon + 1);
    header[3] = (unsigned char)date->tm_mday;
  }
  char orig[8];
  local_text("__cio_origdbfdate", orig, sizeof(orig));
  if (!strcmp(orig, "1")) {
    const char *keys[] = {"year", "month", "day"};
    for (int j = 0; j < 3; j++) {
      char key[80], value[40];
      snprintf(key, sizeof(key), "__cio_dbf_%s", keys[j]);
      local_text(key, value, sizeof(value));
      if (*value)
        header[j + 1] = (unsigned char)atoi(value);
    }
  }
  dbf_p32(header + 4, n);
  dbf_p16(header + 8, 33 + (unsigned)k * 32);
  dbf_p16(header + 10, (unsigned)record);
  if (fwrite(header, 1, 32, fp) != 32 ||
      fwrite(fields, 32, (size_t)k, fp) != (size_t)k || fputc(13, fp) == EOF) {
    rc = 603;
    goto done;
  }
  for (int i = SF_in1(); i <= SF_in2(); i++) {
    if (!SF_ifobs(i))
      continue;
    if ((i & 4095) == 0 && SF_poll()) {
      rc = 1;
      goto done;
    }
    memset(row, ' ', record);
    size_t offset = 1;
    for (int j = 0; j < k; j++) {
      char value[2046] = {0};
      unsigned char *f = fields + 32 * j;
      if (SF_var_is_string(j + 1)) {
        if ((rc = SF_sdata(j + 1, i, value)))
          goto done;
      } else {
        double x;
        if ((rc = (_stata_)->safevdata(j + 1, i, &x)))
          goto done;
        if (SF_is_missing(x)) {
          if (kinds[j] != 1) {
            uint64_t bits, base;
            double mv = SV_missval;
            memcpy(&bits, &x, 8);
            memcpy(&base, &mv, 8);
            int code = (int)((bits - base) >> 40);
            snprintf(value, sizeof(value), code ? ".%c" : ".",
                     code ? 'a' + code - 1 : 0);
          }
        } else if (kinds[j] == 1) {
          int year;
          unsigned month, day;
          days_civil((int)x, &year, &month, &day);
          snprintf(value, sizeof(value), "%04d%02u%02u", year, month, day);
        } else if (kinds[j] >= 2) {
          snprintf(value, sizeof(value), "%.0f", x);
        } else if (!strcmp(datafmt, "1")) {
          char index[40];
          export_meta(j + 1, "fmtvar", index, sizeof(index));
          int var = atoi(index);
          if (var < k + 1 || var > SF_nvars()) {
            rc = 198;
            goto done;
          }
          if ((rc = SF_sdata(var, i, value)))
            goto done;
        } else
          dbf_number(x, f[16], value, sizeof(value));
      }
      size_t len = strlen(value);
      if (len > f[16]) {
        rc = 108;
        goto done;
      }
      memcpy(row + offset, value, len);
      offset += f[16];
    }
    if (fwrite(row, 1, record, fp) != record) {
      rc = 603;
      goto done;
    }
  }
  if (fputc(26, fp) == EOF) {
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
  free(fields);
  free(row);
  free(kinds);
  return rc;
}
#include "cio_xpf.inc"
static int do_export(const char *kind) {
  if (!strcmp(kind, "shp"))
    return shp_export();
  if (!strcmp(kind, "dbase"))
    return export_dbf();
  char path[4096], replace[10], label[1025];
  cexport_output output;
  memset(&output, 0, sizeof(output));
  local_text("__cio_filename", path, sizeof(path));
  local_text("__cio_replace", replace, sizeof(replace));
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
  local_text("__cio_nmeta", metacount, sizeof(metacount));
  int nmeta = atoi(metacount);
  if (nmeta < nvar)
    nmeta = nvar;
  char (*names)[301] = calloc((size_t)nmeta, 301),
       (*sets)[301] = calloc((size_t)nmeta, 301);
  readstat_variable_t **vars = calloc((size_t)nvar, sizeof(*vars));
  int *dates = calloc((size_t)nvar, sizeof(int));
  char *s = NULL;
  size_t stringcap = 0;
  FILE *fp = NULL;
  if (!vars || !dates || !names || !sets) {
    rc = 909;
    goto done;
  }
  if (!strcmp(kind, "sasxport5")) {
    char rename[8];
    local_text("__cio_rename", rename, sizeof(rename));
    if ((rc = xport_names(nmeta, names, sets, !strcmp(rename, "1"))))
      goto done;
  }
  if (strcmp(kind, "sasxport5"))
    for (int j = 0; j < nmeta; j++)
      export_meta(j + 1, "labelset", sets[j], 301);
  readstat_writer_set_error_handler(w, error_cb);
  readstat_set_data_writer(w, write_cb);
  local_text("__cio_datalabel", label, sizeof(label));
  readstat_writer_set_file_label(w, label);
  char table[40];
  local_text("__cio_table", table, sizeof(table));
  if (*table) {
    if (!strcmp(kind, "sasxport5"))
      table[8] = 0;
    readstat_writer_set_table_name(w, table);
  }
  if (!strncmp(kind, "sasxport", 8))
    readstat_writer_set_file_format_version(w,
                                            !strcmp(kind, "sasxport5") ? 5 : 8);
  else if (!strcmp(kind, "spss"))
    readstat_writer_set_compression(w, READSTAT_COMPRESS_ROWS);
  for (int j = 0; j < nvar; j++) {
    char name[301], type[40], format[257], set[301], count[40];
    export_meta(j + 1, "name", name, sizeof(name));
    export_meta(j + 1, "type", type, sizeof(type));
    export_meta(j + 1, "format", format, sizeof(format));
    export_meta(j + 1, "label", label, sizeof(label));
    if (!strcmp(kind, "sasxport5"))
      copy_text(name, sizeof(name), names[j]);
    int isstr = SF_var_is_string(j + 1), width = isstr ? atoi(type + 3) : 0;
    if (SF_var_is_strl(j + 1)) {
      width = 1;
      for (int i = SF_in1(); i <= SF_in2(); i++)
        if (SF_ifobs(i)) {
          if (SF_var_is_binary(j + 1, i)) {
            rc = 109;
            goto done;
          }
          int length = SF_sdatalen(j + 1, i);
          if (length < 0) {
            rc = 109;
            goto done;
          }
          size_t cap = (size_t)length + 1;
          if (cap > stringcap) {
            char *tmp = realloc(s, cap);
            if (!tmp) {
              rc = 909;
              goto done;
            }
            s = tmp;
            stringcap = cap;
          }
          if (SF_strldata(j + 1, i, s, (int)stringcap) != length) {
            rc = 609;
            goto done;
          }
          s[length] = 0;
          size_t logical = strlen(s) + (!strcmp(kind, "sasxport8") ? 1 : 0);
          if (logical > INT_MAX) {
            rc = 109;
            goto done;
          }
          if ((int)logical > width)
            width = (int)logical;
        }
    }
    vars[j] = readstat_add_variable(
        w, name, isstr ? READSTAT_TYPE_STRING : READSTAT_TYPE_DOUBLE,
        (size_t)width);
    if (!vars[j]) {
      rc = 909;
      goto done;
    }
    readstat_variable_set_label(vars[j], label);
    if (isstr && !strcmp(kind, "sasxport8")) {
      char stringformat[80];
      int display = 0;
      sscanf(format, "%%%d", &display);
      if (SF_var_is_strl(j + 1))
        display = 8;
      snprintf(stringformat, sizeof(stringformat), "$%d",
               display > 0 ? display : 8);
      readstat_variable_set_format(vars[j], stringformat);
    }
    if (!isstr) {
      char fmt[80];
      if (!strncmp(format, "%td", 3))
        dates[j] = 1;
      else if (!strncmp(format, "%tc", 3) || !strncmp(format, "%tC", 3))
        dates[j] = 2;
      if (!strcmp(kind, "spss")) {
        if (dates[j] == 1)
          copy_text(fmt, sizeof(fmt), "DATE11");
        else if (dates[j] == 2)
          copy_text(fmt, sizeof(fmt), "DATETIME23.2");
        else {
          int fw = 0, dec = 0;
          sscanf(format, "%%%d.%d", &fw, &dec);
          if (fw < 8)
            fw = 8;
          snprintf(fmt, sizeof(fmt), "F%d.%d", fw + 3, dec);
        }
      } else {
        if (dates[j] == 1)
          copy_text(fmt, sizeof(fmt), "DATE");
        else if (dates[j] == 2)
          copy_text(fmt, sizeof(fmt),
                    !strcmp(kind, "sasxport5") ? "" : "DATETIME");
        else
          copy_text(fmt, sizeof(fmt), "");
      }
      readstat_variable_set_format(vars[j], fmt);
    }
    export_meta(j + 1, "labelset", set, sizeof(set));
    export_meta(j + 1, "lcount", count, sizeof(count));
    if (!strcmp(kind, "sasxport5"))
      copy_text(set, sizeof(set), sets[j]);
    if (*set && !strcmp(kind, "sasxport5")) {
      for (char *p = set; *p; p++)
        *p = (char)toupper((unsigned char)*p);
      readstat_variable_set_format(vars[j], set);
    }
    if (*set && !strcmp(kind, "spss")) {
      readstat_label_set_t *ls = NULL;
      for (long k = 0; k < w->label_sets_count; k++)
        if (!strcmp(w->label_sets[k]->name, set)) {
          ls = w->label_sets[k];
          break;
        }
      int new_set = !ls;
      if (!ls)
        ls = readstat_add_label_set(w, READSTAT_TYPE_DOUBLE, set);
      if (!ls) {
        rc = 909;
        goto done;
      }
      readstat_variable_set_label_set(vars[j], ls);
      int n = atoi(count);
      for (int k = 1; new_set && k <= n; k++) {
        char key[100], val[64], text[32768];
        snprintf(key, sizeof(key), "__cio_lv_%d_%d", j + 1, k);
        local_text(key, val, sizeof(val));
        snprintf(key, sizeof(key), "__cio_lt_%d_%d", j + 1, k);
        local_text(key, text, sizeof(text));
        if (*val != '.')
          readstat_label_double_value(ls, strtod(val, NULL), text);
      }
    }
    if (isstr && (size_t)width + 1 > stringcap)
      stringcap = (size_t)width + 1;
  }
  if (stringcap) {
    char *tmp = realloc(s, stringcap);
    if (!tmp) {
      rc = 909;
      goto done;
    }
    s = tmp;
  }
  /* Use a sibling temporary and atomic publish; failed writes never damage
     an existing destination. Never unlink the destination before writing. */
  rc = cexport_output_prepare(&output, path, !strcmp(replace, "1"));
  if (rc)
    goto done;
  fp = fopen(output.temporary, "wb");
  if (!fp) {
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
  for (int i = SF_in1(); i <= SF_in2(); i++) {
    if (!SF_ifobs(i))
      continue;
    if ((i & 4095) == 0 && SF_poll()) {
      rc = 1;
      goto done;
    }
    if ((err = readstat_begin_row(w)) != READSTAT_OK)
      goto codec_error;
    for (int j = 0; j < nvar; j++) {
      if (SF_var_is_string(j + 1)) {
        if (SF_var_is_strl(j + 1)) {
          int length = SF_sdatalen(j + 1, i);
          if (length < 0 || (size_t)length >= stringcap) {
            rc = 109;
            goto done;
          }
          if (SF_strldata(j + 1, i, s, (int)stringcap) != length) {
            rc = 609;
            goto done;
          }
          s[length] = 0;
        } else if ((rc = SF_sdata(j + 1, i, s)))
          goto done;
        err = readstat_insert_string_value(w, vars[j], s);
      } else {
        double x;
        if ((rc = (_stata_)->safevdata(j + 1, i, &x)))
          goto done;
        if (SF_is_missing(x)) {
          uint64_t bits, base;
          double mv = SV_missval;
          memcpy(&bits, &x, 8);
          memcpy(&base, &mv, 8);
          int n = (int)((bits - base) >> 40);
          if (n && strcmp(kind, "spss"))
            err = readstat_insert_tagged_missing_value(w, vars[j],
                                                       (char)('A' + n - 1));
          else
            err = readstat_insert_missing_value(w, vars[j]);
        } else {
          if (!strcmp(kind, "spss") && dates[j])
            x = ((dates[j] == 1 ? x * 86400000 : x) + 11903760000000.0) / 1000;
          else if (!strcmp(kind, "sasxport8") && dates[j] == 1)
            x *= 86400000;
          else if (!strncmp(kind, "sasxport",
                            8)) { /* Transport preserves Stata numeric units. */
          } else if (dates[j] >= 2)
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
  if (fclose(fp)) {
    fp = NULL;
    rc = 603;
    goto done;
  }
  fp = NULL;
  if (!strcmp(kind, "sasxport5")) {
    char xpf[4096];
    local_text("__cio_xpfoutput", xpf, sizeof(xpf));
    if (*xpf && (rc = export_xpf(nmeta, xpf, !strcmp(replace, "1"), sets)))
      goto done;
  }
  char sascode[4096];
  local_text("__cio_sascode", sascode, sizeof(sascode));
  if (*sascode &&
      (rc = export_sascode(nmeta, sets, sascode, kind, !strcmp(replace, "1"))))
    goto done;
  rc = cexport_output_commit(&output);
  goto done;
codec_error:
  ctools_error("cexport", "%s", readstat_error_message(err));
  rc = 610;
done:
  if (fp)
    fclose(fp);
  cexport_output_cleanup(&output);
  free(s);
  free(vars);
  free(dates);
  free(names);
  free(sets);
  readstat_writer_free(w);
  return rc;
}
ST_retcode cio_main(const char *args) {
  char op[24], kind[24];
  int index = 0;
  if (sscanf(args, "%23s %23s %d", op, kind, &index) < 1)
    return 198;
  if (!strcmp(op, "workbook"))
    return excel_workbook(kind);
  if (!strcmp(op, "blob")) {
    int col, row;
    size_t offset;
    if (sscanf(args, "%*s %d %d %zu", &col, &row, &offset) != 3)
      return 198;
    return blob_info(col, row, offset);
  }
  if (!strcmp(op, "clear")) {
    cio_cleanup();
    return 0;
  }
  if (!strcmp(op, "scan"))
    return scan(kind);
  if (!strcmp(op, "column"))
    return column_info(atoi(kind));
  if (!strcmp(op, "textnumbers")) {
    int j = atoi(kind) - 1;
    if (j < 0 || j >= cache.nvar || !cache.cols[j].format_numbers ||
        SF_nvars() != 1 || SF_var_is_string(1) || SF_nobs() < (int)cache.nobs)
      return 198;
    for (size_t i = 0; i < cache.nobs; i++) {
      int rc = (_stata_)->safestore(1, (int)i + 1, cache.cols[j].numbers[i]);
      if (rc)
        return rc;
    }
    return 0;
  }
  if (!strcmp(op, "label"))
    return label_info(atoi(kind));
  if (!strcmp(op, "load"))
    return load();
  if (!strcmp(op, "export"))
    return do_export(kind);
  return 198;
}

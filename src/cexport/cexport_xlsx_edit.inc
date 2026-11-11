/* ZIP/XML workbook editing. Untouched parts are copied as compressed streams;
   worksheet edits retain raw XML for every cell outside the exported range. */
typedef struct {
  const char *start, *open_end, *close, *end;
  const char *name;
  size_t name_len;
  int empty;
} xedit_node;

/* Lex tags without interpreting text, formulas, comments, or quoted > signs. */
static const char *xedit_tag(const char *p, const char *limit, xedit_node *n,
                             int *closing) {
  while (p < limit) {
    if (*p++ != '<')
      continue;
    const char *start = p - 1;
    if (p == limit)
      return NULL;
    if (*p == '!' || *p == '?') {
      const char *terminator =
          *p == '?'                                       ? "?>"
          : (limit - p >= 3 && !memcmp(p, "!--", 3))      ? "-->"
          : (limit - p >= 8 && !memcmp(p, "![CDATA[", 8)) ? "]]>"
                                                          : ">";
      size_t width = strlen(terminator);
      while (p < limit &&
             ((size_t)(limit - p) < width || memcmp(p, terminator, width)))
        p++;
      if (p == limit)
        return NULL;
      p += width;
      continue;
    }
    *closing = *p == '/';
    if (*closing)
      p++;
    const char *name = p;
    while (p < limit && !isspace((unsigned char)*p) && *p != '/' && *p != '>')
      p++;
    size_t name_len = (size_t)(p - name);
    char quote = 0;
    while (p < limit) {
      if (quote) {
        if (*p == quote)
          quote = 0;
      } else if (*p == '\'' || *p == '"')
        quote = *p;
      else if (*p == '>')
        break;
      p++;
    }
    if (p == limit || !name_len)
      return NULL;
    *n = (xedit_node){
        start, p + 1, p + 1, p + 1, name, name_len, p > start && p[-1] == '/'};
    return p + 1;
  }
  return NULL;
}
static int xedit_is(const xedit_node *n, const char *name) {
  const char *local = n->name;
  size_t len = n->name_len;
  for (size_t j = 0; j < n->name_len; j++)
    if (n->name[j] == ':') {
      local = n->name + j + 1;
      len = n->name_len - j - 1;
    }
  return strlen(name) == len && !memcmp(local, name, len);
}
static int xedit_find(const char *p, const char *limit, const char *name,
                      xedit_node *out) {
  xedit_node n;
  int closing;
  while ((p = xedit_tag(p, limit, &n, &closing))) {
    if (closing || !xedit_is(&n, name))
      continue;
    *out = n;
    if (n.empty)
      return 1;
    int depth = 1;
    while ((p = xedit_tag(p, limit, &n, &closing))) {
      if (closing)
        depth--;
      else if (!n.empty)
        depth++;
      if (!depth) {
        out->close = n.start;
        out->end = n.end;
        return 1;
      }
    }
    return 0;
  }
  return 0;
}
/* Attribute values are decoded only for lookup; preserved XML stays byte exact.
 */
static int xedit_attr_span(const xedit_node *n, const char *key, const char **a,
                           const char **b, const char **v, size_t *len) {
  const char *p = n->name + n->name_len, *end = n->open_end - 1;
  while (p < end) {
    while (p < end && isspace((unsigned char)*p))
      p++;
    const char *start = p;
    while (p < end && !isspace((unsigned char)*p) && *p != '=' && *p != '/' &&
           *p != '>')
      p++;
    size_t k = (size_t)(p - start);
    while (p < end && isspace((unsigned char)*p))
      p++;
    if (p == end || *p != '=')
      break;
    p++;
    while (p < end && isspace((unsigned char)*p))
      p++;
    if (p == end || (*p != '\'' && *p != '"'))
      break;
    char quote = *p++;
    const char *value = p;
    while (p < end && *p != quote)
      p++;
    if (p == end)
      break;
    if (k == strlen(key) && !memcmp(start, key, k)) {
      *a = start;
      *b = p + 1;
      *v = value;
      *len = (size_t)(p - value);
      return 1;
    }
    p++;
  }
  return 0;
}
static char *xedit_attr(const xedit_node *n, const char *key) {
  const char *a, *b, *v;
  size_t len;
  if (!xedit_attr_span(n, key, &a, &b, &v, &len))
    return NULL;
  char *s = malloc(len + 1);
  if (!s)
    return NULL;
  memcpy(s, v, len);
  s[len] = 0;
  len = xlsx_xml_decode_entities(s, len);
  s[len] = 0;
  return s;
}
static char *xedit_set_attr(const xedit_node *n, const char *key,
                            unsigned value) {
  StringBuffer sb;
  if (!strbuf_init(&sb, (size_t)(n->end - n->start) + 64))
    return NULL;
  const char *a, *b, *v;
  size_t len;
  if (xedit_attr_span(n, key, &a, &b, &v, &len)) {
    strbuf_append_n(&sb, n->start, (size_t)(a - n->start));
    strbuf_appendf(&sb, "%s=\"%u\"", key, value);
    strbuf_append_n(&sb, b, (size_t)(n->end - b));
  } else {
    const char *cut = n->open_end - (n->empty ? 2 : 1);
    strbuf_append_n(&sb, n->start, (size_t)(cut - n->start));
    strbuf_appendf(&sb, " %s=\"%u\"", key, value);
    strbuf_append_n(&sb, cut, (size_t)(n->end - cut));
  }
  return sb.data;
}
static char *xedit_replace(const char *xml, const xedit_node *n,
                           const char *replacement) {
  StringBuffer sb;
  if (!strbuf_init(&sb, strlen(xml) + strlen(replacement) + 64))
    return NULL;
  strbuf_append_n(&sb, xml, (size_t)(n->start - xml));
  strbuf_append(&sb, replacement);
  strbuf_append(&sb, n->end);
  return sb.data;
}
static char *xedit_insert(const char *xml, const char *parent,
                          const char *child) {
  xedit_node n;
  if (!xedit_find(xml, xml + strlen(xml), parent, &n))
    return NULL;
  StringBuffer sb;
  if (!strbuf_init(&sb, strlen(xml) + strlen(child) + 64))
    return NULL;
  const char *cut = n.empty ? n.open_end - 2 : n.close;
  strbuf_append_n(&sb, xml, (size_t)(cut - xml));
  if (n.empty)
    strbuf_append(&sb, ">");
  strbuf_append(&sb, child);
  if (n.empty) {
    strbuf_append(&sb, "</");
    strbuf_append_n(&sb, n.name, n.name_len);
    strbuf_append(&sb, ">");
  }
  strbuf_append(&sb, n.empty ? n.end : n.close);
  return sb.data;
}
static char *xedit_zip_text(mz_zip_archive *zip, const char *name) {
  size_t len;
  char *raw = mz_zip_reader_extract_file_to_heap(zip, name, &len, 0);
  if (!raw)
    return NULL;
  char *s = malloc(len + 1);
  if (s) {
    memcpy(s, raw, len);
    s[len] = 0;
  }
  mz_free(raw);
  return s;
}
static int xedit_cell(const xedit_node *n, int *col, int *row) {
  char *r = xedit_attr(n, "r");
  int ok = r && xlsx_xml_parse_cell_ref(r, col, row);
  free(r);
  return ok;
}
static int xedit_row(const xedit_node *n) {
  char *r = xedit_attr(n, "r");
  int i = r ? atoi(r) : 0;
  free(r);
  return i;
}
static void xedit_write_cell(StringBuffer *sb, const xedit_node *fresh,
                             const xedit_node *old, int keep,
                             unsigned style_base) {
  char *style = xedit_attr(fresh, "s");
  unsigned index = style ? (unsigned)strtoul(style, NULL, 10) : 0;
  free(style);
  /* keepcellfmt keeps the old style of every cell but one with a %t
   * format, which keeps its own, as native export does; fresh numbers then
   * carry no number format. */
  if (keep && old && !CEXPORT_STYLE_IS_TFORMAT(index)) {
    style = xedit_attr(old, "s");
    index = style ? (unsigned)strtoul(style, NULL, 10) : 0;
    free(style);
  } else
    index += style_base;
  char *cell = xedit_set_attr(fresh, "s", index);
  if (!cell) {
    free(sb->data);
    sb->data = NULL;
    return;
  }
  strbuf_append(sb, cell);
  free(cell);
}
static int xedit_next_cell(const char **p, const char *end, xedit_node *n,
                           int *col) {
  int row;
  if (!xedit_find(*p, end, "c", n) || !xedit_cell(n, col, &row))
    return 0;
  *p = n->end;
  return 1;
}
static void xedit_merge_row(StringBuffer *sb, const xedit_node *old,
                            const xedit_node *fresh, int keep, unsigned base) {
  if (!fresh) {
    strbuf_append_n(sb, old->start, (size_t)(old->end - old->start));
    return;
  }
  if (old && !old->empty)
    strbuf_append_n(sb, old->start, (size_t)(old->open_end - old->start));
  else
    strbuf_appendf(
        sb,
        "<row r=\"%d\" "
        "xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">",
        xedit_row(fresh));
  const char *op = old && !old->empty ? old->open_end : NULL,
             *fp = fresh->open_end;
  xedit_node oc, fc;
  int oi = 0, fi = 0;
  int haso = op && xedit_next_cell(&op, old->close, &oc, &oi);
  int hasf = xedit_next_cell(&fp, fresh->close, &fc, &fi);
  while (haso || hasf) {
    if (haso && (!hasf || oi < fi)) {
      strbuf_append_n(sb, oc.start, (size_t)(oc.end - oc.start));
      haso = xedit_next_cell(&op, old->close, &oc, &oi);
    } else {
      xedit_write_cell(sb, &fc, haso && oi == fi ? &oc : NULL, keep, base);
      if (haso && oi == fi)
        haso = xedit_next_cell(&op, old->close, &oc, &oi);
      hasf = xedit_next_cell(&fp, fresh->close, &fc, &fi);
    }
  }
  if (old && !old->empty)
    strbuf_append_n(sb, old->close, (size_t)(old->end - old->close));
  else
    strbuf_append(sb, "</row>");
}
/* Give every old row and cell an explicit r, as readers infer it (a row
   without r follows the previous row, a cell without r the previous cell in
   its row), so the merge, the used range and shared formulas see each cell
   where Excel does. Returns xml itself when every position is explicit. */
static char *xedit_explicit_positions(char *xml) {
  xedit_node sd, n;
  int closing, missing = 0;
  const char *a, *b, *v;
  size_t vl, len = strlen(xml);
  if (!xedit_find(xml, xml + len, "sheetData", &sd) || sd.empty)
    return xml;
  for (const char *p = sd.open_end;
       !missing && (p = xedit_tag(p, sd.close, &n, &closing));) {
    int col, row;
    if (!closing && xedit_is(&n, "row"))
      missing = xedit_row(&n) <= 0;
    else if (!closing && xedit_is(&n, "c"))
      missing = !xedit_cell(&n, &col, &row);
  }
  if (!missing)
    return xml;
  StringBuffer sb;
  if (!strbuf_init(&sb, len + len / 4 + 64))
    return NULL;
  const char *copied = xml;
  int row = 0, col = -1;
  for (const char *p = sd.open_end; (p = xedit_tag(p, sd.close, &n, &closing));) {
    char text[32] = "";
    if (closing)
      continue;
    if (xedit_is(&n, "row")) {
      int i = xedit_row(&n);
      row = i > 0 ? i : row + 1;
      col = -1;
      if (i <= 0)
        snprintf(text, sizeof(text), "%d", row);
    } else if (xedit_is(&n, "c")) {
      int c, r;
      if (xedit_cell(&n, &c, &r))
        col = c;
      else {
        char letters[16];
        xlsx_col_to_letters(++col + 1, letters);
        snprintf(text, sizeof(text), "%s%d", letters, row);
      }
    }
    if (!*text)
      continue;
    /* Replace an unusable r, or add one before the end of the start tag. */
    const char *cut = n.open_end - (n.empty ? 2 : 1), *resume = cut;
    int present = xedit_attr_span(&n, "r", &a, &b, &v, &vl);
    if (present) {
      cut = a;
      resume = b;
    }
    strbuf_append_n(&sb, copied, (size_t)(cut - copied));
    strbuf_appendf(&sb, "%sr=\"%s\"", present ? "" : " ", text);
    copied = resume;
  }
  strbuf_append(&sb, copied);
  return sb.data;
}
#include "cexport_xlsx_formula.inc"
static char *xedit_worksheet(const char *oldxml, const char *newxml, int keep,
                             unsigned base) {
  xedit_node od, nd;
  if (!xedit_find(oldxml, oldxml + strlen(oldxml), "sheetData", &od) ||
      !xedit_find(newxml, newxml + strlen(newxml), "sheetData", &nd))
    return NULL;
  StringBuffer sb;
  if (!strbuf_init(&sb, strlen(oldxml) + strlen(newxml) + 256))
    return NULL;
  strbuf_append_n(&sb, oldxml, (size_t)(od.start - oldxml));
  strbuf_append(
      &sb,
      "<sheetData "
      "xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">");
  const char *op = od.open_end, *np = nd.open_end;
  xedit_node o, n;
  int haso = !od.empty && xedit_find(op, od.close, "row", &o);
  int hasn = !nd.empty && xedit_find(np, nd.close, "row", &n);
  while (haso || hasn) {
    int oi = haso ? xedit_row(&o) : INT_MAX,
        ni = hasn ? xedit_row(&n) : INT_MAX;
    if (oi < ni)
      xedit_merge_row(&sb, &o, NULL, keep, base);
    else
      xedit_merge_row(&sb, haso && oi == ni ? &o : NULL, &n, keep, base);
    if (haso && oi <= ni) {
      op = o.end;
      haso = xedit_find(op, od.close, "row", &o);
    }
    if (hasn && ni <= oi) {
      np = n.end;
      hasn = xedit_find(np, nd.close, "row", &n);
    }
  }
  strbuf_append(&sb, "</sheetData>");
  strbuf_append(&sb, od.end);
  if (!sb.data)
    return NULL;
  /* Recompute the used range, including explicitly blank exported cells. */
  xedit_node cell, dim;
  int minc = INT_MAX, minr = INT_MAX, maxc = -1, maxr = -1, col, row;
  const char *p = sb.data, *end = p + sb.len;
  while (xedit_find(p, end, "c", &cell)) {
    if (xedit_cell(&cell, &col, &row)) {
      if (col < minc)
        minc = col;
      if (col > maxc)
        maxc = col;
      if (row < minr)
        minr = row;
      if (row > maxr)
        maxr = row;
    }
    p = cell.end;
  }
  if (maxc >= 0 && xedit_find(sb.data, end, "dimension", &dim)) {
    char a[16], b[16], text[128];
    xlsx_col_to_letters(minc + 1, a);
    xlsx_col_to_letters(maxc + 1, b);
    snprintf(
        text, sizeof(text),
        "<dimension "
        "xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
        "ref=\"%s%d:%s%d\"/>",
        a, minr, b, maxr);
    char *result = xedit_replace(sb.data, &dim, text);
    free(sb.data);
    return result;
  }
  return sb.data;
}
/* The number format of a cell format that sets nothing else (font, fill,
   border and parent style 0, no alignment or protection), or UINT_MAX. */
static unsigned xedit_plain_xf(const xedit_node *xf) {
  static const char *const keys[] = {"numFmtId", "fontId", "fillId",
                                     "borderId", "xfId"};
  unsigned values[5];
  if (!xf->empty)
    return UINT_MAX;
  for (int k = 0; k < 5; k++) {
    char *v = xedit_attr(xf, keys[k]);
    values[k] = v ? (unsigned)strtoul(v, NULL, 10) : 0;
    free(v);
  }
  return values[1] || values[2] || values[3] || values[4] ? UINT_MAX
                                                          : values[0];
}
/* Excel may escape the format code's literals when it saves the workbook. */
static int xedit_same_format(const char *code, const char *ours) {
  for (; *code; code++) {
    if (*code == '\\')
      continue;
    if (*code != *ours++)
      return 0;
  }
  return !*ours;
}
/* Our cell styles, one per CEXPORT_STYLE_* code in order, at *base. A run of
   them already in the workbook (one that ctools wrote, or added by an
   earlier edit) is reused, as is our date-time format; otherwise they are
   appended, retaining every existing style index. Appending on every edit
   would pass Excel's limits on cell formats and custom number formats. */
static char *xedit_styles(char *xml, unsigned *base) {
  static const char ours[] = "yyyy-mm-dd hh:mm:ss.000";
  xedit_node xfs, xf, nfmt;
  if (!xedit_find(xml, xml + strlen(xml), "cellXfs", &xfs))
    return NULL;
  unsigned fmt = 164, have_fmt = 0, found_fmt = 0;
  const char *p = xml, *xml_end = xml + strlen(xml);
  while (xedit_find(p, xml_end, "numFmt", &nfmt)) {
    char *id = xedit_attr(&nfmt, "numFmtId"), *code = xedit_attr(&nfmt, "formatCode");
    if (id && strtoul(id, NULL, 10) >= fmt)
      fmt = (unsigned)strtoul(id, NULL, 10) + 1;
    if (id && code && !have_fmt && xedit_same_format(code, ours)) {
      found_fmt = (unsigned)strtoul(id, NULL, 10);
      have_fmt = 1;
    }
    free(id);
    free(code);
    p = nfmt.end;
  }
  if (have_fmt)
    fmt = found_fmt;
  const unsigned ids[CEXPORT_STYLE_COUNT] = {0, 14, fmt, 22, 19, 1, 1, 2, 3, 4, 11};
  *base = 0;
  size_t cap = 0;
  unsigned *plain = NULL;
  p = xfs.open_end;
  while (!xfs.empty && xedit_find(p, xfs.close, "xf", &xf)) {
    if (have_fmt && *base == cap) {
      size_t grown_cap = cap ? 2 * cap : 64;
      unsigned *grown = realloc(plain, grown_cap * sizeof(*grown));
      if (!grown) {
        free(plain);
        return NULL;
      }
      plain = grown;
      cap = grown_cap;
    }
    if (have_fmt)
      plain[*base] = xedit_plain_xf(&xf);
    (*base)++;
    p = xf.end;
  }
  for (unsigned k = 0; have_fmt && k + CEXPORT_STYLE_COUNT <= *base; k++) {
    int j = 0;
    while (j < CEXPORT_STYLE_COUNT && plain[k + (unsigned)j] == ids[j])
      j++;
    if (j == CEXPORT_STYLE_COUNT) {
      free(plain);
      *base = k;
      return strdup(xml);
    }
  }
  free(plain);
  char child[512];
  snprintf(
      child, sizeof(child),
      "<numFmt "
      "xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
      "numFmtId=\"%u\" formatCode=\"yyyy-mm-dd hh:mm:ss.000\"/>",
      fmt);
  char *result;
  if (have_fmt)
    result = strdup(xml);
  else if (xedit_find(xml, xml + strlen(xml), "numFmts", &nfmt)) {
    char *count = xedit_attr(&nfmt, "count");
    unsigned total = count ? (unsigned)strtoul(count, NULL, 10) : 0;
    free(count);
    char *node = xedit_set_attr(&nfmt, "count", total + 1);
    if (!node)
      return NULL;
    char *updated = xedit_replace(xml, &nfmt, node);
    free(node);
    if (!updated)
      return NULL;
    result = xedit_insert(updated, "numFmts", child);
    free(updated);
  } else {
    StringBuffer sb;
    if (!strbuf_init(&sb, strlen(xml) + 1024))
      return NULL;
    xedit_node root;
    if (!xedit_find(xml, xml + strlen(xml), "styleSheet", &root)) {
      free(sb.data);
      return NULL;
    }
    strbuf_append_n(&sb, xml, (size_t)(root.open_end - xml));
    strbuf_append(
        &sb,
        "<numFmts "
        "xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
        "count=\"1\">");
    strbuf_append(&sb, child);
    strbuf_append(&sb, "</numFmts>");
    strbuf_append(&sb, root.open_end);
    result = sb.data;
  }
  if (!result)
    return NULL;
  /* xfs still points into xml; never apply it to result. */
  if (!xedit_find(result, result + strlen(result), "cellXfs", &xfs)) {
    free(result);
    return NULL;
  }
  char *node = xedit_set_attr(&xfs, "count", *base + CEXPORT_STYLE_COUNT);
  if (!node) {
    free(result);
    return NULL;
  }
  char *updated = xedit_replace(result, &xfs, node);
  free(node);
  free(result);
  if (!updated)
    return NULL;
  StringBuffer cells;
  if (!strbuf_init(&cells, 1536)) {
    free(updated);
    return NULL;
  }
  for (int i = 0; i < CEXPORT_STYLE_COUNT; i++)
    strbuf_appendf(
        &cells,
        "<xf xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
        "numFmtId=\"%u\" fontId=\"0\" fillId=\"0\" borderId=\"0\" xfId=\"0\"%s/>",
        ids[i], ids[i] ? " applyNumberFormat=\"1\"" : "");
  result = cells.data ? xedit_insert(updated, "cellXfs", cells.data) : NULL;
  free(cells.data);
  free(updated);
  return result;
}

/* Simple Unicode case folding for the scripts sheet names mostly use: ASCII,
   Latin-1, Latin Extended-A, Greek and Cyrillic. */
static unsigned xedit_fold(unsigned c) {
  if ((c >= 'A' && c <= 'Z') || (c >= 0xC0 && c <= 0xDE && c != 0xD7) ||
      (c >= 0x391 && c <= 0x3A9 && c != 0x3A2) || (c >= 0x410 && c <= 0x42F))
    return c + 0x20;
  if (c >= 0x400 && c <= 0x40F)
    return c + 0x50;
  if (c == 0x178)
    return 0xFF;
  if (c == 0x17F)
    return 's';
  if (((c >= 0x100 && c <= 0x12F) || (c >= 0x132 && c <= 0x137) ||
       (c >= 0x14A && c <= 0x177)) && !(c & 1))
    return c + 1;
  if (((c >= 0x139 && c <= 0x148) || (c >= 0x179 && c <= 0x17E)) && (c & 1))
    return c + 1;
  return c;
}
static unsigned xedit_next_codepoint(const unsigned char **s) {
  const unsigned char *p = *s;
  unsigned c = *p++;
  int more = c >= 0xF0 ? 3 : c >= 0xE0 ? 2 : c >= 0xC0 ? 1 : 0;
  if (more)
    c &= 0x3F >> more;
  while (more-- && (*p & 0xC0) == 0x80)
    c = (c << 6) | (*p++ & 0x3F);
  *s = p;
  return c;
}
/* Whether two UTF-8 sheet names are equal ignoring case, as Excel compares
   them: it refuses a sheet whose name differs from another only in case. */
static int xedit_same_name(const char *a, const char *b) {
  const unsigned char *p = (const unsigned char *)a, *q = (const unsigned char *)b;
  while (*p && *q)
    if (xedit_fold(xedit_next_codepoint(&p)) != xedit_fold(xedit_next_codepoint(&q)))
      return 0;
  return !*p && !*q;
}
/* The sheet element's relationship id: r:id, or the id attribute under
   another prefix bound to the relationships namespace (ns1:id), as some
   writers declare it. A sheet's only prefixed id attribute is that one. */
static char *xedit_rel_id(const xedit_node *n) {
  char *r = xedit_attr(n, "r:id");
  if (r)
    return r;
  for (const char *p = n->name + n->name_len; p + 3 <= n->open_end; p++) {
    if (memcmp(p, ":id", 3))
      continue;
    const char *q = p + 3;
    while (q < n->open_end && isspace((unsigned char)*q))
      q++;
    if (q >= n->open_end || *q != '=')
      continue;
    const char *s = p;
    while (s > n->start && !isspace((unsigned char)s[-1]))
      s--;
    if (s == p || p - s > 32)
      continue;
    char key[40];
    snprintf(key, sizeof(key), "%.*s:id", (int)(p - s), s);
    return xedit_attr(n, key);
  }
  return NULL;
}
/* Whether the workbook relationship rid is a worksheet rather than a chart,
   dialog or macro sheet, whose parts hold no sheetData. A relationship that
   rels does not list counts as a worksheet: the target lookup reports it. */
static int xedit_rel_worksheet(const char *rels, const char *rid) {
  const char *p = rels, *end = rels + strlen(rels);
  xedit_node node;
  while (xedit_find(p, end, "Relationship", &node)) {
    char *id = xedit_attr(&node, "Id");
    int match = id && !strcmp(id, rid);
    free(id);
    if (match) {
      char *type = xedit_attr(&node, "Type");
      size_t n = type ? strlen(type) : 0;
      int worksheet = !type || (n >= 10 && !strcmp(type + n - 10, "/worksheet"));
      free(type);
      return worksheet;
    }
    p = node.end;
  }
  return 1;
}
/* Strict Open XML packages name the same elements in other namespaces; the
   parts written here (injected elements declare the transitional ones) then
   take the package's own. OPC namespaces (content types, package
   relationships) are the same in both. */
static char *xedit_strict(char *xml) {
  static const char *const from[] = {
      "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "http://schemas.openxmlformats.org/officeDocument/2006/relationships"};
  static const char *const to[] = {
      "http://purl.oclc.org/ooxml/spreadsheetml/main",
      "http://purl.oclc.org/ooxml/officeDocument/relationships"};
  StringBuffer sb;
  if (!strbuf_init(&sb, strlen(xml) + 64))
    return NULL;
  const char *p = xml;
  while (*p) {
    int k = -1;
    const char *hit = NULL;
    for (int j = 0; j < 2; j++) {
      const char *h = strstr(p, from[j]);
      if (h && (!hit || h < hit)) {
        hit = h;
        k = j;
      }
    }
    if (!hit) {
      strbuf_append(&sb, p);
      break;
    }
    strbuf_append_n(&sb, p, (size_t)(hit - p));
    strbuf_append(&sb, to[k]);
    p = hit + strlen(from[k]);
  }
  return sb.data;
}

static ST_retcode xlsx_edit_file(XLSXExportContext *ctx,
                                 ctools_filtered_data *filtered) {
  mz_zip_archive in = {0}, out = {0};
  ST_retcode rc = 610;
  char *book = NULL, *rels = NULL, *types = NULL, *styles = NULL, *sheet = NULL,
       *fresh = NULL, *merged = NULL, *newshared = NULL;
  char *newbook = NULL, *newrels = NULL, *newtypes = NULL, *newstyles = NULL;
  char shared_path[MAX_PATH_LEN] = "xl/sharedStrings.xml";
  char target[MAX_PATH_LEN] = {0}, chain[MAX_PATH_LEN] = {0}, rid[256] = {0};
  /* A file that is not a workbook: native export reports that the file
   * exists (r(602)), or with sheet() that it cannot be loaded. */
  int named = ctx->sheet_action || ctx->sheet_specified;
  ST_retcode unreadable = named ? 603 : 602;
  snprintf(ctx->error_message, sizeof(ctx->error_message),
           named ? "file %.400s could not be loaded\n"
                 : "file %.400s already exists\n",
           ctx->input_filename);
  if (!mz_zip_reader_init_file(&in, ctx->input_filename, 0)) {
    rc = unreadable;
    goto done;
  }
  book = xedit_zip_text(&in, "xl/workbook.xml");
  rels = xedit_zip_text(&in, "xl/_rels/workbook.xml.rels");
  types = xedit_zip_text(&in, "[Content_Types].xml");
  styles = xedit_zip_text(&in, "xl/styles.xml");
  int had_styles = styles != NULL;
  if (!book || !rels || !types) {
    rc = unreadable;
    goto done;
  }
  ctx->error_message[0] = 0;
  if (ctx->invalid_sheet) {
    rc = 9901;
    snprintf(ctx->error_message, sizeof(ctx->error_message),
             "invalid worksheet name %.150s\n", ctx->sheet_name);
    goto done;
  }
  if (!styles)
    styles = strdup(STYLES_XML);
  if (!styles) {
    rc = 920;
    goto done;
  }
  xedit_node sheets, node;
  unsigned maxid = 0;
  if (!xedit_find(book, book + strlen(book), "sheets", &sheets))
    goto done;
  const char *p = sheets.open_end;
  int first = 1;
  while (!sheets.empty && xedit_find(p, sheets.close, "sheet", &node)) {
    char *name = xedit_attr(&node, "name"), *id = xedit_attr(&node, "sheetId");
    if (id && strtoul(id, NULL, 10) > maxid)
      maxid = (unsigned)strtoul(id, NULL, 10);
    /* Without sheet(), native export edits the first worksheet, whatever its
     * name; a later sheet named Sheet1 is not a second match. Chart, dialog
     * and macro sheets are no worksheets: skipped here, refused by name. */
    int named = name && ctx->sheet_specified && !strcasecmp(name, ctx->sheet_name);
    /* A name equal to another only ignoring the case of non-ASCII letters
     * (Données, DONNÉES) is neither that sheet nor a new one: Excel would
     * repair the workbook, and native export fails (r(9901)). */
    if (name && ctx->sheet_specified && !named && xedit_same_name(name, ctx->sheet_name)) {
      rc = 9901;
      snprintf(ctx->error_message, sizeof(ctx->error_message),
               "worksheet name %.150s differs from the existing worksheet "
               "%.150s only in case\n",
               ctx->sheet_name, name);
      free(name);
      free(id);
      goto done;
    }
    if (name && (named || (!ctx->sheet_specified && first))) {
      char *r = xedit_rel_id(&node);
      if (!r) {
        free(name);
        free(id);
        goto done;
      }
      if (xedit_rel_worksheet(rels, r)) {
        snprintf(rid, sizeof(rid), "%s", r);
        first = 0;
      } else if (named) {
        rc = 601;
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "sheet %.150s is a chart, dialog or macro sheet, not a "
                 "worksheet\n",
                 ctx->sheet_name);
        free(r);
        free(name);
        free(id);
        goto done;
      }
      free(r);
    }
    free(name);
    free(id);
    p = node.end;
  }
  if (rid[0] && !ctx->sheet_action) {
    rc = 602;
    snprintf(ctx->error_message, sizeof(ctx->error_message),
             "worksheet already exists; specify sheet(..., modify) or "
             "sheet(..., replace)\n");
    goto done;
  }
  if (rid[0]) {
    p = rels;
    const char *rels_end = rels + strlen(rels);
    while (xedit_find(p, rels_end, "Relationship", &node)) {
      char *id = xedit_attr(&node, "Id"), *t = xedit_attr(&node, "Target");
      if (id && t && !strcmp(id, rid)) {
        if (*t == '/')
          snprintf(target, sizeof(target), "%s", t + 1);
        else
          snprintf(target, sizeof(target), "xl/%s", t);
      }
      free(id);
      free(t);
      p = node.end;
    }
    if (!target[0])
      goto done;
    sheet = xedit_zip_text(&in, target);
    if (!sheet)
      goto done;
    xedit_node data;
    if (!xedit_find(sheet, sheet + strlen(sheet), "sheetData", &data)) {
      snprintf(ctx->error_message, sizeof(ctx->error_message),
               "worksheet part %.300s holds no sheetData\n", target);
      goto done;
    }
    char *positioned = xedit_explicit_positions(sheet);
    if (!positioned) {
      rc = 920;
      goto done;
    }
    if (positioned != sheet) {
      free(sheet);
      sheet = positioned;
    }
  } else {
    unsigned id = maxid + 1;
    do {
      snprintf(target, sizeof(target), "xl/worksheets/sheet%u.xml", id++);
    } while (mz_zip_reader_locate_file(&in, target, NULL, 0) >= 0);
    unsigned relation = 1;
    do {
      snprintf(rid, sizeof(rid), "ctoolsSheet%u", relation++);
    } while (strstr(rels, rid));
    char escaped[MAX_SHEET_NAME * 7 + 1], child[2048];
    xlsx_escape_xml(ctx->sheet_name, escaped, sizeof(escaped));
    snprintf(
        child, sizeof(child),
        "<sheet "
        "xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
        "xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/"
        "relationships\" "
        "name=\"%s\" sheetId=\"%u\" r:id=\"%s\"/>",
        escaped, maxid + 1, rid);
    newbook = xedit_insert(book, "sheets", child);
    snprintf(child, sizeof(child),
             "<Relationship "
             "xmlns=\"http://schemas.openxmlformats.org/package/2006/"
             "relationships\" Id=\"%s\" "
             "Type=\"http://schemas.openxmlformats.org/officeDocument/2006/"
             "relationships/"
             "worksheet\" Target=\"%s\"/>",
             rid, target + 3);
    newrels = xedit_insert(rels, "Relationships", child);
    snprintf(
        child, sizeof(child),
        "<Override "
        "xmlns=\"http://schemas.openxmlformats.org/package/2006/"
        "content-types\" "
        "PartName=\"/%s\" "
        "ContentType=\"application/"
        "vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>",
        target);
    newtypes = xedit_insert(types, "Types", child);
    if (!newbook || !newrels || !newtypes) {
      rc = 920;
      goto done;
    }
  }

  if (ctx->num_shared_strings) {
    const char *rp = newrels ? newrels : rels;
    const char *rend = rp+strlen(rp);
    int found=0;
    while (xedit_find(rp,rend,"Relationship",&node)) {
      char *type=xedit_attr(&node,"Type"), *target_s=xedit_attr(&node,"Target");
      if(type && target_s && strstr(type,"/sharedStrings")) {
        snprintf(shared_path,sizeof(shared_path),target_s[0]=='/' ? "%s" : "xl/%s", target_s[0]=='/' ? target_s+1 : target_s);
        found=1;
      }
      free(type);free(target_s);rp=node.end;
    }
    char *old=found ? xedit_zip_text(&in,shared_path) : strdup(SHARED_STRINGS_XML);
    if(!old || !xlsx_xml_validate_document(old,strlen(old),"sst")) {free(old);goto done;}
    /* The end is found once: rescanning it for every item is quadratic. */
    const char *sp=old, *old_end=old+strlen(old);
    while(xedit_find(sp,old_end,"si",&node)) {ctx->shared_offset++;sp=node.end;}
    char *items=xlsx_shared_items(ctx);
    if(!items) {free(old);rc=920;goto done;}
    newshared=xedit_insert(old,"sst",items);
    free(old);free(items);
    if(!newshared) {rc=920;goto done;}
    xedit_node root;
    if(!xedit_find(newshared,newshared+strlen(newshared),"sst",&root)) goto done;
    char *updated=xedit_set_attr(&root,"uniqueCount",ctx->shared_offset+ctx->num_shared_strings);
    if(!updated) {rc=920;goto done;}
    free(newshared);newshared=updated;
    /* The total reference count is optional. Remove a stale count rather
       than confusing unique strings with references across all worksheets. */
    if (xedit_find(newshared, newshared+strlen(newshared), "sst", &root)) {
      const char *a, *b, *value;
      size_t length;
      if (xedit_attr_span(&root, "count", &a, &b, &value, &length)) {
        StringBuffer clean;
        if (!strbuf_init(&clean, strlen(newshared)+1)) {rc=920;goto done;}
        strbuf_append_n(&clean, root.start, (size_t)(a-root.start));
        strbuf_append_n(&clean, b, (size_t)(root.end-b));
        if (!clean.data) {rc=920;goto done;}
        free(newshared);newshared=clean.data;
      }
    }
    if(!found) {
      char id[64];unsigned serial=1;
      do {snprintf(id,sizeof(id),"ctoolsStrings%u",serial++);} while(strstr(newrels ? newrels : rels,id));
      char relation[1024];
      snprintf(relation,sizeof(relation),"<Relationship xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\" Id=\"%s\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/sharedStrings\" Target=\"sharedStrings.xml\"/>",id);
      updated=xedit_insert(newrels ? newrels : rels,"Relationships",relation);
      if(!updated) {rc=920;goto done;} free(newrels);newrels=updated;
      updated=xedit_insert(newtypes ? newtypes : types,"Types","<Override xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\" PartName=\"/xl/sharedStrings.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml\"/>");
      if(!updated) {rc=920;goto done;} free(newtypes);newtypes=updated;
    }
  }
  /* Dates in 1904 workbooks use a different serial origin. */
  if (xedit_find(book, book + strlen(book), "workbookPr", &node)) {
    char *flag = xedit_attr(&node, "date1904");
    if (flag && (!strcmp(flag, "1") || !strcmp(flag, "true"))) {
      for (int j = 0; j < ctx->nvars; j++)
        if (CEXPORT_STYLE_IS_DATE(ctx->date_cols[j]))
          for (size_t i = 0; i < filtered->data.nobs; i++)
            if (!SF_is_missing(filtered->data.vars[j].data.dbl[i]))
              filtered->data.vars[j].data.dbl[i] -= 1462;
    }
    free(flag);
  }
  unsigned base;
  newstyles = xedit_styles(styles, &base);
  if (!newstyles) {
    rc = 920;
    goto done;
  }
  xlsx_worksheet_segments segs;
  if (!xlsx_gen_worksheet_segments(ctx, filtered, &segs)) {
    snprintf(ctx->error_message, sizeof(ctx->error_message),
             "Memory allocation failed for worksheet\n");
    rc = 920;
    goto done;
  }
  StringBuffer sb;
  if (!strbuf_init(&sb, 4096)) {
    xlsx_segments_free(&segs);
    rc = 920;
    goto done;
  }
  for (size_t i = 0; i < segs.num_segments; i++)
    strbuf_append_n(&sb, segs.segments[i].data, segs.segments[i].len);
  xlsx_segments_free(&segs);
  fresh = sb.data;
  if (!fresh) {
    rc = 920;
    goto done;
  }
  /* A replacement uses a clean sheet; keepcellfmt still maps styles by cell. */
  if (sheet && ctx->sheet_action == 1)
    merged = xedit_worksheet(sheet, fresh, ctx->keepcellfmt, base);
  else {
    static const char empty[] =
        "<worksheet "
        "xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/"
        "main\"><dimension ref=\"A1\"/><sheetData/></worksheet>";
    merged = xedit_worksheet(empty, fresh, 0, base);
  }
  if (!merged) {
    rc = 920;
    goto done;
  }
  if (sheet && ctx->sheet_action == 1) {
    char *updated = xedit_shared_anchors(sheet, merged);
    if (!updated) {
      rc = 920;
      goto done;
    }
    free(merged);
    merged = updated;
  }
  /* Our cell styles need a styles part that the workbook references. */
  if (!had_styles) {
    char *withrel = xedit_insert(
        newrels ? newrels : rels, "Relationships",
        "<Relationship "
        "xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\" "
        "Id=\"ctoolsStyles\" "
        "Type=\"http://schemas.openxmlformats.org/officeDocument/2006/"
        "relationships/styles\" Target=\"styles.xml\"/>");
    char *withtype = xedit_insert(
        newtypes ? newtypes : types, "Types",
        "<Override "
        "xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\" "
        "PartName=\"/xl/styles.xml\" "
        "ContentType=\"application/"
        "vnd.openxmlformats-officedocument.spreadsheetml.styles+xml\"/>");
    if (!withrel || !withtype) {
      free(withrel);
      free(withtype);
      rc = 920;
      goto done;
    }
    free(newrels);
    newrels = withrel;
    free(newtypes);
    newtypes = withtype;
  }
  char *cleanrels =
      xedit_remove_chain(newrels ? newrels : rels, "Relationship", "Type",
                         "/calcChain", chain, sizeof(chain));
  char *cleantypes =
      xedit_remove_chain(newtypes ? newtypes : types, "Override", "ContentType",
                         "calcChain+xml", NULL, 0);
  if (!cleanrels || !cleantypes) {
    free(cleanrels);
    free(cleantypes);
    rc = 920;
    goto done;
  }
  free(newrels);
  newrels = cleanrels;
  free(newtypes);
  newtypes = cleantypes;
  char *recalculate = xedit_recalculate(newbook ? newbook : book);
  if (!recalculate) {
    rc = 920;
    goto done;
  }
  free(newbook);
  newbook = recalculate;
  /* A Strict Open XML workbook: the parts written take its namespaces, or
   * injected rows and cells would leave the sheet's namespace. */
  if (strstr(book, "http://purl.oclc.org/ooxml/spreadsheetml/main")) {
    char **parts[] = {&newbook, &newrels, &newstyles, &merged, &newshared};
    for (int j = 0; j < 5; j++) {
      if (!*parts[j])
        continue;
      char *converted = xedit_strict(*parts[j]);
      if (!converted) {
        rc = 920;
        goto done;
      }
      free(*parts[j]);
      *parts[j] = converted;
    }
  }
  rc = 603;
  if (!mz_zip_writer_init_file(&out, ctx->filename, 0))
    goto done;
  const char *names[] = {"xl/workbook.xml", "xl/_rels/workbook.xml.rels",
                         "[Content_Types].xml", "xl/styles.xml", target};
  const char *values[] = {newbook ? newbook : book, newrels ? newrels : rels,
                          newtypes ? newtypes : types, newstyles, merged};
  for (unsigned i = 0; i < mz_zip_reader_get_num_files(&in); i++) {
    mz_zip_archive_file_stat stat;
    if (!mz_zip_reader_file_stat(&in, i, &stat))
      goto done;
    int changed = (*chain && !strcmp(stat.m_filename, chain)) || (newshared && !strcmp(stat.m_filename,shared_path));
    for (int j = 0; j < 5; j++)
      if (!strcmp(stat.m_filename, names[j]))
        changed = 1;
    if (!changed && !mz_zip_writer_add_from_zip_reader(&out, &in, i))
      goto done;
  }
  for (int j = 0; j < 5; j++)
    if (!mz_zip_writer_add_mem(&out, names[j], values[j], strlen(values[j]),
                               MZ_DEFAULT_COMPRESSION))
      goto done;
  if (newshared && !mz_zip_writer_add_mem(&out,shared_path,newshared,strlen(newshared),MZ_DEFAULT_COMPRESSION)) goto done;
  if (!mz_zip_writer_finalize_archive(&out))
    goto done;
  rc = 0;
done:
  /* Every failure says what failed: a missing or malformed part (r(610)),
   * memory (r(920)) or the output (r(603)). */
  if (rc && !ctx->error_message[0])
    snprintf(ctx->error_message, sizeof(ctx->error_message),
             rc == 920   ? "not enough memory to edit workbook %.400s\n"
             : rc == 603 ? "could not write the edited workbook %.400s\n"
                         : "workbook %.400s has a missing or malformed part "
                           "that cexport excel cannot edit\n",
             ctx->input_filename);
  if (out.m_pState)
    mz_zip_writer_end(&out);
  if (in.m_pState)
    mz_zip_reader_end(&in);
  free(book);
  free(rels);
  free(types);
  free(styles);
  free(sheet);
  free(fresh);
  free(merged);
  free(newbook);
  free(newrels);
  free(newtypes);
  free(newstyles);
  free(newshared);
  return rc;
}

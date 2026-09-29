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
  if (keep && old) {
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
/* Append our three cell styles while retaining every existing style index. */
static char *xedit_styles(char *xml, unsigned *base) {
  xedit_node xfs, xf, nfmt;
  if (!xedit_find(xml, xml + strlen(xml), "cellXfs", &xfs))
    return NULL;
  *base = 0;
  const char *p = xfs.open_end;
  while (!xfs.empty && xedit_find(p, xfs.close, "xf", &xf)) {
    (*base)++;
    p = xf.end;
  }
  unsigned fmt = 164;
  p = xml;
  while (xedit_find(p, xml + strlen(xml), "numFmt", &nfmt)) {
    char *id = xedit_attr(&nfmt, "numFmtId");
    if (id && strtoul(id, NULL, 10) >= fmt)
      fmt = (unsigned)strtoul(id, NULL, 10) + 1;
    free(id);
    p = nfmt.end;
  }
  char child[512];
  snprintf(
      child, sizeof(child),
      "<numFmt "
      "xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
      "numFmtId=\"%u\" formatCode=\"yyyy-mm-dd hh:mm:ss.000\"/>",
      fmt);
  char *result;
  if (xedit_find(xml, xml + strlen(xml), "numFmts", &nfmt)) {
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
  xedit_find(result, result + strlen(result), "cellXfs", &xfs);
  char *node = xedit_set_attr(&xfs, "count", *base + 3);
  if (!node) {
    free(result);
    return NULL;
  }
  char *updated = xedit_replace(result, &xfs, node);
  free(node);
  free(result);
  if (!updated)
    return NULL;
  snprintf(
      child, sizeof(child),
      "<xf xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
      "numFmtId=\"0\" fontId=\"0\" fillId=\"0\" borderId=\"0\" xfId=\"0\"/>"
      "<xf xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
      "numFmtId=\"14\" fontId=\"0\" fillId=\"0\" borderId=\"0\" xfId=\"0\"/>"
      "<xf xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
      "numFmtId=\"%u\" fontId=\"0\" fillId=\"0\" borderId=\"0\" xfId=\"0\"/>",
      fmt);
  result = xedit_insert(updated, "cellXfs", child);
  free(updated);
  return result;
}

static ST_retcode xlsx_edit_file(XLSXExportContext *ctx,
                                 ctools_filtered_data *filtered) {
  mz_zip_archive in = {0}, out = {0};
  ST_retcode rc = 610;
  char *book = NULL, *rels = NULL, *types = NULL, *styles = NULL, *sheet = NULL,
       *fresh = NULL, *merged = NULL;
  char *newbook = NULL, *newrels = NULL, *newtypes = NULL, *newstyles = NULL;
  char target[MAX_PATH_LEN] = {0}, chain[MAX_PATH_LEN] = {0}, rid[256] = {0};
  if (!mz_zip_reader_init_file(&in, ctx->input_filename, 0))
    goto done;
  book = xedit_zip_text(&in, "xl/workbook.xml");
  rels = xedit_zip_text(&in, "xl/_rels/workbook.xml.rels");
  types = xedit_zip_text(&in, "[Content_Types].xml");
  styles = xedit_zip_text(&in, "xl/styles.xml");
  if (!book || !rels || !types)
    goto done;
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
    if (name && ((!ctx->sheet_specified && first) ||
                 !strcasecmp(name, ctx->sheet_name))) {
      char *r = xedit_attr(&node, "r:id");
      if (!r) {
        free(name);
        free(id);
        goto done;
      }
      snprintf(rid, sizeof(rid), "%s", r);
      free(r);
    }
    free(name);
    free(id);
    p = node.end;
    first = 0;
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
    while (xedit_find(p, rels + strlen(rels), "Relationship", &node)) {
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
  } else {
    unsigned id = maxid + 1;
    do {
      snprintf(target, sizeof(target), "xl/worksheets/sheet%u.xml", id++);
    } while (mz_zip_reader_locate_file(&in, target, NULL, 0) >= 0);
    unsigned relation = 1;
    do {
      snprintf(rid, sizeof(rid), "ctoolsSheet%u", relation++);
    } while (strstr(rels, rid));
    char escaped[MAX_SHEET_NAME * 6 + 1], child[2048];
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
  /* Dates in 1904 workbooks use a different serial origin. */
  if (xedit_find(book, book + strlen(book), "workbookPr", &node)) {
    char *flag = xedit_attr(&node, "date1904");
    if (flag && (!strcmp(flag, "1") || !strcmp(flag, "true"))) {
      for (int j = 0; j < ctx->nvars; j++)
        if (ctx->date_cols[j])
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
    int changed = *chain && !strcmp(stat.m_filename, chain);
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
  if (!mz_zip_writer_finalize_archive(&out))
    goto done;
  rc = 0;
done:
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
  return rc;
}

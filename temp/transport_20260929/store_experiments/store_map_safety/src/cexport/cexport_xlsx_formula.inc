/* Move a shared-formula anchor after its original cell is overwritten. */
static int xedit_namechar(unsigned char c) {
  return isalnum(c) || c == '_' || c == '.' || c >= 128;
}
static char *xedit_translate_formula(const char *formula, int dc, int dr) {
  StringBuffer out;
  if (!strbuf_init(&out, strlen(formula) + 64))
    return NULL;
  const char *p = formula;
  while (*p) {
    /* Strings, quoted worksheet names and structured-reference columns. */
    if (*p == '"' || *p == '\'' || *p == '[') {
      char quote = *p == '[' ? ']' : *p;
      const char *start = p++;
      int depth = 1;
      while (*p && depth) {
        if (quote == ']' && *p == '[')
          depth++;
        if (*p++ == quote) {
          if (quote != ']' && *p == quote) {
            p++;
            continue;
          }
          depth--;
        }
      }
      strbuf_append_n(&out, start, (size_t)(p - start));
      continue;
    }
    const char *start = p, *q = p;
    /* Unquoted sheet prefixes (including 3D ranges) are identifiers. */
    const char *sheet = p;
    while (xedit_namechar((unsigned char)*sheet) || *sheet == '$' ||
           *sheet == ':')
      sheet++;
    if (sheet > p && *sheet == '!') {
      strbuf_append_n(&out, p, (size_t)(sheet + 1 - p));
      p = sheet + 1;
      continue;
    }
    int fixed_col = *q == '$';
    if (fixed_col)
      q++;
    const char *letters = q;
    int col = 0;
    while (isalpha((unsigned char)*q) && q - letters < 3) {
      col = col * 26 + toupper((unsigned char)*q++) - 'A' + 1;
    }
    int fixed_row = *q == '$';
    if (fixed_row)
      q++;
    const char *digits = q;
    int row = 0;
    while (isdigit((unsigned char)*q) && q - digits < 7)
      row = row * 10 + *q++ - '0';
    int allowed_before =
        start == formula || !xedit_namechar((unsigned char)start[-1]);
    /* Native shifts A1 cells; whole-row/column ranges retain their spelling. */
    if (allowed_before && q > digits && letters < digits && col >= 1 &&
        col <= 16384 && row >= 1 && row <= 1048576 &&
        !xedit_namechar((unsigned char)*q) && *q != '(') {
      int newcol = fixed_col ? col : col + dc,
          newrow = fixed_row ? row : row + dr;
      if (newcol < 1 || newcol > 16384 || newrow < 1 || newrow > 1048576)
        strbuf_append(&out, "#REF!");
      else {
        char column[16];
        xlsx_col_to_letters(newcol, column);
        if (islower((unsigned char)*letters))
          for (char *s = column; *s; s++)
            *s = (char)tolower((unsigned char)*s);
        strbuf_appendf(&out, "%s%s%s%d", fixed_col ? "$" : "", column,
                       fixed_row ? "$" : "", newrow);
      }
      p = q;
    } else
      strbuf_append_n(&out, p++, 1);
  }
  return out.data;
}
static int xedit_shared_id(const xedit_node *f, const char *id) {
  const char *a, *b, *v;
  size_t length;
  return xedit_attr_span(f, "si", &a, &b, &v, &length) &&
         strlen(id) == length && !memcmp(id, v, length);
}
static char *xedit_shared_anchors(const char *oldxml, const char *merged) {
  char *result = strdup(merged);
  if (!result)
    return NULL;
  const char *p = oldxml, *end = p + strlen(p);
  xedit_node cell, anchor;
  while (xedit_find(p, end, "c", &cell)) {
    p = cell.end;
    if (!xedit_find(cell.open_end, cell.close, "f", &anchor) || anchor.empty ||
        anchor.close == anchor.open_end)
      continue;
    char *type = xedit_attr(&anchor, "t"), *id = xedit_attr(&anchor, "si");
    int shared = type && !strcmp(type, "shared") && id;
    free(type);
    if (!shared) {
      free(id);
      continue;
    }
    int col, row;
    if (!xedit_cell(&cell, &col, &row)) {
      free(id);
      continue;
    }
    const char *q = result, *limit = q + strlen(q);
    xedit_node next, f, first = {0};
    int nc = 0, nr = 0, surviving_anchor = 0;
    while (xedit_find(q, limit, "c", &next)) {
      q = next.end;
      if (!xedit_find(next.open_end, next.close, "f", &f) ||
          !xedit_shared_id(&f, id))
        continue;
      if (!f.empty && f.close > f.open_end) {
        surviving_anchor = 1;
        break;
      }
      if (!first.start && xedit_cell(&next, &nc, &nr))
        first = f;
    }
    free(id);
    if (surviving_anchor || !first.start)
      continue;
    size_t length = (size_t)(anchor.close - anchor.open_end);
    char *formula = malloc(length + 1);
    if (!formula)
      goto failed;
    memcpy(formula, anchor.open_end, length);
    formula[length] = 0;
    length = xlsx_xml_decode_entities(formula, length);
    formula[length] = 0;
    char *translated = xedit_translate_formula(formula, nc - col, nr - row);
    free(formula);
    if (!translated)
      goto failed;
    length = strlen(translated);
    if (length > (SIZE_MAX - 1) / 6) {
      free(translated);
      goto failed;
    }
    char *escaped = malloc(length * 6 + 1);
    if (!escaped) {
      free(translated);
      goto failed;
    }
    xlsx_escape_xml(translated, escaped, length * 6 + 1);
    free(translated);
    StringBuffer node;
    if (!strbuf_init(&node, length * 6 + 256)) {
      free(escaped);
      goto failed;
    }
    strbuf_append_n(&node, anchor.start,
                    (size_t)(anchor.open_end - anchor.start));
    strbuf_append(&node, escaped);
    free(escaped);
    strbuf_append_n(&node, anchor.close, (size_t)(anchor.end - anchor.close));
    char *updated = node.data ? xedit_replace(result, &first, node.data) : NULL;
    free(node.data);
    if (!updated)
      goto failed;
    free(result);
    result = updated;
  }
  return result;
failed:
  free(result);
  return NULL;
}
/* Remove calculation-chain declarations, returning the related ZIP part. */
static char *xedit_remove_chain(const char *xml, const char *element,
                                const char *attribute, const char *suffix,
                                char *target, size_t capacity) {
  StringBuffer out;
  if (!strbuf_init(&out, strlen(xml) + 1))
    return NULL;
  const char *p = xml, *keep = xml, *end = xml + strlen(xml);
  xedit_node node;
  while (xedit_find(p, end, element, &node)) {
    p = node.end;
    char *value = xedit_attr(&node, attribute);
    size_t n = value ? strlen(value) : 0, k = strlen(suffix);
    int remove = n >= k && !strcmp(value + n - k, suffix);
    free(value);
    if (!remove)
      continue;
    if (target) {
      char *part = xedit_attr(&node, "Target");
      if (!part) {
        free(out.data);
        return NULL;
      }
      snprintf(target, capacity, "%s%s", *part == '/' ? "" : "xl/",
               *part == '/' ? part + 1 : part);
      free(part);
    }
    strbuf_append_n(&out, keep, (size_t)(node.start - keep));
    keep = node.end;
  }
  strbuf_append(&out, keep);
  return out.data;
}

static char *xedit_recalculate(const char *book) {
  xedit_node calc;
  if (!xedit_find(book, book + strlen(book), "calcPr", &calc))
    return xedit_insert(
        book, "workbook",
        "<calcPr "
        "xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
        "calcId=\"124519\" fullCalcOnLoad=\"1\"/>");
  char *node = xedit_set_attr(&calc, "fullCalcOnLoad", 1);
  if (!node)
    return NULL;
  char *result = xedit_replace(book, &calc, node);
  free(node);
  return result;
}

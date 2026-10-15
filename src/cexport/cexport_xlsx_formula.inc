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
/* The first qualifying old anchor of each si, and what the merged sheet holds
 * for that si: a surviving anchor, or the first dependent with a cell ref. */
typedef struct {
  xedit_node anchor, first; /* anchor in oldxml; dependent in merged */
  char *id;                 /* decoded si */
  int col, row, nc, nr, survives;
} xedit_shared_group;
typedef struct {
  const char *start, *end;
  char *text;
} xedit_shared_edit;
static int xedit_shared_cmp(const void *x, const void *y) {
  const xedit_shared_group *a = x, *b = y;
  int c = strcmp(a->id, b->id);
  return c ? c : (a->anchor.start > b->anchor.start) -
                     (a->anchor.start < b->anchor.start);
}
static int xedit_shared_edit_cmp(const void *x, const void *y) {
  const xedit_shared_edit *a = x, *b = y;
  return (a->start > b->start) - (a->start < b->start);
}
/* Merged cells match a decoded id by their raw si bytes. */
static xedit_shared_group *xedit_shared_lookup(xedit_shared_group *groups,
                                               size_t count, const char *v,
                                               size_t length) {
  size_t lo = 0, hi = count;
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2, n = strlen(groups[mid].id);
    int c = memcmp(groups[mid].id, v, n < length ? n : length);
    if (!c)
      c = (n > length) - (n < length);
    if (!c)
      return &groups[mid];
    if (c < 0)
      lo = mid + 1;
    else
      hi = mid;
  }
  return NULL;
}
/* One pass over each sheet, then one rebuild: a relocation makes its si's
 * anchor survive, so only the first old anchor of an si can move. */
static char *xedit_shared_anchors(const char *oldxml, const char *merged) {
  const char *mend = merged + strlen(merged);
  xedit_shared_group *groups = NULL;
  xedit_shared_edit *edits = NULL;
  size_t ngroups = 0, nedits = 0, cap = 0, ecap = 0, extra = 0;
  char *result = NULL;
  const char *p = oldxml, *end = p + strlen(p);
  xedit_node cell, anchor, f;
  while (xedit_find(p, end, "c", &cell)) {
    p = cell.end;
    if (!xedit_find(cell.open_end, cell.close, "f", &anchor) || anchor.empty ||
        anchor.close == anchor.open_end)
      continue;
    char *type = xedit_attr(&anchor, "t"), *id = xedit_attr(&anchor, "si");
    int shared = type && !strcmp(type, "shared") && id;
    free(type);
    int col, row;
    if (!shared || !xedit_cell(&cell, &col, &row)) {
      free(id);
      continue;
    }
    if (ngroups == cap) {
      size_t grown = cap ? cap * 2 : 16;
      xedit_shared_group *g = realloc(groups, grown * sizeof(*groups));
      if (!g) {
        free(id);
        goto done;
      }
      groups = g;
      cap = grown;
    }
    groups[ngroups++] =
        (xedit_shared_group){anchor, {0}, id, col, row, 0, 0, 0};
  }
  if (ngroups)
    qsort(groups, ngroups, sizeof(*groups), xedit_shared_cmp);
  size_t kept = 0;
  for (size_t g = 0; g < ngroups; g++) {
    if (kept && !strcmp(groups[kept - 1].id, groups[g].id))
      free(groups[g].id);
    else
      groups[kept++] = groups[g];
  }
  ngroups = kept;
  const char *q = ngroups ? merged : mend;
  while (xedit_find(q, mend, "c", &cell)) {
    q = cell.end;
    const char *a, *b, *v;
    size_t length;
    if (!xedit_find(cell.open_end, cell.close, "f", &f) ||
        !xedit_attr_span(&f, "si", &a, &b, &v, &length))
      continue;
    xedit_shared_group *g = xedit_shared_lookup(groups, ngroups, v, length);
    if (!g || g->survives)
      continue;
    if (!f.empty && f.close > f.open_end)
      g->survives = 1;
    else if (!g->first.start && xedit_cell(&cell, &g->nc, &g->nr))
      g->first = f;
  }
  for (size_t g = 0; g < ngroups; g++) {
    const xedit_shared_group *s = &groups[g];
    if (s->survives || !s->first.start)
      continue;
    size_t length = (size_t)(s->anchor.close - s->anchor.open_end);
    char *formula = malloc(length + 1);
    if (!formula)
      goto done;
    memcpy(formula, s->anchor.open_end, length);
    formula[length] = 0;
    length = xlsx_xml_decode_entities(formula, length);
    formula[length] = 0;
    char *translated =
        xedit_translate_formula(formula, s->nc - s->col, s->nr - s->row);
    free(formula);
    if (!translated)
      goto done;
    length = strlen(translated);
    if (length > (SIZE_MAX - 1) / 6) {
      free(translated);
      goto done;
    }
    char *escaped = malloc(length * 6 + 1);
    if (!escaped) {
      free(translated);
      goto done;
    }
    xlsx_escape_xml(translated, escaped, length * 6 + 1);
    free(translated);
    StringBuffer node;
    if (!strbuf_init(&node, length * 6 + 256)) {
      free(escaped);
      goto done;
    }
    strbuf_append_n(&node, s->anchor.start,
                    (size_t)(s->anchor.open_end - s->anchor.start));
    strbuf_append(&node, escaped);
    free(escaped);
    strbuf_append_n(&node, s->anchor.close,
                    (size_t)(s->anchor.end - s->anchor.close));
    if (!node.data)
      goto done;
    if (nedits == ecap) {
      size_t grown = ecap ? ecap * 2 : 16;
      xedit_shared_edit *e = realloc(edits, grown * sizeof(*edits));
      if (!e) {
        free(node.data);
        goto done;
      }
      edits = e;
      ecap = grown;
    }
    edits[nedits++] =
        (xedit_shared_edit){s->first.start, s->first.end, node.data};
    extra += node.len;
  }
  if (nedits)
    qsort(edits, nedits, sizeof(*edits), xedit_shared_edit_cmp);
  StringBuffer out;
  if (!strbuf_init(&out, (size_t)(mend - merged) + extra + 64))
    goto done;
  const char *keep = merged;
  for (size_t e = 0; e < nedits; e++) {
    strbuf_append_n(&out, keep, (size_t)(edits[e].start - keep));
    strbuf_append(&out, edits[e].text);
    keep = edits[e].end;
  }
  strbuf_append(&out, keep);
  result = out.data; /* NULL after a failed grow */
done:
  for (size_t g = 0; g < ngroups; g++)
    free(groups[g].id);
  free(groups);
  for (size_t e = 0; e < nedits; e++)
    free(edits[e].text);
  free(edits);
  return result;
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

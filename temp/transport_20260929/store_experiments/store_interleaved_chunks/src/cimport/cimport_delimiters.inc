/* Normalize native delimiter sets/literal strings into one unused byte, after
   encoding conversion. Quoted contents remain byte exact; the parallel parser
   and column caches can then use their existing single-byte fast path. */
static int cimport_normalize_delimiters(CImportContext *ctx) {
  char delimiters[4096] = "", collapse_text[8] = "", literal_text[8] = "";
  SF_macro_use("___cimport_delimiters", delimiters, sizeof(delimiters) - 1);
  SF_macro_use("___cimport_collapse", collapse_text, sizeof(collapse_text) - 1);
  SF_macro_use("___cimport_asstring", literal_text, sizeof(literal_text) - 1);
  int collapse = !strcmp(collapse_text, "1"),
      literal = !strcmp(literal_text, "1");
  size_t width = strlen(delimiters);
  if (!width || (width == 1 && !collapse))
    return 0;
  unsigned char used[256] = {0};
  for (size_t i = 0; i < ctx->file_size; i++)
    used[(unsigned char)ctx->file_data[i]] = 1;
  unsigned char separator = 1;
  while (separator < 32 && (used[separator] || separator == '\r' ||
                            separator == '\n' || separator == '\t'))
    separator++;
  if (separator == 32) {
    if (used[127])
      return 198;
    separator = 127;
  }
  char *output = malloc(ctx->file_size + 1);
  if (!output)
    return 920;
  size_t out = 0;
  int quoted = 0, previous = 0;
  for (size_t i = 0; i < ctx->file_size;) {
    char ch = ctx->file_data[i];
    if (ch == '"' && ctx->bindquotes != CIMPORT_BINDQUOTES_NOBIND) {
      if (quoted && i + 1 < ctx->file_size && ctx->file_data[i + 1] == '"') {
        output[out++] = ch;
        output[out++] = ch;
        i += 2;
        previous = 0;
        continue;
      }
      quoted = !quoted;
    }
    if (ch == '\n' && ctx->bindquotes != CIMPORT_BINDQUOTES_STRICT)
      quoted = 0;
    size_t match = 0;
    if (!quoted && ch != '\n' && ch != '\r') {
      if (literal) {
        if (width <= ctx->file_size - i &&
            !memcmp(ctx->file_data + i, delimiters, width))
          match = width;
      } else {
        for (size_t d = 0; d < width;) {
          unsigned char first = (unsigned char)delimiters[d];
          size_t chars = first < 0x80   ? 1
                         : first < 0xe0 ? 2
                         : first < 0xf0 ? 3
                                        : 4;
          if (chars > width - d)
            chars = 1;
          if (chars <= ctx->file_size - i &&
              !memcmp(ctx->file_data + i, delimiters + d, chars)) {
            match = chars;
            break;
          }
          d += chars;
        }
      }
    }
    if (match) {
      if (!collapse || !previous)
        output[out++] = (char)separator;
      i += match;
      previous = 1;
    } else {
      output[out++] = ch;
      i++;
      previous = 0;
    }
  }
  output[out] = 0;
  free(ctx->converted_data);
  ctx->converted_data = output;
  ctx->converted_size = out;
  ctx->file_data = output;
  ctx->file_size = out;
  ctx->delimiter = (char)separator;
  return 0;
}

static int cimport_check_quoted_rows(const CImportContext *ctx) {
  if (ctx->bindquotes != CIMPORT_BINDQUOTES_STRICT || ctx->max_quoted_rows <= 0)
    return 0;
  int quoted = 0, lines = 1;
  for (size_t i = 0; i < ctx->file_size; i++) {
    char ch = ctx->file_data[i];
    if (ch == '"') {
      if (quoted && i + 1 < ctx->file_size && ctx->file_data[i + 1] == '"') {
        i++;
        continue;
      }
      quoted = !quoted;
      lines = 1;
    } else if (ch == '\n' && quoted && ++lines > ctx->max_quoted_rows) {
      cimport_display_error("Unmatched quote exceeded maxquotedrows()\n");
      return 5101;
    }
  }
  return 0;
}

/* File offset where a parsed row ends: its line terminator, or EOF. */
static uint64_t cimport_row_end(const CImportParsedRow *row) {
  uint64_t end = row->offset;
  for (int f = 0; f < row->num_fields; f++)
    end += (row->lengths[f] & 0x7FFFFFFFu) + 1;
  return row->num_fields ? end - 1 : end;
}

/* Select rows before inference so discarded text/outliers cannot determine the
   retained columns' storage. Keep the header reference for names and labels.
   rowrange() arrives as file lines and is applied here, after the header
   decision, the way native reads it: every physical line counts (the header,
   blank lines, strict quoted newlines), and a row is kept when it ends on or
   after the first line and the row before it (or the header) ended before the
   last line. */
static void cimport_select_rows(CImportContext *ctx) {
  size_t first, last;
  cimport_row_range(&first, &last);
  if (!first && !last)
    return;
  char auto_text[8] = "";
  SF_macro_use("___cimport_autoheader", auto_text, sizeof(auto_text) - 1);
  CImportLineCursor cursor = {ctx->file_data, 0};
  size_t prev = 0, data_rows = 0, total = 0;
  for (int c = 0; c < ctx->num_chunks; c++) {
    CImportParsedChunk *chunk = &ctx->chunks[c];
    size_t header = c == 0 && ctx->has_header ? 1 : 0, write = header;
    cimport_init_col_stats(chunk->col_stats, chunk->num_col_stats);
    /* Nothing is read once a row has ended on or after the last line. */
    for (size_t r = 0; r < chunk->num_rows && !(last && prev >= last); r++) {
      size_t line = cimport_line_at(
          &cursor, ctx->file_data + cimport_row_end(chunk->rows[r]));
      int selected = r >= header && cimport_line_in_range(first, last, prev, line);
      /* An explicit header seeds native inference from the first data row. */
      int seed = !selected && r >= header && data_rows == 0 &&
                 ctx->has_header && strcmp(auto_text, "1");
      prev = line;
      if (r < header)
        continue;
      data_rows++;
      if (!selected && !seed)
        continue;
      CImportParsedRow *row = chunk->rows[r];
      uint64_t *values = cimport_row_values(row);
      CImportFieldRef field = {row->offset, 0};
      if (selected)
        chunk->rows[write++] = row;
      for (int f = 0; f < row->num_fields; f++) {
        field.length = row->lengths[f];
        if (f < chunk->num_col_stats) {
          CImportColumnParseStats *stats = &chunk->col_stats[f];
          int len = (int)CIMPORT_FIELD_LENGTH(field);
          if (selected) {
            values[f] = cimport_accumulate_field(ctx, stats, &field, ctx->col_force[f]);
          } else {
            /* The seed row fixes type and width, not numeric extrema. */
            if (len > stats->max_field_len)
              stats->max_field_len = len;
            ctx->columns[f].inference_width = cimport_extract_field_fast(
                ctx->file_data, &field, NULL, INT_MAX, cimport_extract_quote(ctx));
            if (field.length & CIMPORT_FIELD_QUOTED_FLAG)
              stats->has_quotes = true;
            if (len && !cimport_context_numeric(ctx, &field))
              stats->seen_string = true;
          }
        }
        field.offset += CIMPORT_FIELD_LENGTH(field) + 1;
      }
      if (selected)
        total++;
    }
    chunk->num_rows = write;
  }
  ctx->total_rows = total;
}

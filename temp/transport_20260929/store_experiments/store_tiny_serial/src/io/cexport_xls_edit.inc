/* Retain the original compound-file streams and edit only Workbook's BIFF data.
 */
#include "vendor/libxls/include/libxls/xls_endian.h"
#include "vendor/libxls/include/xls.h"
static unsigned xe16(const unsigned char *p) {
  return p[0] | ((unsigned)p[1] << 8);
}
static uint32_t xe32(const unsigned char *p) {
  return (uint32_t)xe16(p) | ((uint32_t)xe16(p + 2) << 16);
}
typedef struct {
  unsigned id, len;
  size_t start, end;
  const unsigned char *data;
} xe_record;
static int xe_next(const unsigned char *data, size_t len, size_t *pos,
                   xe_record *r) {
  if (*pos > len || len - *pos < 4)
    return 0;
  unsigned n = xe16(data + *pos + 2);
  if (n > 8224 || n > len - *pos - 4)
    return 0;
  *r = (xe_record){xe16(data + *pos), n, *pos, *pos + 4 + n, data + *pos + 4};
  *pos = r->end;
  return 1;
}
static int xe_copy(FILE *f, const xe_record *r) {
  return record(f, r->id, r->data, r->len);
}
static int xe_cell_id(unsigned id) {
  return id == 0x203 || id == 0xfd || id == 0x204 || id == 0x27e ||
         id == 0x201 || id == 0x205 || id == 0x6 || id == 0xd6;
}
static int xe_inside(unsigned row, unsigned col, int r0, int c0, size_t n,
                     int k) {
  return row >= (unsigned)r0 && (size_t)row < (size_t)r0 + n &&
         col >= (unsigned)c0 && col < (unsigned)(c0 + k);
}
/* Write a cell from the generated sheet with its SST/style indices relocated.
 */
static int xe_new_cell(FILE *f, const xe_record *r, unsigned strings,
                       unsigned styles, xlsWorkSheet *old, int keep,
                       int date1904) {
  if (r->len < 6)
    return 610;
  unsigned char b[8224];
  memcpy(b, r->data, r->len);
  unsigned style = xe16(b + 4) + styles;
  if (keep && old) {
    xlsCell *cell = xls_cell(old, (XLSWORD)xe16(b), (XLSWORD)xe16(b + 2));
    if (cell && cell->id)
      style = cell->xf;
  }
  if (date1904 && r->id == 0x203 && r->len >= 14 && xe16(b + 4) > 0) {
    uint64_t bits = (uint64_t)xe32(b + 6) | ((uint64_t)xe32(b + 10) << 32);
    double value;
    memcpy(&value, &bits, 8);
    value -= 1462;
    memcpy(&bits, &value, 8);
    put64(b + 6, bits);
  }
  put16(b + 4, style);
  if (r->id == 0xfd && r->len >= 10)
    put32(b + 6, xe32(b + 6) + strings);
  return record(f, r->id, b, r->len);
}
static unsigned char *xe_file(FILE *f, size_t *n) {
  long end = ftell(f);
  if (end < 0 || fseek(f, 0, SEEK_SET))
    return NULL;
  unsigned char *data = malloc((size_t)end);
  if (!data)
    return NULL;
  if (fread(data, 1, (size_t)end, f) != (size_t)end) {
    free(data);
    return NULL;
  }
  *n = (size_t)end;
  return data;
}
static int xe_globals(FILE *f, const unsigned char *old, size_t oldlen,
                      const unsigned char *fresh, size_t freshlen, int add,
                      long **bounds, unsigned *nbounds, unsigned *strings,
                      unsigned *styles, size_t *newsheet) {
  size_t p = 0, np = 0;
  xe_record r, nr, sst = {0}, bound = {0};
  *strings = 0;
  *styles = 0;
  *nbounds = 0;
  *bounds = NULL;
  *newsheet = 0;
  while (xe_next(fresh, freshlen, &np, &nr)) {
    if (nr.id == 0xfc)
      sst = nr;
    if (nr.id == 0x85)
      bound = nr;
    if (nr.id == 0x0a) {
      *newsheet = nr.end;
      break;
    }
  }
  if (!sst.start || sst.len < 8 || !bound.start || !*newsheet)
    return 610;
  int rc = 0, wrote_sst = 0;
  while (xe_next(old, oldlen, &p, &r)) {
    if (r.id == 0x0a)
      break;
    if (r.id == 0xe0)
      (*styles)++;
    if (r.id == 0x85) {
      long *tmp = realloc(*bounds, ((size_t)*nbounds + 1) * sizeof(long));
      if (!tmp)
        return 909;
      *bounds = tmp;
      (*bounds)[(*nbounds)++] = ftell(f) + 4;
    }
    if (r.id == 0xff)
      continue; /* EXTSST indexes change when appending strings. */
    if (r.id == 0xfc) {
      if (r.len < 8)
        return 610;
      unsigned char b[8224];
      memcpy(b, r.data, r.len);
      *strings = xe32(b + 4);
      put32(b, xe32(b) + xe32(sst.data));
      put32(b + 4, *strings + xe32(sst.data + 4));
      if ((rc = record(f, r.id, b, r.len)))
        return rc;
      /* Preserve original SST continuations; the new strings start at a
         string boundary, so their first continuation needs no flag byte. */
      size_t next = p;
      xe_record continuation;
      while (xe_next(old, oldlen, &next, &continuation) &&
             continuation.id == 0x3c) {
        if ((rc = xe_copy(f, &continuation)))
          return rc;
        p = next;
      }
      if ((rc = record(f, 0x3c, sst.data + 8, sst.len - 8)))
        return rc;
      next = sst.end;
      while (xe_next(fresh, freshlen, &next, &continuation) &&
             continuation.id == 0x3c)
        if ((rc = xe_copy(f, &continuation)))
          return rc;
      wrote_sst = 1;
      continue;
    }
    if ((rc = xe_copy(f, &r)))
      return rc;
  }
  if (r.id != 0x0a)
    return 610;
  np = 0;
  while (xe_next(fresh, freshlen, &np, &nr) && nr.id != 0x0a) {
    if (nr.id == 0xe0 || (!wrote_sst && (nr.id == 0xfc || nr.id == 0x3c)))
      if ((rc = xe_copy(f, &nr)))
        return rc;
  }
  if (add) {
    long *tmp = realloc(*bounds, ((size_t)*nbounds + 1) * sizeof(long));
    if (!tmp)
      return 909;
    *bounds = tmp;
    (*bounds)[(*nbounds)++] = ftell(f) + 4;
    if ((rc = xe_copy(f, &bound)))
      return rc;
  }
  return record(f, 0x0a, NULL, 0);
}
static int xe_sheet(FILE *f, const unsigned char *old, size_t oldlen,
                    size_t start, const unsigned char *fresh, size_t freshlen,
                    size_t newstart, unsigned strings, unsigned styles,
                    int action, int keep, int r0, int c0, size_t n, int k,
                    xlsWorkSheet *ws) {
  size_t p = start;
  xe_record r;
  int rc, skip_string = 0;
  while (xe_next(old, oldlen, &p, &r)) {
    if (r.id == 0x0a)
      break;
    if (r.id == 0x20b || r.id == 0xd7)
      continue; /* INDEX/DBCELL byte offsets */
    if (r.id == 0x200 && r.len >= 14) {
      unsigned char b[14];
      memcpy(b, r.data, 14);
      if ((uint32_t)r0 < xe32(b))
        put32(b, (uint32_t)r0);
      if ((size_t)r0 + n > xe32(b + 4))
        put32(b + 4, (uint32_t)((size_t)r0 + n));
      if (c0 < (int)xe16(b + 8))
        put16(b + 8, (unsigned)c0);
      if (c0 + k > (int)xe16(b + 10))
        put16(b + 10, (unsigned)(c0 + k));
      if ((rc = record(f, r.id, b, 14)))
        return rc;
      continue;
    }
    if (skip_string && (r.id == 0x207 || r.id == 0x3c))
      continue;
    skip_string = 0;
    if (xe_cell_id(r.id) && r.len >= 6 &&
        xe_inside(xe16(r.data), xe16(r.data + 2), r0, c0, n, k)) {
      if (r.id == 6)
        skip_string = 1;
      continue;
    }
    if ((r.id == 0xbd || r.id == 0xbe) && r.len >= 6) {
      unsigned row = xe16(r.data), first = xe16(r.data + 2),
               last = xe16(r.data + r.len - 2);
      size_t width = r.id == 0xbd ? 6 : 2;
      if (last < first || (size_t)(last - first + 1) * width != r.len - 6)
        return 610;
      for (unsigned col = first; col <= last; col++) {
        if (xe_inside(row, col, r0, c0, n, k))
          continue;
        unsigned char b[10];
        put16(b, row);
        put16(b + 2, col);
        memcpy(b + 4, r.data + 4 + (col - first) * width, width);
        if ((rc = record(f, r.id == 0xbd ? 0x27e : 0x201, b,
                         (unsigned)(4 + width))))
          return rc;
      }
      continue;
    }
    if ((rc = xe_copy(f, &r)))
      return rc;
  }
  if (r.id != 0x0a)
    return 610;
  p = newstart;
  while (xe_next(fresh, freshlen, &p, &r) && r.id != 0x0a)
    if (xe_cell_id(r.id))
      if ((rc = xe_new_cell(f, &r, strings, styles, ws, keep && action == 1,
                            ws ? ws->workbook->is1904 : 0)))
        return rc;
  return record(f, 0x0a, NULL, 0);
}
/* Keep all original sectors/streams and install a new FAT for the appended
   Workbook. This also preserves embedded objects and document properties. */
static int xe_compound(FILE *out, const char *path, xlsWorkBook *wb, FILE *biff,
                       size_t bytes) {
  FILE *original = fopen(path, "rb");
  if (!original)
    return 601;
  int rc = 610;
  if (fseek(original, 0, SEEK_END)) {
    fclose(original);
    return 603;
  }
  long length = ftell(original);
  if (length < 512 || length % 512 || fseek(original, 0, SEEK_SET)) {
    fclose(original);
    return 610;
  }
  unsigned char *data = malloc((size_t)length);
  if (!data) {
    fclose(original);
    return 909;
  }
  if (fread(data, 1, (size_t)length, original) != (size_t)length) {
    rc = 603;
    goto done;
  }
  OLE2 *ole = wb->olestr->ole;
  if (xe16(data + 30) != 9)
    goto done;
  size_t original_sectors = (size_t)length / 512 - 1,
         stream = bytes < 4096 ? 4096 : bytes;
  size_t n = (stream + 511) / 512, nfat = 1, ndifat = 0;
  for (;;) {
    size_t next = (original_sectors + n + nfat + ndifat + 127) / 128,
           dif = next > 109 ? (next - 109 + 126) / 127 : 0;
    if (next == nfat && dif == ndifat)
      break;
    nfat = next;
    ndifat = dif;
  }
  size_t total = original_sectors + n + nfat + ndifat;
  if (total > UINT32_MAX || nfat > UINT32_MAX) {
    rc = 920;
    goto done;
  }
  uint32_t *fat = malloc(nfat * 512);
  if (!fat) {
    rc = 909;
    goto done;
  }
  for (size_t i = 0; i < nfat * 128; i++)
    fat[i] = UINT32_MAX;
  for (size_t i = 0; i < original_sectors && i < ole->SecIDCount; i++)
    fat[i] = xlsIntVal(ole->SecID[i]);
  size_t dir = ole->dirstart, steps = 0;
  int found = 0;
  while (dir < original_sectors && steps++ < original_sectors) {
    for (size_t i = 0; i < 512; i += 128) {
      unsigned char *entry = data + 512 * (1 + dir) + i;
      if (entry[66] != 2)
        continue;
      unsigned chars = xe16(entry + 64) / 2;
      char name[33] = {0};
      if (chars > 33)
        continue;
      for (unsigned j = 0; j + 1 < chars; j++)
        name[j] = entry[2 * j + 1] ? '?' : (char)entry[2 * j];
      if (!strcmp(name, "Workbook") || !strcmp(name, "Book")) {
        put32(entry + 116, (uint32_t)original_sectors);
        put64(entry + 120, stream);
        found = 1;
      }
    }
    if (dir >= ole->SecIDCount)
      break;
    dir = xlsIntVal(ole->SecID[dir]);
  }
  if (!found) {
    free(fat);
    goto done;
  }
  for (size_t i = 0; i < n; i++)
    fat[original_sectors + i] =
        i + 1 < n ? (uint32_t)(original_sectors + i + 1) : 0xfffffffe;
  for (size_t i = 0; i < nfat; i++)
    fat[original_sectors + n + i] = 0xfffffffd;
  for (size_t i = 0; i < ndifat; i++)
    fat[original_sectors + n + nfat + i] = 0xfffffffc;
  put32(data + 44, (uint32_t)nfat);
  put32(data + 68,
        ndifat ? (uint32_t)(original_sectors + n + nfat) : 0xfffffffe);
  put32(data + 72, (uint32_t)ndifat);
  for (size_t i = 0; i < 109; i++)
    put32(data + 76 + 4 * i,
          i < nfat ? (uint32_t)(original_sectors + n + i) : UINT32_MAX);
  rc = 603;
  if (fwrite(data, 1, (size_t)length, out) != (size_t)length ||
      fseek(biff, 0, SEEK_SET)) {
    free(fat);
    goto done;
  }
  unsigned char sector[512];
  for (size_t i = 0; i < n; i++) {
    memset(sector, 0, 512);
    size_t take = bytes > 512 ? 512 : bytes;
    if (fread(sector, 1, take, biff) != take ||
        fwrite(sector, 1, 512, out) != 512) {
      free(fat);
      goto done;
    }
    bytes -= take;
  }
  for (size_t i = 0; i < nfat; i++) {
    for (size_t j = 0; j < 128; j++)
      put32(sector + 4 * j, fat[i * 128 + j]);
    if (fwrite(sector, 1, 512, out) != 512) {
      free(fat);
      goto done;
    }
  }
  free(fat);
  for (size_t i = 0; i < ndifat; i++) {
    for (size_t j = 0; j < 127; j++) {
      size_t index = 109 + 127 * i + j;
      put32(sector + 4 * j, index < nfat
                                ? (uint32_t)(original_sectors + n + index)
                                : UINT32_MAX);
    }
    put32(sector + 508, i + 1 < ndifat
                            ? (uint32_t)(original_sectors + n + nfat + i + 1)
                            : 0xfffffffe);
    if (fwrite(sector, 1, 512, out) != 512)
      goto done;
  }
  rc = 0;
done:
  free(data);
  fclose(original);
  return rc;
}
static int xe_edit(const char *path, const char *sheet, int action, int keep,
                   int row0, int col0, size_t n, int k, FILE *generated,
                   FILE **result, xlsWorkBook **book) {
  xls_error_t error;
  int rc = 610, which = -1;
  xlsWorkBook *wb = xls_open_file(path, "UTF-8", &error);
  if (!wb)
    return 610;
  *book = wb;
  for (unsigned i = 0; i < wb->sheets.count; i++)
    if (!*sheet ? i == 0 : !strcasecmp(sheet, wb->sheets.sheet[i].name))
      which = (int)i;
  if (which >= 0 && !action)
    return 602;
  if (wb->is5ver)
    return 610;
  size_t oldlen = wb->olestr->size, freshlen;
  unsigned char *old = malloc(oldlen), *fresh = xe_file(generated, &freshlen);
  long *bounds = NULL;
  unsigned nbounds, strings, styles;
  size_t newstart;
  xlsWorkSheet *ws = NULL;
  if (!old || !fresh) {
    rc = 909;
    goto done;
  }
  if (ole2_seek(wb->olestr, 0) < 0 ||
      ole2_read(old, 1, oldlen, wb->olestr) != (ssize_t)oldlen)
    goto done;
  *result = tmpfile();
  if (!*result) {
    rc = 603;
    goto done;
  }
  rc = xe_globals(*result, old, oldlen, fresh, freshlen, which < 0, &bounds,
                  &nbounds, &strings, &styles, &newstart);
  if (rc)
    goto done;
  if (which >= 0 && action == 1) {
    ws = xls_getWorkSheet(wb, which);
    if (!ws || xls_parseWorkSheet(ws) != LIBXLS_OK) {
      rc = 610;
      goto done;
    }
  }
  for (unsigned i = 0; i < nbounds; i++) {
    long offset = ftell(*result);
    if (offset < 0 || offset > UINT32_MAX ||
        fseek(*result, bounds[i], SEEK_SET)) {
      rc = 603;
      goto done;
    }
    unsigned char b[4];
    put32(b, (uint32_t)offset);
    if (fwrite(b, 1, 4, *result) != 4 || fseek(*result, offset, SEEK_SET)) {
      rc = 603;
      goto done;
    }
    if (i == wb->sheets.count || ((int)i == which && action == 2)) {
      size_t p = newstart;
      xe_record r;
      while (xe_next(fresh, freshlen, &p, &r)) {
        rc = xe_cell_id(r.id) ? xe_new_cell(*result, &r, strings, styles, NULL,
                                            0, wb->is1904)
                              : xe_copy(*result, &r);
        if (rc || r.id == 0x0a)
          break;
      }
    } else if ((int)i == which) {
      rc = xe_sheet(*result, old, oldlen, wb->sheets.sheet[i].filepos, fresh,
                    freshlen, newstart, strings, styles, action, keep, row0,
                    col0, n, k, ws);
    } else {
      size_t p = wb->sheets.sheet[i].filepos;
      xe_record r;
      int depth = 0;
      while (xe_next(old, oldlen, &p, &r)) {
        if (r.id == 0x809)
          depth++;
        if ((rc = xe_copy(*result, &r)))
          break;
        if (r.id == 0x0a && --depth == 0)
          break;
      }
    }
    if (rc)
      goto done;
  }
  rc = 0;
done:
  if (ws)
    xls_close_WS(ws);
  free(old);
  free(fresh);
  free(bounds);
  return rc;
}

/* View one XPORT5 member as a standalone library. XPORT is made of 80-byte
   records; adjacent MEMBER/DSCRPTR headers identify member boundaries. */
typedef struct {
  FILE *fp;
  readstat_off_t begin, end, pos;
} cio_xport_view;
static int xport_view_open(const char *path, void *ctx) {
  cio_xport_view *v = ctx;
  v->pos = 0;
  v->fp = fopen(path, "rb");
  return v->fp ? 0 : -1;
}
static int xport_view_close(void *ctx) {
  cio_xport_view *v = ctx;
  if (!v->fp)
    return 0;
  int rc = fclose(v->fp);
  v->fp = NULL;
  return rc;
}
static readstat_off_t xport_view_seek(readstat_off_t offset,
                                      readstat_io_flags_t whence, void *ctx) {
  cio_xport_view *v = ctx;
  readstat_off_t base = whence == READSTAT_SEEK_SET   ? 0
                        : whence == READSTAT_SEEK_CUR ? v->pos
                                                      : 240 + v->end - v->begin;
  if (offset < -base || offset > INT64_MAX - base)
    return -1;
  return v->pos = base + offset;
}
static ssize_t xport_view_read(void *buffer, size_t n, void *ctx) {
  cio_xport_view *v = ctx;
  size_t got = 0;
  readstat_off_t size = 240 + v->end - v->begin;
  if (v->pos >= size)
    return 0;
  if ((uint64_t)n > (uint64_t)(size - v->pos))
    n = (size_t)(size - v->pos);
  while (got < n) {
    size_t take = n - got;
    readstat_off_t offset = v->pos;
    if (offset < 240) {
      if (take > (size_t)(240 - offset))
        take = (size_t)(240 - offset);
    } else
      offset += v->begin - 240;
#ifdef _WIN32
    if (_fseeki64(v->fp, offset, SEEK_SET))
      return -1;
#else
    if (fseeko(v->fp, offset, SEEK_SET))
      return -1;
#endif
    size_t read = fread((char *)buffer + got, 1, take, v->fp);
    v->pos += (readstat_off_t)read;
    got += read;
    if (read < take) {
      if (ferror(v->fp))
        return -1;
      break;
    }
  }
  return (ssize_t)got;
}
static readstat_error_t xport_view_update(long size,
                                          readstat_progress_handler progress,
                                          void *user, void *ctx) {
  cio_xport_view *v = ctx;
  if (progress && progress(size ? (double)v->pos / size : 1, user) ==
                      READSTAT_HANDLER_ABORT)
    return READSTAT_ERROR_USER_ABORT;
  return READSTAT_OK;
}
static int xport_select(const char *path, const char *member,
                        readstat_parser_t *p, cio_xport_view *v) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return 601;
  unsigned char record[240];
  readstat_off_t position = 240, selected = -1, end = 0;
  size_t count = 0;
  char members[32768] = "";
  if (fread(record, 1, 240, fp) != 240) {
    fclose(fp);
    return 459;
  }
  while (fread(record, 1, 80, fp) == 80) {
    if (!memcmp(record, "HEADER RECORD*******MEMBER  ", 28)) {
      if (fread(record + 80, 1, 160, fp) != 160) {
        fclose(fp);
        return 459;
      }
      if (memcmp(record + 80, "HEADER RECORD*******DSCRPTR ", 28)) {
        fclose(fp);
        return 610;
      }
      char name[9];
      memcpy(name, record + 168, 8);
      name[8] = 0;
      for (int j = 7; j >= 0 && name[j] == ' '; j--)
        name[j] = 0;
      if (count && selected >= 0 && !end)
        end = position;
      if ((!*member && !count) || !strcasecmp(member, name))
        selected = position;
      size_t length = strlen(members), extra = strlen(name) + (count ? 1 : 0);
      if (length + extra >= sizeof(members)) {
        fclose(fp);
        return 103;
      }
      if (count)
        strcat(members, " ");
      strcat(members, name);
      count++;
      position += 160;
    }
    position += 80;
    if (SF_poll()) {
      fclose(fp);
      return 1;
    }
  }
  if (ferror(fp)) {
    fclose(fp);
    return 610;
  }
  fclose(fp);
  if (!count)
    return 610;
  if (selected < 0)
    return 111;
  /* A selected member is bounded by the next header, or by EOF. */
  v->begin = selected;
  v->end = end ? end : position;
  v->pos = 0;
  v->fp = NULL;
  int rc = save_number("__cio_nmembers", count);
  if (!rc)
    rc = cio_save("__cio_members", members);
  if (rc)
    return rc;
  readstat_set_open_handler(p, xport_view_open);
  readstat_set_close_handler(p, xport_view_close);
  readstat_set_seek_handler(p, xport_view_seek);
  readstat_set_read_handler(p, xport_view_read);
  readstat_set_update_handler(p, xport_view_update);
  readstat_set_io_ctx(p, v);
  return 0;
}

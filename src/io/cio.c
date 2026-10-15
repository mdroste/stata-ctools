#include "cio_internal.h"

static cio_dataset cache;
void cio_cleanup(void) { cio_dataset_free(&cache); }
ST_retcode cio_main(const char *args) {
  char op[24], kind[24];
  int index = 0;
  if (sscanf(args, "%23s %23s %d", op, kind, &index) < 1)
    return 198;
  if (!strcmp(op, "workbook"))
    return cio_excel_workbook(kind);
  if (!strcmp(op, "blob")) {
    int col, row;
    size_t offset;
    if (sscanf(args, "%*s %d %d %zu", &col, &row, &offset) != 3)
      return 198;
    return cio_blob_info(&cache, col, row, offset);
  }
  if (!strcmp(op, "clear")) {
    cio_dataset_free(&cache);
    return 0;
  }
  if (!strcmp(op, "scan"))
    return cio_scan(&cache, kind);
  if (!strcmp(op, "column"))
    return cio_column_info(&cache, atoi(kind));
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
    return cio_label_info(&cache, atoi(kind));
  if (!strcmp(op, "load"))
    return cio_load(&cache);
  if (!strcmp(op, "export"))
    return cio_do_export(kind);
  return 198;
}

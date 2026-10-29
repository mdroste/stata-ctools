#ifndef CTOOLS_CIO_INTERNAL_H
#define CTOOLS_CIO_INTERNAL_H
#include "cio.h"
#include "formats/ctools_output.h"
#include "formats/excel_address.h"
#include "ctools_config.h"
#include "ctools_data_io.h"
#include "ctools_runtime.h"
#include "ctools_types.h"
#include "vendor/readstat/readstat.h"
#include "vendor/libxls/include/xls.h"
#include "cimport/cimport_xlsx.h"
#include "cio_excel_format.h"
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/* The dispatcher owns a dataset across scan/metadata/load. Readers receive it
 * explicitly; exporters own their output transactions. No reader owns cache. */
typedef struct
{
    char name[301], label[1025], source_format[257], format[80], labelset[301], excel_column[16];
    int is_string, width, date_kind, binary, force_double, format_numbers;
    size_t *lengths;
    double *numbers;
    char **strings;
} cio_column;
typedef struct
{
    char *set, *text;
    double value;
} cio_label;
typedef struct
{
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
typedef struct
{
    FILE *fp;
    readstat_off_t begin, end, pos;
} cio_xport_view;
/* Observations an export writes, resolved once for all of its passes. The
 * selection token is "touse" (a sample marker follows the exported variables),
 * "skipif" (every observation in range) or empty (the plugin call's if). */
typedef struct
{
    ctools_filtered_data marked; /* "touse": the marker's observations */
    perm_idx_t *checked;         /* plain if: observations SF_ifobs selects */
    const perm_idx_t *map;       /* NULL: the range from first */
    size_t n;
    int first, marker;
} cio_rows;
static inline int cio_row(const cio_rows *rows, size_t r)
{
    return rows->map ? (int)rows->map[r] : rows->first + (int)r;
}
static inline int cio_export_nvar(const char *selection)
{
    return SF_nvars() - !strcmp(selection, "touse");
}
int cio_export_rows(cio_rows *rows, const char *selection);
void cio_rows_free(cio_rows *rows);
int cio_local_number(const char *key, size_t v);

uint32_t cio_shp_be32(const unsigned char *p);
void cio_shp_putbe32(unsigned char *p, uint32_t v);
double cio_shp_double(const unsigned char *p);
void cio_shp_putdouble(unsigned char *p, double x);
int cio_shp_dims(unsigned type, int *z, int *m, int *point, int *multi);
int cio_shp_add_row(cio_dataset *dataset, int id, double x, double y, double z, double m, int order,
                    const unsigned char *blob, size_t bytes);
int cio_scan_shp(cio_dataset *dataset, const char *path);
int cio_blob_info(cio_dataset *dataset, int col, int obs, size_t offset);
int cio_shp_export(void);
/* libxls leaves a sheet name NULL when it cannot decode it; read it as "". */
static inline const char *cio_xls_sheet_name(const xlsWorkBook *w, unsigned j)
{
    const char *name = w->sheets.sheet[j].name;
    return name ? name : "";
}
int cio_xls_nonempty(const xlsCell *c);
void cio_xls_formula_record(XLSWORD id, XLSWORD len, XLSBYTE *record);
int cio_xls_result_order(const void *a, const void *b);
void cio_xls_results_reset(void);
int cio_xls_class(const xlsCell *c, const char **text);
int cio_xls_ref(const char *s, int *row, int *col, const char **end);
void cio_xls_column_name(int index, char *name);
uint8_t cio_xls_native_format(xlsWorkBook *w, const xlsCell *cell);
int cio_scan_xls(cio_dataset *dataset, const char *path);
int cio_xport_view_open(const char *path, void *ctx);
int cio_xport_view_close(void *ctx);
readstat_off_t cio_xport_view_seek(readstat_off_t offset, readstat_io_flags_t whence, void *ctx);
ssize_t cio_xport_view_read(void *buffer, size_t n, void *ctx);
readstat_error_t cio_xport_view_update(long size, readstat_progress_handler progress, void *user,
                                       void *ctx);
int cio_xport_select(const char *path, const char *member, readstat_parser_t *p, cio_xport_view *v);
int cio_xlsx_stored(XLSXContext *x, int row, int col);
int cio_xlsx_raw(XLSXContext *x, int row, int col, double *number, const char **text);
uint8_t cio_xlsx_native_format(XLSXContext *x, int row, int col);
int cio_xlsx_numeric(int type);
int cio_scan_xlsx(cio_dataset *dataset, const char *path);
int cio_export_xpf(int nvar, const char *path, int replace, char (*mapped)[301], ctools_output *pending);
int cio_export_sascode(int nvar, char (*mapped)[301], const char *path, const char *kind,
                       int replace, ctools_output *pending);
int cio_excel_sheet_info(int index, const char *name, int row0, int col0, int row1, int col1);
int cio_excel_workbook(const char *kind);
void *cio_realloc(void *p, size_t n, size_t width);
void cio_dataset_free(cio_dataset *d);
void cio_copy_text(char *out, size_t cap, const char *s);
int cio_local_text(const char *key, char *out, int size);
int cio_save(const char *key, const char *value);
int cio_save_number(const char *key, size_t v);
double cio_tagged_missing(char tag);
double cio_number_value(readstat_value_t v);
int cio_reserve_rows(cio_dataset *d, size_t rows);
int cio_civil_days(int y, unsigned m, unsigned d);
void cio_days_civil(int days, int *y, unsigned *m, unsigned *d);
void cio_excel_date_text(double days, uint8_t code, char *out, size_t cap);
int cio_export_meta(int index, const char *field, char *out, int cap);
int cio_metadata_cb(readstat_metadata_t *m, void *ctx);
void cio_import_format(cio_column *c, const char *kind);
int cio_variable_cb(int index, readstat_variable_t *v, const char *set, void *ctx);
int cio_value_cb(int row, readstat_variable_t *v, readstat_value_t value, void *ctx);
int cio_label_cb(const char *set, readstat_value_t value, const char *text, void *ctx);
void cio_error_cb(const char *msg, void *ctx);
int cio_progress_cb(double progress, void *ctx);
int cio_read_xpf(cio_dataset *dataset, const char *path);
int cio_scan(cio_dataset *dataset, const char *kind);
ssize_t cio_write_cb(const void *data, size_t len, void *ctx);
int cio_xport_names(int nvar, char (*names)[301], char (*sets)[301], int rename);
int cio_do_export(const char *kind, const char *selection);
unsigned cio_dbf_u16(const unsigned char *p);
uint32_t cio_dbf_u32(const unsigned char *p);
void cio_dbf_p16(unsigned char *p, unsigned v);
void cio_dbf_p32(unsigned char *p, uint32_t v);
int cio_scan_dbf(cio_dataset *dataset, const char *path);
void cio_dbf_number(double x, int width, char *out, size_t cap);
int cio_export_dbf(const char *selection);
const char *cio_storage(cio_dataset *dataset, cio_column *c, char *out, size_t len);
int cio_column_info(cio_dataset *dataset, int index);
int cio_label_info(cio_dataset *dataset, int index);
int cio_load(cio_dataset *dataset);
#endif

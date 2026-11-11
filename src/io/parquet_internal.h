/* Native Parquet reader for cimport parquet.
 *
 * The reader follows the Apache Parquet format specification: a Thrift
 * compact-protocol footer (FileMetaData), row groups of column chunks, and
 * pages that are compressed (Snappy, Gzip, Zstd, LZ4, Brotli) and encoded (PLAIN,
 * dictionary, RLE/bit-packed hybrid, DELTA_*, BYTE_STREAM_SPLIT). Only flat
 * schemas are decoded: nested (group or repeated) top-level fields are
 * reported with their type name and skipped, as Stata's import parquet does.
 *
 * Types are named and converted the way StataNow's import parquet (which is
 * built on Apache Arrow) presents them, so that both commands create the same
 * variables. */
#ifndef CTOOLS_PARQUET_INTERNAL_H
#define CTOOLS_PARQUET_INTERNAL_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

/* Physical types */
enum
{
    PQ_BOOLEAN = 0,
    PQ_INT32 = 1,
    PQ_INT64 = 2,
    PQ_INT96 = 3,
    PQ_FLOAT = 4,
    PQ_DOUBLE = 5,
    PQ_BYTE_ARRAY = 6,
    PQ_FLBA = 7
};

/* Codecs */
enum
{
    PQ_UNCOMPRESSED = 0,
    PQ_SNAPPY = 1,
    PQ_GZIP = 2,
    PQ_LZO = 3,
    PQ_BROTLI = 4,
    PQ_LZ4 = 5,
    PQ_ZSTD = 6,
    PQ_LZ4_RAW = 7
};

/* Encodings */
enum
{
    PQ_PLAIN = 0,
    PQ_PLAIN_DICTIONARY = 2,
    PQ_RLE = 3,
    PQ_BIT_PACKED = 4,
    PQ_DELTA_BINARY_PACKED = 5,
    PQ_DELTA_LENGTH_BYTE_ARRAY = 6,
    PQ_DELTA_BYTE_ARRAY = 7,
    PQ_RLE_DICTIONARY = 8,
    PQ_BYTE_STREAM_SPLIT = 9
};

/* The value types a column is presented as (Arrow's names, see pq_type_name) */
typedef enum
{
    PA_BOOL,
    PA_INT8,
    PA_INT16,
    PA_INT32,
    PA_INT64,
    PA_UINT8,
    PA_UINT16,
    PA_UINT32,
    PA_UINT64,
    PA_HALF,
    PA_FLOAT,
    PA_DOUBLE,
    PA_DATE32,
    PA_TIME32,
    PA_TIME64,
    PA_TIMESTAMP,
    PA_DECIMAL128,
    PA_DECIMAL256,
    PA_UTF8,
    PA_BINARY,
    PA_FLBA,
    PA_NULL,
    PA_STRUCT,
    PA_LIST,
    PA_MAP,
    PA_INVALID,
    /* Restored from the ARROW:schema metadata of files Arrow writes */
    PA_DICT,
    PA_LARGE_UTF8,
    PA_LARGE_BINARY,
    PA_DURATION
} pq_arrow;

typedef struct
{
    char *name;
    int type; /* physical type, -1 for a group */
    int type_length;
    int repetition; /* 0 required, 1 optional, 2 repeated, -1 unset */
    int num_children;
    int converted; /* -1 unset */
    int scale, precision;
    int logical; /* LogicalType union member id, 0 none */
    int lt_scale, lt_precision, lt_unit, lt_bits, lt_signed;
} pq_schema;

typedef struct
{
    int type, codec, has_path, has_meta;
    int64_t num_values, total_compressed, data_offset, dict_offset;
} pq_chunk;

typedef struct
{
    int64_t num_rows;
    pq_chunk *chunks; /* one per leaf */
} pq_group;

typedef struct
{
    const char *name;
    pq_arrow arrow;
    int leaf;    /* leaf column index, -1 when nested */
    int element; /* schema index */
    int max_def;
    int unit; /* time units: 1 ms, 2 us, 3 ns */
    int scale, precision;
    int physical, type_length;
    int first_leaf, nleaf; /* the leaf columns under this field */
    int origin_type, origin_dict; /* from ARROW:schema; -1 when absent */
} pq_field;

typedef struct
{
    int type; /* 0 data, 1 index, 2 dictionary, 3 data v2 */
    int32_t usize, csize;
    int32_t num_values;
    int encoding, def_enc, rep_enc;
    int32_t num_nulls, num_rows, def_len, rep_len;
    int compressed;
} pq_page;

typedef struct
{
    FILE *fp;
    char *path; /* column readers open their own handles, so they can run in parallel */
    int64_t size;
    uint8_t *footer;
    pq_schema *schema;
    int nschema;
    int64_t num_rows;
    pq_group *groups;
    int ngroups;
    int nleaves;
    pq_field *fields;
    int nfields;
} pq_file;

/* One output column's decoded rows (only rows in the requested range are
 * kept, but strings are measured over every row, as Stata measures them). */
typedef struct
{
    const pq_field *field;
    size_t n;
    uint8_t *valid;
    double *num;
    int64_t *i64;
    uint64_t *off; /* n+1 offsets into bytes */
    uint8_t *bytes;
    size_t nbytes, capbytes;
    size_t max_len; /* longest value in the whole column, in bytes */
    int has_nul;
    int any;      /* a non-missing value in the range */
    double maxv;  /* largest value in the range (NaN ignored) */
    int wide;     /* an int64 value in the range is outside Stata's long */
    /* Dictionary columns: entries in order, unified over row groups */
    uint64_t *doff;
    uint8_t *dbytes;
    size_t dn, dcap, dbn, dbcap, dmax_len;
    uint32_t *dhash;
    size_t dhcap;
} pq_values;

/* Thrift footer (parquet_thrift.c) */
int pq_open(pq_file *f, const char *path, char *err, size_t errcap);
void pq_close(pq_file *f);
const char *pq_type_name(pq_arrow a);
int pq_page_header(const uint8_t *p, size_t n, pq_page *h, size_t *used);

/* Decompression (parquet_codec.c); returns 0 or an error message */
const char *pq_decompress(int codec, const uint8_t *src, size_t n, uint8_t *dst, size_t out);
const char *pq_zstd_decompress(const uint8_t *src, size_t n, uint8_t *dst, size_t out);
const char *pq_brotli_decompress(const uint8_t *src, size_t n, uint8_t *dst, size_t out);

/* Column decoding (parquet_decode.c). first/last are 0-based inclusive row
 * bounds; measure_all 1 decodes every row group where the result depends on
 * the whole column (text widths, dictionaries), 2 decodes everything. */
int pq_read_column(pq_file *f, const pq_field *field, int64_t first, int64_t last, int measure_all,
                   pq_values *out, char *err, size_t errcap);
void pq_values_free(pq_values *v);
int pq_check_field(pq_file *f, const pq_field *field, char *err, size_t errcap);
void pq_decimal_text(const uint8_t *be, size_t n, int scale, char *out, size_t cap);

#endif

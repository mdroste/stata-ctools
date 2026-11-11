#ifndef CTOOLS_DTA_H
#define CTOOLS_DTA_H
/*
 * Stata .dta writer: format 118 (up to 32,767 variables) or 119.
 *
 * The plugin writes a dataset to a file that the ado file loads with
 * `use, clear`. This reaches what the SPI cannot write: strL values of any
 * length (binary included), every storage type, new variables, variable and
 * value labels, formats, characteristics and the data label. See
 * docs/development/DTA_WRITER.md for the format facts and measurements.
 *
 * Indices: columns are 0-based positions in spec->cols, and the file numbers
 * variables from 1 in the same order. Rows are 0-based; the file numbers
 * observations from 1. The writer never calls the SPI, so it may run with
 * data the caller has already copied out of Stata.
 *
 * Layout (what drives the design): records are row-major with fields packed
 * in variable order, and Stata reads the strL section in step with the strL
 * cells, observation by observation. Work is therefore split by blocks of
 * rows, built in parallel and appended to the file in order by one writer
 * thread, the fastest way to write a file on NTFS.
 */
#include <stddef.h>
#include <stdint.h>

typedef enum
{
    CTOOLS_DTA_BYTE = 1,
    CTOOLS_DTA_INT,
    CTOOLS_DTA_LONG,
    CTOOLS_DTA_FLOAT,
    CTOOLS_DTA_DOUBLE,
    CTOOLS_DTA_STR, /* str#: width 1..2045 */
    CTOOLS_DTA_STRL
} ctools_dta_type;

typedef enum
{
    /* const double *: Stata's double encoding, missing values included. For
     * the numeric types; an integer type takes only integers in its range. */
    CTOOLS_DTA_SRC_F64 = 1,
    /* const char *const *: NUL-terminated text; a NULL pointer reads as "". */
    CTOOLS_DTA_SRC_CSTR,
    /* const char *: slots of `stride` bytes; a value ends at its first NUL
     * or at the end of its slot. */
    CTOOLS_DTA_SRC_FLAT,
    /* const ctools_dta_span *: any bytes. */
    CTOOLS_DTA_SRC_SPAN,
    /* Filled a row group at a time by spec->producer. */
    CTOOLS_DTA_SRC_PRODUCER
} ctools_dta_source;

/* A string value. binary: 0 stores text unless the bytes hold a NUL (then
 * binary); 1 stores a binary strL. A str# value must not hold NUL bytes. */
typedef struct
{
    const char *ptr;
    size_t len;
    int binary;
} ctools_dta_span;

/* Column flag: the column's values never hold NUL bytes, so strL values are
 * stored as text without scanning them. */
#define CTOOLS_DTA_COL_TEXT 1u
/* Column flag for a strL column read through rows: cells that read the same
 * source row share one stored value (the first such cell stores it, the
 * others point to it, as Stata's own files do), so the file holds and use
 * loads each value once. */
#define CTOOLS_DTA_COL_SHARE 2u

typedef struct
{
    const char *name; /* UTF-8 Stata name */
    ctools_dta_type type;
    int width;               /* str#: 1..2045, or 0 for the longest value (not PRODUCER) */
    const char *format;      /* NULL: the type's default */
    const char *label;       /* variable label, NULL for none */
    const char *value_label; /* attached value label, NULL for none */
    ctools_dta_source source;
    const void *data;
    size_t stride; /* FLAT: bytes per slot */
    /* Optional, not PRODUCER: output row i reads source row rows[i]; a
     * negative entry gives missing (.) or "". */
    const int64_t *rows;
    unsigned flags;
} ctools_dta_column;

/* A characteristic of the dataset (owner "_dta") or of a column. Notes are
 * characteristics: owner[note0] holds their count, owner[note1]... the text. */
typedef struct
{
    const char *owner;
    const char *name;
    const char *text;
} ctools_dta_char;

/* A value label: n distinct values and their texts, in any order. */
typedef struct
{
    const char *name;
    size_t n;
    const int32_t *values;
    const char *const *texts;
} ctools_dta_value_label;

/* A filled row group: one entry per PRODUCER column, in column order, each
 * of the group's row count: double * for a numeric column, ctools_dta_span *
 * for a string column. The writer allocates the array; fill sets it. */
typedef struct
{
    void **columns;
    void *state; /* the producer's, handed back to release */
} ctools_dta_group;

enum
{
    CTOOLS_DTA_PHASE_CELLS = 0, /* every PRODUCER column is read */
    CTOOLS_DTA_PHASE_STRLS = 1  /* only strL PRODUCER columns are read */
};

/* Rows of PRODUCER columns, a group at a time. fill runs on worker threads,
 * concurrently for different groups, and returns a Stata return code. It
 * may be called once or twice per group (the data come before every strL
 * in the file, so a group is filled again for its strLs unless the writer
 * could keep it); release is called once for every successful fill. */
typedef struct
{
    size_t ngroups;
    const uint64_t *starts; /* ngroups + 1 offsets: 0, ascending, ..., nobs */
    const uint64_t *bytes;  /* optional: memory a filled group holds */
    void *ctx;
    int (*fill)(void *ctx, size_t group, int phase, ctools_dta_group *out);
    void (*release)(void *ctx, ctools_dta_group *group);
} ctools_dta_producer;

typedef struct
{
    uint64_t nobs;
    size_t ncols;
    const ctools_dta_column *cols;
    int release;            /* 0: 118, or 119 past 32,767 variables; 118 or 119 */
    const char *data_label; /* NULL for none */
    const char *timestamp;  /* NULL for none, else "dd Mon yyyy hh:mm" */
    /* 1-based column numbers the rows are sorted by. Stata trusts this
     * without checking it, so set it only for an order the caller made. */
    const int *sortlist;
    size_t nsort;
    const ctools_dta_char *chars;
    size_t nchars;
    const ctools_dta_value_label *labels;
    size_t nlabels;
    const ctools_dta_producer *producer;
    int threads;       /* 0: ctools_get_max_threads(); 1: no worker threads */
    int temporary;     /* nonzero: a short-lived file (FILE_ATTRIBUTE_TEMPORARY) */
    size_t keep_bytes; /* producer groups kept between the passes; 0: 1 GiB */
} ctools_dta_spec;

typedef struct
{
    uint64_t file_bytes;
    uint64_t strls;       /* strL values stored */
    double seconds;       /* the whole write */
    double write_seconds; /* of which in file writes */
    char message[256];    /* why a write failed */
} ctools_dta_result;

void ctools_dta_spec_init(ctools_dta_spec *spec);

/* Write spec to path, replacing any file there. Returns 0, or 198 (invalid
 * spec), 459 (a value that does not fit its storage type), 603 (cannot open
 * the file), 693 (write failed), 920 (memory) or a producer's code. A failed
 * write removes the file. result may be NULL. */
int ctools_dta_write(const char *path, const ctools_dta_spec *spec, ctools_dta_result *result);

/* Internal plugin command "dta" (tests and tools; ctools_dta_plugin.c). */
int ctools_dta_main(const char *args);

#endif /* CTOOLS_DTA_H */

#ifndef CTOOLS_DTA_SCHEMA_H
#define CTOOLS_DTA_SCHEMA_H
/*
 * The description of a Stata dataset, without its observations.
 *
 * A command that rebuilds the dataset must carry over what the SPI cannot
 * reach: storage types, display formats, variable and value labels,
 * characteristics (notes included), the data label and unattached value
 * labels. The ado file saves the dataset with no observations (a file of
 * about 640 bytes per variable), the plugin reads that description here,
 * edits it and writes the description of the result, and the ado file loads
 * it with `use` before the plugin fills in the values through the SPI.
 *
 * Formats 118 and 119, little-endian, as Stata 14 and later write them (see
 * docs/development/DTA_WRITER.md for the layout). Records are kept as the
 * file holds them: fixed-width, NUL-padded fields are copied byte for byte,
 * so a format or label this code does not understand still survives. Value
 * labels stay in their file encoding, one record per label.
 *
 * Indices are zero-based positions in vars[]; the sort list holds zero-based
 * positions too (the file numbers variables from 1).
 */
#include <stddef.h>
#include <stdint.h>

enum
{
    CTOOLS_DTA_CODE_STR_MAX = 2045, /* 1..2045: str# */
    CTOOLS_DTA_CODE_STRL = 32768,
    CTOOLS_DTA_CODE_DOUBLE = 65526,
    CTOOLS_DTA_CODE_FLOAT = 65527,
    CTOOLS_DTA_CODE_LONG = 65528,
    CTOOLS_DTA_CODE_INT = 65529,
    CTOOLS_DTA_CODE_BYTE = 65530
};

#define CTOOLS_DTA_NAME_BYTES 129  /* a name: up to 32 characters, 128 bytes, and a NUL */
#define CTOOLS_DTA_FORMAT_BYTES 57 /* a display format and a NUL */
#define CTOOLS_DTA_LABEL_BYTES 321 /* a variable or data label: 80 characters and a NUL */

typedef struct
{
    uint16_t type; /* a CTOOLS_DTA_CODE_* value or a str# width */
    char name[CTOOLS_DTA_NAME_BYTES];
    char format[CTOOLS_DTA_FORMAT_BYTES];
    char value_label[CTOOLS_DTA_NAME_BYTES];
    char label[CTOOLS_DTA_LABEL_BYTES];
} ctools_dta_schema_var;

/* A characteristic of the dataset (owner "_dta") or of a variable. */
typedef struct
{
    char owner[CTOOLS_DTA_NAME_BYTES];
    char name[CTOOLS_DTA_NAME_BYTES];
    char *text; /* owned, NUL-terminated */
} ctools_dta_schema_char;

/* One value label: its name and its table as the file stores it
 * (u32 count, u32 text bytes, u32 offsets[count], i32 values[count], text). */
typedef struct
{
    char name[CTOOLS_DTA_NAME_BYTES];
    unsigned char *table; /* owned */
    size_t table_bytes;
} ctools_dta_schema_label;

typedef struct
{
    int release; /* 118 or 119 */
    size_t nvars;
    ctools_dta_schema_var *vars; /* owned [nvars] */
    size_t nsort;
    uint32_t *sort; /* owned [nsort]: zero-based positions in vars */
    char data_label[CTOOLS_DTA_LABEL_BYTES];
    size_t nchars, chars_cap;
    ctools_dta_schema_char *chars; /* owned */
    size_t nlabels, labels_cap;
    ctools_dta_schema_label *labels; /* owned */
} ctools_dta_schema;

void ctools_dta_schema_init(ctools_dta_schema *s);
void ctools_dta_schema_free(ctools_dta_schema *s); /* idempotent */

/*
 * Read the description of the dataset in path. The file may hold
 * observations; they are skipped. Returns 0, or a Stata return code with the
 * reason in msg: 601 (cannot open), 610 (not a format 118 or 119 file, or a
 * damaged one), 920 (memory). On failure s is left freed.
 */
int ctools_dta_schema_read(const char *path, ctools_dta_schema *s, char *msg, size_t cap);

/*
 * Write s to path as a dataset with no observations, replacing any file
 * there. The release follows the number of variables (119 past 32,767).
 * Returns 0, 198 (a description Stata could not load: too many variables, a
 * sort entry out of range), 603 (cannot create), 693 (write failed) or 920.
 * A failed write removes the file.
 */
int ctools_dta_schema_write(const char *path, const ctools_dta_schema *s, char *msg, size_t cap);

/* Characteristics. owner is "_dta" or a variable name. */
const char *ctools_dta_schema_char_get(const ctools_dta_schema *s, const char *owner,
                                       const char *name);
/* Set (add or replace) a characteristic; an empty text removes it, as
 * `char define` does. Returns 0 or 920. */
int ctools_dta_schema_char_set(ctools_dta_schema *s, const char *owner, const char *name,
                               const char *text);
void ctools_dta_schema_char_remove(ctools_dta_schema *s, const char *owner, const char *name);
/* Append a characteristic without looking for an earlier one of the same
 * owner and name (the caller knows there is none). Returns 0 or 920. */
int ctools_dta_schema_char_add(ctools_dta_schema *s, const char *owner, const char *name,
                               const char *text);

/* Value labels. find returns the position in labels[] or -1. */
long ctools_dta_schema_label_find(const ctools_dta_schema *s, const char *name);
/* Entries of labels[at]: count, and for entry k its value and text (the text
 * points into the table). Returns 0, or 610 when the table is damaged. */
int ctools_dta_schema_label_count(const ctools_dta_schema *s, size_t at, size_t *count);
int ctools_dta_schema_label_entry(const ctools_dta_schema *s, size_t at, size_t k, int32_t *value,
                                  const char **text);
/* Define a value label from n (value, text) pairs with distinct values, in
 * any order. Returns 0, 110 when the name is already defined, 198 or 920. */
int ctools_dta_schema_label_define(ctools_dta_schema *s, const char *name, size_t n,
                                   const int32_t *values, const char *const *texts);

/* Stata's default display format for a storage type (str#: %#s, at least %9s). */
void ctools_dta_schema_default_format(uint16_t type, char format[CTOOLS_DTA_FORMAT_BYTES]);

#endif /* CTOOLS_DTA_SCHEMA_H */

#ifndef CIMPORT_CACHE_H
#define CIMPORT_CACHE_H

#include <stdatomic.h>
#include "cimport_parse.h"

/* Private, process-local cache. Blocks own their descriptors and payload in one
 * allocation. Source references are relative to file_start; offsets stay 64-bit
 * in the parser. A block spanning more than UINT32_MAX bytes is not cached. */
enum
{
    CIMPORT_CACHE_MISSING,
    CIMPORT_CACHE_CONSTANT,
    CIMPORT_CACHE_U8,
    CIMPORT_CACHE_U16,
    CIMPORT_CACHE_U32,
    CIMPORT_CACHE_FLOAT,
    CIMPORT_CACHE_DOUBLE,
    CIMPORT_CACHE_TEXT,
    CIMPORT_CACHE_FIXED_TEXT,
    CIMPORT_CACHE_INLINE_TEXT
};

typedef struct
{
    double base;
    size_t offset;
    uint32_t length;
    unsigned char kind;
} CImportCachedColumn;

typedef struct CImportCachedBlock
{
    struct CImportCachedBlock *next;
    uint64_t file_start, file_end;
    size_t rows, bytes;
    CImportCachedColumn columns[];
} CImportCachedBlock;

typedef struct
{
    uint64_t limit;
    atomic_uint_fast64_t used;
} CImportCacheBudget;

typedef struct
{
    double *values;   /* Column-major scratch, stride capacity */
    uint64_t *fields; /* Relative offset in high half, length/flags in low half */
    CImportCachedColumn *plan;
    unsigned char *numeric;
    size_t capacity, rows;
    uint64_t file_start, file_end;
    int columns;
    bool wide;
} CImportCacheBuilder;

int cimport_cache_init(CImportCacheBuilder *builder, int columns);
void cimport_cache_destroy(CImportCacheBuilder *builder);
void cimport_cache_append(CImportCacheBuilder *builder, const CImportFieldRef *fields, int count,
                          const double *values, uint64_t start, uint64_t end);
/* NULL means caching is unavailable (budget, allocation, or reference width).
 * The caller retains prior blocks and reparses this block and its suffix later. */
CImportCachedBlock *cimport_cache_pack(CImportCacheBuilder *builder, CImportCacheBudget *budget,
                                       double missing, bool as_float, bool as_double,
                                       const char *file_data, char quote);
void cimport_cache_free(CImportCachedBlock *block, CImportCacheBudget *budget);
double cimport_cache_number(const CImportCachedBlock *block, int column, size_t row,
                            double missing);
CImportFieldRef cimport_cache_field(const CImportCachedBlock *block, int column, size_t row);
const char *cimport_cache_text(const CImportCachedBlock *block, int column, size_t row);

#endif

#ifndef CSQL_SCAN_H
#define CSQL_SCAN_H
/*
 * The machinery both table functions share (__csql_stata and read_dta).
 *
 * DuckDB keeps a query's row order (preserve_insertion_order) by running a
 * plain SELECT, filter or join over a source that reports no batch order on
 * one thread, and the C API has no way to report it. When only one DuckDB
 * thread reads a scan, the scan therefore decodes batches of chunks ahead
 * on csql's own worker threads (ctx->pool) and hands DuckDB one finished
 * chunk per call; when several do (aggregates, joins' build sides, ORDER BY),
 * each call decodes its own chunk. Rows come out in order either way.
 *
 * A scan decodes rows into csql_pre buffers in DuckDB's physical layout:
 * values, a validity bitmap and, for VARCHAR, duckdb_string_t entries that
 * point into memory valid for the whole query (a mapped .dta) or into the
 * batch's arena (copied into the vector when served). Columns flagged as
 * deferred (strL in memory) are read when the chunk is served, through
 * csql_read_strls.
 */
#include "csql/csql_internal.h"

typedef struct
{
    char *p;
    size_t len, cap;
} csql_arena;

typedef struct
{
    size_t n;         /* rows */
    int64_t r0;       /* first row (0-based) */
    void **data;      /* [nproj] physical values */
    uint64_t **valid; /* [nproj] bitmaps, 1 = valid */
    unsigned char *nulls; /* [nproj] any row invalid */
    unsigned char *copy;  /* [nproj] strings in the arena: copy when served */
    csql_arena arena;
    int rc;
} csql_pre;

/* Decode rows [pre->r0, pre->r0 + pre->n) of the projected columns into pre
 * (buffers allocated by the engine). Runs on any thread; returns 0 or a
 * Stata return code. */
typedef int (*csql_decode_fn)(void *bind, const idx_t *proj, size_t nproj, csql_pre *pre);

typedef struct
{
    csql_ctx *ctx;
    void *bind;
    int64_t nrows;
    size_t nsrc;          /* source columns */
    const csql_map *map;  /* [nsrc] */
    const int *deferred;  /* [nsrc] strL read when served; NULL for none */
    int stable;           /* string pointers outlive the query (mapped file) */
    csql_decode_fn decode;
} csql_scan_def;

/* Bitmap helpers */
static inline void csql_valid_set(uint64_t *v, size_t i) { v[i >> 6] |= UINT64_C(1) << (i & 63); }

/* A string at slot s: inlined when it fits, else pointing at p. */
void csql_str_set(duckdb_string_t *s, const char *p, uint32_t len);
/* Make room for bytes more in the arena (call before handing out any of
 * it: the arena never moves what it has handed out); 0 or 920. */
int csql_arena_reserve(csql_arena *a, size_t bytes);
/* Copy n bytes into the reserved room; returns their address. */
char *csql_arena_add(csql_arena *a, const char *p, size_t n);

/* Table function callbacks: the bind data must begin with a csql_scan_def. */
void csql_scan_init(duckdb_init_info info);
void csql_scan_local_init(duckdb_init_info info);
void csql_scan_function(duckdb_function_info info, duckdb_data_chunk output);

/* Physical bytes per value of a mapped type. */
size_t csql_map_width(csql_map m);

/* The worker threads a call created for decoding (after the database is
 * closed). */
void csql_scan_pool_free(csql_ctx *ctx);

#endif /* CSQL_SCAN_H */

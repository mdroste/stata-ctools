#ifndef CSQL_INTERNAL_H
#define CSQL_INTERNAL_H
/*
 * csql internals: one plugin call opens an in-memory (or file) DuckDB
 * database, registers the Stata dataset and the .dta reader as table
 * functions, runs the SQL script, sends the last statement's rows to a .dta
 * file (which the ado loads with use), another file, or the screen, and
 * closes the database again. Nothing is kept between plugin calls, so no
 * DuckDB thread runs while Stata does.
 *
 * Threads: DuckDB runs the query on its own worker threads and on the Stata
 * thread that called the plugin. Worker threads read numeric and str#
 * variables through the SPI, as ctools' loaders do; strL values are read
 * only on the Stata thread, which serves requests from the workers while it
 * waits for the query (csql_service_*).
 */
#include <stddef.h>
#include <stdint.h>
#include <pthread.h>
#include <stdatomic.h>

#include "csql/csql_api.h"
#include "formats/ctools_dta_schema.h"

#define DK(fn) (csql_db.duckdb_##fn)

#define CSQL_NAME_BYTES CTOOLS_DTA_NAME_BYTES
#define CSQL_FORMAT_BYTES CTOOLS_DTA_FORMAT_BYTES
#define CSQL_LABEL_BYTES CTOOLS_DTA_LABEL_BYTES
#define CSQL_STR_MAX 2045

/* Stata missing (.) as stored in a double; every value from 2^1023 up is
 * missing (.a to .z follow .), as is any NaN. */
#define CSQL_MISSING_BITS UINT64_C(0x7FE0000000000000)
static inline int csql_is_missing(double z) { return !(z < 8.98846567431157954e307); }
double csql_missing(void);

/* Stata's date origin is 1960-01-01, DuckDB's 1970-01-01. */
#define CSQL_DAYS_1960_1970 3653
#define CSQL_MS_1960_1970 INT64_C(315619200000)

/* ------------------------------------------------------------ metadata */

/* A value label: n (value, text) pairs. */
typedef struct
{
    char name[CSQL_NAME_BYTES];
    size_t n;
    int32_t *values;
    char **texts;
} csql_vlabel;

/* What Stata knows about a variable beyond its values. */
typedef struct
{
    char name[CSQL_NAME_BYTES];
    /* the column name SQL sees: name, or name_2, name_3 ... for a variable
     * whose name differs from an earlier one only in case (SQL names are
     * not case-sensitive) */
    char sqlname[CSQL_NAME_BYTES + 8];
    uint16_t type; /* CTOOLS_DTA_CODE_* or a str# width */
    char format[CSQL_FORMAT_BYTES];
    char vlabel[CSQL_NAME_BYTES];
    char label[CSQL_LABEL_BYTES];
} csql_meta_var;

/* A dataset whose variables a query may pass through: the data in memory or
 * a .dta file. Results take labels and formats from the first source that
 * was read and has a variable of the same name. */
typedef struct csql_source
{
    char origin[4096]; /* "stata" or the path */
    size_t nvars;
    csql_meta_var *vars;
    size_t nlabels;
    csql_vlabel *labels;
    int used; /* bound by a statement of this call */
    int order; /* position among used sources */
    struct csql_source *next;
} csql_source;

void csql_source_free(csql_source *s);
const csql_vlabel *csql_source_label(const csql_source *s, const char *name);

/* How a Stata numeric variable reads in SQL, from its storage type and
 * display format: DATE for %td, TIMESTAMP for %tc, else by storage type. */
typedef enum
{
    CSQL_MAP_TINYINT = 1,
    CSQL_MAP_SMALLINT,
    CSQL_MAP_INTEGER,
    CSQL_MAP_FLOAT,
    CSQL_MAP_DOUBLE,
    CSQL_MAP_DATE,
    CSQL_MAP_TIMESTAMP,
    CSQL_MAP_VARCHAR
} csql_map;

csql_map csql_map_stata_var(uint16_t type, const char *format, int datetime);
duckdb_logical_type csql_map_logical(csql_map m);
/* Store a Stata double z (not missing) into a DuckDB vector slot. */
void csql_map_store(csql_map m, void *data, size_t row, double z, int *invalid);

/* ------------------------------------------------------------ context */

/* A request from a worker for strL values, served by the Stata thread. */
typedef struct
{
    int var;        /* 1-based plugin variable */
    int64_t row0;   /* first 1-based observation */
    size_t n;       /* observations */
    duckdb_vector out;
    int rc;
    int done;
} csql_strl_request;

typedef struct csql_ctx
{
    duckdb_database db;
    duckdb_connection con;
    uint64_t main_thread;
    int threads;
    int datetime;  /* map %td/%tc variables to DATE/TIMESTAMP */
    int verbose;

    /* The data in memory: plugin variables 1..k, observations 1..N. */
    csql_source *stata;
    int64_t stata_nobs;
    int stata_binds; /* times a statement bound __csql_stata */

    /* Sources registered (the data in memory and every .dta opened). */
    csql_source *sources;
    int nused;

    /* strL service */
    pthread_mutex_t lock;
    pthread_cond_t wake;  /* the Stata thread waits here */
    pthread_cond_t served;
    csql_strl_request *pending;
    int waiting_workers;
    atomic_int cancelled;
    int service_ready;

    /* worker threads for decoding batches (csql_scan.c) */
    void *pool;
    int pool_failed;

    char error[2048];
} csql_ctx;

uint64_t csql_thread_id(void);
int csql_on_main(const csql_ctx *ctx);

/* Ask the Stata thread for strL values (from a worker), or read them
 * directly (on the Stata thread). Returns 0 or a Stata return code. */
int csql_read_strls(csql_ctx *ctx, int var, int64_t row0, size_t n, duckdb_vector out);
/* The Stata thread: serve one pending request, if any. */
void csql_service_once(csql_ctx *ctx);

/* ------------------------------------------------------------ modules */

/* csql_meta.c: the ado's description of the data in memory. */
int csql_meta_read(const char *path, csql_source **out, char *msg, size_t cap);
/* "byte", "str12", "strL", ... as a .dta type code; 0 when invalid. */
uint16_t csql_meta_type_code(const char *type, size_t len);

/* csql_stata_scan.c: table function __csql_stata(). */
int csql_stata_register(csql_ctx *ctx);

/* csql_dta_scan.c: table function read_dta(path) and FROM 'file.dta'. */
int csql_dta_register(csql_ctx *ctx);
/* Unmap every .dta opened by this call (after the database is closed). */
void csql_dta_close_all(void);

/* csql_result.c: types stored as they are (others are cast to VARCHAR). */
int csql_type_supported(duckdb_type t);

/* csql_sql.c */
/* Byte range [*start, *end) of the last statement of script (without a
 * trailing semicolon); 0, or -1 when the script has no statement. */
int csql_sql_last_statement(const char *sql, size_t *start, size_t *end);
/* Append text to a growing buffer; 0 or 920. */
typedef struct
{
    char *p;
    size_t len, cap;
} csql_buf;
int csql_buf_add(csql_buf *b, const char *text, size_t n);
int csql_buf_str(csql_buf *b, const char *text);
int csql_buf_ident(csql_buf *b, const char *name);   /* "name" with "" doubled */
int csql_buf_literal(csql_buf *b, const char *text); /* 'text' with '' doubled */
void csql_buf_free(csql_buf *b);
/* A valid, unique Stata name from a column name (strtoname rules); returns
 * 1 when the name changed. taken: names already used, n of them. */
int csql_stata_name(const char *in, char out[CSQL_NAME_BYTES], char (*taken)[CSQL_NAME_BYTES],
                    size_t n);
int csql_name_reserved(const char *name);
/* Equal names, ignoring the case of ASCII letters (as SQL compares names). */
int csql_name_eq_ci(const char *a, const char *b);
/* Set vars[].sqlname: unique without regard to case. */
int csql_sql_names(csql_meta_var *vars, size_t n);

/* csql_result.c */
typedef enum
{
    CSQL_OUT_DISPLAY = 0,
    CSQL_OUT_DTA,  /* write a .dta (loaded by the ado for clear) */
} csql_out_mode;

typedef struct
{
    csql_out_mode mode;
    const char *path;    /* DTA */
    int temporary;       /* the file is the ado's tempfile */
    int labels;          /* carry variable and value labels */
    int64_t rows;        /* DISPLAY: rows to list */
    int linesize;
    int strl_min;        /* strings longer than this become strL (<= 2045) */
} csql_out_opts;

/* Receive the rows of res (a query result) and send them to opts. */
int csql_result_output(csql_ctx *ctx, duckdb_result *res, int streaming, const csql_out_opts *opts,
                       int64_t *nrows, int64_t *ncols);

/* csql_install.c: verify and unpack a DuckDB release archive. */
int csql_install_extract(const char *archive, const char *sha256_hex, const char *member,
                         const char *dest, char *msg, size_t cap);

#endif /* CSQL_INTERNAL_H */

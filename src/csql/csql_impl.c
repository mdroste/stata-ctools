/*
 * csql_impl.c - run an SQL script in DuckDB from Stata.
 *
 * One plugin call does all of it: load the DuckDB library, open a database
 * (in memory unless database() names a file), register the data in memory
 * (__csql_stata, views `stata` and `data`) and the .dta reader, run every
 * statement of the script, send the rows of the last one to their
 * destination, close the database. The Stata thread runs DuckDB's tasks
 * itself between polls for Break and serves strL reads (csql_internal.h).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <sched.h>
#endif

#include "stplugin.h"
#include "csql/csql_impl.h"
#include "csql/csql_internal.h"
#include "csql/csql_scan.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "ctools_cancel.h"
#include "ctools_parse.h"
#include "ctools_file.h"

/* ------------------------------------------------------------ helpers */

/* A local macro of any length (NULL on failure). */
static char *csql_local(const char *name)
{
    char macro[160];
    snprintf(macro, sizeof(macro), "_%s", name);
    for (size_t cap = (size_t)1 << 12; cap <= ((size_t)1 << 24); cap *= 8) {
        char *buf = calloc(cap, 1);
        if (!buf) return NULL;
        if (SF_macro_use(macro, buf, (int)cap) != 0) buf[0] = 0;
        if (strlen(buf) < cap - 1) return buf;
        free(buf);
    }
    return NULL;
}

/* Print text line by line (SF_error and SF_display take one line a call). */
static void csql_print(int error, const char *prefix, const char *text)
{
    const char *p = text;
    int first = 1;
    while (p && *p) {
        const char *nl = strchr(p, '\n');
        size_t n = nl ? (size_t)(nl - p) : strlen(p);
        while (n && (p[n - 1] == '\r' || p[n - 1] == ' ')) n--;
        if (n || !first) {
            size_t cap = n + (prefix ? strlen(prefix) : 0) + 4;
            char *line = malloc(cap);
            if (!line) return;
            snprintf(line, cap, "%s%.*s\n", first && prefix ? prefix : "", (int)n, p);
            /* SF_error shows SMCL as typed; braces in DuckDB's messages are literal */
            if (error) SF_error(line);
            else SF_display(line);
            free(line);
        }
        first = 0;
        p = nl ? nl + 1 : NULL;
    }
}

static int csql_fail(const char *text, int rc)
{
    csql_print(1, "csql: ", text);
    return rc;
}

static char *read_text(const char *path, char *msg, size_t cap)
{
    FILE *f = ctools_fopen(path, "rb");
    if (!f) {
        snprintf(msg, cap, "could not open %s", path);
        return NULL;
    }
    ctools_fseek(f, 0, SEEK_END);
    long long n = (long long)ctools_ftell(f);
    ctools_fseek(f, 0, SEEK_SET);
    char *buf = n >= 0 ? malloc((size_t)n + 1) : NULL;
    if (!buf || fread(buf, 1, (size_t)n, f) != (size_t)n) {
        fclose(f);
        free(buf);
        snprintf(msg, cap, "could not read %s", path);
        return NULL;
    }
    fclose(f);
    buf[n] = '\0';
    /* a UTF-8 byte order mark */
    if (n >= 3 && (unsigned char)buf[0] == 0xEF && (unsigned char)buf[1] == 0xBB &&
        (unsigned char)buf[2] == 0xBF)
        memmove(buf, buf + 3, (size_t)n - 2);
    return buf;
}

static void csql_yield(void)
{
#ifdef _WIN32
    SwitchToThread();
#else
    sched_yield();
#endif
}

/* Wait up to a millisecond for a strL request. */
static void csql_wait(csql_ctx *ctx)
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    ts.tv_nsec += 1000000;
    if (ts.tv_nsec >= 1000000000) {
        ts.tv_sec++;
        ts.tv_nsec -= 1000000000;
    }
    pthread_mutex_lock(&ctx->lock);
    if (!ctx->pending) pthread_cond_timedwait(&ctx->wake, &ctx->lock, &ts);
    pthread_mutex_unlock(&ctx->lock);
}

static void csql_cancel(csql_ctx *ctx)
{
    if (atomic_exchange(&ctx->cancelled, 1)) return;
    DK(interrupt)(ctx->con);
    pthread_mutex_lock(&ctx->lock);
    pthread_cond_broadcast(&ctx->served);
    pthread_mutex_unlock(&ctx->lock);
}

/*
 * Execute a prepared statement: the Stata thread runs DuckDB's tasks one at
 * a time, serving strL requests and polling Break between them. Returns 0
 * with *res set, 1 on Break, or 198 with the error in ctx->error.
 */
static int csql_execute(csql_ctx *ctx, duckdb_prepared_statement stmt, int streaming, duckdb_result *res)
{
    duckdb_pending_result pending = NULL;
    duckdb_state st = streaming && DK(pending_prepared_streaming) ? DK(pending_prepared_streaming)(stmt, &pending)
                                                                  : DK(pending_prepared)(stmt, &pending);
    memset(res, 0, sizeof(*res));
    if (st != DuckDBSuccess) {
        const char *e = pending ? DK(pending_error)(pending) : NULL;
        snprintf(ctx->error, sizeof(ctx->error), "%s", e && *e ? e : "could not run the statement");
        if (pending) DK(destroy_pending)(&pending);
        return 198;
    }
    int idle = 0;
    for (;;) {
        duckdb_pending_state ps = DK(pending_execute_task)(pending);
        if (ps == DUCKDB_PENDING_RESULT_READY) break;
        if (ps == DUCKDB_PENDING_ERROR) {
            const char *e = DK(pending_error)(pending);
            snprintf(ctx->error, sizeof(ctx->error), "%s", e && *e ? e : "the statement failed");
            DK(destroy_pending)(&pending);
            return atomic_load(&ctx->cancelled) ? 1 : 198;
        }
        csql_service_once(ctx);
        if (ctools_cancel_requested()) csql_cancel(ctx);
        if (ps == DUCKDB_PENDING_NO_TASKS_AVAILABLE) {
            if (++idle < 256) csql_yield();
            else csql_wait(ctx);
        } else {
            idle = 0;
        }
    }
    st = DK(execute_pending)(pending, res);
    DK(destroy_pending)(&pending);
    if (st != DuckDBSuccess) {
        const char *e = DK(result_error)(res);
        snprintf(ctx->error, sizeof(ctx->error), "%s", e && *e ? e : "the statement failed");
        DK(destroy_result)(res);
        return atomic_load(&ctx->cancelled) ? 1 : 198;
    }
    return 0;
}

/* Run one statement of setup SQL and discard its result. */
static int csql_exec_sql(csql_ctx *ctx, const char *sql)
{
    duckdb_prepared_statement stmt = NULL;
    duckdb_result res;
    if (DK(prepare)(ctx->con, sql, &stmt) != DuckDBSuccess) {
        const char *e = stmt ? DK(prepare_error)(stmt) : NULL;
        snprintf(ctx->error, sizeof(ctx->error), "%s", e && *e ? e : "could not prepare a statement");
        DK(destroy_prepare)(&stmt);
        return 198;
    }
    int rc = csql_execute(ctx, stmt, 0, &res);
    DK(destroy_prepare)(&stmt);
    if (!rc) DK(destroy_result)(&res);
    return rc;
}

static int has_suffix(const char *s, const char *suffix)
{
    size_t n = strlen(s), m = strlen(suffix);
    if (m > n) return 0;
    for (size_t i = 0; i < m; i++) {
        char a = s[n - m + i], b = suffix[i];
        if (a >= 'A' && a <= 'Z') a = (char)(a - 'A' + 'a');
        if (a != b) return 0;
    }
    return 1;
}

/* The reader for csql's using file, by extension. */
static int using_reader(csql_buf *b, const char *path)
{
    const char *fn = NULL;
    char base[4096];
    snprintf(base, sizeof(base), "%s", path);
    /* compressed text: judge by the name before .gz/.zst */
    size_t n = strlen(base);
    if (has_suffix(base, ".gz")) base[n - 3] = '\0';
    else if (has_suffix(base, ".zst")) base[n - 4] = '\0';
    if (has_suffix(base, ".parquet") || has_suffix(base, ".parq") || has_suffix(base, ".pq")) fn = "read_parquet";
    else if (has_suffix(base, ".csv") || has_suffix(base, ".tsv") || has_suffix(base, ".txt") ||
             has_suffix(base, ".tab") || has_suffix(base, ".dat"))
        fn = "read_csv";
    else if (has_suffix(base, ".json") || has_suffix(base, ".jsonl") || has_suffix(base, ".ndjson"))
        fn = "read_json";
    else if (has_suffix(base, ".dta")) fn = "read_dta";
    else if (has_suffix(base, ".xlsx")) fn = "read_xlsx";
    int rc = 0;
    if (fn) {
        rc = csql_buf_str(b, fn);
        rc = rc ? rc : csql_buf_str(b, "(");
        rc = rc ? rc : csql_buf_literal(b, path);
        rc = rc ? rc : csql_buf_str(b, ")");
    } else {
        rc = csql_buf_literal(b, path);
    }
    return rc;
}

static void reset_used(csql_ctx *ctx)
{
    for (csql_source *s = ctx->sources; s; s = s->next) s->used = 0, s->order = 0;
    ctx->nused = 0;
}

/* ------------------------------------------------------------ run */

typedef struct
{
    int verbose, labels, datetime, readonly;
    int rows, linesize, strl;
    char copy[16]; /* parquet, csv, json: saving() through COPY */
} run_opts;

static void parse_run_opts(const char *args, run_opts *o)
{
    int v;
    memset(o, 0, sizeof(*o));
    o->labels = 1;
    o->datetime = 1;
    o->rows = 20;
    o->linesize = 80;
    o->strl = 2045;
    o->verbose = ctools_parse_bool_option(args, "verbose");
    o->readonly = ctools_parse_bool_option(args, "readonly");
    if (ctools_parse_int_option(args, "labels", &v) == CTOOLS_PARSE_OK) o->labels = v != 0;
    if (ctools_parse_int_option(args, "datetime", &v) == CTOOLS_PARSE_OK) o->datetime = v != 0;
    if (ctools_parse_int_option(args, "rows", &v) == CTOOLS_PARSE_OK && v >= 0) o->rows = v;
    if (ctools_parse_int_option(args, "linesize", &v) == CTOOLS_PARSE_OK && v > 0) o->linesize = v;
    if (ctools_parse_int_option(args, "strl", &v) == CTOOLS_PARSE_OK && v >= 1 && v <= 2045) o->strl = v;
    ctools_parse_string_option(args, "copy", o->copy, sizeof(o->copy));
}

typedef struct
{
    char *lib, *sqlfile, *meta, *out, *usingf, *database, *memory, *tempdir;
} run_locals;

static void free_locals(run_locals *l)
{
    free(l->lib), free(l->sqlfile), free(l->meta), free(l->out);
    free(l->usingf), free(l->database), free(l->memory), free(l->tempdir);
}

static int open_database(csql_ctx *ctx, const run_locals *l, const run_opts *o)
{
    duckdb_config config = NULL;
    char text[64], *err = NULL;
    if (DK(create_config)(&config) != DuckDBSuccess) return csql_fail("could not configure DuckDB", 920);
    snprintf(text, sizeof(text), "%d", ctx->threads);
    DK(set_config)(config, "threads", text);
    if (l->memory[0] && DK(set_config)(config, "memory_limit", l->memory) != DuckDBSuccess) {
        DK(destroy_config)(&config);
        return csql_fail("invalid memory() limit", 198);
    }
    if (l->tempdir[0]) DK(set_config)(config, "temp_directory", l->tempdir);
    if (o->readonly && l->database[0]) DK(set_config)(config, "access_mode", "READ_ONLY");
    duckdb_state st = DK(open_ext)(l->database[0] ? l->database : NULL, &ctx->db, config, &err);
    DK(destroy_config)(&config);
    if (st != DuckDBSuccess) {
        char msg[2400];
        snprintf(msg, sizeof(msg), "could not open the database%s%s: %s", l->database[0] ? " " : "",
                 l->database, err ? err : "unknown error");
        if (err) DK(free)(err);
        ctx->db = NULL;
        return csql_fail(msg, 603);
    }
    if (DK(connect)(ctx->db, &ctx->con) != DuckDBSuccess) {
        ctx->con = NULL;
        return csql_fail("could not connect to the database", 603);
    }
    return 0;
}

/* SELECT #1, CAST(#2 AS VARCHAR) AS "name", ... FROM (<last>) */
static int wrap_casts(duckdb_prepared_statement stmt, const char *last, size_t len, csql_buf *out)
{
    if (!DK(prepared_statement_column_count) || !DK(prepared_statement_column_type) ||
        !DK(prepared_statement_column_name))
        return 0;
    idx_t n = DK(prepared_statement_column_count)(stmt);
    int need = 0;
    for (idx_t j = 0; j < n && !need; j++)
        need = !csql_type_supported(DK(prepared_statement_column_type)(stmt, j));
    if (!need) return 0;
    int rc = csql_buf_str(out, "SELECT ");
    for (idx_t j = 0; j < n && !rc; j++) {
        char pos[32];
        const char *name = DK(prepared_statement_column_name)(stmt, j);
        int ok = csql_type_supported(DK(prepared_statement_column_type)(stmt, j));
        snprintf(pos, sizeof(pos), "#%llu", (unsigned long long)j + 1);
        rc = csql_buf_str(out, j ? ", " : "");
        if (!ok) rc = rc ? rc : csql_buf_str(out, "CAST(");
        rc = rc ? rc : csql_buf_str(out, pos);
        if (!ok) rc = rc ? rc : csql_buf_str(out, " AS VARCHAR)");
        rc = rc ? rc : csql_buf_str(out, " AS ");
        rc = rc ? rc : csql_buf_ident(out, name ? name : "");
        if (name) DK(free)((void *)name);
    }
    rc = rc ? rc : csql_buf_str(out, " FROM (");
    rc = rc ? rc : csql_buf_add(out, last, len);
    rc = rc ? rc : csql_buf_str(out, "\n) AS __csql_q");
    return rc ? -1 : 1;
}

static int csql_run(const char *args)
{
    csql_ctx ctx;
    run_opts o;
    run_locals l = {0};
    char msg[4400] = "";
    int rc = 0;
    double t0 = ctools_timer_seconds();
    duckdb_extracted_statements ext = NULL;
    char *sql = NULL;
    int64_t nrows = 0, ncols = 0, changed = -1;

    memset(&ctx, 0, sizeof(ctx));
    parse_run_opts(args, &o);
    l.lib = csql_local("__csql_lib");
    l.sqlfile = csql_local("__csql_sqlfile");
    l.meta = csql_local("__csql_meta");
    l.out = csql_local("__csql_out");
    l.usingf = csql_local("__csql_using");
    l.database = csql_local("__csql_database");
    l.memory = csql_local("__csql_memory");
    l.tempdir = csql_local("__csql_tempdir");
    if (!l.lib || !l.sqlfile || !l.meta || !l.out || !l.usingf || !l.database || !l.memory || !l.tempdir) {
        free_locals(&l);
        return csql_fail("out of memory", 920);
    }
    rc = csql_api_load(l.lib, msg, sizeof(msg));
    if (rc) {
        free_locals(&l);
        return csql_fail(msg, rc);
    }
    sql = read_text(l.sqlfile, msg, sizeof(msg));
    if (!sql) {
        free_locals(&l);
        return csql_fail(msg, 601);
    }

    ctx.main_thread = csql_thread_id();
    ctx.threads = ctools_get_max_threads();
    if (ctx.threads < 1) ctx.threads = 1;
    ctx.datetime = o.datetime;
    ctx.verbose = o.verbose;
    atomic_init(&ctx.cancelled, 0);
    pthread_mutex_init(&ctx.lock, NULL);
    pthread_cond_init(&ctx.wake, NULL);
    pthread_cond_init(&ctx.served, NULL);
    ctx.service_ready = 1;

    if (l.meta[0]) {
        rc = csql_meta_read(l.meta, &ctx.stata, msg, sizeof(msg));
        if (rc) {
            csql_fail(msg, rc);
            goto done;
        }
        if ((ST_int)ctx.stata->nvars != SF_nvars()) {
            rc = csql_fail("the description of the data in memory does not match the variables", 198);
            goto done;
        }
        ctx.stata_nobs = (int64_t)SF_nobs();
        ctx.stata->next = ctx.sources;
        ctx.sources = ctx.stata;
    }

    rc = open_database(&ctx, &l, &o);
    if (rc) goto done;
    if ((ctx.stata && (rc = csql_stata_register(&ctx)) != 0) || (rc = csql_dta_register(&ctx)) != 0) {
        csql_fail("could not register csql's table functions", rc);
        goto done;
    }

    /* views: stata (the data in memory) and data (the using file, else stata) */
    {
        csql_buf b = {0};
        if (ctx.stata && ctx.stata->nvars) {
            csql_buf_str(&b, "CREATE OR REPLACE TEMP VIEW stata AS SELECT * FROM __csql_stata()");
            rc = csql_exec_sql(&ctx, b.p);
            b.len = 0;
        }
        if (!rc && (l.usingf[0] || (ctx.stata && ctx.stata->nvars))) {
            csql_buf_str(&b, "CREATE OR REPLACE TEMP VIEW data AS SELECT * FROM ");
            if (l.usingf[0]) using_reader(&b, l.usingf);
            else csql_buf_str(&b, "stata");
            rc = csql_exec_sql(&ctx, b.p);
        }
        csql_buf_free(&b);
        if (rc) {
            char text[2400];
            snprintf(text, sizeof(text), "%s%s", l.usingf[0] ? "could not read the using file: " : "", ctx.error);
            csql_fail(text, rc == 1 ? 1 : l.usingf[0] ? 601 : 198);
            if (rc != 1) rc = l.usingf[0] ? 601 : 198;
            goto done;
        }
        reset_used(&ctx);
    }
    double t_setup = ctools_timer_seconds();

    idx_t nstmt = DK(extract_statements)(ctx.con, sql, &ext);
    if (nstmt == 0) {
        const char *e = ext ? DK(extract_statements_error)(ext) : NULL;
        rc = csql_fail(e && *e ? e : "the query has no SQL statement", 198);
        goto done;
    }
    for (idx_t i = 0; i < nstmt && !rc; i++) {
        int last = i + 1 == nstmt;
        duckdb_prepared_statement stmt = NULL;
        int binds_before = ctx.stata_binds;
        char prefix[64] = "";
        if (nstmt > 1) snprintf(prefix, sizeof(prefix), "statement %llu: ", (unsigned long long)i + 1);
        if (DK(prepare_extracted_statement)(ctx.con, ext, i, &stmt) != DuckDBSuccess) {
            const char *e = stmt ? DK(prepare_error)(stmt) : NULL;
            snprintf(ctx.error, sizeof(ctx.error), "%s%s", prefix, e && *e ? e : "invalid statement");
            DK(destroy_prepare)(&stmt);
            rc = csql_fail(ctx.error, 198);
            break;
        }
        if (!last) {
            duckdb_result res;
            rc = csql_execute(&ctx, stmt, 0, &res);
            DK(destroy_prepare)(&stmt);
            if (rc) {
                char text[2200];
                snprintf(text, sizeof(text), "%s%s", prefix, ctx.error);
                csql_fail(text, rc);
            } else {
                DK(destroy_result)(&res);
            }
            continue;
        }

        /* the last statement: its rows go to the destination */
        duckdb_statement_type stype = DK(prepared_statement_type)(stmt);
        size_t s0 = 0, s1 = 0;
        int have_text = csql_sql_last_statement(sql, &s0, &s1) == 0;
        csql_buf alt = {0};
        if (o.copy[0]) {
            if (stype != DUCKDB_STATEMENT_TYPE_SELECT || !have_text) {
                DK(destroy_prepare)(&stmt);
                rc = csql_fail("saving() needs a query: the last statement must return rows", 198);
                break;
            }
            csql_buf_str(&alt, "COPY (");
            csql_buf_add(&alt, sql + s0, s1 - s0);
            csql_buf_str(&alt, "\n) TO ");
            csql_buf_literal(&alt, l.out);
            csql_buf_str(&alt, " (FORMAT ");
            csql_buf_str(&alt, o.copy);
            csql_buf_str(&alt, ")");
        } else if (stype == DUCKDB_STATEMENT_TYPE_SELECT && have_text) {
            if (wrap_casts(stmt, sql + s0, s1 - s0, &alt) < 0) {
                DK(destroy_prepare)(&stmt);
                rc = csql_fail("out of memory", 920);
                break;
            }
        }
        if (alt.len) {
            DK(destroy_prepare)(&stmt);
            binds_before = ctx.stata_binds;
            if (DK(prepare)(ctx.con, alt.p, &stmt) != DuckDBSuccess) {
                const char *e = stmt ? DK(prepare_error)(stmt) : NULL;
                snprintf(ctx.error, sizeof(ctx.error), "%s%s", prefix, e && *e ? e : "invalid statement");
                DK(destroy_prepare)(&stmt);
                csql_buf_free(&alt);
                rc = csql_fail(ctx.error, 198);
                break;
            }
        }
        csql_buf_free(&alt);
        /* Rows for Stata: a materialized result (built by every thread),
         * read chunk by chunk, so DuckDB's copy is the only one. A listing
         * streams: it keeps only the rows it shows. A streaming result has
         * the Stata thread wait inside DuckDB, so strL reads from the data
         * in memory, which need that thread, rule it out. */
        int strl_mem = 0;
        if (ctx.stata && ctx.stata_binds != binds_before)
            for (size_t j = 0; j < ctx.stata->nvars && !strl_mem; j++)
                strl_mem = ctx.stata->vars[j].type == CTOOLS_DTA_CODE_STRL;
        int random_ok = DK(result_get_chunk) && DK(result_chunk_count);
        int streaming = !o.copy[0] && !strl_mem && (!l.out[0] || !random_ok);
        duckdb_result res;
        rc = csql_execute(&ctx, stmt, streaming, &res);
        DK(destroy_prepare)(&stmt);
        if (rc) {
            char text[2200];
            snprintf(text, sizeof(text), "%s%s", prefix, ctx.error);
            csql_fail(text, rc);
            break;
        }
        duckdb_result_type rt = DK(result_return_type)(res);
        if (o.copy[0] || rt != DUCKDB_RESULT_TYPE_QUERY_RESULT) {
            if (o.copy[0] || rt == DUCKDB_RESULT_TYPE_CHANGED_ROWS) changed = (int64_t)DK(rows_changed)(&res);
            if (!o.copy[0] && l.out[0]) {
                rc = csql_fail("the last statement returns no rows to load", 198);
            } else if (o.copy[0]) {
                nrows = changed; /* the ado reports the file */
            } else if (changed >= 0) {
                char text[128];
                snprintf(text, sizeof(text), "{txt}(%lld row%s changed)\n", (long long)changed, changed == 1 ? "" : "s");
                SF_display(text);
            }
            DK(destroy_result)(&res);
            continue;
        }
        csql_out_opts out = {0};
        out.mode = l.out[0] ? CSQL_OUT_DTA : CSQL_OUT_DISPLAY;
        out.path = l.out;
        out.temporary = ctools_parse_bool_option(args, "temporary");
        out.labels = o.labels;
        out.rows = o.rows;
        out.linesize = o.linesize;
        out.strl_min = o.strl;
        ctx.error[0] = '\0';
        rc = csql_result_output(&ctx, &res, streaming, &out, &nrows, &ncols);
        DK(destroy_result)(&res);
        if (rc) {
            if (rc == 1 && !ctx.error[0]) snprintf(ctx.error, sizeof(ctx.error), "--Break--");
            if (rc != 1 || ctx.error[0]) csql_fail(ctx.error[0] ? ctx.error : "could not store the result", rc);
        }
    }
    if (o.verbose && !rc)
        ctools_msg("csql", "DuckDB %s, %d threads: setup %.3f s, total %.3f s", csql_db.version, ctx.threads,
                   t_setup - t0, ctools_timer_seconds() - t0);

done:
    if (ext) DK(destroy_extracted)(&ext);
    if (ctx.con) DK(disconnect)(&ctx.con);
    if (ctx.db) DK(close)(&ctx.db);
    csql_scan_pool_free(&ctx);
    csql_dta_close_all();
    while (ctx.sources) {
        csql_source *n = ctx.sources->next;
        csql_source_free(ctx.sources);
        ctx.sources = n;
    }
    pthread_cond_destroy(&ctx.wake);
    pthread_cond_destroy(&ctx.served);
    pthread_mutex_destroy(&ctx.lock);
    free(sql);
    free_locals(&l);
    if (!rc) {
        char text[64];
        SF_scal_save("__csql_N", (double)nrows);
        SF_scal_save("__csql_k", (double)ncols);
        SF_scal_save("__csql_changed", changed >= 0 ? (double)changed : SV_missval);
        snprintf(text, sizeof(text), "%s", csql_db.version);
        SF_macro_save("_csql_version", text);
    }
    return rc;
}

/* ------------------------------------------------------------ entry */

ST_retcode csql_main(const char *args)
{
    char op[32] = "";
    char msg[4400] = "";
    sscanf(args ? args : "", "%31s", op);
    if (!strcmp(op, "run")) return csql_run(args + 3);
    if (!strcmp(op, "version")) {
        char *lib = csql_local("__csql_lib");
        int rc = lib ? csql_api_load(lib, msg, sizeof(msg)) : 920;
        free(lib);
        if (rc) return csql_fail(msg[0] ? msg : "out of memory", rc);
        SF_macro_save("_csql_version", csql_db.version);
        return 0;
    }
    if (!strcmp(op, "install")) {
        char *archive = csql_local("__csql_archive"), *sha = csql_local("__csql_sha256");
        char *member = csql_local("__csql_member"), *dest = csql_local("__csql_dest");
        int rc = archive && sha && member && dest ? csql_install_extract(archive, sha, member, dest, msg, sizeof(msg)) : 920;
        free(archive), free(sha), free(member), free(dest);
        if (rc) return csql_fail(msg[0] ? msg : "out of memory", rc);
        return 0;
    }
    return csql_fail("unknown subcommand", 198);
}

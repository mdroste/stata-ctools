/*
 * csql_stata_scan.c - the data in memory as a DuckDB table function.
 *
 * __csql_stata() returns one column per variable passed to the plugin, in
 * dataset order, and one row per observation. The views `stata` (and `data`
 * when csql has no using file) select from it. Numeric missing values read
 * as NULL; strings keep "" (Stata has no NULL string). byte, int, long,
 * float and double read as TINYINT, SMALLINT, INTEGER, FLOAT and DOUBLE;
 * variables with a %td or %tc format read as DATE and TIMESTAMP.
 *
 * Only the columns a query uses are read (projection pushdown). DuckDB calls
 * the scan from several threads; each call takes the next block of rows.
 * strL values go through the Stata thread (csql_read_strls).
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "stplugin.h"
#include "csql/csql_internal.h"
#include "csql/csql_scan.h"

#ifdef _WIN32
#include <windows.h>
#endif

/* ------------------------------------------------------------ mapping */

static int fmt_is(const char *format, const char *code)
{
    /* %td, %-td, %tdD_m_Y ...; %d (and %dD_m_Y) is Stata's old name for %td */
    const char *p = format;
    if (*p != '%') return 0;
    p++;
    if (*p == '-') p++;
    if (strcmp(code, "td") == 0 && p[0] == 'd') return 1;
    return strncmp(p, code, strlen(code)) == 0;
}

csql_map csql_map_stata_var(uint16_t type, const char *format, int datetime)
{
    if (type == CTOOLS_DTA_CODE_STRL || (type >= 1 && type <= CSQL_STR_MAX)) return CSQL_MAP_VARCHAR;
    if (datetime && format) {
        if (fmt_is(format, "td")) return CSQL_MAP_DATE;
        if (fmt_is(format, "tc")) return CSQL_MAP_TIMESTAMP;
    }
    switch (type) {
    case CTOOLS_DTA_CODE_BYTE: return CSQL_MAP_TINYINT;
    case CTOOLS_DTA_CODE_INT: return CSQL_MAP_SMALLINT;
    case CTOOLS_DTA_CODE_LONG: return CSQL_MAP_INTEGER;
    case CTOOLS_DTA_CODE_FLOAT: return CSQL_MAP_FLOAT;
    default: return CSQL_MAP_DOUBLE;
    }
}

duckdb_logical_type csql_map_logical(csql_map m)
{
    duckdb_type t;
    switch (m) {
    case CSQL_MAP_TINYINT: t = DUCKDB_TYPE_TINYINT; break;
    case CSQL_MAP_SMALLINT: t = DUCKDB_TYPE_SMALLINT; break;
    case CSQL_MAP_INTEGER: t = DUCKDB_TYPE_INTEGER; break;
    case CSQL_MAP_FLOAT: t = DUCKDB_TYPE_FLOAT; break;
    case CSQL_MAP_DATE: t = DUCKDB_TYPE_DATE; break;
    case CSQL_MAP_TIMESTAMP: t = DUCKDB_TYPE_TIMESTAMP; break;
    case CSQL_MAP_VARCHAR: t = DUCKDB_TYPE_VARCHAR; break;
    default: t = DUCKDB_TYPE_DOUBLE; break;
    }
    return DK(create_logical_type)(t);
}

void csql_map_store(csql_map m, void *data, size_t row, double z, int *invalid)
{
    *invalid = 0;
    switch (m) {
    case CSQL_MAP_TINYINT: ((int8_t *)data)[row] = (int8_t)z; break;
    case CSQL_MAP_SMALLINT: ((int16_t *)data)[row] = (int16_t)z; break;
    case CSQL_MAP_INTEGER: ((int32_t *)data)[row] = (int32_t)z; break;
    case CSQL_MAP_FLOAT: ((float *)data)[row] = (float)z; break;
    case CSQL_MAP_DATE: {
        double d = floor(z) - CSQL_DAYS_1960_1970;
        if (d < -2147483647.0 || d > 2147483646.0) *invalid = 1;
        else ((duckdb_date *)data)[row].days = (int32_t)d;
        break;
    }
    case CSQL_MAP_TIMESTAMP: {
        double us = (z - (double)CSQL_MS_1960_1970) * 1000.0;
        if (!(us > -9.2e18 && us < 9.2e18)) *invalid = 1;
        else ((duckdb_timestamp *)data)[row].micros = (int64_t)(us < 0 ? us - 0.5 : us + 0.5);
        break;
    }
    default: ((double *)data)[row] = z; break;
    }
}

/* ------------------------------------------------------------ threads */

uint64_t csql_thread_id(void)
{
#ifdef _WIN32
    return (uint64_t)GetCurrentThreadId();
#else
    return (uint64_t)(uintptr_t)pthread_self();
#endif
}

int csql_on_main(const csql_ctx *ctx)
{
    return csql_thread_id() == ctx->main_thread;
}

static int strl_read_direct(int var, int64_t row0, size_t n, duckdb_vector out)
{
    char small[256];
    for (size_t i = 0; i < n; i++) {
        ST_int obs = (ST_int)(row0 + (int64_t)i);
        ST_int len = SF_sdatalen(var, obs);
        if (len < 0) return 3;
        char *buf = (size_t)len + 1 <= sizeof(small) ? small : malloc((size_t)len + 1);
        if (!buf) return 920;
        if (len > 0) SF_strldata(var, obs, buf, len + 1);
        DK(vector_assign_string_element_len)(out, (idx_t)i, buf, (idx_t)len);
        if (buf != small) free(buf);
    }
    return 0;
}

int csql_read_strls(csql_ctx *ctx, int var, int64_t row0, size_t n, duckdb_vector out)
{
    csql_strl_request req;
    if (csql_on_main(ctx)) return strl_read_direct(var, row0, n, out);
    if (!ctx->service_ready) return 198;
    req.var = var;
    req.row0 = row0;
    req.n = n;
    req.out = out;
    req.rc = 0;
    req.done = 0;
    pthread_mutex_lock(&ctx->lock);
    while (ctx->pending && !atomic_load(&ctx->cancelled)) pthread_cond_wait(&ctx->served, &ctx->lock);
    if (atomic_load(&ctx->cancelled)) {
        pthread_mutex_unlock(&ctx->lock);
        return 1;
    }
    ctx->pending = &req;
    pthread_cond_signal(&ctx->wake);
    while (!req.done && !atomic_load(&ctx->cancelled)) pthread_cond_wait(&ctx->served, &ctx->lock);
    if (!req.done) {
        /* cancelled before the Stata thread took the request */
        if (ctx->pending == &req) ctx->pending = NULL;
        pthread_cond_broadcast(&ctx->served);
        pthread_mutex_unlock(&ctx->lock);
        return 1;
    }
    pthread_mutex_unlock(&ctx->lock);
    return req.rc;
}

void csql_service_once(csql_ctx *ctx)
{
    csql_strl_request *req;
    if (!ctx->service_ready) return;
    pthread_mutex_lock(&ctx->lock);
    req = ctx->pending;
    pthread_mutex_unlock(&ctx->lock);
    if (!req) return;
    int rc = strl_read_direct(req->var, req->row0, req->n, req->out);
    pthread_mutex_lock(&ctx->lock);
    req->rc = rc;
    req->done = 1;
    ctx->pending = NULL;
    pthread_cond_broadcast(&ctx->served);
    pthread_mutex_unlock(&ctx->lock);
}

/* ------------------------------------------------------------ scan */

typedef struct
{
    csql_scan_def def; /* first: the engine reads it */
    csql_map *map;     /* per variable */
    int *strl;         /* per variable */
    int *width;        /* per variable: str# width */
} stata_bind;

static void stata_bind_free(void *p)
{
    stata_bind *b = p;
    if (!b) return;
    free(b->map);
    free(b->strl);
    free(b->width);
    free(b);
}

/* Rows [pre->r0, pre->r0 + pre->n) through the SPI (any thread; strL
 * columns are deferred to the engine). */
static int stata_decode(void *arg, const idx_t *proj, size_t nproj, csql_pre *pre)
{
    stata_bind *b = arg;
    size_t n = pre->n, need = 0;
    char sbuf[CSQL_STR_MAX + 2];
    for (size_t c = 0; c < nproj; c++) {
        idx_t j = proj[c];
        if (j < b->def.nsrc && b->map[j] == CSQL_MAP_VARCHAR && !b->strl[j] && b->width[j] > 12)
            need += n * (size_t)b->width[j];
    }
    if (need && csql_arena_reserve(&pre->arena, need) != 0) return 920;
    for (size_t c = 0; c < nproj; c++) {
        idx_t j = proj[c];
        if (j >= b->def.nsrc || b->strl[j]) continue;
        ST_int var = (ST_int)(j + 1);
        ST_int obs = (ST_int)(pre->r0 + 1);
        uint64_t *valid = pre->valid[c];
        if (b->map[j] == CSQL_MAP_VARCHAR) {
            duckdb_string_t *out = pre->data[c];
            for (size_t i = 0; i < n; i++) {
                sbuf[0] = '\0';
                SF_sdata(var, obs + (ST_int)i, sbuf);
                size_t len = strlen(sbuf);
                const char *p = sbuf;
                if (len > 12) {
                    p = csql_arena_add(&pre->arena, sbuf, len);
                    if (!p) return 920;
                    pre->copy[c] = 1;
                }
                csql_str_set(&out[i], p, (uint32_t)len);
                csql_valid_set(valid, i);
            }
            continue;
        }
        void *data = pre->data[c];
        csql_map m = b->map[j];
        for (size_t i = 0; i < n; i++) {
            ST_double z = 0;
            int bad = 1;
            SF_vdata(var, obs + (ST_int)i, &z);
            if (!csql_is_missing(z)) csql_map_store(m, data, i, z, &bad);
            if (!bad) csql_valid_set(valid, i);
        }
    }
    return 0;
}

static void stata_scan_bind(duckdb_bind_info info)
{
    csql_ctx *ctx = DK(bind_get_extra_info)(info);
    csql_source *src = ctx->stata;
    if (!src || src->nvars == 0) {
        DK(bind_set_error)(info, "there are no variables in memory");
        return;
    }
    stata_bind *b = calloc(1, sizeof(*b));
    if (b) {
        b->map = calloc(src->nvars, sizeof(*b->map));
        b->strl = calloc(src->nvars, sizeof(*b->strl));
        b->width = calloc(src->nvars, sizeof(*b->width));
    }
    if (!b || !b->map || !b->strl || !b->width) {
        stata_bind_free(b);
        DK(bind_set_error)(info, "out of memory");
        return;
    }
    for (size_t j = 0; j < src->nvars; j++) {
        const csql_meta_var *v = &src->vars[j];
        b->map[j] = csql_map_stata_var(v->type, v->format, ctx->datetime);
        b->strl[j] = v->type == CTOOLS_DTA_CODE_STRL;
        b->width[j] = v->type <= CSQL_STR_MAX ? v->type : 0;
        duckdb_logical_type t = csql_map_logical(b->map[j]);
        DK(bind_add_result_column)(info, v->sqlname, t);
        DK(destroy_logical_type)(&t);
    }
    b->def = (csql_scan_def){ctx, b, ctx->stata_nobs, src->nvars, b->map, b->strl, 0, stata_decode};
    DK(bind_set_cardinality)(info, (idx_t)ctx->stata_nobs, true);
    DK(bind_set_bind_data)(info, b, stata_bind_free);
    ctx->stata_binds++;
    if (!src->used) {
        src->used = 1;
        src->order = ++ctx->nused;
    }
}

int csql_stata_register(csql_ctx *ctx)
{
    duckdb_table_function f = DK(create_table_function)();
    duckdb_state st;
    if (!f) return 920;
    DK(table_function_set_name)(f, "__csql_stata");
    DK(table_function_set_extra_info)(f, ctx, NULL);
    DK(table_function_set_bind)(f, stata_scan_bind);
    DK(table_function_set_init)(f, csql_scan_init);
    DK(table_function_set_local_init)(f, csql_scan_local_init);
    DK(table_function_set_function)(f, csql_scan_function);
    DK(table_function_supports_projection_pushdown)(f, true);
    st = DK(register_table_function)(ctx->con, f);
    DK(destroy_table_function)(&f);
    return st == DuckDBSuccess ? 0 : 198;
}

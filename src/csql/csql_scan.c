/*
 * csql_scan.c - batching, prefetch and hand-off shared by csql's table
 * functions (see csql_scan.h).
 */
#include <stdlib.h>
#include <string.h>

#include "csql/csql_scan.h"
#include "ctools_threads.h"

#define CSQL_BATCH 32 /* chunks decoded ahead by one reader: 64K rows */

size_t csql_map_width(csql_map m)
{
    switch (m) {
    case CSQL_MAP_TINYINT: return 1;
    case CSQL_MAP_SMALLINT: return 2;
    case CSQL_MAP_INTEGER: case CSQL_MAP_FLOAT: case CSQL_MAP_DATE: return 4;
    case CSQL_MAP_VARCHAR: return sizeof(duckdb_string_t);
    default: return 8;
    }
}

void csql_str_set(duckdb_string_t *s, const char *p, uint32_t len)
{
    memset(s, 0, sizeof(*s));
    s->value.inlined.length = len;
    if (len <= 12) {
        if (len) memcpy(s->value.inlined.inlined, p, len);
    } else {
        memcpy(s->value.pointer.prefix, p, 4);
        s->value.pointer.ptr = (char *)p;
    }
}

int csql_arena_reserve(csql_arena *a, size_t bytes)
{
    if (a->len + bytes <= a->cap) return 0;
    if (a->len) return 920; /* never move strings already handed out */
    size_t cap = a->cap ? a->cap : (size_t)1 << 16;
    while (cap < bytes) cap *= 2;
    char *q = realloc(a->p, cap);
    if (!q) return 920;
    a->p = q;
    a->cap = cap;
    return 0;
}

char *csql_arena_add(csql_arena *a, const char *p, size_t n)
{
    if (a->len + n > a->cap) return NULL;
    char *dst = a->p + a->len;
    memcpy(dst, p, n);
    a->len += n;
    return dst;
}

/* ------------------------------------------------------------ states */

typedef struct
{
    atomic_llong next; /* next row to claim */
    atomic_int nlocal; /* DuckDB threads reading this scan */
    size_t nproj;
    idx_t *proj;
} scan_global;

typedef struct
{
    csql_pre pre[CSQL_BATCH];
    size_t npre, serve;
    size_t cap; /* rows per chunk */
    size_t nproj;
    int registered; /* counted in nlocal (on its first call) */
} scan_local;

static void pre_free(csql_pre *p, size_t nproj)
{
    for (size_t c = 0; c < nproj; c++) {
        if (p->data) free(p->data[c]);
        if (p->valid) free(p->valid[c]);
    }
    free(p->data);
    free(p->valid);
    free(p->nulls);
    free(p->copy);
    free(p->arena.p);
    memset(p, 0, sizeof(*p));
}

static int pre_alloc(csql_pre *p, const csql_scan_def *def, const scan_global *g, size_t cap)
{
    size_t k = g->nproj ? g->nproj : 1, words = (cap + 63) / 64;
    p->data = calloc(k, sizeof(void *));
    p->valid = calloc(k, sizeof(uint64_t *));
    p->nulls = calloc(k, 1);
    p->copy = calloc(k, 1);
    if (!p->data || !p->valid || !p->nulls || !p->copy) return 920;
    for (size_t c = 0; c < g->nproj; c++) {
        idx_t j = g->proj[c];
        size_t w = j < def->nsrc ? csql_map_width(def->map[j]) : 8;
        p->data[c] = malloc(cap * w);
        p->valid[c] = malloc(words * sizeof(uint64_t));
        if (!p->data[c] || !p->valid[c]) return 920;
    }
    return 0;
}

static void global_free(void *p)
{
    scan_global *g = p;
    if (!g) return;
    free(g->proj);
    free(g);
}

static void local_free(void *p)
{
    scan_local *l = p;
    if (!l) return;
    for (int i = 0; i < CSQL_BATCH; i++) pre_free(&l->pre[i], l->nproj);
    free(l);
}

void csql_scan_init(duckdb_init_info info)
{
    const csql_scan_def *def = DK(init_get_bind_data)(info);
    scan_global *g = calloc(1, sizeof(*g));
    idx_t n = DK(init_get_column_count)(info);
    if (g) g->proj = calloc(n ? n : 1, sizeof(idx_t));
    if (!g || !g->proj) {
        global_free(g);
        DK(init_set_error)(info, "out of memory");
        return;
    }
    atomic_init(&g->next, 0);
    atomic_init(&g->nlocal, 0);
    g->nproj = n;
    for (idx_t i = 0; i < n; i++) g->proj[i] = DK(init_get_column_index)(info, i);
    DK(init_set_max_threads)(info, (idx_t)(def->ctx->threads > 0 ? def->ctx->threads : 1));
    DK(init_set_init_data)(info, g, global_free);
}

/* One per DuckDB thread reading the scan; it counts itself among the
 * readers on its first call (local init cannot see the global state). */
void csql_scan_local_init(duckdb_init_info info)
{
    scan_local *l = calloc(1, sizeof(*l));
    if (!l) {
        DK(init_set_error)(info, "out of memory");
        return;
    }
    DK(init_set_init_data)(info, l, local_free);
}

/* ------------------------------------------------------------ decoding */

typedef struct
{
    const csql_scan_def *def;
    const scan_global *g;
    csql_pre *pre;
} decode_job;

static void decode_one(const csql_scan_def *def, const scan_global *g, csql_pre *pre)
{
    size_t words = (pre->n + 63) / 64;
    for (size_t c = 0; c < g->nproj; c++) {
        memset(pre->valid[c], 0, words * sizeof(uint64_t));
        pre->copy[c] = 0;
    }
    pre->arena.len = 0;
    pre->rc = def->decode(def->bind, g->proj, g->nproj, pre);
    for (size_t c = 0; c < g->nproj; c++) {
        size_t valid = 0;
        for (size_t w = 0; w < words; w++) valid += (size_t)__builtin_popcountll(pre->valid[c][w]);
        pre->nulls[c] = valid != pre->n;
    }
}

static void *decode_worker(void *arg)
{
    decode_job *j = arg;
    decode_one(j->def, j->g, j->pre);
    return NULL;
}

static ctools_persistent_pool *scan_pool(csql_ctx *ctx)
{
    pthread_mutex_lock(&ctx->lock);
    if (!ctx->pool && !ctx->pool_failed && ctx->threads > 1) {
        ctools_persistent_pool *p = calloc(1, sizeof(*p));
        if (p && ctools_persistent_pool_init(p, (size_t)ctx->threads) == 0) ctx->pool = p;
        else {
            free(p);
            ctx->pool_failed = 1;
        }
    }
    ctools_persistent_pool *pool = ctx->pool;
    pthread_mutex_unlock(&ctx->lock);
    return pool;
}

void csql_scan_pool_free(csql_ctx *ctx)
{
    if (!ctx->pool) return;
    ctools_persistent_pool_destroy(ctx->pool);
    free(ctx->pool);
    ctx->pool = NULL;
}

/* ------------------------------------------------------------ hand-off */

static int serve(const csql_scan_def *def, const scan_global *g, const csql_pre *pre,
                 duckdb_data_chunk output)
{
    size_t n = pre->n, words = (n + 63) / 64;
    for (size_t c = 0; c < g->nproj; c++) {
        duckdb_vector vec = DK(data_chunk_get_vector)(output, (idx_t)c);
        idx_t j = g->proj[c];
        if (j >= def->nsrc) {
            DK(vector_ensure_validity_writable)(vec);
            uint64_t *v = DK(vector_get_validity)(vec);
            memset(v, 0, words * sizeof(uint64_t));
            continue;
        }
        if (def->deferred && def->deferred[j]) {
            int rc = csql_read_strls(def->ctx, (int)j + 1, pre->r0 + 1, n, vec);
            if (rc) return rc;
            continue;
        }
        if (def->map[j] == CSQL_MAP_VARCHAR) {
            const duckdb_string_t *src = pre->data[c];
            duckdb_string_t *dst = DK(vector_get_data)(vec);
            int direct = def->stable && !pre->copy[c];
            for (size_t i = 0; i < n; i++) {
                if (pre->nulls[c] && !((pre->valid[c][i >> 6] >> (i & 63)) & 1)) continue;
                uint32_t len = src[i].value.inlined.length;
                if (len <= 12 || direct) dst[i] = src[i];
                else DK(vector_assign_string_element_len)(vec, (idx_t)i, src[i].value.pointer.ptr, len);
            }
        } else {
            memcpy(DK(vector_get_data)(vec), pre->data[c], n * csql_map_width(def->map[j]));
        }
        if (pre->nulls[c]) {
            DK(vector_ensure_validity_writable)(vec);
            memcpy(DK(vector_get_validity)(vec), pre->valid[c], words * sizeof(uint64_t));
        }
    }
    DK(data_chunk_set_size)(output, (idx_t)n);
    return 0;
}

void csql_scan_function(duckdb_function_info info, duckdb_data_chunk output)
{
    const csql_scan_def *def = DK(function_get_bind_data)(info);
    scan_global *g = DK(function_get_init_data)(info);
    scan_local *l = DK(function_get_local_init_data)(info);
    csql_ctx *ctx = def->ctx;
    if (!l->registered) {
        l->registered = 1;
        l->nproj = g->nproj;
        l->cap = (size_t)DK(vector_size)();
        atomic_fetch_add(&g->nlocal, 1);
    }
    if (l->serve >= l->npre) {
        ctools_persistent_pool *pool = NULL;
        size_t batch = 1;
        if (atomic_load(&g->nlocal) == 1 && ctx->threads > 1 && (pool = scan_pool(ctx)) != NULL)
            batch = CSQL_BATCH;
        long long r0 = atomic_fetch_add(&g->next, (long long)(batch * l->cap));
        l->npre = l->serve = 0;
        if (r0 >= def->nrows || atomic_load(&ctx->cancelled)) {
            DK(data_chunk_set_size)(output, 0);
            return;
        }
        size_t left = (size_t)(def->nrows - r0), k = 0;
        for (; k < batch && k * l->cap < left; k++) {
            csql_pre *p = &l->pre[k];
            if (!p->data && pre_alloc(p, def, g, l->cap) != 0) {
                DK(function_set_error)(info, "out of memory");
                return;
            }
            p->r0 = r0 + (long long)(k * l->cap);
            p->n = left - k * l->cap < l->cap ? left - k * l->cap : l->cap;
        }
        if (k > 1 && pool) {
            decode_job jobs[CSQL_BATCH];
            for (size_t i = 0; i < k; i++) jobs[i] = (decode_job){def, g, &l->pre[i]};
            if (ctools_persistent_pool_submit_batch(pool, decode_worker, jobs, k, sizeof(jobs[0])) != 0)
                for (size_t i = 0; i < k; i++) decode_one(def, g, &l->pre[i]);
            else
                ctools_persistent_pool_wait(pool);
        } else {
            for (size_t i = 0; i < k; i++) decode_one(def, g, &l->pre[i]);
        }
        for (size_t i = 0; i < k; i++)
            if (l->pre[i].rc) {
                DK(function_set_error)(info, l->pre[i].rc == 920 ? "out of memory" : "could not read the data");
                return;
            }
        l->npre = k;
    }
    int rc = serve(def, g, &l->pre[l->serve++], output);
    if (rc) DK(function_set_error)(info, rc == 1 ? "cancelled" : "could not read a strL value");
}

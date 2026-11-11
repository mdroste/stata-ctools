/* cimport parquet: plugin entry points.
 *
 *   cparquet schema         footer: rows, columns, names and type names
 *   cparquet describe       the same; then every page is decoded, and an
 *                           unreadable file reports -1 rows, as Stata's
 *                           describe does
 *   cparquet scan           decodes the selected columns over the row range
 *                           (in parallel) and reports each one's Stata
 *                           storage type, format and dictionary size
 *   cparquet load           stores numeric and str# columns into the
 *                           variables named on the plugin call
 *   cparquet blobfile j     writes column j's values, length-prefixed, for
 *                           Mata to store (strL and text holding NUL bytes)
 *   cparquet dictfile j     writes column j's dictionary entries, for its
 *                           value labels
 *   cparquet clear          frees the scanned columns
 *
 * Inputs and results are globals: CTOOLS_PQ_file (absolute path),
 * CTOOLS_PQ_cols (0-based field indices), CTOOLS_PQ_first and CTOOLS_PQ_last
 * (1-based rows); the ado file drops them all when it is done. */
#include "cparquet.h"
#include "parquet_internal.h"
#include "ctools_file.h"
#include "ctools_threads.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef enum
{
    ST_BYTE,
    ST_INT,
    ST_LONG,
    ST_FLOAT,
    ST_DOUBLE,
    ST_HEX,
    ST_STR,
    ST_STRL,
    ST_EMPTY /* a dictionary column without entries: str1, all empty */
} pq_storage;

static struct
{
    pq_file file;
    int open;
    pq_values *cols;
    pq_storage *storage;
    int ncols;
} state;

void cparquet_cleanup(void)
{
    for (int j = 0; state.cols && j < state.ncols; j++)
        pq_values_free(&state.cols[j]);
    free(state.cols);
    free(state.storage);
    if (state.open)
        pq_close(&state.file);
    memset(&state, 0, sizeof(state));
}

static int local_text(const char *name, char *out, int cap)
{
    out[0] = 0;
    return SF_macro_use((char *)name, out, cap - 1);
}

static int save_local(const char *name, const char *value)
{
    return SF_macro_save((char *)name, (char *)value);
}

static int save_number(const char *name, long long v)
{
    char s[32];
    snprintf(s, sizeof(s), "%lld", v);
    return save_local(name, s);
}

static int report(int rc, const char *err)
{
    if (err && *err)
    {
        char line[1100];
        snprintf(line, sizeof(line), "%s\n", err);
        SF_error(line);
    }
    return rc;
}

static int open_file(pq_file *f, char *err, size_t cap)
{
    char path[4096];
    if (local_text("CTOOLS_PQ_file", path, sizeof(path)) || !*path)
        return 198;
    return pq_open(f, path, err, cap);
}

static int is_supported(pq_arrow a)
{
    return a != PA_NULL && a != PA_STRUCT && a != PA_LIST && a != PA_MAP && a != PA_INVALID &&
           a != PA_DURATION;
}

static int describe(int check)
{
    pq_file f;
    char err[1024] = "";
    int rc = open_file(&f, err, sizeof(err));
    if (rc)
    {
        pq_close(&f);
        return rc == 603 ? 603 : report(rc, err);
    }
    /* Stata's describe reads the data; data it cannot read show as -1 rows */
    long long rows = f.num_rows;
    for (int i = 0; check && i < f.nfields && rows >= 0; i++)
    {
        const pq_field *fd = &f.fields[i];
        pq_values v;
        int bad;
        if (fd->leaf >= 0 && is_supported(fd->arrow))
        {
            bad = pq_read_column(&f, fd, 0, -1, 2, &v, err, sizeof(err));
            pq_values_free(&v);
        }
        else
            bad = pq_check_field(&f, fd, err, sizeof(err));
        if (bad == 909)
        {
            pq_close(&f);
            return 909;
        }
        if (bad)
            rows = -1;
    }
    rc = save_number("CTOOLS_PQ_N", rows);
    rc = rc ? rc : save_number("CTOOLS_PQ_k", f.nfields);
    char name[64];
    for (int i = 0; i < f.nfields && !rc; i++)
    {
        snprintf(name, sizeof(name), "CTOOLS_PQ_name_%d", i + 1);
        rc = save_local(name, f.fields[i].name);
        snprintf(name, sizeof(name), "CTOOLS_PQ_type_%d", i + 1);
        rc = rc ? rc : save_local(name, pq_type_name(f.fields[i].arrow));
    }
    pq_close(&f);
    return rc;
}

/* The storage type import parquet gives a column, from the values in range
 * (numbers) or in the whole column (text lengths). */
static pq_storage storage_of(const pq_values *v, int *width)
{
    const double maxfloat = 1.70141173319264430e+38;
    *width = 0;
    switch (v->field->arrow)
    {
    case PA_BOOL:
        return ST_BYTE;
    case PA_INT8:
    case PA_UINT8:
        return v->any && v->maxv > 100 ? ST_INT : ST_BYTE;
    case PA_INT16:
    case PA_UINT16:
        return v->any && v->maxv > 32740 ? ST_LONG : ST_INT;
    case PA_INT32:
    case PA_UINT32:
        return v->any && v->maxv > 2147483620.0 ? ST_DOUBLE : ST_LONG;
    case PA_INT64:
    case PA_UINT64:
        if (v->wide)
        {
            *width = 16;
            return ST_HEX;
        }
        return ST_LONG;
    case PA_HALF:
    case PA_FLOAT:
        return v->any && v->maxv > maxfloat ? ST_DOUBLE : ST_FLOAT;
    case PA_DATE32:
        return ST_LONG;
    case PA_DOUBLE:
    case PA_TIME32:
    case PA_TIME64:
    case PA_TIMESTAMP:
        return ST_DOUBLE;
    case PA_DECIMAL128:
    case PA_DECIMAL256:
        *width = v->field->precision + 2;
        return *width > 2045 ? ST_STRL : ST_STR;
    case PA_DICT:
        /* Codes into the column's dictionary, labeled with its entries */
        if (!v->dn)
        {
            *width = 1;
            return ST_EMPTY;
        }
        return v->dn <= 100 ? ST_BYTE : v->dn <= 32740 ? ST_INT : ST_LONG;
    case PA_UTF8:
    case PA_LARGE_UTF8:
        *width = v->max_len ? (int)(v->max_len > 2046 ? 2046 : v->max_len) : 1;
        return v->max_len > 2045 ? ST_STRL : ST_STR;
    default:
        return ST_STRL;
    }
}

typedef struct
{
    pq_file *file;
    const pq_field *field;
    int64_t first, last;
    pq_values *out;
    int rc;
    char err[512];
} pq_task;

static void *decode_task(void *arg)
{
    pq_task *t = arg;
    t->rc = pq_read_column(t->file, t->field, t->first, t->last, 1, t->out, t->err, sizeof(t->err));
    return NULL;
}

static int scan(void)
{
    cparquet_cleanup();
    char err[1024] = "", text[64];
    int rc = open_file(&state.file, err, sizeof(err));
    if (rc)
    {
        pq_close(&state.file);
        return rc == 603 ? 603 : report(rc, err);
    }
    state.open = 1;
    long long first, last;
    if (local_text("CTOOLS_PQ_first", text, sizeof(text)) || sscanf(text, "%lld", &first) != 1 ||
        local_text("CTOOLS_PQ_last", text, sizeof(text)) || sscanf(text, "%lld", &last) != 1)
        return 198;
    /* Reject the requested SPI row count before allocating or decoding columns. */
    if (first < 1 || last < first || last > state.file.num_rows)
        return 198;
    if ((uint64_t)last - (uint64_t)first + 1 > INT32_MAX)
        return 920;
    size_t cap = (size_t)state.file.nfields * 12 + 64;
    char *list = malloc(cap);
    if (!list)
        return 909;
    if (local_text("CTOOLS_PQ_cols", list, (int)(cap > 0x7fffffff ? 0x7fffffff : cap)))
    {
        free(list);
        return 198;
    }
    int n = 0;
    for (char *p = list; *p;)
    {
        while (*p == ' ')
            p++;
        if (!*p)
            break;
        n++;
        while (*p && *p != ' ')
            p++;
    }
    state.cols = calloc(n ? (size_t)n : 1, sizeof(pq_values));
    state.storage = calloc(n ? (size_t)n : 1, sizeof(pq_storage));
    if (!state.cols || !state.storage)
    {
        free(list);
        return 909;
    }
    pq_task *tasks = calloc(n ? (size_t)n : 1, sizeof(pq_task));
    if (!tasks)
    {
        free(list);
        return 909;
    }
    char *p = list;
    for (int j = 0; j < n && !rc; j++)
    {
        int index = (int)strtol(p, &p, 10);
        if (index < 0 || index >= state.file.nfields ||
            state.file.fields[index].leaf < 0 || !is_supported(state.file.fields[index].arrow))
            rc = 198;
        else
            tasks[j] = (pq_task){&state.file, &state.file.fields[index], first - 1, last - 1,
                                 &state.cols[j], 0, ""};
    }
    state.ncols = rc ? 0 : n;
    /* Columns decode independently, each with its own file handle */
    ctools_persistent_pool *pool = !rc && n > 1 ? ctools_get_global_pool() : NULL;
    if (pool && !ctools_persistent_pool_submit_batch(pool, decode_task, tasks, (size_t)n,
                                                     sizeof(pq_task)))
        ctools_persistent_pool_wait(pool);
    else
        for (int j = 0; j < n && !rc; j++)
            decode_task(&tasks[j]);
    /* The first failing column, in order, reports */
    for (int j = 0; j < n && !rc; j++)
        if (tasks[j].rc)
        {
            rc = tasks[j].rc;
            snprintf(err, sizeof(err), "%s", tasks[j].err);
        }
    free(tasks);
    for (int j = 0; j < n && !rc; j++)
    {
        const pq_field *fd = state.cols[j].field;
        int width;
        pq_storage st = storage_of(&state.cols[j], &width);
        state.storage[j] = st;
        static const char *names[] = {"byte", "int", "long", "float", "double"};
        char type[32], key[64];
        if (st <= ST_DOUBLE)
            snprintf(type, sizeof(type), "%s", names[st]);
        else if (st == ST_STRL)
            snprintf(type, sizeof(type), "strL");
        else if (st == ST_EMPTY)
            snprintf(type, sizeof(type), "str1");
        else
            snprintf(type, sizeof(type), "str%d", width);
        snprintf(key, sizeof(key), "CTOOLS_PQ_vtype_%d", j + 1);
        rc = save_local(key, type);
        /* A labeled dictionary column is displayed as wide as its longest
           entry (in bytes, at most 2045), when that exceeds the default */
        char fmt[32] = "";
        if (fd->arrow == PA_DATE32)
            snprintf(fmt, sizeof(fmt), "%%td");
        else if (fd->arrow == PA_TIMESTAMP)
            snprintf(fmt, sizeof(fmt), "%%tc");
        else if (fd->arrow == PA_DICT && st != ST_EMPTY)
        {
            size_t w = state.cols[j].dmax_len > 2045 ? 2045 : state.cols[j].dmax_len;
            if (w > (st == ST_LONG ? 12u : 8u))
                snprintf(fmt, sizeof(fmt), "%%%zu.0g", w);
        }
        snprintf(key, sizeof(key), "CTOOLS_PQ_fmt_%d", j + 1);
        rc = rc ? rc : save_local(key, fmt);
        snprintf(key, sizeof(key), "CTOOLS_PQ_dict_%d", j + 1);
        rc = rc ? rc : save_number(key, fd->arrow == PA_DICT && st != ST_EMPTY ? (long long)state.cols[j].dn : 0);
        /* strL and text with NUL bytes go through Mata, which stores bytes */
        snprintf(key, sizeof(key), "CTOOLS_PQ_blob_%d", j + 1);
        rc = rc ? rc : save_number(key, st == ST_STRL || (st == ST_STR && state.cols[j].has_nul));
    }
    free(list);
    if (rc)
    {
        cparquet_cleanup();
        return rc == 691 || rc == 3900 ? report(rc, err) : rc;
    }
    return 0;
}

static double stored(pq_storage st, double x)
{
    if (!isfinite(x))
        return SV_missval;
    switch (st)
    {
    case ST_BYTE:
        return x >= -127 && x <= 100 ? x : SV_missval;
    case ST_INT:
        return x >= -32767 && x <= 32740 ? x : SV_missval;
    case ST_LONG:
        return x >= -2147483647.0 && x <= 2147483620.0 ? x : SV_missval;
    case ST_FLOAT:
        return fabs(x) <= 1.70141173319264430e+38 ? x : SV_missval;
    default:
        return x;
    }
}

static int load(void)
{
    if (!state.cols || SF_nvars() != state.ncols)
        return 198;
    size_t nobs = state.ncols ? state.cols[0].n : 0;
    if ((size_t)SF_nobs() < nobs)
        return 198;
    char *buf = NULL;
    size_t bufcap = 0;
    int rc = 0;
    for (int j = 0; j < state.ncols && !rc; j++)
    {
        pq_values *v = &state.cols[j];
        pq_storage st = state.storage[j];
        if (st == ST_STRL || (st == ST_STR && v->has_nul))
            continue;
        int var = j + 1;
        if (st == ST_EMPTY)
            continue; /* created empty */
        if ((st >= ST_HEX) != !!SF_var_is_string(var))
        {
            rc = 109;
            break;
        }
        for (size_t i = 0; i < v->n && !rc; i++)
        {
            if ((i & 65535) == 0 && SF_poll())
            {
                rc = 1;
                break;
            }
            int obs = (int)i + 1;
            if (st == ST_HEX)
            {
                char hex[24] = "";
                if (v->valid[i])
                    snprintf(hex, sizeof(hex), "%llx", (unsigned long long)v->i64[i]);
                rc = SF_sstore(var, obs, hex);
            }
            else if (st == ST_STR)
            {
                size_t len = (size_t)(v->off[i + 1] - v->off[i]);
                if (!v->valid[i] || !len)
                {
                    rc = SF_sstore(var, obs, "");
                    continue;
                }
                if (len + 1 > bufcap)
                {
                    free(buf);
                    bufcap = len + 1;
                    buf = malloc(bufcap);
                    if (!buf)
                    {
                        rc = 909;
                        break;
                    }
                }
                memcpy(buf, v->bytes + v->off[i], len);
                buf[len] = 0;
                rc = SF_sstore(var, obs, buf);
            }
            else
            {
                double x = SV_missval;
                if (v->valid[i])
                    x = v->i64 ? (double)v->i64[i] : v->num[i];
                rc = SF_vstore(var, obs, stored(st, x));
            }
        }
    }
    free(buf);
    return rc;
}

/* A dictionary column's entries, length-prefixed, for its value labels */
static int dictfile(const pq_values *v)
{
    char path[4096];
    if (local_text("CTOOLS_PQ_blobpath", path, sizeof(path)) || !*path)
        return 198;
    FILE *out = ctools_fopen(path, "wb");
    if (!out)
        return 603;
    int rc = 0;
    for (size_t e = 0; e < v->dn && !rc; e++)
    {
        uint64_t len = v->doff[e + 1] - v->doff[e];
        uint8_t head[4] = {(uint8_t)len, (uint8_t)(len >> 8), (uint8_t)(len >> 16),
                           (uint8_t)(len >> 24)};
        if (len > UINT32_MAX || fwrite(head, 1, 4, out) != 4 ||
            (len && fwrite(v->dbytes + v->doff[e], 1, (size_t)len, out) != len))
            rc = 693;
    }
    if (fclose(out) && !rc)
        rc = 693;
    return rc;
}

static int blobfile(int j, int dict)
{
    if (!state.cols || j < 1 || j > state.ncols)
        return 198;
    pq_values *v = &state.cols[j - 1];
    if (dict)
        return dictfile(v);
    char path[4096];
    if (local_text("CTOOLS_PQ_blobpath", path, sizeof(path)) || !*path)
        return 198;
    FILE *out = ctools_fopen(path, "wb");
    if (!out)
        return 603;
    int rc = 0;
    for (size_t i = 0; i < v->n && !rc; i++)
    {
        uint64_t len = v->valid[i] ? v->off[i + 1] - v->off[i] : 0;
        if (len > UINT32_MAX)
        {
            rc = 109;
            break;
        }
        uint8_t head[4] = {(uint8_t)len, (uint8_t)(len >> 8), (uint8_t)(len >> 16),
                           (uint8_t)(len >> 24)};
        if (fwrite(head, 1, 4, out) != 4 || (len && fwrite(v->bytes + v->off[i], 1, (size_t)len, out) != len))
            rc = 693;
    }
    if (fclose(out) && !rc)
        rc = 693;
    return rc;
}

ST_retcode cparquet_main(const char *args)
{
    char op[32] = "";
    int j = 0;
    if (!args || sscanf(args, "%31s %d", op, &j) < 1)
        return 198;
    if (!strcmp(op, "describe"))
        return describe(1);
    if (!strcmp(op, "schema"))
        return describe(0);
    if (!strcmp(op, "scan"))
        return scan();
    if (!strcmp(op, "load"))
        return load();
    if (!strcmp(op, "blobfile") || !strcmp(op, "dictfile"))
        return blobfile(j, op[0] == 'd');
    if (!strcmp(op, "clear"))
    {
        cparquet_cleanup();
        return 0;
    }
    return 198;
}

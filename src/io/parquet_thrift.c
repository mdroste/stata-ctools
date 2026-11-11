/* Parquet footer: Thrift compact protocol decoding of FileMetaData and the
 * mapping of its schema to the value types import parquet presents. */
#include "parquet_internal.h"
#include "ctools_file.h"
#include <limits.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
    const uint8_t *p, *end;
    int bad, depth;
} tr;

enum
{
    T_STOP = 0,
    T_TRUE = 1,
    T_FALSE = 2,
    T_BYTE = 3,
    T_I16 = 4,
    T_I32 = 5,
    T_I64 = 6,
    T_DOUBLE = 7,
    T_BINARY = 8,
    T_LIST = 9,
    T_SET = 10,
    T_MAP = 11,
    T_STRUCT = 12
};

static uint64_t tr_varint(tr *t)
{
    uint64_t v = 0;
    for (int shift = 0; shift < 70 && t->p < t->end; shift += 7)
    {
        uint8_t b = *t->p++;
        if (shift < 64)
            v |= (uint64_t)(b & 0x7f) << shift;
        if (!(b & 0x80))
            return v;
    }
    t->bad = 1;
    return 0;
}

static int64_t tr_zigzag(tr *t)
{
    uint64_t u = tr_varint(t);
    return (int64_t)(u >> 1) ^ -(int64_t)(u & 1);
}

/* Any integer field, whatever width the writer chose */
static int64_t tr_int(tr *t, int type)
{
    if (type == T_BYTE)
    {
        if (t->p >= t->end)
        {
            t->bad = 1;
            return 0;
        }
        return (int8_t)*t->p++;
    }
    if (type == T_I16 || type == T_I32 || type == T_I64)
        return tr_zigzag(t);
    t->bad = 1;
    return 0;
}

static int32_t tr_i32(tr *t, int type)
{
    int64_t v = tr_int(t, type);
    if (v < INT32_MIN || v > INT32_MAX)
        t->bad = 1;
    return (int32_t)v;
}

static const uint8_t *tr_binary(tr *t, size_t *n)
{
    uint64_t len = tr_varint(t);
    if (t->bad || len > (uint64_t)(t->end - t->p))
    {
        t->bad = 1;
        *n = 0;
        return NULL;
    }
    const uint8_t *s = t->p;
    t->p += len;
    *n = (size_t)len;
    return s;
}

/* Field header: returns the type (T_STOP at the end of a struct) */
static int tr_field(tr *t, int *id)
{
    if (t->p >= t->end)
    {
        t->bad = 1;
        return T_STOP;
    }
    uint8_t b = *t->p++;
    if (b == 0)
        return T_STOP;
    int delta = b >> 4;
    if (delta)
        *id += delta;
    else
        *id = (int)tr_zigzag(t);
    return b & 0x0f;
}

static uint64_t tr_list(tr *t, int *etype)
{
    if (t->p >= t->end)
    {
        t->bad = 1;
        return 0;
    }
    uint8_t b = *t->p++;
    *etype = b & 0x0f;
    uint64_t n = b >> 4;
    if (n == 15)
        n = tr_varint(t);
    /* Every element takes at least one byte (an empty struct is its stop) */
    if (n > (uint64_t)(t->end - t->p))
        t->bad = 1;
    return t->bad ? 0 : n;
}

static void tr_skip(tr *t, int type);

static void tr_skip_struct(tr *t)
{
    int id = 0;
    if (++t->depth > 64)
    {
        t->bad = 1;
        return;
    }
    for (;;)
    {
        int type = tr_field(t, &id);
        if (t->bad || type == T_STOP)
            break;
        tr_skip(t, type);
    }
    t->depth--;
}

static void tr_skip(tr *t, int type)
{
    size_t n;
    switch (type)
    {
    case T_TRUE:
    case T_FALSE:
        break;
    case T_BYTE:
        if (t->p >= t->end)
            t->bad = 1;
        else
            t->p++;
        break;
    case T_I16:
    case T_I32:
    case T_I64:
        tr_varint(t);
        break;
    case T_DOUBLE:
        if (t->end - t->p < 8)
            t->bad = 1;
        else
            t->p += 8;
        break;
    case T_BINARY:
        tr_binary(t, &n);
        break;
    case T_LIST:
    case T_SET:
    {
        int et;
        uint64_t count = tr_list(t, &et);
        for (uint64_t i = 0; i < count && !t->bad; i++)
            tr_skip(t, et == T_TRUE || et == T_FALSE ? T_BYTE : et);
        break;
    }
    case T_MAP:
    {
        uint64_t count = tr_varint(t);
        if (count && t->p >= t->end)
            t->bad = 1;
        if (t->bad || count > (uint64_t)(t->end - t->p))
        {
            t->bad = 1;
            break;
        }
        int kv = count ? *t->p++ : 0;
        for (uint64_t i = 0; i < count && !t->bad; i++)
        {
            int k = kv >> 4, v = kv & 15;
            tr_skip(t, k == T_TRUE || k == T_FALSE ? T_BYTE : k);
            tr_skip(t, v == T_TRUE || v == T_FALSE ? T_BYTE : v);
        }
        break;
    }
    case T_STRUCT:
        tr_skip_struct(t);
        break;
    default:
        t->bad = 1;
    }
}

/* LogicalType union: the member id names the type; some carry parameters. */
static void read_logical(tr *t, pq_schema *s)
{
    int id = 0;
    for (;;)
    {
        int type = tr_field(t, &id);
        if (t->bad || type == T_STOP)
            break;
        if (type != T_STRUCT)
        {
            tr_skip(t, type);
            continue;
        }
        s->logical = id;
        int fid = 0;
        for (;;)
        {
            int ft = tr_field(t, &fid);
            if (t->bad || ft == T_STOP)
                break;
            if (id == 5 && fid == 1)
                s->lt_scale = tr_i32(t, ft);
            else if (id == 5 && fid == 2)
                s->lt_precision = tr_i32(t, ft);
            else if ((id == 7 || id == 8) && fid == 2 && ft == T_STRUCT)
            {
                int uid = 0;
                for (;;)
                {
                    int ut = tr_field(t, &uid);
                    if (t->bad || ut == T_STOP)
                        break;
                    if (ut == T_STRUCT)
                        s->lt_unit = uid;
                    tr_skip(t, ut);
                }
            }
            else if (id == 10 && fid == 1)
                s->lt_bits = (int)tr_int(t, ft);
            else if (id == 10 && fid == 2 && (ft == T_TRUE || ft == T_FALSE))
                s->lt_signed = ft == T_TRUE;
            else
                tr_skip(t, ft);
        }
    }
}

static int read_schema(tr *t, pq_schema *s)
{
    int id = 0;
    memset(s, 0, sizeof(*s));
    s->type = -1;
    s->repetition = -1;
    s->num_children = -1;
    s->converted = -1;
    s->lt_signed = 1;
    for (;;)
    {
        int type = tr_field(t, &id);
        if (t->bad || type == T_STOP)
            break;
        size_t n;
        const uint8_t *b;
        switch (id)
        {
        case 1:
            s->type = tr_i32(t, type);
            break;
        case 2:
            s->type_length = tr_i32(t, type);
            break;
        case 3:
            s->repetition = tr_i32(t, type);
            break;
        case 4:
            if (type != T_BINARY)
            {
                t->bad = 1;
                break;
            }
            b = tr_binary(t, &n);
            s->name = malloc(n + 1);
            if (!s->name)
                return 909;
            if (n)
                memcpy(s->name, b, n);
            s->name[n] = 0;
            break;
        case 5:
            s->num_children = tr_i32(t, type);
            break;
        case 6:
            s->converted = tr_i32(t, type);
            break;
        case 7:
            s->scale = tr_i32(t, type);
            break;
        case 8:
            s->precision = tr_i32(t, type);
            break;
        case 10:
            if (type == T_STRUCT)
                read_logical(t, s);
            else
                tr_skip(t, type);
            break;
        default:
            tr_skip(t, type);
        }
    }
    if (!s->name)
    {
        s->name = calloc(1, 1);
        if (!s->name)
            return 909;
    }
    return 0;
}

static void read_column_meta(tr *t, pq_chunk *c)
{
    int id = 0;
    c->has_meta = 1;
    for (;;)
    {
        int type = tr_field(t, &id);
        if (t->bad || type == T_STOP)
            break;
        switch (id)
        {
        case 1:
            c->type = tr_i32(t, type);
            break;
        case 4:
            c->codec = tr_i32(t, type);
            break;
        case 5:
            c->num_values = tr_int(t, type);
            break;
        case 7:
            c->total_compressed = tr_int(t, type);
            break;
        case 9:
            c->data_offset = tr_int(t, type);
            break;
        case 11:
            c->dict_offset = tr_int(t, type);
            break;
        default:
            tr_skip(t, type);
        }
    }
}

static int read_row_group(tr *t, pq_group *g, int nleaves)
{
    int id = 0, seen = 0;
    for (;;)
    {
        int type = tr_field(t, &id);
        if (t->bad || type == T_STOP)
            break;
        if (id == 1 && type == T_LIST)
        {
            int et;
            uint64_t n = tr_list(t, &et);
            if (t->bad || et != T_STRUCT || n != (uint64_t)nleaves || seen)
            {
                t->bad = 1;
                break;
            }
            seen = 1;
            g->chunks = calloc(n ? n : 1, sizeof(pq_chunk));
            if (!g->chunks)
                return 909;
            for (uint64_t i = 0; i < n && !t->bad; i++)
            {
                pq_chunk *c = &g->chunks[i];
                c->dict_offset = -1;
                int cid = 0;
                for (;;)
                {
                    int ct = tr_field(t, &cid);
                    if (t->bad || ct == T_STOP)
                        break;
                    if (cid == 1 && ct == T_BINARY)
                    {
                        size_t len;
                        tr_binary(t, &len);
                        c->has_path = 1;
                    }
                    else if (cid == 3 && ct == T_STRUCT)
                        read_column_meta(t, c);
                    else
                        tr_skip(t, ct);
                }
            }
        }
        else if (id == 3)
            g->num_rows = tr_int(t, type);
        else
            tr_skip(t, type);
    }
    if (!seen)
        t->bad = 1;
    return 0;
}

const char *pq_type_name(pq_arrow a)
{
    static const char *names[] = {
        "bool",       "int8",          "int16",       "int32",      "int64",     "uint8",
        "uint16",     "uint32",        "uint64",      "halffloat",  "float",     "double",
        "date32",     "time32",        "time64",      "timestamp",  "decimal128", "decimal256",
        "utf8",       "binary",        "fixed_size_binary", "null", "struct",    "list",
        "map",        "unknown",       "dictionary",  "large_utf8", "large_binary", "duration"};
    return names[a];
}

/* Parquet's logical (or, failing that, converted) annotation of a primitive,
 * as Arrow's reader presents it. */
static pq_arrow primitive_type(const pq_schema *s, pq_field *f)
{
    int lt = s->logical, ct = s->converted;
    f->precision = lt == 5 ? s->lt_precision : s->precision;
    f->scale = lt == 5 ? s->lt_scale : s->scale;
    pq_arrow decimal = f->precision > 38 ? PA_DECIMAL256 : PA_DECIMAL128;
    if (lt == 11)
        return PA_NULL;
    switch (s->type)
    {
    case PQ_BOOLEAN:
        return lt || ct >= 0 ? PA_INVALID : PA_BOOL;
    case PQ_INT32:
        if (lt == 10)
        {
            if (s->lt_bits == 8)
                return s->lt_signed ? PA_INT8 : PA_UINT8;
            if (s->lt_bits == 16)
                return s->lt_signed ? PA_INT16 : PA_UINT16;
            if (s->lt_bits == 32)
                return s->lt_signed ? PA_INT32 : PA_UINT32;
            return PA_INVALID;
        }
        if (lt == 6)
            return PA_DATE32;
        if (lt == 7)
        {
            f->unit = 1;
            return s->lt_unit == 1 ? PA_TIME32 : PA_INVALID;
        }
        if (lt == 5)
            return decimal;
        if (lt)
            return PA_INVALID;
        switch (ct)
        {
        case -1:
        case 17:
            return PA_INT32;
        case 15:
            return PA_INT8;
        case 16:
            return PA_INT16;
        case 11:
            return PA_UINT8;
        case 12:
            return PA_UINT16;
        case 13:
            return PA_UINT32;
        case 6:
            return PA_DATE32;
        case 7:
            f->unit = 1;
            return PA_TIME32;
        case 5:
            return decimal;
        }
        return PA_INVALID;
    case PQ_INT64:
        if (lt == 10)
            return s->lt_bits == 64 ? (s->lt_signed ? PA_INT64 : PA_UINT64) : PA_INVALID;
        if (lt == 7)
        {
            f->unit = s->lt_unit;
            return s->lt_unit == 2 || s->lt_unit == 3 ? PA_TIME64 : PA_INVALID;
        }
        if (lt == 8)
        {
            f->unit = s->lt_unit;
            return s->lt_unit >= 1 && s->lt_unit <= 3 ? PA_TIMESTAMP : PA_INVALID;
        }
        if (lt == 5)
            return decimal;
        if (lt)
            return PA_INVALID;
        switch (ct)
        {
        case -1:
        case 18:
            return PA_INT64;
        case 14:
            return PA_UINT64;
        case 8:
            f->unit = 2;
            return PA_TIME64;
        case 9:
            f->unit = 1;
            return PA_TIMESTAMP;
        case 10:
            f->unit = 2;
            return PA_TIMESTAMP;
        case 5:
            return decimal;
        }
        return PA_INVALID;
    case PQ_INT96:
        f->unit = 3;
        return lt ? PA_INVALID : PA_TIMESTAMP;
    case PQ_FLOAT:
        return lt ? PA_INVALID : PA_FLOAT;
    case PQ_DOUBLE:
        return lt ? PA_INVALID : PA_DOUBLE;
    case PQ_BYTE_ARRAY:
        if (lt == 1 || lt == 12)
            return PA_UTF8;
        if (lt == 4 || lt == 13)
            return PA_BINARY;
        if (lt == 5)
            return decimal;
        if (lt)
            return PA_INVALID;
        if (ct == 0 || ct == 19)
            return PA_UTF8;
        if (ct == 4 || ct == 20 || ct == -1)
            return PA_BINARY;
        if (ct == 5)
            return decimal;
        return PA_INVALID;
    case PQ_FLBA:
        if (s->type_length <= 0)
            return PA_INVALID;
        if (lt == 5)
            return decimal;
        if (lt == 15)
            return s->type_length == 2 ? PA_HALF : PA_INVALID;
        if (lt == 14)
            return s->type_length == 16 ? PA_FLBA : PA_INVALID;
        if (lt)
            return PA_INVALID;
        if (ct == 5)
            return decimal;
        if (ct == -1 || ct == 21)
            return PA_FLBA;
        return PA_INVALID;
    }
    return PA_INVALID;
}

/* Walks the subtree at schema[*pos]; returns 0 or nonzero when malformed */
static int walk(pq_file *f, int *pos, int depth)
{
    if (*pos >= f->nschema || depth > 100)
        return 1;
    pq_schema *s = &f->schema[(*pos)++];
    if (s->num_children > 0)
    {
        for (int c = 0; c < s->num_children; c++)
            if (walk(f, pos, depth + 1))
                return 1;
        return 0;
    }
    if (s->type < 0 || s->type > PQ_FLBA)
        return 1;
    f->nleaves++;
    return 0;
}

static int build_fields(pq_file *f)
{
    if (f->nschema < 1 || f->schema[0].num_children < 0)
        return 1;
    int n = f->schema[0].num_children;
    f->fields = calloc(n ? (size_t)n : 1, sizeof(pq_field));
    if (!f->fields)
        return 909;
    f->nfields = n;
    int pos = 1;
    for (int i = 0; i < n; i++)
    {
        pq_field *fd = &f->fields[i];
        int start = pos, leaf = f->nleaves;
        if (walk(f, &pos, 1))
            return 1;
        fd->first_leaf = leaf;
        fd->nleaf = f->nleaves - leaf;
        pq_schema *s = &f->schema[start];
        fd->name = s->name;
        fd->element = start;
        fd->physical = s->type;
        fd->type_length = s->type_length;
        fd->leaf = -1;
        fd->origin_type = -1;
        if (s->num_children > 0)
        {
            if (s->logical == 3 || s->converted == 3 || s->repetition == 2)
                fd->arrow = PA_LIST;
            else if (s->logical == 2 || s->converted == 1 || s->converted == 2)
                fd->arrow = PA_MAP;
            else
                fd->arrow = PA_STRUCT;
        }
        else if (s->repetition == 2)
            fd->arrow = PA_LIST;
        else
        {
            fd->arrow = primitive_type(s, fd);
            if (fd->arrow == PA_INVALID)
                return 1;
            fd->leaf = leaf;
            fd->max_def = s->repetition == 1;
        }
    }
    return pos == f->nschema ? 0 : 1;
}

static void read_data_header(tr *t, pq_page *h, int v2)
{
    int id = 0;
    for (;;)
    {
        int type = tr_field(t, &id);
        if (t->bad || type == T_STOP)
            break;
        if (type == T_TRUE || type == T_FALSE)
        {
            /* v2 is_compressed; a dictionary page's is_sorted */
            if (v2 && id == 7)
                h->compressed = type == T_TRUE;
            continue;
        }
        if (id == 1)
            h->num_values = tr_i32(t, type);
        else if (!v2 && id == 2)
            h->encoding = tr_i32(t, type);
        else if (!v2 && id == 3)
            h->def_enc = tr_i32(t, type);
        else if (!v2 && id == 4)
            h->rep_enc = tr_i32(t, type);
        else if (v2 && id == 2)
            h->num_nulls = tr_i32(t, type);
        else if (v2 && id == 3)
            h->num_rows = tr_i32(t, type);
        else if (v2 && id == 4)
            h->encoding = tr_i32(t, type);
        else if (v2 && id == 5)
            h->def_len = tr_i32(t, type);
        else if (v2 && id == 6)
            h->rep_len = tr_i32(t, type);
        else
            tr_skip(t, type);
    }
}

int pq_page_header(const uint8_t *p, size_t n, pq_page *h, size_t *used)
{
    tr t = {p, p + n, 0, 0};
    int id = 0, seen = 0;
    memset(h, 0, sizeof(*h));
    h->type = -1;
    h->compressed = 1;
    h->def_enc = PQ_RLE;
    for (;;)
    {
        int type = tr_field(&t, &id);
        if (t.bad || type == T_STOP)
            break;
        if (id == 1)
            h->type = tr_i32(&t, type);
        else if (id == 2)
            h->usize = tr_i32(&t, type);
        else if (id == 3)
            h->csize = tr_i32(&t, type);
        else if ((id == 5 || id == 7 || id == 8) && type == T_STRUCT)
        {
            read_data_header(&t, h, id == 8);
            seen = 1;
        }
        else
            tr_skip(&t, type);
    }
    if (t.bad || h->type < 0 || h->usize < 0 || h->csize < 0 || h->num_values < 0 ||
        h->def_len < 0 || h->rep_len < 0 || (h->type != 1 && !seen))
        return 1;
    *used = (size_t)(t.p - p);
    return 0;
}

void pq_close(pq_file *f)
{
    if (f->fp)
        fclose(f->fp);
    free(f->path);
    for (int i = 0; f->schema && i < f->nschema; i++)
        free(f->schema[i].name);
    for (int g = 0; f->groups && g < f->ngroups; g++)
        free(f->groups[g].chunks);
    free(f->schema);
    free(f->groups);
    free(f->fields);
    free(f->footer);
    memset(f, 0, sizeof(*f));
}

/* ---- ARROW:schema --------------------------------------------------------
 * Files written by Arrow carry their Arrow schema, base64-encoded as an IPC
 * Schema message (a flatbuffer). Arrow's reader uses it to restore some types
 * the Parquet annotations cannot express; import parquet sees those types,
 * so they are read here too. Only what changes the result is used: a
 * dictionary (categorical) field over string or binary data, the large
 * string and binary types, and durations stored as INT64. A value that is
 * not a readable schema fails as it does in Arrow. */
static int base64_value(int c)
{
    if (c >= 'A' && c <= 'Z')
        return c - 'A';
    if (c >= 'a' && c <= 'z')
        return c - 'a' + 26;
    if (c >= '0' && c <= '9')
        return c - '0' + 52;
    if (c == '+')
        return 62;
    if (c == '/')
        return 63;
    return -1;
}

static uint8_t *base64_decode(const uint8_t *s, size_t n, size_t *out)
{
    while (n && (s[n - 1] == '=' || s[n - 1] == '\n' || s[n - 1] == '\r'))
        n--;
    if (n % 4 == 1)
        return NULL;
    uint8_t *d = malloc(n / 4 * 3 + 3);
    if (!d)
        return NULL;
    size_t k = 0;
    uint32_t acc = 0;
    int bits = 0;
    for (size_t i = 0; i < n; i++)
    {
        int v = base64_value(s[i]);
        if (v < 0)
        {
            free(d);
            return NULL;
        }
        acc = (acc << 6) | (uint32_t)v;
        bits += 6;
        if (bits >= 8)
        {
            bits -= 8;
            d[k++] = (uint8_t)(acc >> bits);
        }
    }
    *out = k;
    return d;
}

typedef struct
{
    const uint8_t *b;
    size_t n;
} fbuf;

static int fb_u32(const fbuf *f, size_t pos, uint32_t *v)
{
    if (pos > f->n || f->n - pos < 4)
        return 1;
    *v = (uint32_t)f->b[pos] | (uint32_t)f->b[pos + 1] << 8 | (uint32_t)f->b[pos + 2] << 16 |
         (uint32_t)f->b[pos + 3] << 24;
    return 0;
}

static int fb_u16(const fbuf *f, size_t pos, uint16_t *v)
{
    if (pos > f->n || f->n - pos < 2)
        return 1;
    *v = (uint16_t)(f->b[pos] | f->b[pos + 1] << 8);
    return 0;
}

/* Position of a table's field (0 when absent); -1 when malformed */
static int64_t fb_field(const fbuf *f, size_t table, int slot)
{
    uint32_t so;
    uint16_t vsize, tsize, off;
    if (fb_u32(f, table, &so))
        return -1;
    int64_t vt = (int64_t)table - (int32_t)so;
    if (vt < 0 || fb_u16(f, (size_t)vt, &vsize) || fb_u16(f, (size_t)vt + 2, &tsize) || vsize < 4)
        return -1;
    if ((size_t)(4 + 2 * slot) + 2 > vsize)
        return 0;
    if (fb_u16(f, (size_t)vt + 4 + 2 * (size_t)slot, &off))
        return -1;
    if (!off)
        return 0;
    if (off >= tsize || table + off >= f->n)
        return -1;
    return (int64_t)(table + off);
}

/* The object an offset field refers to */
static int64_t fb_ref(const fbuf *f, int64_t pos)
{
    uint32_t rel;
    if (pos <= 0 || fb_u32(f, (size_t)pos, &rel))
        return -1;
    uint64_t at = (uint64_t)pos + rel;
    return at < f->n ? (int64_t)at : -1;
}

static int apply_arrow_schema(pq_file *pf, const uint8_t *value, size_t len)
{
    size_t n;
    uint8_t *raw = base64_decode(value, len, &n);
    if (!raw)
        return 1;
    fbuf f = {raw, n};
    size_t start = 0;
    uint32_t word;
    int rc = 1;
    /* An encapsulated message: continuation marker and length, or (before
       Arrow 0.15) the length alone */
    if (fb_u32(&f, 0, &word))
        goto done;
    if (word == 0xffffffffu)
        start = 8;
    else
        start = 4;
    if (fb_u32(&f, start - 4, &word) || word > n - start)
        goto done;
    f.b = raw + start;
    f.n = word;
    uint32_t root;
    if (fb_u32(&f, 0, &root) || root >= f.n)
        goto done;
    int64_t ht = fb_field(&f, root, 1), hd = fb_field(&f, root, 2);
    if (ht <= 0 || hd <= 0 || f.b[ht] != 1)
        goto done;
    int64_t schema = fb_ref(&f, hd);
    if (schema < 0)
        goto done;
    int64_t fields = fb_field(&f, (size_t)schema, 1);
    if (fields < 0)
        goto done;
    uint32_t count = 0;
    int64_t vec = fields ? fb_ref(&f, fields) : 0;
    if (fields && (vec < 0 || fb_u32(&f, (size_t)vec, &count) || count > (f.n - (size_t)vec) / 4))
        goto done;
    uint64_t dict_ids[64];
    int dict_types[64];
    uint32_t ndict = 0;
    for (uint32_t i = 0; i < count; i++)
    {
        int64_t field = fb_ref(&f, vec + 4 + 4 * (int64_t)i);
        if (field < 0)
            goto done;
        int64_t tt = fb_field(&f, (size_t)field, 2), dict = fb_field(&f, (size_t)field, 4);
        if (tt < 0 || dict < 0)
            goto done;
        if (dict > 0)
        {
            /* Arrow refuses one dictionary id for two different value types */
            int64_t enc = fb_ref(&f, dict), idf;
            uint64_t id = 0;
            if (enc < 0 || (idf = fb_field(&f, (size_t)enc, 0)) < 0)
                goto done;
            if (idf > 0)
            {
                uint32_t lo, hi;
                if (fb_u32(&f, (size_t)idf, &lo) || fb_u32(&f, (size_t)idf + 4, &hi))
                    goto done;
                id = (uint64_t)hi << 32 | lo;
            }
            int type = tt ? f.b[tt] : 0;
            for (uint32_t k = 0; k < ndict; k++)
                if (dict_ids[k] == id && dict_types[k] != type)
                    goto done;
            if (ndict < 64)
            {
                dict_ids[ndict] = id;
                dict_types[ndict++] = type;
            }
        }
        /* Arrow applies the stored types only when the field counts agree */
        if (count == (uint32_t)pf->nfields)
        {
            pf->fields[i].origin_type = tt ? f.b[tt] : 0;
            pf->fields[i].origin_dict = dict > 0;
        }
    }
    for (int i = 0; count == (uint32_t)pf->nfields && i < pf->nfields; i++)
    {
        pq_field *fd = &pf->fields[i];
        int binary_like = fd->arrow == PA_UTF8 || fd->arrow == PA_BINARY;
        if (fd->origin_dict && binary_like)
            fd->arrow = PA_DICT;
        else if (fd->origin_type == 20 && fd->arrow == PA_UTF8)
            fd->arrow = PA_LARGE_UTF8;
        else if (fd->origin_type == 19 && fd->arrow == PA_BINARY)
            fd->arrow = PA_LARGE_BINARY;
        else if (fd->origin_type == 18 && fd->arrow == PA_INT64)
            fd->arrow = PA_DURATION;
    }
    rc = 0;
done:
    free(raw);
    return rc;
}

static int parse_footer(pq_file *f, const uint8_t *buf, size_t len)
{
    tr t = {buf, buf + len, 0, 0};
    int id = 0, rc = 0, have_schema = 0, have_groups = 0;
    const uint8_t *arrow = NULL;
    size_t arrow_len = 0;
    for (;;)
    {
        int type = tr_field(&t, &id);
        if (t.bad || type == T_STOP)
            break;
        if (id == 2 && type == T_LIST && !have_schema)
        {
            int et;
            uint64_t n = tr_list(&t, &et);
            if (t.bad || et != T_STRUCT || n == 0 || n > INT_MAX / 2)
                return 1;
            f->schema = calloc(n, sizeof(pq_schema));
            if (!f->schema)
                return 909;
            for (uint64_t i = 0; i < n && !t.bad; i++)
            {
                if ((rc = read_schema(&t, &f->schema[i])))
                    return rc;
                f->nschema = (int)i + 1;
            }
            have_schema = 1;
            if (t.bad)
                return 1;
            if ((rc = build_fields(f)))
                return rc;
        }
        else if (id == 3)
            f->num_rows = tr_int(&t, type);
        else if (id == 5 && type == T_LIST)
        {
            /* key_value_metadata: only ARROW:schema matters */
            int et;
            uint64_t n = tr_list(&t, &et);
            if (t.bad || (n && et != T_STRUCT))
                return 1;
            for (uint64_t i = 0; i < n && !t.bad; i++)
            {
                int kid = 0;
                const uint8_t *key = NULL, *value = NULL;
                size_t klen = 0, vlen = 0;
                for (;;)
                {
                    int kt = tr_field(&t, &kid);
                    if (t.bad || kt == T_STOP)
                        break;
                    if (kid == 1 && kt == T_BINARY)
                        key = tr_binary(&t, &klen);
                    else if (kid == 2 && kt == T_BINARY)
                        value = tr_binary(&t, &vlen);
                    else
                        tr_skip(&t, kt);
                }
                if (key && klen == 12 && !memcmp(key, "ARROW:schema", 12) && value)
                {
                    arrow = value;
                    arrow_len = vlen;
                }
            }
        }
        else if (id == 4 && type == T_LIST && have_schema && !have_groups)
        {
            int et;
            uint64_t n = tr_list(&t, &et);
            if (t.bad || (n && et != T_STRUCT) || n > INT_MAX / 2)
                return 1;
            f->groups = calloc(n ? n : 1, sizeof(pq_group));
            if (!f->groups)
                return 909;
            f->ngroups = (int)n;
            for (uint64_t i = 0; i < n && !t.bad; i++)
                if ((rc = read_row_group(&t, &f->groups[i], f->nleaves)))
                    return rc;
            have_groups = 1;
        }
        else
            tr_skip(&t, type);
    }
    if (t.bad || !have_schema || f->num_rows < 0)
        return 1;
    if (arrow)
        return apply_arrow_schema(f, arrow, arrow_len);
    return 0;
}

int pq_open(pq_file *f, const char *path, char *err, size_t errcap)
{
    memset(f, 0, sizeof(*f));
    err[0] = 0;
    f->path = malloc(strlen(path) + 1);
    if (!f->path)
        return 909;
    memcpy(f->path, path, strlen(path) + 1);
    f->fp = ctools_fopen(path, "rb");
    if (!f->fp)
        return 603;
    uint8_t head[4], tail[8];
    if (fread(head, 1, sizeof(head), f->fp) != sizeof(head) || memcmp(head, "PAR1", 4) ||
        ctools_fseek(f->fp, 0, SEEK_END) || (f->size = ctools_ftell(f->fp)) < 12 ||
        ctools_fseek(f->fp, f->size - 8, SEEK_SET) || fread(tail, 1, 8, f->fp) != 8 ||
        memcmp(tail + 4, "PAR1", 4))
    {
        snprintf(err, errcap, "failed to get parquet schema");
        return 3900;
    }
    uint32_t len = (uint32_t)tail[0] | (uint32_t)tail[1] << 8 | (uint32_t)tail[2] << 16 |
                   (uint32_t)tail[3] << 24;
    if (len == 0 || (int64_t)len > f->size - 12)
    {
        snprintf(err, errcap, "failed to get parquet schema");
        return 3900;
    }
    f->footer = malloc(len);
    if (!f->footer)
        return 909;
    if (ctools_fseek(f->fp, f->size - 8 - len, SEEK_SET) || fread(f->footer, 1, len, f->fp) != len)
    {
        snprintf(err, errcap, "failed to get parquet schema");
        return 3900;
    }
    int rc = parse_footer(f, f->footer, len);
    if (rc == 909)
        return 909;
    if (rc)
    {
        snprintf(err, errcap, "failed to get parquet schema");
        return 3900;
    }
    int64_t rows = 0;
    for (int g = 0; g < f->ngroups; g++)
    {
        if (f->groups[g].num_rows < 0 || f->groups[g].num_rows > INT64_MAX - rows)
        {
            snprintf(err, errcap, "failed to get parquet schema");
            return 3900;
        }
        rows += f->groups[g].num_rows;
    }
    /* The footer's row count must agree with its row groups */
    if (rows != f->num_rows)
    {
        snprintf(err, errcap, "failed to get parquet schema");
        return 3900;
    }
    return 0;
}

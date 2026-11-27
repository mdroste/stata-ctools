/*
 * csql_result.c - the rows of a query, to a .dta file or to the screen.
 *
 * The rows arrive as DuckDB chunks (up to 2,048 rows each). They are kept
 * as they are: one pass over them, split across threads, measures what each
 * column needs (integer range, longest string, NUL bytes, ENUM dictionary),
 * which decides the Stata storage type as compress would; then ctools' .dta
 * writer pulls row groups of chunks through its producer interface on its
 * worker threads and the ado loads the file with use. Numeric NULL, NaN and
 * infinities become missing (.); string NULL becomes "".
 *
 * Storage: BOOLEAN byte; integers the narrowest of byte/int/long that holds
 * every value, else double; FLOAT float; DOUBLE and DECIMAL double; DATE long
 * %td; TIMESTAMP (any precision or zone, as UTC) double %tc; TIME double
 * %tcHH:MM:SS; VARCHAR str# up to 2,045 bytes, else strL; BLOB strL; ENUM
 * integer codes 1, 2, ... with a value label of the same name (name_1, ...
 * when an inherited label has that name). FLOAT past Stata's float range
 * (2^127) is missing, as DOUBLE past 2^1023 is.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "csql/csql_internal.h"
#include "formats/ctools_dta.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "ctools_cancel.h"

/* ------------------------------------------------------------ columns */

typedef enum
{
    K_BOOL = 1,
    K_INT,     /* signed integers up to 64 bits */
    K_UINT,    /* unsigned up to 64 bits */
    K_HUGE,    /* HUGEINT */
    K_UHUGE,   /* UHUGEINT */
    K_FLOAT,
    K_DOUBLE,
    K_DECIMAL,
    K_DATE,
    K_TS,      /* timestamps: micros (scale to ms via div) */
    K_TIME,
    K_STRING,
    K_BLOB,
    K_ENUM,
    K_NULL,
    K_OTHER
} col_kind;

typedef struct
{
    char name[CSQL_NAME_BYTES];
    char *orig;      /* DuckDB's column name */
    int renamed;
    duckdb_type type;
    col_kind kind;
    int bytes;       /* K_INT/K_UINT: 1, 2, 4, 8 */
    double div;      /* K_TS: units per millisecond; K_TIME: units per ms; K_DECIMAL: 10^scale */
    duckdb_type internal; /* K_DECIMAL, K_ENUM storage */
    int scale;
    /* stats */
    int64_t imin, imax; /* exact for 64-bit kinds */
    double dmin, dmax;  /* K_HUGE, K_UHUGE, K_UINT above INT64_MAX */
    int any;
    size_t maxlen;
    int has_nul;
    uint64_t imprecise; /* integers beyond 2^53 stored as double */
    uint64_t nonfinite; /* NaN, infinities, out of Stata's range */
    /* enum */
    size_t ndict;
    char **dict;
    /* decision */
    ctools_dta_type stype;
    int width;
    char format[CSQL_FORMAT_BYTES];
    char vlabel[CSQL_NAME_BYTES];
    char label[CSQL_LABEL_BYTES];
} res_col;

/* The rows: chunks held here (fetched from a streaming result), or, with
 * random set, read by index from the materialized result whenever they are
 * needed, so that DuckDB's copy is the only one. */
typedef struct
{
    duckdb_data_chunk *chunks;
    uint64_t *starts; /* nchunks + 1 */
    size_t nchunks, cap;
    uint64_t nrows;
    int random;
    duckdb_result *res;
} res_rows;

static col_kind kind_of(duckdb_logical_type lt, res_col *c)
{
    duckdb_type t = DK(get_type_id)(lt);
    c->type = t;
    c->div = 1;
    switch (t) {
    case DUCKDB_TYPE_BOOLEAN: return K_BOOL;
    case DUCKDB_TYPE_TINYINT: c->bytes = 1; return K_INT;
    case DUCKDB_TYPE_SMALLINT: c->bytes = 2; return K_INT;
    case DUCKDB_TYPE_INTEGER: c->bytes = 4; return K_INT;
    case DUCKDB_TYPE_BIGINT: c->bytes = 8; return K_INT;
    case DUCKDB_TYPE_UTINYINT: c->bytes = 1; return K_UINT;
    case DUCKDB_TYPE_USMALLINT: c->bytes = 2; return K_UINT;
    case DUCKDB_TYPE_UINTEGER: c->bytes = 4; return K_UINT;
    case DUCKDB_TYPE_UBIGINT: c->bytes = 8; return K_UINT;
    case DUCKDB_TYPE_HUGEINT: return K_HUGE;
    case DUCKDB_TYPE_UHUGEINT: return K_UHUGE;
    case DUCKDB_TYPE_FLOAT: return K_FLOAT;
    case DUCKDB_TYPE_DOUBLE: return K_DOUBLE;
    case DUCKDB_TYPE_DECIMAL:
        c->scale = DK(decimal_scale)(lt);
        c->internal = DK(decimal_internal_type)(lt);
        c->div = pow(10.0, c->scale);
        return K_DECIMAL;
    case DUCKDB_TYPE_DATE: return K_DATE;
    case DUCKDB_TYPE_TIMESTAMP:
    case DUCKDB_TYPE_TIMESTAMP_TZ: c->div = 1000; return K_TS;
    case DUCKDB_TYPE_TIMESTAMP_S: c->div = 0.001; return K_TS;
    case DUCKDB_TYPE_TIMESTAMP_MS: c->div = 1; return K_TS;
    case DUCKDB_TYPE_TIMESTAMP_NS:
    case DUCKDB_TYPE_TIMESTAMP_TZ_NS: c->div = 1000000; return K_TS;
    case DUCKDB_TYPE_TIME: c->div = 1000; return K_TIME;
    case DUCKDB_TYPE_TIME_NS: c->div = 1000000; return K_TIME;
    case DUCKDB_TYPE_VARCHAR: return K_STRING;
    case DUCKDB_TYPE_BLOB: return K_BLOB;
    case DUCKDB_TYPE_ENUM: {
        c->internal = DK(enum_internal_type)(lt);
        c->ndict = DK(enum_dictionary_size)(lt);
        c->dict = calloc(c->ndict ? c->ndict : 1, sizeof(char *));
        for (size_t k = 0; c->dict && k < c->ndict; k++) c->dict[k] = DK(enum_dictionary_value)(lt, (idx_t)k);
        return K_ENUM;
    }
    case DUCKDB_TYPE_SQLNULL: return K_NULL;
    default: return K_OTHER;
    }
}

/* Types csql stores directly; the others are cast to VARCHAR in SQL. */
int csql_type_supported(duckdb_type t)
{
    switch (t) {
    case DUCKDB_TYPE_BOOLEAN: case DUCKDB_TYPE_TINYINT: case DUCKDB_TYPE_SMALLINT:
    case DUCKDB_TYPE_INTEGER: case DUCKDB_TYPE_BIGINT: case DUCKDB_TYPE_UTINYINT:
    case DUCKDB_TYPE_USMALLINT: case DUCKDB_TYPE_UINTEGER: case DUCKDB_TYPE_UBIGINT:
    case DUCKDB_TYPE_HUGEINT: case DUCKDB_TYPE_UHUGEINT: case DUCKDB_TYPE_FLOAT:
    case DUCKDB_TYPE_DOUBLE: case DUCKDB_TYPE_DECIMAL: case DUCKDB_TYPE_DATE:
    case DUCKDB_TYPE_TIMESTAMP: case DUCKDB_TYPE_TIMESTAMP_TZ: case DUCKDB_TYPE_TIMESTAMP_S:
    case DUCKDB_TYPE_TIMESTAMP_MS: case DUCKDB_TYPE_TIMESTAMP_NS: case DUCKDB_TYPE_TIMESTAMP_TZ_NS:
    case DUCKDB_TYPE_TIME: case DUCKDB_TYPE_TIME_NS: case DUCKDB_TYPE_VARCHAR: case DUCKDB_TYPE_BLOB:
    case DUCKDB_TYPE_ENUM: case DUCKDB_TYPE_SQLNULL:
        return 1;
    default: return 0;
    }
}

static inline int row_valid(const uint64_t *valid, size_t i)
{
    return !valid || ((valid[i >> 6] >> (i & 63)) & 1);
}

/* The 128-bit unsigned hi * 2^64 + lo, rounded once to the nearest double
 * (ties to even): the top 64 bits, with every lower bit folded into the
 * last one as a sticky bit, round as the whole number would. */
static double u128_to_double(uint64_t hi, uint64_t lo)
{
    if (!hi) return (double)lo;
    int lz = __builtin_clzll(hi);
    int shift = 64 - lz; /* bits below the top 64 */
    uint64_t top, sticky;
    if (lz == 0) {
        top = hi;
        sticky = lo != 0;
    } else {
        top = (hi << lz) | (lo >> shift);
        sticky = (lo & ((UINT64_C(1) << shift) - 1)) != 0;
    }
    return ldexp((double)(top | sticky), shift);
}

/* HUGEINT is two's complement: a negative value converts as its magnitude. */
static inline double huge_to_double(duckdb_hugeint h)
{
    if (h.upper >= 0) return u128_to_double((uint64_t)h.upper, h.lower);
    uint64_t lo = ~h.lower + 1, hi = ~(uint64_t)h.upper + (lo == 0);
    return -u128_to_double(hi, lo);
}

static inline double uhuge_to_double(duckdb_uhugeint h)
{
    return u128_to_double(h.upper, h.lower);
}

/* Signed integer value i of a vector of width bytes. */
static inline int64_t int_at(const void *d, int bytes, size_t i)
{
    switch (bytes) {
    case 1: return ((const int8_t *)d)[i];
    case 2: return ((const int16_t *)d)[i];
    case 4: return ((const int32_t *)d)[i];
    default: return ((const int64_t *)d)[i];
    }
}

static inline uint64_t uint_at(const void *d, int bytes, size_t i)
{
    switch (bytes) {
    case 1: return ((const uint8_t *)d)[i];
    case 2: return ((const uint16_t *)d)[i];
    case 4: return ((const uint32_t *)d)[i];
    default: return ((const uint64_t *)d)[i];
    }
}

static inline uint32_t enum_at(const void *d, duckdb_type internal, size_t i)
{
    switch (internal) {
    case DUCKDB_TYPE_UTINYINT: return ((const uint8_t *)d)[i];
    case DUCKDB_TYPE_USMALLINT: return ((const uint16_t *)d)[i];
    default: return ((const uint32_t *)d)[i];
    }
}

static inline double decimal_at(const void *d, duckdb_type internal, size_t i, double div)
{
    switch (internal) {
    case DUCKDB_TYPE_SMALLINT: return ((const int16_t *)d)[i] / div;
    case DUCKDB_TYPE_INTEGER: return ((const int32_t *)d)[i] / div;
    case DUCKDB_TYPE_BIGINT: return ((const int64_t *)d)[i] / div;
    default: return huge_to_double(((const duckdb_hugeint *)d)[i]) / div;
    }
}

static inline void string_at(const void *d, size_t i, const char **p, size_t *len)
{
    const duckdb_string_t *s = (const duckdb_string_t *)d + i;
    uint32_t n = s->value.inlined.length;
    *len = n;
    *p = n <= 12 ? s->value.inlined.inlined : s->value.pointer.ptr;
}

/* The value of row i as a Stata double (numeric kinds). */
static double value_at(const res_col *c, const void *d, size_t i)
{
    switch (c->kind) {
    case K_BOOL: return ((const bool *)d)[i] ? 1.0 : 0.0;
    case K_INT: return (double)int_at(d, c->bytes, i);
    case K_UINT: return (double)uint_at(d, c->bytes, i);
    case K_HUGE: return huge_to_double(((const duckdb_hugeint *)d)[i]);
    case K_UHUGE: return uhuge_to_double(((const duckdb_uhugeint *)d)[i]);
    case K_FLOAT: return (double)((const float *)d)[i];
    case K_DOUBLE: return ((const double *)d)[i];
    case K_DECIMAL: return decimal_at(d, c->internal, i, c->div);
    case K_DATE: {
        int32_t days = ((const duckdb_date *)d)[i].days;
        /* +-infinity dates */
        if (days == INT32_MAX || days == -INT32_MAX) return csql_missing();
        return (double)days + CSQL_DAYS_1960_1970;
    }
    case K_TS: {
        int64_t v = ((const int64_t *)d)[i];
        if (v == INT64_MAX || v == -INT64_MAX) return csql_missing();
        double ms = c->div == 1 ? (double)v : c->div < 1 ? (double)v * 1000.0
                                                         : (double)(v / (int64_t)c->div) +
                                                               (double)(v % (int64_t)c->div) / c->div;
        return ms + (double)CSQL_MS_1960_1970;
    }
    case K_TIME: {
        int64_t v = ((const int64_t *)d)[i];
        return (double)(v / (int64_t)c->div) + (double)(v % (int64_t)c->div) / c->div;
    }
    case K_ENUM: return (double)enum_at(d, c->internal, i) + 1.0;
    default: return csql_missing();
    }
}

/* ------------------------------------------------------------ stats */

typedef struct
{
    int64_t imin, imax;
    double dmin, dmax;
    int any;
    size_t maxlen;
    int has_nul;
    uint64_t imprecise, nonfinite;
} col_stats;

static void stats_merge(col_stats *a, const col_stats *b)
{
    if (b->any) {
        if (!a->any) {
            a->imin = b->imin, a->imax = b->imax, a->dmin = b->dmin, a->dmax = b->dmax;
        } else {
            if (b->imin < a->imin) a->imin = b->imin;
            if (b->imax > a->imax) a->imax = b->imax;
            if (b->dmin < a->dmin) a->dmin = b->dmin;
            if (b->dmax > a->dmax) a->dmax = b->dmax;
        }
        a->any = 1;
    }
    if (b->maxlen > a->maxlen) a->maxlen = b->maxlen;
    a->has_nul |= b->has_nul;
    a->imprecise += b->imprecise;
    a->nonfinite += b->nonfinite;
}

static void stats_chunk(const res_col *c, duckdb_vector vec, size_t n, col_stats *s)
{
    const void *d = DK(vector_get_data)(vec);
    const uint64_t *valid = DK(vector_get_validity)(vec);
    for (size_t i = 0; i < n; i++) {
        if (!row_valid(valid, i)) continue;
        switch (c->kind) {
        case K_BOOL: case K_ENUM: case K_DATE: break;
        case K_INT: {
            int64_t v = int_at(d, c->bytes, i);
            if (!s->any) s->imin = s->imax = v, s->any = 1;
            else {
                if (v < s->imin) s->imin = v;
                if (v > s->imax) s->imax = v;
            }
            if (v > (INT64_C(1) << 53) || v < -(INT64_C(1) << 53)) s->imprecise++;
            break;
        }
        case K_UINT: {
            uint64_t u = uint_at(d, c->bytes, i);
            int64_t v = u > (uint64_t)INT64_MAX ? INT64_MAX : (int64_t)u;
            if (!s->any) s->imin = s->imax = v, s->any = 1;
            else {
                if (v < s->imin) s->imin = v;
                if (v > s->imax) s->imax = v;
            }
            if (u > (UINT64_C(1) << 53)) s->imprecise++;
            break;
        }
        case K_HUGE: case K_UHUGE: case K_DECIMAL: {
            double v = value_at(c, d, i);
            if (c->kind == K_DECIMAL && c->scale) break;
            int64_t iv = v >= 9.2e18 ? INT64_MAX : v <= -9.2e18 ? INT64_MIN : (int64_t)v;
            if (!s->any) s->imin = s->imax = iv, s->any = 1;
            else {
                if (iv < s->imin) s->imin = iv;
                if (iv > s->imax) s->imax = iv;
            }
            if (fabs(v) > 9007199254740992.0) s->imprecise++;
            break;
        }
        case K_FLOAT: case K_DOUBLE: {
            /* Stata's float range ends below 2^127, its double range below
             * 2^1023; past it (and NaN) a value is missing */
            double v = value_at(c, d, i), lim = c->kind == K_FLOAT ? 0x1p127 : 0x1p1023;
            if (!(v < lim && v > -lim)) s->nonfinite++;
            break;
        }
        case K_STRING: case K_BLOB: {
            const char *p;
            size_t len;
            string_at(d, i, &p, &len);
            if (len > s->maxlen) s->maxlen = len;
            if (!s->has_nul && len && memchr(p, 0, len)) s->has_nul = 1;
            break;
        }
        default: break;
        }
    }
}

/* ------------------------------------------------------------ decisions */

static void decide(res_col *c, int strl_min)
{
    c->format[0] = '\0';
    switch (c->kind) {
    case K_BOOL: c->stype = CTOOLS_DTA_BYTE; break;
    case K_NULL: c->stype = CTOOLS_DTA_BYTE; break;
    case K_INT: case K_UINT: case K_HUGE: case K_UHUGE: case K_DECIMAL:
        if (c->kind == K_DECIMAL && c->scale) {
            c->stype = CTOOLS_DTA_DOUBLE;
            break;
        }
        if (!c->any || (c->imin >= -127 && c->imax <= 100)) c->stype = CTOOLS_DTA_BYTE;
        else if (c->imin >= -32767 && c->imax <= 32740) c->stype = CTOOLS_DTA_INT;
        else if (c->imin >= -2147483647 && c->imax <= 2147483620) c->stype = CTOOLS_DTA_LONG;
        else {
            /* whole numbers past long: show every digit */
            c->stype = CTOOLS_DTA_DOUBLE;
            snprintf(c->format, sizeof(c->format), "%%20.0g");
        }
        break;
    case K_FLOAT: c->stype = CTOOLS_DTA_FLOAT; break;
    case K_DOUBLE: c->stype = CTOOLS_DTA_DOUBLE; break;
    case K_DATE:
        c->stype = CTOOLS_DTA_LONG;
        snprintf(c->format, sizeof(c->format), "%%td");
        break;
    case K_TS:
        c->stype = CTOOLS_DTA_DOUBLE;
        snprintf(c->format, sizeof(c->format), "%%tc");
        break;
    case K_TIME:
        c->stype = CTOOLS_DTA_DOUBLE;
        snprintf(c->format, sizeof(c->format), "%%tcHH:MM:SS");
        break;
    case K_ENUM:
        c->stype = c->ndict <= 100 ? CTOOLS_DTA_BYTE : c->ndict <= 32740 ? CTOOLS_DTA_INT : CTOOLS_DTA_LONG;
        break;
    case K_STRING: case K_BLOB:
        if (c->kind == K_BLOB || c->has_nul || c->maxlen > (size_t)strl_min) {
            c->stype = CTOOLS_DTA_STRL;
        } else {
            c->stype = CTOOLS_DTA_STR;
            c->width = c->maxlen ? (int)c->maxlen : 1;
        }
        break;
    default: c->stype = CTOOLS_DTA_STRL; break;
    }
}

static int numeric_stype(ctools_dta_type t)
{
    return t != CTOOLS_DTA_STR && t != CTOOLS_DTA_STRL;
}

static int format_ok(const char *f, int string)
{
    size_t n = strlen(f);
    if (n < 2 || n > 56 || f[0] != '%') return 0;
    for (size_t i = 0; i < n; i++)
        if ((unsigned char)f[i] < 0x20 || (unsigned char)f[i] > 0x7e) return 0;
    return string ? f[n - 1] == 's' : f[n - 1] != 's';
}

static int fmt_kind(const char *f, const char *code)
{
    if (f[0] != '%') return 0;
    f++;
    if (*f == '-') f++;
    if (!strcmp(code, "td") && f[0] == 'd') return 1;
    return !strncmp(f, code, strlen(code));
}

/* Cut text to at most 80 characters (whole UTF-8 characters). */
static void label_copy(char *dst, const char *src)
{
    size_t o = 0;
    int chars = 0;
    while (*src && chars < 80 && o < CSQL_LABEL_BYTES - 5) {
        unsigned char ch = (unsigned char)*src;
        int len = ch < 0x80 ? 1 : (ch >> 5) == 6 ? 2 : (ch >> 4) == 14 ? 3 : (ch >> 3) == 30 ? 4 : 1;
        int k = 0;
        for (; k < len && src[k]; k++) dst[o++] = src[k];
        src += k;
        chars++;
    }
    dst[o] = '\0';
}

/* Labels and formats from the first source read that has the variable. */
static void inherit(csql_ctx *ctx, res_col *c, int labels, const csql_source **label_src)
{
    const csql_source *best = NULL;
    const csql_meta_var *bv = NULL;
    const char *want = c->orig ? c->orig : c->name;
    *label_src = NULL;
    /* the column as SQL named it; an exact match first, then one that
     * differs only in case (SELECT PRICE reads price) */
    for (int pass = 0; pass < 2 && !bv; pass++)
        for (const csql_source *s = ctx->sources; s; s = s->next) {
            if (!s->used) continue;
            for (size_t j = 0; j < s->nvars; j++)
                if (pass ? csql_name_eq_ci(s->vars[j].sqlname, want) : !strcmp(s->vars[j].sqlname, want)) {
                    if (!best || s->order < best->order) best = s, bv = &s->vars[j];
                    break;
                }
        }
    if (bv) {
        int src_string = bv->type == CTOOLS_DTA_CODE_STRL || bv->type <= CSQL_STR_MAX;
        int dst_string = !numeric_stype(c->stype);
        if (labels && bv->label[0]) label_copy(c->label, bv->label);
        if (!src_string && !dst_string && format_ok(bv->format, 0)) {
            /* a date keeps its date format, a plain number its number format;
             * a %td or %tc source read as DATE or TIMESTAMP gives no format to
             * a plain number of the same name (year(date) AS date, say) */
            int src_date = fmt_kind(bv->format, "td"), src_tc = fmt_kind(bv->format, "tc");
            int copy;
            if (c->kind == K_DATE) copy = src_date;
            else if (c->kind == K_TS) copy = src_tc;
            else if (c->kind == K_TIME || c->kind == K_ENUM || c->kind == K_BOOL) copy = 0;
            else copy = !(ctx->datetime && (src_date || src_tc));
            if (copy) snprintf(c->format, sizeof(c->format), "%s", bv->format);
        }
        if (src_string && dst_string && bv->format[1] == '-' && c->stype == CTOOLS_DTA_STR)
            snprintf(c->format, sizeof(c->format), "%%-%ds", c->width < 9 ? 9 : c->width);
        if (labels && bv->vlabel[0] && !dst_string && c->kind != K_ENUM && c->kind != K_DATE &&
            c->kind != K_TS && c->kind != K_TIME &&
            (c->stype == CTOOLS_DTA_BYTE || c->stype == CTOOLS_DTA_INT || c->stype == CTOOLS_DTA_LONG) &&
            csql_source_label(best, bv->vlabel)) {
            snprintf(c->vlabel, sizeof(c->vlabel), "%s", bv->vlabel);
            *label_src = best;
        }
    }
    if (c->renamed && !c->label[0] && c->orig) label_copy(c->label, c->orig);
}

/* Whether a label written already holds the values and texts of vl. */
static int label_same(const ctools_dta_value_label *a, const csql_vlabel *vl)
{
    if (a->n != vl->n) return 0;
    for (size_t k = 0; k < vl->n; k++)
        if (a->values[k] != vl->values[k] || strcmp(a->texts[k], vl->texts[k])) return 0;
    return 1;
}

/* A value label name for base: base itself when no label written or
 * inherited by a column has it, else base_1, base_2, ... (base cut to whole
 * characters so that the name stays within 32 characters). */
static void label_name_free(char out[CSQL_NAME_BYTES], const char *base, const res_col *cols,
                            const csql_source *const *col_src, size_t ncols,
                            const ctools_dta_value_label *labs, size_t nlabs)
{
    snprintf(out, CSQL_NAME_BYTES, "%s", base);
    for (int k = 1;; k++) {
        int taken = 0;
        for (size_t l = 0; l < nlabs && !taken; l++) taken = !strcmp(labs[l].name, out);
        for (size_t j = 0; j < ncols && !taken; j++)
            taken = col_src[j] && cols[j].kind != K_ENUM && cols[j].vlabel[0] && !strcmp(cols[j].vlabel, out);
        if (!taken) return;
        char suffix[16];
        int sl = snprintf(suffix, sizeof(suffix), "_%d", k);
        size_t o = 0, chars = 0;
        for (const char *p = base; *p;) {
            size_t len = 1;
            while (p[len] && ((unsigned char)p[len] & 0xC0) == 0x80) len++;
            if (chars + 1 + (size_t)sl > 32 || o + len + (size_t)sl >= CSQL_NAME_BYTES) break;
            memcpy(out + o, p, len);
            o += len, p += len, chars++;
        }
        memcpy(out + o, suffix, (size_t)sl + 1);
    }
}

/* ------------------------------------------------------------ producer */

typedef struct
{
    res_col *cols;
    size_t ncols;
    const res_rows *rows;
    size_t chunks_per_group;
} prod_ctx;

typedef struct
{
    void **cols;
    duckdb_data_chunk *owned;
    size_t nowned;
} prod_group;

static void prod_group_free(const prod_ctx *p, prod_group *g)
{
    if (!g) return;
    for (size_t j = 0; g->cols && j < p->ncols; j++) free(g->cols[j]);
    for (size_t k = 0; k < g->nowned; k++) DK(destroy_data_chunk)(&g->owned[k]);
    free(g->owned);
    free(g->cols);
    free(g);
}

/* Chunk ch of the rows: held, or fetched from the materialized result
 * (*owned set: the caller destroys it). */
static duckdb_data_chunk rows_chunk(const res_rows *r, size_t ch, int *owned)
{
    *owned = r->random;
    return r->random ? DK(result_get_chunk)(*r->res, (idx_t)ch) : r->chunks[ch];
}

static int prod_fill(void *arg, size_t group, int phase, ctools_dta_group *out)
{
    prod_ctx *p = arg;
    const res_rows *r = p->rows;
    size_t c0 = group * p->chunks_per_group, c1 = c0 + p->chunks_per_group;
    if (c1 > r->nchunks) c1 = r->nchunks;
    uint64_t n = r->starts[c1] - r->starts[c0];
    int strings = 0;
    prod_group *g = calloc(1, sizeof(*g));
    if (!g || !(g->cols = calloc(p->ncols ? p->ncols : 1, sizeof(void *))) ||
        (r->random && !(g->owned = calloc(c1 - c0 ? c1 - c0 : 1, sizeof(duckdb_data_chunk))))) {
        prod_group_free(p, g);
        return 920;
    }
    for (size_t j = 0; j < p->ncols; j++) {
        const res_col *c = &p->cols[j];
        int string = !numeric_stype(c->stype);
        out->columns[j] = NULL;
        if (phase == CTOOLS_DTA_PHASE_STRLS && c->stype != CTOOLS_DTA_STRL) continue;
        g->cols[j] = malloc((size_t)(n ? n : 1) * (string ? sizeof(ctools_dta_span) : sizeof(double)));
        if (!g->cols[j]) {
            prod_group_free(p, g);
            return 920;
        }
        out->columns[j] = g->cols[j];
        strings |= string;
    }
    size_t at = 0;
    for (size_t ch = c0; ch < c1; ch++) {
        int owned;
        duckdb_data_chunk chunk = rows_chunk(r, ch, &owned);
        size_t m = (size_t)(r->starts[ch + 1] - r->starts[ch]);
        if (!chunk) {
            prod_group_free(p, g);
            return 920;
        }
        for (size_t j = 0; j < p->ncols; j++) {
            const res_col *c = &p->cols[j];
            int string = !numeric_stype(c->stype);
            void *col = g->cols[j];
            if (!col) continue;
            duckdb_vector vec = DK(data_chunk_get_vector)(chunk, (idx_t)j);
            const void *d = DK(vector_get_data)(vec);
            const uint64_t *valid = DK(vector_get_validity)(vec);
            if (string) {
                ctools_dta_span *sp = (ctools_dta_span *)col + at;
                int binary = c->kind == K_BLOB;
                for (size_t i = 0; i < m; i++) {
                    if (c->kind != K_STRING && c->kind != K_BLOB) {
                        sp[i] = (ctools_dta_span){"", 0, 0};
                        continue;
                    }
                    if (!row_valid(valid, i)) {
                        sp[i] = (ctools_dta_span){"", 0, 0};
                        continue;
                    }
                    const char *s;
                    size_t len;
                    string_at(d, i, &s, &len);
                    sp[i] = (ctools_dta_span){s, len, binary};
                }
            } else {
                double *x = (double *)col + at;
                double miss = csql_missing();
                double lim = c->kind == K_FLOAT ? 0x1p127 : 0x1p1023; /* as in stats_chunk */
                for (size_t i = 0; i < m; i++) {
                    if (c->kind == K_NULL || !row_valid(valid, i)) {
                        x[i] = miss;
                        continue;
                    }
                    double v = value_at(c, d, i);
                    if (!(v < lim && v > -lim)) v = miss;
                    x[i] = v;
                }
            }
        }
        at += m;
        /* string spans point into the chunk: keep it until release */
        if (owned) {
            if (strings) g->owned[g->nowned++] = chunk;
            else DK(destroy_data_chunk)(&chunk);
        }
    }
    out->state = g;
    return 0;
}

static void prod_release(void *arg, ctools_dta_group *group)
{
    prod_group_free(arg, group->state);
    group->state = NULL;
}

/* ------------------------------------------------------------ display */

static void civil_from_days(int64_t z, int *y, unsigned *m, unsigned *d)
{
    z += 719468;
    int64_t era = (z >= 0 ? z : z - 146096) / 146097;
    unsigned doe = (unsigned)(z - era * 146097);
    unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    int64_t yy = (int64_t)yoe + era * 400;
    unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    unsigned mp = (5 * doy + 2) / 153;
    *d = doy - (153 * mp + 2) / 5 + 1;
    *m = mp < 10 ? mp + 3 : mp - 9;
    *y = (int)(yy + (*m <= 2));
}

/* Text of one cell for the listing; returns 1 for right alignment. */
static int cell_text(const res_col *c, duckdb_vector vec, size_t i, char *buf, size_t cap)
{
    const void *d = DK(vector_get_data)(vec);
    const uint64_t *valid = DK(vector_get_validity)(vec);
    int numeric = c->kind != K_STRING && c->kind != K_BLOB && c->kind != K_ENUM && c->kind != K_OTHER;
    if (!row_valid(valid, i) || c->kind == K_NULL) {
        snprintf(buf, cap, numeric ? "." : "");
        return numeric;
    }
    switch (c->kind) {
    case K_BOOL: snprintf(buf, cap, "%s", ((const bool *)d)[i] ? "true" : "false"); return 1;
    case K_INT: snprintf(buf, cap, "%lld", (long long)int_at(d, c->bytes, i)); return 1;
    case K_UINT: snprintf(buf, cap, "%llu", (unsigned long long)uint_at(d, c->bytes, i)); return 1;
    case K_HUGE: case K_UHUGE: snprintf(buf, cap, "%.0f", value_at(c, d, i)); return 1;
    case K_FLOAT: snprintf(buf, cap, "%.7g", value_at(c, d, i)); return 1;
    case K_DOUBLE: case K_DECIMAL: {
        double v = value_at(c, d, i);
        if (c->kind == K_DECIMAL && c->scale <= 9) snprintf(buf, cap, "%.*f", c->scale, v);
        else snprintf(buf, cap, "%.9g", v);
        return 1;
    }
    case K_DATE: {
        int32_t days = ((const duckdb_date *)d)[i].days;
        int y;
        unsigned m, dd;
        civil_from_days(days, &y, &m, &dd);
        snprintf(buf, cap, "%04d-%02u-%02u", y, m, dd);
        return 1;
    }
    case K_TS: case K_TIME: {
        int64_t v = ((const int64_t *)d)[i];
        double per_s = c->div * 1000.0; /* units per second */
        int64_t ups = (int64_t)per_s;
        int64_t secs = ups > 0 ? (v >= 0 ? v / ups : -((-v + ups - 1) / ups)) : v * 1000;
        int64_t frac = ups > 0 ? v - secs * ups : 0;
        int64_t days = secs >= 0 ? secs / 86400 : -((-secs + 86399) / 86400);
        int64_t sod = secs - days * 86400;
        char tail[32] = "";
        if (frac) {
            int digits = ups == 1000 ? 3 : ups == 1000000 ? 6 : 9;
            snprintf(tail, sizeof(tail), ".%0*lld", digits, (long long)frac);
        }
        if (c->kind == K_TIME)
            snprintf(buf, cap, "%02lld:%02lld:%02lld%s", (long long)(sod / 3600), (long long)(sod / 60 % 60),
                     (long long)(sod % 60), tail);
        else {
            int y;
            unsigned m, dd;
            civil_from_days(days, &y, &m, &dd);
            snprintf(buf, cap, "%04d-%02u-%02u %02lld:%02lld:%02lld%s", y, m, dd, (long long)(sod / 3600),
                     (long long)(sod / 60 % 60), (long long)(sod % 60), tail);
        }
        return 1;
    }
    case K_ENUM: {
        uint32_t k = enum_at(d, c->internal, i);
        snprintf(buf, cap, "%s", c->dict && k < c->ndict && c->dict[k] ? c->dict[k] : "");
        return 0;
    }
    case K_STRING: case K_BLOB: {
        const char *s;
        size_t len;
        string_at(d, i, &s, &len);
        size_t o = 0;
        for (size_t k = 0; k < len && o + 1 < cap; k++) {
            unsigned char ch = (unsigned char)s[k];
            buf[o++] = ch < 0x20 || ch == 0x7f ? ' ' : (char)ch;
        }
        buf[o] = '\0';
        return 0;
    }
    default: snprintf(buf, cap, "?"); return 0;
    }
}

static size_t utf8_width(const char *s)
{
    size_t w = 0;
    for (; *s; s++)
        if (((unsigned char)*s & 0xC0) != 0x80) w++;
    return w;
}

/* Cut s to at most w characters, marking the cut with "~". */
static void utf8_cut(char *s, size_t w)
{
    size_t chars = 0;
    for (char *p = s; *p; p++) {
        if (((unsigned char)*p & 0xC0) != 0x80) {
            if (chars == w - 1 && utf8_width(p) > 1) {
                p[0] = '~';
                p[1] = '\0';
                return;
            }
            chars++;
        }
    }
}

/* SMCL-safe copy: braces become {c -(} and {c )-}. */
static void smcl_add(csql_buf *b, const char *s)
{
    for (; *s; s++) {
        if (*s == '{') csql_buf_str(b, "{c -(}");
        else if (*s == '}') csql_buf_str(b, "{c )-}");
        else csql_buf_add(b, s, 1);
    }
}

static void pad(csql_buf *b, size_t n)
{
    for (size_t i = 0; i < n; i++) csql_buf_add(b, " ", 1);
}

static int display_rows(res_col *cols, size_t ncols, const res_rows *r, uint64_t shown, uint64_t total,
                        int linesize)
{
    enum { CELL = 128, MAXW = 30 };
    size_t *width = calloc(ncols ? ncols : 1, sizeof(size_t));
    char *cells = calloc((size_t)(shown ? shown : 1) * (ncols ? ncols : 1), CELL);
    unsigned char *right = calloc((size_t)(shown ? shown : 1) * (ncols ? ncols : 1), 1);
    if (!width || !cells || !right) {
        free(width), free(cells), free(right);
        return 920;
    }
    uint64_t row = 0;
    for (size_t ch = 0; ch < r->nchunks && row < shown; ch++) {
        size_t m = (size_t)(r->starts[ch + 1] - r->starts[ch]);
        for (size_t i = 0; i < m && row < shown; i++, row++)
            for (size_t j = 0; j < ncols; j++) {
                char *cell = cells + ((size_t)row * ncols + j) * CELL;
                duckdb_vector vec = DK(data_chunk_get_vector)(r->chunks[ch], (idx_t)j);
                right[(size_t)row * ncols + j] = (unsigned char)cell_text(&cols[j], vec, i, cell, CELL);
                if (utf8_width(cell) > MAXW) utf8_cut(cell, MAXW);
            }
    }
    for (size_t j = 0; j < ncols; j++) {
        width[j] = utf8_width(cols[j].orig ? cols[j].orig : cols[j].name);
        if (width[j] > MAXW) width[j] = MAXW;
        for (uint64_t k = 0; k < shown; k++) {
            size_t w = utf8_width(cells + ((size_t)k * ncols + j) * CELL);
            if (w > width[j]) width[j] = w;
        }
    }
    char num[32];
    int rw = snprintf(num, sizeof(num), "%llu", (unsigned long long)(shown ? shown : 1));
    size_t lead = (size_t)rw + 2;
    if (linesize < 40) linesize = 40;
    csql_buf b = {0};
    SF_display("\n");
    for (size_t j0 = 0; j0 < ncols;) {
        size_t used = lead, j1 = j0;
        while (j1 < ncols && (j1 == j0 || used + 3 + width[j1] <= (size_t)linesize)) used += 3 + width[j1++];
        /* header */
        b.len = 0;
        csql_buf_str(&b, "{txt}");
        pad(&b, lead);
        for (size_t j = j0; j < j1; j++) {
            char name[CELL];
            snprintf(name, sizeof(name), "%s", cols[j].orig ? cols[j].orig : cols[j].name);
            if (utf8_width(name) > width[j]) utf8_cut(name, width[j]);
            size_t w = utf8_width(name);
            int numeric = numeric_stype(cols[j].stype) && cols[j].kind != K_ENUM;
            pad(&b, 3);
            if (numeric) pad(&b, width[j] - w);
            smcl_add(&b, name);
            if (!numeric) pad(&b, width[j] - w);
        }
        csql_buf_str(&b, "\n");
        SF_display(b.p);
        b.len = 0;
        csql_buf_str(&b, "{txt}");
        pad(&b, lead);
        for (size_t j = j0; j < j1; j++) {
            pad(&b, 3);
            char rule[64];
            snprintf(rule, sizeof(rule), "{hline %zu}", width[j]);
            csql_buf_str(&b, rule);
        }
        csql_buf_str(&b, "\n");
        SF_display(b.p);
        for (uint64_t k = 0; k < shown; k++) {
            b.len = 0;
            snprintf(num, sizeof(num), "%*llu.", rw, (unsigned long long)(k + 1));
            csql_buf_str(&b, "{txt}");
            pad(&b, 1);
            csql_buf_str(&b, num);
            csql_buf_str(&b, "{res}");
            for (size_t j = j0; j < j1; j++) {
                const char *cell = cells + ((size_t)k * ncols + j) * CELL;
                size_t w = utf8_width(cell);
                pad(&b, 3);
                if (right[(size_t)k * ncols + j]) pad(&b, width[j] - w);
                smcl_add(&b, cell);
                if (!right[(size_t)k * ncols + j]) pad(&b, width[j] - w);
            }
            csql_buf_str(&b, "\n");
            SF_display(b.p);
        }
        SF_display("\n");
        j0 = j1;
    }
    b.len = 0;
    char foot[160];
    if (shown < total)
        snprintf(foot, sizeof(foot), "{txt}(%llu rows, %zu columns; first %llu shown)\n",
                 (unsigned long long)total, ncols, (unsigned long long)shown);
    else
        snprintf(foot, sizeof(foot), "{txt}(%llu row%s, %zu column%s)\n", (unsigned long long)total,
                 total == 1 ? "" : "s", ncols, ncols == 1 ? "" : "s");
    SF_display(foot);
    csql_buf_free(&b);
    free(width), free(cells), free(right);
    return 0;
}

/* ------------------------------------------------------------ output */

static int rows_push(res_rows *r, duckdb_data_chunk chunk, size_t n)
{
    if (r->nchunks + 1 >= r->cap) {
        size_t cap = r->cap ? r->cap * 2 : 256;
        duckdb_data_chunk *c = realloc(r->chunks, cap * sizeof(*c));
        if (!c) return 920;
        r->chunks = c;
        uint64_t *s = realloc(r->starts, (cap + 1) * sizeof(*s));
        if (!s) return 920;
        r->starts = s;
        r->cap = cap;
    }
    if (!r->nchunks) r->starts[0] = 0;
    r->chunks[r->nchunks++] = chunk;
    r->nrows += n;
    r->starts[r->nchunks] = r->nrows;
    return 0;
}

static void rows_free(res_rows *r)
{
    for (size_t i = 0; r->chunks && i < r->nchunks; i++) DK(destroy_data_chunk)(&r->chunks[i]);
    free(r->chunks);
    free(r->starts);
    memset(r, 0, sizeof(*r));
}

static void cols_free(res_col *cols, size_t n)
{
    for (size_t j = 0; cols && j < n; j++) {
        for (size_t k = 0; cols[j].dict && k < cols[j].ndict; k++) DK(free)(cols[j].dict[k]);
        free(cols[j].dict);
    }
    free(cols);
}

int csql_result_output(csql_ctx *ctx, duckdb_result *res, int streaming, const csql_out_opts *opts,
                       int64_t *nrows_out, int64_t *ncols_out)
{
    (void)streaming;
    size_t ncols = (size_t)DK(column_count)(res);
    res_col *cols = calloc(ncols ? ncols : 1, sizeof(*cols));
    res_rows rows = {0};
    int rc = 0;
    double t0 = ctools_timer_seconds();
    if (!cols) return 920;
    char (*names)[CSQL_NAME_BYTES] = calloc(ncols ? ncols : 1, CSQL_NAME_BYTES);
    if (!names) {
        free(cols);
        return 920;
    }
    for (size_t j = 0; j < ncols; j++) {
        res_col *c = &cols[j];
        duckdb_logical_type lt = DK(column_logical_type)(res, (idx_t)j);
        c->kind = kind_of(lt, c);
        DK(destroy_logical_type)(&lt);
        const char *nm = DK(column_name)(res, (idx_t)j);
        c->orig = (char *)nm;
        c->renamed = csql_stata_name(nm ? nm : "", c->name, names, j);
        memcpy(names[j], c->name, CSQL_NAME_BYTES);
        if (c->kind == K_OTHER && opts->mode == CSQL_OUT_DTA) {
            snprintf(ctx->error, sizeof(ctx->error),
                     "column %s has a type csql cannot store in Stata; cast it in the query "
                     "(for example CAST(%s AS VARCHAR))",
                     nm ? nm : "", nm ? nm : "");
            rc = 459;
        }
    }
    free(names);

    /* rows */
    uint64_t limit = opts->mode == CSQL_OUT_DISPLAY ? (uint64_t)(opts->rows > 0 ? opts->rows : 0) : UINT64_MAX;
    uint64_t total = 0;
    int random = opts->mode == CSQL_OUT_DTA && !streaming && DK(result_get_chunk) && DK(result_chunk_count);
    if (!rc && random) {
        rows.random = 1;
        rows.res = res;
        rows.nchunks = (size_t)DK(result_chunk_count)(*res);
        rows.starts = calloc(rows.nchunks + 1, sizeof(uint64_t));
        if (!rows.starts) rc = 920;
    }
    while (!rc && !random) {
        duckdb_data_chunk chunk = DK(fetch_chunk)(*res);
        if (!chunk) break;
        size_t n = (size_t)DK(data_chunk_get_size)(chunk);
        if (n == 0) {
            DK(destroy_data_chunk)(&chunk);
            continue;
        }
        total += n;
        if (opts->mode == CSQL_OUT_DISPLAY && rows.nrows >= limit) {
            DK(destroy_data_chunk)(&chunk);
        } else if ((rc = rows_push(&rows, chunk, n)) != 0) {
            DK(destroy_data_chunk)(&chunk);
        }
        if (!rc && ctools_cancel_requested()) {
            atomic_store(&ctx->cancelled, 1);
            DK(interrupt)(ctx->con);
            rc = 1;
        }
    }
    if (!rc) {
        const char *err = DK(result_error)(res);
        if (err && *err) {
            snprintf(ctx->error, sizeof(ctx->error), "%s", err);
            rc = 198;
        }
    }
    double t_fetch = ctools_timer_seconds();

    if (!rc && opts->mode == CSQL_OUT_DISPLAY) {
        for (size_t j = 0; j < ncols; j++) decide(&cols[j], CSQL_STR_MAX);
        uint64_t shown = rows.nrows < limit ? rows.nrows : limit;
        if (ncols) rc = display_rows(cols, ncols, &rows, shown, total, opts->linesize);
    } else if (!rc) {
        /* measure, in parallel over chunks */
        col_stats *part = NULL;
        int nt = ctx->threads > 0 ? ctx->threads : 1;
        if (nt > (int)rows.nchunks) nt = rows.nchunks ? (int)rows.nchunks : 1;
        part = calloc((size_t)nt * (ncols ? ncols : 1), sizeof(col_stats));
        if (!part) rc = 920;
        if (!rc) {
            int failed = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(dynamic, 8)
#endif
            for (long ch = 0; ch < (long)rows.nchunks; ch++) {
#ifdef _OPENMP
                int tid = omp_get_thread_num();
#else
                int tid = 0;
#endif
                int owned;
                duckdb_data_chunk chunk = rows_chunk(&rows, (size_t)ch, &owned);
                if (!chunk) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
                    failed = 1;
                    continue;
                }
                size_t m = owned ? (size_t)DK(data_chunk_get_size)(chunk)
                                 : (size_t)(rows.starts[ch + 1] - rows.starts[ch]);
                if (owned) rows.starts[ch + 1] = m; /* sizes; summed below */
                for (size_t j = 0; j < ncols; j++)
                    stats_chunk(&cols[j], DK(data_chunk_get_vector)(chunk, (idx_t)j), m,
                                &part[(size_t)tid * ncols + j]);
                if (owned) DK(destroy_data_chunk)(&chunk);
            }
            if (failed) rc = 920;
            if (random) {
                rows.starts[0] = 0;
                for (size_t ch = 0; ch < rows.nchunks; ch++) rows.starts[ch + 1] += rows.starts[ch];
                rows.nrows = total = rows.starts[rows.nchunks];
            }
            for (size_t j = 0; j < ncols; j++) {
                col_stats s = {0};
                for (int t = 0; t < nt; t++) stats_merge(&s, &part[(size_t)t * ncols + j]);
                cols[j].imin = s.imin, cols[j].imax = s.imax, cols[j].dmin = s.dmin, cols[j].dmax = s.dmax;
                cols[j].any = s.any, cols[j].maxlen = s.maxlen, cols[j].has_nul = s.has_nul;
                cols[j].imprecise = s.imprecise, cols[j].nonfinite = s.nonfinite;
            }
        }
        free(part);
        double t_stats = ctools_timer_seconds();

        /* the dataset's description */
        ctools_dta_column *dcols = calloc(ncols ? ncols : 1, sizeof(*dcols));
        ctools_dta_value_label *labs = calloc(ncols ? ncols : 1, sizeof(*labs));
        size_t nlabs = 0;
        int32_t **enum_vals = calloc(ncols ? ncols : 1, sizeof(int32_t *));
        /* per label written: the source it came from (NULL for an ENUM) and
         * its name there */
        const csql_source **lab_src = calloc(ncols ? ncols : 1, sizeof(*lab_src));
        const char **lab_orig = calloc(ncols ? ncols : 1, sizeof(*lab_orig));
        const csql_source **col_src = calloc(ncols ? ncols : 1, sizeof(*col_src));
        if (!dcols || !labs || !enum_vals || !lab_src || !lab_orig || !col_src) rc = rc ? rc : 920;
        for (size_t j = 0; !rc && j < ncols; j++) {
            decide(&cols[j], opts->strl_min > 0 && opts->strl_min < CSQL_STR_MAX ? opts->strl_min : CSQL_STR_MAX);
            inherit(ctx, &cols[j], opts->labels, &col_src[j]);
        }
        /* Value label names: an inherited label keeps its name; another
         * label of that name (an ENUM's, or one from another source with
         * different contents) gets name_1, name_2, ..., a name no label
         * written or inherited has. */
        for (size_t j = 0; !rc && j < ncols; j++) {
            res_col *c = &cols[j];
            const csql_source *lsrc = col_src[j];
            if (c->kind == K_ENUM && c->ndict) {
                char lname[CSQL_NAME_BYTES];
                label_name_free(lname, c->name, cols, col_src, ncols, labs, nlabs);
                enum_vals[j] = malloc(c->ndict * sizeof(int32_t));
                if (!enum_vals[j]) {
                    rc = 920;
                    break;
                }
                for (size_t k = 0; k < c->ndict; k++) enum_vals[j][k] = (int32_t)(k + 1);
                snprintf(c->vlabel, sizeof(c->vlabel), "%s", lname);
                lab_src[nlabs] = NULL, lab_orig[nlabs] = NULL;
                labs[nlabs++] = (ctools_dta_value_label){c->vlabel, c->ndict, enum_vals[j],
                                                         (const char *const *)c->dict};
            } else if (c->vlabel[0] && lsrc) {
                const csql_vlabel *vl = csql_source_label(lsrc, c->vlabel);
                size_t l = 0;
                /* the same label of the same source, written already */
                for (; l < nlabs; l++)
                    if (lab_src[l] == lsrc && !strcmp(lab_orig[l], vl->name)) break;
                /* else one of the same name: shared when it says the same */
                if (l == nlabs)
                    for (l = 0; l < nlabs; l++)
                        if (!strcmp(labs[l].name, vl->name)) {
                            if (!label_same(&labs[l], vl)) l = nlabs;
                            break;
                        }
                if (l < nlabs) {
                    if (strcmp(labs[l].name, c->vlabel)) snprintf(c->vlabel, sizeof(c->vlabel), "%s", labs[l].name);
                } else {
                    int clash = 0;
                    for (size_t k = 0; k < nlabs && !clash; k++) clash = !strcmp(labs[k].name, vl->name);
                    if (clash) {
                        char lname[CSQL_NAME_BYTES];
                        label_name_free(lname, vl->name, cols, col_src, ncols, labs, nlabs);
                        snprintf(c->vlabel, sizeof(c->vlabel), "%s", lname);
                    }
                    lab_src[nlabs] = lsrc, lab_orig[nlabs] = vl->name;
                    labs[nlabs++] = (ctools_dta_value_label){clash ? c->vlabel : vl->name, vl->n, vl->values,
                                                             (const char *const *)vl->texts};
                }
            }
            dcols[j] = (ctools_dta_column){.name = c->name,
                                           .type = c->stype,
                                           .width = c->stype == CTOOLS_DTA_STR ? c->width : 0,
                                           .format = c->format[0] ? c->format : NULL,
                                           .label = c->label[0] ? c->label : NULL,
                                           .value_label = c->vlabel[0] ? c->vlabel : NULL,
                                           .source = CTOOLS_DTA_SRC_PRODUCER,
                                           .flags = c->kind == K_STRING && !c->has_nul ? CTOOLS_DTA_COL_TEXT : 0};
        }
        /* notes the user should see */
        for (size_t j = 0; !rc && j < ncols; j++) {
            if (cols[j].imprecise)
                ctools_msg("", "{txt}note: %s has %llu integer value%s beyond 2^53, stored as double "
                           "and rounded", cols[j].name, (unsigned long long)cols[j].imprecise,
                           cols[j].imprecise == 1 ? "" : "s");
            if (cols[j].nonfinite)
                ctools_msg("", "{txt}note: %s has %llu value%s that Stata cannot hold (NaN, infinite or "
                           "too large), stored as missing", cols[j].name,
                           (unsigned long long)cols[j].nonfinite, cols[j].nonfinite == 1 ? "" : "s");
        }

        if (!rc) {
            size_t per = 32; /* chunks per row group: about 64K rows */
            size_t ngroups = (rows.nchunks + per - 1) / per;
            uint64_t *starts = malloc((ngroups + 1) * sizeof(uint64_t));
            if (!starts) rc = 920;
            for (size_t g = 0; !rc && g <= ngroups; g++) {
                size_t ch = g * per < rows.nchunks ? g * per : rows.nchunks;
                starts[g] = rows.nchunks ? rows.starts[ch] : 0;
            }
            prod_ctx pc = {cols, ncols, &rows, per};
            ctools_dta_producer pr = {ngroups, starts, NULL, &pc, prod_fill, prod_release};
            ctools_dta_spec spec;
            ctools_dta_result wres;
            ctools_dta_spec_init(&spec);
            spec.nobs = rows.nrows;
            spec.ncols = ncols;
            spec.cols = dcols;
            spec.labels = labs;
            spec.nlabels = nlabs;
            spec.producer = &pr;
            spec.threads = ctx->threads;
            spec.temporary = opts->temporary;
            if (!rc) {
                memset(&wres, 0, sizeof(wres));
                rc = ctools_dta_write(opts->path, &spec, &wres);
                if (rc) snprintf(ctx->error, sizeof(ctx->error), "%s", wres.message[0] ? wres.message : "could not write the result");
            }
            free(starts);
            if (ctx->verbose)
                ctools_msg("csql", "rows %.3f s, types %.3f s, .dta %.3f s (%.0f MB)", t_fetch - t0,
                           t_stats - t_fetch, ctools_timer_seconds() - t_stats,
                           (double)wres.file_bytes / 1e6);
        }
        for (size_t j = 0; enum_vals && j < ncols; j++) free(enum_vals[j]);
        free(enum_vals);
        free(dcols);
        free(labs);
        free(lab_src);
        free(lab_orig);
        free(col_src);
    }
    if (!rc) {
        *nrows_out = (int64_t)total;
        *ncols_out = (int64_t)ncols;
        /* renamed columns, for the ado's report */
        csql_buf b = {0};
        for (size_t j = 0; j < ncols; j++)
            if (cols[j].renamed && opts->mode == CSQL_OUT_DTA) {
                csql_buf_str(&b, b.len ? " " : "");
                csql_buf_str(&b, cols[j].name);
            }
        SF_macro_save("_csql_renamed", b.p ? b.p : "");
        csql_buf_free(&b);
    }
    rows_free(&rows);
    cols_free(cols, ncols);
    return rc;
}

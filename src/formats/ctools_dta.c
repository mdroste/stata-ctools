/*
 * ctools_dta.c - Stata .dta writer (formats 118 and 119)
 *
 * Passes over the rows, each split into units of contiguous rows:
 *   plan   (only when needed) longest str# values and strL bytes per range
 *   cells  the data records, one ~1 MiB block per unit
 *   strls  the GSO records of the strL cells that define their value
 * Units are built in parallel a wave at a time. While the pool builds wave
 * k + 1, one I/O thread appends wave k to the file in order, so building
 * overlaps writing. Header and trailer go out from the calling thread, and
 * the map is patched last, inside the bytes already written.
 */
#include "formats/ctools_dta.h"

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ctools_config.h"
#include "ctools_file.h"
#include "ctools_threads.h"

double ctools_timer_seconds(void); /* ctools_runtime.c */

#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif

#define DTA_BLOCK_BYTES ((size_t)1 << 20)          /* data records per unit */
#define DTA_PLAN_ROWS ((uint64_t)65536)            /* rows per plan unit */
#define DTA_GSO_UNIT_BYTES ((uint64_t)4 << 20)     /* copied strL bytes per unit */
#define DTA_EXTERNAL_BYTES ((size_t)64 << 10)      /* written from their source */
#define DTA_WRITE_PIECE ((size_t)4 << 20)          /* measured sweet spot: 1-8 MiB */
#define DTA_WAVE_BYTES ((uint64_t)256 << 20)       /* built bytes in flight per wave */
/* Units built at once from columns in memory. Encoding outruns the single
 * appending thread, and more encoders only take memory bandwidth from its
 * copy into the file cache (measured: 1.6 GB of doubles appended in 0.49 s
 * with 4 encoders, 0.54 s with 16). Producer passes, and columns gathered
 * through rows (random reads, which wait on memory), use every thread. */
#define DTA_IO_WAVE_UNITS 4
#define DTA_PRODUCER_WAVE_BYTES ((uint64_t)512 << 20)
#define DTA_KEEP_BYTES ((size_t)1 << 30)
#define DTA_GSO_HEAD 20 /* "GSO", v u32, o u64, t u8, len u32 */
#define DTA_MAX_STRL ((uint64_t)2000000000)
#define DTA_MAX_STR 2045

/* Stata's missing values: base codes for ., then +1 step per .a ... .z */
#define DTA_D_MISS UINT64_C(0x7FE0000000000000)
#define DTA_F_MISS UINT32_C(0x7F000000)
#define DTA_L_MISS INT32_C(2147483621)
#define DTA_I_MISS 32741
#define DTA_B_MISS 101

/* ------------------------------------------------------------------------
 * File
 * ------------------------------------------------------------------------ */

typedef struct
{
#ifdef _WIN32
    HANDLE h;
#else
    int fd;
#endif
    uint64_t pos;
    double busy; /* seconds spent in writes */
    int open;
} dta_file;

static int dta_file_open(dta_file *f, const char *path, int temporary)
{
    memset(f, 0, sizeof(*f));
#ifdef _WIN32
    DWORD attr = (temporary ? FILE_ATTRIBUTE_TEMPORARY : FILE_ATTRIBUTE_NORMAL) | FILE_FLAG_SEQUENTIAL_SCAN;
    f->h = ctools_create_file(path, GENERIC_WRITE, 0, CREATE_ALWAYS, attr);
    if (f->h == INVALID_HANDLE_VALUE) return -1;
#else
    (void)temporary;
    f->fd = open(path, O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0666);
    if (f->fd < 0) return -1;
#endif
    f->open = 1;
    return 0;
}

static int dta_file_write(dta_file *f, const void *buf, size_t len)
{
    const char *p = buf;
    const double t0 = ctools_timer_seconds();
    while (len) {
        size_t n = len < DTA_WRITE_PIECE ? len : DTA_WRITE_PIECE;
#ifdef _WIN32
        DWORD done = 0;
        if (!WriteFile(f->h, p, (DWORD)n, &done, NULL) || done == 0) return -1;
        n = done;
#else
        ssize_t done = write(f->fd, p, n);
        if (done < 0 && errno == EINTR) continue;
        if (done <= 0) return -1;
        n = (size_t)done;
#endif
        p += n;
        len -= n;
        f->pos += n;
    }
    f->busy += ctools_timer_seconds() - t0;
    return 0;
}

/* Rewrite bytes already written (the map), leaving the append position. */
static int dta_file_patch(dta_file *f, const void *buf, size_t len, uint64_t at)
{
#ifdef _WIN32
    LARGE_INTEGER off;
    DWORD done = 0;
    off.QuadPart = (LONGLONG)at;
    if (!SetFilePointerEx(f->h, off, NULL, FILE_BEGIN)) return -1;
    if (!WriteFile(f->h, buf, (DWORD)len, &done, NULL) || done != len) return -1;
    off.QuadPart = (LONGLONG)f->pos;
    return SetFilePointerEx(f->h, off, NULL, FILE_BEGIN) ? 0 : -1;
#else
    const char *p = buf;
    while (len) {
        ssize_t done = pwrite(f->fd, p, len, (off_t)at);
        if (done < 0 && errno == EINTR) continue;
        if (done <= 0) return -1;
        p += done;
        len -= (size_t)done;
        at += (uint64_t)done;
    }
    return 0;
#endif
}

static int dta_file_close(dta_file *f)
{
    if (!f->open) return 0;
    f->open = 0;
#ifdef _WIN32
    return CloseHandle(f->h) ? 0 : -1;
#else
    return close(f->fd);
#endif
}

/* ------------------------------------------------------------------------
 * Byte buffer for the header and trailer
 * ------------------------------------------------------------------------ */

typedef struct
{
    unsigned char *p;
    size_t len, cap;
    int failed;
} dta_bytes;

static unsigned char *dta_bytes_grow(dta_bytes *b, size_t n)
{
    if (b->failed) return NULL;
    if (n > SIZE_MAX - b->len) {
        b->failed = 1;
        return NULL;
    }
    if (b->len + n > b->cap) {
        size_t cap = b->cap ? b->cap : 4096;
        while (cap < b->len + n) cap = cap > SIZE_MAX / 2 ? b->len + n : cap * 2;
        unsigned char *p = realloc(b->p, cap);
        if (!p) {
            b->failed = 1;
            return NULL;
        }
        b->p = p;
        b->cap = cap;
    }
    unsigned char *at = b->p + b->len;
    b->len += n;
    return at;
}

static void dta_put(dta_bytes *b, const void *src, size_t n)
{
    unsigned char *at = dta_bytes_grow(b, n);
    if (at) {
        if (src) memcpy(at, src, n);
        else memset(at, 0, n);
    }
}

static void dta_tag(dta_bytes *b, const char *tag) { dta_put(b, tag, strlen(tag)); }

static inline void dta_le(unsigned char *p, uint64_t v, int n)
{
    for (int i = 0; i < n; i++) p[i] = (unsigned char)(v >> (8 * i));
}

static void dta_put_le(dta_bytes *b, uint64_t v, int n)
{
    unsigned char t[8];
    dta_le(t, v, n);
    dta_put(b, t, (size_t)n);
}

/* text in a field of n bytes, zero-padded (the caller checked its length) */
static void dta_put_fixed(dta_bytes *b, const char *text, size_t n)
{
    size_t len = text ? strlen(text) : 0;
    dta_put(b, text, len);
    dta_put(b, NULL, n - len);
}

/* ------------------------------------------------------------------------
 * Writer state
 * ------------------------------------------------------------------------ */

typedef struct
{
    const ctools_dta_column *c;
    uint16_t code;   /* file type code */
    uint32_t width;  /* bytes per cell */
    size_t offset;   /* in the record */
    size_t slot;     /* PRODUCER: entry in ctools_dta_group.columns */
    int measure;     /* str# whose width the plan pass measures */
    /* CTOOLS_DTA_COL_SHARE: for each source row, the first output row that
     * reads it (UINT64_MAX for none); that row's cell defines the value and
     * later ones point to it */
    uint64_t *first;
} dta_col;

typedef struct
{
    uint64_t r0, r1;
    size_t group;   /* producer group */
    uint64_t bytes; /* pacing estimate */
} dta_unit;

typedef struct
{
    size_t at, len;  /* inline bytes buf[at, at + len) ... */
    const char *ext; /* ... or len bytes at ext */
} dta_seg;

typedef struct dta_writer dta_writer;

typedef struct
{
    dta_writer *w;
    const dta_unit *unit;
    size_t slot;  /* fixed position among both sets of tasks */
    size_t index; /* unit number in its pass */
    int phase;
    ctools_dta_group g; /* a fill this task owns */
    int filled;
    int kept; /* uses w->kept[unit->group] */
    void *(*fn)(void *);
    int done; /* built; guarded by the writer's lock */
    unsigned char *buf;
    size_t cap, len;
    dta_seg *segs;
    size_t nsegs, segcap;
    uint64_t strls;
    int rc;
    char msg[200];
} dta_task;

struct dta_writer
{
    const ctools_dta_spec *spec;
    int release;
    size_t K, W;
    dta_col *cols;
    size_t nproducer, nstrl, nmeasure;
    int producer_strl;   /* some strL column is PRODUCER */
    int threads;
    ctools_persistent_pool *pool;
    size_t maxwave;
    dta_task *sets[2];   /* task slots, alternating between waves */
    dta_file file;
    /* plan results, per plan unit */
    uint64_t *plan_bytes;
    size_t *plan_max; /* [unit * nmeasure + j] */
    /* producer groups kept from the cells pass for the strls pass */
    ctools_dta_group *kept;
    unsigned char *keep, *kept_filled;
    void **group_columns; /* per task slot: nproducer entries */
    uint64_t strls, strl_bytes;
    pthread_mutex_t lock; /* task completion, for the appending thread */
    pthread_cond_t cond;
    int sync;
};

#if defined(__GNUC__) || defined(__clang__)
__attribute__((format(printf, 4, 5)))
#endif
static int dta_fail(char *msg, size_t cap, int rc, const char *fmt, ...)
{
    if (msg && cap && !msg[0]) {
        va_list ap;
        va_start(ap, fmt);
        vsnprintf(msg, cap, fmt, ap);
        va_end(ap);
    }
    return rc;
}

/* ------------------------------------------------------------------------
 * Validation
 * ------------------------------------------------------------------------ */

static const char *const dta_reserved[] = {
    "_all", "_b", "byte", "_coef", "_cons", "double", "float", "if", "in", "int", "long",
    "_n", "_N", "_pi", "_pred", "_rc", "_skip", "strL", "using", "with", NULL};

/* A Stata name: 1-32 characters (128 bytes) of letters, digits and _,
 * not starting with a digit, valid UTF-8, not a reserved word or str#. */
static int dta_name_ok(const char *s)
{
    if (!s) return 0;
    size_t bytes = strlen(s), chars = 0;
    if (!bytes || bytes > 128) return 0;
    for (size_t i = 0; i < bytes;) {
        unsigned char c = (unsigned char)s[i];
        size_t n = c < 0x80 ? 1 : (c >> 5) == 6 ? 2 : (c >> 4) == 14 ? 3 : (c >> 3) == 30 ? 4 : 0;
        if (!n || i + n > bytes) return 0;
        for (size_t k = 1; k < n; k++)
            if (((unsigned char)s[i + k] >> 6) != 2) return 0;
        if (n == 1 && !(c == '_' || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                        (i && c >= '0' && c <= '9')))
            return 0;
        i += n;
        chars++;
    }
    if (chars > 32) return 0;
    for (int i = 0; dta_reserved[i]; i++)
        if (!strcmp(s, dta_reserved[i])) return 0;
    if (!strncmp(s, "str", 3) && s[3]) {
        const char *p = s + 3;
        while (*p >= '0' && *p <= '9') p++;
        if (!*p && s[3] != '0' && atoi(s + 3) >= 1 && atoi(s + 3) <= DTA_MAX_STR) return 0;
    }
    return 1;
}

static int dta_text_ok(const char *s, size_t max) { return !s || strlen(s) <= max; }

static int dta_format_ok(const char *f, int string)
{
    if (!f) return 1;
    size_t n = strlen(f);
    if (n < 2 || n > 56 || f[0] != '%') return 0;
    for (size_t i = 0; i < n; i++)
        if ((unsigned char)f[i] < 0x20 || (unsigned char)f[i] > 0x7e) return 0;
    return string ? f[n - 1] == 's' : f[n - 1] != 's';
}

static int dta_cmp_str(const void *a, const void *b)
{
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

static int dta_validate(dta_writer *w, char *msg, size_t cap)
{
    const ctools_dta_spec *s = w->spec;
    if (!s || (s->ncols && !s->cols)) return dta_fail(msg, cap, 198, "invalid dataset description");
    if (s->release && s->release != 118 && s->release != 119)
        return dta_fail(msg, cap, 198, "release %d is not 118 or 119", s->release);
    w->K = s->ncols;
    w->release = s->release ? s->release : (w->K > 32767 ? 119 : 118);
    if (w->release == 118 && w->K > 32767)
        return dta_fail(msg, cap, 198, "format 118 holds at most 32,767 variables");
    if (w->K > 16777215) return dta_fail(msg, cap, 198, "too many variables");
    if (s->nobs >= (w->release == 118 ? (UINT64_C(1) << 48) : (UINT64_C(1) << 40)))
        return dta_fail(msg, cap, 198, "too many observations for format %d", w->release);
    if (!dta_text_ok(s->data_label, 320)) return dta_fail(msg, cap, 198, "data label too long");
    if (s->timestamp && strlen(s->timestamp) > 17) return dta_fail(msg, cap, 198, "invalid timestamp");

    const ctools_dta_producer *pr = s->producer;
    for (size_t j = 0; j < w->K; j++) {
        const ctools_dta_column *c = &s->cols[j];
        int string = c->type == CTOOLS_DTA_STR || c->type == CTOOLS_DTA_STRL;
        if (!dta_name_ok(c->name)) return dta_fail(msg, cap, 198, "invalid variable name %.130s", c->name ? c->name : "");
        if (c->type < CTOOLS_DTA_BYTE || c->type > CTOOLS_DTA_STRL)
            return dta_fail(msg, cap, 198, "variable %s: invalid storage type", c->name);
        if (string ? c->source == CTOOLS_DTA_SRC_F64
                   : (c->source != CTOOLS_DTA_SRC_F64 && c->source != CTOOLS_DTA_SRC_PRODUCER))
            return dta_fail(msg, cap, 198, "variable %s: source does not match its type", c->name);
        if (c->source < CTOOLS_DTA_SRC_F64 || c->source > CTOOLS_DTA_SRC_PRODUCER)
            return dta_fail(msg, cap, 198, "variable %s: invalid source", c->name);
        if (c->source == CTOOLS_DTA_SRC_PRODUCER) {
            if (!pr) return dta_fail(msg, cap, 198, "variable %s: no producer", c->name);
            if (c->rows) return dta_fail(msg, cap, 198, "variable %s: a producer column has no row map", c->name);
        } else if (s->nobs && !c->data) {
            return dta_fail(msg, cap, 198, "variable %s: no data", c->name);
        }
        if (c->source == CTOOLS_DTA_SRC_FLAT && !c->stride)
            return dta_fail(msg, cap, 198, "variable %s: zero stride", c->name);
        if (c->type == CTOOLS_DTA_STR &&
            (c->width < 0 || c->width > DTA_MAX_STR || (!c->width && c->source == CTOOLS_DTA_SRC_PRODUCER)))
            return dta_fail(msg, cap, 198, "variable %s: invalid str width %d", c->name, c->width);
        if (!dta_format_ok(c->format, string)) return dta_fail(msg, cap, 198, "variable %s: invalid format", c->name);
        if (!dta_text_ok(c->label, 320)) return dta_fail(msg, cap, 198, "variable %s: label too long", c->name);
        if (c->value_label && c->value_label[0] && (string || !dta_name_ok(c->value_label)))
            return dta_fail(msg, cap, 198, "variable %s: invalid value label", c->name);
    }
    if (pr) {
        if (!pr->fill || !pr->release || (pr->ngroups && !pr->starts))
            return dta_fail(msg, cap, 198, "invalid producer");
        uint64_t at = 0;
        for (size_t g = 0; g < pr->ngroups; g++) {
            if (pr->starts[g] != at || pr->starts[g + 1] <= at)
                return dta_fail(msg, cap, 198, "producer groups do not cover the rows in order");
            at = pr->starts[g + 1];
        }
        if (at != s->nobs || (!pr->ngroups && s->nobs)) return dta_fail(msg, cap, 198, "producer groups do not cover the rows");
    }

    /* distinct names */
    if (w->K > 1) {
        const char **names = malloc(w->K * sizeof(*names));
        if (!names) return dta_fail(msg, cap, 920, "out of memory");
        for (size_t j = 0; j < w->K; j++) names[j] = s->cols[j].name;
        qsort(names, w->K, sizeof(*names), dta_cmp_str);
        int dup = 0;
        for (size_t j = 1; j < w->K && !dup; j++)
            if (!strcmp(names[j - 1], names[j])) {
                dta_fail(msg, cap, 198, "variable %s defined twice", names[j]);
                dup = 1;
            }
        free(names);
        if (dup) return 198;
    }
    if (s->nsort > w->K || (s->nsort && !s->sortlist)) return dta_fail(msg, cap, 198, "invalid sort list");
    for (size_t i = 0; i < s->nsort; i++) {
        if (s->sortlist[i] < 1 || (size_t)s->sortlist[i] > w->K) return dta_fail(msg, cap, 198, "invalid sort list");
        for (size_t k = 0; k < i; k++)
            if (s->sortlist[k] == s->sortlist[i]) return dta_fail(msg, cap, 198, "sort list repeats a variable");
    }
    for (size_t i = 0; i < s->nchars; i++) {
        const ctools_dta_char *ch = &s->chars[i];
        int found = ch->owner && !strcmp(ch->owner, "_dta");
        for (size_t j = 0; j < w->K && !found && ch->owner; j++) found = !strcmp(ch->owner, s->cols[j].name);
        if (!found || !ch->name || !ch->name[0] || strlen(ch->name) > 128 || !ch->text ||
            strlen(ch->text) >= UINT32_MAX - 260)
            return dta_fail(msg, cap, 198, "invalid characteristic");
    }
    for (size_t i = 0; i < s->nlabels; i++) {
        const ctools_dta_value_label *l = &s->labels[i];
        if (!dta_name_ok(l->name) || (l->n && (!l->values || !l->texts)) || l->n > INT32_MAX / 8)
            return dta_fail(msg, cap, 198, "invalid value label");
        uint64_t text = 0;
        for (size_t k = 0; k < l->n; k++) {
            if (!l->texts[k]) return dta_fail(msg, cap, 198, "value label %s: missing text", l->name);
            text += strlen(l->texts[k]) + 1;
        }
        if (text > INT32_MAX) return dta_fail(msg, cap, 198, "value label %s too long", l->name);
        for (size_t k = 0; k < i; k++)
            if (!strcmp(s->labels[k].name, l->name)) return dta_fail(msg, cap, 198, "value label %s defined twice", l->name);
    }
    return 0;
}

/* ------------------------------------------------------------------------
 * Values
 * ------------------------------------------------------------------------ */

/* Stata's doubles run from -c(maxdouble) to c(maxdouble), below 2^1023 in
 * magnitude; values beyond them in either direction are stored as missing. */
static inline int dta_valid(double x) { return x < 0x1p1023 && x > -0x1p1023; }

/* The missing value a double stands for: 0 for ., 1..26 for .a .. .z. NaN,
 * infinities, values of -2^1023 or less and other patterns at or above
 * 2^1023 are taken as . */
static inline unsigned dta_mcode(double x)
{
    uint64_t bits, d;
    memcpy(&bits, &x, 8);
    d = bits - DTA_D_MISS;
    if (d <= (UINT64_C(26) << 40) && !(d & ((UINT64_C(1) << 40) - 1))) return (unsigned)(d >> 40);
    return 0;
}

static const char *dta_type_name(const dta_col *dc)
{
    switch (dc->c->type) {
    case CTOOLS_DTA_BYTE: return "byte";
    case CTOOLS_DTA_INT: return "int";
    case CTOOLS_DTA_LONG: return "long";
    case CTOOLS_DTA_FLOAT: return "float";
    case CTOOLS_DTA_DOUBLE: return "double";
    case CTOOLS_DTA_STR: return "str#";
    default: return "strL";
    }
}

/* The string at output row `row` (group row i for a PRODUCER column). */
static inline void dta_string(const dta_col *dc, void *const *gcols, uint64_t row, size_t i,
                              const char **p, size_t *n, int *binary)
{
    const ctools_dta_column *c = dc->c;
    *binary = 0;
    if (c->source == CTOOLS_DTA_SRC_PRODUCER) {
        const ctools_dta_span *sp = (const ctools_dta_span *)gcols[dc->slot] + i;
        *p = sp->ptr;
        *n = sp->ptr ? sp->len : 0;
        *binary = sp->binary;
        return;
    }
    int64_t s = c->rows ? c->rows[row] : (int64_t)row;
    if (s < 0) {
        *p = "";
        *n = 0;
        return;
    }
    switch (c->source) {
    case CTOOLS_DTA_SRC_CSTR: {
        const char *v = ((const char *const *)c->data)[s];
        *p = v ? v : "";
        *n = v ? strlen(v) : 0;
        break;
    }
    case CTOOLS_DTA_SRC_FLAT: {
        const char *v = (const char *)c->data + (size_t)s * c->stride;
        const char *z = memchr(v, 0, c->stride);
        *p = v;
        *n = z ? (size_t)(z - v) : c->stride;
        break;
    }
    default: {
        const ctools_dta_span *sp = (const ctools_dta_span *)c->data + s;
        *p = sp->ptr ? sp->ptr : "";
        *n = sp->ptr ? sp->len : 0;
        *binary = sp->binary;
        break;
    }
    }
}

static inline int dta_defines(const dta_col *dc, uint64_t row)
{
    if (!dc->first) return 1;
    int64_t s = dc->c->rows[row];
    return s < 0 || dc->first[s] == row;
}

static int dta_misfit(dta_task *t, const dta_col *dc, uint64_t row, double x)
{
    t->rc = 459;
    snprintf(t->msg, sizeof(t->msg), "variable %s, observation %llu: %.17g does not fit storage type %s",
             dc->c->name, (unsigned long long)row + 1, x, dta_type_name(dc));
    return 459;
}

/* One numeric column of rows [r0, r0 + rows) into records at out (stride W). */
static int dta_encode_numeric(dta_task *t, const dta_col *dc, void *const *gcols, uint64_t r0, size_t rows,
                              unsigned char *out, size_t W)
{
    const ctools_dta_column *c = dc->c;
    const double *v = c->source == CTOOLS_DTA_SRC_PRODUCER ? (const double *)gcols[dc->slot] : (const double *)c->data;
    const int local = c->source == CTOOLS_DTA_SRC_PRODUCER;
    const int64_t *map = c->rows;
#define DTA_X(i) (local ? v[i] : map ? (map[r0 + (i)] < 0 ? 0x1p1023 : v[map[r0 + (i)]]) : v[r0 + (i)])
    switch (c->type) {
    case CTOOLS_DTA_DOUBLE:
        for (size_t i = 0; i < rows; i++) {
            double x = DTA_X(i);
            if (dta_valid(x)) memcpy(out + i * W, &x, 8);
            else {
                uint64_t b = DTA_D_MISS + ((uint64_t)dta_mcode(x) << 40);
                memcpy(out + i * W, &b, 8);
            }
        }
        break;
    case CTOOLS_DTA_FLOAT:
        for (size_t i = 0; i < rows; i++) {
            double x = DTA_X(i);
            if (dta_valid(x)) {
                float f = (float)x;
                if (!(fabsf(f) < 0x1p127f)) return dta_misfit(t, dc, r0 + i, x);
                memcpy(out + i * W, &f, 4);
            } else {
                uint32_t b = DTA_F_MISS + (uint32_t)dta_mcode(x) * 0x800u;
                memcpy(out + i * W, &b, 4);
            }
        }
        break;
    case CTOOLS_DTA_LONG:
        for (size_t i = 0; i < rows; i++) {
            double x = DTA_X(i);
            int32_t b;
            if (dta_valid(x)) {
                if (!(x >= -2147483647.0 && x <= 2147483620.0) || (double)(int32_t)x != x)
                    return dta_misfit(t, dc, r0 + i, x);
                b = (int32_t)x;
            } else {
                b = DTA_L_MISS + (int32_t)dta_mcode(x);
            }
            memcpy(out + i * W, &b, 4);
        }
        break;
    case CTOOLS_DTA_INT:
        for (size_t i = 0; i < rows; i++) {
            double x = DTA_X(i);
            int16_t b;
            if (dta_valid(x)) {
                if (!(x >= -32767.0 && x <= 32740.0) || (double)(int16_t)x != x) return dta_misfit(t, dc, r0 + i, x);
                b = (int16_t)x;
            } else {
                b = (int16_t)(DTA_I_MISS + (int)dta_mcode(x));
            }
            memcpy(out + i * W, &b, 2);
        }
        break;
    default: /* byte */
        for (size_t i = 0; i < rows; i++) {
            double x = DTA_X(i);
            int8_t b;
            if (dta_valid(x)) {
                if (!(x >= -127.0 && x <= 100.0) || (double)(int8_t)x != x) return dta_misfit(t, dc, r0 + i, x);
                b = (int8_t)x;
            } else {
                b = (int8_t)(DTA_B_MISS + (int)dta_mcode(x));
            }
            out[i * W] = (unsigned char)b;
        }
        break;
    }
#undef DTA_X
    return 0;
}

static int dta_encode_string(dta_task *t, const dta_col *dc, size_t col, void *const *gcols, uint64_t r0,
                             size_t rows, unsigned char *out, size_t W)
{
    const int release = t->w->release;
    for (size_t i = 0; i < rows; i++) {
        const char *p;
        size_t n;
        int binary;
        unsigned char *cell = out + i * W;
        dta_string(dc, gcols, r0 + i, i, &p, &n, &binary);
        if (dc->code == 32768) {
            /* an empty cell (no source row included) stores nothing */
            uint64_t o = n && dc->first ? dc->first[dc->c->rows[r0 + i]] + 1 : r0 + i + 1;
            if (!n) {
                memset(cell, 0, 8);
            } else if (release == 118) {
                dta_le(cell, col + 1, 2);
                dta_le(cell + 2, o, 6);
            } else {
                dta_le(cell, col + 1, 3);
                dta_le(cell + 3, o, 5);
            }
        } else {
            if (n > dc->width) {
                t->rc = 459;
                snprintf(t->msg, sizeof(t->msg), "variable %s, observation %llu: %zu bytes do not fit str%u",
                         dc->c->name, (unsigned long long)(r0 + i + 1), n, (unsigned)dc->width);
                return 459;
            }
            memcpy(cell, p, n);
            memset(cell + n, 0, dc->width - n);
        }
    }
    return 0;
}

/* ------------------------------------------------------------------------
 * Producer groups
 * ------------------------------------------------------------------------ */

static int dta_fill(dta_task *t, int phase, const ctools_dta_group **out)
{
    dta_writer *w = t->w;
    const ctools_dta_producer *pr = w->spec->producer;
    size_t g = t->unit->group;
    *out = NULL;
    if (!pr || !w->nproducer) return 0;
    if (t->kept && w->kept_filled[g]) { /* strls pass of a kept group */
        *out = &w->kept[g];
        return 0;
    }
    ctools_dta_group *dst = t->kept ? &w->kept[g] : &t->g;
    void **cols = t->kept ? w->kept[g].columns : w->group_columns + t->slot * w->nproducer;
    memset(cols, 0, w->nproducer * sizeof(*cols));
    dst->columns = cols;
    dst->state = NULL;
    int rc = pr->fill(pr->ctx, g, phase, dst);
    if (rc) {
        snprintf(t->msg, sizeof(t->msg), "producer failed with r(%d)", rc);
        t->rc = rc;
        return rc;
    }
    for (size_t j = 0; j < w->K; j++) {
        const dta_col *dc = &w->cols[j];
        int needed = phase == CTOOLS_DTA_PHASE_CELLS || dc->code == 32768;
        if (dc->c->source == CTOOLS_DTA_SRC_PRODUCER && needed && !cols[dc->slot] &&
            t->unit->r1 > t->unit->r0) {
            pr->release(pr->ctx, dst);
            t->rc = 198;
            snprintf(t->msg, sizeof(t->msg), "producer left variable %s empty", dc->c->name);
            return 198;
        }
    }
    if (t->kept) w->kept_filled[g] = 1;
    else t->filled = 1;
    *out = dst;
    return 0;
}

static void dta_release_task(dta_task *t)
{
    dta_writer *w = t->w;
    const ctools_dta_producer *pr = w->spec->producer;
    if (t->filled) {
        pr->release(pr->ctx, &t->g);
        t->filled = 0;
    }
    if (t->kept && t->phase == CTOOLS_DTA_PHASE_STRLS && w->kept_filled[t->unit->group]) {
        pr->release(pr->ctx, &w->kept[t->unit->group]);
        w->kept_filled[t->unit->group] = 0;
    }
}

static int dta_reserve(dta_task *t, size_t need)
{
    if (need <= t->cap) return 0;
    size_t cap = t->cap ? t->cap : 4096;
    while (cap < need) cap = cap > SIZE_MAX / 2 ? need : cap * 2;
    unsigned char *p = realloc(t->buf, cap);
    if (!p) {
        t->rc = 920;
        snprintf(t->msg, sizeof(t->msg), "out of memory");
        return 920;
    }
    t->buf = p;
    t->cap = cap;
    return 0;
}

static int dta_seg_push(dta_task *t, size_t at, size_t len, const char *ext)
{
    if (!len) return 0;
    if (!ext && t->nsegs && !t->segs[t->nsegs - 1].ext &&
        t->segs[t->nsegs - 1].at + t->segs[t->nsegs - 1].len == at) {
        t->segs[t->nsegs - 1].len += len;
        return 0;
    }
    if (t->nsegs == t->segcap) {
        size_t cap = t->segcap ? 2 * t->segcap : 16;
        dta_seg *s = realloc(t->segs, cap * sizeof(*s));
        if (!s) {
            t->rc = 920;
            snprintf(t->msg, sizeof(t->msg), "out of memory");
            return 920;
        }
        t->segs = s;
        t->segcap = cap;
    }
    t->segs[t->nsegs++] = (dta_seg){at, len, ext};
    return 0;
}

/* ------------------------------------------------------------------------
 * Workers
 * ------------------------------------------------------------------------ */

static void *dta_plan_worker(void *arg)
{
    dta_task *t = arg;
    dta_writer *w = t->w;
    uint64_t bytes = 0;
    size_t *max = w->plan_max ? w->plan_max + t->index * w->nmeasure : NULL;
    size_t m = 0;
    for (size_t j = 0; j < w->K; j++) {
        const dta_col *dc = &w->cols[j];
        int strl = dc->code == 32768;
        if (dc->c->source == CTOOLS_DTA_SRC_PRODUCER || (!strl && !dc->measure)) continue;
        size_t longest = 0;
        for (uint64_t r = t->unit->r0; r < t->unit->r1; r++) {
            const char *p;
            size_t n;
            int binary;
            dta_string(dc, NULL, r, 0, &p, &n, &binary);
            if (strl) {
                if (n && dta_defines(dc, r)) bytes += DTA_GSO_HEAD + 1 + (n < DTA_EXTERNAL_BYTES ? n : 0);
            } else if (n > longest) {
                longest = n;
            }
        }
        if (dc->measure) max[m++] = longest;
    }
    w->plan_bytes[t->index] = bytes;
    return NULL;
}

static void *dta_cells_worker(void *arg)
{
    dta_task *t = arg;
    dta_writer *w = t->w;
    const ctools_dta_group *g;
    size_t rows = (size_t)(t->unit->r1 - t->unit->r0);
    if (dta_fill(t, CTOOLS_DTA_PHASE_CELLS, &g)) return NULL;
    if (dta_reserve(t, rows * w->W)) goto done;
    void *const *gcols = g ? g->columns : NULL;
    for (size_t j = 0; j < w->K && !t->rc; j++) {
        const dta_col *dc = &w->cols[j];
        if (dc->code > 32768) dta_encode_numeric(t, dc, gcols, t->unit->r0, rows, t->buf + dc->offset, w->W);
        else dta_encode_string(t, dc, j, gcols, t->unit->r0, rows, t->buf + dc->offset, w->W);
    }
    t->len = rows * w->W;
done:
    /* the records hold copies; only a kept group lives on for the strls */
    if (t->filled) {
        w->spec->producer->release(w->spec->producer->ctx, &t->g);
        t->filled = 0;
    }
    return NULL;
}

static void *dta_strls_worker(void *arg)
{
    dta_task *t = arg;
    dta_writer *w = t->w;
    const ctools_dta_group *g;
    const uint64_t r0 = t->unit->r0;
    const size_t rows = (size_t)(t->unit->r1 - r0);
    if (w->producer_strl && dta_fill(t, CTOOLS_DTA_PHASE_STRLS, &g)) return NULL;
    if (!w->producer_strl) g = NULL;
    void *const *gcols = g ? g->columns : NULL;

    /* bytes copied: headers, short payloads, terminators */
    size_t need = 0;
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < w->K; j++) {
            const dta_col *dc = &w->cols[j];
            if (dc->code != 32768) continue;
            const char *p;
            size_t n;
            int binary;
            dta_string(dc, gcols, r0 + i, i, &p, &n, &binary);
            if (!n || !dta_defines(dc, r0 + i)) continue;
            if (n > DTA_MAX_STRL) {
                t->rc = 459;
                snprintf(t->msg, sizeof(t->msg), "variable %s, observation %llu: strL longer than 2,000,000,000 bytes",
                         dc->c->name, (unsigned long long)(r0 + i + 1));
                return NULL;
            }
            need += DTA_GSO_HEAD + 1 + (n < DTA_EXTERNAL_BYTES ? n : 0);
        }
    if (dta_reserve(t, need ? need : 1)) return NULL;

    size_t len = 0, start = 0;
    for (size_t i = 0; i < rows; i++) {
        const uint64_t o = r0 + i + 1;
        for (size_t j = 0; j < w->K; j++) {
            const dta_col *dc = &w->cols[j];
            if (dc->code != 32768) continue;
            const char *p;
            size_t n;
            int binary;
            dta_string(dc, gcols, r0 + i, i, &p, &n, &binary);
            if (!n || !dta_defines(dc, r0 + i)) continue;
            if (!binary && !(dc->c->flags & CTOOLS_DTA_COL_TEXT) && memchr(p, 0, n)) binary = 1;
            unsigned char *h = t->buf + len;
            memcpy(h, "GSO", 3);
            dta_le(h + 3, j + 1, 4);
            dta_le(h + 7, o, 8);
            h[15] = binary ? 129 : 130; /* 130: text stored with its terminating NUL */
            dta_le(h + 16, binary ? n : n + 1, 4);
            len += DTA_GSO_HEAD;
            if (n < DTA_EXTERNAL_BYTES) {
                memcpy(t->buf + len, p, n);
                len += n;
            } else {
                if (dta_seg_push(t, start, len - start, NULL) || dta_seg_push(t, 0, n, p)) return NULL;
                start = len;
            }
            if (!binary) t->buf[len++] = 0;
            t->strls++;
        }
    }
    if (dta_seg_push(t, start, len - start, NULL)) return NULL;
    t->len = len;
    return NULL;
}

/* ------------------------------------------------------------------------
 * Waves
 * ------------------------------------------------------------------------ */

/* Append a built unit to the file. */
static int dta_write_task(dta_writer *w, int phase, dta_task *t)
{
    if (phase == CTOOLS_DTA_PHASE_CELLS) return dta_file_write(&w->file, t->buf, t->len) ? 693 : 0;
    for (size_t k = 0; k < t->nsegs; k++) {
        const dta_seg *s = &t->segs[k];
        if (dta_file_write(&w->file, s->ext ? (const void *)s->ext : t->buf + s->at, s->len)) return 693;
        w->strl_bytes += s->len;
    }
    w->strls += t->strls;
    return 0;
}

/* Wait until a worker has built t (each signals as it finishes). */
static void dta_wait_task(dta_writer *w, dta_task *t)
{
    pthread_mutex_lock(&w->lock);
    while (!t->done) pthread_cond_wait(&w->cond, &w->lock);
    pthread_mutex_unlock(&w->lock);
}

/* The appending thread of a pass takes the waves in order and writes each
 * unit as soon as it is built. At most two waves are pushed and not yet
 * written, one per set of task slots. After a failed write, or at a unit
 * that failed, it stops writing but still marks waves written, so the pass
 * never waits on it forever. Shared fields are guarded by w->lock. */
typedef struct
{
    dta_writer *w;
    int phase;
    dta_task *tasks[2];
    size_t n[2];
    size_t pushed, written; /* waves */
    int closed;
    int rc;
} dta_appender;

static void *dta_appender_main(void *arg)
{
    dta_appender *a = arg;
    dta_writer *w = a->w;
    int stop = 0;
    for (size_t k = 0;; k++) {
        pthread_mutex_lock(&w->lock);
        while (a->pushed <= k && !a->closed) pthread_cond_wait(&w->cond, &w->lock);
        int more = a->pushed > k;
        dta_task *tasks = a->tasks[k & 1];
        size_t n = a->n[k & 1];
        pthread_mutex_unlock(&w->lock);
        if (!more) break;
        for (size_t i = 0; i < n; i++) {
            dta_task *t = &tasks[i];
            dta_wait_task(w, t);
            if (stop || t->rc) {
                stop = 1;
                continue;
            }
            int rc = dta_write_task(w, a->phase, t);
            if (rc) {
                pthread_mutex_lock(&w->lock);
                a->rc = rc;
                pthread_mutex_unlock(&w->lock);
                stop = 1;
            }
        }
        pthread_mutex_lock(&w->lock);
        a->written = k + 1;
        pthread_cond_broadcast(&w->cond);
        pthread_mutex_unlock(&w->lock);
    }
    return NULL;
}

static void *dta_task_main(void *arg)
{
    dta_task *t = arg;
    t->fn(t);
    pthread_mutex_lock(&t->w->lock);
    t->done = 1;
    pthread_cond_broadcast(&t->w->cond);
    pthread_mutex_unlock(&t->w->lock);
    return NULL;
}

/* Start building tasks: 1 when they run on the pool (finish with
 * ctools_persistent_pool_wait), 0 when they were built here. */
static int dta_start(dta_writer *w, dta_task *tasks, size_t n)
{
    if (n > 1 && w->pool &&
        ctools_persistent_pool_submit_batch(w->pool, dta_task_main, tasks, n, sizeof(*tasks)) == 0)
        return 1;
    for (size_t i = 0; i < n; i++) dta_task_main(&tasks[i]);
    return 0;
}

static void dta_finish_wave(dta_task *tasks, size_t n)
{
    for (size_t i = 0; tasks && i < n; i++) dta_release_task(&tasks[i]);
}

/* Build units with fn a wave at a time, at most cap_units per wave. With
 * io, one thread appends each unit in order as soon as it is built, while
 * the pool builds the next wave. Returns the first failure in unit order. */
static int dta_pass(dta_writer *w, int phase, void *(*fn)(void *), const dta_unit *units, size_t nunits,
                    uint64_t budget, size_t cap_units, int io, char *msg, size_t cap)
{
    int rc = 0, running = 0;
    dta_appender a;
    memset(&a, 0, sizeof(a));
    a.w = w;
    a.phase = phase;
    pthread_t thread;
    if (io && w->pool && nunits > 1) {
        pthread_attr_t attr;
        if (ctools_worker_attr_init(&attr) == 0) {
            running = pthread_create(&thread, &attr, dta_appender_main, &a) == 0;
            pthread_attr_destroy(&attr);
        }
    }
    size_t built[2] = {0, 0}; /* units of the wave last built in each set */
    if (cap_units < 1) cap_units = 1;
    if (cap_units > w->maxwave) cap_units = w->maxwave;
    for (size_t u = 0, wave = 0; u < nunits && !rc; wave++) {
        dta_task *set = w->sets[wave & 1];
        if (wave >= 2) { /* this set's previous wave must be out first */
            if (running) {
                pthread_mutex_lock(&w->lock);
                while (a.written < wave - 1) pthread_cond_wait(&w->cond, &w->lock);
                int failed = a.rc;
                pthread_mutex_unlock(&w->lock);
                if (failed) rc = dta_fail(msg, cap, failed, "cannot write the file");
            }
            dta_finish_wave(set, built[wave & 1]);
            built[wave & 1] = 0;
            if (rc) break;
        }
        size_t n = 0;
        uint64_t bytes = 0;
        while (u + n < nunits && n < cap_units && (n == 0 || bytes + units[u + n].bytes <= budget))
            bytes += units[u + n++].bytes;
        for (size_t i = 0; i < n; i++) {
            dta_task *t = &set[i];
            t->w = w;
            t->fn = fn;
            t->unit = &units[u + i];
            t->index = u + i;
            t->phase = phase;
            t->filled = 0;
            t->kept = w->keep && w->keep[units[u + i].group];
            t->len = 0;
            t->nsegs = 0;
            t->strls = 0;
            t->rc = 0;
            t->done = 0;
            t->msg[0] = 0;
        }
        built[wave & 1] = n;
        int pooled = dta_start(w, set, n);
        if (running) {
            pthread_mutex_lock(&w->lock);
            a.tasks[wave & 1] = set;
            a.n[wave & 1] = n;
            a.pushed = wave + 1;
            pthread_cond_broadcast(&w->cond);
            pthread_mutex_unlock(&w->lock);
        }
        if (pooled) ctools_persistent_pool_wait(w->pool);
        for (size_t i = 0; i < n && !rc; i++)
            if (set[i].rc) rc = dta_fail(msg, cap, set[i].rc, "%s", set[i].msg);
        if (io && !running)
            for (size_t i = 0; i < n && !rc; i++)
                if (dta_write_task(w, phase, &set[i])) rc = dta_fail(msg, cap, 693, "cannot write the file");
        u += n;
    }
    if (running) {
        pthread_mutex_lock(&w->lock);
        a.closed = 1;
        pthread_cond_broadcast(&w->cond);
        pthread_mutex_unlock(&w->lock);
        pthread_join(thread, NULL);
        if (!rc && a.rc) rc = dta_fail(msg, cap, a.rc, "cannot write the file");
    }
    dta_finish_wave(w->sets[0], built[0]);
    dta_finish_wave(w->sets[1], built[1]);
    return rc;
}

/* ------------------------------------------------------------------------
 * Header and trailer
 * ------------------------------------------------------------------------ */

static const char *dta_default_format(const dta_col *dc, char *buf, size_t cap)
{
    switch (dc->c->type) {
    case CTOOLS_DTA_BYTE:
    case CTOOLS_DTA_INT: return "%8.0g";
    case CTOOLS_DTA_LONG: return "%12.0g";
    case CTOOLS_DTA_FLOAT: return "%9.0g";
    case CTOOLS_DTA_DOUBLE: return "%10.0g";
    case CTOOLS_DTA_STR:
        snprintf(buf, cap, "%%%us", (unsigned)(dc->width < 9 ? 9 : dc->width));
        return buf;
    default: return "%9s";
    }
}

static void dta_header(dta_writer *w, dta_bytes *b, uint64_t *map, size_t *map_at)
{
    const ctools_dta_spec *s = w->spec;
    const int k = w->release == 118 ? 2 : 4;
    char fmt[16];
    map[0] = 0;
    dta_tag(b, w->release == 118 ? "<stata_dta><header><release>118</release><byteorder>LSF</byteorder><K>"
                                 : "<stata_dta><header><release>119</release><byteorder>LSF</byteorder><K>");
    dta_put_le(b, w->K, k);
    dta_tag(b, "</K><N>");
    dta_put_le(b, s->nobs, 8);
    dta_tag(b, "</N><label>");
    size_t label = s->data_label ? strlen(s->data_label) : 0;
    dta_put_le(b, label, 2);
    dta_put(b, s->data_label, label);
    dta_tag(b, "</label><timestamp>");
    size_t stamp = s->timestamp ? strlen(s->timestamp) : 0;
    dta_put_le(b, stamp, 1);
    dta_put(b, s->timestamp, stamp);
    dta_tag(b, "</timestamp></header>");
    map[1] = b->len;
    dta_tag(b, "<map>");
    *map_at = b->len;
    dta_put(b, NULL, 14 * 8);
    dta_tag(b, "</map>");
    map[2] = b->len;
    dta_tag(b, "<variable_types>");
    for (size_t j = 0; j < w->K; j++) dta_put_le(b, w->cols[j].code, 2);
    dta_tag(b, "</variable_types>");
    map[3] = b->len;
    dta_tag(b, "<varnames>");
    for (size_t j = 0; j < w->K; j++) dta_put_fixed(b, w->cols[j].c->name, 129);
    dta_tag(b, "</varnames>");
    map[4] = b->len;
    dta_tag(b, "<sortlist>");
    for (size_t i = 0; i <= w->K; i++) dta_put_le(b, i < s->nsort ? (uint64_t)s->sortlist[i] : 0, k);
    dta_tag(b, "</sortlist>");
    map[5] = b->len;
    dta_tag(b, "<formats>");
    for (size_t j = 0; j < w->K; j++) {
        const char *f = w->cols[j].c->format;
        dta_put_fixed(b, f ? f : dta_default_format(&w->cols[j], fmt, sizeof(fmt)), 57);
    }
    dta_tag(b, "</formats>");
    map[6] = b->len;
    dta_tag(b, "<value_label_names>");
    for (size_t j = 0; j < w->K; j++) dta_put_fixed(b, w->cols[j].c->value_label, 129);
    dta_tag(b, "</value_label_names>");
    map[7] = b->len;
    dta_tag(b, "<variable_labels>");
    for (size_t j = 0; j < w->K; j++) dta_put_fixed(b, w->cols[j].c->label, 321);
    dta_tag(b, "</variable_labels>");
    map[8] = b->len;
    dta_tag(b, "<characteristics>");
    for (size_t i = 0; i < s->nchars; i++) {
        const ctools_dta_char *ch = &s->chars[i];
        size_t text = strlen(ch->text) + 1;
        dta_tag(b, "<ch>");
        dta_put_le(b, 129 + 129 + text, 4);
        dta_put_fixed(b, ch->owner, 129);
        dta_put_fixed(b, ch->name, 129);
        dta_put(b, ch->text, text);
        dta_tag(b, "</ch>");
    }
    dta_tag(b, "</characteristics>");
    map[9] = b->len;
    dta_tag(b, "<data>");
}

typedef struct
{
    int32_t value;
    const char *text;
} dta_entry;

static int dta_cmp_entry(const void *a, const void *b)
{
    int32_t x = ((const dta_entry *)a)->value, y = ((const dta_entry *)b)->value;
    return (x > y) - (x < y);
}

static int dta_value_labels(dta_writer *w, dta_bytes *b, char *msg, size_t cap)
{
    const ctools_dta_spec *s = w->spec;
    dta_tag(b, "<value_labels>");
    for (size_t i = 0; i < s->nlabels; i++) {
        const ctools_dta_value_label *l = &s->labels[i];
        dta_entry *e = l->n ? malloc(l->n * sizeof(*e)) : NULL;
        if (l->n && !e) return dta_fail(msg, cap, 920, "out of memory");
        size_t text = 0;
        for (size_t k = 0; k < l->n; k++) {
            e[k] = (dta_entry){l->values[k], l->texts[k]};
            text += strlen(l->texts[k]) + 1;
        }
        if (l->n) qsort(e, l->n, sizeof(*e), dta_cmp_entry);
        for (size_t k = 1; k < l->n; k++)
            if (e[k].value == e[k - 1].value) {
                free(e);
                return dta_fail(msg, cap, 198, "value label %s labels %d twice", l->name, (int)e[k].value);
            }
        dta_tag(b, "<lbl>");
        dta_put_le(b, 8 + 8 * (uint64_t)l->n + text, 4);
        dta_put_fixed(b, l->name, 129);
        dta_put(b, NULL, 3);
        dta_put_le(b, l->n, 4);
        dta_put_le(b, text, 4);
        size_t at = 0;
        for (size_t k = 0; k < l->n; k++) {
            dta_put_le(b, at, 4);
            at += strlen(e[k].text) + 1;
        }
        for (size_t k = 0; k < l->n; k++) dta_put_le(b, (uint32_t)e[k].value, 4);
        for (size_t k = 0; k < l->n; k++) dta_put(b, e[k].text, strlen(e[k].text) + 1);
        dta_tag(b, "</lbl>");
        free(e);
    }
    dta_tag(b, "</value_labels>");
    return 0;
}

/* ------------------------------------------------------------------------
 * Entry point
 * ------------------------------------------------------------------------ */

void ctools_dta_spec_init(ctools_dta_spec *spec)
{
    if (spec) memset(spec, 0, sizeof(*spec));
}

static void dta_writer_free(dta_writer *w)
{
    const ctools_dta_producer *pr = w->spec ? w->spec->producer : NULL;
    if (w->kept && w->kept_filled && pr)
        for (size_t g = 0; g < pr->ngroups; g++)
            if (w->kept_filled[g]) pr->release(pr->ctx, &w->kept[g]);
    for (int s = 0; s < 2; s++) {
        if (!w->sets[s]) continue;
        for (size_t i = 0; i < w->maxwave; i++) {
            free(w->sets[s][i].buf);
            free(w->sets[s][i].segs);
        }
        free(w->sets[s]);
    }
    if (w->kept) {
        for (size_t g = 0; pr && g < pr->ngroups; g++) free(w->kept[g].columns);
        free(w->kept);
    }
    free(w->keep);
    free(w->kept_filled);
    free(w->group_columns);
    for (size_t j = 0; w->cols && j < w->K; j++) free(w->cols[j].first);
    free(w->cols);
    free(w->plan_bytes);
    free(w->plan_max);
    if (w->sync) {
        pthread_cond_destroy(&w->cond);
        pthread_mutex_destroy(&w->lock);
    }
}

int ctools_dta_write(const char *path, const ctools_dta_spec *spec, ctools_dta_result *result)
{
    ctools_dta_result local;
    if (!result) result = &local;
    memset(result, 0, sizeof(*result));
    char *msg = result->message;
    const size_t cap = sizeof(result->message);
    const uint16_t one = 1;
    if (*(const unsigned char *)&one != 1) return dta_fail(msg, cap, 198, "big-endian hosts are not supported");
    if (!path || !*path) return dta_fail(msg, cap, 198, "no file name");
    const double started = ctools_timer_seconds();

    dta_writer w;
    memset(&w, 0, sizeof(w));
    w.spec = spec;
    int rc = dta_validate(&w, msg, cap);
    dta_unit *units = NULL, *gso_units = NULL, *plan_units = NULL;
    dta_bytes b = {0};
    int opened = 0;
    if (rc) goto done;
    const ctools_dta_producer *pr = spec->producer;
    const uint64_t N = spec->nobs;

    /* columns */
    w.cols = calloc(w.K ? w.K : 1, sizeof(*w.cols));
    if (!w.cols) {
        rc = dta_fail(msg, cap, 920, "out of memory");
        goto done;
    }
    for (size_t j = 0; j < w.K; j++) {
        dta_col *dc = &w.cols[j];
        dc->c = &spec->cols[j];
        switch (dc->c->type) {
        case CTOOLS_DTA_BYTE: dc->code = 65530, dc->width = 1; break;
        case CTOOLS_DTA_INT: dc->code = 65529, dc->width = 2; break;
        case CTOOLS_DTA_LONG: dc->code = 65528, dc->width = 4; break;
        case CTOOLS_DTA_FLOAT: dc->code = 65527, dc->width = 4; break;
        case CTOOLS_DTA_DOUBLE: dc->code = 65526, dc->width = 8; break;
        case CTOOLS_DTA_STRL: dc->code = 32768, dc->width = 8; break;
        default:
            dc->width = (uint32_t)dc->c->width;
            dc->measure = !dc->width;
            if (dc->measure) w.nmeasure++;
            break;
        }
        if (dc->c->source == CTOOLS_DTA_SRC_PRODUCER) {
            dc->slot = w.nproducer++;
            if (dc->code == 32768) w.producer_strl = 1;
        }
        if (dc->code == 32768) w.nstrl++;
    }
    w.threads = spec->threads > 0 ? spec->threads : ctools_get_max_threads();
    if (w.threads < 1) w.threads = 1;
    w.pool = w.threads > 1 ? ctools_get_global_pool() : NULL;
    w.maxwave = 2 * (size_t)w.threads;
    if (w.maxwave < 2) w.maxwave = 2;
    if (pthread_mutex_init(&w.lock, NULL)) {
        rc = dta_fail(msg, cap, 920, "cannot create a lock");
        goto done;
    }
    if (pthread_cond_init(&w.cond, NULL)) {
        pthread_mutex_destroy(&w.lock);
        rc = dta_fail(msg, cap, 920, "cannot create a lock");
        goto done;
    }
    w.sync = 1;
    for (int s = 0; s < 2; s++) {
        if (!(w.sets[s] = calloc(w.maxwave, sizeof(dta_task)))) {
            rc = dta_fail(msg, cap, 920, "out of memory");
            goto done;
        }
        for (size_t i = 0; i < w.maxwave; i++) w.sets[s][i].slot = (size_t)s * w.maxwave + i;
    }
    if (w.nproducer && !(w.group_columns = calloc(2 * w.maxwave * w.nproducer, sizeof(void *)))) {
        rc = dta_fail(msg, cap, 920, "out of memory");
        goto done;
    }

    /* shared strL columns: the first output row reading each source row */
    for (size_t j = 0; j < w.K && N; j++) {
        dta_col *dc = &w.cols[j];
        const ctools_dta_column *c = dc->c;
        if (!(c->flags & CTOOLS_DTA_COL_SHARE) || dc->code != 32768 || !c->rows ||
            c->source == CTOOLS_DTA_SRC_PRODUCER)
            continue;
        int64_t top = -1;
        for (uint64_t r = 0; r < N; r++)
            if (c->rows[r] > top) top = c->rows[r];
        if (top < 0) continue;
        dc->first = malloc(((size_t)top + 1) * sizeof(uint64_t));
        if (!dc->first) {
            rc = dta_fail(msg, cap, 920, "out of memory");
            goto done;
        }
        memset(dc->first, 0xff, ((size_t)top + 1) * sizeof(uint64_t));
        for (uint64_t r = 0; r < N; r++) {
            int64_t src = c->rows[r];
            if (src >= 0 && dc->first[src] == UINT64_MAX) dc->first[src] = r;
        }
    }

    /* plan: longest str# values and strL bytes, for columns held in memory */
    int need_plan = w.nmeasure > 0 || (w.nstrl && !w.producer_strl);
    size_t nplan = need_plan && N ? (size_t)((N + DTA_PLAN_ROWS - 1) / DTA_PLAN_ROWS) : 0;
    if (nplan) {
        plan_units = calloc(nplan, sizeof(*plan_units));
        w.plan_bytes = calloc(nplan, sizeof(*w.plan_bytes));
        w.plan_max = w.nmeasure ? calloc(nplan * w.nmeasure, sizeof(*w.plan_max)) : NULL;
        if (!plan_units || !w.plan_bytes || (w.nmeasure && !w.plan_max)) {
            rc = dta_fail(msg, cap, 920, "out of memory");
            goto done;
        }
        for (size_t u = 0; u < nplan; u++) {
            plan_units[u].r0 = (uint64_t)u * DTA_PLAN_ROWS;
            plan_units[u].r1 = plan_units[u].r0 + DTA_PLAN_ROWS < N ? plan_units[u].r0 + DTA_PLAN_ROWS : N;
        }
        rc = dta_pass(&w, CTOOLS_DTA_PHASE_CELLS, dta_plan_worker, plan_units, nplan, UINT64_MAX, w.maxwave, 0,
                      msg, cap);
        if (rc) goto done;
    }
    for (size_t j = 0, m = 0; j < w.K; j++) {
        dta_col *dc = &w.cols[j];
        if (!dc->measure) continue;
        size_t longest = 0;
        for (size_t u = 0; u < nplan; u++)
            if (w.plan_max[u * w.nmeasure + m] > longest) longest = w.plan_max[u * w.nmeasure + m];
        if (longest > DTA_MAX_STR) {
            rc = dta_fail(msg, cap, 459, "variable %s: a value of %zu bytes does not fit str#; use strL",
                          dc->c->name, longest);
            goto done;
        }
        dc->width = (uint32_t)(longest ? longest : 1);
        dc->code = (uint16_t)dc->width;
        m++;
    }
    for (size_t j = 0; j < w.K; j++) {
        if (w.cols[j].c->type == CTOOLS_DTA_STR) w.cols[j].code = (uint16_t)w.cols[j].width;
        w.cols[j].offset = w.W;
        w.W += w.cols[j].width;
    }

    /* units: data blocks (producer groups when there is a producer) */
    size_t nunits = 0;
    if (pr) {
        nunits = pr->ngroups;
        units = calloc(nunits ? nunits : 1, sizeof(*units));
        if (!units) {
            rc = dta_fail(msg, cap, 920, "out of memory");
            goto done;
        }
        for (size_t g = 0; g < nunits; g++)
            units[g] = (dta_unit){pr->starts[g], pr->starts[g + 1], g, pr->bytes ? pr->bytes[g] : 0};
    } else if (N && w.W) {
        uint64_t R = DTA_BLOCK_BYTES / w.W ? DTA_BLOCK_BYTES / w.W : 1;
        nunits = (size_t)((N + R - 1) / R);
        units = calloc(nunits, sizeof(*units));
        if (!units) {
            rc = dta_fail(msg, cap, 920, "out of memory");
            goto done;
        }
        for (size_t u = 0; u < nunits; u++) {
            uint64_t r0 = (uint64_t)u * R, r1 = r0 + R < N ? r0 + R : N;
            units[u] = (dta_unit){r0, r1, 0, (r1 - r0) * w.W};
        }
    }
    /* producer groups kept from the cells pass, so their strLs need no refill */
    if (pr && w.producer_strl && pr->bytes && pr->ngroups) {
        size_t budget = spec->keep_bytes ? spec->keep_bytes : DTA_KEEP_BYTES, used = 0;
        w.keep = calloc(pr->ngroups, 1);
        w.kept_filled = calloc(pr->ngroups, 1);
        w.kept = calloc(pr->ngroups, sizeof(*w.kept));
        if (!w.keep || !w.kept_filled || !w.kept) {
            rc = dta_fail(msg, cap, 920, "out of memory");
            goto done;
        }
        for (size_t g = 0; g < pr->ngroups; g++) {
            if (pr->bytes[g] > budget - used) continue;
            if (!(w.kept[g].columns = calloc(w.nproducer, sizeof(void *)))) {
                rc = dta_fail(msg, cap, 920, "out of memory");
                goto done;
            }
            w.keep[g] = 1;
            used += (size_t)pr->bytes[g];
        }
    }

    /* header */
    uint64_t map[14] = {0};
    size_t map_at = 0;
    dta_header(&w, &b, map, &map_at);
    if (b.failed) {
        rc = dta_fail(msg, cap, 920, "out of memory");
        goto done;
    }
    if (dta_file_open(&w.file, path, spec->temporary)) {
        rc = dta_fail(msg, cap, 603, "cannot create %.200s", path);
        goto done;
    }
    opened = 1;
    if (dta_file_write(&w.file, b.p, b.len)) {
        rc = dta_fail(msg, cap, 693, "cannot write the file");
        goto done;
    }

    /* data */
    int gathers = 0;
    for (size_t j = 0; j < w.K; j++) gathers |= spec->cols[j].rows != NULL;
    const size_t wave_units = pr || gathers ? w.maxwave : DTA_IO_WAVE_UNITS;
    rc = dta_pass(&w, CTOOLS_DTA_PHASE_CELLS, dta_cells_worker, units, nunits,
                  pr ? DTA_PRODUCER_WAVE_BYTES : DTA_WAVE_BYTES, wave_units, 1, msg, cap);
    if (rc) goto done;
    if (w.file.pos != map[9] + 6 + N * w.W) {
        rc = dta_fail(msg, cap, 459, "data section has the wrong size");
        goto done;
    }
    static const char mid[] = "</data><strls>";
    if (dta_file_write(&w.file, mid, sizeof(mid) - 1)) {
        rc = dta_fail(msg, cap, 693, "cannot write the file");
        goto done;
    }
    map[10] = map[9] + 6 + N * w.W + 7;

    /* strLs: producer groups, or plan ranges cut to about DTA_GSO_UNIT_BYTES */
    if (w.nstrl && N) {
        const dta_unit *su = units;
        size_t nsu = nunits;
        if (!pr) {
            size_t count = 0;
            for (size_t u = 0; u < nplan; u++) {
                uint64_t rows = plan_units[u].r1 - plan_units[u].r0;
                uint64_t k = (w.plan_bytes[u] + DTA_GSO_UNIT_BYTES - 1) / DTA_GSO_UNIT_BYTES;
                count += (size_t)(k < 1 ? 1 : k > rows ? rows : k);
            }
            gso_units = calloc(count ? count : 1, sizeof(*gso_units));
            if (!gso_units) {
                rc = dta_fail(msg, cap, 920, "out of memory");
                goto done;
            }
            nsu = 0;
            for (size_t u = 0; u < nplan; u++) {
                uint64_t r0 = plan_units[u].r0, rows = plan_units[u].r1 - r0;
                uint64_t k = (w.plan_bytes[u] + DTA_GSO_UNIT_BYTES - 1) / DTA_GSO_UNIT_BYTES;
                k = k < 1 ? 1 : k > rows ? rows : k;
                for (uint64_t i = 0; i < k; i++)
                    gso_units[nsu++] = (dta_unit){r0 + rows * i / k, r0 + rows * (i + 1) / k, 0, w.plan_bytes[u] / k};
            }
            su = gso_units;
        }
        rc = dta_pass(&w, CTOOLS_DTA_PHASE_STRLS, dta_strls_worker, su, nsu,
                      pr ? DTA_PRODUCER_WAVE_BYTES : DTA_WAVE_BYTES, wave_units, 1, msg, cap);
        if (rc) goto done;
    }
    result->strls = w.strls;

    /* trailer and map */
    b.len = 0;
    dta_tag(&b, "</strls>");
    size_t labels_at = b.len;
    if ((rc = dta_value_labels(&w, &b, msg, cap)) != 0) goto done;
    size_t end_at = b.len;
    dta_tag(&b, "</stata_dta>");
    if (b.failed) {
        rc = dta_fail(msg, cap, 920, "out of memory");
        goto done;
    }
    uint64_t trailer_at = w.file.pos;
    if (trailer_at != map[10] + 7 + w.strl_bytes) {
        rc = dta_fail(msg, cap, 459, "strL section has the wrong size");
        goto done;
    }
    map[11] = trailer_at + labels_at;
    map[12] = trailer_at + end_at;
    map[13] = trailer_at + b.len;
    unsigned char m[14 * 8];
    for (int i = 0; i < 14; i++) dta_le(m + 8 * i, map[i], 8);
    if (dta_file_write(&w.file, b.p, b.len) || dta_file_patch(&w.file, m, sizeof(m), map_at)) {
        rc = dta_fail(msg, cap, 693, "cannot write the file");
        goto done;
    }
    result->file_bytes = w.file.pos;
    result->write_seconds = w.file.busy;

done:
    if (opened && dta_file_close(&w.file) && !rc) rc = dta_fail(msg, cap, 693, "cannot close the file");
    if (opened && rc) ctools_remove(path);
    dta_writer_free(&w);
    free(units);
    free(gso_units);
    free(plan_units);
    free(b.p);
    result->seconds = ctools_timer_seconds() - started;
    return rc;
}

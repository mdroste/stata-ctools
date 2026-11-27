/*
 * csql_dta_scan.c - Stata .dta files as a DuckDB table function.
 *
 *   SELECT ... FROM read_dta('file.dta')
 *   SELECT ... FROM 'file.dta'              (replacement scan)
 *   read_dta('file.dta', decode := true)    (value-labeled columns as text)
 *
 * Formats 113 to 119 (Stata 8 and later), either byte order. The file is
 * mapped into memory and read in place: only the columns a query uses are
 * decoded, by several threads at once, so a file larger than memory scans
 * like any other. Missing values (. and .a to .z) read as NULL; strL values
 * are found through an index of the file's strL section, built the first
 * time a query reads a strL column. Variable names, labels, formats and
 * value labels are recorded so that results keep them.
 *
 * Strings are passed on as stored: UTF-8 in formats 118 and 119, the file's
 * own encoding (often Latin-1) in older ones, as Stata's use does.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#include "csql/csql_internal.h"
#include "csql/csql_scan.h"

/* ------------------------------------------------------------ mapping */

typedef struct
{
    const unsigned char *base;
    uint64_t size;
#ifdef _WIN32
    HANDLE file, mapping;
#else
    int fd;
#endif
} dta_map;

static int dta_map_open(dta_map *m, const char *path, char *msg, size_t cap)
{
    memset(m, 0, sizeof(*m));
#ifdef _WIN32
    int n = MultiByteToWideChar(CP_UTF8, 0, path, -1, NULL, 0);
    wchar_t *w = n > 0 ? malloc((size_t)n * sizeof(wchar_t)) : NULL;
    if (!w || MultiByteToWideChar(CP_UTF8, 0, path, -1, w, n) <= 0) {
        free(w);
        snprintf(msg, cap, "file %s not found", path);
        return 601;
    }
    m->file = CreateFileW(w, GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, OPEN_EXISTING,
                          FILE_ATTRIBUTE_NORMAL, NULL);
    free(w);
    if (m->file == INVALID_HANDLE_VALUE) {
        m->file = NULL;
        snprintf(msg, cap, "file %s not found", path);
        return 601;
    }
    LARGE_INTEGER sz;
    if (!GetFileSizeEx(m->file, &sz)) {
        CloseHandle(m->file);
        snprintf(msg, cap, "could not read %s", path);
        return 692;
    }
    m->size = (uint64_t)sz.QuadPart;
    if (m->size == 0) {
        CloseHandle(m->file);
        snprintf(msg, cap, "%s is empty", path);
        return 610;
    }
    m->mapping = CreateFileMappingW(m->file, NULL, PAGE_READONLY, 0, 0, NULL);
    if (m->mapping) m->base = MapViewOfFile(m->mapping, FILE_MAP_READ, 0, 0, 0);
    if (!m->base) {
        if (m->mapping) CloseHandle(m->mapping);
        CloseHandle(m->file);
        snprintf(msg, cap, "could not map %s into memory", path);
        return 692;
    }
#else
    struct stat st;
    m->fd = open(path, O_RDONLY);
    if (m->fd < 0) {
        snprintf(msg, cap, "file %s not found", path);
        return 601;
    }
    if (fstat(m->fd, &st) != 0 || st.st_size <= 0) {
        close(m->fd);
        snprintf(msg, cap, "could not read %s", path);
        return 610;
    }
    m->size = (uint64_t)st.st_size;
    void *p = mmap(NULL, (size_t)m->size, PROT_READ, MAP_SHARED, m->fd, 0);
    if (p == MAP_FAILED) {
        close(m->fd);
        snprintf(msg, cap, "could not map %s into memory", path);
        return 692;
    }
    m->base = p;
#endif
    return 0;
}

static void dta_map_close(dta_map *m)
{
    if (!m->base) return;
#ifdef _WIN32
    UnmapViewOfFile((void *)m->base);
    CloseHandle(m->mapping);
    CloseHandle(m->file);
#else
    munmap((void *)m->base, (size_t)m->size);
    close(m->fd);
#endif
    m->base = NULL;
}

/* ------------------------------------------------------------ file */

typedef struct
{
    uint64_t key; /* v << 40 | o; 0 = empty slot */
    uint64_t off; /* contents */
    uint32_t len;
} gso_slot;

typedef struct dta_file
{
    char path[4096];
    dta_map map;
    int release, msf;
    uint64_t nobs;
    uint32_t nvar;
    uint16_t *types;    /* 117+ codes */
    uint32_t *offs;     /* field offsets in a record */
    uint64_t width;     /* record bytes */
    uint64_t data_off;  /* first record */
    uint64_t strls_off, strls_end;
    csql_source *meta;  /* owned by the ctx source list */
    /* strL index, built once */
    pthread_mutex_t gso_lock;
    int gso_state; /* 0 not built, 1 built, -1 failed */
    gso_slot *gso;
    size_t gso_cap;
    struct dta_file *next;
} dta_file;

static uint64_t rd_u(const dta_file *f, const unsigned char *p, int n)
{
    uint64_t v = 0;
    if (f->msf)
        for (int i = 0; i < n; i++) v = (v << 8) | p[i];
    else
        for (int i = n - 1; i >= 0; i--) v = (v << 8) | p[i];
    return v;
}

static const unsigned char *find_tag(const unsigned char *p, const unsigned char *end, const char *tag)
{
    size_t n = strlen(tag);
    while (p + n <= end) {
        const unsigned char *q = memchr(p, tag[0], (size_t)(end - p));
        if (!q || q + n > end) return NULL;
        if (memcmp(q, tag, n) == 0) return q;
        p = q + 1;
    }
    return NULL;
}

static void copy_field(char *dst, size_t cap, const unsigned char *src, size_t width)
{
    size_t n = 0;
    while (n < width && src[n]) n++;
    if (n >= cap) n = cap - 1;
    memcpy(dst, src, n);
    dst[n] = '\0';
}

static int dta_type_width(uint16_t t)
{
    if (t >= 1 && t <= CSQL_STR_MAX) return t;
    switch (t) {
    case CTOOLS_DTA_CODE_STRL: return 8;
    case CTOOLS_DTA_CODE_DOUBLE: return 8;
    case CTOOLS_DTA_CODE_FLOAT: return 4;
    case CTOOLS_DTA_CODE_LONG: return 4;
    case CTOOLS_DTA_CODE_INT: return 2;
    case CTOOLS_DTA_CODE_BYTE: return 1;
    default: return 0;
    }
}

/* Value-label tables: n, txtlen, off[n], val[n], txt. */
static int add_vlabel(dta_file *f, csql_source *s, const char *name, const unsigned char *t,
                      uint64_t len)
{
    if (len < 8) return 610;
    uint64_t n = rd_u(f, t, 4), txtlen = rd_u(f, t + 4, 4);
    if (8 + 8 * n + txtlen > len) return 610;
    csql_vlabel *x = realloc(s->labels, (s->nlabels + 1) * sizeof(*x));
    if (!x) return 920;
    s->labels = x;
    csql_vlabel *lab = &s->labels[s->nlabels];
    memset(lab, 0, sizeof(*lab));
    snprintf(lab->name, sizeof(lab->name), "%s", name);
    lab->values = malloc((n ? n : 1) * sizeof(int32_t));
    lab->texts = calloc(n ? n : 1, sizeof(char *));
    if (!lab->values || !lab->texts) {
        free(lab->values);
        free(lab->texts);
        return 920;
    }
    s->nlabels++;
    const unsigned char *offs = t + 8, *vals = t + 8 + 4 * n, *txt = t + 8 + 8 * n;
    for (uint64_t k = 0; k < n; k++) {
        uint64_t o = rd_u(f, offs + 4 * k, 4);
        if (o >= txtlen) return 610;
        size_t l = 0;
        while (o + l < txtlen && txt[o + l]) l++;
        char *text = malloc(l + 1);
        if (!text) return 920;
        memcpy(text, txt + o, l);
        text[l] = '\0';
        lab->values[k] = (int32_t)(uint32_t)rd_u(f, vals + 4 * k, 4);
        lab->texts[k] = text;
        lab->n++;
    }
    return 0;
}

static int dta_parse_new(dta_file *f, csql_source *s, char *msg, size_t cap)
{
    const unsigned char *b = f->map.base, *end = b + f->map.size;
    const unsigned char *p = find_tag(b, end, "<release>");
    if (!p || p + 12 > end) return 610;
    f->release = (p[9] - '0') * 100 + (p[10] - '0') * 10 + (p[11] - '0');
    if (f->release < 117 || f->release > 119) {
        snprintf(msg, cap, "%s: .dta format %d is not supported", f->path, f->release);
        return 610;
    }
    p = find_tag(p, end, "<byteorder>");
    if (!p || p + 14 > end) return 610;
    f->msf = memcmp(p + 11, "MSF", 3) == 0;
    p = find_tag(p, end, "<K>");
    if (!p) return 610;
    int kbytes = f->release == 119 ? 4 : 2;
    f->nvar = (uint32_t)rd_u(f, p + 3, kbytes);
    p = find_tag(p, end, "<N>");
    if (!p) return 610;
    int nbytes = f->release == 117 ? 4 : 8;
    f->nobs = rd_u(f, p + 3, nbytes);
    /* </N><label>, the data label with its length, </label><timestamp>,
     * the timestamp with its length, </timestamp></header><map> */
    p += 3 + nbytes;
    if (p + 11 > end || memcmp(p, "</N><label>", 11)) return 610;
    p += 11;
    int lbytes = f->release == 117 ? 1 : 2;
    if (p + lbytes > end) return 610;
    p += lbytes + rd_u(f, p, lbytes);
    if (p + 19 > end || memcmp(p, "</label><timestamp>", 19)) return 610;
    p += 19;
    if (p + 1 > end) return 610;
    p += 1 + p[0];
    if (p + 26 > end || memcmp(p, "</timestamp></header><map>", 26)) return 610;
    p += 21;
    if (p + 5 + 14 * 8 > end) return 610;
    uint64_t map[14];
    for (int i = 0; i < 14; i++) map[i] = rd_u(f, p + 5 + 8 * i, 8);
    for (int i = 2; i < 14; i++)
        if (map[i] > f->map.size) return 610;

    size_t namew = f->release == 117 ? 33 : 129, fmtw = f->release == 117 ? 49 : 57;
    size_t labw = f->release == 117 ? 81 : 321;
    uint32_t K = f->nvar;
    const unsigned char *types = b + map[2] + strlen("<variable_types>");
    const unsigned char *names = b + map[3] + strlen("<varnames>");
    const unsigned char *fmts = b + map[5] + strlen("<formats>");
    const unsigned char *vlns = b + map[6] + strlen("<value_label_names>");
    const unsigned char *labs = b + map[7] + strlen("<variable_labels>");
    if (types + 2 * (uint64_t)K > end || names + namew * K > end || fmts + fmtw * K > end ||
        vlns + namew * K > end || labs + labw * K > end)
        return 610;
    f->types = calloc(K ? K : 1, sizeof(uint16_t));
    f->offs = calloc(K ? K : 1, sizeof(uint32_t));
    s->vars = calloc(K ? K : 1, sizeof(csql_meta_var));
    if (!f->types || !f->offs || !s->vars) return 920;
    s->nvars = K;
    uint64_t w = 0;
    for (uint32_t j = 0; j < K; j++) {
        uint16_t t = (uint16_t)rd_u(f, types + 2 * j, 2);
        int tw = dta_type_width(t);
        if (!tw) return 610;
        f->types[j] = t;
        f->offs[j] = (uint32_t)w;
        w += (uint64_t)tw;
        csql_meta_var *v = &s->vars[j];
        v->type = t;
        copy_field(v->name, sizeof(v->name), names + namew * j, namew);
        copy_field(v->format, sizeof(v->format), fmts + fmtw * j, fmtw);
        copy_field(v->vlabel, sizeof(v->vlabel), vlns + namew * j, namew);
        copy_field(v->label, sizeof(v->label), labs + labw * j, labw);
    }
    f->width = w;
    f->data_off = map[9] + strlen("<data>");
    if (f->data_off + f->width * f->nobs > f->map.size) return 610;
    f->strls_off = map[10] + strlen("<strls>");
    f->strls_end = map[11];

    /* value labels */
    const unsigned char *q = b + map[11] + strlen("<value_labels>"), *vend = b + map[12];
    if (vend > end) vend = end;
    while (q + 5 <= vend && memcmp(q, "<lbl>", 5) == 0) {
        q += 5;
        if (q + 4 + namew + 3 > vend) return 610;
        uint64_t len = rd_u(f, q, 4);
        char name[CSQL_NAME_BYTES];
        copy_field(name, sizeof(name), q + 4, namew);
        const unsigned char *t = q + 4 + namew + 3;
        if (t + len > vend) return 610;
        int rc = add_vlabel(f, s, name, t, len);
        if (rc) return rc;
        q = t + len;
        if (q + 6 <= vend && memcmp(q, "</lbl>", 6) == 0) q += 6;
    }
    return 0;
}

static int dta_parse_old(dta_file *f, csql_source *s, char *msg, size_t cap)
{
    const unsigned char *b = f->map.base, *end = b + f->map.size;
    f->release = b[0];
    if (f->release < 113 || f->release > 115) {
        snprintf(msg, cap, "%s: .dta format %d is not supported (use and save it in Stata first)",
                 f->path, f->release);
        return 610;
    }
    if (f->map.size < 109) return 610;
    f->msf = b[1] == 1;
    f->nvar = (uint32_t)rd_u(f, b + 4, 2);
    f->nobs = rd_u(f, b + 6, 4);
    uint32_t K = f->nvar;
    size_t fmtw = f->release == 113 ? 12 : 49;
    const unsigned char *p = b + 109;
    const unsigned char *typ = p, *names = typ + K, *srt = names + 33 * (uint64_t)K;
    const unsigned char *fmts = srt + 2 * ((uint64_t)K + 1), *vlns = fmts + fmtw * K;
    const unsigned char *labs = vlns + 33 * (uint64_t)K, *ext = labs + 81 * (uint64_t)K;
    if (ext > end) return 610;
    f->types = calloc(K ? K : 1, sizeof(uint16_t));
    f->offs = calloc(K ? K : 1, sizeof(uint32_t));
    s->vars = calloc(K ? K : 1, sizeof(csql_meta_var));
    if (!f->types || !f->offs || !s->vars) return 920;
    s->nvars = K;
    uint64_t w = 0;
    for (uint32_t j = 0; j < K; j++) {
        unsigned t = typ[j];
        uint16_t code;
        if (t >= 1 && t <= 244) code = (uint16_t)t;
        else if (t == 251) code = CTOOLS_DTA_CODE_BYTE;
        else if (t == 252) code = CTOOLS_DTA_CODE_INT;
        else if (t == 253) code = CTOOLS_DTA_CODE_LONG;
        else if (t == 254) code = CTOOLS_DTA_CODE_FLOAT;
        else if (t == 255) code = CTOOLS_DTA_CODE_DOUBLE;
        else return 610;
        f->types[j] = code;
        f->offs[j] = (uint32_t)w;
        w += (uint64_t)dta_type_width(code);
        csql_meta_var *v = &s->vars[j];
        v->type = code;
        copy_field(v->name, sizeof(v->name), names + 33 * j, 33);
        copy_field(v->format, sizeof(v->format), fmts + fmtw * j, fmtw);
        copy_field(v->vlabel, sizeof(v->vlabel), vlns + 33 * j, 33);
        copy_field(v->label, sizeof(v->label), labs + 81 * j, 81);
    }
    f->width = w;
    /* expansion fields: type (1), length (4), contents; ends with 0, 0 */
    p = ext;
    for (;;) {
        if (p + 5 > end) return 610;
        unsigned type = p[0];
        uint64_t len = rd_u(f, p + 1, 4);
        p += 5;
        if (type == 0 && len == 0) break;
        if (p + len > end) return 610;
        p += len;
    }
    f->data_off = (uint64_t)(p - b);
    if (f->data_off + f->width * f->nobs > f->map.size) return 610;
    f->strls_off = f->strls_end = 0;
    const unsigned char *q = b + f->data_off + f->width * f->nobs;
    while (q + 4 + 33 + 3 <= end) {
        uint64_t len = rd_u(f, q, 4);
        char name[CSQL_NAME_BYTES];
        copy_field(name, sizeof(name), q + 4, 33);
        const unsigned char *t = q + 40;
        if (t + len > end) break;
        int rc = add_vlabel(f, s, name, t, len);
        if (rc) return rc;
        q = t + len;
    }
    return 0;
}

static void dta_file_free(dta_file *f)
{
    if (!f) return;
    dta_map_close(&f->map);
    free(f->types);
    free(f->offs);
    free(f->gso);
    pthread_mutex_destroy(&f->gso_lock);
    free(f);
}

/* Files opened by this call; freed by csql_dta_close_all. */
static dta_file *g_files;

void csql_dta_close_all(void)
{
    while (g_files) {
        dta_file *n = g_files->next;
        dta_file_free(g_files);
        g_files = n;
    }
}

static int dta_open(csql_ctx *ctx, const char *path, dta_file **out, char *msg, size_t cap)
{
    for (dta_file *f = g_files; f; f = f->next)
        if (strcmp(f->path, path) == 0) {
            *out = f;
            return 0;
        }
    dta_file *f = calloc(1, sizeof(*f));
    csql_source *s = calloc(1, sizeof(*s));
    if (!f || !s) {
        free(f);
        free(s);
        snprintf(msg, cap, "out of memory");
        return 920;
    }
    pthread_mutex_init(&f->gso_lock, NULL);
    snprintf(f->path, sizeof(f->path), "%s", path);
    snprintf(s->origin, sizeof(s->origin), "%s", path);
    int rc = dta_map_open(&f->map, path, msg, cap);
    if (!rc) {
        msg[0] = '\0';
        rc = f->map.base[0] == '<' ? dta_parse_new(f, s, msg, cap) : dta_parse_old(f, s, msg, cap);
        if (!rc) rc = csql_sql_names(s->vars, s->nvars);
        if (rc && !msg[0])
            snprintf(msg, cap, rc == 920 ? "out of memory" : "%s is not a Stata dataset or is damaged",
                     path);
    }
    if (rc) {
        csql_source_free(s);
        dta_file_free(f);
        return rc;
    }
    f->meta = s;
    s->next = ctx->sources;
    ctx->sources = s;
    f->next = g_files;
    g_files = f;
    *out = f;
    return 0;
}

/* ------------------------------------------------------------ strL index */

static uint64_t gso_hash(uint64_t k)
{
    k ^= k >> 33;
    k *= UINT64_C(0xff51afd7ed558ccd);
    k ^= k >> 33;
    return k;
}

static int gso_build(dta_file *f)
{
    pthread_mutex_lock(&f->gso_lock);
    if (f->gso_state) {
        int st = f->gso_state;
        pthread_mutex_unlock(&f->gso_lock);
        return st > 0 ? 0 : 610;
    }
    const unsigned char *b = f->map.base, *p = b + f->strls_off, *end = b + f->strls_end;
    int obytes = f->release == 117 ? 4 : 8;
    size_t count = 0;
    /* count, then insert */
    for (const unsigned char *q = p; q + 3 + 4 + obytes + 5 <= end && !memcmp(q, "GSO", 3);) {
        uint64_t len = rd_u(f, q + 3 + 4 + obytes + 1, 4);
        q += 3 + 4 + obytes + 5 + len;
        count++;
    }
    size_t cap = 16;
    while (cap < count * 2) cap *= 2;
    gso_slot *t = calloc(cap, sizeof(*t));
    if (!t) {
        f->gso_state = -1;
        pthread_mutex_unlock(&f->gso_lock);
        return 920;
    }
    for (const unsigned char *q = p; q + 3 + 4 + obytes + 5 <= end && !memcmp(q, "GSO", 3);) {
        uint64_t v = rd_u(f, q + 3, 4), o = rd_u(f, q + 7, obytes);
        unsigned type = q[3 + 4 + obytes];
        uint64_t len = rd_u(f, q + 3 + 4 + obytes + 1, 4);
        const unsigned char *data = q + 3 + 4 + obytes + 5;
        if (data + len > end) break;
        uint64_t key = (v << 40) | o;
        size_t h = (size_t)gso_hash(key) & (cap - 1);
        while (t[h].key && t[h].key != key) h = (h + 1) & (cap - 1);
        t[h].key = key;
        t[h].off = (uint64_t)(data - b);
        /* text strLs (130) store their terminating NUL */
        t[h].len = (uint32_t)(type == 130 && len > 0 && data[len - 1] == 0 ? len - 1 : len);
        q = data + len;
    }
    f->gso = t;
    f->gso_cap = cap;
    f->gso_state = 1;
    pthread_mutex_unlock(&f->gso_lock);
    return 0;
}

static int gso_find(const dta_file *f, uint64_t v, uint64_t o, const char **ptr, uint32_t *len)
{
    uint64_t key = (v << 40) | o;
    size_t h = (size_t)gso_hash(key) & (f->gso_cap - 1);
    while (f->gso[h].key) {
        if (f->gso[h].key == key) {
            *ptr = (const char *)f->map.base + f->gso[h].off;
            *len = f->gso[h].len;
            return 1;
        }
        h = (h + 1) & (f->gso_cap - 1);
    }
    return 0;
}

/* ------------------------------------------------------------ decoding */

/* A field as a Stata double (missing values as such). */
static double dta_number(const dta_file *f, uint16_t t, const unsigned char *p)
{
    switch (t) {
    case CTOOLS_DTA_CODE_BYTE: {
        int v = (int8_t)p[0];
        return v > 100 ? csql_missing() : (double)v;
    }
    case CTOOLS_DTA_CODE_INT: {
        int v = (int16_t)(uint16_t)rd_u(f, p, 2);
        return v > 32740 ? csql_missing() : (double)v;
    }
    case CTOOLS_DTA_CODE_LONG: {
        int32_t v = (int32_t)(uint32_t)rd_u(f, p, 4);
        return v > 2147483620 ? csql_missing() : (double)v;
    }
    case CTOOLS_DTA_CODE_FLOAT: {
        uint32_t bits = (uint32_t)rd_u(f, p, 4);
        float x;
        memcpy(&x, &bits, 4);
        /* float missing values start at 2^127 */
        return (bits & 0x7fffffffu) >= 0x7f000000u && !(bits & 0x80000000u) ? csql_missing()
                                                                              : (double)x;
    }
    default: {
        uint64_t bits = rd_u(f, p, 8);
        double x;
        memcpy(&x, &bits, 8);
        return x;
    }
    }
}

typedef struct
{
    csql_scan_def def; /* first: the engine reads it */
    csql_ctx *ctx;
    dta_file *file;
    csql_map *map;
    int *decode; /* per variable: decode its value labels */
    const csql_vlabel **labels;
} dta_bind;

static void dta_bind_free(void *p)
{
    dta_bind *b = p;
    if (!b) return;
    free(b->map);
    free(b->decode);
    free(b->labels);
    free(b);
}

static const char *label_text(const csql_vlabel *lab, double z)
{
    if (z != (double)(int32_t)z) return NULL;
    int32_t v = (int32_t)z;
    for (size_t k = 0; k < lab->n; k++)
        if (lab->values[k] == v) return lab->texts[k];
    return NULL;
}

/* Rows [pre->r0, pre->r0 + pre->n) from the mapped file (any thread).
 * Strings point into the mapping, which outlives the query. */
static int dta_decode(void *arg, const idx_t *proj, size_t nproj, csql_pre *pre)
{
    dta_bind *b = arg;
    const dta_file *f = b->file;
    size_t n = pre->n, need = 0;
    const unsigned char *rows = f->map.base + f->data_off + (uint64_t)pre->r0 * f->width;
    for (size_t c = 0; c < nproj; c++)
        if (proj[c] < f->nvar && b->decode[proj[c]]) need += n * 32;
    if (need && csql_arena_reserve(&pre->arena, need) != 0) return 920;
    for (size_t c = 0; c < nproj; c++) {
        idx_t j = proj[c];
        if (j >= f->nvar) continue;
        uint16_t t = f->types[j];
        const unsigned char *cell = rows + f->offs[j];
        uint64_t *valid = pre->valid[c];
        if (t >= 1 && t <= CSQL_STR_MAX) {
            duckdb_string_t *out = pre->data[c];
            for (size_t i = 0; i < n; i++, cell += f->width) {
                const unsigned char *z = memchr(cell, 0, t);
                csql_str_set(&out[i], (const char *)cell, (uint32_t)(z ? (size_t)(z - cell) : t));
                csql_valid_set(valid, i);
            }
        } else if (t == CTOOLS_DTA_CODE_STRL) {
            duckdb_string_t *out = pre->data[c];
            for (size_t i = 0; i < n; i++, cell += f->width) {
                uint64_t v, o;
                if (f->release == 117) {
                    v = rd_u(f, cell, 4);
                    o = rd_u(f, cell + 4, 4);
                } else if (f->release == 118) {
                    v = rd_u(f, cell, 2);
                    o = rd_u(f, cell + 2, 6);
                } else {
                    v = rd_u(f, cell, 3);
                    o = rd_u(f, cell + 3, 5);
                }
                const char *ptr = "";
                uint32_t len = 0;
                if ((v || o) && !gso_find(f, v, o, &ptr, &len)) {
                    ptr = "";
                    len = 0;
                }
                csql_str_set(&out[i], ptr, len);
                csql_valid_set(valid, i);
            }
        } else if (b->decode[j]) {
            duckdb_string_t *out = pre->data[c];
            char num[48];
            for (size_t i = 0; i < n; i++, cell += f->width) {
                double z = dta_number(f, t, cell);
                if (csql_is_missing(z)) continue;
                const char *text = label_text(b->labels[j], z);
                size_t len;
                if (text) {
                    len = strlen(text); /* the label's text lives as long as the query */
                } else {
                    len = (size_t)snprintf(num, sizeof(num), "%.15g", z);
                    text = num;
                    if (len > 12) {
                        text = csql_arena_add(&pre->arena, num, len);
                        if (!text) return 920;
                        pre->copy[c] = 1;
                    }
                }
                csql_str_set(&out[i], text, (uint32_t)len);
                csql_valid_set(valid, i);
            }
        } else {
            void *data = pre->data[c];
            csql_map m = b->map[j];
            for (size_t i = 0; i < n; i++, cell += f->width) {
                double z = dta_number(f, t, cell);
                int bad = 1;
                if (!csql_is_missing(z)) csql_map_store(m, data, i, z, &bad);
                if (!bad) csql_valid_set(valid, i);
            }
        }
    }
    return 0;
}

static void dta_scan_bind(duckdb_bind_info info)
{
    csql_ctx *ctx = DK(bind_get_extra_info)(info);
    duckdb_value pv = DK(bind_get_parameter)(info, 0);
    char *path = pv ? DK(get_varchar)(pv) : NULL;
    char msg[4400] = "";
    dta_file *f = NULL;
    int decode = 0;
    if (pv) DK(destroy_value)(&pv);
    duckdb_value dv = DK(bind_get_named_parameter)(info, "decode");
    if (dv) {
        decode = DK(get_bool)(dv);
        DK(destroy_value)(&dv);
    }
    if (!path) {
        DK(bind_set_error)(info, "read_dta: no file name");
        return;
    }
    int rc = dta_open(ctx, path, &f, msg, sizeof(msg));
    DK(free)(path);
    if (rc) {
        DK(bind_set_error)(info, msg);
        return;
    }
    if (f->nvar == 0) {
        DK(bind_set_error)(info, "read_dta: the dataset has no variables");
        return;
    }
    dta_bind *b = calloc(1, sizeof(*b));
    if (b) {
        b->map = calloc(f->nvar, sizeof(*b->map));
        b->decode = calloc(f->nvar, sizeof(int));
        b->labels = calloc(f->nvar, sizeof(*b->labels));
    }
    if (!b || !b->map || !b->decode || !b->labels) {
        dta_bind_free(b);
        DK(bind_set_error)(info, "out of memory");
        return;
    }
    b->ctx = ctx;
    b->file = f;
    for (uint32_t j = 0; j < f->nvar; j++) {
        const csql_meta_var *v = &f->meta->vars[j];
        b->map[j] = csql_map_stata_var(v->type, v->format, ctx->datetime);
        if (decode && v->vlabel[0] && b->map[j] != CSQL_MAP_VARCHAR) {
            b->labels[j] = csql_source_label(f->meta, v->vlabel);
            if (b->labels[j]) {
                b->decode[j] = 1;
                b->map[j] = CSQL_MAP_VARCHAR;
            }
        }
        duckdb_logical_type t = csql_map_logical(b->map[j]);
        DK(bind_add_result_column)(info, v->sqlname, t);
        DK(destroy_logical_type)(&t);
    }
    b->def = (csql_scan_def){ctx, b, (int64_t)f->nobs, f->nvar, b->map, NULL, 1, dta_decode};
    DK(bind_set_cardinality)(info, (idx_t)f->nobs, true);
    DK(bind_set_bind_data)(info, b, dta_bind_free);
    if (!f->meta->used) {
        f->meta->used = 1;
        f->meta->order = ++ctx->nused;
    }
}

/* strL columns need the file's strL index: build it before any thread reads. */
static void dta_scan_init(duckdb_init_info info)
{
    dta_bind *b = DK(init_get_bind_data)(info);
    idx_t n = DK(init_get_column_count)(info);
    for (idx_t i = 0; i < n; i++) {
        idx_t j = DK(init_get_column_index)(info, i);
        if (j < b->file->nvar && b->file->types[j] == CTOOLS_DTA_CODE_STRL) {
            if (gso_build(b->file) != 0) {
                DK(init_set_error)(info, "read_dta: could not read the strL values");
                return;
            }
            break;
        }
    }
    csql_scan_init(info);
}

/* FROM 'file.dta' */
static void dta_replacement(duckdb_replacement_scan_info info, const char *table_name, void *data)
{
    (void)data;
    size_t n = strlen(table_name);
    if (n < 5) return;
    const char *ext = table_name + n - 4;
    if (!(ext[0] == '.' && (ext[1] == 'd' || ext[1] == 'D') && (ext[2] == 't' || ext[2] == 'T') &&
          (ext[3] == 'a' || ext[3] == 'A')))
        return;
    duckdb_value v = DK(create_varchar)(table_name);
    DK(replacement_scan_set_function_name)(info, "read_dta");
    DK(replacement_scan_add_parameter)(info, v);
    DK(destroy_value)(&v);
}

int csql_dta_register(csql_ctx *ctx)
{
    duckdb_table_function f = DK(create_table_function)();
    duckdb_logical_type vt, bt;
    duckdb_state st;
    if (!f) return 920;
    vt = DK(create_logical_type)(DUCKDB_TYPE_VARCHAR);
    bt = DK(create_logical_type)(DUCKDB_TYPE_BOOLEAN);
    DK(table_function_set_name)(f, "read_dta");
    DK(table_function_add_parameter)(f, vt);
    DK(table_function_add_named_parameter)(f, "decode", bt);
    DK(table_function_set_extra_info)(f, ctx, NULL);
    DK(table_function_set_bind)(f, dta_scan_bind);
    DK(table_function_set_init)(f, dta_scan_init);
    DK(table_function_set_local_init)(f, csql_scan_local_init);
    DK(table_function_set_function)(f, csql_scan_function);
    DK(table_function_supports_projection_pushdown)(f, true);
    st = DK(register_table_function)(ctx->con, f);
    DK(destroy_table_function)(&f);
    DK(destroy_logical_type)(&vt);
    DK(destroy_logical_type)(&bt);
    if (st != DuckDBSuccess) return 198;
    DK(add_replacement_scan)(ctx->db, dta_replacement, NULL, NULL);
    return 0;
}

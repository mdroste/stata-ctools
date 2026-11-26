/*
 * ctools_dta_schema.c - read and write the description of a Stata dataset
 * (format 118 or 119, little-endian) without its observations.
 *
 * See ctools_dta_schema.h for the contract and docs/development/DTA_WRITER.md
 * for the file layout. The reader bounds every field against the bytes it
 * holds; the writer emits the sections in file order and fills in the map.
 */
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ctools_dta_schema.h"
#include "ctools_file.h"

#define SCHEMA_MAP_ENTRIES 14
#define SCHEMA_MAX_VARS 16777215 /* a strL cell numbers variables in 3 bytes */

#if defined(__GNUC__) || defined(__clang__)
__attribute__((format(printf, 4, 5)))
#endif
static int schema_fail(char *msg, size_t cap, int rc, const char *fmt, ...)
{
    if (msg && cap)
    {
        va_list ap;
        va_start(ap, fmt);
        vsnprintf(msg, cap, fmt, ap);
        va_end(ap);
    }
    return rc;
}

void ctools_dta_schema_init(ctools_dta_schema *s)
{
    if (s)
        memset(s, 0, sizeof(*s));
}

void ctools_dta_schema_free(ctools_dta_schema *s)
{
    if (!s)
        return;
    free(s->vars);
    free(s->sort);
    for (size_t i = 0; i < s->nchars; i++)
        free(s->chars[i].text);
    free(s->chars);
    for (size_t i = 0; i < s->nlabels; i++)
        free(s->labels[i].table);
    free(s->labels);
    ctools_dta_schema_init(s);
}

/* ------------------------------------------------------------------------
 * Reading
 * ------------------------------------------------------------------------ */

typedef struct
{
    const unsigned char *p;
    size_t len, at;
    int bad;
} schema_cursor;

static uint64_t schema_le(const unsigned char *p, int n)
{
    uint64_t v = 0;
    for (int i = n - 1; i >= 0; i--)
        v = (v << 8) | p[i];
    return v;
}

static const unsigned char *schema_take(schema_cursor *c, size_t n)
{
    if (c->bad || n > c->len - c->at)
    {
        c->bad = 1;
        return NULL;
    }
    const unsigned char *p = c->p + c->at;
    c->at += n;
    return p;
}

static uint64_t schema_uint(schema_cursor *c, int n)
{
    const unsigned char *p = schema_take(c, (size_t)n);
    return p ? schema_le(p, n) : 0;
}

static void schema_tag(schema_cursor *c, const char *tag)
{
    size_t n = strlen(tag);
    const unsigned char *p = schema_take(c, n);
    if (p && memcmp(p, tag, n))
        c->bad = 1;
}

static int schema_peek(const schema_cursor *c, const char *tag)
{
    size_t n = strlen(tag);
    return !c->bad && n <= c->len - c->at && !memcmp(c->p + c->at, tag, n);
}

/* A fixed-width field into dst[n]: the bytes up to the first NUL, and a NUL
 * always (a damaged field without one is cut at n - 1). */
static void schema_field(schema_cursor *c, char *dst, size_t n)
{
    const unsigned char *p = schema_take(c, n);
    if (!p)
        return;
    memcpy(dst, p, n);
    dst[n - 1] = '\0';
}

static int schema_read_at(FILE *fp, uint64_t at, unsigned char *buf, size_t n)
{
    if (ctools_fseek(fp, (long long)at, SEEK_SET) != 0)
        return -1;
    return fread(buf, 1, n, fp) == n ? 0 : -1;
}

static int schema_push_char(ctools_dta_schema *s, const char *owner, const char *name,
                            const char *text, size_t text_len)
{
    if (s->nchars == s->chars_cap)
    {
        size_t cap = s->chars_cap ? s->chars_cap * 2 : 32;
        ctools_dta_schema_char *grown = realloc(s->chars, cap * sizeof(*grown));
        if (!grown)
            return 920;
        s->chars = grown;
        s->chars_cap = cap;
    }
    ctools_dta_schema_char *ch = &s->chars[s->nchars];
    memset(ch, 0, sizeof(*ch));
    ch->text = malloc(text_len + 1);
    if (!ch->text)
        return 920;
    memcpy(ch->text, text, text_len);
    ch->text[text_len] = '\0';
    snprintf(ch->owner, sizeof(ch->owner), "%s", owner);
    snprintf(ch->name, sizeof(ch->name), "%s", name);
    s->nchars++;
    return 0;
}

static int schema_push_label(ctools_dta_schema *s, const char *name, const unsigned char *table,
                             size_t table_bytes)
{
    if (s->nlabels == s->labels_cap)
    {
        size_t cap = s->labels_cap ? s->labels_cap * 2 : 16;
        ctools_dta_schema_label *grown = realloc(s->labels, cap * sizeof(*grown));
        if (!grown)
            return 920;
        s->labels = grown;
        s->labels_cap = cap;
    }
    ctools_dta_schema_label *l = &s->labels[s->nlabels];
    memset(l, 0, sizeof(*l));
    l->table = malloc(table_bytes ? table_bytes : 1);
    if (!l->table)
        return 920;
    memcpy(l->table, table, table_bytes);
    l->table_bytes = table_bytes;
    snprintf(l->name, sizeof(l->name), "%s", name);
    s->nlabels++;
    return 0;
}

int ctools_dta_schema_read(const char *path, ctools_dta_schema *s, char *msg, size_t cap)
{
    if (!s || !path || !*path)
        return schema_fail(msg, cap, 198, "no file name");
    ctools_dta_schema_init(s);
    FILE *fp = ctools_fopen(path, "rb");
    if (!fp)
        return schema_fail(msg, cap, 601, "cannot open the dataset description");

    unsigned char *meta = NULL, *lbl = NULL;
    int rc = 0;

    /* The header up to the map says where the sections are. */
    unsigned char head[512];
    size_t got = fread(head, 1, sizeof(head), fp);
    schema_cursor c = {head, got, 0, 0};
    schema_tag(&c, "<stata_dta><header><release>");
    const unsigned char *rel = schema_take(&c, 3);
    int release = rel && !memcmp(rel, "118", 3) ? 118 : rel && !memcmp(rel, "119", 3) ? 119 : 0;
    schema_tag(&c, "</release><byteorder>");
    const unsigned char *order = schema_take(&c, 3);
    if (c.bad || !release || memcmp(order, "LSF", 3))
    {
        rc = schema_fail(msg, cap, 610, "the dataset description is not a format 118 or 119 file");
        goto done;
    }
    const int kbytes = release == 118 ? 2 : 4;
    schema_tag(&c, "</byteorder><K>");
    uint64_t K = schema_uint(&c, kbytes);
    schema_tag(&c, "</K><N>");
    (void)schema_uint(&c, 8);
    schema_tag(&c, "</N><label>");
    size_t label_len = (size_t)schema_uint(&c, 2);
    const unsigned char *label = schema_take(&c, label_len);
    schema_tag(&c, "</label><timestamp>");
    size_t stamp_len = (size_t)schema_uint(&c, 1);
    (void)schema_take(&c, stamp_len);
    schema_tag(&c, "</timestamp></header><map>");
    uint64_t map[SCHEMA_MAP_ENTRIES];
    for (int i = 0; i < SCHEMA_MAP_ENTRIES; i++)
        map[i] = schema_uint(&c, 8);
    if (c.bad || label_len >= CTOOLS_DTA_LABEL_BYTES || K > SCHEMA_MAX_VARS || map[9] < c.at ||
        map[9] > ((uint64_t)1 << 40) || map[11] < map[9] || map[12] < map[11] ||
        map[12] - map[11] > ((uint64_t)1 << 40))
    {
        rc = schema_fail(msg, cap, 610, "the dataset description has a damaged header");
        goto done;
    }
    s->release = release;
    memcpy(s->data_label, label, label_len);
    s->data_label[label_len] = '\0';

    /* Everything before the observations, then the value labels. */
    const size_t meta_len = (size_t)map[9], lbl_len = (size_t)(map[12] - map[11]);
    meta = malloc(meta_len ? meta_len : 1);
    lbl = malloc(lbl_len ? lbl_len : 1);
    if (!meta || !lbl)
    {
        rc = schema_fail(msg, cap, 920, "out of memory");
        goto done;
    }
    if (schema_read_at(fp, 0, meta, meta_len) || schema_read_at(fp, map[11], lbl, lbl_len))
    {
        rc = schema_fail(msg, cap, 610, "the dataset description is cut short");
        goto done;
    }

    const size_t nvars = (size_t)K;
    s->vars = calloc(nvars ? nvars : 1, sizeof(*s->vars));
    s->sort = calloc(nvars ? nvars : 1, sizeof(*s->sort));
    if (!s->vars || !s->sort)
    {
        rc = schema_fail(msg, cap, 920, "out of memory");
        goto done;
    }
    s->nvars = nvars;

    c = (schema_cursor){meta, meta_len, (size_t)map[2], 0};
    schema_tag(&c, "<variable_types>");
    for (size_t j = 0; j < nvars; j++)
        s->vars[j].type = (uint16_t)schema_uint(&c, 2);
    schema_tag(&c, "</variable_types><varnames>");
    for (size_t j = 0; j < nvars; j++)
        schema_field(&c, s->vars[j].name, CTOOLS_DTA_NAME_BYTES);
    schema_tag(&c, "</varnames><sortlist>");
    int sorting = 1;
    for (size_t j = 0; j <= nvars; j++)
    {
        uint64_t v = schema_uint(&c, kbytes);
        if (sorting && v >= 1 && v <= nvars && j < nvars)
            s->sort[s->nsort++] = (uint32_t)(v - 1);
        else
            sorting = 0;
    }
    schema_tag(&c, "</sortlist><formats>");
    for (size_t j = 0; j < nvars; j++)
        schema_field(&c, s->vars[j].format, CTOOLS_DTA_FORMAT_BYTES);
    schema_tag(&c, "</formats><value_label_names>");
    for (size_t j = 0; j < nvars; j++)
        schema_field(&c, s->vars[j].value_label, CTOOLS_DTA_NAME_BYTES);
    schema_tag(&c, "</value_label_names><variable_labels>");
    for (size_t j = 0; j < nvars; j++)
        schema_field(&c, s->vars[j].label, CTOOLS_DTA_LABEL_BYTES);
    schema_tag(&c, "</variable_labels><characteristics>");
    while (schema_peek(&c, "<ch>"))
    {
        schema_tag(&c, "<ch>");
        size_t n = (size_t)schema_uint(&c, 4);
        char owner[CTOOLS_DTA_NAME_BYTES] = "", name[CTOOLS_DTA_NAME_BYTES] = "";
        if (n < 2 * CTOOLS_DTA_NAME_BYTES)
            c.bad = 1;
        schema_field(&c, owner, CTOOLS_DTA_NAME_BYTES);
        schema_field(&c, name, CTOOLS_DTA_NAME_BYTES);
        const unsigned char *text = schema_take(&c, c.bad ? 0 : n - 2 * CTOOLS_DTA_NAME_BYTES);
        schema_tag(&c, "</ch>");
        if (c.bad)
            break;
        size_t text_len = n - 2 * CTOOLS_DTA_NAME_BYTES;
        text_len = strnlen((const char *)text, text_len);
        if ((rc = schema_push_char(s, owner, name, (const char *)text, text_len)))
        {
            schema_fail(msg, cap, rc, "out of memory");
            goto done;
        }
    }
    schema_tag(&c, "</characteristics>");
    if (c.bad)
    {
        rc = schema_fail(msg, cap, 610, "the dataset description is damaged");
        goto done;
    }

    c = (schema_cursor){lbl, lbl_len, 0, 0};
    schema_tag(&c, "<value_labels>");
    while (schema_peek(&c, "<lbl>"))
    {
        schema_tag(&c, "<lbl>");
        size_t n = (size_t)schema_uint(&c, 4);
        char name[CTOOLS_DTA_NAME_BYTES] = "";
        schema_field(&c, name, CTOOLS_DTA_NAME_BYTES);
        (void)schema_take(&c, 3);
        const unsigned char *table = schema_take(&c, n);
        schema_tag(&c, "</lbl>");
        if (c.bad)
            break;
        if ((rc = schema_push_label(s, name, table, n)))
        {
            schema_fail(msg, cap, rc, "out of memory");
            goto done;
        }
    }
    schema_tag(&c, "</value_labels>");
    if (c.bad)
        rc = schema_fail(msg, cap, 610, "the dataset description has damaged value labels");

done:
    fclose(fp);
    free(meta);
    free(lbl);
    if (rc)
        ctools_dta_schema_free(s);
    return rc;
}

/* ------------------------------------------------------------------------
 * Writing
 * ------------------------------------------------------------------------ */

typedef struct
{
    unsigned char *p;
    size_t len, cap;
    int failed;
} schema_bytes;

static unsigned char *schema_grow(schema_bytes *b, size_t n)
{
    if (b->failed)
        return NULL;
    if (n > b->cap - b->len)
    {
        size_t cap = b->cap ? b->cap : 4096;
        while (cap - b->len < n)
        {
            if (cap > SIZE_MAX / 2)
            {
                b->failed = 1;
                return NULL;
            }
            cap *= 2;
        }
        unsigned char *grown = realloc(b->p, cap);
        if (!grown)
        {
            b->failed = 1;
            return NULL;
        }
        b->p = grown;
        b->cap = cap;
    }
    unsigned char *at = b->p + b->len;
    b->len += n;
    return at;
}

static void schema_put(schema_bytes *b, const void *src, size_t n)
{
    unsigned char *at = schema_grow(b, n);
    if (at && n)
        memcpy(at, src, n);
}

static void schema_put_tag(schema_bytes *b, const char *tag)
{
    schema_put(b, tag, strlen(tag));
}

static void schema_put_le(schema_bytes *b, uint64_t v, int n)
{
    unsigned char *at = schema_grow(b, (size_t)n);
    for (int i = 0; at && i < n; i++)
        at[i] = (unsigned char)(v >> (8 * i));
}

/* text into a NUL-padded field of n bytes (cut to n - 1 bytes) */
static void schema_put_field(schema_bytes *b, const char *text, size_t n)
{
    unsigned char *at = schema_grow(b, n);
    if (!at)
        return;
    size_t len = strnlen(text, n - 1);
    memcpy(at, text, len);
    memset(at + len, 0, n - len);
}

int ctools_dta_schema_write(const char *path, const ctools_dta_schema *s, char *msg, size_t cap)
{
    if (!s || !path || !*path)
        return schema_fail(msg, cap, 198, "no file name");
    if (s->nvars > SCHEMA_MAX_VARS || s->nsort > s->nvars)
        return schema_fail(msg, cap, 198, "invalid dataset description");
    for (size_t i = 0; i < s->nsort; i++)
        if (s->sort[i] >= s->nvars)
            return schema_fail(msg, cap, 198, "invalid sort list");
    const int release = s->nvars > 32767 ? 119 : 118;
    const int kbytes = release == 118 ? 2 : 4;

    schema_bytes b = {0};
    uint64_t map[SCHEMA_MAP_ENTRIES] = {0};
    schema_put_tag(&b,
                   release == 118
                       ? "<stata_dta><header><release>118</release><byteorder>LSF</byteorder><K>"
                       : "<stata_dta><header><release>119</release><byteorder>LSF</byteorder><K>");
    schema_put_le(&b, s->nvars, kbytes);
    schema_put_tag(&b, "</K><N>");
    schema_put_le(&b, 0, 8);
    schema_put_tag(&b, "</N><label>");
    size_t label_len = strnlen(s->data_label, CTOOLS_DTA_LABEL_BYTES - 1);
    schema_put_le(&b, label_len, 2);
    schema_put(&b, s->data_label, label_len);
    schema_put_tag(&b, "</label><timestamp>");
    schema_put_le(&b, 0, 1);
    schema_put_tag(&b, "</timestamp></header>");
    map[1] = b.len;
    schema_put_tag(&b, "<map>");
    const size_t map_at = b.len;
    for (int i = 0; i < SCHEMA_MAP_ENTRIES; i++)
        schema_put_le(&b, 0, 8);
    schema_put_tag(&b, "</map>");
    map[2] = b.len;
    schema_put_tag(&b, "<variable_types>");
    for (size_t j = 0; j < s->nvars; j++)
        schema_put_le(&b, s->vars[j].type, 2);
    schema_put_tag(&b, "</variable_types>");
    map[3] = b.len;
    schema_put_tag(&b, "<varnames>");
    for (size_t j = 0; j < s->nvars; j++)
        schema_put_field(&b, s->vars[j].name, CTOOLS_DTA_NAME_BYTES);
    schema_put_tag(&b, "</varnames>");
    map[4] = b.len;
    schema_put_tag(&b, "<sortlist>");
    for (size_t j = 0; j <= s->nvars; j++)
        schema_put_le(&b, j < s->nsort ? (uint64_t)s->sort[j] + 1 : 0, kbytes);
    schema_put_tag(&b, "</sortlist>");
    map[5] = b.len;
    schema_put_tag(&b, "<formats>");
    for (size_t j = 0; j < s->nvars; j++)
        schema_put_field(&b, s->vars[j].format, CTOOLS_DTA_FORMAT_BYTES);
    schema_put_tag(&b, "</formats>");
    map[6] = b.len;
    schema_put_tag(&b, "<value_label_names>");
    for (size_t j = 0; j < s->nvars; j++)
        schema_put_field(&b, s->vars[j].value_label, CTOOLS_DTA_NAME_BYTES);
    schema_put_tag(&b, "</value_label_names>");
    map[7] = b.len;
    schema_put_tag(&b, "<variable_labels>");
    for (size_t j = 0; j < s->nvars; j++)
        schema_put_field(&b, s->vars[j].label, CTOOLS_DTA_LABEL_BYTES);
    schema_put_tag(&b, "</variable_labels>");
    map[8] = b.len;
    schema_put_tag(&b, "<characteristics>");
    for (size_t i = 0; i < s->nchars; i++)
    {
        const ctools_dta_schema_char *ch = &s->chars[i];
        size_t text = strlen(ch->text) + 1;
        if (text > UINT32_MAX - 2 * CTOOLS_DTA_NAME_BYTES)
        {
            free(b.p);
            return schema_fail(msg, cap, 198, "characteristic %s[%s] is too long", ch->owner,
                               ch->name);
        }
        schema_put_tag(&b, "<ch>");
        schema_put_le(&b, 2 * CTOOLS_DTA_NAME_BYTES + text, 4);
        schema_put_field(&b, ch->owner, CTOOLS_DTA_NAME_BYTES);
        schema_put_field(&b, ch->name, CTOOLS_DTA_NAME_BYTES);
        schema_put(&b, ch->text, text);
        schema_put_tag(&b, "</ch>");
    }
    schema_put_tag(&b, "</characteristics>");
    map[9] = b.len;
    schema_put_tag(&b, "<data></data>");
    map[10] = b.len;
    schema_put_tag(&b, "<strls></strls>");
    map[11] = b.len;
    schema_put_tag(&b, "<value_labels>");
    for (size_t i = 0; i < s->nlabels; i++)
    {
        const ctools_dta_schema_label *l = &s->labels[i];
        schema_put_tag(&b, "<lbl>");
        schema_put_le(&b, l->table_bytes, 4);
        schema_put_field(&b, l->name, CTOOLS_DTA_NAME_BYTES);
        schema_put_le(&b, 0, 3);
        schema_put(&b, l->table, l->table_bytes);
        schema_put_tag(&b, "</lbl>");
    }
    schema_put_tag(&b, "</value_labels>");
    map[12] = b.len;
    schema_put_tag(&b, "</stata_dta>");
    map[13] = b.len;
    if (b.failed)
    {
        free(b.p);
        return schema_fail(msg, cap, 920, "out of memory");
    }
    for (int i = 0; i < SCHEMA_MAP_ENTRIES; i++)
        for (int k = 0; k < 8; k++)
            b.p[map_at + (size_t)i * 8 + (size_t)k] = (unsigned char)(map[i] >> (8 * k));

    FILE *fp = ctools_fopen(path, "wb");
    if (!fp)
    {
        free(b.p);
        return schema_fail(msg, cap, 603, "cannot create the dataset description");
    }
    int failed = fwrite(b.p, 1, b.len, fp) != b.len;
    failed |= fclose(fp) != 0;
    free(b.p);
    if (failed)
    {
        ctools_remove(path);
        return schema_fail(msg, cap, 693, "cannot write the dataset description");
    }
    return 0;
}

/* ------------------------------------------------------------------------
 * Characteristics
 * ------------------------------------------------------------------------ */

static long schema_char_at(const ctools_dta_schema *s, const char *owner, const char *name)
{
    for (size_t i = 0; i < s->nchars; i++)
        if (!strcmp(s->chars[i].name, name) && !strcmp(s->chars[i].owner, owner))
            return (long)i;
    return -1;
}

const char *ctools_dta_schema_char_get(const ctools_dta_schema *s, const char *owner,
                                       const char *name)
{
    long at = schema_char_at(s, owner, name);
    return at < 0 ? NULL : s->chars[at].text;
}

void ctools_dta_schema_char_remove(ctools_dta_schema *s, const char *owner, const char *name)
{
    long at = schema_char_at(s, owner, name);
    if (at < 0)
        return;
    free(s->chars[at].text);
    memmove(&s->chars[at], &s->chars[at + 1], (s->nchars - (size_t)at - 1) * sizeof(*s->chars));
    s->nchars--;
}

int ctools_dta_schema_char_add(ctools_dta_schema *s, const char *owner, const char *name,
                               const char *text)
{
    return schema_push_char(s, owner, name, text, strlen(text));
}

int ctools_dta_schema_char_set(ctools_dta_schema *s, const char *owner, const char *name,
                               const char *text)
{
    if (!text || !*text)
    {
        ctools_dta_schema_char_remove(s, owner, name);
        return 0;
    }
    long at = schema_char_at(s, owner, name);
    if (at < 0)
        return schema_push_char(s, owner, name, text, strlen(text));
    char *copy = malloc(strlen(text) + 1);
    if (!copy)
        return 920;
    strcpy(copy, text);
    free(s->chars[at].text);
    s->chars[at].text = copy;
    return 0;
}

/* ------------------------------------------------------------------------
 * Value labels
 * ------------------------------------------------------------------------ */

long ctools_dta_schema_label_find(const ctools_dta_schema *s, const char *name)
{
    for (size_t i = 0; i < s->nlabels; i++)
        if (!strcmp(s->labels[i].name, name))
            return (long)i;
    return -1;
}

int ctools_dta_schema_label_count(const ctools_dta_schema *s, size_t at, size_t *count)
{
    *count = 0;
    if (at >= s->nlabels)
        return 610;
    const ctools_dta_schema_label *l = &s->labels[at];
    if (l->table_bytes < 8)
        return 610;
    uint64_t n = schema_le(l->table, 4), text = schema_le(l->table + 4, 4);
    if (n > (l->table_bytes - 8) / 8 || text != l->table_bytes - 8 - 8 * n)
        return 610;
    *count = (size_t)n;
    return 0;
}

int ctools_dta_schema_label_entry(const ctools_dta_schema *s, size_t at, size_t k, int32_t *value,
                                  const char **text)
{
    size_t n;
    if (ctools_dta_schema_label_count(s, at, &n) || k >= n)
        return 610;
    const ctools_dta_schema_label *l = &s->labels[at];
    const size_t text_bytes = l->table_bytes - 8 - 8 * n;
    const unsigned char *texts = l->table + 8 + 8 * n;
    uint64_t off = schema_le(l->table + 8 + 4 * k, 4);
    if (off >= text_bytes || !memchr(texts + off, 0, text_bytes - (size_t)off))
        return 610;
    *value = (int32_t)(uint32_t)schema_le(l->table + 8 + 4 * n + 4 * k, 4);
    *text = (const char *)texts + off;
    return 0;
}

typedef struct
{
    int32_t value;
    const char *text;
} schema_entry;

static int schema_cmp_entry(const void *a, const void *b)
{
    int32_t x = ((const schema_entry *)a)->value, y = ((const schema_entry *)b)->value;
    return (x > y) - (x < y);
}

int ctools_dta_schema_label_define(ctools_dta_schema *s, const char *name, size_t n,
                                   const int32_t *values, const char *const *texts)
{
    if (!name || !*name || strlen(name) >= CTOOLS_DTA_NAME_BYTES || n > (size_t)INT32_MAX / 8)
        return 198;
    if (ctools_dta_schema_label_find(s, name) >= 0)
        return 110;
    schema_entry *e = malloc((n ? n : 1) * sizeof(*e));
    if (!e)
        return 920;
    size_t text_bytes = 0;
    for (size_t k = 0; k < n; k++)
    {
        e[k] = (schema_entry){values[k], texts[k]};
        text_bytes += strlen(texts[k]) + 1;
    }
    qsort(e, n, sizeof(*e), schema_cmp_entry); /* Stata writes entries by value */
    for (size_t k = 1; k < n; k++)
        if (e[k].value == e[k - 1].value)
        {
            free(e);
            return 198;
        }
    if (text_bytes > (size_t)INT32_MAX)
    {
        free(e);
        return 198;
    }
    const size_t table_bytes = 8 + 8 * n + text_bytes;
    unsigned char *table = malloc(table_bytes);
    if (!table)
    {
        free(e);
        return 920;
    }
    schema_bytes b = {table, 0, table_bytes, 0};
    schema_put_le(&b, n, 4);
    schema_put_le(&b, text_bytes, 4);
    size_t off = 0;
    for (size_t k = 0; k < n; k++)
    {
        schema_put_le(&b, off, 4);
        off += strlen(e[k].text) + 1;
    }
    for (size_t k = 0; k < n; k++)
        schema_put_le(&b, (uint32_t)e[k].value, 4);
    for (size_t k = 0; k < n; k++)
        schema_put(&b, e[k].text, strlen(e[k].text) + 1);
    free(e);
    int rc = schema_push_label(s, name, table, table_bytes);
    free(table);
    return rc;
}

void ctools_dta_schema_default_format(uint16_t type, char format[CTOOLS_DTA_FORMAT_BYTES])
{
    switch (type)
    {
    case CTOOLS_DTA_CODE_BYTE:
    case CTOOLS_DTA_CODE_INT:
        strcpy(format, "%8.0g");
        break;
    case CTOOLS_DTA_CODE_LONG:
        strcpy(format, "%12.0g");
        break;
    case CTOOLS_DTA_CODE_FLOAT:
        strcpy(format, "%9.0g");
        break;
    case CTOOLS_DTA_CODE_DOUBLE:
        strcpy(format, "%10.0g");
        break;
    case CTOOLS_DTA_CODE_STRL:
        strcpy(format, "%9s");
        break;
    default:
        snprintf(format, CTOOLS_DTA_FORMAT_BYTES, "%%%us", (unsigned)(type < 9 ? 9 : type));
        break;
    }
}

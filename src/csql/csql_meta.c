/*
 * csql_meta.c - the description of the data in memory, written by csql.ado.
 *
 * The SPI does not report storage types, display formats or labels, so the
 * ado writes them to a file the plugin reads here. Each record is one
 * letter, a space and length-prefixed fields ("<bytes>:<text>"), so names
 * and labels may hold any bytes:
 *
 *   V name type format value-label variable-label   (in variable order)
 *   L label-name value text                          (one per entry)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "csql/csql_internal.h"
#include "ctools_file.h"

double csql_missing(void)
{
    uint64_t bits = CSQL_MISSING_BITS;
    double z;
    memcpy(&z, &bits, sizeof(z));
    return z;
}

void csql_source_free(csql_source *s)
{
    if (!s) return;
    for (size_t i = 0; i < s->nlabels; i++) {
        for (size_t k = 0; k < s->labels[i].n; k++) free(s->labels[i].texts[k]);
        free(s->labels[i].texts);
        free(s->labels[i].values);
    }
    free(s->labels);
    free(s->vars);
    free(s);
}

const csql_vlabel *csql_source_label(const csql_source *s, const char *name)
{
    for (size_t i = 0; s && i < s->nlabels; i++)
        if (strcmp(s->labels[i].name, name) == 0) return &s->labels[i];
    return NULL;
}

/* One "<len>:<bytes>" field; returns its text (not terminated) and length. */
static int meta_field(const char **p, const char *end, const char **text, size_t *len)
{
    size_t v = 0;
    const char *q = *p;
    if (q >= end || *q < '0' || *q > '9') return -1;
    while (q < end && *q >= '0' && *q <= '9') v = v * 10 + (size_t)(*q++ - '0');
    if (q >= end || *q != ':' || (size_t)(end - q - 1) < v) return -1;
    *text = q + 1;
    *len = v;
    *p = q + 1 + v;
    return 0;
}

static void meta_copy(char *dst, size_t cap, const char *src, size_t len)
{
    if (len >= cap) len = cap - 1;
    memcpy(dst, src, len);
    dst[len] = '\0';
}

uint16_t csql_meta_type_code(const char *type, size_t len)
{
    if (len == 4 && !memcmp(type, "byte", 4)) return CTOOLS_DTA_CODE_BYTE;
    if (len == 3 && !memcmp(type, "int", 3)) return CTOOLS_DTA_CODE_INT;
    if (len == 4 && !memcmp(type, "long", 4)) return CTOOLS_DTA_CODE_LONG;
    if (len == 5 && !memcmp(type, "float", 5)) return CTOOLS_DTA_CODE_FLOAT;
    if (len == 6 && !memcmp(type, "double", 6)) return CTOOLS_DTA_CODE_DOUBLE;
    if (len == 4 && !memcmp(type, "strL", 4)) return CTOOLS_DTA_CODE_STRL;
    if (len > 3 && !memcmp(type, "str", 3)) {
        unsigned v = 0;
        for (size_t i = 3; i < len; i++) {
            if (type[i] < '0' || type[i] > '9') return 0;
            v = v * 10 + (unsigned)(type[i] - '0');
            if (v > CSQL_STR_MAX) return 0;
        }
        return v >= 1 ? (uint16_t)v : 0;
    }
    return 0;
}

int csql_meta_read(const char *path, csql_source **out, char *msg, size_t cap)
{
    FILE *f = ctools_fopen(path, "rb");
    char *buf = NULL;
    long size;
    csql_source *s = NULL;
    size_t vcap = 0, lcap = 0;
    int rc = 0;
    *out = NULL;
    if (!f) {
        snprintf(msg, cap, "could not open the description of the data in memory (%s)", path);
        return 601;
    }
    if (ctools_fseek(f, 0, SEEK_END) != 0 || (size = (long)ctools_ftell(f)) < 0 ||
        ctools_fseek(f, 0, SEEK_SET) != 0) {
        fclose(f);
        snprintf(msg, cap, "could not read %s", path);
        return 692;
    }
    buf = malloc((size_t)size + 1);
    s = calloc(1, sizeof(*s));
    if (!buf || !s || fread(buf, 1, (size_t)size, f) != (size_t)size) {
        fclose(f);
        free(buf);
        free(s);
        snprintf(msg, cap, buf && s ? "could not read %s" : "out of memory", path);
        return buf && s ? 692 : 920;
    }
    fclose(f);
    snprintf(s->origin, sizeof(s->origin), "stata");

    const char *p = buf, *end = buf + size;
    while (p < end && !rc) {
        char kind = *p;
        if (kind == '\n' || kind == '\r') {
            p++;
            continue;
        }
        if (p + 2 > end || p[1] != ' ') {
            rc = 610;
            break;
        }
        p += 2;
        const char *t[5];
        size_t l[5];
        int nf = kind == 'V' ? 5 : kind == 'L' ? 3 : 0;
        if (!nf) {
            rc = 610;
            break;
        }
        for (int k = 0; k < nf && !rc; k++)
            if (meta_field(&p, end, &t[k], &l[k])) rc = 610;
        if (rc) break;
        if (kind == 'V') {
            if (s->nvars == vcap) {
                vcap = vcap ? vcap * 2 : 64;
                csql_meta_var *v = realloc(s->vars, vcap * sizeof(*v));
                if (!v) {
                    rc = 920;
                    break;
                }
                s->vars = v;
            }
            csql_meta_var *v = &s->vars[s->nvars++];
            memset(v, 0, sizeof(*v));
            meta_copy(v->name, sizeof(v->name), t[0], l[0]);
            v->type = csql_meta_type_code(t[1], l[1]);
            meta_copy(v->format, sizeof(v->format), t[2], l[2]);
            meta_copy(v->vlabel, sizeof(v->vlabel), t[3], l[3]);
            meta_copy(v->label, sizeof(v->label), t[4], l[4]);
            if (!v->type) rc = 610;
        } else {
            char name[CSQL_NAME_BYTES], num[32];
            meta_copy(name, sizeof(name), t[0], l[0]);
            meta_copy(num, sizeof(num), t[1], l[1]);
            csql_vlabel *lab = NULL;
            if (s->nlabels && strcmp(s->labels[s->nlabels - 1].name, name) == 0)
                lab = &s->labels[s->nlabels - 1];
            if (!lab) {
                if (s->nlabels == lcap) {
                    lcap = lcap ? lcap * 2 : 16;
                    csql_vlabel *x = realloc(s->labels, lcap * sizeof(*x));
                    if (!x) {
                        rc = 920;
                        break;
                    }
                    s->labels = x;
                }
                lab = &s->labels[s->nlabels++];
                memset(lab, 0, sizeof(*lab));
                snprintf(lab->name, sizeof(lab->name), "%s", name);
            }
            int32_t *vals = realloc(lab->values, (lab->n + 1) * sizeof(int32_t));
            if (vals) lab->values = vals;
            char **texts = realloc(lab->texts, (lab->n + 1) * sizeof(char *));
            if (texts) lab->texts = texts;
            char *text = malloc(l[2] + 1);
            if (!vals || !texts || !text) {
                free(text);
                rc = 920;
                break;
            }
            memcpy(text, t[2], l[2]);
            text[l[2]] = '\0';
            lab->values[lab->n] = (int32_t)strtol(num, NULL, 10);
            lab->texts[lab->n] = text;
            lab->n++;
        }
    }
    free(buf);
    if (!rc) rc = csql_sql_names(s->vars, s->nvars);
    if (rc) {
        csql_source_free(s);
        snprintf(msg, cap, rc == 920 ? "out of memory" : "the description of the data in memory is damaged");
        return rc;
    }
    *out = s;
    return 0;
}

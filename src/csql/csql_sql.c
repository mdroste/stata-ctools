/*
 * csql_sql.c - SQL text and Stata names.
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "csql/csql_internal.h"

/* ------------------------------------------------------------ buffers */

int csql_buf_add(csql_buf *b, const char *text, size_t n)
{
    if (b->len + n + 1 > b->cap) {
        size_t cap = b->cap ? b->cap : 256;
        while (cap < b->len + n + 1) cap *= 2;
        char *p = realloc(b->p, cap);
        if (!p) return 920;
        b->p = p;
        b->cap = cap;
    }
    memcpy(b->p + b->len, text, n);
    b->len += n;
    b->p[b->len] = '\0';
    return 0;
}

int csql_buf_str(csql_buf *b, const char *text)
{
    return csql_buf_add(b, text, strlen(text));
}

static int csql_buf_quoted(csql_buf *b, const char *text, char q)
{
    int rc = csql_buf_add(b, &q, 1);
    for (const char *p = text; !rc && *p; p++) {
        rc = csql_buf_add(b, p, 1);
        if (!rc && *p == q) rc = csql_buf_add(b, &q, 1);
    }
    return rc ? rc : csql_buf_add(b, &q, 1);
}

int csql_buf_ident(csql_buf *b, const char *name)
{
    return csql_buf_quoted(b, name, '"');
}

int csql_buf_literal(csql_buf *b, const char *text)
{
    return csql_buf_quoted(b, text, '\'');
}

void csql_buf_free(csql_buf *b)
{
    free(b->p);
    b->p = NULL;
    b->len = b->cap = 0;
}

/* ------------------------------------------------------------ statements */

/*
 * The parser's lexical rules that decide where a statement ends: quoted
 * strings ('', E'' with backslash escapes, $tag$ ... $tag$), quoted names
 * (""), line comments (--) and block comments (nested, as in PostgreSQL).
 */
int csql_sql_last_statement(const char *sql, size_t *start, size_t *end)
{
    size_t n = strlen(sql), i = 0, seg = 0;
    size_t best_start = 0, best_end = 0;
    int have = 0, content = 0;
    size_t content_end = 0;
    while (i <= n) {
        char c = i < n ? sql[i] : ';';
        if (c == ';') {
            if (content) {
                best_start = seg;
                best_end = content_end;
                have = 1;
            }
            seg = i + 1;
            content = 0;
            i++;
            continue;
        }
        if (c == '-' && i + 1 < n && sql[i + 1] == '-') {
            while (i < n && sql[i] != '\n') i++;
            continue;
        }
        if (c == '/' && i + 1 < n && sql[i + 1] == '*') {
            int depth = 1;
            i += 2;
            while (i < n && depth) {
                if (sql[i] == '/' && i + 1 < n && sql[i + 1] == '*') depth++, i += 2;
                else if (sql[i] == '*' && i + 1 < n && sql[i + 1] == '/') depth--, i += 2;
                else i++;
            }
            continue;
        }
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f' || c == '\v') {
            i++;
            continue;
        }
        if (!content) {
            /* the statement starts at its first token (after comments) */
            seg = i;
            content = 1;
        }
        if (c == '\'' || c == '"') {
            int backslash = c == '\'' && i > 0 && (sql[i - 1] == 'E' || sql[i - 1] == 'e') &&
                            (i < 2 || !((sql[i - 2] >= 'A' && sql[i - 2] <= 'Z') ||
                                        (sql[i - 2] >= 'a' && sql[i - 2] <= 'z') ||
                                        (sql[i - 2] >= '0' && sql[i - 2] <= '9') || sql[i - 2] == '_'));
            i++;
            while (i < n) {
                if (backslash && sql[i] == '\\' && i + 1 < n) {
                    i += 2;
                    continue;
                }
                if (sql[i] == c) {
                    if (i + 1 < n && sql[i + 1] == c) {
                        i += 2;
                        continue;
                    }
                    break;
                }
                i++;
            }
            i++;
            content_end = i > n ? n : i;
            continue;
        }
        if (c == '$') {
            /* $tag$ ... $tag$; a $ followed by digits is a parameter */
            size_t j = i + 1;
            while (j < n && (sql[j] == '_' || (sql[j] >= 'A' && sql[j] <= 'Z') ||
                             (sql[j] >= 'a' && sql[j] <= 'z') || (unsigned char)sql[j] >= 0x80 ||
                             (j > i + 1 && sql[j] >= '0' && sql[j] <= '9')))
                j++;
            if (j < n && sql[j] == '$') {
                size_t taglen = j - i + 1;
                const char *tag = sql + i;
                size_t k = j + 1;
                while (k < n && !(k + taglen <= n && memcmp(sql + k, tag, taglen) == 0)) k++;
                i = k + taglen <= n ? k + taglen : n;
                content_end = i;
                continue;
            }
        }
        i++;
        content_end = i;
    }
    if (!have) return -1;
    *start = best_start;
    *end = best_end;
    return 0;
}

/* ------------------------------------------------------------ names */

static unsigned long name_hash_ci(const char *s)
{
    unsigned long h = 1469598103UL;
    for (; *s; s++) {
        unsigned char c = (unsigned char)*s;
        if (c >= 'A' && c <= 'Z') c = (unsigned char)(c - 'A' + 'a');
        h = (h ^ c) * 16777619UL;
    }
    return h;
}

int csql_name_eq_ci(const char *a, const char *b)
{
    for (;; a++, b++) {
        unsigned char x = (unsigned char)*a, y = (unsigned char)*b;
        if (x >= 'A' && x <= 'Z') x = (unsigned char)(x - 'A' + 'a');
        if (y >= 'A' && y <= 'Z') y = (unsigned char)(y - 'A' + 'a');
        if (x != y) return 0;
        if (!x) return 1;
    }
}

int csql_sql_names(csql_meta_var *vars, size_t n)
{
    size_t cap = 16;
    while (cap < 2 * (n + 1)) cap *= 2;
    const char **slot = calloc(cap, sizeof(*slot));
    if (!slot) return 920;
    for (size_t j = 0; j < n; j++) {
        snprintf(vars[j].sqlname, sizeof(vars[j].sqlname), "%s", vars[j].name);
        for (unsigned k = 2;; k++) {
            size_t h = name_hash_ci(vars[j].sqlname) & (cap - 1);
            int taken = 0;
            while (slot[h]) {
                if (csql_name_eq_ci(slot[h], vars[j].sqlname)) {
                    taken = 1;
                    break;
                }
                h = (h + 1) & (cap - 1);
            }
            if (!taken) {
                slot[h] = vars[j].sqlname;
                break;
            }
            snprintf(vars[j].sqlname, sizeof(vars[j].sqlname), "%s_%u", vars[j].name, k);
        }
    }
    free(slot);
    return 0;
}

static const char *const csql_reserved[] = {
    "_all", "_b", "byte", "_coef", "_cons", "double", "float", "if", "in", "int", "long",
    "_n", "_N", "_pi", "_pred", "_rc", "_se", "_skip", "strL", "using", "with", NULL};

int csql_name_reserved(const char *name)
{
    for (int i = 0; csql_reserved[i]; i++)
        if (strcmp(name, csql_reserved[i]) == 0) return 1;
    if (strncmp(name, "str", 3) == 0 && name[3] >= '1' && name[3] <= '9') {
        const char *p = name + 3;
        long v = 0;
        while (*p >= '0' && *p <= '9' && v < 100000) v = v * 10 + (*p++ - '0');
        if (!*p && v >= 1 && v <= 2045) return 1;
    }
    return 0;
}

/* Decode one UTF-8 character; returns its length, 0 when invalid. */
static int csql_utf8_char(const unsigned char *s, unsigned *cp)
{
    if (s[0] < 0x80) {
        *cp = s[0];
        return 1;
    }
    if ((s[0] & 0xE0) == 0xC0 && (s[1] & 0xC0) == 0x80) {
        *cp = ((unsigned)(s[0] & 0x1F) << 6) | (s[1] & 0x3F);
        return *cp >= 0x80 ? 2 : 0;
    }
    if ((s[0] & 0xF0) == 0xE0 && (s[1] & 0xC0) == 0x80 && (s[2] & 0xC0) == 0x80) {
        *cp = ((unsigned)(s[0] & 0x0F) << 12) | ((unsigned)(s[1] & 0x3F) << 6) | (s[2] & 0x3F);
        return *cp >= 0x800 ? 3 : 0;
    }
    if ((s[0] & 0xF8) == 0xF0 && (s[1] & 0xC0) == 0x80 && (s[2] & 0xC0) == 0x80 &&
        (s[3] & 0xC0) == 0x80) {
        *cp = ((unsigned)(s[0] & 0x07) << 18) | ((unsigned)(s[1] & 0x3F) << 12) |
              ((unsigned)(s[2] & 0x3F) << 6) | (s[3] & 0x3F);
        return *cp >= 0x10000 && *cp <= 0x10FFFF ? 4 : 0;
    }
    return 0;
}

/* Letters outside ASCII that Stata accepts in names (the common scripts;
 * anything else becomes an underscore). */
static int csql_unicode_letter(unsigned cp)
{
    if (cp >= 0xC0 && cp <= 0x24F) return cp != 0xD7 && cp != 0xF7;
    return (cp >= 0x370 && cp <= 0x3FF && cp != 0x37E && cp != 0x387) || /* Greek */
           (cp >= 0x400 && cp <= 0x52F) ||                               /* Cyrillic */
           (cp >= 0x5D0 && cp <= 0x5EA) ||                               /* Hebrew */
           (cp >= 0x620 && cp <= 0x64A) ||                               /* Arabic */
           (cp >= 0x1E00 && cp <= 0x1FFF) ||                             /* Latin, Greek ext. */
           (cp >= 0x3040 && cp <= 0x30FF) ||                             /* kana */
           (cp >= 0x3400 && cp <= 0x4DBF) || (cp >= 0x4E00 && cp <= 0x9FFF) ||
           (cp >= 0xAC00 && cp <= 0xD7A3); /* Hangul */
}

static int csql_name_taken(const char *name, char (*taken)[CSQL_NAME_BYTES], size_t n)
{
    for (size_t i = 0; i < n; i++)
        if (strcmp(taken[i], name) == 0) return 1;
    return 0;
}

/* Copy at most max characters of src (valid UTF-8) into dst. */
static void csql_name_cut(char *dst, const char *src, int max)
{
    const unsigned char *s = (const unsigned char *)src;
    size_t o = 0;
    int chars = 0;
    while (*s && chars < max) {
        unsigned cp;
        int len = csql_utf8_char(s, &cp);
        if (!len) len = 1;
        if (o + (size_t)len >= CSQL_NAME_BYTES) break;
        memcpy(dst + o, s, (size_t)len);
        o += (size_t)len;
        s += len;
        chars++;
    }
    dst[o] = '\0';
}

int csql_stata_name(const char *in, char out[CSQL_NAME_BYTES], char (*taken)[CSQL_NAME_BYTES],
                    size_t n)
{
    char base[4 * 32 + 8];
    size_t o = 0;
    int chars = 0;
    const unsigned char *s = (const unsigned char *)in;
    while (*s && chars < 32) {
        unsigned cp = 0;
        int len = csql_utf8_char(s, &cp);
        if (len == 0) {
            base[o++] = '_';
            s++;
        } else if (len == 1) {
            char c = (char)cp;
            int ok = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') ||
                     c == '_';
            if (chars == 0 && c >= '0' && c <= '9') {
                base[o++] = '_';
                chars++;
                if (chars >= 32) break;
            }
            base[o++] = ok ? c : '_';
            s++;
        } else if (csql_unicode_letter(cp)) {
            memcpy(base + o, s, (size_t)len);
            o += (size_t)len;
            s += len;
        } else {
            base[o++] = '_';
            s += len;
        }
        chars++;
    }
    base[o] = '\0';
    if (o == 0) snprintf(base, sizeof(base), "v");
    if (csql_name_reserved(base)) {
        char tmp[sizeof(base) + 1];
        snprintf(tmp, sizeof(tmp), "_%s", base);
        csql_name_cut(base, tmp, 32);
    }
    if (!csql_name_taken(base, taken, n)) {
        snprintf(out, CSQL_NAME_BYTES, "%s", base);
        return strcmp(out, in) != 0;
    }
    for (unsigned long k = 1;; k++) {
        char suffix[24], cut[CSQL_NAME_BYTES], cand[CSQL_NAME_BYTES + 24];
        int slen = snprintf(suffix, sizeof(suffix), "_%lu", k);
        csql_name_cut(cut, base, 32 - slen);
        snprintf(cand, sizeof(cand), "%s%s", cut, suffix);
        if (!csql_name_taken(cand, taken, n) && !csql_name_reserved(cand)) {
            snprintf(out, CSQL_NAME_BYTES, "%s", cand);
            return 1;
        }
    }
}

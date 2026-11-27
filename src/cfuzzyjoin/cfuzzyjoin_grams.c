/*
    cfuzzyjoin_grams.c
    matchit's similarity functions: a string becomes a list of grams

    Each function follows the Mata function of the same name in matchit.ado
    1.5.2 and freqindex.ado 1.3.1 (J.D. Raffo), byte by byte: substr() and
    strlen() count bytes, tokens() splits on spaces (not tabs) and strips
    quotes, strupper() and strtrim() touch ASCII letters and spaces only, and
    strpos(s, "") is 1. Grams are listed with repetitions; the vocabulary
    turns them into counts (or presence, under matchit's unitflag).
*/

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "cfuzzyjoin_internal.h"

/* ----------------------------------------------------------------------------
 * Method names
 * -------------------------------------------------------------------------- */

static const char *skip_blank(const char *p)
{
    while (*p == ' ' || *p == '\t')
        p++;
    return p;
}

/* One argument: a number, or a string in double quotes (or bare). */
static const char *parse_arg(const char *p, char *out, size_t size)
{
    size_t n = 0;
    p = skip_blank(p);
    if (*p == '"')
    {
        p++;
        while (*p && *p != '"')
        {
            if (n + 1 < size)
                out[n++] = *p;
            p++;
        }
        if (*p == '"')
            p++;
    }
    else
    {
        while (*p && *p != ',')
        {
            if (n + 1 < size)
                out[n++] = *p;
            p++;
        }
        while (n > 0 && (out[n - 1] == ' ' || out[n - 1] == '\t'))
            n--;
    }
    out[n] = '\0';
    return skip_blank(p);
}

static int parse_number(const char *text, double *out)
{
    char *end;
    if (!*text)
        return -1;
    double v = strtod(text, &end);
    end = (char *)skip_blank(end);
    if (*end)
        return -1;
    *out = v;
    return 0;
}

static int parse_size(const char *text, int *out)
{
    double v;
    if (parse_number(text, &v) || v < 1 || v > 2045 || v != (double)(int)v)
        return -1;
    *out = (int)v;
    return 0;
}

int cfj_method_parse(const char *text, cfj_method *m, const char **what)
{
    char name[64], args[3][64];
    int nargs = 0;
    size_t n = 0;
    const char *p = skip_blank(text);
    memset(m, 0, sizeof(*m));
    *what = "similarity function not found";
    while (*p && *p != ',' && *p != ' ' && *p != '\t')
    {
        if (n + 1 >= sizeof(name))
            return -1;
        name[n++] = *p++;
    }
    name[n] = '\0';
    p = skip_blank(p);
    while (*p == ',')
    {
        if (nargs == 3)
        {
            *what = "too many arguments to the similarity function";
            return -1;
        }
        p = parse_arg(p + 1, args[nargs], sizeof(args[nargs]));
        nargs++;
    }
    if (*p)
    {
        *what = "invalid similarity function arguments";
        return -1;
    }
    if (!n || !strcmp(name, "bigram"))
        m->kind = CFJ_SIM_BIGRAM;
    else if (!strcmp(name, "ngram"))
        m->kind = CFJ_SIM_NGRAM;
    else if (!strcmp(name, "ngram_circ"))
        m->kind = CFJ_SIM_NGRAM_CIRC;
    else if (!strcmp(name, "token"))
        m->kind = CFJ_SIM_TOKEN;
    else if (!strcmp(name, "cotoken"))
        m->kind = CFJ_SIM_COTOKEN;
    else if (!strcmp(name, "scotoken"))
        m->kind = CFJ_SIM_SCOTOKEN;
    else if (!strcmp(name, "firstgram"))
        m->kind = CFJ_SIM_FIRSTGRAM;
    else if (!strcmp(name, "tkngram"))
        m->kind = CFJ_SIM_TKNGRAM;
    else if (!strcmp(name, "soundex"))
        m->kind = CFJ_SIM_SOUNDEX;
    else if (!strcmp(name, "soundex_nara"))
        m->kind = CFJ_SIM_SOUNDEX_NARA;
    else if (!strcmp(name, "soundex_fk"))
        m->kind = CFJ_SIM_SOUNDEX_FK;
    else if (!strcmp(name, "soundex_ext"))
        m->kind = CFJ_SIM_SOUNDEX_EXT;
    else if (!strcmp(name, "nysiis_fk"))
        m->kind = CFJ_SIM_NYSIIS_FK;
    else if (!strcmp(name, "token_soundex"))
        m->kind = CFJ_SIM_TOKEN_SOUNDEX;
    else if (!strcmp(name, "tokenwrap"))
        m->kind = CFJ_SIM_TOKENWRAP;
    else
        return -1;

    /* The arguments each function takes: (size) and (unitflag) */
    int sized = 0, unitarg = 0;
    switch (m->kind)
    {
    case CFJ_SIM_BIGRAM:
    case CFJ_SIM_TOKEN:
    case CFJ_SIM_TOKEN_SOUNDEX:
        unitarg = 1;
        break;
    case CFJ_SIM_NGRAM:
    case CFJ_SIM_NGRAM_CIRC:
    case CFJ_SIM_TKNGRAM:
        sized = 1;
        unitarg = 1;
        break;
    case CFJ_SIM_FIRSTGRAM:
        sized = 1;
        break;
    case CFJ_SIM_TOKENWRAP:
        if (nargs != 1)
        {
            *what = "tokenwrap needs one argument: the function applied to each token";
            return -1;
        }
        if (!strcmp(args[0], "soundex"))
            m->wrapped = CFJ_SIM_SOUNDEX;
        else if (!strcmp(args[0], "soundex_nara"))
            m->wrapped = CFJ_SIM_SOUNDEX_NARA;
        else if (!strcmp(args[0], "soundex_fk"))
            m->wrapped = CFJ_SIM_SOUNDEX_FK;
        else if (!strcmp(args[0], "soundex_ext"))
            m->wrapped = CFJ_SIM_SOUNDEX_EXT;
        else if (!strcmp(args[0], "nysiis_fk"))
            m->wrapped = CFJ_SIM_NYSIIS_FK;
        else
        {
            *what = "tokenwrap applies soundex, soundex_nara, soundex_fk, soundex_ext or nysiis_fk";
            return -1;
        }
        return 0;
    default:
        break;
    }
    int at = 0;
    if (sized)
    {
        if (nargs < 1 || parse_size(args[0], &m->n))
        {
            *what = "the similarity function needs a gram size: a positive integer";
            return -1;
        }
        at = 1;
    }
    if (nargs > at + unitarg)
    {
        *what = "too many arguments to the similarity function";
        return -1;
    }
    if (unitarg && nargs > at)
    {
        double flag;
        if (parse_number(args[at], &flag))
        {
            *what = "the unit flag of the similarity function must be a number";
            return -1;
        }
        m->unit = flag == 1;
    }
    return 0;
}

/* ----------------------------------------------------------------------------
 * Gram lists
 * -------------------------------------------------------------------------- */

void cfj_grams_init(cfj_grams *g)
{
    memset(g, 0, sizeof(*g));
}

void cfj_grams_free(cfj_grams *g)
{
    free(g->buf);
    free(g->off);
    free(g->len);
    free(g->work);
    free(g->tstart);
    free(g->tlen);
    memset(g, 0, sizeof(*g));
}

static int grow(void **p, size_t *cap, size_t need, size_t elem)
{
    if (need <= *cap)
        return 0;
    size_t c = *cap ? *cap : 16;
    while (c < need)
        c *= 2;
    void *q = realloc(*p, c * elem);
    if (!q)
        return -1;
    *p = q;
    *cap = c;
    return 0;
}

static int add(cfj_grams *g, const char *s, size_t len)
{
    if (g->n + 1 > g->ncap)
    {
        size_t cap = g->ncap ? 2 * g->ncap : 32;
        uint32_t *off = (uint32_t *)realloc(g->off, cap * sizeof(uint32_t));
        if (!off)
            return -1;
        g->off = off;
        uint32_t *l = (uint32_t *)realloc(g->len, cap * sizeof(uint32_t));
        if (!l)
            return -1;
        g->len = l;
        g->ncap = cap;
    }
    if (grow((void **)&g->buf, &g->cap, g->used + len + 1, 1))
        return -1;
    memcpy(g->buf + g->used, s, len);
    g->off[g->n] = (uint32_t)g->used;
    g->len[g->n] = (uint32_t)len;
    g->used += len;
    g->n++;
    return 0;
}

/* Two pieces joined by one space, as invtokens((a, b)) */
static int add_pair(cfj_grams *g, const char *a, size_t la, const char *b, size_t lb)
{
    if (grow((void **)&g->work, &g->wcap, la + lb + 2, 1))
        return -1;
    memcpy(g->work, a, la);
    g->work[la] = ' ';
    memcpy(g->work + la + 1, b, lb);
    return add(g, g->work, la + 1 + lb);
}

/* ----------------------------------------------------------------------------
 * Mata's tokens(s): space-separated; "..." and `"..."' are one token without
 * their quotes (an unterminated quote runs to the end); a closing quote ends
 * the token even when no space follows. Other characters, tabs included, are
 * part of a token.
 * -------------------------------------------------------------------------- */

static int tokenize(cfj_grams *g, const char *s, size_t len)
{
    size_t i = 0;
    g->ntok = 0;
    while (i < len)
    {
        while (i < len && s[i] == ' ')
            i++;
        if (i >= len)
            break;
        size_t start, end;
        if (s[i] == '"')
        {
            start = i + 1;
            size_t q = start;
            while (q < len && s[q] != '"')
                q++;
            end = q;
            i = q < len ? q + 1 : len;
        }
        else if (s[i] == '`' && i + 1 < len && s[i + 1] == '"')
        {
            int depth = 1;
            size_t q = i + 2;
            start = q;
            end = len;
            while (q < len)
            {
                if (s[q] == '`' && q + 1 < len && s[q + 1] == '"')
                {
                    depth++;
                    q += 2;
                }
                else if (s[q] == '"' && q + 1 < len && s[q + 1] == '\'')
                {
                    if (--depth == 0)
                    {
                        end = q;
                        break;
                    }
                    q += 2;
                }
                else
                    q++;
            }
            i = end < len ? end + 2 : len;
        }
        else
        {
            start = i;
            while (i < len && s[i] != ' ')
                i++;
            end = i;
        }
        if (g->ntok + 1 > g->tcap)
        {
            size_t cap = g->tcap ? 2 * g->tcap : 16;
            uint32_t *a = (uint32_t *)realloc(g->tstart, cap * sizeof(uint32_t));
            if (!a)
                return -1;
            g->tstart = a;
            uint32_t *b = (uint32_t *)realloc(g->tlen, cap * sizeof(uint32_t));
            if (!b)
                return -1;
            g->tlen = b;
            g->tcap = cap;
        }
        g->tstart[g->ntok] = (uint32_t)start;
        g->tlen[g->ntok] = (uint32_t)(end - start);
        g->ntok++;
    }
    return 0;
}

/* Bytewise comparison, as Mata compares strings */
static int str_less(const char *a, size_t la, const char *b, size_t lb)
{
    size_t n = la < lb ? la : lb;
    int c = memcmp(a, b, n);
    if (c)
        return c < 0;
    return la < lb;
}

/* ----------------------------------------------------------------------------
 * Phonetic keys
 * -------------------------------------------------------------------------- */

static int is_alpha(unsigned char c)
{
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}

static unsigned char up(unsigned char c)
{
    return (c >= 'a' && c <= 'z') ? (unsigned char)(c - 32) : c;
}

static int sx_code(unsigned char c)
{
    switch (c)
    {
    case 'B': case 'F': case 'P': case 'V':
        return 1;
    case 'C': case 'G': case 'J': case 'K': case 'Q': case 'S': case 'X': case 'Z':
        return 2;
    case 'D': case 'T':
        return 3;
    case 'L':
        return 4;
    case 'M': case 'N':
        return 5;
    case 'R':
        return 6;
    default:
        return 0;
    }
}

/* Stata's soundex(): non-letters are skipped entirely; vowels and Y (and,
 * except under NARA's rules, H and W) separate letters of the same code. */
size_t cfj_soundex(const char *s, size_t len, int nara, char *out)
{
    size_t i = 0;
    while (i < len && !is_alpha((unsigned char)s[i]))
        i++;
    if (i == len)
    {
        out[0] = '\0';
        return 0;
    }
    unsigned char first = up((unsigned char)s[i]);
    out[0] = (char)first;
    int last = sx_code(first);
    size_t n = 1;
    for (size_t j = i + 1; j < len && n < 4; j++)
    {
        unsigned char c = (unsigned char)s[j];
        if (!is_alpha(c))
            continue;
        c = up(c);
        if (c == 'A' || c == 'E' || c == 'I' || c == 'O' || c == 'U' || c == 'Y')
        {
            last = 0;
            continue;
        }
        if (c == 'H' || c == 'W')
        {
            if (!nara)
                last = 0;
            continue;
        }
        int d = sx_code(c);
        if (d != last)
            out[n++] = (char)('0' + d);
        last = d;
    }
    while (n < 4)
        out[n++] = '0';
    out[4] = '\0';
    return 4;
}

/* strupper(strtrim(s)) into w; returns the length */
static size_t upper_trim(const char *s, size_t len, char *w)
{
    size_t a = 0, b = len;
    while (a < b && s[a] == ' ')
        a++;
    while (b > a && s[b - 1] == ' ')
        b--;
    for (size_t i = a; i < b; i++)
        w[i - a] = (char)up((unsigned char)s[i]);
    return b - a;
}

static int fk_digit(unsigned char c)
{
    /* matchit's simf_soundex_fk classes (no vowels, H, W or Y) */
    return sx_code(c);
}

/* simf_soundex_fk / simf_soundex_ext. ext adds each key as it grows. */
static int soundex_fk(cfj_grams *g, const char *s, size_t len, int ext)
{
    if (grow((void **)&g->work, &g->wcap, 2 * len + 2, 1))
        return -1;
    char *u = g->work;
    size_t L = upper_trim(s, len, u);
    if (!L)
        return 0;
    char *key = g->work + len + 1;
    size_t k = 0;
    key[k++] = u[0];
    int prev = -1;
    for (size_t j = 1; j < L; j++)
    {
        unsigned char c = (unsigned char)u[j];
        if ((int)c != prev)
        {
            int d = fk_digit(c);
            if (d > 0)
            {
                key[k++] = (char)('0' + d);
                if (ext && add(g, key, k))
                    return -1;
            }
        }
        prev = c;
    }
    if (!ext)
        return add(g, key, k);
    return 0;
}

static int is_vowel(unsigned char c)
{
    return c == 'A' || c == 'E' || c == 'I' || c == 'O' || c == 'U';
}

/* simf_nysiis_fk */
static int nysiis_fk(cfj_grams *g, const char *s, size_t len)
{
    if (grow((void **)&g->work, &g->wcap, 3 * len + 8, 1))
        return -1;
    char *u = g->work;
    size_t L = upper_trim(s, len, u);
    if (!L)
        return 0;
    char first[3] = {0};
    size_t nfirst = 0;
    const char *e = u;
    size_t el = L;
    if (L >= 3 && !memcmp(u, "MAC", 3))
    {
        memcpy(first, "MC", 2), nfirst = 2;
        e += 3, el -= 3;
    }
    else if (L >= 2 && !memcmp(u, "KN", 2))
    {
        first[0] = 'N', nfirst = 1;
        e += 2, el -= 2;
    }
    else if (u[0] == 'K')
    {
        first[0] = 'C', nfirst = 1;
        e += 1, el -= 1;
    }
    else if (L >= 2 && (!memcmp(u, "PH", 2) || !memcmp(u, "PF", 2)))
    {
        memcpy(first, "FF", 2), nfirst = 2;
        e += 2, el -= 2;
    }
    else if (L >= 3 && !memcmp(u, "SCH", 3))
    {
        memcpy(first, "SS", 2), nfirst = 2;
        e += 3, el -= 3;
    }
    char lastkey = 0;
    if (el >= 2)
    {
        const char *t = e + el - 2;
        if (!memcmp(t, "EE", 2) || !memcmp(t, "IE", 2))
            lastkey = 'Y', el -= 2;
        else if (!memcmp(t, "DT", 2) || !memcmp(t, "RT", 2) || !memcmp(t, "RD", 2) ||
                 !memcmp(t, "NT", 2) || !memcmp(t, "ND", 2))
            lastkey = 'D', el -= 2;
    }
    if (!nfirst)
    {
        if (el >= 1)
        {
            first[0] = e[0], nfirst = 1;
            e += 1, el -= 1;
        }
    }
    char *fin = g->work + len + 1;
    size_t nf = 0;
    memcpy(fin, first, nfirst);
    nf = nfirst;
    char prev[3], cur[3];
    size_t nprev = nfirst, ncur = 0;
    memcpy(prev, first, nfirst);
    for (size_t j = 0; j < el; j++)
    {
        unsigned char c = (unsigned char)e[j];
        int has_next = j + 1 < el;
        unsigned char next = has_next ? (unsigned char)e[j + 1] : 0;
        if (c == 'E' && next == 'V')
        {
            memcpy(cur, "AF", 2), ncur = 2;
            j++;
        }
        else if (is_vowel(c))
            cur[0] = 'A', ncur = 1;
        else if (c == 'Q')
            cur[0] = 'G', ncur = 1;
        else if (c == 'Z')
            cur[0] = 'S', ncur = 1;
        else if (c == 'M')
            cur[0] = 'N', ncur = 1;
        else if (c == 'K')
        {
            if (next == 'N')
            {
                cur[0] = 'N', ncur = 1;
                j++;
            }
            else
                cur[0] = 'C', ncur = 1;
        }
        else if (c == 'S' && j + 2 < el && e[j + 1] == 'C' && e[j + 2] == 'H')
        {
            cur[0] = 'S', ncur = 1;
            j += 2;
        }
        else if (c == 'P' && next == 'H')
        {
            memcpy(cur, "FF", 2), ncur = 2;
            j++;
        }
        else if (c == 'H')
        {
            ncur = 0;
            if (j > 0)
                cur[0] = e[j - 1], ncur = 1;
            /* strpos("AEIOU", "") is 1: an H at the end counts as before a
             * vowel */
            if (!has_next || is_vowel(next))
                ncur = 0;
            if (j > 0 && is_vowel((unsigned char)e[j - 1]))
                ncur = 0;
        }
        else if (c == 'W')
        {
            cur[0] = 'W', ncur = 1;
            if (j > 0 && is_vowel((unsigned char)e[j - 1]))
                ncur = 0;
        }
        else if (c >= 'A' && c <= 'Z')
            cur[0] = (char)c, ncur = 1;
        /* anything else leaves curkey as it was */
        if (ncur)
        {
            if (ncur != nprev || memcmp(cur, prev, ncur))
            {
                memcpy(fin + nf, cur, ncur);
                nf += ncur;
            }
            memcpy(prev, cur, ncur);
            nprev = ncur;
        }
    }
    size_t nlast = lastkey ? 1 : 0;
    if (nprev != nlast || (nlast && prev[0] != lastkey))
    {
        if (lastkey)
            fin[nf++] = lastkey;
    }
    if (nf && fin[nf - 1] == 'S')
        nf--;
    if (nf >= 2 && fin[nf - 2] == 'A' && fin[nf - 1] == 'Y')
    {
        fin[nf - 2] = 'Y';
        nf--;
    }
    if (nf && fin[nf - 1] == 'A')
        nf--;
    return add(g, fin, nf);
}

/* One phonetic function applied to a piece of text (simf_soundex & co.) */
static int phonetic(cfj_grams *g, cfj_sim_kind kind, const char *s, size_t len)
{
    char key[8];
    switch (kind)
    {
    case CFJ_SIM_SOUNDEX:
    case CFJ_SIM_SOUNDEX_NARA:
    {
        size_t k = cfj_soundex(s, len, kind == CFJ_SIM_SOUNDEX_NARA, key);
        return add(g, key, k);
    }
    case CFJ_SIM_SOUNDEX_FK:
        return soundex_fk(g, s, len, 0);
    case CFJ_SIM_SOUNDEX_EXT:
        return soundex_fk(g, s, len, 1);
    case CFJ_SIM_NYSIIS_FK:
        return nysiis_fk(g, s, len);
    default:
        return 0;
    }
}

/* ----------------------------------------------------------------------------
 * Grams
 * -------------------------------------------------------------------------- */

/* substr windows of size n; the whole string when fewer than two fit */
static int ngrams(cfj_grams *g, const char *s, size_t len, int n, int whole)
{
    long long t = (long long)len - (n - 1);
    if (t > 1)
    {
        for (long long j = 0; j < t; j++)
            if (add(g, s + j, (size_t)n))
                return -1;
        return 0;
    }
    return whole ? add(g, s, len) : 0;
}

int cfj_extract(const cfj_method *m, const char *s, size_t len, cfj_grams *g)
{
    g->used = 0;
    g->n = 0;
    if (!s)
        s = "", len = 0;
    switch (m->kind)
    {
    case CFJ_SIM_BIGRAM:
        return ngrams(g, s, len, 2, 1);
    case CFJ_SIM_NGRAM:
        return ngrams(g, s, len, m->n, 1);
    case CFJ_SIM_NGRAM_CIRC:
    {
        long long t = (long long)len - (m->n - 1);
        if (t <= 1)
            return add(g, s, len);
        size_t n = (size_t)m->n, total = len + 1 + (n - 1);
        if (grow((void **)&g->work, &g->wcap, total + 1, 1))
            return -1;
        memcpy(g->work, s, len);
        g->work[len] = ' ';
        memcpy(g->work + len + 1, s, n - 1);
        /* add() copies into buf and leaves work alone */
        for (size_t j = 0; j + n <= total; j++)
            if (add(g, g->work + j, n))
                return -1;
        return 0;
    }
    case CFJ_SIM_TOKEN:
        if (tokenize(g, s, len))
            return -1;
        for (size_t i = 0; i < g->ntok; i++)
            if (add(g, s + g->tstart[i], g->tlen[i]))
                return -1;
        return 0;
    case CFJ_SIM_COTOKEN:
    case CFJ_SIM_SCOTOKEN:
    {
        if (tokenize(g, s, len))
            return -1;
        size_t nt = g->ntok;
        if (m->kind == CFJ_SIM_SCOTOKEN && nt == 1)
            return add(g, s + g->tstart[0], g->tlen[0]);
        for (size_t i = 1; i < nt; i++)
        {
            const char *a = s + g->tstart[i - 1], *b = s + g->tstart[i];
            size_t la = g->tlen[i - 1], lb = g->tlen[i];
            if (m->kind == CFJ_SIM_SCOTOKEN && (add(g, a, la) || add(g, b, lb)))
                return -1;
            int rc = str_less(a, la, b, lb) ? add_pair(g, a, la, b, lb) : add_pair(g, b, lb, a, la);
            if (rc)
                return -1;
        }
        return 0;
    }
    case CFJ_SIM_FIRSTGRAM:
        if (tokenize(g, s, len))
            return -1;
        for (size_t i = 0; i < g->ntok; i++)
        {
            size_t l = g->tlen[i] < (uint32_t)m->n ? g->tlen[i] : (size_t)m->n;
            if (add(g, s + g->tstart[i], l))
                return -1;
        }
        return 0;
    case CFJ_SIM_TKNGRAM:
        if (tokenize(g, s, len))
            return -1;
        for (size_t i = 0; i < g->ntok; i++)
            if (add(g, s + g->tstart[i], g->tlen[i]))
                return -1;
        return ngrams(g, s, len, m->n, 0);
    case CFJ_SIM_SOUNDEX:
    case CFJ_SIM_SOUNDEX_NARA:
    case CFJ_SIM_SOUNDEX_FK:
    case CFJ_SIM_SOUNDEX_EXT:
    case CFJ_SIM_NYSIIS_FK:
        return phonetic(g, m->kind, s, len);
    case CFJ_SIM_TOKEN_SOUNDEX:
    case CFJ_SIM_TOKENWRAP:
    {
        if (tokenize(g, s, len))
            return -1;
        cfj_sim_kind k = m->kind == CFJ_SIM_TOKEN_SOUNDEX ? CFJ_SIM_SOUNDEX : m->wrapped;
        for (size_t i = 0; i < g->ntok; i++)
            if (phonetic(g, k, s + g->tstart[i], g->tlen[i]))
                return -1;
        return 0;
    }
    }
    return 0;
}

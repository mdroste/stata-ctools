/* A Brotli decoder written from RFC 7932, for Parquet pages.
 *
 * A page's uncompressed size is known in advance, so the stream decodes
 * straight into the caller's buffer and backward distances reach into it;
 * the window size only decides which distances are static dictionary
 * references. All of RFC 7932 is read: uncompressed and metadata meta-blocks,
 * simple and complex prefix codes, block switching, context modeling,
 * NPOSTFIX/NDIRECT distance codes and the static dictionary with its 121
 * transforms. Every read is bounds-checked; malformed input yields an error,
 * never a read outside the input or a write outside the output.
 *
 * As in libbrotli (which Arrow uses), padding bits must be zero and bytes
 * after the last meta-block are ignored; the stream must produce exactly the
 * page's declared size. */
#include "parquet_internal.h"
#include <stdlib.h>
#include <string.h>

#define BBAD "Corrupt brotli compressed data."

/* brotli_dictionary[122784] */
#include "parquet_brotli_dict.inc"

/* Words of length l start at dict_offset[l]; there are 1 << dict_ndbits[l] */
static const uint32_t dict_offset[25] = {0,     0,     0,     0,      0,      4096,   9216,
                                         21504, 35840, 44032, 53248,  63488,  74752,  87040,
                                         93696, 100864, 104704, 106752, 108928, 113536, 115968,
                                         118528, 119872, 121280, 122016};
static const uint8_t dict_ndbits[25] = {0,  0,  0,  0,  10, 10, 11, 11, 10, 10, 10, 10, 10,
                                        9,  9,  8,  7,  7,  8,  7,  7,  6,  6,  5,  5};

/* ---- Bit reader (least significant bit first) ---------------------------- */
typedef struct
{
    const uint8_t *start, *p, *end;
    uint64_t val; /* unread bits, lowest first; bits above nb are zero */
    int nb;
    int bad; /* a read went past the end of the input */
} bbits;

static inline void bfill(bbits *b)
{
    while (b->nb <= 56 && b->p < b->end)
    {
        b->val |= (uint64_t)*b->p++ << b->nb;
        b->nb += 8;
    }
}

static void bpast(bbits *b)
{
    b->bad = 1;
    b->val = 0;
    b->nb = 0;
    b->p = b->end;
}

/* Reads n <= 32 bits; past the end of the input it reads zeros and sets bad */
static inline uint32_t bget(bbits *b, int n)
{
    if (b->nb < n)
    {
        bfill(b);
        if (b->nb < n)
        {
            bpast(b);
            return 0;
        }
    }
    uint32_t v = (uint32_t)(b->val & ((1ull << n) - 1));
    b->val >>= n;
    b->nb -= n;
    return v;
}

/* Skips to the next byte boundary; the skipped bits must be zero */
static int balign(bbits *b)
{
    int k = b->nb & 7;
    if (k && bget(b, k))
        return 1;
    return b->bad;
}

/* The byte offset of the next unread byte (at a byte boundary) */
static size_t bpos(const bbits *b) { return (size_t)(b->p - b->start) - (size_t)(b->nb >> 3); }

static void bseek(bbits *b, size_t at)
{
    b->p = b->start + at;
    b->val = 0;
    b->nb = 0;
}

/* ---- Prefix codes --------------------------------------------------------- */
#define ROOT 8

typedef struct
{
    uint16_t value; /* symbol, or a sub-table's offset from the table start */
    uint8_t bits;   /* code length */
    uint8_t sub;    /* index bits of the sub-table, 0 for a symbol */
} hent;

/* Tables live in one growing pool and are named by their offset in it */
typedef struct
{
    hent *e;
    size_t n, cap;
} hpool;

static hent *hreserve(hpool *pl, size_t k)
{
    if (pl->cap - pl->n < k)
    {
        size_t cap = pl->cap ? pl->cap : 8192;
        while (cap - pl->n < k)
            cap *= 2;
        hent *e = realloc(pl->e, cap * sizeof(hent));
        if (!e)
            return NULL;
        pl->e = e;
        pl->cap = cap;
    }
    return pl->e + pl->n;
}

static uint32_t reverse(uint32_t code, int len)
{
    uint32_t r = 0;
    for (int k = 0; k < len; k++, code >>= 1)
        r = (r << 1) | (code & 1);
    return r;
}

/* Builds the table of a complete canonical prefix code from code lengths
 * (0 for unused symbols). Returns 0, 1 (not a complete code) or 2 (memory). */
static int hbuild(hpool *pl, const uint8_t *len, int n, uint32_t *at)
{
    int count[16] = {0};
    for (int s = 0; s < n; s++)
        count[len[s]]++;
    count[0] = 0;
    int left = 1;
    for (int l = 1; l < 16; l++)
    {
        left = (left << 1) - count[l];
        if (left < 0)
            return 1;
    }
    if (left)
        return 1;
    uint32_t first[16], next[16], code = 0;
    for (int l = 1; l < 16; l++)
    {
        code = (code + (uint32_t)count[l - 1]) << 1;
        first[l] = code;
    }
    /* Sub-tables: one per root index of codes longer than ROOT bits, sized for
     * the longest of them */
    uint8_t subbits[1 << ROOT];
    uint16_t suboff[1 << ROOT];
    memset(subbits, 0, sizeof(subbits));
    memcpy(next, first, sizeof(next));
    for (int s = 0; s < n; s++)
        if (len[s] > ROOT)
        {
            uint32_t r = reverse(next[len[s]]++, len[s]) & ((1u << ROOT) - 1);
            if (len[s] - ROOT > subbits[r])
                subbits[r] = (uint8_t)(len[s] - ROOT);
        }
    size_t total = 1u << ROOT;
    for (int r = 0; r < (1 << ROOT); r++)
        if (subbits[r])
        {
            suboff[r] = (uint16_t)total;
            total += (size_t)1 << subbits[r];
        }
    hent *t = hreserve(pl, total);
    if (!t)
        return 2;
    for (int r = 0; r < (1 << ROOT); r++)
        if (subbits[r])
        {
            t[r].value = suboff[r];
            t[r].bits = ROOT;
            t[r].sub = subbits[r];
        }
    memcpy(next, first, sizeof(next));
    for (int s = 0; s < n; s++)
    {
        int l = len[s];
        if (!l)
            continue;
        uint32_t rv = reverse(next[l]++, l);
        hent e = {(uint16_t)s, (uint8_t)l, 0};
        if (l <= ROOT)
            for (uint32_t i = rv; i < (1u << ROOT); i += 1u << l)
                t[i] = e;
        else
        {
            uint32_t r = rv & ((1u << ROOT) - 1);
            hent *sub = t + suboff[r];
            for (uint32_t i = rv >> ROOT; i < (1u << subbits[r]); i += 1u << (l - ROOT))
                sub[i] = e;
        }
    }
    *at = (uint32_t)pl->n;
    pl->n += total;
    return 0;
}

/* A code of one symbol, which takes no bits */
static int hsingle(hpool *pl, int sym, uint32_t *at)
{
    hent *t = hreserve(pl, 1u << ROOT);
    if (!t)
        return 2;
    for (int i = 0; i < (1 << ROOT); i++)
    {
        t[i].value = (uint16_t)sym;
        t[i].bits = 0;
        t[i].sub = 0;
    }
    *at = (uint32_t)pl->n;
    pl->n += 1u << ROOT;
    return 0;
}

static inline int hread(bbits *b, const hent *t)
{
    if (b->nb < 15)
        bfill(b);
    uint32_t v = (uint32_t)b->val;
    hent e = t[v & ((1u << ROOT) - 1)];
    if (e.sub)
        e = t[e.value + ((v >> ROOT) & ((1u << e.sub) - 1))];
    if (e.bits > b->nb)
    {
        bpast(b);
        return 0;
    }
    b->val >>= e.bits;
    b->nb -= e.bits;
    return e.value;
}

/* ---- Decoder state ------------------------------------------------------- */
typedef struct
{
    int n;                  /* block types */
    uint32_t types, counts; /* tables of the block type and count codes (n >= 2) */
    int cur, prev;
    uint32_t left; /* symbols left in the current block */
} bswitch;

typedef struct
{
    bbits br;
    hpool pool;
    uint8_t *dst;
    size_t pos, out, window;
    int32_t dist[4]; /* last distances, dist[(didx - 1) & 3] the latest */
    unsigned didx;
    bswitch sw[3]; /* literals, insert-and-copy commands, distances */
    uint8_t cmode[256];
    uint8_t cmapl[64 * 256];
    uint8_t cmapd[4 * 256];
    uint32_t tl[256], ti[256], td[256];
} bdec;

/* Section 3.4 and 3.5: a simple or complex prefix code over alphabet symbols */
static const uint8_t cl_order[18] = {1, 2, 3, 4, 0, 5, 17, 6, 16, 7, 8, 9, 10, 11, 12, 13, 14, 15};
/* The fixed code of code length code lengths, indexed by the next 4 bits */
static const uint8_t cl_plen[16] = {2, 2, 2, 3, 2, 2, 2, 4, 2, 2, 2, 3, 2, 2, 2, 4};
static const uint8_t cl_pval[16] = {0, 4, 3, 2, 0, 4, 3, 1, 0, 4, 3, 2, 0, 4, 3, 5};

static int read_code(bdec *d, int alphabet, uint32_t *at)
{
    bbits *b = &d->br;
    uint8_t len[704];
    memset(len, 0, (size_t)alphabet);
    uint32_t hskip = bget(b, 2);
    if (hskip == 1)
    {
        static const uint8_t shape[5][4] = {{0}, {0}, {1, 1}, {1, 2, 2}, {2, 2, 2, 2}};
        static const uint8_t shape4[4] = {1, 2, 3, 3};
        int abits = 0;
        while ((1 << abits) < alphabet)
            abits++;
        int nsym = (int)bget(b, 2) + 1, sym[4];
        for (int i = 0; i < nsym; i++)
        {
            sym[i] = (int)bget(b, abits);
            if (sym[i] >= alphabet)
                return 1;
            for (int j = 0; j < i; j++)
                if (sym[j] == sym[i])
                    return 1;
        }
        if (b->bad)
            return 1;
        if (nsym == 1)
            return hsingle(&d->pool, sym[0], at);
        int tree = nsym == 4 ? (int)bget(b, 1) : 0;
        for (int i = 0; i < nsym; i++)
            len[sym[i]] = tree ? shape4[i] : shape[nsym][i];
        return b->bad ? 1 : hbuild(&d->pool, len, alphabet, at);
    }
    uint8_t cl[18] = {0};
    int space = 32, num = 0, only = 0;
    for (int i = (int)hskip; i < 18 && space > 0; i++)
    {
        if (b->nb < 4)
            bfill(b);
        uint32_t v = (uint32_t)b->val & 15;
        if (cl_plen[v] > b->nb)
        {
            bpast(b);
            return 1;
        }
        b->val >>= cl_plen[v];
        b->nb -= cl_plen[v];
        int c = cl_pval[v];
        cl[cl_order[i]] = (uint8_t)c;
        if (c)
        {
            space -= 32 >> c;
            num++;
            only = cl_order[i];
        }
    }
    if (b->bad || (num != 1 && space != 0))
        return 1;
    size_t mark = d->pool.n;
    uint32_t clt;
    int rc = num == 1 ? hsingle(&d->pool, only, &clt) : hbuild(&d->pool, cl, 18, &clt);
    if (rc)
        return rc;
    const hent *t = d->pool.e + clt;
    int sym = 0, prev = 8, repeat = 0, rlen = 0;
    int32_t left = 32768;
    while (sym < alphabet && left > 0)
    {
        int c = hread(b, t);
        if (b->bad)
            return 1;
        if (c < 16)
        {
            repeat = 0;
            len[sym++] = (uint8_t)c;
            if (c)
            {
                prev = c;
                left -= 32768 >> c;
            }
            continue;
        }
        /* 16 repeats the previous nonzero length, 17 repeats zero; repeat
         * codes in a row combine their counts */
        int extra = c == 16 ? 2 : 3, nl = c == 16 ? prev : 0;
        if (rlen != nl)
        {
            repeat = 0;
            rlen = nl;
        }
        int old = repeat;
        if (repeat > 0)
            repeat = (repeat - 2) << extra;
        repeat += (int)bget(b, extra) + 3;
        int delta = repeat - old;
        if (b->bad || delta > alphabet - sym)
            return 1;
        memset(len + sym, rlen, (size_t)delta);
        sym += delta;
        if (rlen)
            left -= delta << (15 - rlen);
    }
    d->pool.n = mark; /* the code length code is no longer needed */
    if (left != 0)
        return 1;
    return hbuild(&d->pool, len, alphabet, at);
}

/* NBLTYPES and NTREES: 1..256 */
static int read_count256(bbits *b)
{
    if (!bget(b, 1))
        return 1;
    int n = (int)bget(b, 3);
    return 1 + (1 << n) + (int)bget(b, n);
}

/* Block count codes (section 6) */
static const uint16_t bc_base[26] = {1,   5,   9,   13,  17,  25,  33,   41,   49,   65,   81,   97,   113,
                                     145, 177, 209, 241, 305, 369, 497, 753, 1265, 2289, 4337, 8433, 16625};
static const uint8_t bc_extra[26] = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5,
                                     5, 5, 5, 6, 6, 7, 8, 9, 10, 11, 12, 13, 24};

static uint32_t read_block_count(bbits *b, const hent *t)
{
    int c = hread(b, t);
    return bc_base[c] + bget(b, bc_extra[c]);
}

static int read_switch(bdec *d, bswitch *s)
{
    bbits *b = &d->br;
    s->n = read_count256(b);
    s->cur = 0;
    s->prev = 1;
    s->left = UINT32_MAX; /* one block type: never switched */
    if (s->n >= 2)
    {
        int rc = read_code(d, s->n + 2, &s->types);
        rc = rc ? rc : read_code(d, 26, &s->counts);
        if (rc)
            return rc;
        s->left = read_block_count(b, d->pool.e + s->counts);
    }
    return b->bad;
}

static void block_switch(bdec *d, bswitch *s)
{
    if (s->n < 2)
    {
        s->left = UINT32_MAX;
        return;
    }
    int sym = hread(&d->br, d->pool.e + s->types);
    int t = sym == 0 ? s->prev : sym == 1 ? s->cur + 1 : sym - 2;
    if (t >= s->n)
        t -= s->n;
    s->prev = s->cur;
    s->cur = t;
    s->left = read_block_count(&d->br, d->pool.e + s->counts);
}

/* Section 7.3: a context map of size entries over ntrees prefix codes */
static int read_cmap(bdec *d, int ntrees, uint8_t *map, size_t size)
{
    bbits *b = &d->br;
    if (ntrees < 2)
    {
        memset(map, 0, size);
        return b->bad;
    }
    int rlemax = bget(b, 1) ? (int)bget(b, 4) + 1 : 0;
    size_t mark = d->pool.n;
    uint32_t at;
    int rc = read_code(d, ntrees + rlemax, &at);
    if (rc)
        return rc;
    const hent *t = d->pool.e + at;
    for (size_t i = 0; i < size;)
    {
        int sym = hread(b, t);
        if (b->bad)
            return 1;
        if (sym == 0)
            map[i++] = 0;
        else if (sym <= rlemax)
        {
            size_t reps = ((size_t)1 << sym) + bget(b, sym);
            if (reps > size - i)
                return 1;
            memset(map + i, 0, reps);
            i += reps;
        }
        else
            map[i++] = (uint8_t)(sym - rlemax);
    }
    if (bget(b, 1))
    {
        /* Inverse move-to-front transform */
        uint8_t mtf[256];
        for (int k = 0; k < 256; k++)
            mtf[k] = (uint8_t)k;
        for (size_t i = 0; i < size; i++)
        {
            int idx = map[i];
            uint8_t v = mtf[idx];
            map[i] = v;
            memmove(mtf + 1, mtf, (size_t)idx);
            mtf[0] = v;
        }
    }
    d->pool.n = mark;
    for (size_t i = 0; i < size; i++)
        if (map[i] >= ntrees)
            return 1;
    return b->bad;
}

/* ---- Literal context (section 7.1) ---------------------------------------- */
/* Lut0 and Lut1 of the UTF8 mode (lut0[p1] | lut1[p2]); both tables were
 * checked against libbrotli's choice of literal code for all 65536 (p1, p2) */
static const uint8_t lut0[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 0, 0, 4, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    8, 12, 16, 12, 12, 20, 12, 16, 24, 28, 12, 12, 32, 12, 36, 12,
    44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 32, 32, 24, 40, 28, 12,
    12, 48, 52, 52, 52, 48, 52, 52, 52, 48, 52, 52, 52, 52, 52, 48,
    52, 52, 52, 52, 52, 48, 52, 52, 52, 52, 52, 24, 12, 28, 12, 12,
    12, 56, 60, 60, 60, 56, 60, 60, 60, 56, 60, 60, 60, 60, 60, 56,
    60, 60, 60, 60, 60, 56, 60, 60, 60, 60, 60, 24, 12, 28, 12, 0,
    0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
    0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
    0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
    0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
    2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3,
    2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3,
    2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3,
    2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3,
};
static const uint8_t lut1[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
    1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1,
    1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
};

/* Signed mode: 0, 1-15, 16-63, 64-127, 128-191, 192-239, 240-254, 255 */
static inline int lut2(uint8_t v)
{
    return v == 0 ? 0 : v < 16 ? 1 : v < 64 ? 2 : v < 128 ? 3 : v < 192 ? 4 : v < 240 ? 5 : v < 255 ? 6 : 7;
}

static inline int literal_context(int mode, uint8_t p1, uint8_t p2)
{
    switch (mode)
    {
    case 0: /* LSB6 */
        return p1 & 0x3f;
    case 1: /* MSB6 */
        return p1 >> 2;
    case 2: /* UTF8 */
        return lut0[p1] | lut1[p2];
    default: /* Signed */
        return (lut2(p1) << 3) | lut2(p2);
    }
}

/* ---- Static dictionary transforms (appendix B) ----------------------------- */
/* Transform types: 0 identity, 1-9 omit the last 1-9 bytes, 10 uppercase the
 * first character, 11 uppercase all, 12-20 omit the first 1-9 bytes */
#define UF 10
#define UA 11
#define OF(n) (11 + (n))
typedef struct
{
    const char *prefix;
    uint8_t type;
    const char *suffix;
} btransform;

static const btransform transforms[] = {
    {"", 0, ""},          {"", 0, " "},          {" ", 0, " "},          {"", OF(1), ""},
    {"", UF, " "},        {"", 0, " the "},      {" ", 0, ""},           {"s ", 0, " "},
    {"", 0, " of "},      {"", UF, ""},          {"", 0, " and "},       {"", OF(2), ""},
    {"", 1, ""},          {", ", 0, " "},        {"", 0, ", "},          {" ", UF, " "},
    {"", 0, " in "},      {"", 0, " to "},       {"e ", 0, " "},         {"", 0, "\""},
    {"", 0, "."},         {"", 0, "\">"},        {"", 0, "\n"},          {"", 3, ""},
    {"", 0, "]"},         {"", 0, " for "},      {"", OF(3), ""},        {"", 2, ""},
    {"", 0, " a "},       {"", 0, " that "},     {" ", UF, ""},          {"", 0, ". "},
    {".", 0, ""},         {" ", 0, ", "},        {"", OF(4), ""},        {"", 0, " with "},
    {"", 0, "'"},         {"", 0, " from "},     {"", 0, " by "},        {"", OF(5), ""},
    {"", OF(6), ""},      {" the ", 0, ""},      {"", 4, ""},            {"", 0, ". The "},
    {"", UA, ""},         {"", 0, " on "},       {"", 0, " as "},        {"", 0, " is "},
    {"", 7, ""},          {"", 1, "ing "},       {"", 0, "\n\t"},        {"", 0, ":"},
    {" ", 0, ". "},       {"", 0, "ed "},        {"", OF(9), ""},        {"", OF(7), ""},
    {"", 6, ""},          {"", 0, "("},          {"", UF, ", "},         {"", 8, ""},
    {"", 0, " at "},      {"", 0, "ly "},        {" the ", 0, " of "},   {"", 5, ""},
    {"", 9, ""},          {" ", UF, ", "},       {"", UF, "\""},         {".", 0, "("},
    {"", UA, " "},        {"", UF, "\">"},       {"", 0, "=\""},         {" ", 0, "."},
    {".com/", 0, ""},     {" the ", 0, " of the "}, {"", UF, "'"},       {"", 0, ". This "},
    {"", 0, ","},         {".", 0, " "},         {"", UF, "("},          {"", UF, "."},
    {"", 0, " not "},     {" ", 0, "=\""},       {"", 0, "er "},         {" ", UA, " "},
    {"", 0, "al "},       {" ", UA, ""},         {"", 0, "='"},          {"", UA, "\""},
    {"", UF, ". "},       {" ", 0, "("},         {"", 0, "ful "},        {" ", UF, ". "},
    {"", 0, "ive "},      {"", 0, "less "},      {"", UA, "'"},          {"", 0, "est "},
    {" ", UF, "."},       {"", UA, "\">"},       {" ", 0, "='"},         {"", UF, ","},
    {"", 0, "ize "},      {"", UA, "."},         {"\xc2\xa0", 0, ""},    {" ", 0, ","},
    {"", UF, "=\""},      {"", UA, "=\""},       {"", 0, "ous "},        {"", UA, ", "},
    {"", UF, "='"},       {" ", UF, ","},        {" ", UA, "=\""},       {" ", UA, ", "},
    {"", UA, ","},        {"", UA, "("},         {"", UA, ". "},         {" ", UA, "."},
    {"", UA, "='"},       {" ", UA, ". "},       {" ", UF, "=\""},       {" ", UA, "='"},
    {" ", UF, "='"}};
typedef char transforms_count_check[sizeof(transforms) / sizeof(transforms[0]) == 121 ? 1 : -1];

/* RFC 7932's simplified UTF-8 uppercasing of the character at p; returns
 * its length. p must have two bytes of room after the word. */
static int to_upper(uint8_t *p)
{
    if (p[0] < 0xc0)
    {
        if (p[0] >= 'a' && p[0] <= 'z')
            p[0] ^= 32;
        return 1;
    }
    if (p[0] < 0xe0)
    {
        p[1] ^= 32;
        return 2;
    }
    p[2] ^= 5;
    return 3;
}

/* Writes word idx of length len under transform t to `to`, which has room
 * bytes; returns the length written, or 0 with *bad set when it does not fit */
static size_t dict_word(uint8_t *to, size_t room, int len, uint32_t idx, int t, int *bad)
{
    const btransform *tr = &transforms[t];
    const uint8_t *w = brotli_dictionary + dict_offset[len] + (size_t)idx * (size_t)len;
    uint8_t buf[32];
    int skip = 0, n = len;
    if (tr->type >= OF(1))
        skip = tr->type - OF(0);
    else if (tr->type >= 1 && tr->type <= 9)
        n -= tr->type;
    if (n < 0)
        n = 0;
    if (skip > n)
        skip = n;
    n -= skip;
    memcpy(buf, w + skip, (size_t)n);
    memset(buf + n, 0, 4);
    if (tr->type == UF)
        to_upper(buf);
    else if (tr->type == UA)
        for (int i = 0; i < n;)
            i += to_upper(buf + i);
    size_t pl = strlen(tr->prefix), sl = strlen(tr->suffix), total = pl + (size_t)n + sl;
    if (total > room)
    {
        *bad = 1;
        return 0;
    }
    memcpy(to, tr->prefix, pl);
    memcpy(to + pl, buf, (size_t)n);
    memcpy(to + pl + n, tr->suffix, sl);
    return total;
}

/* ---- Compressed meta-blocks ------------------------------------------------ */
/* Insert-and-copy length codes (section 5) */
static const uint32_t ins_base[24] = {0,  1,  2,   3,   4,   5,   6,   8,    10,   14,   18,   26,
                                      34, 50, 66,  98,  130, 194, 322, 578,  1090, 2114, 6210, 22594};
static const uint8_t ins_extra[24] = {0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3,
                                      4, 4, 5, 5, 6, 7, 8, 9, 10, 12, 14, 24};
static const uint32_t copy_base[24] = {2,  3,  4,  5,  6,   7,   8,   9,   10,  12,   14,   18,
                                       22, 30, 38, 54, 70, 102, 134, 198, 326, 582, 1094, 2118};
static const uint8_t copy_extra[24] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2,
                                       3, 3, 4, 4, 5, 5, 6, 7, 8, 9, 10, 24};
/* The insert and copy code ranges of each group of 64 command codes */
static const uint8_t cell_ins[11] = {0, 0, 0, 0, 8, 8, 0, 16, 8, 16, 16};
static const uint8_t cell_copy[11] = {0, 8, 0, 8, 0, 8, 16, 0, 16, 8, 16};
/* Distance codes 0-15: which last distance, and the change to it */
static const uint8_t short_which[16] = {0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1};
static const int8_t short_delta[16] = {0, 0, 0, 0, -1, 1, -2, 2, -3, 3, -1, 1, -2, 2, -3, 3};

static int metablock(bdec *d, size_t mlen)
{
    bbits *b = &d->br;
    int rc;
    d->pool.n = 0;
    for (int k = 0; k < 3; k++)
        if ((rc = read_switch(d, &d->sw[k])))
            return rc;
    int npostfix = (int)bget(b, 2);
    int ndirect = (int)bget(b, 4) << npostfix;
    for (int i = 0; i < d->sw[0].n; i++)
        d->cmode[i] = (uint8_t)bget(b, 2);
    int ntl = read_count256(b);
    if ((rc = read_cmap(d, ntl, d->cmapl, (size_t)64 * (size_t)d->sw[0].n)))
        return rc;
    int ntd = read_count256(b);
    if ((rc = read_cmap(d, ntd, d->cmapd, (size_t)4 * (size_t)d->sw[2].n)))
        return rc;
    for (int i = 0; i < ntl; i++)
        if ((rc = read_code(d, 256, &d->tl[i])))
            return rc;
    for (int i = 0; i < d->sw[1].n; i++)
        if ((rc = read_code(d, 704, &d->ti[i])))
            return rc;
    int dalpha = 16 + ndirect + (48 << npostfix);
    for (int i = 0; i < ntd; i++)
        if ((rc = read_code(d, dalpha, &d->td[i])))
            return rc;
    if (b->bad)
        return 1;

    const hent *T = d->pool.e; /* the pool no longer grows */
    uint8_t *dst = d->dst;
    size_t pos = d->pos, left = mlen;
    bswitch *L = &d->sw[0], *I = &d->sw[1], *D = &d->sw[2];
    int mode = d->cmode[0];
    const uint8_t *lmap = d->cmapl;
    const uint32_t postmask = (1u << npostfix) - 1;
    while (left)
    {
        if (b->bad)
            return 1;
        if (!I->left)
            block_switch(d, I);
        I->left--;
        int cmd = hread(b, T + d->ti[I->cur]);
        int cell = cmd >> 6;
        int ic = cell_ins[cell] + ((cmd >> 3) & 7), cc = cell_copy[cell] + (cmd & 7);
        size_t ins = ins_base[ic] + bget(b, ins_extra[ic]);
        size_t copy = copy_base[cc] + bget(b, copy_extra[cc]);
        if (ins > left)
            return 1;
        left -= ins;
        for (size_t k = 0; k < ins; k++)
        {
            if (!L->left)
            {
                block_switch(d, L);
                mode = d->cmode[L->cur];
                lmap = d->cmapl + 64 * L->cur;
            }
            L->left--;
            uint8_t p1 = pos ? dst[pos - 1] : 0, p2 = pos > 1 ? dst[pos - 2] : 0;
            dst[pos++] = (uint8_t)hread(b, T + d->tl[lmap[literal_context(mode, p1, p2)]]);
            if (b->bad)
                return 1;
        }
        /* The meta-block may end after a command's literals */
        if (!left)
            break;
        int dcode = 0;
        if (cmd >= 128)
        {
            if (!D->left)
                block_switch(d, D);
            D->left--;
            int dctx = copy > 4 ? 3 : (int)copy - 2;
            dcode = hread(b, T + d->td[d->cmapd[4 * D->cur + dctx]]);
        }
        int64_t distance;
        if (dcode < 16)
        {
            distance = (int64_t)d->dist[(d->didx - 1u - short_which[dcode]) & 3] + short_delta[dcode];
            if (distance <= 0)
                return 1;
        }
        else if (dcode < 16 + ndirect)
            distance = dcode - 15;
        else
        {
            int x = dcode - ndirect - 16;
            int nbits = 1 + (x >> (npostfix + 1));
            int64_t offset = ((int64_t)(2 + ((x >> npostfix) & 1)) << nbits) - 4;
            distance = ((offset + bget(b, nbits)) << npostfix) + (x & postmask) + ndirect + 1;
        }
        if (b->bad)
            return 1;
        size_t maxd = pos < d->window ? pos : d->window;
        if ((uint64_t)distance > maxd)
        {
            /* A static dictionary word; the last distances are unchanged */
            if (copy < 4 || copy > 24)
                return 1;
            uint64_t id = (uint64_t)distance - maxd - 1;
            int nb = dict_ndbits[copy];
            if ((id >> nb) >= 121)
                return 1;
            int bad = 0;
            size_t n = dict_word(dst + pos, left, (int)copy, (uint32_t)(id & ((1u << nb) - 1)),
                                 (int)(id >> nb), &bad);
            if (bad)
                return 1;
            pos += n;
            left -= n;
            continue;
        }
        if (copy > left)
            return 1;
        const uint8_t *from = dst + pos - (size_t)distance;
        uint8_t *to = dst + pos;
        if ((size_t)distance >= copy)
            memcpy(to, from, copy);
        else
            for (size_t k = 0; k < copy; k++)
                to[k] = from[k];
        pos += copy;
        left -= copy;
        if (dcode)
        {
            d->dist[d->didx & 3] = (int32_t)distance;
            d->didx++;
        }
    }
    d->pos = pos;
    return b->bad;
}

/* ---- Stream ---------------------------------------------------------------- */
static int stream(bdec *d, size_t n)
{
    bbits *b = &d->br;
    int wbits = 16;
    if (bget(b, 1))
    {
        int k = (int)bget(b, 3);
        if (k)
            wbits = 17 + k;
        else
        {
            k = (int)bget(b, 3);
            if (k == 1) /* the large-window marker is not RFC 7932 */
                return 1;
            wbits = k ? 8 + k : 17;
        }
    }
    d->window = ((size_t)1 << wbits) - 16;
    for (;;)
    {
        int last = (int)bget(b, 1);
        if (last && bget(b, 1)) /* ISLASTEMPTY */
            break;
        int nibbles = (int)bget(b, 2);
        if (b->bad)
            return 1;
        if (nibbles == 3)
        {
            /* Metadata, skipped; libbrotli also accepts it as the last meta-block */
            if (bget(b, 1))
                return 1;
            int bytes = (int)bget(b, 2);
            size_t skip = 0;
            for (int k = 0; k < bytes; k++)
            {
                uint32_t v = bget(b, 8);
                if (k == bytes - 1 && bytes > 1 && !v)
                    return 1;
                skip |= (size_t)v << (8 * k);
            }
            if (bytes)
                skip++;
            if (balign(b))
                return 1;
            size_t at = bpos(b);
            if (skip > n - at)
                return 1;
            bseek(b, at + skip);
            if (last)
                break;
            continue;
        }
        nibbles += 4;
        size_t mlen = 0;
        for (int k = 0; k < nibbles; k++)
        {
            uint32_t v = bget(b, 4);
            if (k == nibbles - 1 && nibbles > 4 && !v)
                return 1;
            mlen |= (size_t)v << (4 * k);
        }
        mlen++;
        int raw = last ? 0 : (int)bget(b, 1);
        if (b->bad || mlen > d->out - d->pos)
            return 1;
        if (raw)
        {
            if (balign(b))
                return 1;
            size_t at = bpos(b);
            if (mlen > n - at)
                return 1;
            memcpy(d->dst + d->pos, b->start + at, mlen);
            d->pos += mlen;
            bseek(b, at + mlen);
        }
        else
        {
            int rc = metablock(d, mlen);
            if (rc)
                return rc;
        }
        if (last)
            break;
    }
    if (balign(b))
        return 1;
    return d->pos == d->out ? 0 : 1;
}

const char *pq_brotli_decompress(const uint8_t *src, size_t n, uint8_t *dst, size_t out)
{
    bdec *d = malloc(sizeof(bdec));
    if (!d)
        return "Brotli decompression failed: out of memory";
    d->br.start = d->br.p = src;
    d->br.end = src + n;
    d->br.val = 0;
    d->br.nb = 0;
    d->br.bad = 0;
    d->pool.e = NULL;
    d->pool.n = d->pool.cap = 0;
    d->dst = dst;
    d->pos = 0;
    d->out = out;
    d->dist[0] = 16;
    d->dist[1] = 15;
    d->dist[2] = 11;
    d->dist[3] = 4;
    d->didx = 0;
    int rc = stream(d, n);
    free(d->pool.e);
    free(d);
    if (rc == 2)
        return "Brotli decompression failed: out of memory";
    return rc ? BBAD : NULL;
}

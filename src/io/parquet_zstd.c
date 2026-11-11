/* A compact Zstandard decoder written from RFC 8878, for Parquet pages.
 *
 * A page's uncompressed size is known in advance, so frames decode straight
 * into the caller's buffer and matches may reach back to any earlier byte of
 * it. Dictionaries are not supported (Parquet writers do not use them).
 * Every read is bounds-checked; malformed input yields an error, never a read
 * outside the input or a write outside the output. */
#include "parquet_internal.h"
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define ZBAD "Corrupt ZSTD compressed data."

/* ---- Backward bit stream (FSE and Huffman streams) ---------------------- */
typedef struct
{
    const uint8_t *s;
    size_t n;
    int64_t bit; /* bits left to read; negative once past the start */
} zback;

static int zback_init(zback *b, const uint8_t *s, size_t n)
{
    if (n == 0 || s[n - 1] == 0)
        return 1;
    int top = 7;
    while (!(s[n - 1] >> top))
        top--;
    b->s = s;
    b->n = n;
    b->bit = (int64_t)(n - 1) * 8 + top;
    return 0;
}

/* Reads n (<= 57) bits; bits before the stream start read as zero */
static uint64_t zback_read(zback *b, int n)
{
    if (n == 0)
        return 0;
    b->bit -= n;
    int64_t lo = b->bit;
    uint64_t v = 0;
    for (int k = n - 1; k >= 0; k--)
    {
        int64_t at = lo + k;
        int bitv = at >= 0 ? (b->s[at >> 3] >> (at & 7)) & 1 : 0;
        v = (v << 1) | (uint64_t)bitv;
    }
    return v;
}

/* Faster path when the requested bits are all inside the stream */
static inline uint64_t zback_bits(zback *b, int n)
{
    if (n == 0)
        return 0;
    if (b->bit - n < 0 || n > 56)
        return zback_read(b, n);
    b->bit -= n;
    int64_t lo = b->bit;
    size_t byte = (size_t)(lo >> 3);
    int shift = (int)(lo & 7);
    uint64_t w = 0;
    if (byte + 8 <= b->n)
    {
        /* One little-endian word: every bit asked for lies in it (n <= 56) */
        const uint8_t *q = b->s + byte;
        w = (uint64_t)q[0] | (uint64_t)q[1] << 8 | (uint64_t)q[2] << 16 | (uint64_t)q[3] << 24 |
            (uint64_t)q[4] << 32 | (uint64_t)q[5] << 40 | (uint64_t)q[6] << 48 |
            (uint64_t)q[7] << 56;
    }
    else
    {
        int need = (shift + n + 7) >> 3;
        for (int k = 0; k < need; k++)
            w |= (uint64_t)b->s[byte + k] << (8 * k);
    }
    return (w >> shift) & ((1ull << n) - 1);
}

/* ---- FSE ----------------------------------------------------------------- */
typedef struct
{
    int log;
    uint16_t size;
    uint8_t sym[512];
    uint8_t nbits[512];
    uint16_t base[512];
} zfse;

static int highbit(uint32_t v)
{
    int n = -1;
    while (v)
    {
        v >>= 1;
        n++;
    }
    return n;
}

static int fse_build(zfse *t, const int16_t *norm, int nsym, int log)
{
    if (log > 9 || nsym > 256)
        return 1;
    int size = 1 << log, high = size - 1;
    uint16_t next[256];
    t->log = log;
    t->size = (uint16_t)size;
    for (int s = 0; s < nsym; s++)
    {
        if (norm[s] == -1)
        {
            t->sym[high--] = (uint8_t)s;
            next[s] = 1;
        }
        else
            next[s] = (uint16_t)norm[s];
    }
    int step = (size >> 1) + (size >> 3) + 3, mask = size - 1, pos = 0;
    for (int s = 0; s < nsym; s++)
        for (int i = 0; i < norm[s]; i++)
        {
            t->sym[pos] = (uint8_t)s;
            do
                pos = (pos + step) & mask;
            while (pos > high);
        }
    if (pos != 0)
        return 1;
    for (int i = 0; i < size; i++)
    {
        int s = t->sym[i];
        uint32_t ns = next[s]++;
        if (ns == 0)
            return 1;
        int nb = log - highbit(ns);
        t->nbits[i] = (uint8_t)nb;
        t->base[i] = (uint16_t)((ns << nb) - (uint32_t)size);
    }
    return 0;
}

/* Forward little-endian bit reader for FSE table headers */
typedef struct
{
    const uint8_t *s;
    size_t n, bit;
} zfwd;

static int zfwd_read(zfwd *r, int n, uint32_t *v)
{
    uint32_t out = 0;
    for (int k = 0; k < n; k++, r->bit++)
    {
        size_t byte = r->bit >> 3;
        if (byte >= r->n)
            return 1;
        out |= (uint32_t)((r->s[byte] >> (r->bit & 7)) & 1) << k;
    }
    *v = out;
    return 0;
}

/* Reads a table description; returns bytes consumed or 0 on error */
static size_t fse_header(zfse *t, const uint8_t *s, size_t n, int maxsym, int maxlog)
{
    zfwd r = {s, n, 0};
    uint32_t v;
    if (zfwd_read(&r, 4, &v))
        return 0;
    int log = (int)v + 5;
    if (log > maxlog)
        return 0;
    int remaining = (1 << log) + 1, sym = 0;
    int16_t norm[256];
    memset(norm, 0, sizeof(norm));
    while (remaining > 1 && sym <= maxsym)
    {
        int bits = highbit((uint32_t)remaining) + 1;
        if (zfwd_read(&r, bits, &v))
            return 0;
        uint32_t lower = (1u << (bits - 1)) - 1, threshold = (1u << bits) - 1 - (uint32_t)remaining;
        if ((v & lower) < threshold)
        {
            r.bit--;
            v &= lower;
        }
        else if (v > lower)
            v -= threshold;
        int p = (int)v - 1;
        remaining -= p < 0 ? -p : p;
        norm[sym++] = (int16_t)p;
        if (p == 0)
        {
            uint32_t rep;
            for (;;)
            {
                if (zfwd_read(&r, 2, &rep))
                    return 0;
                for (uint32_t i = 0; i < rep; i++)
                {
                    if (sym > maxsym)
                        return 0;
                    norm[sym++] = 0;
                }
                if (rep != 3)
                    break;
            }
        }
    }
    if (remaining != 1 || sym > maxsym + 1)
        return 0;
    if (fse_build(t, norm, sym, log))
        return 0;
    return (r.bit + 7) >> 3;
}

static void fse_rle(zfse *t, uint8_t sym)
{
    t->log = 0;
    t->size = 1;
    t->sym[0] = sym;
    t->nbits[0] = 0;
    t->base[0] = 0;
}

/* ---- Huffman literals ------------------------------------------------------ */
typedef struct
{
    int maxbits;
    uint8_t sym[1 << 11];
    uint8_t nbits[1 << 11];
} zhuf;

static int huf_from_weights(zhuf *h, const uint8_t *w, int n)
{
    if (n < 1 || n > 255)
        return 1;
    uint32_t total = 0;
    for (int i = 0; i < n; i++)
    {
        if (w[i] > 11)
            return 1;
        if (w[i])
            total += 1u << (w[i] - 1);
    }
    if (total == 0)
        return 1;
    int maxbits = highbit(total) + 1;
    if (maxbits > 11)
        return 1;
    uint32_t rest = (1u << maxbits) - total;
    if (rest & (rest - 1))
        return 1;
    uint8_t bits[256];
    memset(bits, 0, sizeof(bits));
    for (int i = 0; i < n; i++)
        bits[i] = w[i] ? (uint8_t)(maxbits + 1 - w[i]) : 0;
    bits[n] = (uint8_t)(maxbits + 1 - (highbit(rest) + 1));
    int count = n + 1;
    int rank[13] = {0}, idx[13];
    for (int i = 0; i < count; i++)
        rank[bits[i]]++;
    idx[maxbits] = 0;
    for (int b = maxbits; b >= 1; b--)
    {
        idx[b - 1] = idx[b] + rank[b] * (1 << (maxbits - b));
        for (int k = idx[b]; k < idx[b - 1]; k++)
            h->nbits[k] = (uint8_t)b;
    }
    if (idx[0] != (1 << maxbits))
        return 1;
    for (int i = 0; i < count; i++)
        if (bits[i])
        {
            int len = 1 << (maxbits - bits[i]);
            memset(&h->sym[idx[bits[i]]], i, (size_t)len);
            idx[bits[i]] += len;
        }
    h->maxbits = maxbits;
    return 0;
}

/* Tree description; returns bytes consumed or 0 */
static size_t huf_header(zhuf *h, const uint8_t *s, size_t n)
{
    if (n < 1)
        return 0;
    uint8_t w[256];
    int count = 0;
    size_t used;
    if (s[0] >= 128)
    {
        count = s[0] - 127;
        used = 1 + (size_t)(count + 1) / 2;
        if (used > n)
            return 0;
        for (int i = 0; i < count; i++)
            w[i] = (uint8_t)(i & 1 ? s[1 + i / 2] & 15 : s[1 + i / 2] >> 4);
    }
    else
    {
        size_t size = s[0];
        if (size == 0 || 1 + size > n)
            return 0;
        zfse local;
        size_t hdr = fse_header(&local, s + 1, size, 255, 6);
        if (!hdr || hdr >= size)
            return 0;
        zback b;
        if (zback_init(&b, s + 1 + hdr, size - hdr))
            return 0;
        uint32_t s1 = (uint32_t)zback_bits(&b, local.log), s2 = (uint32_t)zback_bits(&b, local.log);
        for (;;)
        {
            if (count >= 255)
                return 0;
            w[count++] = local.sym[s1];
            s1 = local.base[s1] + (uint32_t)zback_bits(&b, local.nbits[s1]);
            if (b.bit < 0)
            {
                if (count >= 255)
                    return 0;
                w[count++] = local.sym[s2];
                break;
            }
            if (count >= 255)
                return 0;
            w[count++] = local.sym[s2];
            s2 = local.base[s2] + (uint32_t)zback_bits(&b, local.nbits[s2]);
            if (b.bit < 0)
            {
                if (count >= 255)
                    return 0;
                w[count++] = local.sym[s1];
                break;
            }
        }
        used = 1 + size;
    }
    if (huf_from_weights(h, w, count))
        return 0;
    return used;
}

static int huf_stream(const zhuf *h, const uint8_t *s, size_t n, uint8_t *out, size_t count)
{
    zback b;
    if (zback_init(&b, s, n))
        return 1;
    int mb = h->maxbits;
    uint32_t mask = (1u << mb) - 1;
    uint32_t state = (uint32_t)zback_bits(&b, mb);
    for (size_t i = 0; i < count; i++)
    {
        out[i] = h->sym[state];
        int nb = h->nbits[state];
        state = ((state << nb) | (uint32_t)zback_bits(&b, nb)) & mask;
    }
    return b.bit != -mb;
}

/* ---- Sequences ------------------------------------------------------------ */
static const int16_t LL_DEF[36] = {4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2,
                                   2, 2, 2, 2, 2, 2, 2, 3, 2, 1, 1, 1, 1, 1, -1, -1, -1, -1};
static const int16_t ML_DEF[53] = {1, 4, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1};
static const int16_t OF_DEF[29] = {1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1,
                                   1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1};
static const uint32_t LL_BASE[36] = {0,  1,  2,   3,   4,   5,    6,    7,    8,    9,     10,    11,
                                     12, 13, 14,  15,  16,  18,   20,   22,   24,   28,    32,    40,
                                     48, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536};
static const uint8_t LL_BITS[36] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  1,
                                    1, 1, 2, 2, 3, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
static const uint32_t ML_BASE[53] = {
    3,  4,  5,  6,  7,  8,  9,  10, 11, 12,  13,  14,  15,  16,   17,   18,   19,   20,
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,  31,  32,  33,  34,   35,   37,   39,   41,
    43, 47, 51, 59, 67, 83, 99, 131, 259, 515, 1027, 2051, 4099, 8195, 16387, 32771, 65539};
static const uint8_t ML_BITS[53] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                                    2, 2, 3, 3, 4, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

typedef struct
{
    zfse ll, of, ml;
    int have_ll, have_of, have_ml;
    zhuf huf;
    int have_huf;
    uint32_t rep[3];
    uint8_t lit[1 << 17];
} zstate;

static size_t seq_table(zfse *t, int *have, int mode, const uint8_t *s, size_t n, const int16_t *def,
                        int ndef, int deflog, int maxsym, int maxlog)
{
    if (mode == 0)
    {
        if (fse_build(t, def, ndef, deflog))
            return (size_t)-1;
        *have = 1;
        return 0;
    }
    if (mode == 1)
    {
        if (n < 1 || s[0] > maxsym)
            return (size_t)-1;
        fse_rle(t, s[0]);
        *have = 1;
        return 1;
    }
    if (mode == 2)
    {
        size_t used = fse_header(t, s, n, maxsym, maxlog);
        if (!used)
            return (size_t)-1;
        *have = 1;
        return used;
    }
    return *have ? 0 : (size_t)-1;
}

static const char *block(zstate *z, const uint8_t *s, size_t n, uint8_t *dst, size_t *pos, size_t out)
{
    /* Literals section */
    if (n < 1)
        return ZBAD;
    int ltype = s[0] & 3, sf = (s[0] >> 2) & 3;
    size_t regen, csize = 0, hdr;
    const uint8_t *lit;
    if (ltype <= 1)
    {
        if (sf == 0 || sf == 2)
        {
            regen = s[0] >> 3;
            hdr = 1;
        }
        else if (sf == 1)
        {
            if (n < 2)
                return ZBAD;
            regen = (s[0] >> 4) | (size_t)s[1] << 4;
            hdr = 2;
        }
        else
        {
            if (n < 3)
                return ZBAD;
            regen = (s[0] >> 4) | (size_t)s[1] << 4 | (size_t)s[2] << 12;
            hdr = 3;
        }
        if (regen > sizeof(z->lit))
            return ZBAD;
        if (ltype == 0)
        {
            if (n - hdr < regen)
                return ZBAD;
            lit = s + hdr;
            hdr += regen;
        }
        else
        {
            if (n - hdr < 1)
                return ZBAD;
            memset(z->lit, s[hdr], regen);
            lit = z->lit;
            hdr += 1;
        }
    }
    else
    {
        int streams = sf == 0 ? 1 : 4;
        if (sf <= 1)
        {
            if (n < 3)
                return ZBAD;
            uint32_t v = s[0] | (uint32_t)s[1] << 8 | (uint32_t)s[2] << 16;
            regen = (v >> 4) & 0x3ff;
            csize = (v >> 14) & 0x3ff;
            hdr = 3;
        }
        else if (sf == 2)
        {
            if (n < 4)
                return ZBAD;
            uint32_t v = s[0] | (uint32_t)s[1] << 8 | (uint32_t)s[2] << 16 | (uint32_t)s[3] << 24;
            regen = (v >> 4) & 0x3fff;
            csize = (v >> 18) & 0x3fff;
            hdr = 4;
        }
        else
        {
            if (n < 5)
                return ZBAD;
            uint64_t v = s[0] | (uint64_t)s[1] << 8 | (uint64_t)s[2] << 16 | (uint64_t)s[3] << 24 |
                         (uint64_t)s[4] << 32;
            regen = (size_t)((v >> 4) & 0x3ffff);
            csize = (size_t)((v >> 22) & 0x3ffff);
            hdr = 5;
        }
        if (regen > sizeof(z->lit) || n - hdr < csize)
            return ZBAD;
        const uint8_t *c = s + hdr;
        size_t cn = csize;
        if (ltype == 2)
        {
            size_t used = huf_header(&z->huf, c, cn);
            if (!used)
                return ZBAD;
            z->have_huf = 1;
            c += used;
            cn -= used;
        }
        else if (!z->have_huf)
            return ZBAD;
        if (streams == 1)
        {
            if (huf_stream(&z->huf, c, cn, z->lit, regen))
                return ZBAD;
        }
        else
        {
            if (cn < 6)
                return ZBAD;
            size_t l1 = c[0] | (size_t)c[1] << 8, l2 = c[2] | (size_t)c[3] << 8,
                   l3 = c[4] | (size_t)c[5] << 8;
            if (l1 + l2 + l3 > cn - 6)
                return ZBAD;
            size_t l4 = cn - 6 - l1 - l2 - l3, seg = (regen + 3) / 4;
            if (regen < 3 * seg)
                return ZBAD;
            const uint8_t *p = c + 6;
            if (huf_stream(&z->huf, p, l1, z->lit, seg) ||
                huf_stream(&z->huf, p + l1, l2, z->lit + seg, seg) ||
                huf_stream(&z->huf, p + l1 + l2, l3, z->lit + 2 * seg, seg) ||
                huf_stream(&z->huf, p + l1 + l2 + l3, l4, z->lit + 3 * seg, regen - 3 * seg))
                return ZBAD;
        }
        lit = z->lit;
        hdr += csize;
    }
    s += hdr;
    n -= hdr;

    /* Sequences section */
    if (n < 1)
        return ZBAD;
    size_t nseq = s[0], used = 1;
    if (nseq >= 128)
    {
        if (nseq < 255)
        {
            if (n < 2)
                return ZBAD;
            nseq = ((nseq - 128) << 8) + s[1];
            used = 2;
        }
        else
        {
            if (n < 3)
                return ZBAD;
            nseq = s[1] + ((size_t)s[2] << 8) + 0x7f00;
            used = 3;
        }
    }
    s += used;
    n -= used;
    size_t lp = 0;
    if (nseq)
    {
        if (n < 1)
            return ZBAD;
        int modes = s[0];
        if (modes & 3)
            return ZBAD;
        s++;
        n--;
        size_t k = seq_table(&z->ll, &z->have_ll, (modes >> 6) & 3, s, n, LL_DEF, 36, 6, 35, 9);
        if (k == (size_t)-1)
            return ZBAD;
        s += k;
        n -= k;
        k = seq_table(&z->of, &z->have_of, (modes >> 4) & 3, s, n, OF_DEF, 29, 5, 31, 8);
        if (k == (size_t)-1)
            return ZBAD;
        s += k;
        n -= k;
        k = seq_table(&z->ml, &z->have_ml, (modes >> 2) & 3, s, n, ML_DEF, 53, 6, 52, 9);
        if (k == (size_t)-1)
            return ZBAD;
        s += k;
        n -= k;
        zback b;
        if (zback_init(&b, s, n))
            return ZBAD;
        uint32_t sl = (uint32_t)zback_bits(&b, z->ll.log), so = (uint32_t)zback_bits(&b, z->of.log),
                 sm = (uint32_t)zback_bits(&b, z->ml.log);
        for (size_t q = 0; q < nseq; q++)
        {
            int lc = z->ll.sym[sl], oc = z->of.sym[so], mc = z->ml.sym[sm];
            if (lc > 35 || mc > 52 || oc > 31)
                return ZBAD;
            uint64_t ofv = (1ull << oc) + zback_bits(&b, oc);
            size_t ml = ML_BASE[mc] + (size_t)zback_bits(&b, ML_BITS[mc]);
            size_t ll = LL_BASE[lc] + (size_t)zback_bits(&b, LL_BITS[lc]);
            uint64_t offset;
            if (ofv > 3)
            {
                offset = ofv - 3;
                z->rep[2] = z->rep[1];
                z->rep[1] = z->rep[0];
                z->rep[0] = (uint32_t)offset;
            }
            else
            {
                int idx = (int)ofv - 1 + (ll == 0);
                if (idx == 0)
                    offset = z->rep[0];
                else
                {
                    offset = idx == 3 ? (uint64_t)z->rep[0] - 1 : z->rep[idx];
                    if (idx > 1)
                        z->rep[2] = z->rep[1];
                    z->rep[1] = z->rep[0];
                    z->rep[0] = (uint32_t)offset;
                }
            }
            if (q + 1 < nseq)
            {
                sl = z->ll.base[sl] + (uint32_t)zback_bits(&b, z->ll.nbits[sl]);
                sm = z->ml.base[sm] + (uint32_t)zback_bits(&b, z->ml.nbits[sm]);
                so = z->of.base[so] + (uint32_t)zback_bits(&b, z->of.nbits[so]);
            }
            if (ll > regen - lp || ll > out - *pos)
                return ZBAD;
            memcpy(dst + *pos, lit + lp, ll);
            lp += ll;
            *pos += ll;
            if (offset == 0 || offset > *pos || ml > out - *pos)
                return ZBAD;
            uint8_t *d = dst + *pos;
            const uint8_t *from = d - offset;
            if (offset >= ml)
                memcpy(d, from, ml);
            else
                for (size_t i = 0; i < ml; i++)
                    d[i] = from[i];
            *pos += ml;
        }
        if (b.bit != 0)
            return ZBAD;
    }
    else if (n != 0)
        return ZBAD;
    if (regen - lp > out - *pos)
        return ZBAD;
    memcpy(dst + *pos, lit + lp, regen - lp);
    *pos += regen - lp;
    return NULL;
}

const char *pq_zstd_decompress(const uint8_t *s, size_t n, uint8_t *dst, size_t out)
{
    zstate *z = malloc(sizeof(zstate));
    if (!z)
        return "ZSTD decompression failed: out of memory";
    size_t pos = 0;
    const char *message = NULL;
    while (n > 0 && !message)
    {
        if (n < 4)
        {
            message = ZBAD;
            break;
        }
        uint32_t magic = s[0] | (uint32_t)s[1] << 8 | (uint32_t)s[2] << 16 | (uint32_t)s[3] << 24;
        if ((magic & 0xfffffff0u) == 0x184d2a50u)
        {
            if (n < 8)
            {
                message = ZBAD;
                break;
            }
            uint32_t skip = s[4] | (uint32_t)s[5] << 8 | (uint32_t)s[6] << 16 | (uint32_t)s[7] << 24;
            if (skip > n - 8)
            {
                message = ZBAD;
                break;
            }
            s += 8 + skip;
            n -= 8 + skip;
            continue;
        }
        if (magic != 0xfd2fb528u || n < 5)
        {
            message = ZBAD;
            break;
        }
        uint8_t fhd = s[4];
        int fcs_flag = fhd >> 6, single = (fhd >> 5) & 1, checksum = (fhd >> 2) & 1, did = fhd & 3;
        if (fhd & 8)
        {
            message = ZBAD;
            break;
        }
        size_t h = 5 + (single ? 0 : 1);
        static const int did_size[4] = {0, 1, 2, 4};
        static const int fcs_size[4] = {0, 2, 4, 8};
        int fcs_n = fcs_flag == 0 && single ? 1 : fcs_size[fcs_flag];
        if (did)
        {
            /* A dictionary id of zero means none */
            uint32_t id = 0;
            if (n < h + (size_t)did_size[did])
            {
                message = ZBAD;
                break;
            }
            for (int k = 0; k < did_size[did]; k++)
                id |= (uint32_t)s[h + k] << (8 * k);
            if (id)
            {
                message = "ZSTD decompression failed: dictionaries are not supported";
                break;
            }
        }
        h += (size_t)did_size[did];
        if (n < h + (size_t)fcs_n)
        {
            message = ZBAD;
            break;
        }
        uint64_t fcs = 0;
        int have_fcs = fcs_n > 0;
        for (int k = 0; k < fcs_n; k++)
            fcs |= (uint64_t)s[h + k] << (8 * k);
        if (fcs_n == 2)
            fcs += 256;
        h += (size_t)fcs_n;
        s += h;
        n -= h;
        size_t start = pos;
        memset(z, 0, offsetof(zstate, lit));
        z->rep[0] = 1;
        z->rep[1] = 4;
        z->rep[2] = 8;
        for (;;)
        {
            if (n < 3)
            {
                message = ZBAD;
                break;
            }
            uint32_t bh = s[0] | (uint32_t)s[1] << 8 | (uint32_t)s[2] << 16;
            s += 3;
            n -= 3;
            int last = bh & 1, type = (bh >> 1) & 3;
            size_t size = bh >> 3;
            if (type == 0)
            {
                if (size > n || size > out - pos)
                {
                    message = ZBAD;
                    break;
                }
                memcpy(dst + pos, s, size);
                pos += size;
                s += size;
                n -= size;
            }
            else if (type == 1)
            {
                if (n < 1 || size > out - pos)
                {
                    message = ZBAD;
                    break;
                }
                memset(dst + pos, s[0], size);
                pos += size;
                s += 1;
                n -= 1;
            }
            else if (type == 2)
            {
                if (size > n || size > (1 << 17))
                {
                    message = ZBAD;
                    break;
                }
                message = block(z, s, size, dst, &pos, out);
                if (message)
                    break;
                s += size;
                n -= size;
            }
            else
            {
                message = ZBAD;
                break;
            }
            if (last)
                break;
        }
        if (message)
            break;
        if (have_fcs && pos - start != fcs)
        {
            message = ZBAD;
            break;
        }
        if (checksum)
        {
            if (n < 4)
            {
                message = ZBAD;
                break;
            }
            s += 4;
            n -= 4;
        }
    }
    free(z);
    if (!message && pos != out)
        message = ZBAD;
    return message;
}

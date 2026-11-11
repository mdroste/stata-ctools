/* Column chunk decoding: pages, definition levels and value encodings, and
 * the conversion of each value to what import parquet stores. */
#include "parquet_internal.h"
#include "ctools_file.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
    const uint8_t *p;
    uint32_t n;
} pq_str;

/* Decoded values of one page (or of a dictionary) */
typedef struct
{
    size_t count;
    int width;        /* fixed-width physical types (FLBA: its length) */
    uint8_t *fixed;   /* count * width bytes */
    pq_str *strs;     /* BYTE_ARRAY values */
    uint8_t *scratch; /* bytes rebuilt by DELTA_BYTE_ARRAY */
} pq_vals;

static void vals_free(pq_vals *v)
{
    free(v->fixed);
    free(v->strs);
    free(v->scratch);
    memset(v, 0, sizeof(*v));
}

static uint32_t le32(const uint8_t *p)
{
    return p[0] | (uint32_t)p[1] << 8 | (uint32_t)p[2] << 16 | (uint32_t)p[3] << 24;
}

static int uleb(const uint8_t **p, const uint8_t *end, uint64_t *v)
{
    uint64_t out = 0;
    for (int shift = 0; shift < 70; shift += 7)
    {
        if (*p >= end)
            return 1;
        uint8_t b = *(*p)++;
        if (shift < 64)
            out |= (uint64_t)(b & 0x7f) << shift;
        if (!(b & 0x80))
        {
            *v = out;
            return 0;
        }
    }
    return 1;
}

/* LSB-first unpacking of count values of width bits */
static void unpack(const uint8_t *s, int width, size_t count, uint32_t *out)
{
    uint64_t mask = width >= 32 ? 0xffffffffull : ((1ull << width) - 1);
    size_t bit = 0;
    for (size_t i = 0; i < count; i++, bit += (size_t)width)
    {
        size_t byte = bit >> 3;
        int shift = (int)(bit & 7);
        uint64_t w = 0;
        int need = (shift + width + 7) >> 3;
        for (int k = 0; k < need; k++)
            w |= (uint64_t)s[byte + k] << (8 * k);
        out[i] = (uint32_t)((w >> shift) & mask);
    }
}

/* RLE / bit-packed hybrid: exactly count values (a final bit-packed group
 * may hold more, which are ignored). Returns nonzero when data run out. */
static int rle_decode(const uint8_t *p, size_t n, int width, uint32_t *out, size_t count)
{
    const uint8_t *end = p + n;
    size_t got = 0;
    if (width < 0 || width > 32)
        return 1;
    int nbytes = (width + 7) / 8;
    while (got < count)
    {
        uint64_t header;
        if (uleb(&p, end, &header))
            return 1;
        if (header & 1)
        {
            uint64_t groups = header >> 1;
            if (groups == 0 || (width && groups > (uint64_t)(end - p) / (uint64_t)width))
                return 1;
            size_t bytes = (size_t)groups * (size_t)width;
            size_t left = count - got;
            size_t take = groups >= left / 8 + 1 ? left : (size_t)groups * 8;
            if (width)
                unpack(p, width, take, out + got);
            else
                memset(out + got, 0, take * sizeof(uint32_t));
            got += take;
            p += bytes;
        }
        else
        {
            uint64_t run = header >> 1;
            if ((size_t)(end - p) < (size_t)nbytes)
                return 1;
            uint32_t v = 0;
            for (int k = 0; k < nbytes; k++)
                v |= (uint32_t)p[k] << (8 * k);
            p += nbytes;
            if (width < 32 && (v >> width))
                return 1;
            if (run == 0)
                return 1;
            size_t take = run < count - got ? (size_t)run : count - got;
            for (size_t i = 0; i < take; i++)
                out[got + i] = v;
            got += take;
        }
    }
    return 0;
}

/* DELTA_BINARY_PACKED; *used receives the bytes consumed. bits is 32 or 64:
 * arithmetic wraps at that width, as the encoding requires. */
static int delta_decode(const uint8_t *p, size_t n, int bits, int64_t *out, size_t count,
                        size_t *total, size_t *used)
{
    const uint8_t *start = p, *end = p + n;
    uint64_t block, minis, all, first;
    if (uleb(&p, end, &block) || uleb(&p, end, &minis) || uleb(&p, end, &all) ||
        uleb(&p, end, &first))
        return 1;
    if (block == 0 || block % 128 || minis == 0 || block % minis || (block / minis) % 32 ||
        block > 1u << 20)
        return 1;
    if (all > count)
        return 1;
    *total = (size_t)all;
    uint64_t mask = bits == 32 ? 0xffffffffull : ~0ull;
    uint64_t value = ((first >> 1) ^ -(first & 1)) & mask;
    size_t got = 0, per = (size_t)(block / minis);
    if (all)
        out[got++] = bits == 32 ? (int64_t)(int32_t)(uint32_t)value : (int64_t)value;
    while (got < all)
    {
        uint64_t zmin;
        if (uleb(&p, end, &zmin) || (uint64_t)(end - p) < minis)
            return 1;
        uint64_t mindelta = (zmin >> 1) ^ -(zmin & 1);
        const uint8_t *widths = p;
        p += minis;
        for (uint64_t m = 0; m < minis && got < all; m++)
        {
            int w = widths[m];
            size_t bytes = per * (size_t)w / 8;
            if (w > bits || bytes > (size_t)(end - p))
                return 1;
            size_t take = per < all - got ? per : (size_t)(all - got);
            for (size_t i = 0; i < take; i++)
            {
                uint64_t d;
                if (w <= 32)
                {
                    uint32_t one;
                    if (w)
                    {
                        size_t bit = i * (size_t)w;
                        uint64_t word = 0;
                        int need = ((int)(bit & 7) + w + 7) >> 3;
                        for (int k = 0; k < need; k++)
                            word |= (uint64_t)p[(bit >> 3) + k] << (8 * k);
                        one = (uint32_t)((word >> (bit & 7)) & (w == 32 ? 0xffffffffull : ((1ull << w) - 1)));
                    }
                    else
                        one = 0;
                    d = one;
                }
                else
                {
                    size_t bit = i * (size_t)w;
                    d = 0;
                    for (int k = 0; k < w; k++, bit++)
                        d |= (uint64_t)((p[bit >> 3] >> (bit & 7)) & 1) << k;
                }
                value = (value + mindelta + d) & mask;
                out[got++] = bits == 32 ? (int64_t)(int32_t)(uint32_t)value : (int64_t)value;
            }
            p += bytes;
        }
    }
    *used = (size_t)(p - start);
    return 0;
}

static int fixed_width(int physical, int type_length)
{
    switch (physical)
    {
    case PQ_BOOLEAN:
        return 1;
    case PQ_INT32:
    case PQ_FLOAT:
        return 4;
    case PQ_INT64:
    case PQ_DOUBLE:
        return 8;
    case PQ_INT96:
        return 12;
    case PQ_FLBA:
        return type_length;
    }
    return 0;
}

static int alloc_vals(pq_vals *v, int physical, int width, size_t count)
{
    v->count = count;
    v->width = width;
    if (physical == PQ_BYTE_ARRAY)
        v->strs = malloc((count ? count : 1) * sizeof(pq_str));
    else
    {
        if (width && count > SIZE_MAX / (size_t)width)
            return 909;
        v->fixed = malloc(count * (size_t)width + 1);
    }
    return (v->strs || v->fixed) ? 0 : 909;
}

/* PLAIN values */
static int plain(const uint8_t *p, size_t n, int physical, int width, pq_vals *v, size_t count)
{
    if (alloc_vals(v, physical, width, count))
        return 909;
    if (physical == PQ_BOOLEAN)
    {
        if ((count + 7) / 8 > n)
            return 1;
        for (size_t i = 0; i < count; i++)
            v->fixed[i] = (p[i >> 3] >> (i & 7)) & 1;
        return 0;
    }
    if (physical == PQ_BYTE_ARRAY)
    {
        const uint8_t *end = p + n;
        for (size_t i = 0; i < count; i++)
        {
            if (end - p < 4)
                return 1;
            uint32_t len = le32(p);
            p += 4;
            if (len > (size_t)(end - p))
                return 1;
            v->strs[i].p = p;
            v->strs[i].n = len;
            p += len;
        }
        return 0;
    }
    if (count * (size_t)width > n)
        return 1;
    memcpy(v->fixed, p, count * (size_t)width);
    return 0;
}

static int delta_strings(const uint8_t *p, size_t n, pq_vals *v, size_t count, int prefixed, int width)
{
    int64_t *prefix = NULL, *lens = malloc((count ? count : 1) * sizeof(int64_t));
    size_t total, used = 0, used2 = 0;
    int rc = 1;
    if (!lens)
        return 909;
    if (prefixed)
    {
        prefix = malloc((count ? count : 1) * sizeof(int64_t));
        if (!prefix)
        {
            free(lens);
            return 909;
        }
        if ((rc = delta_decode(p, n, 32, prefix, count, &total, &used)) || total != count)
        {
            rc = rc ? rc : 1;
            goto done;
        }
        rc = 1;
    }
    if ((rc = delta_decode(p + used, n - used, 32, lens, count, &total, &used2)) || total != count)
    {
        rc = rc ? rc : 1;
        goto done;
    }
    rc = 1;
    const uint8_t *data = p + used + used2, *end = p + n;
    size_t need = 0;
    for (size_t i = 0; i < count; i++)
    {
        if (lens[i] < 0 || (prefixed && prefix[i] < 0))
            goto done;
        need += (size_t)lens[i] + (prefixed ? (size_t)prefix[i] : 0);
        if (need > ((size_t)1 << 40))
            goto done;
    }
    if (width)
        v->fixed = malloc(count * (size_t)width + 1);
    else
        v->strs = malloc((count ? count : 1) * sizeof(pq_str));
    v->count = count;
    v->width = width;
    if (prefixed)
        v->scratch = malloc(need + 1);
    if ((width ? !v->fixed : !v->strs) || (prefixed && !v->scratch))
    {
        rc = 909;
        goto done;
    }
    size_t at = 0;
    const uint8_t *prev = NULL;
    size_t prevn = 0;
    for (size_t i = 0; i < count; i++)
    {
        size_t len = (size_t)lens[i];
        if (len > (size_t)(end - data))
            goto done;
        const uint8_t *value = data;
        size_t vn = len;
        if (prefixed)
        {
            size_t k = (size_t)prefix[i];
            if (k > prevn)
                goto done;
            uint8_t *dst = v->scratch + at;
            if (k)
                memmove(dst, prev, k);
            memcpy(dst + k, data, len);
            value = dst;
            vn = k + len;
            at += vn;
        }
        data += len;
        prev = value;
        prevn = vn;
        if (width)
        {
            if (vn != (size_t)width)
                goto done;
            memcpy(v->fixed + i * (size_t)width, value, vn);
        }
        else
        {
            v->strs[i].p = value;
            v->strs[i].n = (uint32_t)vn;
        }
    }
    rc = 0;
done:
    free(prefix);
    free(lens);
    return rc;
}

/* Decodes count non-null values of a data page's value section */
static int decode_values(const uint8_t *p, size_t n, int encoding, const pq_field *fd,
                         const pq_vals *dict, pq_vals *v, size_t count)
{
    int physical = fd->physical, width = fixed_width(physical, fd->type_length);
    switch (encoding)
    {
    case PQ_PLAIN:
        return plain(p, n, physical, width, v, count);
    case PQ_PLAIN_DICTIONARY:
    case PQ_RLE_DICTIONARY:
    {
        if (!dict)
            return 1;
        if (n < 1 && count)
            return 1;
        uint32_t *idx = malloc((count ? count : 1) * sizeof(uint32_t));
        if (!idx)
            return 909;
        if (count && rle_decode(p + 1, n - 1, p[0], idx, count))
        {
            free(idx);
            return 1;
        }
        if (alloc_vals(v, physical, width, count))
        {
            free(idx);
            return 909;
        }
        for (size_t i = 0; i < count; i++)
        {
            if (idx[i] >= dict->count)
            {
                free(idx);
                return 1;
            }
            if (physical == PQ_BYTE_ARRAY)
                v->strs[i] = dict->strs[idx[i]];
            else
                memcpy(v->fixed + i * (size_t)width, dict->fixed + (size_t)idx[i] * (size_t)width,
                       (size_t)width);
        }
        free(idx);
        return 0;
    }
    case PQ_RLE:
    {
        if (physical != PQ_BOOLEAN || n < 4)
            return 1;
        uint32_t len = le32(p);
        if (len > n - 4 || alloc_vals(v, physical, 1, count))
            return 1;
        uint32_t *bits = malloc((count ? count : 1) * sizeof(uint32_t));
        if (!bits)
            return 909;
        int bad = rle_decode(p + 4, len, 1, bits, count);
        for (size_t i = 0; !bad && i < count; i++)
            v->fixed[i] = (uint8_t)bits[i];
        free(bits);
        return bad;
    }
    case PQ_DELTA_BINARY_PACKED:
    {
        if (physical != PQ_INT32 && physical != PQ_INT64)
            return 1;
        int64_t *tmp = malloc((count ? count : 1) * sizeof(int64_t));
        size_t total, used;
        if (!tmp)
            return 909;
        int rc = delta_decode(p, n, physical == PQ_INT32 ? 32 : 64, tmp, count, &total, &used);
        if (!rc && total != count)
            rc = 1;
        if (!rc && alloc_vals(v, physical, width, count))
            rc = 909;
        for (size_t i = 0; !rc && i < count; i++)
        {
            if (physical == PQ_INT32)
            {
                int32_t x = (int32_t)tmp[i];
                memcpy(v->fixed + 4 * i, &x, 4);
            }
            else
                memcpy(v->fixed + 8 * i, &tmp[i], 8);
        }
        free(tmp);
        return rc;
    }
    case PQ_DELTA_LENGTH_BYTE_ARRAY:
        if (physical != PQ_BYTE_ARRAY)
            return 1;
        return delta_strings(p, n, v, count, 0, 0);
    case PQ_DELTA_BYTE_ARRAY:
        if (physical != PQ_BYTE_ARRAY && physical != PQ_FLBA)
            return 1;
        return delta_strings(p, n, v, count, 1, physical == PQ_FLBA ? width : 0);
    case PQ_BYTE_STREAM_SPLIT:
    {
        if (physical != PQ_FLOAT && physical != PQ_DOUBLE && physical != PQ_INT32 &&
            physical != PQ_INT64 && physical != PQ_FLBA)
            return 1;
        if (count * (size_t)width > n || alloc_vals(v, physical, width, count))
            return 1;
        for (int k = 0; k < width; k++)
            for (size_t i = 0; i < count; i++)
                v->fixed[i * (size_t)width + (size_t)k] = p[(size_t)k * count + i];
        return 0;
    }
    }
    return 2; /* unsupported encoding */
}

/* ---- Conversion to Stata values ------------------------------------------- */
static double half_to_double(uint16_t h)
{
    int e = (h >> 10) & 31, m = h & 1023;
    double v;
    if (e == 0)
        v = ldexp((double)m, -24);
    else if (e == 31)
        v = m ? NAN : INFINITY;
    else
        v = ldexp((double)(m | 1024), e - 25);
    return h & 0x8000 ? -v : v;
}

/* Arrow's Decimal ToString: plain digits with the point inserted, or
 * scientific notation when the scale is negative or the adjusted exponent
 * is below -6. be holds a big-endian two's complement integer. */
void pq_decimal_text(const uint8_t *be, size_t n, int scale, char *out, size_t cap)
{
    uint32_t limb[9] = {0}; /* little-endian 32-bit limbs, up to 288 bits */
    int neg = n && (be[0] & 0x80);
    size_t nl = 0;
    if (n > 32)
    {
        snprintf(out, cap, "%s", "");
        return;
    }
    for (size_t i = 0; i < n; i++)
    {
        size_t bit = (n - 1 - i) * 8;
        limb[bit / 32] |= (uint32_t)be[i] << (bit % 32);
    }
    nl = (n * 8 + 31) / 32;
    if (neg)
    {
        /* sign-extend to whole limbs, then negate */
        size_t bits = n * 8;
        if (bits % 32)
            limb[nl - 1] |= ~0u << (bits % 32);
        uint64_t carry = 1;
        for (size_t i = 0; i < nl; i++)
        {
            uint64_t x = (uint64_t)(uint32_t)~limb[i] + carry;
            limb[i] = (uint32_t)x;
            carry = x >> 32;
        }
    }
    char digits[100];
    int nd = 0;
    for (;;)
    {
        int nonzero = 0;
        for (size_t i = 0; i < nl; i++)
            if (limb[i])
                nonzero = 1;
        if (!nonzero)
            break;
        uint64_t rem = 0;
        for (size_t i = nl; i-- > 0;)
        {
            uint64_t cur = (rem << 32) | limb[i];
            limb[i] = (uint32_t)(cur / 1000000000u);
            rem = cur % 1000000000u;
        }
        int more = 0;
        for (size_t i = 0; i < nl; i++)
            if (limb[i])
                more = 1;
        for (int k = 0; k < 9 && (more || rem); k++)
        {
            digits[nd++] = (char)('0' + rem % 10);
            rem /= 10;
        }
        if (nd > 90)
            break;
    }
    if (nd == 0)
        digits[nd++] = '0';
    /* digits are least significant first */
    char d[100];
    for (int i = 0; i < nd; i++)
        d[i] = digits[nd - 1 - i];
    d[nd] = 0;
    char buf[200];
    int pos = 0;
    if (neg)
        buf[pos++] = '-';
    int adjusted = nd - 1 - scale;
    if (scale < 0 || adjusted < -6)
    {
        buf[pos++] = d[0];
        if (nd > 1)
        {
            buf[pos++] = '.';
            memcpy(buf + pos, d + 1, (size_t)nd - 1);
            pos += nd - 1;
        }
        pos += snprintf(buf + pos, sizeof(buf) - (size_t)pos, "E%s%d", adjusted >= 0 ? "+" : "",
                        adjusted);
    }
    else if (scale == 0)
    {
        memcpy(buf + pos, d, (size_t)nd);
        pos += nd;
    }
    else if (nd > scale)
    {
        memcpy(buf + pos, d, (size_t)(nd - scale));
        pos += nd - scale;
        buf[pos++] = '.';
        memcpy(buf + pos, d + nd - scale, (size_t)scale);
        pos += scale;
    }
    else
    {
        buf[pos++] = '0';
        buf[pos++] = '.';
        for (int i = 0; i < scale - nd; i++)
            buf[pos++] = '0';
        memcpy(buf + pos, d, (size_t)nd);
        pos += nd;
    }
    buf[pos] = 0;
    snprintf(out, cap, "%s", buf);
}

static int is_text(pq_arrow a)
{
    return a == PA_UTF8 || a == PA_BINARY || a == PA_FLBA || a == PA_DECIMAL128 ||
           a == PA_DECIMAL256 || a == PA_LARGE_UTF8 || a == PA_LARGE_BINARY;
}

/* Columns whose result depends on every row: text widths and dictionaries */
static int whole_column(pq_arrow a)
{
    return is_text(a) || a == PA_DICT;
}

static uint64_t hash_bytes(const uint8_t *p, size_t n)
{
    uint64_t h = 1469598103934665603ull;
    for (size_t i = 0; i < n; i++)
        h = (h ^ p[i]) * 1099511628211ull;
    return h;
}

/* The code of a dictionary entry, appended when new (entries keep the order
   in which Arrow unifies them: each dictionary page in turn, then values of
   pages that are not dictionary encoded, in order of appearance) */
static int dict_code(pq_values *v, const uint8_t *p, size_t n, uint32_t *code)
{
    if ((v->dn + 1) * 2 > v->dhcap)
    {
        size_t cap = v->dhcap ? v->dhcap * 2 : 1024;
        uint32_t *h = calloc(cap, sizeof(uint32_t));
        if (!h)
            return 909;
        for (size_t e = 0; e < v->dn; e++)
        {
            size_t at = (size_t)hash_bytes(v->dbytes + v->doff[e], (size_t)(v->doff[e + 1] - v->doff[e])) & (cap - 1);
            while (h[at])
                at = (at + 1) & (cap - 1);
            h[at] = (uint32_t)e + 1;
        }
        free(v->dhash);
        v->dhash = h;
        v->dhcap = cap;
    }
    size_t at = (size_t)hash_bytes(p, n) & (v->dhcap - 1);
    while (v->dhash[at])
    {
        size_t e = v->dhash[at] - 1, len = (size_t)(v->doff[e + 1] - v->doff[e]);
        if (len == n && (!n || !memcmp(v->dbytes + v->doff[e], p, n)))
        {
            *code = (uint32_t)e;
            return 0;
        }
        at = (at + 1) & (v->dhcap - 1);
    }
    if (v->dn >= UINT32_MAX - 1)
        return 1;
    if (v->dn + 2 > v->dcap)
    {
        size_t cap = v->dcap ? v->dcap * 2 : 256;
        uint64_t *o = realloc(v->doff, cap * sizeof(uint64_t));
        if (!o)
            return 909;
        v->doff = o;
        v->dcap = cap;
        if (!v->dn)
            v->doff[0] = 0;
    }
    if (v->dbcap - v->dbn < n + 1)
    {
        size_t cap = v->dbcap ? v->dbcap : 4096;
        while (cap - v->dbn < n + 1)
            cap *= 2;
        uint8_t *b = realloc(v->dbytes, cap);
        if (!b)
            return 909;
        v->dbytes = b;
        v->dbcap = cap;
    }
    memcpy(v->dbytes + v->dbn, p, n);
    v->dbn += n;
    v->doff[v->dn + 1] = v->dbn;
    v->dhash[at] = (uint32_t)v->dn + 1;
    *code = (uint32_t)v->dn++;
    if (n > v->dmax_len)
        v->dmax_len = n;
    return 0;
}

static int put_bytes(pq_values *out, size_t row, const uint8_t *p, size_t n)
{
    if (out->capbytes - out->nbytes < n)
    {
        size_t cap = out->capbytes ? out->capbytes : 4096;
        while (cap - out->nbytes < n)
        {
            if (cap > SIZE_MAX / 2)
                return 909;
            cap *= 2;
        }
        uint8_t *b = realloc(out->bytes, cap);
        if (!b)
            return 909;
        out->bytes = b;
        out->capbytes = cap;
    }
    if (n)
        memcpy(out->bytes + out->nbytes, p, n);
    out->nbytes += n;
    out->off[row + 1] = out->nbytes;
    if (n && memchr(p, 0, n))
        out->has_nul = 1;
    return 0;
}

/* Stores the value of one row; idx is the row's position in out (or
 * SIZE_MAX when the row is only measured). */
static int emit(pq_values *out, const pq_field *fd, size_t idx, const pq_vals *v, size_t k)
{
    /* No successful decode may expose an absent or undersized value buffer. */
    if (k >= v->count || (idx != SIZE_MAX && idx >= out->n))
        return 1;
    if (fd->physical == PQ_BYTE_ARRAY) {
        if (!v->strs || (!v->strs[k].p && v->strs[k].n)) return 1;
    } else if (!v->fixed || v->width <= 0) {
        return 1;
    }
    pq_arrow a = fd->arrow;
    if (a == PA_DICT)
    {
        uint32_t code;
        int rc = dict_code(out, v->strs[k].p, v->strs[k].n, &code);
        if (rc || idx == SIZE_MAX)
            return rc;
        out->num[idx] = code;
        out->valid[idx] = 1;
        return 0;
    }
    if (is_text(a))
    {
        uint8_t tmp[64];
        const uint8_t *p;
        size_t n;
        char text[120];
        if (a == PA_DECIMAL128 || a == PA_DECIMAL256)
        {
            if (fd->physical == PQ_INT32 || fd->physical == PQ_INT64)
            {
                int w = fd->physical == PQ_INT32 ? 4 : 8;
                for (int b = 0; b < w; b++)
                    tmp[b] = v->fixed[k * (size_t)w + (size_t)(w - 1 - b)];
                p = tmp;
                n = (size_t)w;
            }
            else if (fd->physical == PQ_BYTE_ARRAY)
            {
                p = v->strs[k].p;
                n = v->strs[k].n;
            }
            else
            {
                p = v->fixed + k * (size_t)v->width;
                n = (size_t)v->width;
            }
            if (n > 32)
                return 1;
            if (idx == SIZE_MAX)
                return 0;
            pq_decimal_text(p, n, fd->scale, text, sizeof(text));
            p = (const uint8_t *)text;
            n = strlen(text);
        }
        else if (fd->physical == PQ_BYTE_ARRAY)
        {
            p = v->strs[k].p;
            n = v->strs[k].n;
        }
        else
        {
            p = v->fixed + k * (size_t)v->width;
            n = (size_t)v->width;
        }
        if (n > out->max_len)
            out->max_len = n;
        if (idx == SIZE_MAX)
            return 0;
        out->valid[idx] = 1;
        return put_bytes(out, idx, p, n);
    }
    if (idx == SIZE_MAX)
        return 0;
    const uint8_t *f = v->fixed + k * (size_t)v->width;
    int32_t i32 = 0;
    int64_t i64 = 0;
    if (v->width == 4)
        memcpy(&i32, f, 4);
    if (v->width == 8)
        memcpy(&i64, f, 8);
    double x = 0;
    switch (a)
    {
    case PA_BOOL:
        x = f[0] ? 1 : 0;
        break;
    case PA_INT8:
        x = (int8_t)i32;
        break;
    case PA_INT16:
        x = (int16_t)i32;
        break;
    case PA_INT32:
        x = i32;
        break;
    case PA_UINT8:
        x = (uint8_t)i32;
        break;
    case PA_UINT16:
        x = (uint16_t)i32;
        break;
    case PA_UINT32:
        x = (uint32_t)i32;
        break;
    case PA_INT64:
    case PA_UINT64:
        out->i64[idx] = i64;
        out->valid[idx] = 1;
        if (a == PA_INT64 ? (i64 < -2147483647LL || i64 > 2147483620LL)
                          : (uint64_t)i64 > 2147483620ULL)
            out->wide = 1;
        out->any = 1;
        return 0;
    case PA_HALF:
    {
        uint16_t h = (uint16_t)(f[0] | f[1] << 8);
        x = half_to_double(h);
        break;
    }
    case PA_FLOAT:
    {
        float fl;
        memcpy(&fl, f, 4);
        x = fl;
        break;
    }
    case PA_DOUBLE:
        memcpy(&x, f, 8);
        break;
    case PA_DATE32:
        x = (double)i32 + 3653.0;
        break;
    case PA_TIME32:
        x = i32;
        break;
    case PA_TIME64:
        x = (double)(i64 / (fd->unit == 2 ? 1000 : 1000000));
        break;
    case PA_TIMESTAMP:
        if (fd->physical == PQ_INT96)
        {
            uint64_t nanos;
            uint32_t julian;
            memcpy(&nanos, f, 8);
            memcpy(&julian, f + 8, 4);
            /* Arrow's conversion: unsigned arithmetic that wraps, and a
               Julian day of zero reads as the epoch */
            int64_t ns = julian ? (int64_t)(((uint64_t)julian - 2440588u) * 86400000000000ull + nanos)
                                : 0;
            x = (double)(ns / 1000000) + 315619200000.0;
        }
        else
        {
            int64_t ms = fd->unit == 1 ? i64 : fd->unit == 2 ? i64 / 1000 : i64 / 1000000;
            x = (double)ms + 315619200000.0;
        }
        break;
    default:
        return 1;
    }
    out->num[idx] = x;
    out->valid[idx] = 1;
    if (!isnan(x) && (!out->any || x > out->maxv))
        out->maxv = x;
    if (!isnan(x))
        out->any = 1;
    return 0;
}

/* ---- Pages ------------------------------------------------------------------ */
typedef struct
{
    pq_file *f;
    const pq_field *fd;
    pq_values *out;
    int64_t row;   /* first row of the group */
    int64_t first, last;
    int measure;   /* text lengths over every row */
    char *err;
    size_t errcap;
    FILE *fp; /* this reader's own handle */
} pq_ctx;

static int fail(pq_ctx *c, const char *message)
{
    snprintf(c->err, c->errcap, "%s", message);
    return 691;
}

static int read_chunk(pq_file *f, FILE *fp, const pq_chunk *ch, uint8_t **buf, size_t *len)
{
    int64_t start = ch->data_offset;
    if (ch->dict_offset > 0 && ch->dict_offset < start)
        start = ch->dict_offset;
    int64_t n = ch->total_compressed;
    if (start < 4 || n <= 0 || start > f->size || n > f->size - start)
        return 1;
    *buf = malloc((size_t)n);
    if (!*buf)
        return 909;
    if (ctools_fseek(fp, start, SEEK_SET) || fread(*buf, 1, (size_t)n, fp) != (size_t)n)
        return 1;
    *len = (size_t)n;
    return 0;
}

/* Decodes (check_only: decompresses) one column chunk of a row group */
static int read_group(pq_ctx *c, const pq_chunk *ch, int64_t nrows, int leaf_only)
{
    const pq_field *fd = c->fd;
    uint8_t *buf = NULL, *page = NULL;
    size_t len = 0, pagecap = 0;
    pq_vals dict = {0}, vals = {0};
    uint32_t *defs = NULL;
    size_t defcap = 0;
    int have_dict = 0, rc = 0;
    if (!ch->has_meta || ch->has_path)
        return fail(c, "column chunk metadata is missing or refers to another file");
    if (ch->num_values == 0 && nrows == 0)
        return 0;
    if ((rc = read_chunk(c->f, c->fp, ch, &buf, &len)))
    {
        free(buf);
        return rc == 909 ? 909 : fail(c, "column chunk lies outside the file");
    }
    size_t pos = 0;
    int64_t done = 0, values_seen = 0;
    int64_t want = leaf_only ? ch->num_values : nrows;
    while ((leaf_only ? values_seen : done) < want)
    {
        pq_page h;
        size_t used;
        if (pos >= len || pq_page_header(buf + pos, len - pos, &h, &used))
        {
            rc = fail(c, "Couldn't deserialize thrift: page header is corrupt");
            break;
        }
        pos += used;
        if ((size_t)h.csize > len - pos)
        {
            rc = fail(c, "page lies outside its column chunk");
            break;
        }
        const uint8_t *body = buf + pos;
        pos += (size_t)h.csize;
        if (h.type == 1)
            continue;
        if (h.type != 0 && h.type != 2 && h.type != 3)
            continue;
        /* Levels of a v2 page are stored uncompressed before the values */
        size_t lev = h.type == 3 ? (size_t)h.def_len + (size_t)h.rep_len : 0;
        if (lev > (size_t)h.csize || lev > (size_t)h.usize)
        {
            rc = fail(c, "page level sizes are inconsistent");
            break;
        }
        size_t usize = (size_t)h.usize;
        /* Reject an impossible uncompressed page before trusting its declared
         * allocation size. V2 levels are included in both size fields. */
        int codec = h.type == 3 && !h.compressed ? PQ_UNCOMPRESSED : ch->codec;
        if (codec == PQ_UNCOMPRESSED && h.csize != h.usize) {
            rc = fail(c, "uncompressed page sizes are inconsistent");
            break;
        }
        if (usize + 1 > pagecap)
        {
            free(page);
            pagecap = usize + 1;
            page = malloc(pagecap);
            if (!page)
            {
                rc = 909;
                break;
            }
        }
        memcpy(page, body, lev);
        /* An uncompressed v2 page is the UNCOMPRESSED codec's copy */
        const char *msg =
            pq_decompress(codec, body + lev, (size_t)h.csize - lev, page + lev, usize - lev);
        if (msg)
        {
            rc = fail(c, msg);
            break;
        }
        if (h.type == 2)
        {
            if (leaf_only)
                continue;
            vals_free(&dict);
            int width = fixed_width(fd->physical, fd->type_length);
            if ((h.encoding != PQ_PLAIN && h.encoding != PQ_PLAIN_DICTIONARY) ||
                (rc = plain(page, usize, fd->physical, width, &dict, (size_t)h.num_values)))
            {
                rc = rc == 909 ? 909 : fail(c, "dictionary page is corrupt");
                break;
            }
            /* String dictionary values point into this page: keep it */
            if (fd->physical == PQ_BYTE_ARRAY)
            {
                dict.scratch = page;
                page = NULL;
                pagecap = 0;
            }
            have_dict = 1;
            /* A dictionary column takes every entry, used or not, in order */
            for (size_t e = 0; fd->arrow == PA_DICT && e < dict.count && !rc; e++)
            {
                uint32_t code;
                rc = dict_code(c->out, dict.strs[e].p, dict.strs[e].n, &code);
                if (rc && rc != 909)
                    rc = fail(c, "too many dictionary entries");
            }
            if (rc)
                break;
            continue;
        }
        values_seen += h.num_values;
        if (leaf_only)
            continue;
        size_t nv = (size_t)h.num_values;
        if (nv > (size_t)(nrows - done))
        {
            rc = fail(c, "a page holds more values than its row group has rows");
            break;
        }
        /* Definition levels (flat columns: 0 null, 1 present) */
        const uint8_t *vp = page;
        size_t vn = usize;
        if (nv + 1 > defcap)
        {
            free(defs);
            defcap = nv + 1;
            defs = malloc(defcap * sizeof(uint32_t));
            if (!defs)
            {
                rc = 909;
                break;
            }
        }
        size_t present = nv;
        if (fd->max_def)
        {
            int bad = 0;
            if (h.type == 3)
            {
                bad = rle_decode(page + h.rep_len, (size_t)h.def_len, 1, defs, nv);
                vp = page + lev;
                vn = usize - lev;
            }
            else if (h.def_enc == PQ_BIT_PACKED)
            {
                /* The format specifies most significant bit first; Arrow,
                   and so import parquet, reads least significant first. */
                size_t bytes = (nv + 7) / 8;
                if (bytes > usize)
                    bad = 1;
                else
                    for (size_t i = 0; i < nv; i++)
                        defs[i] = (page[i >> 3] >> (i & 7)) & 1;
                vp = page + bytes;
                vn = usize - bytes;
            }
            else
            {
                if (usize < 4 || le32(page) > usize - 4)
                    bad = 1;
                else
                {
                    uint32_t dl = le32(page);
                    bad = rle_decode(page + 4, dl, 1, defs, nv);
                    vp = page + 4 + dl;
                    vn = usize - 4 - dl;
                }
            }
            if (bad)
            {
                rc = fail(c, "definition levels are corrupt");
                break;
            }
            present = 0;
            for (size_t i = 0; i < nv; i++)
                present += defs[i];
        }
        else if (h.type == 3)
        {
            vp = page + lev;
            vn = usize - lev;
        }
        /* Rows of this page that are needed at all */
        int64_t p0 = c->row + done, p1 = p0 + (int64_t)nv - 1;
        int in_range = !(p1 < c->first || p0 > c->last);
        if (!in_range && !(c->measure == 2 || (c->measure && whole_column(fd->arrow))))
        {
            done += (int64_t)nv;
            continue;
        }
        vals_free(&vals);
        int d = decode_values(vp, vn, h.encoding, fd, have_dict ? &dict : NULL, &vals, present);
        if (d)
        {
            rc = d == 909 ? 909
                 : d == 2 ? fail(c, "the page uses an encoding cimport parquet does not read")
                          : fail(c, "page values are corrupt");
            break;
        }
        size_t k = 0;
        for (size_t i = 0; i < nv && !rc; i++)
        {
            int64_t r = p0 + (int64_t)i;
            size_t idx = (r >= c->first && r <= c->last) ? (size_t)(r - c->first) : SIZE_MAX;
            if (fd->max_def && !defs[i])
            {
                if (idx != SIZE_MAX && c->out->off)
                    c->out->off[idx + 1] = c->out->nbytes;
                continue;
            }
            int e = emit(c->out, fd, idx, &vals, k++);
            if (e)
                rc = e == 909 ? 909 : fail(c, "a value cannot be converted");
        }
        done += (int64_t)nv;
        if (rc)
            break;
        /* Past the requested rows: only text needs the rest, for lengths */
        if (c->row + done > c->last && !(c->measure == 2 || (c->measure && whole_column(fd->arrow))))
            break;
    }
    if (!rc && leaf_only && values_seen < want)
        rc = fail(c, "column chunk ends early");
    vals_free(&vals);
    vals_free(&dict);
    free(defs);
    free(page);
    free(buf);
    return rc;
}

void pq_values_free(pq_values *v)
{
    free(v->valid);
    free(v->num);
    free(v->i64);
    free(v->off);
    free(v->bytes);
    free(v->doff);
    free(v->dbytes);
    free(v->dhash);
    memset(v, 0, sizeof(*v));
}

int pq_read_column(pq_file *f, const pq_field *field, int64_t first, int64_t last, int measure_all,
                   pq_values *out, char *err, size_t errcap)
{
    memset(out, 0, sizeof(*out));
    out->field = field;
    size_t n = last >= first ? (size_t)(last - first + 1) : 0;
    out->n = n;
    out->valid = calloc(n ? n : 1, 1);
    if (!out->valid)
        return 909;
    if (is_text(field->arrow))
    {
        out->off = calloc(n + 1, sizeof(uint64_t));
        if (!out->off)
            return 909;
    }
    else if (field->arrow == PA_INT64 || field->arrow == PA_UINT64)
    {
        out->i64 = calloc(n ? n : 1, sizeof(int64_t));
        if (!out->i64)
            return 909;
    }
    else
    {
        out->num = malloc((n ? n : 1) * sizeof(double));
        if (!out->num)
            return 909;
    }
    pq_ctx c = {f, field, out, 0, first, last, measure_all, err, errcap, NULL};
    c.fp = ctools_fopen(f->path, "rb");
    if (!c.fp)
    {
        snprintf(err, errcap, "the file could not be read");
        return 691;
    }
    int64_t row = 0;
    for (int g = 0; g < f->ngroups; g++)
    {
        int64_t nr = f->groups[g].num_rows;
        c.row = row;
        int needed = !(row + nr - 1 < first || row > last);
        if (nr > 0 && (needed || measure_all == 2 || (measure_all && whole_column(field->arrow))))
        {
            int rc = read_group(&c, &f->groups[g].chunks[field->leaf], nr, 0);
            if (rc)
            {
                fclose(c.fp);
                return rc;
            }
        }
        row += nr;
    }
    fclose(c.fp);
    /* Rows never reached keep their offsets continuous */
    if (out->off)
        for (size_t i = 0; i < n; i++)
            if (out->off[i + 1] < out->off[i])
                out->off[i + 1] = out->off[i];
    return 0;
}

/* Decompresses every page of every leaf under a (nested) field */
int pq_check_field(pq_file *f, const pq_field *field, char *err, size_t errcap)
{
    pq_values dummy;
    memset(&dummy, 0, sizeof(dummy));
    pq_field leaf = *field;
    pq_ctx c = {f, &leaf, &dummy, 0, 0, -1, 0, err, errcap, NULL};
    c.fp = ctools_fopen(f->path, "rb");
    if (!c.fp)
    {
        snprintf(err, errcap, "the file could not be read");
        return 691;
    }
    int rc = 0;
    for (int g = 0; g < f->ngroups && !rc; g++)
        for (int l = field->first_leaf; l < field->first_leaf + field->nleaf && !rc; l++)
            rc = read_group(&c, &f->groups[g].chunks[l], f->groups[g].num_rows, 1);
    fclose(c.fp);
    return rc;
}

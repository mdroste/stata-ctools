/* Page decompression: Snappy and LZ4 decoders written here from their format
 * descriptions, Gzip/zlib through the bundled libdeflate, Zstd in
 * parquet_zstd.c, Brotli in parquet_brotli.c. Every decoder must produce
 * exactly the page's declared uncompressed size; anything else is reported as
 * corrupt data. */
#include "parquet_internal.h"
#include "vendor/libdeflate/libdeflate.h"
#include <string.h>

static const char *snappy(const uint8_t *s, size_t n, uint8_t *d, size_t out)
{
    const uint8_t *end = s + n;
    uint64_t len = 0;
    for (int shift = 0;; shift += 7)
    {
        if (s >= end || shift > 35)
            return "Corrupt snappy compressed data.";
        uint8_t b = *s++;
        len |= (uint64_t)(b & 0x7f) << shift;
        if (!(b & 0x80))
            break;
    }
    if (len != out)
        return "Corrupt snappy compressed data.";
    size_t pos = 0;
    while (s < end)
    {
        uint8_t tag = *s++;
        size_t length, offset;
        if ((tag & 3) == 0)
        {
            length = (tag >> 2) + 1;
            if (length > 60)
            {
                int bytes = (int)length - 60;
                if (end - s < bytes)
                    return "Corrupt snappy compressed data.";
                length = 0;
                for (int k = 0; k < bytes; k++)
                    length |= (size_t)s[k] << (8 * k);
                length += 1;
                s += bytes;
            }
            if ((size_t)(end - s) < length || out - pos < length)
                return "Corrupt snappy compressed data.";
            memcpy(d + pos, s, length);
            s += length;
            pos += length;
            continue;
        }
        if ((tag & 3) == 1)
        {
            if (s >= end)
                return "Corrupt snappy compressed data.";
            length = ((tag >> 2) & 7) + 4;
            offset = ((size_t)(tag >> 5) << 8) | *s++;
        }
        else if ((tag & 3) == 2)
        {
            if (end - s < 2)
                return "Corrupt snappy compressed data.";
            length = (tag >> 2) + 1;
            offset = s[0] | (size_t)s[1] << 8;
            s += 2;
        }
        else
        {
            if (end - s < 4)
                return "Corrupt snappy compressed data.";
            length = (tag >> 2) + 1;
            offset = s[0] | (size_t)s[1] << 8 | (size_t)s[2] << 16 | (size_t)s[3] << 24;
            s += 4;
        }
        if (offset == 0 || offset > pos || out - pos < length)
            return "Corrupt snappy compressed data.";
        /* Overlapping copies repeat the last offset bytes */
        for (size_t k = 0; k < length; k++, pos++)
            d[pos] = d[pos - offset];
    }
    return pos == out ? NULL : "Corrupt snappy compressed data.";
}

/* One LZ4 block; *produced receives its decompressed length */
static int lz4_block(const uint8_t *s, size_t n, uint8_t *d, size_t cap, size_t *produced)
{
    const uint8_t *end = s + n;
    size_t pos = 0;
    while (s < end)
    {
        uint8_t token = *s++;
        size_t lit = token >> 4;
        if (lit == 15)
        {
            uint8_t b;
            do
            {
                if (s >= end)
                    return 1;
                b = *s++;
                lit += b;
            } while (b == 255);
        }
        if ((size_t)(end - s) < lit || cap - pos < lit)
            return 1;
        memcpy(d + pos, s, lit);
        s += lit;
        pos += lit;
        if (s == end)
            break; /* the last sequence has literals only */
        if (end - s < 2)
            return 1;
        size_t offset = s[0] | (size_t)s[1] << 8;
        s += 2;
        size_t match = (token & 15) + 4;
        if ((token & 15) == 15)
        {
            uint8_t b;
            do
            {
                if (s >= end)
                    return 1;
                b = *s++;
                match += b;
            } while (b == 255);
        }
        if (offset == 0 || offset > pos || cap - pos < match)
            return 1;
        for (size_t k = 0; k < match; k++, pos++)
            d[pos] = d[pos - offset];
    }
    *produced = pos;
    return 0;
}

/* Parquet's LZ4 codec is the Hadoop framing (big-endian decompressed and
 * compressed sizes before each block); older writers stored a bare block,
 * which Arrow accepts when the framing does not fit, and so does this. */
static const char *lz4_hadoop(const uint8_t *s, size_t n, uint8_t *d, size_t out)
{
    const uint8_t *p = s, *end = s + n;
    size_t pos = 0;
    int ok = 1;
    while (end - p >= 8)
    {
        uint32_t raw = (uint32_t)p[0] << 24 | (uint32_t)p[1] << 16 | (uint32_t)p[2] << 8 | p[3];
        uint32_t packed = (uint32_t)p[4] << 24 | (uint32_t)p[5] << 16 | (uint32_t)p[6] << 8 | p[7];
        p += 8;
        size_t produced;
        if (packed > (size_t)(end - p) || raw > out - pos ||
            lz4_block(p, packed, d + pos, raw, &produced) || produced != raw)
        {
            ok = 0;
            break;
        }
        p += packed;
        pos += raw;
    }
    if (ok && p == end && pos == out)
        return NULL;
    size_t produced;
    if (!lz4_block(s, n, d, out, &produced) && produced == out)
        return NULL;
    return "Corrupt Lz4 compressed data.";
}

static const char *lz4_raw(const uint8_t *s, size_t n, uint8_t *d, size_t out)
{
    size_t produced;
    if (!lz4_block(s, n, d, out, &produced) && produced == out)
        return NULL;
    return "Corrupt Lz4 compressed data.";
}

/* Gzip members (RFC 1952), one after another, or a zlib stream (RFC 1950) */
static const char *gzip(const uint8_t *s, size_t n, uint8_t *d, size_t out)
{
    static const char *bad = "GZipCodec failed: corrupt gzip compressed data.";
    struct libdeflate_decompressor *z = libdeflate_alloc_decompressor();
    if (!z)
        return "GZipCodec failed: out of memory.";
    size_t pos = 0;
    const char *message = NULL;
    if (n >= 2 && s[0] != 0x1f && (s[0] & 0x0f) == 8 && ((s[0] << 8) | s[1]) % 31 == 0 &&
        !(s[1] & 0x20))
    {
        size_t in = 0, got = 0;
        if (libdeflate_deflate_decompress_ex(z, s + 2, n - 2, d, out, &in, &got) !=
                LIBDEFLATE_SUCCESS ||
            got != out || n - 2 - in < 4)
            message = bad;
        else
        {
            /* Adler-32 of the output */
            uint32_t a = 1, b = 0;
            for (size_t i = 0; i < out; i++)
            {
                a = (a + d[i]) % 65521;
                b = (b + a) % 65521;
            }
            const uint8_t *t = s + 2 + in;
            uint32_t want = (uint32_t)t[0] << 24 | (uint32_t)t[1] << 16 | (uint32_t)t[2] << 8 | t[3];
            if (want != (b << 16 | a))
                message = bad;
        }
        libdeflate_free_decompressor(z);
        return message;
    }
    while (n > 0)
    {
        if (n < 18 || s[0] != 0x1f || s[1] != 0x8b || s[2] != 8)
        {
            message = bad;
            break;
        }
        uint8_t flags = s[3];
        size_t h = 10;
        if (flags & 4)
        {
            if (n < h + 2)
            {
                message = bad;
                break;
            }
            h += 2 + (s[h] | (size_t)s[h + 1] << 8);
        }
        for (int field = 8; field <= 16; field <<= 1)
            if (flags & field)
            {
                while (h < n && s[h])
                    h++;
                h++;
            }
        if (flags & 2)
            h += 2;
        if (h + 8 > n)
        {
            message = bad;
            break;
        }
        size_t in = 0, got = 0;
        if (libdeflate_deflate_decompress_ex(z, s + h, n - h - 8, d + pos, out - pos, &in, &got) !=
            LIBDEFLATE_SUCCESS)
        {
            message = bad;
            break;
        }
        if (h + in + 8 > n)
        {
            message = bad;
            break;
        }
        const uint8_t *t = s + h + in;
        uint32_t crc = t[0] | (uint32_t)t[1] << 8 | (uint32_t)t[2] << 16 | (uint32_t)t[3] << 24;
        uint32_t size = t[4] | (uint32_t)t[5] << 8 | (uint32_t)t[6] << 16 | (uint32_t)t[7] << 24;
        if (crc != libdeflate_crc32(0, d + pos, got) || size != (uint32_t)got)
        {
            message = bad;
            break;
        }
        pos += got;
        s += h + in + 8;
        n -= h + in + 8;
        /* Trailing zero padding after the last member */
        while (n && !*s)
        {
            s++;
            n--;
        }
    }
    libdeflate_free_decompressor(z);
    if (!message && pos != out)
        message = bad;
    return message;
}

const char *pq_decompress(int codec, const uint8_t *src, size_t n, uint8_t *dst, size_t out)
{
    switch (codec)
    {
    case PQ_UNCOMPRESSED:
        if (n != out)
            return "page size does not match its declared size";
        memcpy(dst, src, n);
        return NULL;
    case PQ_SNAPPY:
        return snappy(src, n, dst, out);
    case PQ_GZIP:
        return gzip(src, n, dst, out);
    case PQ_ZSTD:
        return pq_zstd_decompress(src, n, dst, out);
    case PQ_LZ4:
        return lz4_hadoop(src, n, dst, out);
    case PQ_LZ4_RAW:
        return lz4_raw(src, n, dst, out);
    case PQ_LZO:
        return "While LZO compression is supported by the Parquet format in general, it is "
               "currently not supported by the C++ implementation.";
    case PQ_BROTLI:
        return pq_brotli_decompress(src, n, dst, out);
    }
    return "File has unsupported compression scheme";
}

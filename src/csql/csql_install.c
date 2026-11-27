/*
 * csql_install.c - unpack the DuckDB library from an official release
 * archive (duckdb-shared-libs-<platform>.tar.gz) that `csql, install` has
 * downloaded with Stata's copy.
 *
 * The archive must match the SHA-256 checksum the ado pins for the release;
 * the library is written next to its final name and moved into place only
 * when it is complete, so an interrupted install leaves no partial file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "csql/csql_internal.h"
#include "ctools_file.h"
#include "vendor/libdeflate/libdeflate.h"

/* ------------------------------------------------------------ SHA-256 */

typedef struct
{
    uint32_t h[8];
    uint64_t bytes;
    unsigned char block[64];
    size_t used;
} sha256_ctx;

static const uint32_t sha256_k[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

#define ROR32(x, n) (((x) >> (n)) | ((x) << (32 - (n))))

static void sha256_block(sha256_ctx *s, const unsigned char *p)
{
    uint32_t w[64], a, b, c, d, e, f, g, h;
    for (int i = 0; i < 16; i++)
        w[i] = (uint32_t)p[4 * i] << 24 | (uint32_t)p[4 * i + 1] << 16 | (uint32_t)p[4 * i + 2] << 8 |
               p[4 * i + 3];
    for (int i = 16; i < 64; i++) {
        uint32_t s0 = ROR32(w[i - 15], 7) ^ ROR32(w[i - 15], 18) ^ (w[i - 15] >> 3);
        uint32_t s1 = ROR32(w[i - 2], 17) ^ ROR32(w[i - 2], 19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    a = s->h[0], b = s->h[1], c = s->h[2], d = s->h[3], e = s->h[4], f = s->h[5], g = s->h[6], h = s->h[7];
    for (int i = 0; i < 64; i++) {
        uint32_t S1 = ROR32(e, 6) ^ ROR32(e, 11) ^ ROR32(e, 25);
        uint32_t ch = (e & f) ^ (~e & g);
        uint32_t t1 = h + S1 + ch + sha256_k[i] + w[i];
        uint32_t S0 = ROR32(a, 2) ^ ROR32(a, 13) ^ ROR32(a, 22);
        uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t t2 = S0 + maj;
        h = g, g = f, f = e, e = d + t1, d = c, c = b, b = a, a = t1 + t2;
    }
    s->h[0] += a, s->h[1] += b, s->h[2] += c, s->h[3] += d;
    s->h[4] += e, s->h[5] += f, s->h[6] += g, s->h[7] += h;
}

static void sha256_init(sha256_ctx *s)
{
    static const uint32_t h0[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                                   0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    memcpy(s->h, h0, sizeof(h0));
    s->bytes = 0;
    s->used = 0;
}

static void sha256_update(sha256_ctx *s, const unsigned char *p, size_t n)
{
    s->bytes += n;
    while (n) {
        size_t take = 64 - s->used < n ? 64 - s->used : n;
        memcpy(s->block + s->used, p, take);
        s->used += take, p += take, n -= take;
        if (s->used == 64) {
            sha256_block(s, s->block);
            s->used = 0;
        }
    }
}

static void sha256_final(sha256_ctx *s, char hex[65])
{
    uint64_t bits = s->bytes * 8;
    unsigned char pad = 0x80, zero = 0, len[8];
    sha256_update(s, &pad, 1);
    while (s->used != 56) sha256_update(s, &zero, 1);
    for (int i = 0; i < 8; i++) len[i] = (unsigned char)(bits >> (56 - 8 * i));
    sha256_update(s, len, 8);
    for (int i = 0; i < 8; i++) snprintf(hex + 8 * i, 9, "%08x", s->h[i]);
}

/* ------------------------------------------------------------ archive */

static int read_file(const char *path, unsigned char **out, size_t *n, char *msg, size_t cap)
{
    FILE *f = ctools_fopen(path, "rb");
    if (!f) {
        snprintf(msg, cap, "could not open %s", path);
        return 601;
    }
    ctools_fseek(f, 0, SEEK_END);
    long long size = (long long)ctools_ftell(f);
    ctools_fseek(f, 0, SEEK_SET);
    unsigned char *buf = size > 0 ? malloc((size_t)size) : NULL;
    if (!buf || fread(buf, 1, (size_t)size, f) != (size_t)size) {
        fclose(f);
        free(buf);
        snprintf(msg, cap, size > 0 ? "could not read %s" : "%s is empty", path);
        return 692;
    }
    fclose(f);
    *out = buf;
    *n = (size_t)size;
    return 0;
}

static int gunzip(const unsigned char *in, size_t n, unsigned char **out, size_t *outn, char *msg,
                  size_t cap)
{
    if (n < 18 || in[0] != 0x1f || in[1] != 0x8b || in[2] != 8) {
        snprintf(msg, cap, "the archive is not a gzip file");
        return 610;
    }
    unsigned flg = in[3];
    size_t p = 10;
    if (flg & 4) {
        if (p + 2 > n) return 610;
        p += 2 + (size_t)(in[p] | in[p + 1] << 8);
    }
    if (flg & 8)
        while (p < n && in[p++]) {}
    if (flg & 16)
        while (p < n && in[p++]) {}
    if (flg & 2) p += 2;
    if (p >= n - 8) return 610;
    uint32_t isize = (uint32_t)in[n - 4] | (uint32_t)in[n - 3] << 8 | (uint32_t)in[n - 2] << 16 |
                     (uint32_t)in[n - 1] << 24;
    uint32_t crc = (uint32_t)in[n - 8] | (uint32_t)in[n - 7] << 8 | (uint32_t)in[n - 6] << 16 |
                   (uint32_t)in[n - 5] << 24;
    unsigned char *buf = malloc(isize ? isize : 1);
    struct libdeflate_decompressor *d = libdeflate_alloc_decompressor();
    if (!buf || !d) {
        free(buf);
        if (d) libdeflate_free_decompressor(d);
        snprintf(msg, cap, "out of memory");
        return 920;
    }
    size_t used_in = 0, got = 0;
    enum libdeflate_result r = libdeflate_deflate_decompress_ex(d, in + p, n - p - 8, buf, isize, &used_in, &got);
    libdeflate_free_decompressor(d);
    if (r != LIBDEFLATE_SUCCESS || got != isize || libdeflate_crc32(0, buf, got) != crc) {
        free(buf);
        snprintf(msg, cap, "the archive is damaged");
        return 610;
    }
    *out = buf;
    *outn = got;
    return 0;
}

static uint64_t tar_octal(const unsigned char *p, size_t n)
{
    uint64_t v = 0;
    size_t i = 0;
    while (i < n && (p[i] == ' ' || p[i] == 0)) i++;
    for (; i < n && p[i] >= '0' && p[i] <= '7'; i++) v = v * 8 + (uint64_t)(p[i] - '0');
    return v;
}

static const char *base_name(const char *path)
{
    const char *b = path;
    for (const char *p = path; *p; p++)
        if (*p == '/' || *p == '\\') b = p + 1;
    return b;
}

/* The contents of the regular file whose base name is member. */
static int tar_find(const unsigned char *t, size_t n, const char *member, const unsigned char **data,
                    size_t *len)
{
    char longname[4096] = "";
    size_t p = 0;
    while (p + 512 <= n) {
        const unsigned char *h = t + p;
        if (!h[0]) break;
        uint64_t size = tar_octal(h + 124, 12);
        char type = (char)h[156];
        char name[4096];
        if (longname[0]) {
            snprintf(name, sizeof(name), "%s", longname);
            longname[0] = '\0';
        } else if (!memcmp(h + 257, "ustar", 5) && h[345]) {
            snprintf(name, sizeof(name), "%.155s/%.100s", (const char *)h + 345, (const char *)h);
        } else {
            snprintf(name, sizeof(name), "%.100s", (const char *)h);
        }
        size_t body = p + 512;
        if (body + size > n) return 610;
        if (type == 'L') {
            size_t l = size < sizeof(longname) - 1 ? (size_t)size : sizeof(longname) - 1;
            memcpy(longname, t + body, l);
            longname[l] = '\0';
        } else if (type == 'x') {
            /* pax: "<len> path=<name>\n" records */
            const char *rec = (const char *)t + body, *end = rec + size;
            while (rec < end) {
                long l = strtol(rec, NULL, 10);
                const char *kv = memchr(rec, ' ', (size_t)(end - rec));
                if (l <= 0 || !kv || rec + l > end) break;
                if (!strncmp(kv + 1, "path=", 5)) {
                    size_t vl = (size_t)(rec + l - (kv + 6) - 1);
                    if (vl >= sizeof(longname)) vl = sizeof(longname) - 1;
                    memcpy(longname, kv + 6, vl);
                    longname[vl] = '\0';
                }
                rec += l;
            }
        } else if ((type == '0' || type == 0) && !strcmp(base_name(name), member)) {
            *data = t + body;
            *len = (size_t)size;
            return 0;
        }
        p = body + (size_t)((size + 511) / 512 * 512);
    }
    return 601;
}

int csql_install_extract(const char *archive, const char *sha256_hex, const char *member,
                         const char *dest, char *msg, size_t cap)
{
    unsigned char *gz = NULL, *tar = NULL;
    size_t gzn = 0, tarn = 0;
    int rc = read_file(archive, &gz, &gzn, msg, cap);
    if (rc) return rc;
    if (sha256_hex && *sha256_hex) {
        sha256_ctx s;
        char hex[65];
        sha256_init(&s);
        sha256_update(&s, gz, gzn);
        sha256_final(&s, hex);
        if (strcmp(hex, sha256_hex) != 0) {
            free(gz);
            snprintf(msg, cap, "the archive's SHA-256 checksum is %s, not the expected %s", hex, sha256_hex);
            return 610;
        }
    }
    rc = gunzip(gz, gzn, &tar, &tarn, msg, cap);
    free(gz);
    if (rc) return rc;
    const unsigned char *data = NULL;
    size_t len = 0;
    if (tar_find(tar, tarn, member, &data, &len) != 0) {
        free(tar);
        snprintf(msg, cap, "the archive has no file %s", member);
        return 601;
    }
    char part[4200];
    snprintf(part, sizeof(part), "%s.part", dest);
    FILE *f = ctools_fopen(part, "wb");
    if (!f) {
        free(tar);
        snprintf(msg, cap, "could not create %s", part);
        return 603;
    }
    int ok = fwrite(data, 1, len, f) == len;
    ok = fclose(f) == 0 && ok;
    free(tar);
    if (!ok) {
        ctools_remove(part);
        snprintf(msg, cap, "could not write %s", part);
        return 693;
    }
#ifdef _WIN32
    if (ctools_move_file(part, dest, 1) != 0) {
#else
    if (rename(part, dest) != 0) {
#endif
        ctools_remove(part);
        snprintf(msg, cap, "could not replace %s (is it in use by another Stata?)", dest);
        return 608;
    }
#ifndef _WIN32
    chmod(dest, 0755);
#endif
    return 0;
}

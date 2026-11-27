/*
 * cexport_exp.h
 * printf("%.*e") of a double without the C library, for cexport.
 *
 * The result is the correctly rounded one (round half to even on the exact
 * binary value) that MinGW's, glibc's and macOS's printf all give. It is
 * computed exactly in 128-bit integers when the decimal scaling fits, which
 * covers values from about 1e-22 to 3.4e38 at up to 17 significant digits;
 * otherwise the caller formats with snprintf. MinGW's printf allocates its
 * big numbers from a list guarded by one process-wide lock, so threads that
 * format numbers in parallel queue on it (16 threads ran ten times slower
 * than one).
 */
#ifndef CEXPORT_EXP_H
#define CEXPORT_EXP_H

#include <stdint.h>
#include <string.h>

typedef unsigned __int128 cexport_u128;

static const uint64_t cexport_exp_pow5[28] = {
    1ULL, 5ULL, 25ULL, 125ULL, 625ULL, 3125ULL, 15625ULL, 78125ULL, 390625ULL,
    1953125ULL, 9765625ULL, 48828125ULL, 244140625ULL, 1220703125ULL,
    6103515625ULL, 30517578125ULL, 152587890625ULL, 762939453125ULL,
    3814697265625ULL, 19073486328125ULL, 95367431640625ULL, 476837158203125ULL,
    2384185791015625ULL, 11920928955078125ULL, 59604644775390625ULL,
    298023223876953125ULL, 1490116119384765625ULL, 7450580596923828125ULL
};
static const uint64_t cexport_exp_pow10[20] = {
    1ULL, 10ULL, 100ULL, 1000ULL, 10000ULL, 100000ULL, 1000000ULL, 10000000ULL,
    100000000ULL, 1000000000ULL, 10000000000ULL, 100000000000ULL,
    1000000000000ULL, 10000000000000ULL, 100000000000000ULL,
    1000000000000000ULL, 10000000000000000ULL, 100000000000000000ULL,
    1000000000000000000ULL, 10000000000000000000ULL
};

/* 10^r for r <= 38 */
static inline cexport_u128 cexport_exp_ten(int r)
{
    return r <= 19 ? (cexport_u128)cexport_exp_pow10[r]
                   : (cexport_u128)cexport_exp_pow10[19] * cexport_exp_pow10[r - 19];
}

/* trunc(m * 2^e2 * 10^q) and how the dropped fraction compares with one
 * half (-1 below, 0 equal, 1 above). Returns 0 when the scaling does not fit
 * 128 bits. */
static inline int cexport_exp_scale(uint64_t m, int e2, int q, cexport_u128 *trunc, int *half)
{
    if (q >= 0) {
        if (q > 32) return 0;
        cexport_u128 f = q <= 27 ? (cexport_u128)cexport_exp_pow5[q]
                                 : (cexport_u128)cexport_exp_pow5[27] * cexport_exp_pow5[q - 27];
        cexport_u128 p = (cexport_u128)m * f; /* < 2^53 * 5^32 < 2^128 */
        int sh = e2 + q;
        if (sh >= 0) {
            if (sh >= 64 || (p >> (127 - sh)) != 0) return 0;
            *trunc = p << sh;
            *half = -1;
            return 1;
        }
        int s = -sh;
        if (s >= 127) return 0;
        cexport_u128 rem = p & ((((cexport_u128)1) << s) - 1), hf = ((cexport_u128)1) << (s - 1);
        *trunc = p >> s;
        *half = rem > hf ? 1 : rem == hf ? 0 : -1;
        return 1;
    }
    int r = -q;
    if (r > 38 || e2 < 0 || e2 > 127 - 53) return 0;
    cexport_u128 v = (cexport_u128)m << e2, d = cexport_exp_ten(r);
    cexport_u128 t = v / d, rem = v - t * d;
    *trunc = t;
    *half = rem > d - rem ? 1 : rem == d - rem ? 0 : -1;
    return 1;
}

/* Writes printf("%.*e", decimals, val) to buf (at least 32 bytes) for a
 * finite nonzero normal val and decimals 0..16; estimate is a decimal
 * exponent of val within one of floor(log10(|val|)). Returns the length, or
 * 0 when the value must be formatted by snprintf. */
static inline int cexport_exp_format(double val, int decimals, int estimate, char *buf)
{
    uint64_t bits;
    memcpy(&bits, &val, sizeof(bits));
    int biased = (int)((bits >> 52) & 0x7ff);
    if (biased == 0 || biased == 0x7ff || decimals < 0 || decimals > 16) return 0;
    uint64_t m = (bits & UINT64_C(0x000fffffffffffff)) | (UINT64_C(1) << 52);
    int e2 = biased - 1075;
    const uint64_t lo = cexport_exp_pow10[decimals], hi = cexport_exp_pow10[decimals + 1];
    int k = estimate;
    cexport_u128 t = 0;
    int half = 0, tries = 0;
    for (;;) {
        if (++tries > 3 || !cexport_exp_scale(m, e2, decimals - k, &t, &half)) return 0;
        if (t < lo) k--;            /* val < 10^k */
        else if (t >= hi) k++;      /* val >= 10^(k+1) */
        else break;
    }
    uint64_t n = (uint64_t)t;
    if (half > 0 || (half == 0 && (n & 1))) n++;
    if (n == hi) {                  /* rounded up to the next power of ten */
        n = lo;
        k++;
    }
    char *p = buf;
    if (bits >> 63) *p++ = '-';
    uint64_t lead = n / lo, rest = n - lead * lo;
    *p++ = (char)('0' + lead);
    if (decimals) {
        *p++ = '.';
        for (int i = decimals - 1; i >= 0; i--) {
            p[i] = (char)('0' + rest % 10);
            rest /= 10;
        }
        p += decimals;
    }
    *p++ = 'e';
    *p++ = k < 0 ? '-' : '+';
    int a = k < 0 ? -k : k;
    if (a >= 100) *p++ = (char)('0' + a / 100);
    *p++ = (char)('0' + a / 10 % 10);
    *p++ = (char)('0' + a % 10);
    *p = '\0';
    return (int)(p - buf);
}

#endif /* CEXPORT_EXP_H */

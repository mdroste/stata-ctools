/*
 * ctools_eisel_lemire.h
 *
 * Fast float reconstruction from parsed mantissa + exponent.
 * Uses the Clinger fast path when the result is exact:
 *   - mantissa fits in 53 bits (< 2^53)
 *   - |exponent| <= 22 (exact powers of 10 in double)
 * and otherwise the Eisel-Lemire algorithm, which either returns the correctly
 * rounded double or reports that the caller must fall back to strtod.
 */

#ifndef CTOOLS_EISEL_LEMIRE_H
#define CTOOLS_EISEL_LEMIRE_H

#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include "ctools_pow10_128.h"

/*
 * Eisel-Lemire: correctly rounded mantissa * 10^exp10 for any 64-bit
 * mantissa, or false when the 128-bit product cannot decide the rounding
 * (or the result is subnormal/infinite). Port of Go's strconv eiselLemire64.
 */
static inline bool ctools_eisel_lemire64(uint64_t man, int exp10, bool negative, double *result)
{
    if (man == 0) {
        *result = negative ? -0.0 : 0.0;
        return true;
    }
    if (exp10 < CTOOLS_POW10_128_MIN || exp10 > CTOOLS_POW10_128_MAX) return false;

    /* floor(exp10 * log2(10)) without relying on signed right shifts */
    int scaled = 217706 * exp10;
    int log2_10 = scaled >= 0 ? scaled >> 16 : -((-scaled + 65535) >> 16);
    int clz = __builtin_clzll(man);
    man <<= clz;
    uint64_t ret_exp2 = (uint64_t)(log2_10 + 64 + 1023) - (uint64_t)clz;

    const uint64_t *pow10 = ctools_pow10_128[exp10 - CTOOLS_POW10_128_MIN];
    __uint128_t x = (__uint128_t)man * pow10[1];
    uint64_t x_hi = (uint64_t)(x >> 64), x_lo = (uint64_t)x;

    /* Wider approximation when the truncated product is ambiguous */
    if ((x_hi & 0x1FF) == 0x1FF && x_lo + man < man) {
        __uint128_t y = (__uint128_t)man * pow10[0];
        uint64_t y_hi = (uint64_t)(y >> 64), y_lo = (uint64_t)y;
        uint64_t merged_hi = x_hi, merged_lo = x_lo + y_hi;
        if (merged_lo < x_lo) merged_hi++;
        if ((merged_hi & 0x1FF) == 0x1FF && merged_lo + 1 == 0 && y_lo + man < man) return false;
        x_hi = merged_hi;
        x_lo = merged_lo;
    }

    uint64_t msb = x_hi >> 63;
    uint64_t ret_man = x_hi >> (msb + 9);
    ret_exp2 -= 1 ^ msb;

    /* Exact halfway case: needs the full decimal input to break the tie */
    if (x_lo == 0 && (x_hi & 0x1FF) == 0 && (ret_man & 3) == 1) return false;

    ret_man += ret_man & 1;
    ret_man >>= 1;
    if (ret_man >> 53 > 0) {
        ret_man >>= 1;
        ret_exp2 += 1;
    }
    /* Subnormal, zero, infinite or NaN results are left to strtod */
    if (ret_exp2 - 1 >= 0x7FF - 1) return false;

    uint64_t bits = ret_exp2 << 52 | (ret_man & 0x000FFFFFFFFFFFFFULL);
    if (negative) bits |= 0x8000000000000000ULL;
    memcpy(result, &bits, sizeof(bits));
    return true;
}

/* Exact powers of 10 representable as doubles (10^0 through 10^22).
 * These are all exactly representable in IEEE 754 double precision. */
static const double exact_pow10[23] = {
    1e0,  1e1,  1e2,  1e3,  1e4,  1e5,  1e6,  1e7,  1e8,  1e9,
    1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19,
    1e20, 1e21, 1e22
};

/*
 * Fast mantissa-to-double conversion using Clinger's algorithm.
 *
 * Returns true on success, false if caller should fall back to strtod.
 * Handles ~99% of typical CSV numeric fields (integers, simple decimals).
 *
 * The mantissa must be nonzero and fit in 64 bits.
 * The exp10 is the decimal exponent: value = mantissa * 10^exp10.
 */
static inline bool eisel_lemire_compute(uint64_t mantissa, int exp10,
                                         bool negative, double *result)
{
    /* Max exact integer in double: 2^53 = 9007199254740992 */
    static const uint64_t MAX_EXACT_INT = (1ULL << 53);

    double value;

    if (mantissa < MAX_EXACT_INT) {
        /* Mantissa is exactly representable as double */
        value = (double)mantissa;

        if (exp10 == 0) {
            /* No scaling needed — most common case for integers */
            *result = negative ? -value : value;
            return true;
        }

        if (exp10 > 0) {
            if (exp10 <= 22) {
                /* Single exact multiply: mantissa * 10^exp */
                value *= exact_pow10[exp10];
                *result = negative ? -value : value;
                return true;
            }
            if (exp10 <= 22 + 15) {
                /* Two-step: first scale mantissa into range, then multiply.
                 * mantissa * 10^(exp10-22) must still be exact (< 2^53),
                 * then multiply by 10^22. */
                int first_exp = exp10 - 22;
                double scaled = value * exact_pow10[first_exp];
                /* Check if scaled is still exactly representable */
                if (scaled < (double)MAX_EXACT_INT) {
                    *result = negative ? -(scaled * exact_pow10[22]) : (scaled * exact_pow10[22]);
                    return true;
                }
            }
            return ctools_eisel_lemire64(mantissa, exp10, negative, result);
        } else {
            /* exp10 < 0: division */
            int neg_exp = -exp10;
            if (neg_exp <= 22) {
                /* Single exact division: mantissa / 10^|exp| */
                value /= exact_pow10[neg_exp];
                *result = negative ? -value : value;
                return true;
            }
            return ctools_eisel_lemire64(mantissa, exp10, negative, result);
        }
    }

    /* mantissa >= 2^53: not exactly representable.
     * If exponent is 0, we can still get the right answer since
     * (double)mantissa rounds correctly. */
    if (exp10 == 0) {
        value = (double)mantissa;
        *result = negative ? -value : value;
        return true;
    }

    return ctools_eisel_lemire64(mantissa, exp10, negative, result);
}

#endif /* CTOOLS_EISEL_LEMIRE_H */

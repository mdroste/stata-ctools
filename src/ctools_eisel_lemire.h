/*
 * ctools_eisel_lemire.h
 *
 * Fast float reconstruction from parsed mantissa + exponent.
 * Uses the Clinger fast path for exact results when possible:
 *   - mantissa fits in 53 bits (< 2^53)
 *   - |exponent| <= 22 (exact powers of 10 in double)
 *
 * For larger exponents, uses a two-step multiply approach:
 *   mantissa * 10^22 * 10^(exp-22) or mantissa / 10^22 / 10^(-exp-22)
 *
 * Falls back to strtod for remaining edge cases.
 */

#ifndef CTOOLS_EISEL_LEMIRE_H
#define CTOOLS_EISEL_LEMIRE_H

#include <stdint.h>
#include <stdbool.h>
#include <string.h>

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
            /* Fall through to strtod */
            return false;
        } else {
            /* exp10 < 0: division */
            int neg_exp = -exp10;
            if (neg_exp <= 22) {
                /* Single exact division: mantissa / 10^|exp| */
                value /= exact_pow10[neg_exp];
                *result = negative ? -value : value;
                return true;
            }
            /* For larger negative exponents: two-step division.
             * value / 10^22 / 10^(neg_exp - 22) */
            if (neg_exp <= 22 + 22) {
                value /= exact_pow10[22];
                value /= exact_pow10[neg_exp - 22];
                *result = negative ? -value : value;
                return true;
            }
            return false;
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

    /* For large mantissa with nonzero exponent, fall back to strtod
     * to ensure correct rounding. */
    return false;
}

#endif /* CTOOLS_EISEL_LEMIRE_H */

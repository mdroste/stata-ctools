/* Shared xoshiro256** stream. Keep seeding and draw order stable.
 * uniform(0) returns zero without consuming a draw. Each group owns its state. */
#ifndef CTOOLS_RNG_H
#define CTOOLS_RNG_H
#include <math.h>
#include <stdint.h>
#include <stddef.h>

typedef struct
{
    uint64_t s[4];
} ctools_rng;

static inline uint64_t ctools_rng_rotl64(const uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

static inline uint64_t ctools_rng_next(ctools_rng *state)
{
    const uint64_t result = ctools_rng_rotl64(state->s[1] * 5, 7) * 9;
    const uint64_t t = state->s[1] << 17;

    state->s[2] ^= state->s[0];
    state->s[3] ^= state->s[1];
    state->s[1] ^= state->s[2];
    state->s[0] ^= state->s[3];

    state->s[2] ^= t;
    state->s[3] = ctools_rng_rotl64(state->s[3], 45);

    return result;
}

/* SplitMix64 for seeding */
static inline uint64_t ctools_rng_splitmix64(uint64_t *state)
{
    uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static inline void ctools_rng_seed(ctools_rng *state, uint64_t seed)
{
    uint64_t sm_state = seed;
    state->s[0] = ctools_rng_splitmix64(&sm_state);
    state->s[1] = ctools_rng_splitmix64(&sm_state);
    state->s[2] = ctools_rng_splitmix64(&sm_state);
    state->s[3] = ctools_rng_splitmix64(&sm_state);
}

/* Generate uniform random in [0, n) */
static inline size_t ctools_rng_uniform(ctools_rng *state, size_t n)
{
    if (n == 0)
        return 0;
    /* Rejection sampling to avoid modulo bias */
    uint64_t threshold = (UINT64_MAX - n + 1) % n;
    uint64_t r;
    do
    {
        r = ctools_rng_next(state);
    } while (r < threshold);
    return (size_t)(r % n);
}

/* Uniform double in [0, 1) with 53 random bits. */
static inline double ctools_rng_double(ctools_rng *state)
{
    return (double)(ctools_rng_next(state) >> 11) * 0x1.0p-53;
}

/* Stirling-series correction log(k!) - [(k + 1/2) log(k + 1) - (k + 1) +
 * log(2 pi)/2], as in Hormann's BTRD. */
static inline double ctools_rng_stirling_correction(double k)
{
    static const double table[10] = {
        0.08106146679532726, 0.04134069595540929, 0.02767792568499834,
        0.02079067210376509, 0.01664469118982119, 0.01387612882307075,
        0.01189670994589177, 0.01041126526197209, 0.009255462182712733,
        0.008330563433362871};
    if (k < 10.0)
        return table[(int)k];
    double r = 1.0 / (k + 1.0), r2 = r * r;
    return (1.0 / 12.0 - (1.0 / 360.0 - r2 / 1260.0) * r2) * r;
}

/* Binomial(n, p) draw, exact in distribution. Inversion when n*min(p, 1-p)
 * is below 10, otherwise BTRD (W. Hormann, "The generation of binomial
 * random variates", J. Statist. Comput. Simul. 46, 1993). Returns 0 for
 * n <= 0 or p <= 0 and n for p >= 1. */
static inline int64_t ctools_rng_binomial(ctools_rng *state, int64_t n, double p)
{
    if (n <= 0 || !(p > 0.0))
        return 0;
    if (p >= 1.0)
        return n;
    int flip = p > 0.5;
    double pp = flip ? 1.0 - p : p;
    double dn = (double)n;
    int64_t k;
    if (dn * pp < 10.0)
    {
        /* Inversion by sequential search from 0 */
        double q = 1.0 - pp, s = pp / q, a = (dn + 1.0) * s;
        for (;;)
        {
            double r = exp(dn * log1p(-pp)), u = ctools_rng_double(state);
            k = 0;
            while (u > r)
            {
                u -= r;
                k++;
                if (k > n)
                    break;
                r *= a / (double)k - s;
            }
            if (k <= n)
                break;
        }
    }
    else
    {
        const double spq = sqrt(dn * pp * (1.0 - pp));
        const double b = 1.15 + 2.53 * spq;
        const double a = -0.0873 + 0.0248 * b + 0.01 * pp;
        const double c = dn * pp + 0.5;
        const double alpha = (2.83 + 5.1 / b) * spq;
        const double vr = 0.92 - 4.2 / b;
        const double urvr = 0.86 * vr;
        const double m = floor((dn + 1.0) * pp);
        const double r = pp / (1.0 - pp);
        const double nr = (dn + 1.0) * r;
        const double npq = dn * pp * (1.0 - pp);
        for (;;)
        {
            double u, v = ctools_rng_double(state);
            if (v <= urvr)
            {
                u = v / vr - 0.43;
                double kd = floor((2.0 * a / (0.5 - fabs(u)) + b) * u + c);
                if (kd < 0.0 || kd > dn)
                    continue;
                k = (int64_t)kd;
                break;
            }
            if (v >= vr)
            {
                u = ctools_rng_double(state) - 0.5;
            }
            else
            {
                u = v / vr - 0.93;
                u = (u < 0.0 ? -0.5 : 0.5) - u;
                v = ctools_rng_double(state) * vr;
            }
            double us = 0.5 - fabs(u);
            double kd = floor((2.0 * a / us + b) * u + c);
            if (kd < 0.0 || kd > dn)
                continue;
            k = (int64_t)kd;
            v = v * alpha / (a / (us * us) + b);
            double km = fabs(kd - m);
            if (km <= 15.0)
            {
                /* Recursive evaluation of f(k) / f(m) */
                double f = 1.0;
                if (m < kd)
                {
                    for (double i = m + 1.0; i <= kd; i += 1.0)
                        f *= nr / i - r;
                }
                else if (m > kd)
                {
                    for (double i = kd + 1.0; i <= m; i += 1.0)
                        v *= nr / i - r;
                }
                if (v <= f)
                    break;
                continue;
            }
            /* Squeeze, then the exact log-density comparison */
            v = log(v);
            double rho = (km / npq) * (((km / 3.0 + 0.625) * km + 1.0 / 6.0) / npq + 0.5);
            double t = -km * km / (2.0 * npq);
            if (v < t - rho)
                break;
            if (v > t + rho)
                continue;
            double nm = dn - m + 1.0, nk = dn - kd + 1.0;
            double h = (m + 0.5) * log((m + 1.0) / (r * nm)) +
                       ctools_rng_stirling_correction(m) +
                       ctools_rng_stirling_correction(dn - m);
            if (v <= h + (dn + 1.0) * log(nm / nk) + (kd + 0.5) * log(nk * r / (kd + 1.0)) -
                          ctools_rng_stirling_correction(kd) -
                          ctools_rng_stirling_correction(dn - kd))
                break;
        }
    }
    return flip ? n - k : k;
}
#endif

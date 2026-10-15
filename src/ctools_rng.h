/* Shared xoshiro256** stream. Keep seeding and draw order stable.
 * uniform(0) returns zero without consuming a draw. Each group owns its state. */
#ifndef CTOOLS_RNG_H
#define CTOOLS_RNG_H
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
#endif

#include "cpsmatch_internal.h"
/* Radius weights are double-precision 1/k, 1 <= k <= 2^31-1 (SPI limit).
 * Their lowest possible significand bit is 2^-83. A fixed 84-bit fractional
 * part therefore stores every contribution exactly. At most 2^31-1 treated
 * rows contribute, so two portable uint64_t words have ample headroom.
 * Unsigned carry/borrow also represents negative endpoint totals modulo 2^128. */
typedef struct
{
    uint64_t hi, lo;
} cpsmatch_radius_sum;

static cpsmatch_radius_sum cpsmatch_radius_encode(double weight)
{
    int exponent;
    double fraction = frexp(weight, &exponent);
    uint64_t significand = (uint64_t)ldexp(fraction, 53);
    unsigned shift = (unsigned)(exponent + 31); /* 1..32 for the supported k */
    return (cpsmatch_radius_sum){significand >> (64 - shift), significand << shift};
}

static cpsmatch_radius_sum cpsmatch_radius_add(cpsmatch_radius_sum a, cpsmatch_radius_sum b)
{
    uint64_t lo = a.lo + b.lo;
    return (cpsmatch_radius_sum){a.hi + b.hi + (lo < a.lo), lo};
}

static cpsmatch_radius_sum cpsmatch_radius_sub(cpsmatch_radius_sum a, cpsmatch_radius_sum b)
{
    return (cpsmatch_radius_sum){a.hi - b.hi - (a.lo < b.lo), a.lo - b.lo};
}

static double cpsmatch_radius_value(cpsmatch_radius_sum value)
{
    /* The nonnegative high word is below 2^51 and converts exactly to double.
     * Only conversion of the low word and the final addition can round. */
    return ldexp((double)value.hi, -20) + ldexp((double)value.lo, -84);
}

int cpsmatch_match_radius(const sorted_control_t *controls, size_t n_controls,
                          const double *treated_pscore, const size_t *treated_idx, size_t n_treated,
                          const size_t *control_idx, const perm_idx_t *obs_map, double radius,
                          int common, double common_min, double common_max, double *out_weight,
                          double *out_match, double *out_nn, double *out_support, int *n_matched,
                          int *n_off)
{
    /* One start position per treated row; its existing _nn result gives length.
     * The endpoint array is in sorted-control order, not filtered-row order. */
    perm_idx_t *starts = ctools_safe_malloc2(n_treated, sizeof(*starts));
    cpsmatch_radius_sum *endpoints = ctools_safe_calloc2(n_controls + 1, sizeof(*endpoints));
    if (!starts || !endpoints)
    {
        free(starts);
        free(endpoints);
        return 920;
    }

    int matched = 0, off = 0;
#pragma omp parallel for schedule(static) num_threads(ctools_get_max_threads())                    \
    reduction(+ : matched, off)
    for (size_t t = 0; t < n_treated; t++)
    {
        size_t row = treated_idx[t];
        double ps = treated_pscore[t];
        out_nn[row] = 0.0;
        if (common && (ps < common_min || ps > common_max))
        {
            out_support[row] = 0.0;
            off++;
            continue;
        }
        radius_match_info range;
        cpsmatch_find_radius_range_fast(ps, controls, n_controls, radius, &range);
        if (range.n_matches > 0)
        {
            starts[t] = (perm_idx_t)range.start;
            out_nn[row] = (double)range.n_matches;
            out_match[row] = (double)obs_map[control_idx[controls[range.start].orig_idx]];
            out_support[row] = 1.0;
            out_weight[row] = 1.0;
            matched++;
        }
        else
        {
            out_support[row] = 0.0;
            off++;
        }
    }

    /* Every eligible interval contributes the same 1/_nn to each contained
     * control. Accumulate its endpoints in treated input order, then prefix-sum
     * once instead of visiting all treated/control pairs. */
    for (size_t t = 0; t < n_treated; t++)
    {
        size_t row = treated_idx[t];
        if (out_nn[row] <= 0.0)
            continue;
        size_t start = starts[t], end = start + (size_t)out_nn[row];
        cpsmatch_radius_sum weight = cpsmatch_radius_encode(1.0 / out_nn[row]);
        endpoints[start] = cpsmatch_radius_add(endpoints[start], weight);
        endpoints[end] = cpsmatch_radius_sub(endpoints[end], weight);
    }
    cpsmatch_radius_sum running = {0, 0};
    for (size_t j = 0; j < n_controls; j++)
    {
        running = cpsmatch_radius_add(running, endpoints[j]);
        out_weight[control_idx[controls[j].orig_idx]] = cpsmatch_radius_value(running);
    }
    *n_matched += matched;
    *n_off += off;
    free(starts);
    free(endpoints);
    return 0;
}

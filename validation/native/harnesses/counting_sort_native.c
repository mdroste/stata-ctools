/* Counting eligibility and LSD fallback against an independent stable merge
 * sort. The input permutation is shuffled, so equal keys must preserve its
 * order rather than row-number order. Invalid values at row 1 are deliberately
 * outside the bounded sample on larger cases. */
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"

static ST_plugin host;
ST_plugin *_stata_ = &host;
static int requested_threads = 1;
int ctools_get_max_threads(void)
{
    return requested_threads;
}

static ST_boolean missing(ST_double value)
{
    /* Exercise the actual team at counting's SPI callback, including the
     * 50,000--100,000 row band that previously exceeded threads(1). */
    assert(omp_get_num_threads() <= requested_threads);
    return value >= host.missval;
}

#include "ctools_sort_counting.c"
#include "ctools_sort_radix_lsd.c"

static uint64_t random_state = UINT64_C(0x123456789abcdef);
static uint64_t next_random(void)
{
    random_state ^= random_state << 13;
    random_state ^= random_state >> 7;
    random_state ^= random_state << 17;
    return random_state;
}

static const stata_variable *reference_vars;
static const int *reference_keys;
static size_t reference_nkeys;

static int compare_rows(perm_idx_t a, perm_idx_t b)
{
    for (size_t k = 0; k < reference_nkeys; k++)
    {
        const double *values = reference_vars[reference_keys[k] - 1].data.dbl;
        uint64_t x = ctools_double_to_sortable(values[a]);
        uint64_t y = ctools_double_to_sortable(values[b]);
        if (x != y)
            return x < y ? -1 : 1;
    }
    return 0;
}

static void reference_sort(perm_idx_t *rows, perm_idx_t *scratch, size_t n)
{
    if (n < 2)
        return;
    size_t middle = n / 2, a = 0, b = middle, out = 0;
    reference_sort(rows, scratch, middle);
    reference_sort(rows + middle, scratch, n - middle);
    while (a < middle && b < n)
    {
        scratch[out++] = compare_rows(rows[b], rows[a]) < 0 ? rows[b++] : rows[a++];
    }
    while (a < middle)
        scratch[out++] = rows[a++];
    while (b < n)
        scratch[out++] = rows[b++];
    memcpy(rows, scratch, n * sizeof(*rows));
}

enum fixture
{
    SMALL_INTEGERS,
    ALL_MISSING,
    CONTINUOUS,
    UNSAMPLED_FRACTION,
    UNSAMPLED_EXTENDED_MISSING,
    UNSAMPLED_WIDE_RANGE,
    SAMPLED_WIDE_RANGE,
    EXACT_MAGNITUDE_LIMIT,
    RANGE_LIMIT_INSIDE,
    RANGE_LIMIT_OUTSIDE
};

static void run_case(size_t n, size_t nkeys, int team, enum fixture fixture, size_t key_position)
{
    stata_variable vars[3] = {0};
    int keys[3] = {2, 1, 3};
    stata_data data = {.nobs = n, .nvars = 3, .vars = vars};
    perm_idx_t *expect = malloc(n * sizeof(*expect));
    perm_idx_t *scratch = malloc(n * sizeof(*scratch));
    data.sort_order = malloc(n * sizeof(*data.sort_order));
    assert(expect && scratch && data.sort_order);
    for (size_t j = 0; j < 3; j++)
    {
        vars[j].type = STATA_TYPE_DOUBLE;
        vars[j].nobs = n;
        vars[j].data.dbl = malloc(n * sizeof(double));
        assert(vars[j].data.dbl);
        for (size_t i = 0; i < n; i++)
        {
            double x = (double)(next_random() % 63) - 31;
            if (i % 29 == 0)
                x = host.missval;
            if (i % 37 == 0)
                x = -0.0;
            vars[j].data.dbl[i] = x;
        }
    }
    double *values = vars[keys[key_position] - 1].data.dbl;
    if (fixture == ALL_MISSING || fixture == EXACT_MAGNITUDE_LIMIT)
    {
        for (size_t i = 0; i < n; i++)
        {
            values[i] = fixture == ALL_MISSING ? host.missval : 9007199254740992.0;
        }
    }
    else if (fixture == CONTINUOUS)
    {
        for (size_t i = 0; i < n; i++)
            values[i] = (double)(i % 127) + 0.25;
    }
    else if (fixture == UNSAMPLED_FRACTION)
    {
        values[1] = 0.25;
    }
    else if (fixture == UNSAMPLED_EXTENDED_MISSING)
    {
        uint64_t bits = UINT64_C(0x7fe0010000000000);
        memcpy(values + 1, &bits, sizeof(bits));
    }
    else if (fixture == UNSAMPLED_WIDE_RANGE)
    {
        values[1] = 1000000.0;
    }
    else if (fixture == SAMPLED_WIDE_RANGE || fixture == RANGE_LIMIT_INSIDE ||
             fixture == RANGE_LIMIT_OUTSIDE)
    {
        for (size_t i = 0; i < n; i++)
            values[i] = 0;
        size_t at = fixture == SAMPLED_WIDE_RANGE ? n / 64 + 1 : 1;
        values[at] = fixture == RANGE_LIMIT_INSIDE ? 999999.0 : 1000000.0;
    }

    for (size_t i = 0; i < n; i++)
        data.sort_order[i] = (perm_idx_t)i;
    for (size_t i = n; i > 1; i--)
    {
        size_t j = next_random() % i;
        perm_idx_t row = data.sort_order[i - 1];
        data.sort_order[i - 1] = data.sort_order[j];
        data.sort_order[j] = row;
    }
    memcpy(expect, data.sort_order, n * sizeof(*expect));
    reference_vars = vars;
    reference_keys = keys;
    reference_nkeys = nkeys;
    reference_sort(expect, scratch, n);

    requested_threads = team;
    int eligible = fixture == SMALL_INTEGERS || fixture == ALL_MISSING ||
                   fixture == EXACT_MAGNITUDE_LIMIT || fixture == RANGE_LIMIT_INSIDE;
    stata_retcode rc = ctools_sort_counting_order_only(&data, keys, nkeys);
    assert(rc == (eligible ? STATA_OK : STATA_ERR_UNSUPPORTED_TYPE));
    if (rc == STATA_ERR_UNSUPPORTED_TYPE)
    {
        assert(ctools_sort_radix_lsd_order_only(&data, keys, nkeys) == STATA_OK);
    }
    if (memcmp(expect, data.sort_order, n * sizeof(*expect)) != 0)
    {
        fprintf(stderr, "counting mismatch: n=%zu keys=%zu team=%d fixture=%d position=%zu\n", n,
                nkeys, team, fixture, key_position);
        abort();
    }
    for (size_t j = 0; j < 3; j++)
        free(vars[j].data.dbl);
    free(data.sort_order);
    free(expect);
    free(scratch);
}

int main(void)
{
    uint64_t system_missing = UINT64_C(0x7fe0000000000000);
    memcpy(&host.missval, &system_missing, sizeof(system_missing));
    host.ismissing = missing;
    const size_t sizes[] = {7, 49999, 50000, 99999, 100000, 200001};
    const int teams[] = {1, 2, 4};
    for (size_t s = 0; s < sizeof(sizes) / sizeof(*sizes); s++)
    {
        for (size_t t = 0; t < sizeof(teams) / sizeof(*teams); t++)
        {
            for (int fixture = SMALL_INTEGERS; fixture <= EXACT_MAGNITUDE_LIMIT; fixture++)
            {
                run_case(sizes[s], 1, teams[t], (enum fixture)fixture, 0);
                run_case(sizes[s], 3, teams[t], (enum fixture)fixture, 2);
            }
        }
    }
    run_case(257, 1, 1, RANGE_LIMIT_INSIDE, 0);
    run_case(257, 1, 1, RANGE_LIMIT_OUTSIDE, 0);
    run_case(200001, 3, 4, UNSAMPLED_FRACTION, 0);
    return 0;
}

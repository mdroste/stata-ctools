/* Sorting tests. Select one CTOOLS_CASE_* per compilation.
 * Separate selectors preserve allocator/SPI interposition and independent main functions.
 */

#if defined(CTOOLS_CASE_counting_sort_native)
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
#endif

#if defined(CTOOLS_CASE_csort_ordered_native)
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "csort/csort_impl.c"

static int dispatch_calls;
static int dispatch_permute;
static stata_retcode dispatch_status;
static int test_thread_limit = 1;

int ctools_get_max_threads(void)
{
    return test_thread_limit;
}

/* The serial screens may cover the whole input or leave a short middle.
 * Every possible one-pair inversion must be found at those size boundaries. */
static void check_screen_boundaries(void)
{
    double numbers[1090];
    for (size_t i = 0; i < 1090; i++) numbers[i] = (double)i;
    stata_variable var = {.type = STATA_TYPE_DOUBLE, .data.dbl = numbers};
    stata_data data = {.nvars = 1, .vars = &var};
    int key = 1;
    const size_t sizes[] = {0, 1, 1024, 1025, 1026, 1088, 1089, 1090};
    for (size_t j = 0; j < sizeof(sizes) / sizeof(sizes[0]); j++) {
        data.nobs = sizes[j];
        assert(csort_keys_are_ordered(&data, &key, 1));
        for (size_t i = 1; i < data.nobs; i++) {
            numbers[i] = numbers[i - 1] - 0.5;
            assert(!csort_keys_are_ordered(&data, &key, 1));
            numbers[i] = (double)i;
        }
    }
}

/* More than 200k rows enter the parallel policy. All leading strings tie,
 * so every comparison must inspect the secondary numeric key. Check both a
 * full requested team and a request constrained by the minimum work policy. */
static void check_large_order(size_t nobs, int requested_threads)
{
    double *numbers = malloc(nobs * sizeof(*numbers));
    char **strings = malloc(nobs * sizeof(*strings));
    assert(numbers && strings);
    for (size_t i = 0; i < nobs; i++) {
        numbers[i] = (double)i;
        strings[i] = "shared prefix";
    }
    stata_variable vars[3] = {0};
    vars[1].type = STATA_TYPE_DOUBLE;
    vars[1].data.dbl = numbers;
    vars[2].type = STATA_TYPE_STRING;
    vars[2].data.str = strings;
    stata_data data = {.nobs = nobs, .nvars = 3, .vars = vars};
    int keys[] = {3, 2};
    test_thread_limit = requested_threads;
    assert(csort_keys_are_ordered(&data, keys, 2));

    /* Prefix and tail seams, plus interior and final-row inversions. */
    const size_t positions[] = {1, 1024, 1025, 1026, nobs / 2,
                               nobs - 65, nobs - 64, nobs - 63, nobs - 1};
    for (size_t j = 0; j < sizeof(positions) / sizeof(positions[0]); j++) {
        size_t i = positions[j];
        numbers[i] = numbers[i - 1] - 0.5;
        assert(!csort_keys_are_ordered(&data, keys, 2));
        numbers[i] = (double)i;
    }

    /* Place an inversion at each worker's first pair and on either side.
     * The expected result is independent of which worker detects it. */
    size_t remaining = nobs - 1025 - 64;
    size_t team = remaining / MIN_OBS_PER_THREAD;
    if (team > (size_t)requested_threads) team = (size_t)requested_threads;
    for (size_t worker = 1; worker < team; worker++) {
        size_t base = remaining / team, extra = remaining % team;
        size_t boundary = 1025 + worker * base + (worker < extra ? worker : extra);
        for (size_t i = boundary - 1; i <= boundary + 1; i++) {
            numbers[i] = numbers[i - 1] - 0.5;
            assert(!csort_keys_are_ordered(&data, keys, 2));
            numbers[i] = (double)i;
        }
    }
    assert(csort_keys_are_ordered(&data, keys, 2));
    free(strings);
    free(numbers);
}

/* Isolate the csort adapter's decision from the engines, which have their
 * own stable-reference suites. An explicit algorithm must always reach here. */
stata_retcode ctools_sort_dispatch(stata_data *data, int *sort_vars, size_t nsort,
                                  sort_algorithm_t algorithm)
{
    (void)sort_vars;
    (void)nsort;
    (void)algorithm;
    dispatch_calls++;
    if (dispatch_status != STATA_OK)
        return dispatch_status;
    if (dispatch_permute) {
        data->sort_order[0] = 1;
        data->sort_order[1] = 0;
    }
    return STATA_OK;
}

int main(void)
{
    /* Stata's finite missing encodings: . < .a < .z. */
    double number[] = {-2, -1, -0.0, 0.0, 5, 0x1p1023, 0x1.001p1023, 0x1.01ap1023};
    char *text[] = {"z", "a", "", "\x80", "a", "a", "a", "a"};
    perm_idx_t order[] = {0, 1, 2, 3, 4, 5, 6, 7};
    stata_variable vars[3] = {0};
    vars[1].type = STATA_TYPE_STRING;
    vars[1].data.str = text;
    vars[2].type = STATA_TYPE_DOUBLE;
    vars[2].data.dbl = number;
    stata_data data = {.nobs = 8, .nvars = 3, .vars = vars, .sort_order = order};
    int keys[] = {3, 2};
    int identity = 0;
    unsigned char before[sizeof(number)];
    memcpy(before, number, sizeof(number));

    assert(csort_compute_order(&data, keys, 2, SORT_ALG_AUTO, &identity) == STATA_OK);
    assert(identity && dispatch_calls == 0);
    assert(memcmp(before, number, sizeof(number)) == 0); /* retain signed zero bits */
    for (int alg = SORT_ALG_LSD; alg <= SORT_ALG_IPS4O; alg++) {
        int calls_before = dispatch_calls;
        assert(csort_compute_order(&data, keys, 2, (sort_algorithm_t)alg, &identity) == STATA_OK);
        assert(identity && dispatch_calls == calls_before + 1);
    }

    /* A secondary-key inversion inside equal signed zeros must not skip sort. */
    text[2] = "\xff";
    text[3] = "\x80";
    assert(!csort_keys_are_ordered(&data, keys, 2));
    dispatch_permute = 1;
    assert(csort_compute_order(&data, keys, 2, SORT_ALG_AUTO, &identity) == STATA_OK);
    assert(!identity);
    text[2] = "";
    order[0] = 0;
    order[1] = 1;

    /* Late inversions, distinct missings, and decreasing primary keys. */
    number[7] = 0x1p1023;
    assert(!csort_keys_are_ordered(&data, keys, 2));
    number[7] = 0x1.01ap1023;
    number[1] = -3;
    assert(!csort_keys_are_ordered(&data, keys, 2));
    number[1] = -1;

    /* Exact prefix and unsigned-byte order for a string-only key. */
    char *ordered_text[] = {NULL, "", "a", "a", "aa", "b", "\x80", "\xff"};
    vars[1].data.str = ordered_text;
    int string_key[] = {2};
    assert(csort_keys_are_ordered(&data, string_key, 1));
    ordered_text[7] = "\x7f";
    assert(!csort_keys_are_ordered(&data, string_key, 1));

    /* An engine failure must not inspect an unavailable order or claim identity. */
    dispatch_status = STATA_ERR_MEMORY;
    data.sort_order = NULL;
    assert(csort_compute_order(&data, string_key, 1, SORT_ALG_AUTO, &identity) == STATA_ERR_MEMORY);
    assert(!identity);
    data.nobs = 0;
    assert(csort_keys_are_ordered(&data, keys, 2));
    data.nobs = 1;
    assert(csort_keys_are_ordered(&data, keys, 2));
    check_screen_boundaries();
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    const int thread_limits[] = {1, 2, 8};
    for (size_t j = 0; j < sizeof(thread_limits) / sizeof(thread_limits[0]); j++) {
        check_large_order(400127, thread_limits[j]);
        check_large_order(1000003, thread_limits[j]);
    }
    puts("PASS csort ordered-input adapter");
    return 0;
}
#endif

#if defined(CTOOLS_CASE_msd_sort_native)
/* Reuse the independent stable reference and randomized mixed-key fixtures.
 * The original LSD suite is unchanged; this translation unit runs them on MSD. */
#include "ctools_config.h"
#define ctools_sort_radix_lsd_order_only ctools_sort_radix_msd_order_only
#define main randomized_main
#include "sort_reference.h"
#undef main
#undef ctools_sort_radix_lsd_order_only

static void long_strings(size_t n, int shape, int team, int numeric_first)
{
    stata_variable vars[2] = {0};
    int keys[2] = {numeric_first ? 1 : 2, numeric_first ? 2 : 1};
    stata_data data = {.nobs = n, .nvars = 2, .vars = vars};
    perm_idx_t *expect = malloc(n * sizeof(*expect));
    perm_idx_t *scratch = malloc(n * sizeof(*scratch));
    data.sort_order = malloc(n * sizeof(*data.sort_order));
    vars[0].type = STATA_TYPE_DOUBLE;
    vars[0].data.dbl = malloc(n * sizeof(double));
    vars[1].type = STATA_TYPE_STRING;
    vars[1].str_maxlen = shape == 2 ? 576 : 64;
    vars[1].data.str = calloc(n, sizeof(char *));
    assert(expect && scratch && data.sort_order && vars[0].data.dbl && vars[1].data.str);

    for (size_t i = 0; i < n; i++)
    {
        vars[0].data.dbl[i] = i % 29 == 0 ? -0.0 : (double)(i % 17) - 8;
        if (i % 37 == 0)
        {
            uint64_t bits = UINT64_C(0x7fe0000000000000) + ((uint64_t)(i % 27) << 40);
            memcpy(vars[0].data.dbl + i, &bits, sizeof(bits));
        }
        data.sort_order[i] = (perm_idx_t)i;
        if (i % 113 == 0)
            continue;
        size_t len = vars[1].str_maxlen;
        if (shape == 3 && i % 2 == 0)
            len = (i / 2) % 40;
        else if (i % 79 == 0)
            len = i % 13;
        char *s = malloc(len + 1);
        assert(s);
        for (size_t b = 0; b < len; b++)
        {
            if ((shape == 2 && b < 512) || (shape == 3 && b < 40))
            {
                s[b] = 'P';
            }
            else
            {
                /* A small deterministic pool supplies duplicates and high
                 * bytes; each allocation still ends exactly at its NUL. */
                s[b] = (char)(1 + ((i % 127) * 19 + b * 7) % 255);
            }
        }
        if (len && shape < 2)
            s[0] = (char)('A' + i % (shape == 0 ? 26 : 2));
        s[len] = 0;
        vars[1].data.str[i] = s;
    }
    for (size_t i = n; i > 1; i--)
    {
        size_t j = next() % i;
        perm_idx_t row = data.sort_order[i - 1];
        data.sort_order[i - 1] = data.sort_order[j];
        data.sort_order[j] = row;
    }
    memcpy(expect, data.sort_order, n * sizeof(*expect));
    ref_vars = vars;
    ref_keys = keys;
    ref_nkeys = 2;
    merge_sort(expect, scratch, n);
    threads = team;
    omp_set_num_threads(team);
    assert(ctools_sort_radix_msd_order_only(&data, keys, 2) == STATA_OK);
    if (memcmp(expect, data.sort_order, n * sizeof(*expect)) != 0)
    {
        fprintf(stderr, "MSD mismatch n=%zu shape=%d team=%d numeric_first=%d\n", n, shape, team,
                numeric_first);
        abort();
    }
    for (size_t i = 0; i < n; i++)
        free(vars[1].data.str[i]);
    free(vars[1].data.str);
    free(vars[0].data.dbl);
    free(data.sort_order);
    free(expect);
    free(scratch);
}

int main(void)
{
    omp_set_num_threads(4);
    randomized_main();
    const size_t boundaries[] = {1023, 1024, 1025, 9999, 10000, 10001, 100003, 200001};
    for (size_t i = 0; i < sizeof(boundaries) / sizeof(*boundaries); i++)
    {
        long_strings(boundaries[i], (int)(i % 2), 4, (int)(i % 2));
    }
    for (int team = 1; team <= 4; team *= 2)
    {
        long_strings(10001, 2, team, 0); /* 512-byte common prefix */
        long_strings(20001, 3, team, 1); /* depth-limited merge fallback */
    }
    return 0;
}
#endif

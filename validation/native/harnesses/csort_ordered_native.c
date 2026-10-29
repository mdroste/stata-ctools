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

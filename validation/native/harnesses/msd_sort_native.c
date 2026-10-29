/* Reuse the independent stable reference and randomized mixed-key fixtures.
 * The original LSD suite is unchanged; this translation unit runs them on MSD. */
#include "ctools_config.h"
#define ctools_sort_radix_lsd_order_only ctools_sort_radix_msd_order_only
#define main randomized_main
#include "sort_keys_native.c"
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

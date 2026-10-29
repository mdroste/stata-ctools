#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ctools_types.h"
#include "cpsmatch/cpsmatch_stats.c"

static int fail_boundary;
static void *boundary_malloc(size_t n)
{
    return fail_boundary ? NULL : malloc(n);
}
#define malloc boundary_malloc
#include "ctools_order.c"
#undef malloc

static ST_plugin host;
ST_plugin *_stata_ = &host;
int ctools_get_max_threads(void)
{
    return 4;
}
static int scalar_failure, scalar_calls;
static ST_retcode scalar_save(char *name, double value)
{
    (void)name;
    (void)value;
    return ++scalar_calls == scalar_failure ? 459 : 0;
}

static void statistics(void)
{
    double treatment[] = {0, 0, 0, 1, 1, 1};
    double y[] = {1, 2, 4, 3, 6, 9};
    double weight[] = {2, 1, .5, 1, 1, 1};
    double support[] = {1, 1, 1, 1, 1, 0};
    cpsmatch_statistics s;
    cpsmatch_compute_statistics(6, treatment, y, weight, support, &s);
    assert(s.y1 == 4.5 && s.n1 == 2 && s.var1 == 4.5);
    assert(fabs(s.y0 - 12.0 / 7) < 1e-15 && s.n0 == 3);
    assert(fabs(s.var0 - 7.0 / 3) < 1e-15 && s.weight_squared == 5.25);
    assert(s.y1_u == 6 && s.var1_u == 9);
    assert(s.n0_on == 3 && s.n1_on == 2 && s.n1_off == 1 && s.n0_off == 0);
    for (scalar_failure = 1; scalar_failure <= 17; scalar_failure++)
    {
        scalar_calls = 0;
        assert(cpsmatch_save_statistics(&s) == 459);
        assert(scalar_calls == scalar_failure);
    }
    scalar_failure = 0;
    scalar_calls = 0;
    assert(cpsmatch_save_statistics(&s) == 0 && scalar_calls == 17);
    cpsmatch_compute_statistics(6, treatment, NULL, weight, support, &s);
    assert(s.y1 == host.missval && s.var1 == host.missval && s.n1 == 0);
    assert(s.n1_on == 2 && s.n1_off == 1);
    cpsmatch_compute_statistics(0, NULL, NULL, NULL, NULL, &s);
    assert(s.y0 == host.missval && s.y1 == host.missval && s.n0_on == 0);
    cpsmatch_compute_statistics(1, treatment, y, weight, support, &s);
    assert(s.var0 == host.missval && s.var0_u == host.missval);
}

static void grouping(void)
{
    const size_t n = 70001;
    stata_variable columns[2] = {0};
    stata_data d = {0};
    d.nobs = n;
    d.nvars = 2;
    d.vars = columns;
    columns[0].type = STATA_TYPE_STRING;
    columns[1].type = STATA_TYPE_DOUBLE;
    char **strings = malloc(n * sizeof(*strings));
    double *numbers = malloc(n * sizeof(*numbers));
    perm_idx_t *order = malloc(n * sizeof(*order));
    size_t *expected = malloc((n + 1) * sizeof(*expected)),
           *actual = malloc((n + 1) * sizeof(*actual));
    assert(strings && numbers && order && expected && actual);
    columns[0].data.str = strings;
    columns[1].data.dbl = numbers;
    int keys[] = {0, 1};
    for (int shape = 0; shape < 4; shape++)
    {
        for (size_t i = 0; i < n; i++)
        {
            strings[i] = shape == 0 ? "same" : i < 32768 ? "prefix_a" : "prefix_b";
            numbers[i] = shape < 2 ? 0 : shape == 2 ? (double)(i / 31) : (double)i;
            if (shape == 2 && i >= 65536)
                numbers[i] = host.missval;
            order[i] = (perm_idx_t)(n - i - 1);
        }
        for (int indirect = 0; indirect < 2; indirect++)
        {
            const perm_idx_t *p = indirect ? order : NULL;
            size_t groups = ctools_group_starts(&d, p, keys, 2, expected);
            for (fail_boundary = 0; fail_boundary < 2; fail_boundary++)
            {
                size_t got = ctools_group_starts_parallel(&d, p, keys, 2, actual);
                assert(got == groups && !memcmp(expected, actual, (groups + 1) * sizeof(*actual)));
            }
        }
    }
    fail_boundary = 0;
    free(strings);
    free(numbers);
    free(order);
    free(expected);
    free(actual);
}

int main(void)
{
    host.missval = 8.98846567431158e307;
    host.scalsave = scalar_save;
    statistics();
    grouping();
    puts("MATCH_WINSOR_OPTIMIZATIONS_NATIVE_COMPLETE");
    return 0;
}

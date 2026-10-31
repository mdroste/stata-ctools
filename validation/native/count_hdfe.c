/* Counts, projection identities, and allocation/SPI failure recovery. */
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_simd.h"
#include "alloc.h"
#include <stdarg.h>
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static int nr = 24, nv, write_error, macro_error, request_group;
static double numbers[32][24];
static double saved_n, saved_distinct;
static ST_int nobs(void)
{
    return nr;
}
static ST_int nvars(void)
{
    return nv;
}
static ST_boolean isstr(ST_int v)
{
    assert(v >= 1 && v <= nv);
    return 0;
}
static ST_boolean missing(double x)
{
    return x >= mock.missval;
}
static ST_retcode putnum(ST_int v, ST_int i, double x)
{
    assert(v >= 1 && v <= nv && i >= 1 && i <= nr);
    if (write_error)
        return 459;
    numbers[v - 1][i - 1] = x;
    return 0;
}
static ST_retcode save_macro(char *key, char *value)
{
    if (macro_error)
        return 459;
    if (!strcmp(key, "_ctools_count_N"))
        saved_n = atof(value);
    if (!strcmp(key, "_ctools_count_distinct"))
        saved_distinct = atof(value);
    return 0;
}
static ST_retcode message(char *s)
{
    (void)s;
    return 0;
}
int ctools_get_max_threads(void)
{
    return 1;
}
double ctools_timer_seconds(void)
{
    return 0;
}
void ctools_msg(const char *module, const char *fmt, ...)
{
    (void)module;
    (void)fmt;
}
void ctools_error(const char *module, const char *fmt, ...)
{
    (void)module;
    (void)fmt;
}
void ctools_filtered_data_init(ctools_filtered_data *fd)
{
    memset(fd, 0, sizeof(*fd));
}
void ctools_filtered_data_free(ctools_filtered_data *fd)
{
    for (size_t j = 0; j < fd->data.nvars; j++)
        free(fd->data.vars[j].data.dbl);
    free(fd->data.vars);
    free(fd->data.sort_order);
    free(fd->obs_map);
    memset(fd, 0, sizeof(*fd));
}
stata_retcode ctools_data_load(ctools_filtered_data *fd, int *vars, size_t n, size_t start,
                               size_t end, int flags)
{
    (void)start;
    (void)end;
    fd->data.vars = calloc(n, sizeof(*fd->data.vars));
    if (!fd->data.vars)
        return STATA_ERR_MEMORY;
    fd->data.nvars = n;
    size_t count = 0;
    for (int i = 0; i < nr; i++)
        if (!(flags & CTOOLS_LOAD_MARKER) || (flags & CTOOLS_LOAD_SKIP_IF) || numbers[nv - 1][i])
            count++;
    fd->data.nobs = count;
    if (!(flags & CTOOLS_LOAD_NO_SORT_ORDER))
        fd->data.sort_order = malloc((size_t)nr * sizeof(perm_idx_t));
    fd->obs_map = malloc((size_t)nr * sizeof(perm_idx_t));
    if ((!(flags & CTOOLS_LOAD_NO_SORT_ORDER) && !fd->data.sort_order) || !fd->obs_map)
        return STATA_ERR_MEMORY;
    size_t row = 0;
    for (int i = 0; i < nr; i++)
    {
        if ((flags & CTOOLS_LOAD_MARKER) && !(flags & CTOOLS_LOAD_SKIP_IF) && !numbers[nv - 1][i])
            continue;
        if (fd->data.sort_order)
            fd->data.sort_order[row] = row;
        fd->obs_map[row++] = i + 1;
    }
    for (size_t j = 0; j < n; j++)
    {
        stata_variable *v = &fd->data.vars[j];
        v->type = STATA_TYPE_DOUBLE;
        v->nobs = count;
        v->data.dbl = malloc((size_t)nr * sizeof(double));
        if (!v->data.dbl)
            return STATA_ERR_MEMORY;
        for (size_t i = 0; i < count; i++)
            v->data.dbl[i] = numbers[vars[j] - 1][fd->obs_map[i] - 1];
    }
    return STATA_OK;
}
stata_retcode ctools_store_filtered(double *v, size_t n, int col, perm_idx_t *map)
{
    for (size_t i = 0; i < n; i++)
        if (putnum(col, map[i], v[i]))
            return STATA_ERR_STATA_WRITE;
    return STATA_OK;
}
#include "ctools_parse.c"
#include "ctools_order.c"
#include "ctools_hdfe_utils.c"
#include "hdfe/hdfe_remap.c"
#include "hdfe/hdfe_sample.c"
#include "hdfe/hdfe_solver.c"
#include "cdistinct/cdistinct_impl.c"
#include "chdfe/chdfe_impl.c"

static int distinct_case(void)
{
    nv = 2;
    for (int i = 0; i < nr; i++)
    {
        numbers[0][i] = i % 3;
        numbers[1][i] = i < 12;
    }
    int rc = cdistinct_main("nkeys=1");
    if (!rc)
        assert(saved_n == 12 && saved_distinct == 3);
    return rc;
}

static int distinct_fallback_case(void)
{
    /* Descending unique keys exhaust the bounded trial table. */
    double x[5000];
    for (int i = 0; i < 5000; i++)
        x[i] = 5000 - i;
    stata_variable var = {.type = STATA_TYPE_DOUBLE, .nobs = 5000};
    var.data.dbl = x;
    stata_data data = {.vars = &var, .nvars = 1, .nobs = 5000};
    int key = 0;
    size_t total, distinct;
    int rc = cdistinct_fast_count(&data, &key, 1, 0, &total, &distinct);
    if (!rc)
        assert(total == 5000 && distinct == 5000);
    free(data.sort_order);
    return rc;
}

static void distinct_keys(void)
{
    double x[] = {0, -0.0, mock.missval, 1, 1, -0.0};
    char *s[] = {"", "", "", "a", "ab", ""};
    stata_variable vars[2] = {{.type = STATA_TYPE_DOUBLE, .nobs = 6},
                              {.type = STATA_TYPE_STRING, .nobs = 6}};
    vars[0].data.dbl = x;
    vars[1].data.str = s;
    stata_data data = {.vars = vars, .nvars = 2, .nobs = 6};
    int keys[] = {0, 1};
    size_t total, distinct;
    assert(!cdistinct_fast_count(&data, keys, 2, 1, &total, &distinct));
    assert(total == 6 && distinct == 4);
    assert(!cdistinct_fast_count(&data, keys, 2, 0, &total, &distinct));
    assert(total == 2 && distinct == 2);
    data.nobs = 0;
    assert(!cdistinct_fast_count(&data, keys, 2, 0, &total, &distinct));
    assert(total == 0 && distinct == 0);
    assert(alive == 0);
}

static int hdfe_case(int trans, int accel, int slope)
{
    /* y, x, category, category, optional slope, two outputs, sample, marker. */
    nv = 8 + slope + request_group;
    for (int i = 0; i < nr; i++)
    {
        double a = i % 3, b = (i / 3) % 2;
        /* Balanced cells have residual +1/-1, orthogonal to both categories. */
        double resid = i < 12 ? 1 : -1;
        numbers[0][i] = 4 * a - 2 * b + resid;
        numbers[1][i] = 3 * a + b + 2 * resid;
        numbers[2][i] = a;
        numbers[3][i] = b;
        if (slope)
            numbers[4][i] = 1 + a;
        for (int j = 4 + slope; j < nv; j++)
            numbers[j][i] = j == nv - 1 ? 1 : 0;
    }
    char args[200];
    snprintf(args, sizeof(args), "2 2 0 0 0 10000 13 1 %d %d %d 0 %d tol=1e-12 verbose storeids",
             trans, accel, request_group, slope);
    int rc = chdfe_main(args);
    if (!rc)
    {
        for (int i = 0; i < nr; i++)
        {
            assert(numbers[6 + slope][i] == 1);
            if (!slope)
            {
                assert(fabs(numbers[4][i] - (i < 12 ? 1 : -1)) < 1e-10);
                assert(fabs(numbers[5][i] - (i < 12 ? 2 : -2)) < 1e-10);
            }
        }
    }
    return rc;
}

int main(void)
{
    mock.nobs = nobs;
    mock.nvars = nvars;
    mock.nvar = nvars;
    mock.isstr = isstr;
    mock.ismissing = missing;
    mock.missval = 8.98846567431158e307;
    mock.safestore = putnum;
    mock.macresave = save_macro;
    mock.spouterr = message;
    mock.spoutsml = message;
    int allocations;
    int (*distinct_cases[])(void) = {distinct_case, distinct_fallback_case};
    for (size_t j = 0; j < sizeof(distinct_cases) / sizeof(*distinct_cases); j++)
    {
        calls = 0;
        assert(distinct_cases[j]() == 0 && alive == 0);
        allocations = calls;
        for (int i = 1; i <= allocations; i++)
        {
            calls = 0;
            fail_at = i;
            assert(distinct_cases[j]() == 920 && alive == 0);
        }
        fail_at = 0;
    }
    distinct_keys();
    macro_error = 1;
    assert(distinct_case() == 459 && alive == 0);
    macro_error = 0;
    for (int trans = 0; trans <= 2; trans++)
        for (int accel = 0; accel <= 3; accel++)
        {
            if (trans == 1 && accel == 1)
                continue;
            assert(hdfe_case(trans, accel, 0) == 0 && alive == 0);
        }
    calls = 0;
    assert(hdfe_case(0, 1, 1) == 0 && alive == 0);
    allocations = calls;
    for (int i = 1; i <= allocations; i++)
    {
        calls = 0;
        fail_at = i;
        int rc = hdfe_case(0, 1, 1);
        assert((rc == 920 || rc == 0) && alive == 0);
    }
    fail_at = 0;
    request_group = 1;
    calls = 0;
    assert(hdfe_case(0, 1, 0) == 0 && alive == 0);
    allocations = calls;
    for (int i = 1; i <= allocations; i++)
    {
        calls = 0;
        fail_at = i;
        int rc = hdfe_case(0, 1, 0);
        assert((rc == 920 || rc == 0) && alive == 0);
    }
    fail_at = 0;
    request_group = 0;
    assert(cdistinct_main("nkeys=2147483647") == 198);
    write_error = 1;
    assert(hdfe_case(0, 1, 0) != 0 && alive == 0);
    write_error = 0;
    macro_error = 1;
    assert(hdfe_case(0, 1, 0) == 459 && alive == 0);
    macro_error = 0;
    assert(distinct_case() == 0 && hdfe_case(0, 1, 0) == 0 && alive == 0);
    puts("PASS counting and projection identities, solver choices, allocation failures, SPI "
         "failures, recovery");
    return 0;
}

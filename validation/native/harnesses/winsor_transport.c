/* Production grouping, quantiles, gather and stores; mocked host and loader.
 * The tracked allocator from memory_safety/alloc.h is prepended by the runner. */
#include "ctools_types.h"
#include "ctools_parse.c"
#include "ctools_order.c"
#include "cwinsor/cwinsor_impl.c"

static ST_plugin host;
ST_plugin *_stata_ = &host;
static int targets, byvars, store_at, store_calls, return_bounds, load_calls, read_at, read_calls;
static int host_rows = 256;

static ST_int host_nobs(void) { return host_rows; }

static ST_int host_nvars(void)
{
    return targets + byvars;
}
static ST_retcode save_scalar(char *name, double value)
{
    return 0;
}
static double input_value(int v, int obs);
static ST_retcode read_macro(char *name, char *out, ST_int length)
{
    assert(strcmp(name, "_cwinsor_bounds") == 0);
    snprintf(out, length, "bounds");
    return 0;
}
static ST_retcode store_bound(char *name, ST_int row, ST_int col, double value)
{
    assert(strcmp(name, "bounds") == 0);
    assert(value == input_value(row, col == 1 ? 1 : 255));
    return ++store_calls == store_at ? 459 : 0;
}
static ST_retcode read_bound(char *name, ST_int row, ST_int col, double *value)
{
    assert(strcmp(name, "bounds") == 0);
    *value = input_value(row, col == 1 ? 1 : 255);
    return ++read_calls == read_at ? 459 : 0;
}
int ctools_get_max_threads(void)
{
    return 4;
}
double ctools_timer_seconds(void)
{
    return 0;
}
void ctools_error(const char *a, const char *b, ...) {}
void ctools_error_alloc(const char *a) {}

void ctools_filtered_data_init(ctools_filtered_data *d)
{
    memset(d, 0, sizeof(*d));
}
void ctools_filtered_data_free(ctools_filtered_data *d)
{
    for (size_t j = 0; j < d->data.nvars; j++)
        free(d->data.vars[j].data.dbl);
    free(d->data.vars);
    free(d->data.sort_order);
    free(d->obs_map);
    memset(d, 0, sizeof(*d));
}

static double input_value(int v, int obs)
{
    return obs % 7 == 0 ? host.missval : obs + 1000.0 * v;
}

stata_retcode ctools_data_load(ctools_filtered_data *d, int *vars, size_t n, size_t start,
                               size_t end, int flags)
{
    load_calls++;
    const int count = 128;
    d->data.vars = calloc(n, sizeof(stata_variable));
    if (!d->data.vars)
        return STATA_ERR_MEMORY;
    d->data.nvars = n;
    d->data.nobs = count;
    d->obs_map = malloc(count * sizeof(perm_idx_t));
    if (!d->obs_map)
        goto failure;
    if (!(flags & CTOOLS_LOAD_NO_SORT_ORDER))
    {
        d->data.sort_order = malloc(count * sizeof(perm_idx_t));
        if (!d->data.sort_order)
            goto failure;
    }
    for (int i = 0; i < count; i++)
    {
        d->obs_map[i] = 2 * i + 1;
        if (d->data.sort_order)
            d->data.sort_order[i] = i;
    }
    for (size_t j = 0; j < n; j++)
    {
        d->data.vars[j].type = STATA_TYPE_DOUBLE;
        d->data.vars[j].nobs = count;
        double *x = d->data.vars[j].data.dbl = malloc(count * sizeof(double));
        if (!x)
            goto failure;
        for (int i = 0; i < count; i++)
            x[i] = vars[j] <= targets ? input_value(vars[j], 2 * i + 1) : i % 4;
    }
    return STATA_OK;
failure:
    /* The real loader frees partial allocations before returning an error. */
    ctools_filtered_data_free(d);
    return STATA_ERR_MEMORY;
}

stata_retcode ctools_sort_dispatch(stata_data *d, int *keys, size_t n, sort_algorithm_t alg)
{
    /* Replace the permutation allocation, as a sorting backend is allowed to. */
    perm_idx_t *order = malloc(d->nobs * sizeof(perm_idx_t));
    if (!order)
        return STATA_ERR_MEMORY;
    size_t next = 0;
    for (size_t group = 0; group < 4; group++)
        for (size_t i = group; i < d->nobs; i += 4)
            order[next++] = i;
    free(d->sort_order);
    d->sort_order = order;
    return STATA_OK;
}

stata_retcode ctools_store_filtered(double *values, size_t n, int var, perm_idx_t *map)
{
    for (size_t i = 0; i < n; i++)
    {
        assert(values[i] == input_value(var, map[i]));
        if (__atomic_add_fetch(&store_calls, 1, __ATOMIC_RELAXED) == store_at)
            return 459;
    }
    return STATA_OK;
}
stata_retcode ctools_store_filtered_rowpar(double *v, size_t n, int var, perm_idx_t *map)
{
    return ctools_store_filtered(v, n, var, map);
}

static int run_command(void)
{
    char args[128];
    snprintf(args, sizeof(args), "%d %d %s p=0 q=100 %s", targets, byvars,
             targets == 1 ? "1" : "1 2 3",
             return_bounds == 2 ? "deferbounds" : return_bounds ? "trim returnbounds" : "");
    store_calls = read_calls = load_calls = 0;
    int rc = cwinsor_main(args);
    if (!rc && return_bounds == 2) {
        assert(pending.nvars == (size_t)targets);
        assert(store_calls == 2*targets); /* Only the cutoff matrix was written. */
        rc = cwinsor_main("apply");
        assert(!pending.nvars);
        assert(load_calls == 1);
    }
    cwinsor_cleanup_state();
    leftovers();
    assert(alive == 0);
    return rc;
}

int main(void)
{
    host.nvars = host_nvars;
    host.nobs = host_nobs;
    host.scalsave = save_scalar;
    host.macuse = read_macro;
    host.matstore = store_bound;
    host.safematstore = store_bound;
    host.matel = read_bound;
    host.safematel = read_bound;
    host.missval = 8.98846567431158e307;
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    int checked = 0;
    for (return_bounds = 0; return_bounds <= 2; return_bounds++)
        for (targets = 1; targets <= 3; targets += 2)
            for (byvars = 0; byvars <= 2; byvars += 2)
            {
                if (return_bounds && byvars)
                    continue;
                calls = fail_at = store_at = 0;
                assert(run_command() == 0);
                int allocations = calls;
                for (int failure = 1; failure <= allocations; failure++)
                {
                    calls = 0;
                    fail_at = failure;
                    int rc = run_command();
                    assert(rc == 0 || rc == 920);
                    fail_at = 0;
                    assert(run_command() == 0);
                    checked++;
                }
                for (store_at = 1; store_at <= (return_bounds == 1 ? 2 * targets : 128); store_at *= 2)
                {
                    assert(run_command() == 459);
                    checked++;
                }
                store_at = 0;
                if (return_bounds == 2) {
                    for (read_at = 1; read_at <= 2*targets; read_at++) {
                        assert(run_command() == 459);
                        assert(store_calls == 2*targets); /* No data stores on a cutoff read failure. */
                        checked++;
                    }
                    read_at = 0;
                }
            }
    targets=1; byvars=0; fail_at=store_at=read_at=0;
    assert(cwinsor_main("apply") == 198);
    assert(cwinsor_main("1 0 1 p=0 q=100 deferbounds") == 0);
    host_rows++;
    assert(cwinsor_main("apply") == 198);
    host_rows--;
    assert(!pending.nvars);
    assert(cwinsor_main("1 0 1 p=0 q=100 deferbounds") == 0);
    assert(cwinsor_main("commit") == 0);
    assert(cwinsor_main("apply") == 198);
    leftovers(); assert(alive == 0);
    printf("PASS winsor transport: %d allocation/store injections, no leaks\n", checked);
    return 0;
}

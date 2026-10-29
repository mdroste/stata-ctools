#include "cpsmatch/cpsmatch_internal.h"
#include "cpsmatch/cpsmatch_kernels.c"
#include "cpsmatch/cpsmatch_options.c"
#include "cpsmatch/cpsmatch_stats.c"
#include "cpsmatch/cpsmatch_radius.c"
#include "cpsmatch/cpsmatch_impl.c"

static ST_plugin host;
ST_plugin *_stata_ = &host;
enum
{
    N = 192
};
static double treatment[N], score[N];
static perm_idx_t observations[N];
static stata_variable columns[2];
static int stores;
static ST_int first(void)
{
    return 1;
}
static ST_int last(void)
{
    return N;
}
static ST_boolean missing(double x)
{
    return x >= host.missval;
}
static ST_retcode store(ST_int v, ST_int row, double value)
{
    (void)v;
    (void)row;
    (void)value;
    stores++;
    return 0;
}
static ST_retcode scalar(char *name, double value)
{
    (void)name;
    (void)value;
    return 0;
}
static ST_retcode message(char *text)
{
    (void)text;
    return 0;
}
double ctools_timer_seconds(void)
{
    return 0;
}
int ctools_get_max_threads(void)
{
    return 4;
}
void ctools_filtered_data_init(ctools_filtered_data *d)
{
    memset(d, 0, sizeof(*d));
}
void ctools_filtered_data_free(ctools_filtered_data *d)
{
    (void)d;
}
stata_retcode ctools_data_load(ctools_filtered_data *d, int *v, size_t nv, size_t start, size_t end,
                               int flags)
{
    (void)v;
    (void)nv;
    (void)start;
    (void)end;
    (void)flags;
    columns[0].data.dbl = treatment;
    columns[1].data.dbl = score;
    d->data.nobs = N;
    d->data.nvars = 2;
    d->data.vars = columns;
    d->obs_map = observations;
    return STATA_OK;
}
int main(void)
{
    host.missval = 8.98846567431158e307;
    host.nobs1 = first;
    host.nobs2 = last;
    host.ismissing = missing;
    host.safestore = store;
    host.store = store;
    host.scalsave = scalar;
    host.spouterr = message;
    for (int i = 0; i < N; i++)
    {
        observations[i] = (perm_idx_t)(i + 1);
        treatment[i] = i % 2;
        score[i] = (double)(i % 17) / 17;
    }
    const char *args[] = {"1 2 0 3 4 5 192 0 1 0 0 0 .06 0 0 0",
                          "1 2 0 3 4 5 192 1 1 .1 0 0 .06 1 0 0"};
    for (int method = 0; method < 2; method++)
    {
        calls = fail_at = stores = 0;
        assert(cpsmatch_main(args[method]) == 0 && stores > 0 && alive == 0);
        int total = calls;
        for (int fail = 1; fail <= total; fail++)
        {
            calls = stores = 0;
            fail_at = fail;
            assert(cpsmatch_main(args[method]) == 920);
            assert(stores == 0 && alive == 0);
        }
        printf("path=%d allocation failures checked=%d\n", method, total);
    }
    puts("CPSMATCH_PATHS_OOM_COMPLETE");
}

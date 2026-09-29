/* Full command with mocked Stata loading/output; numerical code is production code. */
#include "alloc.h"
#include "ctools_types.h"
#include "ctools_parse.c"
#include "stplugin.h"
#include <math.h>
static int store_at, store_calls;
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static ST_int one(void) { return 1; }
static ST_int nvars(void) { return 6; }
static ST_int obs(void) { return 12; }
static ST_boolean selected(int i) { return 1; }
static ST_retcode save(char *k, double v) { return 0; }
static ST_retcode store(ST_int v, ST_int i, double x) {
    assert(x == (i - 1) % 10);
    return __atomic_add_fetch(&store_calls, 1, __ATOMIC_RELAXED) == store_at ? 459 : 0;
}
int ctools_get_max_threads(void) { return 8; }
double ctools_timer_seconds(void) { return 0; }
void ctools_error(const char *a, const char *b, ...) {}
void ctools_error_alloc(const char *a) {}
void ctools_msg(const char *a, const char *b, ...) {}
void ctools_filtered_data_init(ctools_filtered_data *d) { memset(d, 0, sizeof(*d)); }
void ctools_filtered_data_free(ctools_filtered_data *d) {
    for (size_t j = 0; j < d->data.nvars; j++)
        free(d->data.vars[j].data.dbl);
    free(d->data.vars);
    free(d->data.sort_order);
    free(d->obs_map);
    memset(d, 0, sizeof(*d));
}
stata_retcode ctools_data_load(ctools_filtered_data *d, int *v, size_t n, size_t s, size_t e,
                               int f) {
    int N = 512;
    d->data.vars = calloc(n, sizeof(stata_variable));
    if (!d->data.vars)
        return STATA_ERR_MEMORY;
    d->data.nvars = n;
    d->data.nobs = N;
    d->data.sort_order = malloc(N * sizeof(perm_idx_t));
    d->obs_map = malloc(N * sizeof(perm_idx_t));
    if (!d->obs_map || !d->data.sort_order)
        return STATA_ERR_MEMORY;
    for (int i = 0; i < N; i++) {
        d->obs_map[i] = i + 1;
        d->data.sort_order[i] = i;
    }
    for (size_t j = 0; j < n; j++) {
        stata_variable *q = d->data.vars + j;
        q->type = STATA_TYPE_DOUBLE;
        q->nobs = N;
        q->data.dbl = malloc(N * sizeof(double));
        if (!q->data.dbl)
            return STATA_ERR_MEMORY;
        for (int i = 0; i < N; i++)
            q->data.dbl[i] = j == 0 ? i / 64 : j == 1 ? i % 64 : i % 10;
    }
    return STATA_OK;
}
stata_retcode ctools_sort_dispatch(stata_data *d, int *v, size_t n, sort_algorithm_t a) {
    return STATA_OK;
}
stata_retcode ctools_apply_permutation(stata_data *d) { return STATA_OK; }

#include "crangestat/crangestat_impl.c"
#include "ctools_order.c"
int main(int argc, char **argv) {
    store_at = argc > 2 ? atoi(argv[2]) : 0;
    fail_at = argc > 1 ? atoi(argv[1]) : 0;
    mock.nobs1 = one;
    mock.nobs2 = obs;
    mock.nobs = obs;
    mock.nvars = nvars;
    mock.selobs = selected;
    mock.scalsave = save;
    mock.safestore = store;
    mock.store = store;
    mock.missval = 8.98846567431158e307;
    int rc = crangestat_main("2 2 1 2 3 4 5 6 1 3 4 0 1 low=0 high=0");
    int first_rc = rc, first_calls = calls;
    assert(alive == 0);
    fail_at = 0;
    store_at = 0;
    calls = 0;
    assert(crangestat_main("2 2 1 2 3 4 5 6 1 3 4 0 1 low=0 high=0") == 0);
    assert(alive == 0);
    rc = first_rc;
    calls = first_calls;
    printf("rc=%d calls=%d alive=%d\n", rc, calls, alive);
    leftovers();
    assert(alive == 0);
    return 0;
}

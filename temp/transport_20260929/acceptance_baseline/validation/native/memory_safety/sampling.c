/* Full command with mocked Stata loading/output; numerical code is production code. */
#include "alloc.h"
#include "ctools_types.h"
#include "stplugin.h"
#include <math.h>
static int store_at, store_calls;
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static ST_int one(void) { return 1; }
static ST_int obs(void) { return 12; }
static ST_boolean selected(int i) { return 1; }
static ST_retcode save(char *k, double v) { return 0; }
static ST_retcode store(ST_int v, ST_int i, double x) {
    return ++store_calls == store_at ? 459 : 0;
}
int ctools_get_max_threads(void) { return 8; }
double ctools_timer_seconds(void) { return 0; }
void ctools_error(const char *a, const char *b, ...) {}
void ctools_error_alloc(const char *a) {}
void ctools_msg(const char *a, const char *b, ...) {}
void ctools_filtered_data_init(ctools_filtered_data *d) { memset(d, 0, sizeof(*d)); }
void ctools_filtered_data_free(ctools_filtered_data *d) {}
stata_retcode ctools_data_load(ctools_filtered_data *d, int *v, size_t n, size_t s, size_t e,
                               int f) {
    abort();
}
stata_retcode ctools_sort_dispatch(stata_data *d, int *v, size_t n, sort_algorithm_t a) { abort(); }
stata_retcode ctools_apply_permutation(stata_data *d) { abort(); }
#ifdef BOOTSTRAP
#include "cbsample/cbsample_impl.c"
#else
#include "csample/csample_impl.c"
#endif
int main(int argc, char **argv) {
    store_at = argc > 2 ? atoi(argv[2]) : 0;
    fail_at = argc > 1 ? atoi(argv[1]) : 0;
    mock.nobs1 = one;
    mock.nobs2 = obs;
    mock.nobs = obs;
    mock.selobs = selected;
    mock.scalsave = save;
    mock.safestore = store;
    mock.store = store;
#ifdef BOOTSTRAP
    int rc = cbsample_main("1 0 0 n=6 seedhi=1 seedlo=2");
#else
    int rc = csample_main("1 0 count=6 seedhi=1 seedlo=2");
#endif
    assert(alive == 0);
    /* Even with eight configured workers, one group gets one RNG/workspace. */
#ifdef BOOTSTRAP
    assert(peak_bytes <= 700);
#else
    assert(peak_bytes <= 512);
#endif
    if (store_at)
        assert(rc == 459);
    printf("rc=%d allocations=%d alive=%d\n", rc, calls, alive);
    leftovers();
    assert(alive == 0);
    return 0;
}

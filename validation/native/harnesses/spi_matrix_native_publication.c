
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "ctools_spi.h"
static ST_plugin host;
ST_plugin *_stata_ = &host;
static int failed_alloc, write_calls, fail_at, host_K = 3;
static double got_b[3], got_V[9];
static void *map_malloc(size_t size) { return failed_alloc ? NULL : malloc(size); }
#define malloc map_malloc
#include "civreghdfe/civreghdfe_impl.c"
#undef malloc
#include "cqreg/cqreg_regress.c"
static ST_retcode message(char *s) { (void)s; return 0; }
static ST_retcode scalar_save(char *s, double v) { (void)s; (void)v; return 0; }
static ST_retcode unsafe_matrix(char *s, ST_int row, ST_int col, double v) {
    (void)s; (void)row; (void)col; (void)v;
    assert(!"unchecked matrix callback"); return 999;
}
static ST_retcode checked_matrix(char *s, ST_int row, ST_int col, double v) {
    write_calls++;
    if (write_calls == fail_at) return 459;
    int b = !strcmp(s, "__civreghdfe_b");
    if (row < 1 || row > (b ? 1 : host_K) || col < 1 || col > host_K) return 503;
    if (b) got_b[col - 1] = v;
    else got_V[(row - 1)*host_K + col - 1] = v;
    return 0;
}
int main(void) {
    host.matstore = unsafe_matrix;
    host.safematstore = checked_matrix;
    host.spouterr = message;
    host.scalsave = scalar_save;
    /* Internal [exog, omitted exog, endog] -> posted [endog, exog, omitted]. */
    double b[] = {2, 3}, V[] = {4, 1, 1, 9};
    int omitted[] = {0, 1, 0};
    assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 0);
    const double expected_b[] = {3, 2, 0}, expected_V[] = {9,1,0, 1,4,0, 0,0,0};
    assert(!memcmp(got_b, expected_b, sizeof got_b));
    assert(!memcmp(got_V, expected_V, sizeof got_V));
    for (int f = 1; f <= 12; f++) {
        write_calls = 0; fail_at = f;
        assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 459);
        assert(write_calls == f);
    }
    fail_at = 0; failed_alloc = 1; write_calls = 0;
    assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 920);
    assert(write_calls == 0);
    failed_alloc = 0; host_K = 2;
    assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 503);
    host_K = 3;
    assert(store_iv_matrices(b, V, 2, 1, 1, omitted) == 503);
    assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 0);
    cqreg_state state = {0};
    state.K = 2; state.N = 12; state.beta = b; state.V = V;
    for (int f = 1; f <= 4; f++) {
        write_calls = 0; fail_at = f;
        assert(store_results(&state) == 459);
        assert(write_calls == f);
    }
    fail_at = 0; host_K = 1;
    assert(store_results(&state) == 503);
    host_K = 2;
    assert(store_results(&state) == 0);
    assert(!memcmp(got_V, V, sizeof V));
    return 0;
}

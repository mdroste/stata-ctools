/* Exercise the production publication boundary with failed SPI writes. */
#include <assert.h>
#include <stdio.h>
#include "cqreg/cqreg_regress.c"

static ST_plugin host;
ST_plugin *_stata_ = &host;
static int scalar_calls, matrix_calls, fail_scalar, fail_matrix;
static ST_retcode scalar_save(char *name, double value)
{
    (void)name; (void)value;
    return ++scalar_calls == fail_scalar ? 459 : 0;
}
ST_retcode ctools_mat_store(const char *name, ST_int row, ST_int col, ST_double value)
{
    (void)name; (void)row; (void)col; (void)value;
    return ++matrix_calls == fail_matrix ? 460 : 0;
}
int main(void)
{
    double b[] = {2, 3}, V[] = {1, 0, 0, 1};
    cqreg_state state = {0};
    state.N = 120; state.K = 2; state.beta = b; state.V = V;
    state.vce_type = CQREG_VCE_CLUSTER;
    state.density_method = CQREG_DENSITY_KERNEL;
    host.scalsave = scalar_save;
    assert(store_results(&state) == 0);
    int nscalar = scalar_calls, nmatrix = matrix_calls;
    assert(nscalar == 13 && nmatrix == 4);
    for (int i = 1; i <= nscalar; i++) {
        scalar_calls = matrix_calls = 0; fail_scalar = i;
        assert(store_results(&state) == 459);
        assert(scalar_calls == i && matrix_calls == 0);
    }
    fail_scalar = 0;
    for (int i = 1; i <= nmatrix; i++) {
        scalar_calls = matrix_calls = 0; fail_matrix = i;
        assert(store_results(&state) == 460);
        assert(scalar_calls == nscalar && matrix_calls == i);
    }
    scalar_calls = matrix_calls = fail_matrix = 0;
    assert(store_results(&state) == 0);
    puts("PASS cqreg required scalar/matrix publication failures and recovery");
    return 0;
}

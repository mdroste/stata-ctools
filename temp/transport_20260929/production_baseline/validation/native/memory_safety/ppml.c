/* Full command with mocked Stata loading/output; numerical code is production code. */
#include "alloc.h"
#include "ctools_types.h"
#include "stplugin.h"
#include <math.h>
static int store_at, store_calls;
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static int clustered = 0, weighted = 0, mode = 0;
static ST_retcode message(char *s) { return 0; }
static ST_int one(void) { return 1; }
static ST_int obs(void) { return 12; }
static ST_int nvars(void) { return 4 + clustered + weighted; }
static ST_retcode scalar(char *k, double *v) {
    *v = 0;
    if (!strcmp(k, "__cpplmhdfe_K"))
        *v = 2;
    else if (!strcmp(k, "__cpplmhdfe_G"))
        *v = 1;
    else if (!strcmp(k, "__cpplmhdfe_maxiter"))
        *v = 10000;
    else if (!strcmp(k, "__cpplmhdfe_tolerance"))
        *v = 1e-10;
    else if (!strcmp(k, "__cpplmhdfe_standardize"))
        *v = 1;
    else if (!strcmp(k, "__cpplmhdfe_vce_type"))
        *v = clustered ? 2 : 0;
    else if (!strcmp(k, "__cpplmhdfe_cluster_terms"))
        *v = clustered;
    else if (!strcmp(k, "__cpplmhdfe_irls_maxiter"))
        *v = 100;
    else if (!strcmp(k, "__cpplmhdfe_irls_tol"))
        *v = 1e-8;
    else if (!strcmp(k, "__cpplmhdfe_has_weights"))
        *v = weighted;
    else if (!strcmp(k, "__cpplmhdfe_weight_type"))
        *v = 1;
    return 0;
}
static ST_retcode save(char *k, double v) { return 0; }
static ST_retcode store(ST_int v, ST_int i, double x) {
    return ++store_calls == store_at ? 459 : 0;
}
static ST_retcode mat(char *k, ST_int i, ST_int j, double x) { return 0; }
static ST_boolean missing(double x) { return x >= 8.98846567431158e307; }
int ctools_get_max_threads(void) { return 1; }
int ctools_get_max_threads_used(void) { return 1; }
int ctools_get_num_procs(void) { return 1; }
double ctools_timer_seconds(void) { return 0; }
void ctools_filtered_data_init(ctools_filtered_data *d) { memset(d, 0, sizeof(*d)); }
void stata_data_free(stata_data *d) {
    for (size_t j = 0; j < d->nvars; j++)
        free(d->vars[j].data.dbl);
    free(d->vars);
    memset(d, 0, sizeof(*d));
}
void ctools_filtered_data_free(ctools_filtered_data *d) {
    stata_data_free(&d->data);
    free(d->obs_map);
    d->obs_map = NULL;
}
stata_retcode ctools_data_load(ctools_filtered_data *d, int *v, size_t n, size_t s, size_t e,
                               int f) {
    d->data.vars = calloc(n, sizeof(stata_variable));
    if (!d->data.vars)
        return STATA_ERR_MEMORY;
    d->data.nvars = n;
    d->data.nobs = 12;
    d->obs_map = malloc(12 * sizeof(perm_idx_t));
    if (!d->obs_map)
        return STATA_ERR_MEMORY;
    for (int i = 0; i < 12; i++)
        d->obs_map[i] = i + 1;
    for (size_t j = 0; j < n; j++) {
        d->data.vars[j].nobs = 12;
        d->data.vars[j].type = STATA_TYPE_DOUBLE;
        double *p = malloc(12 * sizeof(double));
        d->data.vars[j].data.dbl = p;
        if (!p)
            return STATA_ERR_MEMORY;
        for (int i = 0; i < 12; i++)
            p[i] = j == 0                ? (mode == 1 ? 7 : 3 + 2 * i + sin(i))
                   : j == 1              ? i
                   : j == 2              ? i % 3
                   : j == 3 && clustered ? i % 4
                                         : 1;
    }
    return STATA_OK;
}
#include "cpplmhdfe/cpplmhdfe_irls.c"
#include "cpplmhdfe/cpplmhdfe_separation.c"
#include "creghdfe/creghdfe_hdfe.c"
#include "creghdfe/creghdfe_solver.c"
#include "creghdfe/creghdfe_utils.c"
#include "creghdfe/creghdfe_vce.c"
#include "ctools_hdfe_utils.c"
#include "ctools_matrix.c"
#include "ctools_ols.c"
#include "ctools_spi.c"
int main(int argc, char **argv) {
    store_at = argc > 5 ? atoi(argv[5]) : 0;
    fail_at = argc > 1 ? atoi(argv[1]) : 0;
    clustered = argc > 2 ? atoi(argv[2]) : 1;
    weighted = argc > 3 ? atoi(argv[3]) : 1;
    mode = argc > 4 ? atoi(argv[4]) : 0;
    mock.nobs1 = one;
    mock.nobs2 = obs;
    mock.nobs = obs;
    mock.nvars = nvars;
    mock.scalaruse = scalar;
    mock.scalsave = save;
    mock.spouterr = message;
    mock.spoutsml = message;
    mock.safestore = store;
    mock.matstore = mat;
    mock.safematstore = mat;
    mock.ismissing = missing;
    mock.missval = 8.98846567431158e307;
    int rc = do_ppml_regression(0, NULL);
    cleanup_ppml_state();
    int first_rc = rc, first_calls = calls;
    assert(alive == 0);
    fail_at = 0;
    store_at = 0;
    calls = 0;
    assert(do_ppml_regression(0, NULL) == 0);
    assert(alive == 0);
    rc = first_rc;
    calls = first_calls;
    printf("rc=%d allocations=%d alive=%d\n", rc, calls, alive);
    leftovers();
    assert(alive == 0);
    return 0;
}

/* Full IV command: injected allocations must fail cleanly or preserve estimates. */
#include "alloc.h"
#include "ctools_types.h"
#include "stplugin.h"
#include <math.h>

static ST_plugin host;
ST_plugin *_stata_ = &host;
enum
{
    OBS = 40
};
static int weighted, two_fe, singleton;
static double coefficient, variance, sample_size;
static ST_retcode message(char *s)
{
    (void)s;
    return 0;
}
static ST_int one(void)
{
    return 1;
}
static ST_int obs(void)
{
    return OBS;
}
static ST_retcode scalar(char *key, double *value)
{
    *value = 0;
    if (!strcmp(key, "__civreghdfe_K_endog") || !strcmp(key, "__civreghdfe_K_iv"))
        *value = 1;
    else if (!strcmp(key, "__civreghdfe_G"))
        *value = 1 + two_fe;
    else if (!strcmp(key, "__civreghdfe_has_weights"))
        *value = weighted;
    else if (!strcmp(key, "__civreghdfe_weight_type"))
        *value = 1;
    else if (!strcmp(key, "__civreghdfe_maxiter"))
        *value = 10000;
    else if (!strcmp(key, "__civreghdfe_tolerance"))
        *value = 1e-10;
    return 0;
}
static ST_retcode save(char *key, double value)
{
    if (!strcmp(key, "__civreghdfe_N"))
        sample_size = value;
    return 0;
}
static ST_retcode matrix(char *key, ST_int row, ST_int col, double value)
{
    (void)row;
    (void)col;
    if (!strcmp(key, "__civreghdfe_b"))
        coefficient = value;
    if (!strcmp(key, "__civreghdfe_V"))
        variance = value;
    return 0;
}
static ST_boolean missing(double x)
{
    return x >= 8.98846567431158e307;
}
int ctools_get_max_threads(void)
{
    return 1;
}
int ctools_get_max_threads_used(void)
{
    return 1;
}
int ctools_get_num_procs(void)
{
    return 1;
}
double ctools_timer_seconds(void)
{
    return 0;
}
void ctools_filtered_data_init(ctools_filtered_data *data)
{
    memset(data, 0, sizeof(*data));
}
void ctools_filtered_data_free(ctools_filtered_data *data)
{
    for (size_t j = 0; j < data->data.nvars; j++)
        free(data->data.vars[j].data.dbl);
    free(data->data.vars);
    free(data->obs_map);
    memset(data, 0, sizeof(*data));
}
stata_retcode ctools_data_load(ctools_filtered_data *data, int *indices, size_t n, size_t start,
                               size_t end, int flags)
{
    (void)indices;
    (void)start;
    (void)end;
    (void)flags;
    data->data.vars = calloc(n, sizeof(stata_variable));
    if (!data->data.vars)
        return STATA_ERR_MEMORY;
    data->data.nvars = n;
    data->data.nobs = OBS;
    data->obs_map = malloc(OBS * sizeof(perm_idx_t));
    if (!data->obs_map)
        return STATA_ERR_MEMORY;
    for (int i = 0; i < OBS; i++)
        data->obs_map[i] = i + 1;
    for (size_t j = 0; j < n; j++)
    {
        double *p = malloc(OBS * sizeof(double));
        data->data.vars[j].data.dbl = p;
        if (!p)
            return STATA_ERR_MEMORY;
        data->data.vars[j].nobs = OBS;
        data->data.vars[j].type = STATA_TYPE_DOUBLE;
        for (int i = 0; i < OBS; i++)
        {
            double z = sin(0.7 * i) + 0.1 * i;
            double x = z + cos(1.3 * i);
            p[i] = j == 0             ? 2 * x + 3 * (i % 4) + sin(2.1 * i)
                   : j == 1           ? x
                   : j == 2           ? z
                   : j == 3           ? (singleton && (i == 3 || i == 18) ? 100 + i : i % 4)
                   : j == 4 && two_fe ? (i / 3) % 3
                                      : 1 + (i % 3);
        }
    }
    return STATA_OK;
}

#include "civreghdfe/civreghdfe_impl.c"
#include "civreghdfe/civreghdfe_estimate.c"
#include "civreghdfe/civreghdfe_matrix.c"
#include "civreghdfe/civreghdfe_omega.c"
#include "civreghdfe/civreghdfe_tests.c"
#include "civreghdfe/civreghdfe_vce.c"
#include "hdfe/hdfe_solver.c"
#include "hdfe/hdfe_remap.c"
#include "hdfe/hdfe_sample.c"
#include "ctools_hdfe_utils.c"
#include "ctools_matrix.c"
#include "ctools_ols.c"
#include "ctools_spi.c"

int main(int argc, char **argv)
{
    int failure = argc > 1 ? atoi(argv[1]) : 0;
    weighted = argc > 2 ? atoi(argv[2]) : 0;
    two_fe = argc > 3 ? atoi(argv[3]) : 0;
    singleton = argc > 4 ? atoi(argv[4]) : 0;
    host.nobs1 = one;
    host.nobs2 = obs;
    host.nobs = obs;
    host.scalaruse = scalar;
    host.scalsave = save;
    host.spouterr = message;
    host.spoutsml = message;
    host.matstore = matrix;
    host.safematstore = matrix;
    host.ismissing = missing;
    host.missval = 8.98846567431158e307;
    assert(do_iv_regression() == 0);
    assert(alive == 0 && iv_state == NULL);
    double expected_b = coefficient, expected_v = variance, expected_n = sample_size;
    int allocations = calls;
    fail_at = failure;
    calls = 0;
    coefficient = variance = sample_size = NAN;
    int rc = do_iv_regression();
    printf("failure=%d rc=%d allocations=%d alive=%d b=%.17g expected=%.17g\n", failure, rc,
           allocations, alive, coefficient, expected_b);
    leftovers();
    fflush(stdout);
    assert(alive == 0 && iv_state == NULL);
    if (rc == 0)
    {
        assert(fabs(coefficient - expected_b) < 1e-10);
        assert(fabs(variance - expected_v) < 1e-10);
        assert(sample_size == expected_n);
    }
    else
        assert(rc == 920 || rc == 430);
    /* A failed command must not contaminate the next command in the same host. */
    fail_at = 0;
    assert(do_iv_regression() == 0);
    assert(alive == 0 && iv_state == NULL);
    assert(fabs(coefficient - expected_b) < 1e-10);
    return 0;
}

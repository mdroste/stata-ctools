#if defined(CTOOLS_CASE_csort_recovery)
#include "spi_mock.h"
#define ctools_sort_dispatch unused_production_dispatch
#include "ctools_types.c"
#undef ctools_sort_dispatch
stata_retcode ctools_sort_dispatch(stata_data *, int *, size_t, sort_algorithm_t);
#include "ctools_rowstage.h"
static int force_columns;
static stata_retcode recovery_load(ctools_rowstage *, const int *, size_t,
                                   const int *, size_t, size_t);
#define ctools_rowstage_load recovery_load
#include "csort/csort_impl.c"
#undef ctools_rowstage_load
#include "ctools_data_io.c"
#include "ctools_transport.c"
#include "ctools_data_plan.c"
#include "ctools_data_store.c"
#include "ctools_rowstage.c"

static double values[2][70000];
static atomic_int failed;
static stata_retcode recovery_load(ctools_rowstage *rs, const int *vars, size_t nv,
                                   const int *widths, size_t first, size_t nr)
{
    /* Exercise the column fallback after a staging allocation failure. */
    if (force_columns) return STATA_ERR_MEMORY;
    return ctools_rowstage_load(rs, vars, nv, widths, first, nr);
}
stata_retcode csort_stream_sort(int *keys, size_t nk, int *vars, size_t nv,
                                sort_algorithm_t alg, size_t block, int batch,
                                const int *widths, csort_stream_timings *timings)
{
    (void)keys; (void)nk; (void)vars; (void)nv; (void)alg;
    (void)block; (void)batch; (void)widths; (void)timings;
    abort(); /* The test explicitly selects nostream. */
}
static ST_retcode read_value(ST_int v, ST_int row, double *out)
{
    *out = values[v - 1][row - 1];
    return 0;
}
static ST_retcode write_value(ST_int v, ST_int row, double x)
{
    /* Fail once after earlier key writes. The identity rollback must repair
     * every row, including cells written concurrently by other workers. */
    if (v == 1 && row == 3 && !atomic_exchange(&failed, 1)) return 459;
    values[v - 1][row - 1] = x;
    return 0;
}
double ctools_timer_seconds(void) { return 0; }
stata_retcode ctools_sort_dispatch(stata_data *d, int *keys, size_t nk, sort_algorithm_t a)
{
    (void)keys; (void)nk; (void)a;
    /* Isolate publication from the independently tested sort engines. */
    for (size_t i = 0; i < d->nobs; i++) d->sort_order[i] = (perm_idx_t)(d->nobs - 1 - i);
    return STATA_OK;
}
int main(void)
{
    setup();
    mock.safevdata = mock.vdata = read_value;
    mock.safestore = mock.store = write_value;
    for (force_columns = 0; force_columns < 2; force_columns++) {
        for (int parallel = 0; parallel < 2; parallel++) {
            ctools_set_max_threads(parallel ? 8 : 1);
            nrows = parallel ? 70000 : 8;
            atomic_store(&failed, 0);
            for (int i = 0; i < nrows; i++) {
                values[0][i] = nrows - i;
                values[1][i] = 10 * (nrows - i);
            }
            ST_retcode rc = csort_main("1 nostream");
            assert(rc == ctools_stata_rc(STATA_ERR_STATA_WRITE) && atomic_load(&failed));
            for (int i = 0; i < nrows; i++) {
                assert(values[0][i] == nrows - i);
                assert(values[1][i] == 10 * (nrows - i));
            }
            ctools_destroy_global_pool();
        }
    }
    puts("PASS csort row-staged and column partial writes restore original rows");
}
#endif

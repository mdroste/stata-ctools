/* Native missing screening: exact stable compaction, sample publication,
 * default parallelism, all-valid allocation cost, and failure recovery. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdatomic.h>
#include "stplugin.h"
#include "ctools_config.h"

static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static int thread_limit = 4;
int ctools_get_max_threads(void) { return thread_limit; }
static double *marker_values;
static int store_error;
static atomic_uint writers;
static int store(ST_int var, ST_int row, ST_double value)
{
    assert(var == 7 && row >= 1 && value == 0);
    atomic_fetch_or(&writers, 1u << omp_get_thread_num());
    if (store_error) return store_error;
    marker_values[row] = value;
    return 0;
}

static int allocations, fail_at, live;
static size_t allocated_bytes;
static void *allocate(size_t n)
{
    if (++allocations == fail_at) return NULL;
    void *p = malloc(n);
    if (p) { live++; allocated_bytes += n; }
    return p;
}
static void *zero_allocate(size_t n, size_t size)
{
    void *p = allocate(n * size);
    if (p) memset(p, 0, n * size);
    return p;
}
static void release(void *p) { if (p) live--; free(p); }
#define malloc allocate
#define calloc zero_allocate
#define free release
#include "hdfe/hdfe_missing.c"
#undef malloc
#undef calloc
#undef free

static void run(size_t n, int pattern, int failure)
{
    ctools_filtered_data fd = {0};
    const size_t k = 6;
    fd.data.nobs = n; fd.data.nvars = k;
    fd.data.vars = calloc(k, sizeof(stata_variable));
    fd.obs_map = malloc((n + 1) * sizeof(perm_idx_t));
    marker_values = calloc(2*n + 4, sizeof(double));
    unsigned char *expected = malloc(n + 1);
    size_t kept = 0;
    for (size_t i = 0; i < n; i++) {
        fd.obs_map[i] = (perm_idx_t)(2*i + 1);
        marker_values[2*i+1] = 1;
        expected[i] = !(pattern == 3 || (pattern == 2 && i%3 == 0) ||
                        (pattern == 1 && (i%4096 == 0 || i%4096 == 4095 || i+1 == n)));
        kept += expected[i];
    }
    for (size_t j = 0; j < k; j++) {
        stata_variable *v = fd.data.vars + j;
        v->nobs = n; v->type = STATA_TYPE_DOUBLE;
        v->data.dbl = malloc((n + 1) * sizeof(double));
        for (size_t i = 0; i < n; i++) {
            v->data.dbl[i] = (double)(i*k + j);
            if (!expected[i] && i%k == j)
                v->data.dbl[i] = i%3 == 0 ? NAN : i%3 == 1 ? mock.missval : nextafter(mock.missval, INFINITY);
        }
        v->data.dbl[n] = -12345;
    }
    allocations = 0; allocated_bytes = 0; fail_at = failure;
    atomic_store(&writers, 0);
    int rc = hdfe_screen_missing(&fd, 7);
    assert(live == 0);
    if (failure || store_error) {
        assert(rc == (failure ? 920 : store_error));
        assert(fd.data.nobs == n);
        for (size_t i = 0; i < n; i++) assert(fd.obs_map[i] == 2*i+1);
        for (size_t i = 0; i < n; i++) assert(marker_values[2*i+1] == 1);
    } else {
        assert(rc == 0 && fd.data.nobs == kept);
        assert(fd.was_filtered == (kept < n));
        if (kept == n) assert(allocated_bytes <= (n ? 1+(n-1)/HDFE_SCREEN_BLOCK : 0));
        size_t out = 0;
        for (size_t i = 0; i < n; i++) {
            assert(marker_values[2*i] == 0 && marker_values[2*i+1] == expected[i]);
            if (expected[i]) {
                assert(fd.obs_map[out] == 2*i+1);
                for (size_t j = 0; j < k; j++) assert(fd.data.vars[j].data.dbl[out] == (double)(i*k+j));
                out++;
            }
        }
#ifdef _OPENMP
        if (thread_limit > 1 && n > 65536 && pattern == 2) assert((atomic_load(&writers) & ~1u) != 0);
#endif
        if (thread_limit == 1) assert((atomic_load(&writers) & ~1u) == 0);
    }
    for (size_t j = 0; j < k; j++) {
        assert(fd.data.vars[j].data.dbl[n] == -12345);
        free(fd.data.vars[j].data.dbl);
    }
    free(expected); free(marker_values); free(fd.obs_map); free(fd.data.vars);
}

static ST_int stream_first, stream_last, stream_n, stream_badtype;
static ST_int stream_error_row;
static double *stream_columns[6];
static ST_int seven(void) { return 7; }
static ST_int first_row(void) { return stream_first; }
static ST_int last_row(void) { return stream_last; }
static ST_int total_rows(void) { return stream_n; }
static ST_boolean string_column(ST_int j) { return j == stream_badtype; }
static ST_retcode stream_read(ST_int j, ST_int row, ST_double *v)
{
    assert(j >= 1 && j <= 7 && row >= stream_first && row <= stream_last);
    if (row == stream_error_row) return 459;
    if (j == 7) *v = marker_values[row];
    else {
        assert(row % 7 != 0); /* never read inputs outside the initial sample */
        *v = stream_columns[j-1][row];
    }
    return 0;
}

static void streaming(size_t n, int pattern)
{
    stream_first = 11; stream_last = (ST_int)n + 10; stream_n = (ST_int)n + 20;
    marker_values = malloc(((size_t)stream_n + 1) * sizeof(double));
    for (ST_int row = 0; row <= stream_n; row++) marker_values[row] = row % 7 != 0;
    for (int j = 0; j < 6; j++) {
        stream_columns[j] = malloc(((size_t)stream_n + 1) * sizeof(double));
        for (ST_int row = 0; row <= stream_n; row++) {
            double v = row + j;
            if (j == row % 6 && (pattern == 3 || (pattern == 1 && row % 4096 == 0) || (pattern == 2 && row % 3 == 0)))
                v = row % 2 ? NAN : nextafter(mock.missval, INFINITY);
            stream_columns[j][row] = v;
        }
    }
    allocations = 0; allocated_bytes = 0; fail_at = 1; /* must need no heap */
    atomic_store(&writers, 0);
    assert(hdfe_markout_main("") == 0);
    assert(allocations == 0 && live == 0);
    for (ST_int row = 0; row <= stream_n; row++) {
        int expected = row % 7 != 0;
        if (row >= stream_first && row <= stream_last)
            for (int j = 0; j < 6; j++) expected &= stream_columns[j][row] < mock.missval;
        assert(marker_values[row] == expected);
    }
#ifdef _OPENMP
    if (thread_limit > 1 && n > 65536 && pattern == 2) assert((atomic_load(&writers) & ~1u) != 0);
#endif
    if (n) {
        stream_error_row = stream_first;
        assert(hdfe_markout_main("") == 459);
        stream_error_row = 0;
        stream_badtype = 3;
        assert(hdfe_markout_main("") == 109);
        stream_badtype = 0;
        stream_first = 0;
        assert(hdfe_markout_main("") == 198);
        stream_first = 11;
        if (pattern) {
            for (ST_int row = stream_first; row <= stream_last; row++) marker_values[row] = row % 7 != 0;
            store_error = 459;
            if (n > 4096 || pattern == 3) assert(hdfe_markout_main("") == 459);
            store_error = 0;
        }
    }
    for (int j = 0; j < 6; j++) free(stream_columns[j]);
    free(marker_values);
    fail_at = 0;
}

int main(void)
{
    mock.missval = 8.98846567431158e307;
    mock.safestore = store;
    mock.nvars = seven;
    mock.nobs1 = first_row;
    mock.nobs2 = last_row;
    mock.nobs = total_rows;
    mock.isstr = string_column;
    mock.vdata = stream_read;
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    const size_t sizes[] = {0,1,4095,4096,4097,32767,32768,65539,262151};
    for (int t = 1; t <= 4; t *= 4) {
        thread_limit = t;
        for (size_t n = 0; n < sizeof(sizes)/sizeof(*sizes); n++)
            for (int pattern = 0; pattern < 4; pattern++) run(sizes[n], pattern, 0);
        run(65539, 1, 1);
        run(65539, 1, 2);
        for (int error = 1; error <= 2; error++) {
            store_error = error == 1 ? 1 : 459;
            run(65539, 2, 0);
        }
        store_error = 0;
        for (size_t n = 0; n < sizeof(sizes)/sizeof(*sizes); n++)
            for (int pattern = 0; pattern < 4; pattern++) streaming(sizes[n], pattern);
    }
    assert(hdfe_screen_missing(NULL, 7) == 198);
    puts("PASS native HDFE missing screening");
}

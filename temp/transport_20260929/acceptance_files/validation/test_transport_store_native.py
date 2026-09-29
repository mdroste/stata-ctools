"""Numeric filtered-store safety against an independent SPI model.

Exercises ordered, reversed, and repeated destination maps above and below the
parallel threshold; repeated writes must preserve input-order last-write-wins.
Set CTOOLS_TEST_SOURCE to test an isolated source snapshot. Runs without Stata.
"""
from pathlib import Path
import os
import subprocess
import tempfile
from test_p2_native import SPI, ROOT

SOURCE = SPI + r'''
#include <stdint.h>
#include "ctools_config.h"
#undef MIN_OBS_PER_THREAD
#define MIN_OBS_PER_THREAD 256
static int fail_bitmap;
static atomic_int bitmap_attempts;
static void *test_map_calloc(size_t n, size_t bytes) {
    if (bytes == 1) atomic_fetch_add(&bitmap_attempts, 1);
    if (fail_bitmap && bytes == 1) return NULL;
    return calloc(n, bytes);
}
#define calloc test_map_calloc
#include "ctools_data_io.c"
#undef calloc

#define ROWS 6007
#define SPARSE_ROWS 100000
static double stored[2][SPARSE_ROWS + 1];
static pthread_t caller;
static int require_caller;
static atomic_int wrong_thread, worker_used;
static ST_retcode numeric_write(ST_int v, ST_int row, ST_double value) {
    assert(v >= 1 && v <= 2 && row >= 1 && row <= nrows);
    if (require_caller && !pthread_equal(caller, pthread_self())) {
        atomic_store(&wrong_thread, 1);
        return 459; /* Do not itself execute a racy write in a failing build. */
    }
    if (!pthread_equal(caller, pthread_self())) atomic_store(&worker_used, 1);
    atomic_fetch_add(&writes, 1);
    if (row == bad_row) return 459;
    stored[v - 1][row] = value;
    return 0;
}
static ST_retcode unchecked_write(ST_int v, ST_int row, ST_double value) {
    (void)v; (void)row; (void)value;
    assert(!"numeric stores must use the checked SPI callback");
    return 459;
}

static void destinations(size_t n) {
    double *values = malloc(n * sizeof(*values));
    perm_idx_t *map = malloc(n * sizeof(*map));
    double expected[ROWS + 1];
    assert(values && map);
    for (size_t i = 0; i < n; i++) {
        /* Include signed zero and a missing value, and check exact bits. */
        values[i] = i % 23 == 0 ? mock.missval :
                    i % 17 == 0 ? -0.0 : (double)i / 8.0 + 0.125;
    }
    for (int pattern = 0; pattern < 6; pattern++) {
        memset(stored, 0, sizeof(stored));
        memset(expected, 0, sizeof(expected));
        for (size_t i = 0; i < n; i++) {
            map[i] = pattern == 0 ? i + 2 : /* increasing */
                     pattern == 1 ? n - i : /* unique reverse */
                     pattern == 2 ? i / 2 + 2 : /* adjacent duplicates */
                     pattern == 3 ? i % 17 + 2 : /* distant duplicates */
                     pattern == 4 ? 2 : /* every write same destination */
                     i % 2 ? i + 1 : i + 1 < n ? i + 3 : i + 2;
                     /* pattern 5 swaps adjacent destinations: unique,
                      * neither increasing nor decreasing. */
            expected[map[i]] = values[i];
        }
        require_caller = (pattern >= 2 && pattern <= 4) ||
                         (pattern == 5 && fail_bitmap && n > 2);
        writes = wrong_thread = worker_used = 0;
        assert(ctools_store_filtered_rowpar(values,n,2,map) == STATA_OK);
        assert(!wrong_thread && (size_t)writes == n);
        #ifdef _OPENMP
        if (!require_caller && n >= 512 && ctools_get_max_threads() > 1)
            assert(worker_used); /* Unordered unique maps retain parallelism. */
        #endif
        assert(memcmp(stored[1],expected,sizeof(expected)) == 0);
        /* Other columns must never be redirected to by the store callback. */
        for (int row = 1; row <= ROWS; row++) assert(stored[0][row] == 0.0);
        bad_row = (int)map[n / 2];
        assert(ctools_store_filtered_rowpar(values,n,2,map) == STATA_ERR_STATA_WRITE);
        bad_row = 0;
        map[n-1] = ROWS + 1; writes = wrong_thread = 0;
        assert(ctools_store_filtered_rowpar(values,n,2,map) == STATA_ERR_INVALID_INPUT);
        assert(!wrong_thread && !writes); /* Entire map validated before writes. */
        map[n-1] = 0; writes = 0;
        assert(ctools_store_filtered_rowpar(values,n,2,map) == STATA_ERR_INVALID_INPUT && !writes);
    }
    free(map); free(values);
}
/* An unordered sparse range can require more bitmap bytes than the input
 * map itself. Such a map must remain sequential, without attempting scratch
 * allocation; check both unique destinations and distant repeated writes. */
static void sparse_destinations(void) {
    const size_t n = 513;
    double values[513];
    perm_idx_t map[513];
    double *expected = calloc(SPARSE_ROWS + 1, sizeof(*expected));
    assert(expected);
    nrows = SPARSE_ROWS;
    fail_bitmap = 0;
    require_caller = 1;
    for (int duplicates = 0; duplicates <= 1; duplicates++) {
        memset(stored, 0, sizeof(stored));
        memset(expected, 0, (SPARSE_ROWS + 1) * sizeof(*expected));
        for (size_t i = 0; i < n; i++) {
            values[i] = i % 17 == 0 ? -0.0 : (double)i + 0.125;
            if (duplicates) {
                map[i] = i % 11 == 0 || i == n - 1 ? SPARSE_ROWS : i % 7 + 2;
            } else {
                /* Unique but neither increasing nor decreasing. */
                map[i] = i == n - 1 ? SPARSE_ROWS : i % 2 ? i + 1 : i + 3;
            }
            expected[map[i]] = values[i];
        }
        writes = wrong_thread = worker_used = bitmap_attempts = 0;
        assert(ctools_store_filtered_rowpar(values,n,2,map) == STATA_OK);
        assert((size_t)writes == n && !wrong_thread && !worker_used);
        assert(bitmap_attempts == 0);
        assert(memcmp(stored[1],expected,(SPARSE_ROWS + 1)*sizeof(*expected)) == 0);
        if (duplicates) assert(stored[1][SPARSE_ROWS] == values[n - 1]);
        map[n - 1] = SPARSE_ROWS + 1; writes = 0;
        assert(ctools_store_filtered_rowpar(values,n,2,map) == STATA_ERR_INVALID_INPUT && !writes);
    }
    free(expected);
    nrows = ROWS;
    require_caller = 0;
}

int main(void) {
    setup(); nrows = ROWS; string_var = 0; caller = pthread_self();
    mock.safestore = numeric_write; mock.store = unchecked_write;
    assert(ctools_store_filtered_rowpar(NULL,0,1,NULL) == STATA_OK);
    double value = 1; perm_idx_t map = 1;
    assert(ctools_store_filtered_rowpar(NULL,1,1,&map) == STATA_ERR_INVALID_INPUT);
    assert(ctools_store_filtered_rowpar(&value,1,1,NULL) == STATA_ERR_INVALID_INPUT);
    assert(ctools_store_filtered_rowpar(&value,1,3,&map) == STATA_ERR_INVALID_INPUT);
    int thread_counts[] = {1,4,12};
    size_t sizes[] = {1,31,511,512,513,5003};
    for (size_t t = 0; t < sizeof(thread_counts)/sizeof(*thread_counts); t++) {
        ctools_set_max_threads(thread_counts[t]);
        for (fail_bitmap = 0; fail_bitmap <= 1; fail_bitmap++)
            for (size_t n = 0; n < sizeof(sizes)/sizeof(*sizes); n++) destinations(sizes[n]);
        sparse_destinations();
        ctools_destroy_global_pool();
    }
    return 0;
}
'''


def main():
    source_dir = Path(os.environ.get('CTOOLS_TEST_SOURCE', ROOT / 'src'))
    with tempfile.TemporaryDirectory(prefix='ctools-store-native-') as tmp:
        directory = Path(tmp)
        source = directory / 'test.c'
        source.write_text(SOURCE)
        for omp in (False, True):
            binary = directory / ('test-' + str(omp))
            flags = ['-std=c11', '-O2', '-DSD_FASTMODE', '-fno-fast-math',
                     '-fsanitize=' + os.environ.get('CTOOLS_SANITIZERS', 'undefined'),
                     '-fno-omit-frame-pointer', '-pthread', '-I' + str(source_dir)]
            if os.uname().sysname == 'Darwin':
                flags += ['-DSYSTEM=APPLEMAC', '-D_DARWIN_C_SOURCE',
                          '-Wl,-dead_strip', '-Wl,-undefined,dynamic_lookup']
                if omp:
                    prefix = os.environ.get('LIBOMP_PREFIX', '/opt/homebrew/opt/libomp')
                    flags += ['-Xpreprocessor', '-fopenmp', '-I' + prefix + '/include',
                              prefix + '/lib/libomp.a']
            else:
                flags += ['-DSYSTEM=OPUNIX', '-D_GNU_SOURCE', '-ffunction-sections',
                          '-fdata-sections', '-Wl,--gc-sections', '-lm']
                if omp:
                    flags += ['-fopenmp']
            subprocess.run([os.environ.get('CC', 'clang'), *flags, str(source),
                            *[str(source_dir / f) for f in
                              ['ctools_types.c', 'ctools_arena.c', 'ctools_threads.c']],
                            '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True, timeout=120,
                           env=dict(os.environ, UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1'))
            print(f'PASS numeric-store openmp={omp} source={source_dir}', flush=True)


if __name__ == '__main__':
    main()

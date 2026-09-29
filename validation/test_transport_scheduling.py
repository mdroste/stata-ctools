"""Transfer scheduling regressions against a deterministic, independent SPI model.

Optional CTOOLS_TEST_SOURCE selects an experimental source directory. Test both
serial and OpenMP builds; every stored cell, selective variable mapping, sorted
gather, duplicate destination, range check, and SPI failure is checked.
"""
from pathlib import Path
import os
import subprocess
import tempfile
from test_p2_native import SPI, ROOT

SOURCE = SPI + r'''
#include <stdio.h>
#include "ctools_config.h"
#undef MIN_OBS_PER_THREAD
#define MIN_OBS_PER_THREAD 256
#define IO_CHUNK_MIN_NUMERIC_ROWS 512
static int fail_flat;
static void *test_flat_calloc(size_t n, size_t size) {
    if (fail_flat && size == 65 && n > 512) return NULL;
    return calloc(n,size);
}
#define calloc test_flat_calloc
#include "ctools_data_io.c"
#undef calloc

#define COLS 20
#define ROWS 6007
#define WIDTH 64
static char stored[COLS][ROWS + 1][WIDTH + 1];
static double numeric[COLS][ROWS + 1];
static pthread_t caller;
static int serial_required;
static ST_int many_vars(void) { return COLS; }
static ST_boolean string_column(ST_int v) { return v % 2 == 0; }
static ST_boolean never_strl(ST_int v) { (void)v; return 0; }
static void text_value(int v, int row, char *out) {
    if (row % 17 == 0) { out[0] = 0; return; }
    memset(out, 'a' + row % 26, WIDTH);
    char prefix[32];
    int len = snprintf(prefix, sizeof(prefix), "%02d_%06d_", v, row);
    memcpy(out, prefix, (size_t)len);
    out[WIDTH - 2] = (char)0xc3; out[WIDTH - 1] = (char)0xa9;
    out[WIDTH] = 0;
}
static double number_value(int v, int row) {
    return row % 23 == 0 ? mock.missval : row + v * 0.125;
}
static ST_retcode read_text(ST_int v, ST_int row, char *out) {
    assert(v >= 1 && v <= COLS && string_column(v));
    assert(row >= 1 && row <= nrows);
    if (row == bad_row) return 459;
    text_value(v, row, out); return 0;
}
static ST_retcode read_number(ST_int v, ST_int row, double *out) {
    assert(v >= 1 && v <= COLS && !string_column(v));
    assert(row >= 1 && row <= nrows);
    if (row == bad_row) return 459;
    *out = number_value(v, row); return 0;
}
static ST_retcode write_text(ST_int v, ST_int row, char *value) {
    assert(v >= 1 && v <= COLS && string_column(v));
    assert(row >= 1 && row <= nrows && strlen(value) <= WIDTH);
    if (serial_required) assert(pthread_equal(pthread_self(), caller));
    atomic_fetch_add(&writes, 1);
    if (row == bad_row) return 459;
    strcpy(stored[v - 1][row], value); return 0;
}
static ST_retcode write_number(ST_int v, ST_int row, double value) {
    assert(v >= 1 && v <= COLS && !string_column(v));
    assert(row >= 1 && row <= nrows);
    atomic_fetch_add(&writes, 1);
    if (row == bad_row) return 459;
    numeric[v - 1][row] = value; return 0;
}

/* All-numeric multi-column loads go through the chunked scheduler once
 * threads and rows allow several chunks per column; verify both identity
 * and filtered chunk maps, plus the shared error injection. */
static void numeric_chunks(void) {
    int indices[10];
    for (int j = 0; j < 10; j++) indices[j] = 2*j + 1;
    for (filter = 0; filter <= 1; filter++) {
        ctools_filtered_data fd;
        assert(ctools_data_load_ex(&fd, indices, 10, 1, nrows,
               filter ? CTOOLS_LOAD_CHECK_IF : CTOOLS_LOAD_SKIP_IF, NULL) == STATA_OK);
        size_t expected = 0;
        for (int row = 1; row <= nrows; row++) if (!filter || selected(row)) {
            for (int j = 0; j < 10; j++)
                assert(fd.data.vars[j].data.dbl[expected] == number_value(indices[j], row));
            expected++;
        }
        assert(expected == fd.data.nobs);
        ctools_filtered_data_free(&fd);
        bad_row = 777;
        assert(ctools_data_load_ex(&fd, indices, 10, 1, nrows,
               filter ? CTOOLS_LOAD_CHECK_IF : CTOOLS_LOAD_SKIP_IF, NULL) == STATA_ERR_STATA_READ);
        ctools_filtered_data_free(&fd);
        bad_row = 0;
    }
    filter = 0;
}

static void transfers(void) {
    int widths[COLS], indices[COLS];
    for (int j = 0; j < COLS; j++) widths[j] = string_column(j+1) ? WIDTH : 0;
    for (filter = 0; filter <= 1; filter++)
    for (int reverse = 0; reverse <= 1; reverse++)
    for (int hints = 0; hints <= 1; hints++) {
        for (int j = 0; j < COLS; j++) indices[j] = reverse ? COLS-j : j+1;
        ctools_filtered_data fd;
        assert(ctools_data_load_ex(&fd, indices, COLS, 3, nrows-2,
               CTOOLS_LOAD_CHECK_IF, hints ? widths : NULL) == STATA_OK);
        size_t expected = 0;
        for (int row = 3; row <= nrows-2; row++) if (selected(row)) {
            assert(fd.obs_map[expected] == (perm_idx_t)row);
            for (int j = 0; j < COLS; j++) {
                int v = indices[j];
                if (string_column(v)) {
                    char value[WIDTH+1]; text_value(v,row,value);
                    assert(!strcmp(fd.data.vars[j].data.str[expected],value));
                } else assert(fd.data.vars[j].data.dbl[expected] == number_value(v,row));
            }
            expected++;
        }
        assert(expected == fd.data.nobs);
        writes = 0;
        assert(ctools_data_store_selective(&fd.data, indices, COLS, 2) == STATA_OK);
        assert((size_t)writes == expected * COLS);
        for (int sorted = 0; sorted <= !reverse; sorted++) {
            if (sorted) {
                for (size_t i = 0; i < expected; i++) fd.data.sort_order[i] = expected-i-1;
                assert(ctools_data_store_sorted(&fd.data,2) == STATA_OK);
            }
            for (size_t i = 0; i < expected; i++) {
                int source = (int)fd.obs_map[sorted ? expected-i-1 : i];
                for (int j = 0; j < COLS; j++) {
                    int v = indices[j];
                    if (string_column(v)) {
                        char value[WIDTH+1]; text_value(v,source,value);
                        assert(!strcmp(stored[v-1][i+2],value));
                    } else assert(numeric[v-1][i+2] == number_value(v,source));
                }
            }
        }
        bad_row = 777;
        assert(ctools_data_store_selective(&fd.data,indices,COLS,2) == STATA_ERR_STATA_WRITE);
        bad_row = 0;
        ctools_filtered_data_free(&fd);
        bad_row = 777;
        assert(ctools_data_load_ex(&fd,indices,COLS,3,nrows-2,
            CTOOLS_LOAD_CHECK_IF,hints ? widths : NULL) == STATA_ERR_STATA_READ);
        ctools_filtered_data_free(&fd); bad_row = 0;
    }
    filter = 0;
    widths[1] = WIDTH-1;
    assert(ctools_data_load_ex(&(ctools_filtered_data){0},NULL,0,1,nrows,
        CTOOLS_LOAD_SKIP_IF,widths) == STATA_ERR_STATA_READ);
}

static void destinations(void) {
    const size_t n = 5003;
    char **values = malloc(n * sizeof(*values));
    perm_idx_t *map = malloc(n * sizeof(*map));
    char (*expected)[WIDTH+1] = calloc(ROWS+1,WIDTH+1);
    assert(values && map && expected);
    for (size_t i = 0; i < n; i++) values[i] = i % 7 ? "later" : "earlier";
    values[11] = NULL;
    for (int pattern = 0; pattern < 3; pattern++) {
        memset(stored[1],0,sizeof(stored[1]));
        memset(expected,0,(ROWS+1)*(WIDTH+1));
        for (size_t i = 0; i < n; i++) {
            map[i] = pattern == 0 ? i+2 : pattern == 1 ? n-i : i/2+2;
            strcpy(expected[map[i]],values[i] ? values[i] : "");
        }
        serial_required = pattern != 0;
        writes = 0;
        assert(ctools_store_filtered_str(values,n,2,map) == STATA_OK);
        assert((size_t)writes == n);
        for (size_t row = 1; row <= ROWS; row++) assert(!strcmp(expected[row],stored[1][row]));
        bad_row = 777;
        assert(ctools_store_filtered_str(values,n,2,map) == STATA_ERR_STATA_WRITE);
        bad_row = 0;
        map[n-1] = ROWS+1; writes = 0;
        assert(ctools_store_filtered_str(values,n,2,map) == STATA_ERR_INVALID_INPUT && writes == 0);
    }
    serial_required = 0;
    free(expected); free(map); free(values);
}

int main(void) {
    setup(); nrows = ROWS; caller = pthread_self();
    mock.nvars = mock.nvar = many_vars;
    mock.isstr = string_column; mock.isstrl = never_strl;
    mock.safevdata = read_number; mock.safestore = write_number;
    mock.vdata = read_number; mock.store = write_number;
    mock.sdata = read_text; mock.sstore = write_text;
    int thread_cases[] = {1, 4, 12};
    for (int tc = 0; tc < 3; tc++) { int threads = thread_cases[tc];
        ctools_set_max_threads(threads); transfers(); destinations(); numeric_chunks();
        fail_flat = 1; transfers(); fail_flat = 0;
        ctools_destroy_global_pool();
    }
    return 0;
}
'''


def main():
    source_dir = Path(os.environ.get('CTOOLS_TEST_SOURCE', ROOT / 'src'))
    with tempfile.TemporaryDirectory(prefix='ctools-scheduling-') as tmp:
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
            print(f'PASS scheduling openmp={omp} source={source_dir}', flush=True)


if __name__ == '__main__':
    main()

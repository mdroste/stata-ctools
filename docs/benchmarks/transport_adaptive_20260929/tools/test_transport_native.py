"""Transport value/ownership regressions, predicate counts, and scratch OOM.

Runs serial and OpenMP variants under UBSan (set CTOOLS_SANITIZERS=address,undefined
to add ASan on supported systems). Existing SPI fault tests remain
separate and unchanged. No Stata executable is needed for these mocked-SPI tests.
"""
from pathlib import Path
import os
import subprocess
import tempfile
from test_p2_native import SPI, ROOT

SOURCE = SPI + r'''
#include "ctools_config.h"
#undef MIN_OBS_PER_THREAD
#define MIN_OBS_PER_THREAD 256
static int fail_scratch;
static void *test_calloc(size_t n, size_t size) {
    if (fail_scratch && size == sizeof(uint64_t)) return NULL;
    return calloc(n, size);
}
#define calloc test_calloc
#include "ctools_data_io.c"
#undef calloc

static int pattern, cutoff, select_calls, content_width;
static atomic_int type_calls;
static pthread_t calling_thread;
static int keep(int row) {
    switch (pattern) {
        case 0: return 1;
        case 1: return 0;
        case 2: return row != cutoff;
        case 3: return row == cutoff;
        case 4: return row % 3 != 0;
        default: return row > cutoff;
    }
}
static ST_boolean selection(ST_int row) { select_calls++; return keep(row); }
static ST_boolean string_type(ST_int v) { atomic_fetch_add(&type_calls,1); return v == strl; }
static int length_at(int row) { return row % 7 == 0 ? 0 : content_width; }
static ST_retcode varied_string(ST_int v, ST_int row, char *out) {
    assert(v == 2 && row >= 1 && row <= nrows);
    if (row == bad_row) return 459;
    int len = length_at(row);
    memset(out, 'A' + row % 26, len);
    if (len >= 2) { out[0] = (char)0xc3; out[1] = (char)0xa9; }
    out[len] = 0;
    return 0;
}
static ST_int varied_length(ST_int v, ST_int row) { (void)v; return length_at(row); }
static ST_retcode varied_strl(ST_int v, ST_int row, char *out, ST_int capacity) {
    assert(pthread_equal(pthread_self(), calling_thread));
    assert(capacity == length_at(row) + 1);
    return varied_string(v, row, out) ? -1 : length_at(row);
}

static void maps(void) {
    const int sizes[] = {1, 2, 31, 32, 33, 63, 64, 65, 127, 128, 129, 6000};
    for (size_t k = 0; k < sizeof(sizes)/sizeof(*sizes); k++) {
        nrows = sizes[k] + 10;
        for (int start = 1; start <= 3; start += 2)
        for (pattern = 0; pattern <= 5; pattern++)
        for (int where = 0; where < 3; where++)
        for (fail_scratch = 0; fail_scratch <= 1; fail_scratch++) {
            int end = start + sizes[k] - 1;
            cutoff = where == 0 ? start : (where == 1 ? (start+end)/2 : end);
            perm_idx_t *map = NULL;
            size_t count, range, expected = 0;
            int filtered;
            select_calls = 0;
            assert(build_obs_map(start, end, CTOOLS_LOAD_CHECK_IF,
                                 &map, &count, &range, &filtered) == 0);
            for (int row = start; row <= end; row++) if (keep(row)) {
                assert(map && map[expected] == (perm_idx_t)row); expected++;
            }
            assert(count == expected && range == (size_t)sizes[k]);
            assert(filtered == (expected != range));
            if (!fail_scratch && sizes[k] > 32) assert(select_calls == sizes[k]);
            ctools_aligned_free(map);
        }
    }
    fail_scratch = 0;
}

static void values(void) {
    nrows = 6000; string_var = 2;
    mock.sdata = varied_string; mock.isstrl = string_type;
    mock.sdatalen = varied_length; mock.strldata = varied_strl;
    const int widths[] = {0, 1, 7, 31, 244, 2044, 2045};
    int indices[] = {2, 1};
    for (size_t w = 0; w < sizeof(widths)/sizeof(*widths); w++) {
        content_width = widths[w];
        for (strl = 0; strl <= 2; strl += 2)
        for (int hints = 0; hints <= 1; hints++)
        for (size_t columns = 1; columns <= 2; columns++)
        for (pattern = 0; pattern <= 4; pattern += 4) {
            ctools_filtered_data fd;
            int widths_by_spi_index[] = {0, content_width ? content_width : 1};
            type_calls = 0;
            assert(ctools_data_load_ex(&fd, indices, columns, 3, nrows-2,
                CTOOLS_LOAD_CHECK_IF, hints ? widths_by_spi_index : NULL) == 0);
            size_t expected = 0;
            for (int row = 3; row <= nrows-2; row++) if (keep(row)) {
                char value[2046];
                assert(varied_string(2,row,value)==0);
                assert(fd.obs_map[expected] == (perm_idx_t)row);
                assert(!strcmp(fd.data.vars[0].data.str[expected], value));
                if (columns==2) assert(fd.data.vars[1].data.dbl[expected] == row);
                expected++;
            }
            assert(fd.data.nobs == expected);
            assert(type_calls == 1); /* Constant column metadata, not per cell. */
            ctools_filtered_data_free(&fd);
        }
    }
    strl = 0; content_width = 31; pattern = 0;
    int stale[] = {0, 7};
    ctools_filtered_data fd;
    for (size_t columns = 1; columns <= 2; columns++) {
        assert(ctools_data_load_ex(&fd, indices, columns, 1, nrows, 0, stale) == STATA_ERR_STATA_READ);
        ctools_filtered_data_free(&fd);
        bad_row = 333;
        assert(ctools_data_load_ex(&fd, indices, columns, 1, nrows, 0, NULL) == STATA_ERR_STATA_READ);
        ctools_filtered_data_free(&fd); bad_row = 0;
    }
    /* Exhaustion and allocation fallback preserve ownership of private arenas. */
    ctools_string_arena *arena = ctools_string_arena_create(8, CTOOLS_STRING_ARENA_STRDUP_FALLBACK);
    assert(arena);
    char *a = copy_private_string(arena, "1234567");
    char *b = copy_private_string(arena, "fallback");
    char *c = copy_private_string(NULL, "no arena");
    assert(a == arena->base && arena->used == 8 && arena->has_fallback);
    assert(!strcmp(b,"fallback") && !ctools_string_arena_owns(arena,b));
    assert(!strcmp(c,"no arena"));
    free(b); free(c); ctools_string_arena_free(arena);
}

static char **store_expected;
static ST_retcode checked_string_store(ST_int v, ST_int row, char *value) {
    assert(v == 2 && row >= 2 && row <= nrows);
    const char *expected = store_expected[row - 2];
    assert(!strcmp(value, expected ? expected : ""));
    atomic_fetch_add(&writes, 1);
    return row == bad_row ? 459 : 0;
}
static void stores(void) {
    nrows = 6001; string_var = 2; strl = 0;
    size_t n = nrows - 1;
    char **strings = malloc(n * sizeof(char *));
    double *numbers = calloc(n, sizeof(double));
    assert(strings && numbers);
    for (size_t i = 0; i < n; i++) strings[i] = i%3 == 0 ? NULL : (i%3 == 1 ? "" : "text");
    store_expected = strings;
    mock.sstore = checked_string_store;
    stata_variable cols[2] = {
        {.type=STATA_TYPE_DOUBLE, .nobs=n, .data.dbl=numbers},
        {.type=STATA_TYPE_STRING, .nobs=n, .data.str=strings}
    };
    stata_data data = {.nobs=n, .nvars=2, .vars=cols};
    const size_t widths[] = {0, 4, 244, 2044, 2045};
    for (size_t i = 0; i < sizeof(widths)/sizeof(*widths); i++) {
        cols[1].str_maxlen = widths[i];
        assert(ctools_data_store(&data, 2) == 0);
        bad_row = 333;
        assert(ctools_data_store(&data, 2) == STATA_ERR_STATA_WRITE);
        bad_row = 0;
    }
    cols[1].str_maxlen = 3;
    assert(ctools_data_store(&data, 2) == STATA_ERR_STATA_WRITE);
    /* Above the former repack limit the legacy direct fallback is retained.
     * Fail the first callback to test its dispatch without allocating 256 MiB. */
    strings[0] = "text"; bad_row = 2; writes = 0;
    assert(store_single_variable(cols + 1, 2, 2, 256ULL*1024*1024/4 + 1, NULL) == STATA_ERR_STATA_WRITE);
    assert(writes == 1); bad_row = 0;
    free(strings); free(numbers);
}

static ST_retcode cancelled_read(ST_int v,ST_int i,ST_double *out) {(void)v;(void)i;(void)out;return 1;}
static ST_retcode cancelled_store(ST_int v,ST_int i,ST_double value) {(void)v;(void)i;(void)value;return 1;}
static void cancellation(void) {
    setup();nrows=4;string_var=0;strl=0;bad_row=0;pattern=0;
    mock.vdata=cancelled_read;mock.safevdata=cancelled_read;
    ctools_filtered_data fd;int col=1;
    assert(ctools_data_load(&fd,&col,1,1,4,CTOOLS_LOAD_SKIP_IF)==STATA_ERR_CANCELLED);
    ctools_filtered_data_free(&fd);
    double values[]={1,2,3,4};perm_idx_t map[]={1,2,3,4};
    mock.safestore=cancelled_store;mock.store=cancelled_store;
    assert(ctools_store_filtered(values,4,1,map)==STATA_ERR_CANCELLED);
}
int main(void) {
    calling_thread = pthread_self();
    setup(); mock.selobs = selection;
    maps();
    for (int threads = 1; threads <= 4; threads += 3) {
        ctools_set_max_threads(threads);
        values();
        stores();
        ctools_destroy_global_pool();
    }
    cancellation();
    return 0;
}
'''


def main():
    source_dir = Path(os.environ.get('CTOOLS_TEST_SOURCE', ROOT / 'src'))
    with tempfile.TemporaryDirectory(prefix='ctools-transport-native-') as tmp:
        directory = Path(tmp)
        source = directory/'test.c'
        source.write_text(SOURCE)
        for omp in (False, True):
            binary = directory/('test-'+str(omp))
            flags = ['-std=c11', '-O2', '-DSD_FASTMODE', '-fno-fast-math',
                     '-fsanitize='+os.environ.get('CTOOLS_SANITIZERS', 'undefined'), '-fno-omit-frame-pointer',
                     '-pthread', '-I'+str(source_dir)]
            if os.uname().sysname == 'Darwin':
                flags += ['-DSYSTEM=APPLEMAC', '-D_DARWIN_C_SOURCE',
                          '-Wl,-dead_strip', '-Wl,-undefined,dynamic_lookup']
                if omp:
                    prefix = os.environ.get('LIBOMP_PREFIX', '/opt/homebrew/opt/libomp')
                    flags += ['-Xpreprocessor', '-fopenmp', '-I'+prefix+'/include', prefix+'/lib/libomp.a']
            else:
                flags += ['-DSYSTEM=OPUNIX', '-D_GNU_SOURCE', '-ffunction-sections',
                          '-fdata-sections', '-Wl,--gc-sections', '-lm']
                if omp: flags += ['-fopenmp']
            subprocess.run([os.environ.get('CC', 'clang'), *flags, str(source),
                            *[str(source_dir/f) for f in ['ctools_types.c', 'ctools_arena.c', 'ctools_threads.c']],
                            '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True, timeout=120,
                           env=dict(os.environ, UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1"))
            print(f'PASS transport sanitizers={os.environ.get("CTOOLS_SANITIZERS", "undefined")} openmp={omp}', flush=True)


if __name__ == '__main__': main()

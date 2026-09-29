"""Adaptive transport regressions against an independent deterministic SPI model.

Requires the adaptive-load changes: opt-in CTOOLS_LOAD_NO_SORT_ORDER, numeric
row tiles for wide data, and packed maximum-width string arenas with retry.
An older source snapshot is expected to fail rather than silently skip them.

CTOOLS_TEST_SOURCE selects an isolated source tree. Serial and OpenMP builds
run under UBSan by default; CTOOLS_SANITIZERS=address,undefined also enables
ASan on systems with a working runtime. No Stata executable is invoked.
"""
from pathlib import Path
import os
import subprocess
import tempfile
from test_p2_native import SPI, ROOT

SOURCE = SPI + r'''
#include <stdio.h>
#include "ctools_config.h"
#include "ctools_arena.h"

/* Small string fixtures must exercise the same allocation retry and ownership
 * paths as large production loads. Numeric cell/row gates retain real sizes. */
#undef MIN_OBS_PER_THREAD
#define MIN_OBS_PER_THREAD 256

static int fail_aligned;
static atomic_int fail_flat;
static atomic_int fail_large_arena, failed_large_seen, retry_seen;
static atomic_int fail_all_arenas, arena_failures;
static size_t expected_retry;

static void *test_aligned2(size_t count, size_t size)
{
    if (count > 512 && size == sizeof(double) &&
        fail_aligned > 0 && --fail_aligned == 0) return NULL;
    return ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, count, size);
}

static void *test_calloc(size_t count, size_t size)
{
    if (count > 512 && size == 65 && atomic_exchange(&fail_flat, 0)) {
        return NULL;
    }
    return calloc(count, size);
}

static ctools_string_arena *test_arena_create(size_t capacity,
                                              ctools_string_arena_mode mode)
{
    if (expected_retry && capacity == expected_retry)
        atomic_fetch_add(&retry_seen, 1);
    if (atomic_load(&fail_all_arenas)) {
        atomic_fetch_add(&arena_failures, 1);
        return NULL;
    }
    if (capacity > 1000000 && atomic_exchange(&fail_large_arena, 0)) {
        atomic_fetch_add(&failed_large_seen, 1);
        return NULL;
    }
    return ctools_string_arena_create(capacity, mode);
}

#undef ctools_safe_cacheline_alloc2
#define ctools_safe_cacheline_alloc2(count, size) test_aligned2((count), (size))
#define calloc test_calloc
#define ctools_string_arena_create test_arena_create
#include "ctools_data_io.c"
#undef ctools_string_arena_create
#undef calloc

#ifndef CTOOLS_LOAD_NO_SORT_ORDER
#error "Adaptive transport tests require CTOOLS_LOAD_NO_SORT_ORDER"
#endif

static int ncols, shape, width = 2045, content_width = -1, empty_selection, bad_var;
static const perm_idx_t *store_sources;

static ST_int variable_count(void) { return ncols; }
static ST_boolean string_column(ST_int v)
{
    return shape == 2 || (shape == 1 && v % 2 == 0);
}
static ST_boolean fixed_column(ST_int v) { (void)v; return 0; }
static ST_boolean selection(ST_int row)
{
    return !empty_selection && (!filter || row % 2 != 0);
}
static double number_value(int v, int row)
{
    return row + 0.125 * v;
}
static void text_value(int v, int row, char *out)
{
    if (row % 17 == 0) { out[0] = '\0'; return; }
    int length = content_width >= 0 ? content_width : width;
    memset(out, 'a' + (row + v) % 26, (size_t)length);
    if (length >= 2) { out[0] = (char)0xc3; out[1] = (char)0xa9; }
    out[length] = '\0';
}
static ST_retcode read_number(ST_int v, ST_int row, double *out)
{
    assert(v >= 1 && v <= ncols && row >= 1 && row <= nrows);
    if (row == bad_row && (!bad_var || v == bad_var)) return 459;
    *out = number_value(v, row);
    return 0;
}
static ST_retcode read_text(ST_int v, ST_int row, char *out)
{
    assert(v >= 1 && v <= ncols && row >= 1 && row <= nrows);
    if (row == bad_row && (!bad_var || v == bad_var)) return 459;
    text_value(v, row, out);
    return 0;
}
static ST_retcode write_number(ST_int v, ST_int row, double value)
{
    assert(store_sources && v >= 1 && v <= ncols);
    assert(row >= 1 && row <= nrows);
    assert(value == number_value(v, (int)store_sources[row - 1]));
    atomic_fetch_add(&writes, 1);
    return 0;
}
static ST_retcode write_text(ST_int v, ST_int row, char *value)
{
    assert(store_sources && v >= 1 && v <= ncols);
    assert(row >= 1 && row <= nrows);
    char expected[2046];
    text_value(v, (int)store_sources[row - 1], expected);
    assert(!strcmp(value, expected));
    atomic_fetch_add(&writes, 1);
    return 0;
}
static ST_retcode width_macro(char *name, char *out, ST_int capacity)
{
    (void)name;
    size_t used = 0;
    for (int j = 1; j <= ncols; j++) {
        int bytes = snprintf(out + used, (size_t)capacity - used,
                             "%s%d", j > 1 ? "," : "", string_column(j) ? width : 0);
        assert(bytes > 0 && used + (size_t)bytes < (size_t)capacity);
        used += (size_t)bytes;
    }
    return 0;
}

static size_t count_selected(int first, int last)
{
    size_t count = 0;
    for (int row = first; row <= last; row++) count += selection(row) != 0;
    return count;
}

static void check_values(const ctools_filtered_data *fd, int first, int last,
                         int no_order)
{
    size_t out = 0;
    assert(fd->data.nvars == (size_t)ncols);
    if (no_order) assert(fd->data.sort_order == NULL);
    for (int row = first; row <= last; row++) if (selection(row)) {
        assert(fd->obs_map[out] == (perm_idx_t)row);
        if (!no_order) assert(fd->data.sort_order[out] == out);
        for (int j = 1; j <= ncols; j++) {
            const stata_variable *v = fd->data.vars + j - 1;
            if (string_column(j)) {
                char value[2046];
                text_value(j, row, value);
                assert(v->type == STATA_TYPE_STRING && !strcmp(v->data.str[out], value));
            } else {
                assert(v->type == STATA_TYPE_DOUBLE && v->data.dbl[out] == number_value(j, row));
            }
        }
        out++;
    }
    assert(fd->data.nobs == out);
}

/* Real numeric gate boundaries: 128 columns with exactly 50,000 selected
 * rows, and 256 columns with 60,000 selected rows, plus filtered equivalents. Fault injection proves that the speculative tiled allocation
 * was reached in OpenMP builds and falls back without changing any value. */
static void numeric_tiles(void)
{
    const int rows[] = {50004, 100003, 60004, 120003};
    const int columns[] = {128, 128, 256, 256};
    shape = 0;
    for (int c = 0; c < 4; c++) {
        nrows = rows[c]; ncols = columns[c]; filter = c % 2;
        for (int no_order = 0; no_order <= 1; no_order++) {
            int flags = no_order ? CTOOLS_LOAD_NO_SORT_ORDER : CTOOLS_LOAD_CHECK_IF;
            for (int fail = 0; fail <= 3; fail++) {
                fail_aligned = fail;
                ctools_filtered_data fd;
                assert(ctools_data_load(&fd, NULL, 0, 3, nrows - 2, flags) == STATA_OK);
                check_values(&fd, 3, nrows - 2, no_order);
                #ifdef _OPENMP
                assert(fail_aligned == 0);
                #endif
                ctools_filtered_data_free(&fd);
            }
            fail_aligned = 0;
            bad_row = 777;
            ctools_filtered_data fd;
            assert(ctools_data_load(&fd, NULL, 0, 3, nrows - 2, flags) == STATA_ERR_STATA_READ);
            assert(!fd.data.vars && !fd.data.sort_order && !fd.obs_map);
            ctools_filtered_data_free(&fd);
            bad_row = 0;
        }
    }
    fail_aligned = 0;
}

/* Below each numeric gate, no speculative double-column allocation should
 * occur. The interposer targets count>512 with sizeof(double): vars metadata
 * has a different element size, maps/orders use perm_idx_t, and ordinary
 * column loads allocate bytes through ctools_cacheline_alloc instead. These
 * fixtures are also too short to enable the separate chunked-load scheduler. */
static void numeric_below_tile_gates(void)
{
    const int rows[] = {50004, 100003, 50003, 100001, 10004, 20003};
    const int columns[] = {127, 127, 128, 128, 2000, 2000};
    const size_t selected_rows[] = {50000, 50000, 49999, 49999, 10000, 10000};
    shape = 0;
    for (int c = 0; c < 6; c++) {
        nrows = rows[c]; ncols = columns[c]; filter = c % 2;
        for (int no_order = 0; no_order <= 1; no_order++) {
            int flags = no_order ? CTOOLS_LOAD_NO_SORT_ORDER : CTOOLS_LOAD_CHECK_IF;
            ctools_filtered_data fd;
            fail_aligned = 1;
            assert(ctools_data_load(&fd, NULL, 0, 3, nrows - 2, flags) == STATA_OK);
            assert(fail_aligned == 1);
            assert(fd.data.nobs == selected_rows[c]);
            check_values(&fd, 3, nrows - 2, no_order);
            ctools_filtered_data_free(&fd);
            assert(!fd.data.vars && !fd.obs_map && !fd.data.sort_order);
            fail_aligned = 0;
        }
    }
}

/* Full/short/empty maximum-width data remain densely owned by an arena.
 * Rejecting one large reservation must retry the original 64*N estimate,
 * including row-parallel and ordinary/filtered column paths. */
static void packed_strings(void)
{
    const int thread_cases[] = {1, 4, 12};
    nrows = 1109; shape = 2; width = 2045;
    for (size_t t = 0; t < sizeof(thread_cases) / sizeof(*thread_cases); t++) {
        ctools_set_max_threads(thread_cases[t]);
        for (ncols = 1; ncols <= 4; ncols += 3)
        for (filter = 0; filter <= 1; filter++) {
            ctools_filtered_data fd;
            expected_retry = count_selected(3, nrows - 2) * 64;
            fail_large_arena = 1;
            failed_large_seen = retry_seen = 0;
            assert(ctools_data_load(&fd, NULL, 0, 3, nrows - 2, 0) == STATA_OK);
            assert(failed_large_seen == 1 && retry_seen == 1);
            check_values(&fd, 3, nrows - 2, 0);
            ctools_filtered_data_free(&fd);
            ctools_filtered_data_free(&fd);
            expected_retry = 0;
        }
        ctools_destroy_global_pool();
    }
    /* A short actual value must occupy its actual bytes inside the packed
     * arena, even though the declared maximum width remains 2045. */
    const int lengths[] = {0, 8, 64, 2045};
    ctools_set_max_threads(4); ncols = 4;
    for (size_t length = 0; length < sizeof(lengths) / sizeof(*lengths); length++)
    for (filter = 0; filter <= 1; filter++) {
        content_width = lengths[length];
        size_t used = 0;
        for (int row = 3; row <= nrows - 2; row++) if (selection(row))
            used += row % 17 == 0 ? 1 : (size_t)content_width + 1;
        ctools_filtered_data fd;
        assert(ctools_data_load(&fd, NULL, 0, 3, nrows - 2, 0) == STATA_OK);
        check_values(&fd, 3, nrows - 2, 0);
        for (int j = 0; j < ncols; j++) {
            ctools_string_arena *arena = fd.data.vars[j]._arena;
            assert(arena && arena->used == used && !arena->has_fallback);
        }
        ctools_filtered_data_free(&fd);
    }
    content_width = -1;
    /* If both reservations fail, every string is individually owned. Check
     * normal release and release of a partly read column after an SPI error.
     * Leak-enabled ASan builds additionally check for leaked strings. */
    for (ncols = 1; ncols <= 4; ncols += 3)
    for (filter = 0; filter <= 1; filter++) {
        ctools_filtered_data fd;
        expected_retry = count_selected(3, nrows - 2) * 64;
        fail_all_arenas = 1; arena_failures = retry_seen = 0;
        assert(ctools_data_load(&fd, NULL, 0, 3, nrows - 2,
                               CTOOLS_LOAD_NO_SORT_ORDER) == STATA_OK);
        assert(arena_failures == 2 * ncols && retry_seen == ncols);
        for (int j = 0; j < ncols; j++) assert(fd.data.vars[j]._arena == NULL);
        check_values(&fd, 3, nrows - 2, 1);
        ctools_filtered_data_free(&fd);
        ctools_filtered_data_free(&fd);
        assert(!fd.data.vars && !fd.obs_map && !fd.data.sort_order);
        bad_row = 777; arena_failures = retry_seen = 0;
        assert(ctools_data_load(&fd, NULL, 0, 3, nrows - 2,
                               CTOOLS_LOAD_NO_SORT_ORDER) == STATA_ERR_STATA_READ);
        assert(arena_failures >= 2 && retry_seen >= 1);
        assert(!fd.data.vars && !fd.obs_map && !fd.data.sort_order);
        ctools_filtered_data_free(&fd);
        bad_row = 0; fail_all_arenas = 0; expected_retry = 0;
    }
    ctools_destroy_global_pool();
}

/* Default callers retain identity order. Opt-out callers retain correct
 * values/maps, allow ordinary stores, and reject sorted stores before writes.
 * Ownership also survives empty selections and the tiled-memory fallback. */
static void optional_order(void)
{
    ctools_set_max_threads(4);
    nrows = 1109; ncols = 2; width = 8;
    for (shape = 0; shape <= 1; shape++)
    for (filter = 0; filter <= 1; filter++)
    for (int no_order = 0; no_order <= 1; no_order++) {
        ctools_filtered_data fd;
        int flags = no_order ? CTOOLS_LOAD_NO_SORT_ORDER : 0;
        assert(ctools_data_load(&fd, NULL, 0, 3, nrows - 2, flags) == STATA_OK);
        check_values(&fd, 3, nrows - 2, no_order);
        store_sources = fd.obs_map;
        writes = 0;
        assert(ctools_data_store(&fd.data, 1) == STATA_OK);
        assert((size_t)writes == fd.data.nobs * fd.data.nvars);
        writes = 0;
        assert(ctools_data_store_sorted(&fd.data, 1) ==
               (no_order ? STATA_ERR_INVALID_INPUT : STATA_OK));
        assert((size_t)writes == (no_order ? 0 : fd.data.nobs * fd.data.nvars));
        store_sources = NULL;
        ctools_filtered_data_free(&fd);
        ctools_filtered_data_free(&fd);
        assert(!fd.data.vars && !fd.data.sort_order && !fd.obs_map);
    }
    ctools_filtered_data fd;
    empty_selection = 1;
    assert(ctools_data_load(&fd, NULL, 0, 1, nrows, CTOOLS_LOAD_NO_SORT_ORDER) == STATA_OK);
    assert(!fd.data.nobs && !fd.data.sort_order && !fd.obs_map);
    ctools_filtered_data_free(&fd);
    empty_selection = 0; nrows = 0;
    assert(ctools_data_load(&fd, NULL, 0, 0, 0, CTOOLS_LOAD_NO_SORT_ORDER) == STATA_OK);
    assert(!fd.data.nobs && !fd.data.sort_order && !fd.obs_map);
    ctools_filtered_data_free(&fd);

    nrows = 1109; ncols = 4; shape = 2; filter = 0; width = 64;
    fail_flat = 1;
    assert(ctools_data_load(&fd, NULL, 0, 1, nrows,
           CTOOLS_LOAD_SKIP_IF | CTOOLS_LOAD_NO_SORT_ORDER) == STATA_OK);
    assert(!fail_flat);
    check_values(&fd, 1, nrows, 1);
    ctools_filtered_data_free(&fd);
    ctools_destroy_global_pool();
}


/* Reordered mixed varlists must preserve output positions when the strL
 * column is removed from the worker batch and read on the calling thread. */
static pthread_t loading_thread;
static int long_reads;
static ST_boolean reordered_string(ST_int v) { return v == 2 || v == 3; }
static ST_boolean reordered_long(ST_int v) { return v == 3; }
static ST_boolean long_binary(ST_int v, ST_int row)
{
    assert(v == 3 && row >= 1 && row <= nrows);
    assert(pthread_equal(pthread_self(), loading_thread));
    return 0;
}
static ST_int long_length(ST_int v, ST_int row)
{
    assert(v == 3 && row >= 1 && row <= nrows);
    assert(pthread_equal(pthread_self(), loading_thread));
    return row % 17 == 0 ? 0 : width;
}
static ST_retcode read_long(ST_int v, ST_int row, char *out, ST_int capacity)
{
    assert(v == 3 && pthread_equal(pthread_self(), loading_thread));
    int length = long_length(v, row);
    assert(capacity == length + 1);
    long_reads++;
    return read_text(v, row, out) ? -1 : length;
}

static void reordered_long_without_order(void)
{
    const int orders[][4] = {{3, 4, 2, 1}, {2, 3, 1, 4}, {4, 2, 1, 3}};
    const int threads[] = {1, 4, 12};
    const int widths[] = {0, 64, 0, 0}; /* Indexed by SPI variable, not load order. */
    ncols = 4; nrows = 6007; width = 64; content_width = -1;
    loading_thread = pthread_self();
    mock.isstr = reordered_string; mock.isstrl = reordered_long;
    mock.isbinary = long_binary; mock.sdatalen = long_length; mock.strldata = read_long;
    for (size_t t = 0; t < sizeof(threads) / sizeof(*threads); t++) {
        ctools_set_max_threads(threads[t]);
        for (size_t order = 0; order < sizeof(orders) / sizeof(*orders); order++)
        for (filter = 0; filter <= 1; filter++) {
            int indices[4];
            memcpy(indices, orders[order], sizeof(indices));
            ctools_filtered_data fd;
            long_reads = 0;
            assert(ctools_data_load_ex(&fd, indices, 4, 3, nrows - 2,
                                      CTOOLS_LOAD_NO_SORT_ORDER, widths) == STATA_OK);
            assert(!fd.data.sort_order && fd.data.nvars == 4);
            assert((size_t)long_reads == fd.data.nobs);
            size_t out = 0;
            for (int row = 3; row <= nrows - 2; row++) if (selection(row)) {
                assert(fd.obs_map[out] == (perm_idx_t)row);
                for (int j = 0; j < 4; j++) {
                    int v = indices[j];
                    if (reordered_string(v)) {
                        char value[2046];
                        text_value(v, row, value);
                        assert(fd.data.vars[j].type == STATA_TYPE_STRING);
                        assert(!strcmp(fd.data.vars[j].data.str[out], value));
                    } else {
                        assert(fd.data.vars[j].type == STATA_TYPE_DOUBLE);
                        assert(fd.data.vars[j].data.dbl[out] == number_value(v, row));
                    }
                }
                out++;
            }
            assert(out == fd.data.nobs);
            ctools_filtered_data_free(&fd);
            for (bad_var = 1; bad_var <= 4; bad_var++) {
                bad_row = 777;
                assert(ctools_data_load_ex(&fd, indices, 4, 3, nrows - 2,
                                          CTOOLS_LOAD_NO_SORT_ORDER, widths) == STATA_ERR_STATA_READ);
                assert(!fd.data.vars && !fd.obs_map && !fd.data.sort_order);
                ctools_filtered_data_free(&fd);
                bad_row = 0;
            }
            bad_var = 0;
        }
        ctools_destroy_global_pool();
    }
    mock.isstr = string_column; mock.isstrl = fixed_column;
}

int main(void)
{
    setup();
    mock.nvars = mock.nvar = variable_count;
    mock.isstr = string_column; mock.isstrl = fixed_column;
    mock.selobs = selection;
    mock.safevdata = mock.vdata = read_number;
    mock.sdata = read_text; mock.macuse = width_macro;
    mock.safestore = write_number; mock.sstore = write_text;
    ctools_set_max_threads(4);
    numeric_tiles();
    numeric_below_tile_gates();
    ctools_destroy_global_pool();
    packed_strings();
    optional_order();
    reordered_long_without_order();
    return 0;
}
'''


def main():
    source_dir = Path(os.environ.get('CTOOLS_TEST_SOURCE', ROOT / 'src'))
    sanitizers = os.environ.get('CTOOLS_SANITIZERS', 'undefined')
    with tempfile.TemporaryDirectory(prefix='ctools-adaptive-') as tmp:
        directory = Path(tmp)
        source = directory / 'test.c'
        source.write_text(SOURCE)
        for omp in (False, True):
            binary = directory / ('test-' + str(omp))
            flags = ['-std=c11', '-O2', '-DSD_FASTMODE', '-fno-fast-math',
                     '-fsanitize=' + sanitizers, '-fno-omit-frame-pointer',
                     '-pthread', '-I' + str(source_dir)]
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
                            *[str(source_dir / name) for name in
                              ['ctools_types.c', 'ctools_arena.c', 'ctools_threads.c']],
                            '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True, timeout=120,
                           env=dict(os.environ, UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1',
                                    KMP_BLOCKTIME='0', OMP_WAIT_POLICY='PASSIVE'))
            print(f'PASS adaptive transport sanitizers={sanitizers} '
                  f'openmp={omp} source={source_dir}', flush=True)


if __name__ == '__main__':
    main()

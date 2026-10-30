/* Transport tests. Select one CTOOLS_CASE_* per compilation.
 * Separate selectors preserve allocator/SPI interposition and independent main functions.
 */

#if defined(CTOOLS_CASE_bigdata_native_1)
#include "ctools_config.h"
#undef MIN_OBS_PER_THREAD
#define MIN_OBS_PER_THREAD 256
#include "ctools_data_io.c"
#include "ctools_transport.c"
#include "ctools_data_plan.c"
#include "ctools_data_store.c"
static int offset;
static ST_retcode ordered_number(ST_int v, ST_int i, ST_double value) {
    assert(value == nrows - i + offset - 1);
    return putnum(v, i, value);
}
static ST_retcode ordered_string(ST_int v, ST_int i, char *value) {
    size_t row = nrows - i + offset - 2;
    assert(!strcmp(value, row % 3 == 0 ? "text" : ""));
    return putstr(v, i, value);
}
int main(void) {
    setup();
    const size_t n = 6000;
    nrows = n + 1;
    double *numbers = malloc(n * sizeof(double));
    char **strings = malloc(n * sizeof(char *));
    perm_idx_t *order = malloc(n * sizeof(perm_idx_t));
    assert(numbers && strings && order);
    for (size_t i = 0; i < n; i++) {
        numbers[i] = i + 1;
        strings[i] = i % 3 == 0 ? "text" : (i % 3 == 1 ? NULL : "");
        order[i] = n - i - 1;
    }
    stata_variable columns[2] = {
        {.type=STATA_TYPE_DOUBLE, .nobs=n, .data.dbl=numbers},
        /* A wide declared width makes the string gather cross block boundaries. */
        {.type=STATA_TYPE_STRING, .nobs=n, .data.str=strings, .str_maxlen=128}
    };
    stata_data data = {.nobs=n, .nvars=2, .vars=columns, .sort_order=order};
    mock.safestore = ordered_number;
    mock.sstore = ordered_string;
    for (int pool = 0; pool <= 1; pool++) {
        ctools_set_max_threads(pool ? 4 : 1);
        for (int kind = 0; kind < 3; kind++) {
            data.vars = kind == 1 ? columns + 1 : columns;
            data.nvars = kind == 2 ? 2 : 1;
            string_var = kind == 0 ? 0 : (kind == 1 ? 1 : 2);
            for (offset = 1; offset <= 2; offset++) {
                bad_row = 0; writes = 0;
                assert(ctools_data_store_sorted(&data, offset) == 0);
                assert(writes == (int)(n * data.nvars));
                for (size_t i = 0; i < n; i++) {
                    assert(order[i] == n - i - 1 && numbers[i] == i + 1);
                    if (i % 3 == 1) assert(strings[i] == NULL);
                    else assert(!strcmp(strings[i], i % 3 == 0 ? "text" : ""));
                }
                bad_row = 3;
                assert(ctools_data_store_sorted(&data, offset) == STATA_ERR_STATA_WRITE);
                bad_row = 0; writes = 0; order[n - 1] = n;
                assert(ctools_data_store_sorted(&data, offset) == STATA_ERR_INVALID_INPUT);
                assert(writes == 0); order[n - 1] = 0;
                data.sort_order = NULL;
                assert(ctools_data_store_sorted(&data, offset) == STATA_ERR_INVALID_INPUT);
                assert(writes == 0); data.sort_order = order;
                assert(ctools_data_store_sorted(&data, 3) == STATA_ERR_INVALID_INPUT);
                assert(writes == 0);
            }
        }
        ctools_destroy_global_pool();
    }
    free(numbers); free(strings); free(order);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_rowstage_native)
#include <stdio.h>
#include "ctools_config.h"
#undef MIN_OBS_PER_THREAD
#define MIN_OBS_PER_THREAD 256
/* Let a few thousand observations reach the multi-worker paths. */
#define ROWSTAGE_CELLS_PER_THREAD 512
static int fail_rows;
static atomic_int row_buffers;
static void *test_rows_malloc(size_t bytes) {
    /* Only the full-range staged-row buffer is this large in these cases. */
    if (bytes > 500000) {
        atomic_fetch_add(&row_buffers, 1);
        if (fail_rows) return NULL;
    }
    return malloc(bytes);
}
#define malloc test_rows_malloc
#include "ctools_data_io.c"
#include "ctools_transport.c"
#include "ctools_data_plan.c"
#include "ctools_data_store.c"
#include "ctools_rowstage.c"
#undef malloc

/* Seven plugin variables: numeric, narrow string (bounced reads), numeric,
 * a string wide enough for direct slot reads, numeric, string, strL. */
#define COLS 7
#define ROWS 6007
#define MAXW 64
static const int width_of[COLS] = {0, 5, 0, 40, 0, MAXW, 0};
static char stored_text[COLS][ROWS + 1][MAXW + 1];
static double stored_number[COLS][ROWS + 1];
static atomic_int text_writes, number_writes;
static int long_row, cancel_row, strl_var;
static ST_int many_vars(void) { return COLS; }
static ST_boolean string_column(ST_int v) { return width_of[v - 1] > 0 || v == strl_var; }
static ST_boolean strl_column(ST_int v) { return v == strl_var; }
static void text_value(int v, int row, char *out) {
    int width = width_of[v - 1];
    if (row % 17 == 0) { out[0] = 0; return; }
    /* Lengths cycle from 1 to the full width, ending in a UTF-8 pair. */
    int len = 1 + (row + v) % width;
    memset(out, 'a' + (row + v) % 26, (size_t)len);
    if (len >= 2) { out[len - 2] = (char)0xc3; out[len - 1] = (char)0xa9; }
    out[len] = 0;
}
static double number_value(int v, int row) {
    return row % 23 == 0 ? mock.missval : row % 29 == 0 ? -0.0 : row + v * 0.125;
}
static ST_retcode read_text(ST_int v, ST_int row, char *out) {
    assert(v >= 1 && v <= COLS && width_of[v - 1] > 0);
    assert(row >= 1 && row <= nrows);
    if (row == bad_row) return 459;
    if (row == cancel_row) return 1;
    if (row == long_row) { memset(out, 'L', (size_t)width_of[v - 1] + 3); out[width_of[v - 1] + 3] = 0; return 0; }
    text_value(v, row, out); return 0;
}
static ST_retcode read_number(ST_int v, ST_int row, double *out) {
    assert(v >= 1 && v <= COLS && width_of[v - 1] == 0 && v != strl_var);
    assert(row >= 1 && row <= nrows);
    if (row == bad_row) return 459;
    if (row == cancel_row) return 1;
    *out = number_value(v, row); return 0;
}
static ST_retcode write_text(ST_int v, ST_int row, char *value) {
    assert(v >= 1 && v <= COLS && width_of[v - 1] > 0);
    assert(row >= 1 && row <= nrows && strlen(value) <= (size_t)width_of[v - 1]);
    atomic_fetch_add(&text_writes, 1);
    if (row == bad_row) return 459;
    if (row == cancel_row) return 1;
    strcpy(stored_text[v - 1][row], value); return 0;
}
static ST_retcode write_number(ST_int v, ST_int row, double value) {
    assert(v >= 1 && v <= COLS && width_of[v - 1] == 0);
    assert(row >= 1 && row <= nrows);
    atomic_fetch_add(&number_writes, 1);
    if (row == bad_row) return 459;
    if (row == cancel_row) return 1;
    stored_number[v - 1][row] = value; return 0;
}
static ST_retcode unchecked_write(ST_int v, ST_int row, double value) {
    (void)v; (void)row; (void)value;
    assert(!"numeric stores must use the checked SPI callback");
    return 459;
}

static int widths[COLS];
static const int staged_vars[] = {6, 1, 4, 3, 2, 5};  /* reordered, without the last variable */
#define NSTAGED 6

static void clear_stored(void) {
    memset(stored_text, 0x7f, sizeof(stored_text));
    for (int v = 0; v < COLS; v++)
        for (int row = 0; row <= ROWS; row++) stored_number[v][row] = -777.0;
    text_writes = number_writes = 0;
}

/* Destination row dest (one-based) must hold source observation src (one-based), or blanks. */
static void check_row(int dest, int src) {
    for (int j = 0; j < NSTAGED; j++) {
        int v = staged_vars[j];
        if (width_of[v - 1]) {
            char expected[MAXW + 1] = "";
            if (src > 0) text_value(v, src, expected);
            assert(!strcmp(stored_text[v - 1][dest], expected));
        } else {
            double expected = src > 0 ? number_value(v, src) : mock.missval;
            assert(memcmp(&stored_number[v - 1][dest], &expected, sizeof(double)) == 0);
        }
    }
    /* The variable that was not staged is never written. */
    assert(stored_number[6][dest] == -777.0);
}

static void eligibility(void) {
    int numeric_only[] = {1, 3, 5}, with_strl[] = {1, 7}, bad_index[] = {1, COLS + 1};
    int no_hint[COLS], wide_hint[COLS];
    memcpy(no_hint, widths, sizeof(widths)); no_hint[3] = 0;
    memcpy(wide_hint, widths, sizeof(widths)); wide_hint[3] = STATA_STR_MAXLEN;
    size_t n = (size_t)nrows;
    assert(ctools_rowstage_supported(staged_vars, NSTAGED, widths, n));
    assert(ctools_rowstage_supported(numeric_only, 3, NULL, n));
    /* Without the caller's width metadata, strings cannot be staged. */
    assert(!ctools_rowstage_supported(staged_vars, NSTAGED, NULL, n));
    assert(!ctools_rowstage_supported(staged_vars, NSTAGED, no_hint, n));
    assert(!ctools_rowstage_supported(staged_vars, NSTAGED, wide_hint, n));
    assert(!ctools_rowstage_supported(bad_index, 2, widths, n));
    assert(!ctools_rowstage_supported(staged_vars, 0, widths, n));
    /* A string column whose declared-width slots would pass 2 GiB stays with
     * the column transport; numeric columns have no such bound. The widest
     * staged string is MAXW bytes plus its terminator. */
    size_t cap_rows = ((size_t)2048 * 1024 * 1024) / (MAXW + 1);
    assert(ctools_rowstage_supported(staged_vars, NSTAGED, widths, cap_rows));
    assert(!ctools_rowstage_supported(staged_vars, NSTAGED, widths, cap_rows + 1));
    assert(ctools_rowstage_supported(numeric_only, 3, widths, cap_rows + 1));
    strl_var = 7;
    assert(!ctools_rowstage_supported(with_strl, 2, widths, n));
    ctools_rowstage rs;
    assert(ctools_rowstage_load(&rs, with_strl, 2, widths, 1, 100) == STATA_ERR_UNSUPPORTED_TYPE);
    assert(!rs.rows && !rs.var_idx && !rs.width && !rs.offset);
    strl_var = 0;
}

static void round_trips(size_t first, size_t nobs) {
    ctools_rowstage rs;
    assert(ctools_rowstage_load(&rs, staged_vars, NSTAGED, widths, first, nobs) == STATA_OK);
    assert(rs.nobs == nobs && rs.nvars == NSTAGED && rs.stride % 8 == 0);

    /* Columns come back contiguous, in the loader's layout, zero padded. */
    for (size_t j = 0; j < NSTAGED; j++) {
        int v = staged_vars[j];
        stata_variable column;
        assert(ctools_rowstage_column(&rs, j, &column) == STATA_OK && column.nobs == nobs);
        for (size_t i = 0; i < nobs; i++) {
            int row = (int)(first + i);
            if (width_of[v - 1]) {
                char expected[MAXW + 2] = {0};
                text_value(v, row, expected);
                assert(column.type == STATA_TYPE_STRING);
                assert(column.str_maxlen == (size_t)width_of[v - 1]);
                assert(memcmp(column.data.str[i], expected, (size_t)width_of[v - 1] + 1) == 0);
            } else {
                double expected = number_value(v, row);
                assert(column.type == STATA_TYPE_DOUBLE);
                assert(memcmp(&column.data.dbl[i], &expected, sizeof(double)) == 0);
            }
        }
        stata_variable_free(&column);
    }
    stata_variable none;
    assert(ctools_rowstage_column(&rs, NSTAGED, &none) == STATA_ERR_INVALID_INPUT);

    /* In order, to a shifted destination where the range allows one. */
    size_t dest = nobs < (size_t)nrows ? 2 : 1;
    clear_stored();
    assert(ctools_rowstage_store(&rs, NULL, nobs, dest) == STATA_OK);
    for (size_t i = 0; i < nobs; i++) check_row((int)(dest + i), (int)(first + i));
    assert((size_t)number_writes == 3 * nobs && (size_t)text_writes == 3 * nobs);
    assert(stored_number[0][dest - 1] == -777.0);
    if (dest + nobs <= (size_t)nrows) assert(stored_number[0][dest + nobs] == -777.0);
    assert(ctools_rowstage_store(&rs, NULL, nobs + 1, 1) == STATA_ERR_INVALID_INPUT);

    /* Reordered, repeated and absent sources; more output rows than staged rows. */
    size_t n_out = (size_t)nrows - 1;
    perm_idx_t *source = malloc(n_out * sizeof(*source));
    assert(source);
    for (int pattern = 0; pattern < 4; pattern++) {
        for (size_t d = 0; d < n_out; d++) {
            source[d] = pattern == 0 ? (perm_idx_t)(nobs - 1 - d % nobs) :      /* reversed, wrapping */
                        pattern == 1 ? (perm_idx_t)((d * 2654435761u) % nobs) : /* scattered */
                        pattern == 2 ? (perm_idx_t)(d / 3 % nobs) :             /* runs of duplicates */
                        d % 5 == 0 ? CTOOLS_ROW_NONE : (perm_idx_t)(d % nobs);  /* missing sources */
        }
        clear_stored();
        assert(ctools_rowstage_store(&rs, source, n_out, 2) == STATA_OK);
        for (size_t d = 0; d < n_out; d++)
            check_row((int)(2 + d), source[d] == CTOOLS_ROW_NONE ? 0 : (int)(first + source[d]));
        assert((size_t)number_writes == 3 * n_out && (size_t)text_writes == 3 * n_out);
    }

    /* Rejected before anything is written. */
    clear_stored();
    source[n_out / 2] = (perm_idx_t)nobs;
    assert(ctools_rowstage_store(&rs, source, n_out, 2) == STATA_ERR_INVALID_INPUT);
    source[n_out / 2] = 0;
    assert(ctools_rowstage_store(&rs, source, n_out, 3) == STATA_ERR_INVALID_INPUT);
    assert(ctools_rowstage_store(&rs, source, n_out, 0) == STATA_ERR_INVALID_INPUT);
    assert(number_writes == 0 && text_writes == 0);
    assert(ctools_rowstage_store(&rs, source, 0, 1) == STATA_OK);

    /* Write failures and cancellation keep their status. */
    bad_row = 2 + (int)n_out / 2;
    assert(ctools_rowstage_store(&rs, source, n_out, 2) == STATA_ERR_STATA_WRITE);
    bad_row = 0; cancel_row = 2 + (int)n_out / 3;
    assert(ctools_rowstage_store(&rs, source, n_out, 2) == STATA_ERR_CANCELLED);
    cancel_row = 0;

    free(source);
    ctools_rowstage_free(&rs);
    ctools_rowstage_free(&rs);  /* idempotent */
    assert(!rs.rows);
}

static void load_failures(void) {
    ctools_rowstage rs;
    int stale[COLS];
    bad_row = 4000;
    assert(ctools_rowstage_load(&rs, staged_vars, NSTAGED, widths, 1, (size_t)nrows) == STATA_ERR_STATA_READ);
    assert(!rs.rows && !rs.var_idx);
    bad_row = 0; cancel_row = 123;
    assert(ctools_rowstage_load(&rs, staged_vars, NSTAGED, widths, 1, (size_t)nrows) == STATA_ERR_CANCELLED);
    assert(!rs.rows);
    cancel_row = 0;
    /* A value longer than its hint fails the load, for bounced and direct reads,
     * including in the last rows, where direct reads fall back to bounded ones. */
    int long_rows[] = {1, 2999, ROWS - 1, ROWS};
    for (int k = 0; k < 4; k++) {
        long_row = long_rows[k];
        assert(ctools_rowstage_load(&rs, staged_vars, NSTAGED, widths, 1, (size_t)nrows) == STATA_ERR_STATA_READ);
        assert(!rs.rows);
    }
    long_row = 0;
    /* A hint that is too wide is harmless; one that is too narrow is caught. */
    memcpy(stale, widths, sizeof(widths)); stale[1] = 9;
    assert(ctools_rowstage_load(&rs, staged_vars, NSTAGED, stale, 1, 500) == STATA_OK);
    ctools_rowstage_free(&rs);
    stale[1] = 2;
    assert(ctools_rowstage_load(&rs, staged_vars, NSTAGED, stale, 1, 500) == STATA_ERR_STATA_READ);
    /* Ranges and empty requests. */
    assert(ctools_rowstage_load(&rs, staged_vars, NSTAGED, widths, 2, (size_t)nrows) == STATA_ERR_INVALID_INPUT);
    assert(ctools_rowstage_load(&rs, staged_vars, NSTAGED, widths, 0, 10) == STATA_ERR_INVALID_INPUT);
    assert(ctools_rowstage_load(&rs, staged_vars, NSTAGED, widths, 1, 0) == STATA_ERR_INVALID_INPUT);
    assert(ctools_rowstage_load(&rs, NULL, NSTAGED, widths, 1, 10) == STATA_ERR_INVALID_INPUT);
    /* No memory for the staged rows: reported before any read. */
    reads = 0; row_buffers = 0; fail_rows = 1;
    assert(ctools_rowstage_load(&rs, staged_vars, NSTAGED, widths, 1, (size_t)nrows) == STATA_ERR_MEMORY);
    assert(row_buffers == 1 && !rs.rows && !rs.var_idx);
    fail_rows = 0;
}

int main(void) {
    setup(); nrows = ROWS;
    mock.nvars = mock.nvar = many_vars;
    mock.isstr = string_column; mock.isstrl = strl_column;
    mock.safevdata = read_number; mock.vdata = read_number;
    mock.safestore = write_number; mock.store = unchecked_write;
    mock.sdata = read_text; mock.sstore = write_text;
    for (int j = 0; j < COLS; j++) widths[j] = width_of[j];
    int thread_cases[] = {1, 4, 12};
    for (int tc = 0; tc < 3; tc++) {
        ctools_set_max_threads(thread_cases[tc]);
        eligibility();
        round_trips(1, (size_t)nrows);       /* whole range */
        round_trips(5, 3001);                /* interior range */
        round_trips(ROWS, 1);                /* one observation */
        round_trips(1, 31);                  /* below one tile */
        load_failures();
        ctools_destroy_global_pool();
    }
    return 0;
}
#endif

#if defined(CTOOLS_CASE_transport_native_1)
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
#include "ctools_transport.c"
#include "ctools_data_plan.c"
#include "ctools_data_store.c"
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
            assert(ctools_io_build_obs_map(start, end, CTOOLS_LOAD_CHECK_IF,
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
    char *a = ctools_io_copy_private_string(arena, "1234567");
    char *b = ctools_io_copy_private_string(arena, "fallback");
    char *c = ctools_io_copy_private_string(NULL, "no arena");
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
    assert(ctools_io_store_single_variable(cols + 1, 2, 2, 256ULL*1024*1024/4 + 1, NULL) == STATA_ERR_STATA_WRITE);
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
#endif

#if defined(CTOOLS_CASE_transport_scheduling_1)
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
#include "ctools_transport.c"
#include "ctools_data_plan.c"
#include "ctools_data_store.c"
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
#endif

#if defined(CTOOLS_CASE_transport_adaptive_native_1)
#include <stdio.h>
#include "ctools_config.h"
#include "ctools_arena.h"

/* Small string fixtures must exercise the same allocation retry and ownership
 * paths as large production loads. Numeric column/row gates retain real sizes. */
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
#include "ctools_transport.c"
#include "ctools_data_plan.c"
#include "ctools_data_store.c"
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
 * rows, and 256 columns with 60,000 selected rows, plus filtered equivalents.
 * Fault injection proves that the speculative tiled allocation
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
#endif

#if defined(CTOOLS_CASE_transport_store_native_1)
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
#include "ctools_transport.c"
#include "ctools_data_plan.c"
#include "ctools_data_store.c"
#undef calloc

#define ROWS 6007
#define SPARSE_ROWS 100000
static double stored[2][SPARSE_ROWS + 1];
static pthread_t caller;
static int require_caller, cancel_row;
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
    if (row == cancel_row) return 1;
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
        cancel_row = (int)map[n / 2];
        assert(ctools_store_filtered_rowpar(values,n,2,map) == STATA_ERR_CANCELLED);
        cancel_row = 0;
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
        cancel_row = SPARSE_ROWS;
        assert(ctools_store_filtered_rowpar(values,n,2,map) == STATA_ERR_CANCELLED);
        cancel_row = 0;
        map[n - 1] = SPARSE_ROWS + 1; writes = 0;
        assert(ctools_store_filtered_rowpar(values,n,2,map) == STATA_ERR_INVALID_INPUT && !writes);
    }
    free(expected);
    nrows = ROWS;
    require_caller = 0;
}

/* Both sides of the 2 MiB source-size gate must propagate native SPI r(1)
 * as cancellation. Use short output views of the full source columns so the
 * tiled/chunked store cannot accidentally classify by output count. */
static void sorted_gather_cancellation(void) {
    const size_t source_sizes[] = {31,513,262144,262145};
    nrows = ROWS; require_caller = 0;
    for (size_t c = 0; c < sizeof(source_sizes)/sizeof(*source_sizes); c++) {
        size_t source_rows = source_sizes[c];
        size_t n = source_rows < 513 ? source_rows : 513;
        double *values = malloc(source_rows * sizeof(*values));
        perm_idx_t *order = malloc(n * sizeof(*order));
        assert(values && order);
        for (size_t i = 0; i < source_rows; i++) values[i] = (double)i / 8.0 + 0.125;
        for (size_t i = 0; i < n; i++) order[i] = (perm_idx_t)(n - i - 1);
        stata_variable vars[2] = {
            {.type=STATA_TYPE_DOUBLE, .nobs=source_rows, .data.dbl=values},
            {.type=STATA_TYPE_DOUBLE, .nobs=source_rows, .data.dbl=values}
        };
        stata_data data = {.nobs=n, .nvars=2, .vars=vars, .sort_order=order};
        writes = 0;
        assert(ctools_data_store_sorted(&data,1) == STATA_OK);
        assert((size_t)writes == 2*n);
        for (size_t i = 0; i < n; i++)
            for (int j = 0; j < 2; j++)
                assert(stored[j][i+1] == values[order[i]]);
        cancel_row = 17;
        assert(ctools_data_store_sorted(&data,1) == STATA_ERR_CANCELLED);
        cancel_row = 0;
        bad_row = 17;
        assert(ctools_data_store_sorted(&data,1) == STATA_ERR_STATA_WRITE);
        bad_row = 0;
        free(order); free(values);
    }
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
        sorted_gather_cancellation();
        ctools_destroy_global_pool();
    }
    return 0;
}
#endif

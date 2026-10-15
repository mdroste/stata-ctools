
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

"""Regression checks for fused sort stores: order, ownership, and checked SPI.

Run with python3 validation/test_bigdata_native.py. Uses the existing native
test compiler helper and tests both serial and OpenMP/thread-pool execution.
"""
from pathlib import Path
import tempfile

from test_p2_native import SPI, run_test


SOURCE = SPI + r'''
#include "ctools_config.h"
#undef MIN_OBS_PER_THREAD
#define MIN_OBS_PER_THREAD 256
#include "ctools_data_io.c"
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
'''


def main():
    with tempfile.TemporaryDirectory(prefix="ctools-bigdata-native-") as tmp:
        for omp in (False, True):
            run_test(Path(tmp), "sorted_store_" + str(omp), SOURCE,
                     ["ctools_types.c", "ctools_arena.c", "ctools_threads.c"], openmp=omp)


if __name__ == "__main__":
    main()


#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "ctools_spi.h"
static ST_plugin host;
ST_plugin *_stata_ = &host;
static double *matrix;
static int checked_calls, unchecked_calls, fail_write;
static ST_retcode message(char *s) { (void)s; return 0; }
static ST_retcode save(char *s, double v) { (void)s; (void)v; return 0; }
static ST_retcode unchecked(char *s, ST_int row, ST_int col, double value) {
    (void)s;
    unchecked_calls++;
    matrix[(row - 1) * 2 + col - 1] = value;
    return 0;
}
static ST_retcode checked(char *s, ST_int row, ST_int col, double value) {
    checked_calls++;
    if (strcmp(s, "result")) return 111;
    if (row < 1 || row > 2 || col < 1 || col > 2) return 503;
    if (fail_write) return 459;
    matrix[(row - 1) * 2 + col - 1] = value;
    return 0;
}
int main(void) {
    matrix = calloc(4, sizeof(double));
    assert(matrix);
    host.matstore = unchecked;
    host.safematstore = checked;
    host.spouterr = message;
    host.scalsave = save;
    /* Under SD_FASTMODE the old wrapper overran this buffer (ASan). */
    assert(ctools_mat_store("result", 3, 1, 42) == 503);
    assert(ctools_mat_store("result", 1, 3, 42) == 503);
    assert(ctools_mat_store("result", 0, 1, 42) == 503);
    assert(ctools_mat_store("result", 1, 0, 42) == 503);
    assert(ctools_mat_store("absent", 1, 1, 42) == 111);
    fail_write = 1;
    assert(ctools_mat_store("result", 1, 1, 42) == 459);
    for (int i = 0; i < 4; i++) assert(matrix[i] == 0);
    fail_write = 0;
    assert(ctools_mat_store("result", 2, 2, 42) == 0);
    assert(matrix[3] == 42 && checked_calls == 7 && unchecked_calls == 0);
    free(matrix);
    return 0;
}

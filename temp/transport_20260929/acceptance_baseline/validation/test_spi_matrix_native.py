"""Release-mode matrix writes must stay inside host-owned buffers."""
from pathlib import Path
import os
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = r'''
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
'''


PUBLICATION = r"""
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "ctools_spi.h"
static ST_plugin host;
ST_plugin *_stata_ = &host;
static int failed_alloc, write_calls, fail_at, host_K = 3;
static double got_b[3], got_V[9];
static void *map_malloc(size_t size) { return failed_alloc ? NULL : malloc(size); }
#define malloc map_malloc
#include "civreghdfe/civreghdfe_impl.c"
#undef malloc
#include "cqreg/cqreg_regress.c"
static ST_retcode message(char *s) { (void)s; return 0; }
static ST_retcode scalar_save(char *s, double v) { (void)s; (void)v; return 0; }
static ST_retcode unsafe_matrix(char *s, ST_int row, ST_int col, double v) {
    (void)s; (void)row; (void)col; (void)v;
    assert(!"unchecked matrix callback"); return 999;
}
static ST_retcode checked_matrix(char *s, ST_int row, ST_int col, double v) {
    write_calls++;
    if (write_calls == fail_at) return 459;
    int b = !strcmp(s, "__civreghdfe_b");
    if (row < 1 || row > (b ? 1 : host_K) || col < 1 || col > host_K) return 503;
    if (b) got_b[col - 1] = v;
    else got_V[(row - 1)*host_K + col - 1] = v;
    return 0;
}
int main(void) {
    host.matstore = unsafe_matrix;
    host.safematstore = checked_matrix;
    host.spouterr = message;
    host.scalsave = scalar_save;
    /* Internal [exog, omitted exog, endog] -> posted [endog, exog, omitted]. */
    double b[] = {2, 3}, V[] = {4, 1, 1, 9};
    int omitted[] = {0, 1, 0};
    assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 0);
    const double expected_b[] = {3, 2, 0}, expected_V[] = {9,1,0, 1,4,0, 0,0,0};
    assert(!memcmp(got_b, expected_b, sizeof got_b));
    assert(!memcmp(got_V, expected_V, sizeof got_V));
    for (int f = 1; f <= 12; f++) {
        write_calls = 0; fail_at = f;
        assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 459);
        assert(write_calls == f);
    }
    fail_at = 0; failed_alloc = 1; write_calls = 0;
    assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 920);
    assert(write_calls == 0);
    failed_alloc = 0; host_K = 2;
    assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 503);
    host_K = 3;
    assert(store_iv_matrices(b, V, 2, 1, 1, omitted) == 503);
    assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 0);
    cqreg_state state = {0};
    state.K = 2; state.N = 12; state.beta = b; state.V = V;
    for (int f = 1; f <= 4; f++) {
        write_calls = 0; fail_at = f;
        assert(store_results(&state) == 459);
        assert(write_calls == f);
    }
    fail_at = 0; host_K = 1;
    assert(store_results(&state) == 503);
    host_K = 2;
    assert(store_results(&state) == 0);
    assert(!memcmp(got_V, V, sizeof V));
    return 0;
}
"""


def main():
    system = "APPLEMAC" if platform.system() == "Darwin" else "OPUNIX"
    cc = os.environ.get("CC", "/usr/bin/clang" if system == "APPLEMAC" else "cc")
    with tempfile.TemporaryDirectory(prefix="ctools-spi-matrix-") as tmp:
        source = Path(tmp) / "matrix.c"
        source.write_text(SOURCE)
        for mode in ("SD_FASTMODE", "SD_SAFEMODE"):
            binary = Path(tmp) / mode
            subprocess.run([cc, "-std=c11", "-g", "-O1", "-fsanitize=address,undefined",
                            "-fno-omit-frame-pointer", *(["-DSD_FASTMODE"] if mode == "SD_FASTMODE" else []), f"-DSYSTEM={system}",
                            "-I", str(ROOT / "src"), str(source),
                            str(ROOT / "src/ctools_spi.c"), "-o", str(binary)], check=True)
            subprocess.run([str(binary)], check=True, timeout=30)
            print(f"PASS {mode}: invalid matrix writes, injected error, and recovery")
        source.write_text(PUBLICATION)
        binary = Path(tmp) / "publication"
        link = ["-Wl,-dead_strip", "-Wl,-undefined,dynamic_lookup"] if system == "APPLEMAC" else ["-Wl,--gc-sections", "-lm"]
        subprocess.run([cc, "-std=gnu11", "-g", "-O1", "-fsanitize=address,undefined",
                        "-ffunction-sections", "-fdata-sections", "-DSD_FASTMODE", f"-DSYSTEM={system}",
                        "-I", str(ROOT / "src"), str(source), str(ROOT / "src/ctools_spi.c"),
                        *link, "-o", str(binary)], check=True)
        subprocess.run([str(binary)], check=True, timeout=30)
        print("PASS IV/quantile publication: allocation failure, 16 write failures, shapes, recovery")


if __name__ == "__main__":
    main()

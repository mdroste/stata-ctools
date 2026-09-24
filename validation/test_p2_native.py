"""P2 hash, encoding, SPI fault-injection, and compiler-selection regressions."""
import os
from pathlib import Path
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CC = os.environ.get("CC", "clang")


def run_test(directory, name, source, modules=(), openmp=False):
    cfile = directory / f"{name}.c"
    binary = directory / name
    cfile.write_text(source)
    flags = ["-std=c11", "-O2", "-DSD_FASTMODE", "-fno-fast-math", "-ffp-contract=off",
             "-fsanitize=undefined", "-fno-omit-frame-pointer", "-pthread"]
    for path in [ROOT / "src", ROOT / "src/cimport"]:
        flags += ["-I", str(path)]
    if platform.system() == "Darwin":
        flags += ["-DSYSTEM=APPLEMAC", "-D_DARWIN_C_SOURCE", "-Wl,-dead_strip", "-Wl,-undefined,dynamic_lookup"]
        if openmp:
            prefix = os.environ.get("LIBOMP_PREFIX")
            if not prefix:
                prefix = subprocess.check_output(["brew", "--prefix", "libomp"], text=True).strip()
            flags += ["-Xpreprocessor", "-fopenmp", "-I", prefix + "/include", prefix + "/lib/libomp.a"]
    else:
        flags += ["-DSYSTEM=STUNIX", "-D_GNU_SOURCE", "-ffunction-sections", "-fdata-sections", "-Wl,--gc-sections"]
        if openmp:
            flags += ["-fopenmp"]
    libraries = [] if platform.system() == "Darwin" else ["-lm"]
    subprocess.run([CC, *flags, str(cfile), *[str(ROOT / "src" / m) for m in modules], *libraries, "-o", str(binary)], check=True)
    subprocess.run([str(binary)], check=True, timeout=60)
    print(f"PASS {name}", flush=True)


SPI = r'''
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdatomic.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static int string_var, strl, bad_row, filter, binary_value, strl_length = 4, nrows = 6000;
static atomic_int reads, writes;
static ST_int obs(void) { return nrows; }
static ST_int vars(void) { return 2; }
static ST_int first(void) { return 1; }
static ST_boolean selected(ST_int i) { return !filter || (i % 2); }
static ST_boolean is_string(ST_int v) { assert(v >= 1 && v <= 2); return v == string_var; }
static ST_boolean is_strl(ST_int v) { return v == strl; }
static ST_boolean is_binary(ST_int v, ST_int i) { (void)v; (void)i; return binary_value; }
static ST_int string_length(ST_int v, ST_int i) { (void)v; (void)i; return strl_length; }
static ST_retcode getstr(ST_int v, ST_int i, char *out);
static ST_retcode get_strl(ST_int v, ST_int i, char *out, ST_int capacity) {
    assert(capacity >= strl_length + 1 && strl_length == 4);
    if (getstr(v, i, out)) return -1;
    return strl_length;
}
static ST_retcode message(char *s) { (void)s; return 0; }
static ST_retcode macro(char *name, char *buf, ST_int len) { (void)name; (void)len; buf[0] = 0; return 0; }
static ST_retcode getnum(ST_int v, ST_int i, ST_double *out) {
    assert(v >= 1 && v <= 2 && i >= 1 && i <= nrows);
    atomic_fetch_add(&reads, 1); if (i == bad_row) return 459;
    *out = i; return 0;
}
static ST_retcode putnum(ST_int v, ST_int i, ST_double value) {
    (void)value; assert(v >= 1 && v <= 2 && i >= 1 && i <= nrows);
    atomic_fetch_add(&writes, 1); return i == bad_row ? 459 : 0;
}
static ST_retcode getstr(ST_int v, ST_int i, char *out) {
    assert(v >= 1 && v <= 2 && i >= 1 && i <= nrows);
    atomic_fetch_add(&reads, 1); if (i == bad_row) return 459;
    strcpy(out, "text"); return 0;
}
static ST_retcode putstr(ST_int v, ST_int i, char *value) { return putnum(v, i, value != NULL); }
static void setup(void) {
    mock.nobs = obs; mock.nvars = vars; mock.nvar = vars;
    mock.nobs1 = first; mock.nobs2 = obs; mock.selobs = selected;
    mock.isstr = is_string; mock.isstrl = is_strl;
    mock.isbinary = is_binary; mock.sdatalen = string_length; mock.strldata = get_strl;
    mock.safevdata = getnum; mock.safestore = putnum;
    /* Raw numeric callbacks deliberately NULL: SD_FASTMODE must not bypass checks. */
    mock.sdata = getstr; mock.sstore = putstr; mock.macuse = macro;
    mock.spouterr = message; mock.spoutsml = message;
    mock.missval = 8.98846567431158e307;
}
'''


def main():
    with tempfile.TemporaryDirectory(prefix="ctools-p2-native-") as tmp:
        directory = Path(tmp)
        run_test(directory, "signed_labels", r'''
#include <assert.h>
#include "ctools_hash.c"
int main(void) {
    ctools_str_hash_table ht;
    assert(ctools_str_hash_init(&ht, 16) == 0);
    const char *names[] = {"negative", "zero", "positive"};
    int codes[] = {-10, 0, 7};
    for (int i = 0; i < 3; i++) {
        assert(ctools_str_hash_insert_value(&ht, names[i], codes[i]) == 0);
        int value = 999;
        assert(ctools_str_hash_lookup(&ht, names[i], &value) == 1 && value == codes[i]);
    }
    int value = 999;
    assert(ctools_str_hash_lookup(&ht, "absent", &value) == 0 && value == 999);
    ctools_str_hash_free(&ht);
    return 0;
}
''', ["ctools_arena.c"])
        run_test(directory, "encoding_bom", r'''
#include <assert.h>
#include "cimport/cimport_encoding.c"
int main(void) {
    const char le[] = "\xff\xfe\0\0a\0\0\0";
    const char be[] = "\0\0\xfe\xff\0\0\0a";
    CImportEncodingDetection a = cimport_detect_encoding(le, sizeof(le)-1);
    CImportEncodingDetection b = cimport_detect_encoding(be, sizeof(be)-1);
    assert(a.encoding == CIMPORT_ENC_UTF32LE && a.bom_length == 4);
    assert(b.encoding == CIMPORT_ENC_UTF32BE && b.bom_length == 4);
    return 0;
}
''')
        run_test(directory, "spi_failures", SPI + r'''
#include "ctools_config.h"
#undef MIN_OBS_PER_THREAD
#define MIN_OBS_PER_THREAD 256
#include "ctools_data_io.c"
int main(void) {
    setup();
    int indices[] = {1, 2}, widths[] = {8, 8};
    ctools_filtered_data fd;
    for (int pool = 0; pool <= 1; pool++) {
        ctools_set_max_threads(pool ? 4 : 1);
        for (int columns = 1; columns <= 2; columns++)
        for (int strings = 0; strings <= 1; strings++)
        for (int hints = 0; hints <= 1; hints++)
        for (filter = 0; filter <= 1; filter++) {
            string_var = strings ? 1 : 0;
            bad_row = 3;
            assert(ctools_data_load_ex(&fd, indices, columns, 1, nrows, 0, hints ? widths : NULL) == STATA_ERR_STATA_READ);
            assert(fd.data.vars == NULL && fd.obs_map == NULL);
            bad_row = 0;
            assert(ctools_data_load_ex(&fd, indices, columns, 1, nrows, 0, hints ? widths : NULL) == 0);
            bad_row = 3;
            assert(ctools_data_store_ex(&fd.data, indices, columns, 1) == STATA_ERR_STATA_WRITE);
            if (!strings) {
                assert(ctools_store_filtered(fd.data.vars[0].data.dbl, fd.data.nobs, 1, fd.obs_map) == STATA_ERR_STATA_WRITE);
                assert(ctools_store_filtered_rowpar(fd.data.vars[0].data.dbl, fd.data.nobs, 1, fd.obs_map) == STATA_ERR_STATA_WRITE);
            } else {
                assert(ctools_store_filtered_str(fd.data.vars[0].data.str, fd.data.nobs, 1, fd.obs_map) == STATA_ERR_STATA_WRITE);
            }
            ctools_filtered_data_free(&fd);
        }
        ctools_destroy_global_pool();
    }
    bad_row = filter = 0; string_var = 1; strl = 1; reads = 0;
    assert(ctools_data_load(&fd, indices, 1, 1, 10, 0) == 0 && reads == 10);
    ctools_filtered_data_free(&fd);
    reads = 0; strl_length = 2046;
    assert(ctools_data_load(&fd, indices, 1, 1, 10, 0) == STATA_ERR_UNSUPPORTED_TYPE && reads == 0);
    strl_length = 4; binary_value = 1;
    assert(ctools_data_load(&fd, indices, 1, 1, 10, 0) == STATA_ERR_UNSUPPORTED_TYPE && reads == 0);
    binary_value = 0; strl = 0; widths[0] = 1;
    assert(ctools_data_load_ex(&fd, indices, 1, 1, 10, 0, widths) == STATA_ERR_STATA_READ);
    indices[0] = 3;
    assert(ctools_data_load(&fd, indices, 1, 1, 10, 0) == STATA_ERR_INVALID_INPUT);
    indices[0] = 1;
    assert(ctools_data_load(&fd, indices, 1, 10, 1, 0) == STATA_ERR_INVALID_INPUT);
    assert(ctools_data_load(&fd, indices, 1, 1, nrows+1, 0) == STATA_ERR_INVALID_INPUT);
    double values[] = {1, 2}; perm_idx_t map[] = {1, 6001}; writes = 0; string_var = 0;
    assert(ctools_store_filtered_rowpar(values, 2, 1, map) == STATA_ERR_INVALID_INPUT && writes == 0);
    int64_t rows[] = {0, 2, 1}; bad_row = 3; writes = 0;
    assert(ctools_stream_var_permuted(1, rows, 3, 1) == STATA_ERR_STATA_READ && writes == 0);
    string_var = 1; writes = 0;
    assert(ctools_stream_var_permuted(1, rows, 3, 1) == STATA_ERR_STATA_READ && writes == 0);
    strl = 1; reads = 0;
    assert(ctools_stream_var_permuted(1, rows, 3, 1) == STATA_ERR_UNSUPPORTED_TYPE && reads == 0 && writes == 0);
    ctools_destroy_global_pool();
    return 0;
}
''', ["ctools_types.c", "ctools_arena.c", "ctools_threads.c"], openmp=True)

        run_test(directory, "csv_store_failures", SPI + r'''
#include "cimport/cimport_impl.c"
int main(void) {
    setup(); nrows = 3; string_var = 1; bad_row = 2;
    char *strings[] = {"a", "b", "c"};
    CImportColumnCache cache = {.count = 3, .string_data = strings};
    CImportSPIStoreTask task = {.cache = &cache, .var = 1};
    assert(cimport_spi_store_worker(&task) != NULL);
    assert(task.error == 459 && task.rows_stored == 1);
    task.var = 3; writes = 0;
    assert(cimport_spi_store_worker(&task) != NULL && writes == 0);
    string_var = 0;
    CImportParsedRow *rows[3];
    for (int i = 0; i < 3; i++) {
        rows[i] = calloc(1, sizeof(CImportParsedRow) + sizeof(CImportFieldRef));
        rows[i]->num_fields = 1;
        rows[i]->fields[0].offset = i;
        rows[i]->fields[0].length = 1;
    }
    CImportParsedChunk chunk = {.rows = rows, .num_rows = 3};
    CImportContext ctx = {.chunks = &chunk, .num_chunks = 1, .total_rows = 3,
                          .file_data = "123", .decimal_separator = '.'};
    task.var = 1; task.ctx = &ctx; task.col_idx = 0;
    assert(cimport_spi_store_numeric_direct(&task) != NULL);
    assert(task.error == 459 && task.rows_stored == 1);
    for (int i = 0; i < 3; i++) free(rows[i]);
    return 0;
}
''', ["ctools_types.c"])
        run_test(directory, "xlsx_store_failures", SPI + r'''
#include "cimport/cimport_xlsx.c"
int main(void) {
    setup(); nrows = 3; string_var = 1; bad_row = 2;
    char *strings[] = {"a", "b", "c"}; double values[] = {1, 2, 3};
    CImportColumnInfo col = {.type = CIMPORT_COL_STRING};
    CImportColumnCache cache = {.count = 3, .string_data = strings, .numeric_data = values};
    xlsx_store_task task = {.col = &col, .cache = &cache, .var = 1, .stata_nobs = 3};
    assert(xlsx_store_worker(&task) != NULL && task.error == 459);
    string_var = 0; col.type = CIMPORT_COL_NUMERIC;
    assert(xlsx_store_worker(&task) != NULL && task.error == 459);
    task.stata_nobs = 2; writes = 0;
    assert(xlsx_store_worker(&task) != NULL && task.error == 198 && writes == 0);
    return 0;
}
''')
        run_test(directory, "matrix_store_failures", SPI + r'''
#include "cbinscatter/cbinscatter_impl.c"
static ST_retcode scalars(char *s, ST_double d) { (void)s; (void)d; return 0; }
static ST_retcode matrix(char *s, ST_int row, ST_int col, ST_double d) {
    (void)s; (void)row; (void)d; return col == 4 ? 503 : 0;
}
int main(void) {
    setup(); mock.scalsave = scalars; mock.safematstore = matrix;
    BinStats bin = {.n_obs = 3};
    ByGroupResult group = {.num_bins = 1, .bins = &bin};
    BinscatterResults results = {.num_by_groups = 1, .groups = &group};
    BinscatterConfig config = {0};
    assert(store_results(&results, &config) == 503);
    return 0;
}
''')


if __name__ == "__main__":
    main()

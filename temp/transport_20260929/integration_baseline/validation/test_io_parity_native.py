"""UBSan checks for CSV formatting and parsing boundaries, independent of Stata."""
from pathlib import Path
import tempfile
from test_p2_native import run_test

SOURCE = r'''
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static ST_boolean missing(double x) { return x >= mock.missval; }
#include "cexport/cexport_format.c"
#include "cimport/cimport_parse.c"
int main(void) {
    uint64_t bits = UINT64_C(0x7fe0000000000000);
    memcpy(&mock.missval, &bits, 8);
    mock.ismissing = missing;
    char buffer[8192], source[2046];
    for (int i=0; i<2045; i++) source[i] = i % 2 ? 'x' : '"';
    source[2045] = 0;
    int n = cexport_write_quoted_string(source, buffer, sizeof(buffer));
    assert(n == 2045 + 1023 + 2);
    assert(buffer[0] == '"' && buffer[n-1] == '"');
    assert(cexport_write_quoted_string(source, buffer, 256) == -1);
    for (int i=1; i<=26; i++) {
        uint64_t encoded = bits + ((uint64_t)i << 40);
        double x, parsed;
        memcpy(&x, &encoded, 8);
        assert(cexport_double_to_str(x, buffer, sizeof(buffer), false, VARTYPE_DOUBLE) == 2);
        assert(buffer[0] == '.' && buffer[1] == 'a'+i-1);
        assert(cimport_field_looks_numeric_sep(buffer, 2, '.', 0));
        assert(cimport_parse_number(buffer, 2, &parsed, mock.missval, '.', 0));
        assert(parsed == x);
    }
    assert(!cimport_field_looks_numeric_sep("NA", 2, '.', 0));
    assert(!cimport_field_looks_numeric_sep("NaN", 3, '.', 0));
    double number;
    assert(cimport_parse_number("\"1.25\"", 6, &number, mock.missval, '.', 0) && number == 1.25);
    const char *raw = "1e-2";
    CImportFieldRef field = {0, 4};
    bool integer;
    assert(cimport_analyze_numeric_with_sep(raw, &field, '"', '.', 0, &number, &integer));
    assert(!integer && number == 0.01);
    assert(cexport_double_to_str(1e100, buffer, sizeof(buffer), false, VARTYPE_DOUBLE) > 0);
    assert(!strcmp(buffer,"1.0000000000e+100"));
    assert(cexport_double_to_str(1e-5, buffer, sizeof(buffer), false, VARTYPE_DOUBLE) > 0);
    assert(!strcmp(buffer,".00001"));
    assert(cexport_double_to_str(-27915931.612694785, buffer, sizeof(buffer), false, VARTYPE_DOUBLE) > 0);
    assert(!strcmp(buffer, "-27915931.61269479"));
    assert(cexport_double_to_str(136948012240.34935, buffer, sizeof(buffer), false, VARTYPE_DOUBLE) > 0);
    assert(!strcmp(buffer, "136948012240.3494"));
    /* A long quoted field must also survive the mixed-type row writer. */
    stata_variable var = {0}; char *strings[] = {source};
    var.type = STATA_TYPE_STRING; var.data.str = strings;
    vartype_t type = VARTYPE_STRING;
    cexport_context ctx = {0};
    ctx.filtered.data.nvars = 1; ctx.filtered.data.vars = &var;
    ctx.vartypes = &type; ctx.delimiter = ','; ctx.quote_if_needed = true;
    ctx.line_ending[0] = '\n'; ctx.line_ending_len = 1;
    assert(cexport_format_row(&ctx, 0, buffer, sizeof(buffer)) == n + 1);
    assert(buffer[n] == '\n');
    assert(cexport_format_row(&ctx, 0, buffer, 256) == -1);
    return 0;
}
'''

if __name__ == '__main__':
    with tempfile.TemporaryDirectory(prefix='ctools-io-parity-') as directory:
        run_test(Path(directory), 'io_parity', SOURCE, modules=('ctools_types.c',))

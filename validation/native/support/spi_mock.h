
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
    /* Transport reads use the unchecked accessor after validating variable
     * indices and observation ranges up front; writes stay checked (the real
     * unchecked store corrupts columns under concurrent writes). Provide
     * both pairs exactly like Stata does. */
    mock.vdata = getnum; mock.store = putnum;
    mock.sdata = getstr; mock.sstore = putstr; mock.macuse = macro;
    mock.spouterr = message; mock.spoutsml = message;
    mock.missval = 8.98846567431158e307;
}

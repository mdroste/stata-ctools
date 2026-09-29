
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

#include <stdio.h>
#include <sys/resource.h>
#include <mach/mach.h>
#include <time.h>
#include "ctools_data_io.c"
static int ncols=4,content_len=8;
static ST_int wide_vars(void){return ncols;}
static ST_boolean text_col(ST_int v){(void)v;return 1;}
static ST_boolean no_strl(ST_int v){(void)v;return 0;}
static ST_retcode text(ST_int v,ST_int row,char *out){(void)v;memset(out,'a'+row%26,content_len);out[content_len]=0;return 0;}
static ST_retcode width_macro(char *name,char *out,ST_int size){(void)name;(void)size;strcpy(out,"2045,2045,2045,2045");return 0;}
static double now(void){struct timespec ts;clock_gettime(CLOCK_MONOTONIC,&ts);return ts.tv_sec+1e-9*ts.tv_nsec;}
int main(int argc,char **argv){
 setup();nrows=100000;content_len=argc>1?atoi(argv[1]):8;
 mock.nvars=mock.nvar=wide_vars;mock.isstr=text_col;mock.isstrl=no_strl;mock.sdata=text;mock.macuse=width_macro;
 ctools_set_max_threads(4);ctools_filtered_data fd;double start=now();
 assert(ctools_data_load(&fd,NULL,0,1,nrows,CTOOLS_LOAD_SKIP_IF)==0);
 double elapsed=now()-start;struct rusage usage;getrusage(RUSAGE_SELF,&usage);
 mach_task_basic_info_data_t mem;mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;
 assert(task_info(mach_task_self(),MACH_TASK_BASIC_INFO,(task_info_t)&mem,&count)==KERN_SUCCESS);
 for(size_t j=0;j<fd.data.nvars;j++)for(size_t i=0;i<fd.data.nobs;i++) {
    assert(strlen(fd.data.vars[j].data.str[i])==(size_t)content_len);
    if(content_len)assert(fd.data.vars[j].data.str[i][0]=='a'+(i+1)%26);
 }
 start=now();ctools_filtered_data_free(&fd);double cleanup=now()-start;
 printf("{\"length\":%d,\"load_s\":%.6f,\"cleanup_s\":%.6f,\"maxrss_bytes\":%ld,\"resident_bytes\":%llu,\"virtual_bytes\":%llu}\n",content_len,elapsed,cleanup,usage.ru_maxrss,(unsigned long long)mem.resident_size,(unsigned long long)mem.virtual_size);
 ctools_destroy_global_pool();return 0;
}

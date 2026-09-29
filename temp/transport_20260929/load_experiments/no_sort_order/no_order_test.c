
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

#include "ctools_config.h"
#undef MIN_OBS_PER_THREAD
#define MIN_OBS_PER_THREAD 256
static int fail_flat;
static void *flat_calloc(size_t n,size_t size){
    if(fail_flat && size==65 && n>512){fail_flat=0;return NULL;}
    return calloc(n,size);
}
#define calloc flat_calloc
#include "ctools_data_io.c"
#undef calloc
static ST_boolean select_none(ST_int row){(void)row;return 0;}
static ST_int four_vars(void){return 4;}
static ST_boolean all_text(ST_int v){(void)v;return 1;}
static ST_retcode four_text(ST_int v,ST_int row,char *out){
    assert(v>=1 && v<=4 && row>=1 && row<=nrows);strcpy(out,"correct");return 0;
}
int main(void){
 setup();ctools_set_max_threads(4);nrows=6000;
 for(string_var=0;string_var<=2;string_var+=2)for(filter=0;filter<=1;filter++)for(int opt=0;opt<=1;opt++){
  ctools_filtered_data fd;int widths[]={0,8};
  assert(ctools_data_load_ex(&fd,NULL,0,3,nrows-2,opt?CTOOLS_LOAD_NO_SORT_ORDER:0,widths)==STATA_OK);
  assert(fd.data.nobs>0);
  if(opt)assert(!fd.data.sort_order);else for(size_t i=0;i<fd.data.nobs;i++)assert(fd.data.sort_order[i]==i);
  size_t out=0;for(int row=3;row<=nrows-2;row++)if(selected(row)){
   assert(fd.obs_map[out]==(perm_idx_t)row);
   assert(fd.data.vars[0].data.dbl[out]==row);
   if(string_var)assert(!strcmp(fd.data.vars[1].data.str[out],"text"));
   else assert(fd.data.vars[1].data.dbl[out]==row);out++;
  }
  assert(out==fd.data.nobs);writes=0;assert(ctools_data_store(&fd.data,1)==STATA_OK);assert(writes==2*out);
  writes=0;assert(ctools_data_store_sorted(&fd.data,1)==(opt?STATA_ERR_INVALID_INPUT:STATA_OK));assert(writes==(opt?0:2*out));
  ctools_filtered_data_free(&fd);ctools_filtered_data_free(&fd);assert(!fd.data.vars && !fd.data.sort_order && !fd.obs_map);
 }
 mock.selobs=select_none;
 ctools_filtered_data fd;
 assert(ctools_data_load(&fd,NULL,0,1,nrows,CTOOLS_LOAD_NO_SORT_ORDER)==STATA_OK);
 assert(!fd.data.nobs && !fd.data.sort_order);ctools_filtered_data_free(&fd);
 nrows=0;assert(ctools_data_load(&fd,NULL,0,0,0,CTOOLS_LOAD_NO_SORT_ORDER)==STATA_OK);
 assert(!fd.data.nobs && !fd.data.sort_order);ctools_filtered_data_free(&fd);
 mock.selobs=selected;filter=0;nrows=6000;mock.nvars=mock.nvar=four_vars;mock.isstr=all_text;mock.sdata=four_text;
 int widths[]={64,64,64,64};fail_flat=1;
 assert(ctools_data_load_ex(&fd,NULL,0,1,nrows,CTOOLS_LOAD_SKIP_IF|CTOOLS_LOAD_NO_SORT_ORDER,widths)==STATA_OK);
 assert(!fd.data.sort_order && !fail_flat);
 for(int j=0;j<4;j++)for(int i=0;i<nrows;i++)assert(!strcmp(fd.data.vars[j].data.str[i],"correct"));
 ctools_filtered_data_free(&fd);ctools_destroy_global_pool();return 0;
}

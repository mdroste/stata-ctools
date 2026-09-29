
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
#include "ctools_config.h"
#undef MIN_OBS_PER_THREAD
#define MIN_OBS_PER_THREAD 256
static int fail_allocation;
static void *test_alloc2(size_t count,size_t size) {
    if (count>512 && size==sizeof(double) && fail_allocation>0 && --fail_allocation==0) return NULL;
    return ctools_safe_aligned_alloc2(CACHE_LINE_SIZE,count,size);
}
#undef ctools_safe_cacheline_alloc2
#define ctools_safe_cacheline_alloc2(count,size) test_alloc2((count),(size))
#include "ctools_arena.h"
static atomic_int fail_large_arena, failed_large_seen, retry_seen;
static size_t expected_retry;
static ctools_string_arena *test_arena_create(size_t capacity,ctools_string_arena_mode mode) {
    if(capacity==expected_retry) atomic_fetch_add(&retry_seen,1);
    if(capacity>1000000 && atomic_exchange(&fail_large_arena,0)) {
        atomic_fetch_add(&failed_large_seen,1);return NULL;
    }
    return ctools_string_arena_create(capacity,mode);
}
#define ctools_string_arena_create test_arena_create
#include "ctools_data_io.c"
#undef ctools_string_arena_create
static int ncols, shape;
static ST_int wide_vars(void) { return ncols; }
static ST_boolean text_col(ST_int v) { return shape==2 || (shape==1 && v%2==0); }
static ST_boolean no_strl(ST_int v) { (void)v; return 0; }
static ST_retcode wide_num(ST_int v,ST_int row,double *out) {
    assert(v>=1 && v<=ncols && row>=1 && row<=nrows);
    if(row==bad_row) return 459;
    *out=(double)row+0.125*v; return 0;
}
static ST_retcode wide_text(ST_int v,ST_int row,char *out) {
    assert(v>=1 && v<=ncols && row>=1 && row<=nrows);
    if(row==bad_row) return 459;
    if(row%17==0) { out[0]=0; return 0; }
    memset(out,'a'+(row+v)%26,2045);
    out[0]=(char)0xc3;out[1]=(char)0xa9;out[2045]=0;return 0;
}
static ST_retcode wide_macro(char *name,char *out,ST_int capacity) {
    (void)name;
    size_t used=0;
    for(int j=1;j<=ncols;j++) {
        int bytes=snprintf(out+used,(size_t)capacity-used,"%s%d",j>1?",":"",text_col(j)?2045:0);
        assert(bytes>0 && used+(size_t)bytes<(size_t)capacity);used+=(size_t)bytes;
    }
    return 0;
}
int main(void) {
    setup();nrows=1109;
    mock.nvars=mock.nvar=wide_vars;mock.isstr=text_col;mock.isstrl=no_strl;
    mock.vdata=mock.safevdata=wide_num;mock.sdata=wide_text;mock.macuse=wide_macro;
    int cases[]={1,2,3,4,20,64};int thread_cases[]={1,4,12};
    for(int t=0;t<3;t++) {
      ctools_set_max_threads(thread_cases[t]);
      for(int c=0;c<6;c++) {ncols=cases[c];
        for(shape=0;shape<=2;shape++) for(filter=0;filter<=1;filter++) {
          ctools_filtered_data fd;
          assert(ctools_data_load(&fd,NULL,0,3,nrows-2,CTOOLS_LOAD_CHECK_IF)==STATA_OK);
          size_t out=0;
          for(int row=3;row<=nrows-2;row++) if(selected(row)) {
            assert(fd.obs_map[out]==(perm_idx_t)row && fd.data.sort_order[out]==out);
            for(int j=1;j<=ncols;j++) {
              if(text_col(j)) {char value[2046];wide_text(j,row,value);assert(!strcmp(fd.data.vars[j-1].data.str[out],value));}
              else assert(fd.data.vars[j-1].data.dbl[out]==row+0.125*j);
            }
            out++;
          }
          assert(fd.data.nobs==out);ctools_filtered_data_free(&fd);
          bad_row=777;
          assert(ctools_data_load(&fd,NULL,0,3,nrows-2,CTOOLS_LOAD_CHECK_IF)==STATA_ERR_STATA_READ);
          ctools_filtered_data_free(&fd);bad_row=0;
        }
      }
      ctools_destroy_global_pool();
    }
    ncols=64;shape=0;filter=0;
    ctools_filtered_data fd;int invalid=65;
    assert(ctools_data_load(&fd,&invalid,1,1,nrows,0)==STATA_ERR_INVALID_INPUT);
    assert(ctools_data_load(&fd,NULL,0,1,nrows+1,0)==STATA_ERR_INVALID_INPUT);
    assert(ctools_data_load(&fd,NULL,0,99,33,0)==STATA_ERR_INVALID_INPUT);
    for(int fail=1;fail<=3;fail++) {
      fail_allocation=fail;
      assert(ctools_data_load(&fd,NULL,0,3,nrows-2,CTOOLS_LOAD_SKIP_IF)==STATA_OK);
      for(size_t row=0;row<fd.data.nobs;row++) for(int j=1;j<=ncols;j++)
        assert(fd.data.vars[j-1].data.dbl[row]==row+3+0.125*j);
      ctools_filtered_data_free(&fd);
    }
    ctools_destroy_global_pool();
    shape=2;
    for(ncols=1;ncols<=4;ncols+=3) for(filter=0;filter<=1;filter++) {
      size_t selected_count=0;
      for(int row=3;row<=nrows-2;row++) selected_count+=selected(row)!=0;
      expected_retry=selected_count*64;
      fail_large_arena=1;failed_large_seen=0;retry_seen=0;
      assert(ctools_data_load(&fd,NULL,0,3,nrows-2,CTOOLS_LOAD_CHECK_IF)==STATA_OK);
      assert(failed_large_seen==1 && retry_seen==1);
      size_t out=0;
      for(int row=3;row<=nrows-2;row++) if(selected(row)) {
        char value[2046];
        for(int j=1;j<=ncols;j++) {
          wide_text(j,row,value);
          assert(!strcmp(fd.data.vars[j-1].data.str[out],value));
        }
        out++;
      }
      assert(out==fd.data.nobs);ctools_filtered_data_free(&fd);
    }
    ctools_destroy_global_pool();
    return 0;
}

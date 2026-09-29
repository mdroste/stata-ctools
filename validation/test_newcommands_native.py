"""UBSan and allocation-failure checks for the new command kernels.

The mock loader owns real heap copies; every allocation in each command is
failed in turn. This checks cleanup of staged caches and recovery on later calls.
Stata reference tests separately exercise the real SPI, syntax, and metadata.
Set CTOOLS_SANITIZERS=address,undefined to also run ASan (enabled in Linux CI).
"""
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
#include <stdio.h>
#include <math.h>
#include "stplugin.h"
#include "ctools_types.h"
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static int nr, nv, types[16], calls, fail_at, alive, write_error, load_error;
static double numbers[16][16];
static char strings[16][16][2046], widths[16000], nout[64];
static void *owned[1000];
static void track(void *p) { if(p) {for(int i=0;i<1000;i++)if(!owned[i]){owned[i]=p;alive++;return;}abort();} }
static void *test_malloc(size_t n) {if(++calls==fail_at)return NULL;void *p=malloc(n);track(p);return p;}
static void *test_calloc(size_t n,size_t s) {if(++calls==fail_at)return NULL;void *p=calloc(n,s);track(p);return p;}
static void test_free(void *p) {if(!p)return;for(int i=0;i<1000;i++)if(owned[i]==p){owned[i]=NULL;alive--;free(p);return;}assert(!"unowned free");}
static ST_int nobs(void){return nr;}
static ST_int nvars(void){return nv;}
static ST_boolean isstr(ST_int v){assert(v>=1&&v<=nv);return types[v-1];}
static ST_boolean isstrl(ST_int v){(void)v;return 0;}
static ST_retcode putnum(ST_int v,ST_int r,ST_double x){assert(v>=1&&v<=nv&&r>=1&&r<=nr);if(write_error)return 459;numbers[v-1][r-1]=x;return 0;}
static ST_retcode putstr(ST_int v,ST_int r,char *x){assert(v>=1&&v<=nv&&r>=1&&r<=nr);assert(strlen(x)<=2045);if(write_error)return 459;strcpy(strings[v-1][r-1],x);return 0;}
static ST_retcode macro_save(char *key,char *value){if(!strcmp(key,"_csplit_widths"))strcpy(widths,value);else if(!strcmp(key,"_crangejoin_nout"))strcpy(nout,value);else abort();return 0;}
static ST_retcode macro_get(char *key,char *value,ST_int size){const char *s=!strcmp(key,"_csplit_delim1")?"::":":";assert(size>=(int)strlen(s)+1);strcpy(value,s);return 0;}
int ctools_get_max_threads(void){return 1;}
double ctools_timer_seconds(void){return 0;}
void ctools_verbose(const char *m,int v,const char *f,...){(void)m;(void)v;(void)f;}
void ctools_error(const char *m,const char *f,...){(void)m;(void)f;}
#define malloc test_malloc
#define calloc test_calloc
#define free test_free
void ctools_filtered_data_init(ctools_filtered_data *fd){memset(fd,0,sizeof(*fd));}
void ctools_filtered_data_free(ctools_filtered_data *fd){
 for(size_t k=0;k<fd->data.nvars;k++){
  stata_variable *v=fd->data.vars+k;
  if(v->type==STATA_TYPE_STRING){if(v->data.str)for(size_t i=0;i<v->nobs;i++)free(v->data.str[i]);free(v->data.str);}
  else free(v->data.dbl);
 }
 free(fd->data.vars);free(fd->data.sort_order);free(fd->obs_map);memset(fd,0,sizeof(*fd));
}
stata_retcode ctools_data_load(ctools_filtered_data *fd,int *vars,size_t n,size_t start,size_t end,int flags){
 (void)start;(void)end;(void)flags;if(load_error)return load_error;if(!vars)n=nv;
 fd->data.vars=calloc(n,sizeof(stata_variable));if(!fd->data.vars)return STATA_ERR_MEMORY;
 fd->data.nvars=n;fd->data.nobs=nr;
 fd->data.sort_order=malloc(nr*sizeof(perm_idx_t));fd->obs_map=malloc(nr*sizeof(perm_idx_t));
 if(!fd->data.sort_order||!fd->obs_map)return STATA_ERR_MEMORY;
 for(int i=0;i<nr;i++){fd->data.sort_order[i]=i;fd->obs_map[i]=i+1;}
 for(size_t k=0;k<n;k++){
  int source=vars?vars[k]-1:(int)k;stata_variable *v=fd->data.vars+k;v->nobs=nr;v->type=types[source]?STATA_TYPE_STRING:STATA_TYPE_DOUBLE;
  if(types[source]){
   v->data.str=calloc(nr,sizeof(char*));if(!v->data.str)return STATA_ERR_MEMORY;
   for(int i=0;i<nr;i++){size_t len=strlen(strings[source][i])+1;v->data.str[i]=malloc(len);if(!v->data.str[i])return STATA_ERR_MEMORY;memcpy(v->data.str[i],strings[source][i],len);}
  }else{v->data.dbl=malloc(nr*sizeof(double));if(!v->data.dbl)return STATA_ERR_MEMORY;memcpy(v->data.dbl,numbers[source],nr*sizeof(double));}
 }
 return STATA_OK;
}
stata_retcode ctools_store_filtered(double *v,size_t n,int col,perm_idx_t *map){for(size_t i=0;i<n;i++)if(putnum(col,map[i],v[i]))return STATA_ERR_STATA_WRITE;return STATA_OK;}
#include "ctools_order.c"
#include "cipolate/cipolate_impl.c"
#define state split_state
#include "csplit/csplit_impl.c"
#undef state
#define state join_state
#include "crangejoin/crangejoin_impl.c"
#undef state

static void reset(void){assert(alive==0);memset(numbers,0,sizeof(numbers));memset(strings,0,sizeof(strings));memset(types,0,sizeof(types));nr=nv=0;calls=0;write_error=0;}
static int interpolation(void){
 reset();nr=4;nv=3;numbers[0][0]=numbers[0][2]=mock.missval;numbers[0][1]=10;numbers[0][3]=30;
 for(int i=0;i<nr;i++)numbers[1][i]=i;
 int rc=cipolate_main("0 1 0");if(!rc)for(int i=0;i<nr;i++)assert(numbers[2][i]==10*i);
 assert(alive==0);return rc;
}
static int splitting(void){
 reset();nr=2;nv=1;types[0]=1;strcpy(strings[0][0],"a::b:c");strcpy(strings[0][1],"é:東京");
 int rc=csplit_main("scan 2 2046 1 0");
 if(!rc){assert(!strcmp(widths,"2 6 1 1 "));nv=4;for(int k=0;k<nv;k++)types[k]=1;
  rc=csplit_main("write");if(!rc){assert(!strcmp(strings[0][0],"a"));assert(!strcmp(strings[1][0],""));assert(!strcmp(strings[2][0],"b"));assert(!strcmp(strings[3][0],"c"));}}
 csplit_cleanup_cache();assert(alive==0);return rc;
}
static int joining(void){
 reset();nr=3;nv=2;numbers[0][0]=3;numbers[0][1]=1;numbers[0][2]=1;types[1]=1;
 strcpy(strings[1][0],"three");strcpy(strings[1][1],"first");strcpy(strings[1][2],"second");
 int rc=crangejoin_main("using 0 1 0");
 if(!rc){nr=2;nv=5;memset(types,0,sizeof(types));types[4]=1;numbers[0][0]=100;numbers[0][1]=200;numbers[1][0]=1;numbers[1][1]=9;numbers[2][0]=3;numbers[2][1]=10;
  rc=crangejoin_main("prepare 1");
  if(!rc){assert(!strcmp(nout,"4"));nr=4;rc=crangejoin_main("write");
   if(!rc){assert(numbers[0][0]==100&&numbers[0][2]==100&&numbers[0][3]==200);assert(numbers[3][0]==1&&numbers[3][2]==3&&numbers[3][3]==mock.missval);assert(!strcmp(strings[4][0],"first")&&!strcmp(strings[4][1],"second")&&!strcmp(strings[4][3],""));}}}
 crangejoin_cleanup_cache();assert(alive==0);return rc;
}
int main(void){
 mock.nobs=nobs;mock.nvars=nvars;mock.nvar=nvars;mock.isstr=isstr;mock.isstrl=isstrl;
 mock.store=putnum;mock.safestore=putnum;mock.sstore=putstr;mock.macresave=macro_save;mock.macuse=macro_get;mock.missval=8.98846567431158e307;
 int (*tests[])(void)={interpolation,splitting,joining};
 for(int k=0;k<3;k++){fail_at=0;assert(tests[k]()==0);int allocations=calls;
  for(int f=1;f<=allocations;f++){fail_at=f;assert(tests[k]()!=0);}
  fail_at=0;assert(tests[k]()==0);
 }
 assert(csplit_main("write")==198);assert(crangejoin_main("write")==198);assert(alive==0);
 for(int e=STATA_ERR_MEMORY;e<=STATA_ERR_CANCELLED;e++) {
  reset();nr=4;nv=3;load_error=e;fail_at=0;
  assert(cipolate_main("0 1 0")==ctools_stata_rc(e));assert(alive==0);
  for(int i=0;i<nr;i++)assert(numbers[2][i]==0);
 }
 load_error=0;
 /* Exact missing-key grouping and stable ties, independent of Stata's sorter. */
 reset();nr=4;nv=1;numbers[0][0]=mock.missval;numbers[0][1]=nextafter(mock.missval,INFINITY);numbers[0][2]=mock.missval;numbers[0][3]=2;
 ctools_filtered_data fd={0};assert(ctools_data_load(&fd,NULL,0,0,0,0)==0);int keys[]={0};assert(ctools_order_stable(&fd.data,keys,1)==0);
 assert(fd.data.sort_order[0]==3&&fd.data.sort_order[1]==0&&fd.data.sort_order[2]==2&&fd.data.sort_order[3]==1);ctools_filtered_data_free(&fd);assert(alive==0);
 puts("PASS new command kernels, allocation failures, stale-cache rejection, exact key ordering");return 0;
}
'''


def main():
    with tempfile.TemporaryDirectory(prefix='ctools-new-native-') as tmp:
        cfile = Path(tmp) / 'test.c'
        binary = Path(tmp) / 'test'
        cfile.write_text(SOURCE)
        flags = ['-std=c11', '-O1', '-g', '-fsanitize=' + os.environ.get('CTOOLS_SANITIZERS', 'undefined'),
                 '-fno-omit-frame-pointer', '-fno-fast-math', '-ffp-contract=off',
                 '-Wno-unknown-pragmas', '-I', str(ROOT / 'src')]
        flags += ['-DSYSTEM=APPLEMAC', '-D_DARWIN_C_SOURCE'] if platform.system() == 'Darwin' else ['-DSYSTEM=STUNIX', '-D_GNU_SOURCE']
        subprocess.run([os.environ.get('CC','clang'), *flags, str(cfile),
                        str(ROOT/'src/ctools_parse.c'), '-lm', '-o', str(binary)], check=True)
        subprocess.run([str(binary)], check=True, timeout=60)

if __name__ == '__main__':
    main()

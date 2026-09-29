"""Regression contracts for parsing, command sessions, ownership, and statuses."""
from pathlib import Path
import os
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = r'''
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "stplugin.h"
#include "ctools_config.h"
#include "ctools_parse.h"
#include "ctools_runtime.h"
#include "ctools_types.h"
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static int freed[32];
void cmerge_cleanup_cache(void) { freed[CTOOLS_CMD_CMERGE]++; }
void cimport_cleanup_cache(void) { freed[CTOOLS_CMD_CIMPORT]++; }
void cio_cleanup(void) { freed[CTOOLS_CMD_CIO]++; }
void csplit_cleanup_cache(void) { freed[CTOOLS_CMD_CSPLIT]++; }
void crangejoin_cleanup_cache(void) { freed[CTOOLS_CMD_CRANGEJOIN]++; }
void creghdfe_cleanup_state(void) { freed[CTOOLS_CMD_CREGHDFE]++; }
void civreghdfe_cleanup_state(void) { freed[CTOOLS_CMD_CIVREGHDFE]++; }
void cpplmhdfe_cleanup_state(void) { freed[CTOOLS_CMD_CPPLMHDFE]++; }
void cexport_cleanup_state(void) { freed[CTOOLS_CMD_CEXPORT]++; }
static const ctools_command_descriptor merge={"cmerge",CTOOLS_CMD_CMERGE,NULL,cmerge_cleanup_cache};
static const ctools_command_descriptor import={"cimport",CTOOLS_CMD_CIMPORT,NULL,cimport_cleanup_cache};
static const ctools_command_descriptor io={"cio",CTOOLS_CMD_CIO,NULL,cio_cleanup};
static const ctools_command_descriptor join={"crangejoin",CTOOLS_CMD_CRANGEJOIN,NULL,crangejoin_cleanup_cache};
static void invoke(const ctools_command_descriptor *d,const char *phase,int rc) {
 assert(ctools_command_begin(d,phase)==0); ctools_command_finish(rc);
}
static void lifecycle(void) {
 invoke(&merge,"load_using",0); invoke(&merge,"execute",0);
 assert(freed[CTOOLS_CMD_CMERGE]==1);
 invoke(&merge,"load_using",0); invoke(&import,"scan",0); invoke(&merge,"load_using",0);
 assert(freed[CTOOLS_CMD_CIMPORT]==1); // Regression: A -> B -> A.
 invoke(&import,"scan",0); invoke(&io,"scan",0);
 assert(freed[CTOOLS_CMD_CIMPORT]==2); // A -> B -> C, distinct CSV and generic IO.
 invoke(&io,"column 1",0); invoke(&io,"load",0); invoke(&io,"blob 1 1 0",459);
 assert(freed[CTOOLS_CMD_CIO]==1);
 assert(ctools_command_begin(&io,"column 1")==198);
 invoke(&join,"using 1",0); invoke(&join,"prepare 1",920);
 assert(ctools_command_begin(&join,"write")==198);
 invoke(&join,"using 1",0); invoke(&join,"prepare 1",0); invoke(&join,"write",0);
 invoke(&import,"scan",459); assert(ctools_command_begin(&import,"blob 1")==198);
 invoke(&import,"scan",0); invoke(&import,"scan",0); // Interrupted/restarted scan frees old state.
 int count=freed[CTOOLS_CMD_CIMPORT];
 ctools_command_reset(); ctools_command_reset(); assert(freed[CTOOLS_CMD_CIMPORT]==count+1);
 invoke(&import,"scan",0); assert(ctools_command_begin(NULL,"")==198);
 assert(freed[CTOOLS_CMD_CIMPORT]==count+2);
}
static void parsing(void) {
 char b[32]="untouched";int n=47;uint64_t u=47;double d=47;
 assert(ctools_parse_bool_option("notverbose verbose","verbose"));
 assert(!ctools_parse_bool_option("notverbose","verbose"));
 assert(ctools_parse_string_option("nolabel=wrong label=right","label",b,sizeof b)==1&&!strcmp(b,"right"));
 assert(ctools_parse_string_option("label=\"with spaces\"","label",b,sizeof b)==1&&!strcmp(b,"with spaces"));
 assert(ctools_parse_string_option("label=\"unterminated","label",b,sizeof b)==-1);
 assert(ctools_parse_string_option("label=x label=y","label",b,sizeof b)==-1);
 assert(ctools_parse_string_option("label=1234","label",b,4)==-1);
 assert(ctools_parse_int_option("n=2147483647","n",&n)==1&&n==INT_MAX);
 assert(ctools_parse_int_option("n=-2147483648","n",&n)==1&&n==INT_MIN);
 for(size_t i=0;i<6;i++){
  const char *bad[]={"n=2147483648","n=-2147483649","n=1x","n=","n","n=1 n=2"};
  n=47;assert(ctools_parse_int_option(bad[i],"n",&n)==-1&&n==47);
 }
 assert(ctools_parse_int_option("other=3","n",&n)==0&&n==47);
 assert(ctools_parse_u64_option("n=18446744073709551615","n",&u)==1&&u==UINT64_MAX);
 assert(ctools_parse_u64_option("n=18446744073709551616","n",&u)==-1);
 assert(ctools_parse_u64_option("n=-1","n",&u)==-1);
 assert(ctools_parse_double_checked("d=nan","d",&d)==-1&&d==47);
 assert(ctools_parse_double_checked("d=1e999","d",&d)==-1);
 assert(ctools_parse_double_checked("d=1.5junk","d",&d)==-1);
 assert(ctools_parse_seed_option("seedhi=2147483647 seedlo=2147483647",&u)==1&&u==((UINT64_C(1)<<62)-1));
 assert(ctools_parse_seed_option("seedhi=1",&u)==-1);
 assert(ctools_parse_seed_option("seedhi=2147483648 seedlo=0",&u)==-1);
 const char *p=" \t2147483648",*saved=p;assert(ctools_parse_next_int(&p,&n)==-1&&p==saved);
 const char bounded[]={'a','=','1',' ','b'};ctools_arg_cursor c;ctools_arg_token t;
 ctools_args_init(&c,bounded,3);assert(ctools_args_next(&c,&t)==1&&t.length==3);assert(ctools_args_next(&c,&t)==0);
}
/* Model Windows allocator separation on every host: a wrong deallocator aborts.
 * Real Windows builds use _aligned_free through the same explicit storage mode. */
static void *aligned[32];static int live;
static void *alloc_aligned(size_t n) {void *p=ctools_cacheline_alloc(n);assert(p);aligned[live++]=p;return p;}
static void release_aligned(void *p) {if(!p)return;int found=0;for(int i=0;i<live;i++)if(aligned[i]==p){aligned[i]=aligned[--live];found=1;break;}assert(found);ctools_aligned_free(p);}
static void release_heap(void *p) {for(int i=0;i<live;i++)assert(aligned[i]!=p);free(p);}
#define ctools_aligned_free release_aligned
#define free release_heap
#include "ctools_types.c"
#undef free
#undef ctools_aligned_free
static void ownership(void) {
 for(int stage=0;stage<5;stage++){
  stata_data d={.nobs=2,.nvars=2};
  if(stage>=1)d.sort_order=alloc_aligned(2*sizeof(perm_idx_t));
  if(stage>=2){d.vars=alloc_aligned(2*sizeof(stata_variable));memset(d.vars,0,2*sizeof(stata_variable));}
  if(stage>=3){stata_variable *v=d.vars;v->nobs=2;v->data.dbl=alloc_aligned(2*sizeof(double));}
  if(stage>=4){stata_variable *v=d.vars+1;v->type=STATA_TYPE_STRING;v->nobs=2;v->buffer_storage=CTOOLS_BUFFER_HEAP;v->string_storage=CTOOLS_STRINGS_BORROWED;v->data.str=calloc(2,sizeof(char*));v->data.str[0]="borrowed";v->data.str[1]="literal";}
  stata_data_free(&d);stata_data_free(&d);assert(live==0&&!d.vars&&!d.sort_order);
 }
 char *owner=strdup("survives"); char *view[]={owner};
 stata_variable v={.type=STATA_TYPE_STRING,.nobs=1,.data.str=view,.buffer_storage=CTOOLS_BUFFER_BORROWED,.string_storage=CTOOLS_STRINGS_BORROWED};
 stata_variable_free(&v);assert(!strcmp(owner,"survives"));free(owner);
}
int main(void) {
 lifecycle();parsing();ownership();
 assert(ctools_stata_rc(STATA_ERR_MEMORY)==920);
 assert(ctools_stata_rc(STATA_ERR_INVALID_INPUT)==198);
 assert(ctools_stata_rc(STATA_ERR_STATA_READ)==459);
 assert(ctools_stata_rc(STATA_ERR_STATA_WRITE)==459);
 assert(ctools_stata_rc(STATA_ERR_CANCELLED)==1);
 assert(ctools_stata_rc(STATA_ERR_UNSUPPORTED_TYPE)==109);
 puts("PASS parser, lifecycle, partial ownership, allocator pairing, and status contracts");
}
'''


def main():
    with tempfile.TemporaryDirectory(prefix='ctools-contracts-') as tmp:
        c=Path(tmp)/'test.c';exe=Path(tmp)/'test';c.write_text(SOURCE)
        flags=['-std=c11','-O1','-g','-fsanitize='+os.environ.get('CTOOLS_SANITIZERS','undefined'),'-I'+str(ROOT/'src'),'-pthread','-ffunction-sections','-fdata-sections']
        flags += ['-DSYSTEM=APPLEMAC','-D_DARWIN_C_SOURCE'] if platform.system()=='Darwin' else ['-DSYSTEM=STUNIX','-D_GNU_SOURCE']
        flags += ['-Wl,-dead_strip'] if platform.system()=='Darwin' else ['-Wl,--gc-sections']
        subprocess.run([os.environ.get('CC','clang'),*flags,str(c),*[str(ROOT/'src'/p) for p in ['ctools_parse.c','ctools_runtime.c','ctools_threads.c','ctools_arena.c']],'-lm','-o',str(exe)],check=True)
        subprocess.run([str(exe)],check=True)


if __name__=='__main__':main()

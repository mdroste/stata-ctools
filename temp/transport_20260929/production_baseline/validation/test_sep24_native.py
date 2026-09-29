"""Native regressions for SEP24: allocation failure, sorting, pool shutdown, file publication.

Optional CTOOLS_TEST_OPENMP_PREFIX builds the same sort tests with OpenMP and
constrained teams. All binaries run under UBSan with bounded timeouts; set CTOOLS_SANITIZERS=address,undefined to also enable ASan.
"""
import os
from pathlib import Path
import platform
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
ALLOC=r'''
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
static int calls, fail_at, live_count;
static void *live[10000];
static void record(void *p) { if(p) { for(int i=0;i<10000;i++) if(!live[i]) {live[i]=p;live_count++;return;} abort(); } }
static void *a_malloc(size_t n) { if(++calls==fail_at)return NULL;void *p=malloc(n);if(p)memset(p,0xa5,n);record(p);return p; }
static void *a_calloc(size_t n,size_t s) {if(++calls==fail_at)return NULL;void *p=calloc(n,s);record(p);return p;}
static int a_align(void **p,size_t a,size_t n) {if(++calls==fail_at)return ENOMEM;int rc=posix_memalign(p,a,n);if(!rc)record(*p);return rc;}
static void a_free(void *p) {if(!p)return;for(int i=0;i<10000;i++)if(live[i]==p){live[i]=NULL;live_count--;free(p);return;}assert(!"invalid free");}
#define malloc a_malloc
#define calloc a_calloc
#define posix_memalign a_align
#define free a_free
'''

VCE=ALLOC+r'''
#include "ctools_matrix.c"
int main(void) {
 double X[]={1,2,3,4},D[]={1.0/30},e[]={1,-1,2,1},w[]={1,2,1,2},V[1]; int ids[]={0,0,1,1};
 ctools_vce_data d={.X_eff=X,.D=D,.resid=e,.weights=w,.N=4,.K=1,.weight_type=1,.normalize_weights=1};
 for(int cluster=0;cluster<2;cluster++) {
  calls=0;fail_at=0;V[0]=123;
  assert((cluster?ctools_vce_cluster(&d,ids,2,1,V):ctools_vce_robust(&d,1,V))==0);
  int allocations=calls;assert(V[0]>0 && V[0]!=123 && live_count==0);
  for(int f=1;f<=allocations;f++) {calls=0;fail_at=f;V[0]=123;
   assert((cluster?ctools_vce_cluster(&d,ids,2,1,V):ctools_vce_robust(&d,1,V))==920);
   assert(V[0]==123 && live_count==0);
  }
 }
 fail_at=0; memset(e,0,sizeof(e));assert(ctools_vce_robust(&d,1,V)==0 && V[0]==0);
 V[0]=123;assert(ctools_vce_cluster(&d,ids,1,1,V)==498 && V[0]==123);
 e[0]=1e308;V[0]=123;assert(ctools_vce_robust(&d,1,V)==498 && V[0]==123);
 assert(live_count==0);return 0;
}
'''

ARENA=ALLOC+r'''
#include "ctools_arena.c"
int main(void) {
 ctools_arena a;ctools_arena_init(&a,32);
 size_t capacity=0;
 for(int run=0;run<50;run++) {
  assert(ctools_arena_alloc(&a,24));assert(ctools_arena_alloc(&a,64));
  void *aligned=ctools_arena_alloc_aligned(&a,53,64);assert(aligned && (uintptr_t)aligned%64==0);
  assert(!ctools_arena_alloc(&a,SIZE_MAX));assert(!ctools_arena_alloc_aligned(&a,SIZE_MAX,64));
  if(run==0)capacity=live_count;else assert(capacity==(size_t)live_count);
  ctools_arena_reset(&a);
 }
 ctools_arena_free(&a);assert(live_count==0);return 0;
}
'''

POOL=r'''
#include <assert.h>
#include <errno.h>
#include <pthread.h>
#include "ctools_threads.h"
static int calls, fail_at;
static int create(pthread_t *t,const pthread_attr_t *a,void *(*f)(void*),void *p) {
 if(++calls==fail_at)return EAGAIN;return pthread_create(t,a,f,p);
}
#define pthread_create create
#include "ctools_threads.c"
int main(void) {
 for(int run=0;run<100;run++)for(int f=1;f<=4;f++) {
  ctools_persistent_pool pool;calls=0;fail_at=f;
  assert(ctools_persistent_pool_init(&pool,4)!=0);
  assert(!pool.initialized && !pool.workers && pool.num_workers==0);
  fail_at=0;assert(ctools_persistent_pool_init(&pool,2)==0);
  ctools_persistent_pool_destroy(&pool);
 }
 return 0;
}
'''

SORT_FAULT=ALLOC+r'''
#include "ctools_types.h"
int ctools_get_max_threads(void) {return 4;}
#include "ctools_sort_ALGORITHM.c"
int main(void) {
 perm_idx_t order[128];uint64_t keys[128];char *strings[128];
 for(int i=0;i<128;i++){keys[i]=128-i;strings[i]=i%2?"az":"ba";}
 for(int string=0;string<2;string++) {
  calls=0;fail_at=0;for(int i=0;i<128;i++)order[i]=i;
  assert((string?STRING_CALL:NUMERIC_CALL)==0);int allocations=calls;assert(live_count==0);
  for(int f=1;f<=allocations;f++) {
   calls=0;fail_at=f;for(int i=0;i<128;i++)order[i]=i;
   assert((string?STRING_CALL:NUMERIC_CALL)!=0);assert(live_count==0);
  }
 }
 return 0;
}
'''
SORT=r'''
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
#include "ctools_types.h"
static ST_plugin mock;ST_plugin *_stata_=&mock;
static ST_boolean missing(ST_double x){return x>=8.98846567431158e307;}
int ctools_get_max_threads(void) {return 4;}
#include "ctools_sort_ALGORITHM.c"
int main(void) {
 mock.ismissing=missing;mock.missval=8.98846567431158e307;
 const size_t N=400001;int key=1;
 stata_variable var={.type=STATA_TYPE_STRING,.nobs=N};
 char **s=malloc(N*sizeof(char*));double *x=malloc(N*sizeof(double));
 perm_idx_t *order=malloc(N*sizeof(perm_idx_t));unsigned char *seen=malloc(N);
 stata_data d={.nobs=N,.nvars=1,.vars=&var,.sort_order=order};
 for(int string=0;string<2;string++) {
  for(size_t i=0;i<N;i++){s[i]=i%7==0?"":(i%3==0?"a":(i%2?"az":"ba"));x[i]=(double)((N-i)%2345);order[i]=i;}
  var.type=string?STATA_TYPE_STRING:STATA_TYPE_DOUBLE;if(string)var.data.str=s;else var.data.dbl=x;
  assert(ctools_sort_ALGORITHM_order_only(&d,&key,1)==0);memset(seen,0,N);
  for(size_t i=0;i<N;i++){assert(order[i]<N && !seen[order[i]]);seen[order[i]]=1;if(i){if(string)assert(strcmp(s[order[i-1]],s[order[i]])<=0);else assert(x[order[i-1]]<=x[order[i]]);}}
 }
 free(s);free(x);free(order);free(seen);return 0;
}
'''

ATOMIC=r'''
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
static ST_plugin mock;ST_plugin *_stata_=&mock;
#include "cexport/cexport_parse.c"
static void put(const char *p,const char *s){FILE *f=fopen(p,"w");assert(f);assert(fputs(s,f)>=0);assert(!fclose(f));}
static void check(const char *p,const char *s){char buf[30]={0};FILE *f=fopen(p,"r");assert(f);fread(buf,1,29,f);fclose(f);assert(!strcmp(buf,s));}
int main(int argc,char **argv){
 assert(argc==2);cexport_output out;const char *p=argv[1];put(p,"original");
 assert(!cexport_output_prepare(&out,p,1));put(out.temporary,"partial");cexport_output_cleanup(&out);check(p,"original");
 assert(!cexport_output_prepare(&out,p,0));put(out.temporary,"new");assert(cexport_output_commit(&out)==602);cexport_output_cleanup(&out);check(p,"original");
 assert(!cexport_output_prepare(&out,p,1));put(out.temporary,"complete");assert(!cexport_output_commit(&out));cexport_output_cleanup(&out);check(p,"complete");
 remove(p);return 0;
}
'''

KEY_GUARD=ALLOC+r'''
#include "ctools_types.h"
#include "cmerge/cmerge_keys.c"
#include "cmerge/cmerge_group_search.c"
#include "cmerge/cmerge_join.c"
int main(void) {
 double number[]={1};char *text[]={"1"};
 stata_variable a[2]={{.type=STATA_TYPE_DOUBLE,.data.dbl=number},{.type=STATA_TYPE_STRING,.data.str=text}};
 stata_variable b[2]={{.type=STATA_TYPE_DOUBLE,.data.dbl=number},{.type=STATA_TYPE_DOUBLE,.data.dbl=number}};
 stata_data da={.nobs=1,.nvars=2,.vars=a},db={.nobs=1,.nvars=2,.vars=b};
 cmerge_output_spec_t *output=(void*)1;
 assert(cmerge_sorted_join(&da,&db,2,MERGE_1_1,&output)==-4 && !output);
 assert(cmerge_sorted_join(&db,&da,2,MERGE_M_M,&output)==-4 && !output);
 assert(live_count==0);return 0;
}
'''

RANK_SOLVE=ALLOC+r'''
#include "cbinscatter/cbinscatter_resid.c"
int main(void) {
 double gram[9]={0},rhs[3]={0},beta[3];
 for(int i=0;i<5;i++){double x[]={1,1e-9*i,2e-9*i},y=3+2*i;
  for(int j=0;j<3;j++){rhs[j]+=x[j]*y;for(int k=0;k<3;k++)gram[j*3+k]+=x[j]*x[k];}}
 assert(cbinscatter_solve_with_collinearity(gram,rhs,3,beta)==0);
 int allocations=calls;assert(live_count==0);
 for(int i=0;i<5;i++)assert(fabs(beta[0]+1e-9*i*beta[1]+2e-9*i*beta[2]-(3+2*i))<1e-10);
 for(int f=1;f<=allocations;f++) {calls=0;fail_at=f;beta[0]=beta[1]=beta[2]=123;
  assert(cbinscatter_solve_with_collinearity(gram,rhs,3,beta)==920);
  assert(beta[0]==123 && beta[1]==123 && beta[2]==123 && live_count==0);}
 fail_at=0;gram[0]=-1;
 assert(cbinscatter_solve_with_collinearity(gram,rhs,3,beta)==199);
 assert(live_count==0);return 0;
}
'''

def compile_run(tmp,name,source,flags=(),env=None,args=()):
    cfile=tmp/(name+'.c');binary=tmp/name;cfile.write_text(source)
    base=['-std=c11','-O1','-g','-fno-fast-math','-fno-omit-frame-pointer','-fsanitize='+os.environ.get('CTOOLS_SANITIZERS','undefined'),'-I',str(ROOT/'src')]
    if platform.system()=='Darwin':base+=['-DSYSTEM=APPLEMAC','-D_DARWIN_C_SOURCE','-Wl,-dead_strip','-Wl,-undefined,dynamic_lookup']
    else:base+=['-DSYSTEM=STUNIX','-D_POSIX_C_SOURCE=200809L','-ffunction-sections','-Wl,--gc-sections']
    subprocess.run([os.environ.get('CC','/usr/bin/clang' if platform.system()=='Darwin' else 'clang'),*base,str(cfile),*flags,'-pthread','-lm','-o',str(binary)],check=True,capture_output=True,text=True)
    run_env=dict(os.environ,ASAN_OPTIONS='detect_leaks=0:symbolize=0:abort_on_error=1',UBSAN_OPTIONS='halt_on_error=1:symbolize=0',**(env or {}))
    subprocess.run([str(binary),*map(str,args)],check=True,env=run_env,timeout=45,capture_output=True,text=True)
    print('PASS',name,flush=True)

def main():
    with tempfile.TemporaryDirectory(prefix='ctools-sep24-native-') as directory:
        tmp=Path(directory)
        for name,source in [('key_type_guard',KEY_GUARD),('rank_solve',RANK_SOLVE),('vce_faults',VCE),('arena_reset',ARENA),('pool_create_failure',POOL),('atomic_output',ATOMIC)]:
            compile_run(tmp,name,source,args=[tmp/'literal threads(2) path'] if name=='atomic_output' else [])
        for alg,num,string in [('sample','sample_sort_numeric_impl(order,keys,128,4)','sample_sort_string_impl(order,strings,128,4,2)'),('merge','parallel_merge_sort_numeric(order,keys,128,4)','parallel_merge_sort_string(order,strings,128,4,2)')]:
            compile_run(tmp,alg+'_faults',SORT_FAULT.replace('ALGORITHM',alg).replace('NUMERIC_CALL',num).replace('STRING_CALL',string))
        for alg in ['sample','merge','radix_lsd','counting','ips4o']:
            # Public LSD name omits "radix_".
            source=SORT.replace('ALGORITHM',alg)
            if alg == 'counting': source=source.replace('string<2','string<1')
            extra=[str(ROOT/'src/ctools_arena.c')] if alg in ('counting','ips4o') else []
            compile_run(tmp,alg+'_noomp',source,extra)
            prefix=os.environ.get('CTOOLS_TEST_OPENMP_PREFIX')
            if prefix:
                flags=['-Xpreprocessor','-fopenmp','-I',str(Path(prefix)/'include'),str(Path(prefix)/'lib/libomp.a')]
                for limit in ['1','2']:
                    compile_run(tmp,alg+'_team'+limit,source,flags+extra,{'OMP_THREAD_LIMIT':limit,'OMP_DYNAMIC':'TRUE'})

if __name__=='__main__':
    try:main()
    except subprocess.CalledProcessError as e:
        print(e.stdout or '',e.stderr or '',flush=True);raise

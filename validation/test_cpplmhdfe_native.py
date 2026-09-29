"""ASan/UBSan tests for native ReLU separation and heterogeneous projections."""
import os
from pathlib import Path
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = r'''
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
#include "ctools_hdfe_utils.c"
#include "ctools_matrix.c"
#include "ctools_ols.c"
#include "creghdfe/creghdfe_solver.c"
#include "cpplmhdfe/cpplmhdfe_irls.c"
static int message(char *s) { (void)s; return 0; }
static void init(HDFE_State *s, int n, int groups) {
 memset(s,0,sizeof(*s));s->N=n;s->G=groups;s->num_threads=1;
 s->tolerance=1e-12;s->maxiter=10000;s->has_weights=1;
 s->weights=calloc(n,sizeof(double));s->factors=calloc(groups,sizeof(FE_Factor));
 for(int g=0;g<groups;g++) {
  FE_Factor *f=&s->factors[g];f->num_levels=2;f->has_intercept=1;
  f->levels=calloc(n,sizeof(int));f->counts=calloc(2,sizeof(double));
  f->weighted_counts=calloc(2,sizeof(double));
  for(int i=0;i<n;i++){f->levels[i]=1+i%2;f->counts[i%2]++;}
 }
 assert(!ctools_hdfe_alloc_buffers(s,0,4));
}
static void destroy(HDFE_State *s) {free(s->weights);ctools_hdfe_state_cleanup(s);}
int main(void) {
 mock.spouterr=message;mock.spoutsml=message;
 for(int run=0;run<30;run++) {
  HDFE_State s;init(&s,60,1);
  double y[60],x[120];int sep[60];
  for(int i=0;i<60;i++){y[i]=i<15?0:1;x[i]=sin(i);x[60+i]=x[i]+(i<15);}
  assert(ppml_relu(&s,y,x,60,2,sep,1)==15);
  for(int i=0;i<60;i++)assert(sep[i]==(i<15));
  /* A tiny positive outcome must impose the equality constraint too. */
  y[0]=1e-200;assert(ppml_relu(&s,y,x,60,2,sep,1)==0);
  /* Zeros distributed over a full-rank interior are not separation. */
  for(int i=0;i<60;i++){y[i]=(i%5)?1:0;x[60+i]=cos(i);}
  assert(ppml_relu(&s,y,x,60,2,sep,1)==0);
  destroy(&s);
 }
 HDFE_State s;init(&s,60,2);
 double y[60],x[60];int sep[60];
 for(int i=0;i<60;i++) {
  s.factors[0].levels[i]=1+(i>=20);
  s.factors[1].levels[i]=1+(i>=40);
  y[i]=(i>=20 && i<40)?0:1;x[i]=sin(i);
 }
 assert(ppml_relu(&s,y,x,60,1,sep,1)==20);
 for(int i=0;i<60;i++)assert(sep[i]==(i>=20&&i<40));
 destroy(&s);
 init(&s,60,2);
 s.factors[1].slope=calloc(60,sizeof(double));
 s.factors[1].slope_center=calloc(2,sizeof(double));
 for(int i=0;i<60;i++) {
  s.factors[1].slope[i]=i/60.0;
  s.weights[i]=1+(i%3);
  y[i]=(i%2?3:2)+(i%2?-2:4)*s.factors[1].slope[i];
 }
 ppml_update_fe_weights(&s,s.weights,60);
 HDFE_SolveResult r=partial_out_columns(&s,y,60,1,1);
 assert(!r.status);
 for(int i=0;i<60;i++)assert(fabs(y[i])<1e-9);
 destroy(&s);
 {
  double a[]={1,0, 0,1, 1,1};int mask[4]={0};
  assert(ppml_simplex_tableau(a,3,2,100,mask)==3);
  double b[]={1,0, -1,0, 0,1, 0,-1};
  assert(ppml_simplex_tableau(b,4,2,100,mask)==0);
  double c[]={1,1, 1,-1, 0,1, 0,-1};
  assert(ppml_simplex_tableau(c,4,2,100,mask)==2);
  assert(mask[0] && mask[1] && !mask[2] && !mask[3]);
  double duplicate[]={1,1, 2,2, -1,-1};
  assert(ppml_simplex_tableau(duplicate,3,2,100,mask)==0);
 }
 init(&s,60,2);
 for(int i=0;i<60;i++) {
  s.factors[0].levels[i]=1+(i%2);
  s.factors[1].levels[i]=1+(i>=30);
  s.weights[i]=1+(i%3);y[i]=3*(i%2)+7*(i>=30)+sin(i);
 }
 ppml_update_fe_weights(&s,s.weights,60);
 double baseline[60];memcpy(baseline,y,sizeof(y));
 assert(!partial_out_columns(&s,baseline,60,1,1).status);
 ppml_options.btol=1e-12;ppml_options.conlim=1e8;
 assert(!ppml_lsmr(&s,y,60,1).status);
 for(int i=0;i<60;i++)assert(fabs(y[i]-baseline[i])<1e-9);
 for(int method=0;method<5;method++)for(int transform=0;transform<3;transform++) {
  if(transform==1 && (method==0 || method==4))continue;
  for(int i=0;i<60;i++)y[i]=3*(i%2)+7*(i>=30)+sin(i);
  ppml_options.acceleration=method;ppml_options.transform=transform;
  assert(!ppml_map_options(&s,y,60,1).status);
  for(int i=0;i<60;i++)assert(fabs(y[i]-baseline[i])<1e-8);
 }
 destroy(&s);
 return 0;
}
'''

def main():
    with tempfile.TemporaryDirectory(prefix='ctools-ppml-native-') as tmp:
        src=Path(tmp)/'test.c'; binary=Path(tmp)/'test'
        src.write_text(SOURCE)
        flags=['-std=c11','-O1','-g','-fno-fast-math','-fno-omit-frame-pointer',
               '-fsanitize=address,undefined','-I',str(ROOT/'src')]
        if platform.system()=='Darwin':
            flags += ['-DSYSTEM=APPLEMAC','-Wl,-dead_strip']
        else:
            flags += ['-DSYSTEM=STUNIX','-ffunction-sections','-fdata-sections','-Wl,--gc-sections']
        subprocess.run([os.environ.get('CC','/usr/bin/clang' if platform.system()=='Darwin' else 'clang'),*flags,str(src),'-lm','-o',str(binary)],check=True)
        subprocess.run([str(binary)],check=True,timeout=60)
        print('PASS: native regressor/FE separation, overlap, tiny positive outcomes, heterogeneous slopes, simplex, LSMR (ASan/UBSan)')

if __name__=='__main__': main()

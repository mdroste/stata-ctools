"""Exercise stable numeric radix and tiled merges against an independent qsort.

Uses the real OpenMP runtime (including teams smaller than the requested count).
Set CTOOLS_LIBOMP_PREFIX on macOS. --source-root can test a frozen baseline.
"""
from pathlib import Path
import argparse
import os
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = r'''
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ctools_order.h"
#include "ctools_config.h"
static int threads = 1, keys[] = {0,1};
static stata_data data;
static size_t *rank_in;
int ctools_get_max_threads(void) { return threads; }
static int reference(const void *pa, const void *pb) {
    size_t a=*(const perm_idx_t*)pa, b=*(const perm_idx_t*)pb;
    for(int k=0;k<2;k++) {
        int c;
        if(data.vars[k].type==STATA_TYPE_STRING) c=strcmp(data.vars[k].data.str[a],data.vars[k].data.str[b]);
        else {double x=data.vars[k].data.dbl[a],y=data.vars[k].data.dbl[b];c=(x>y)-(x<y);}
        if(c)return c;
    }
    return (rank_in[a]>rank_in[b])-(rank_in[a]<rank_in[b]);
}
static unsigned long long seed=923842;
static unsigned rnd(void) { seed=seed*6364136223846793005ULL+1;return seed>>32; }
static void check(size_t n,int strings) {
    data.nobs=n;data.nvars=2;data.vars=calloc(2,sizeof(*data.vars));
    data.sort_order=malloc((n+1)*sizeof(perm_idx_t));
    rank_in=malloc((n+1)*sizeof(*rank_in));
    perm_idx_t *expected=malloc((n+1)*sizeof(*expected)), *input=malloc((n+1)*sizeof(*input));
    char (*text)[24]=malloc((n+1)*24);
    for(int k=0;k<2;k++) {data.vars[k].type=STATA_TYPE_DOUBLE;data.vars[k].data.dbl=malloc((n+1)*sizeof(double));}
    if(strings){free(data.vars[0].data.dbl);data.vars[0].type=STATA_TYPE_STRING;data.vars[0].data.str=malloc((n+1)*sizeof(char*));}
    for(size_t i=0;i<n;i++) {
        double value=(int)(rnd()%1024)-512;
        if(i%17==0)value=0;
        if(i%19==0)value=-0.0;
        if(i%23==0) {value=8.98846567431158e307;for(unsigned j=0;j<i%27;j++)value=nextafter(value,INFINITY);}
        if(strings){snprintf(text[i],24,"group%u",rnd()%101);data.vars[0].data.str[i]=text[i];}
        else data.vars[0].data.dbl[i]=value;
        data.vars[1].data.dbl[i]=(int)(rnd()%13)-6;
        input[i]=(perm_idx_t)i;
    }
    for(size_t i=n;i>1;i--){size_t j=rnd()%i;perm_idx_t t=input[i-1];input[i-1]=input[j];input[j]=t;}
    for(size_t i=0;i<n;i++)rank_in[input[i]]=i;
    memcpy(expected,input,n*sizeof(*expected));qsort(expected,n,sizeof(*expected),reference);
    int counts[]={1,2,12};
    for(int t=0;t<3;t++) {
        threads=counts[t];memcpy(data.sort_order,input,n*sizeof(*input));
        assert(ctools_order_stable(&data,keys,2)==STATA_OK);
        assert(!memcmp(data.sort_order,expected,n*sizeof(*expected)));
        assert(ctools_order_stable(&data,keys,2)==STATA_OK);
        assert(!memcmp(data.sort_order,expected,n*sizeof(*expected)));
    }
    free(text);free(expected);free(input);free(rank_in);free(data.sort_order);
    free(data.vars[0].data.dbl);free(data.vars[1].data.dbl);free(data.vars);
}
int main(void) {
    size_t sizes[]={0,1,257,32767,32768,32769,200003};
    for(size_t s=0;s<sizeof(sizes)/sizeof(*sizes);s++)for(int strings=0;strings<2;strings++)check(sizes[s],strings);
    puts("PASS stable parallel order: numeric/string keys, signed zeros, extended missings, partial tiles, incoming tie order, threads 1/2/12");
}
'''

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source-root',type=Path,default=ROOT)
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='ctools-order-') as tmp:
        cfile=Path(tmp)/'test.c'; binary=Path(tmp)/'test'
        cfile.write_text(SOURCE)
        flags=['-std=c11','-O2','-g','-fsanitize=undefined','-fno-fast-math','-ffp-contract=off',
               '-Wall','-Wextra','-Werror','-I',str(args.source_root/'src')]
        if platform.system()=='Darwin':
            prefix=Path(os.environ.get('CTOOLS_LIBOMP_PREFIX','/private/tmp/ctools-p2-fixes/libomp-arm21'))
            flags+=['-DSYSTEM=APPLEMAC','-D_DARWIN_C_SOURCE','-Xpreprocessor','-fopenmp','-I',str(prefix/'include')]
            libs=[str(prefix/'lib/libomp.a')]
        else:
            flags+=['-DSYSTEM=STUNIX','-D_GNU_SOURCE','-fopenmp'];libs=[]
        subprocess.run([os.environ.get('CC','clang'),*flags,str(cfile),str(args.source_root/'src/ctools_order.c'),*libs,'-lm','-o',str(binary)],check=True)
        subprocess.run([str(binary)],check=True,timeout=120)
        subprocess.run([str(binary)],env={**os.environ,'OMP_THREAD_LIMIT':'2'},check=True,timeout=120)

if __name__=='__main__':
    main()

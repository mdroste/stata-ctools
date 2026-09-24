"""ASan/UBSan checks of sample selection, cluster identity, and merge contracts."""
import os
from pathlib import Path
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]

SOURCE = r'''
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static int calls, fail_at, live;
static ST_boolean ismissing(ST_double x) { return x >= 8.98846567431158e307; }
static void *alloc(size_t n) {
    if (++calls == fail_at) return NULL;
    void *p = malloc(n); if (p) live++; return p;
}
static void *zero(size_t n, size_t k) {
    if (++calls == fail_at) return NULL;
    void *p = calloc(n,k); if (p) live++; return p;
}
static void release(void *p) { if (p) live--; free(p); }
#define malloc alloc
#define calloc zero
#define free release
#include "ctools_hdfe_utils.c"
#include "cpplmhdfe/cpplmhdfe_separation.c"
int main(void) {
    mock.missval = 8.98846567431158e307;
    mock.ismissing = ismissing;
    double codes[] = {-1, -.1, .1, .2, 1e12, 1e12+1, -.1, 8.98846567431158e307};
    int ids[8], groups;
    assert(ctools_numeric_to_cluster_ids(codes,8,ids,&groups) == 0);
    assert(groups == 6 && ids[7] == -1 && ids[1] == ids[6]);
    for (int i=0; i<6; i++) for (int j=0; j<i; j++) assert(ids[i] != ids[j]);
    assert(live == 0);
    calls=0; fail_at=1;
    assert(ctools_numeric_to_cluster_ids(codes,8,ids,&groups) < 0 && live == 0);
    fail_at=0;
    /* Removed blocks at start/middle/end, with an initial excluded row and
       a singleton. Allocations are exact-sized, so ASan sees boundary errors. */
    for (int block=0; block<3; block++) {
        int *level = malloc(14*sizeof(int));
        int *mask = malloc(14*sizeof(int));
        double *y = malloc(14*sizeof(double));
        for (int i=0; i<14; i++) { level[i]=i/4+1; mask[i]=1; y[i]=1; }
        mask[13]=0;
        for (int i=block*4; i<block*4+4; i++) y[i]=0;
        int *levels[]={level}, nlevels[]={4}, singletons, separated;
        int base=live;
        calls=0;
        assert(ppml_select_fe_sample(y,levels,nlevels,1,14,mask,&singletons,&separated)==0);
        int allocations=calls;
        assert(singletons==1 && separated==4 && live==base);
        for (int i=0; i<14; i++) assert(mask[i] == (i<12 && i/4!=block));
        for (int failure=1; failure<=allocations; failure++) {
            for(int i=0;i<14;i++) mask[i]=i!=13;
            calls=0; fail_at=failure;
            assert(ppml_select_fe_sample(y,levels,nlevels,1,14,mask,&singletons,&separated)<0);
            assert(live==base);
        }
        fail_at=0;
        free(level); free(mask); free(y);
        assert(live==0);
    }
    /* Removing a zero group creates a singleton in a crossing FE. */
    int a[]={1,1,2,2,2}, b[]={1,1,1,2,2}, mask[]={1,1,1,1,1};
    double y[]={0,0,1,1,1};
    int *levels[]={a,b}, counts[]={2,2}, ns, nz;
    assert(ppml_select_fe_sample(y,levels,counts,2,5,mask,&ns,&nz)==0);
    assert(nz==2 && ns==1 && !mask[2] && mask[3] && mask[4]);
    assert(live==0);
    return 0;
}
'''

def main():
    with tempfile.TemporaryDirectory(prefix='ctools-sep22-native-') as tmp:
        cfile=Path(tmp)/'test.c'; binary=Path(tmp)/'test'
        cfile.write_text(SOURCE)
        flags=['-std=c11','-O1','-g','-fno-fast-math','-fno-omit-frame-pointer',
               '-fsanitize=address,undefined','-I',str(ROOT/'src')]
        if platform.system() == 'Darwin':
            flags += ['-DSYSTEM=APPLEMAC','-Wl,-dead_strip','-Wl,-undefined,dynamic_lookup']
        else:
            flags += ['-DSYSTEM=STUNIX','-ffunction-sections','-fdata-sections','-Wl,--gc-sections']
        subprocess.run([os.environ.get('CC', '/usr/bin/clang' if platform.system() == 'Darwin' else 'clang'),*flags,str(cfile),'-lm','-o',str(binary)],check=True)
        subprocess.run([str(binary)],check=True,timeout=60)
        print('PASS: full-double clusters; original-row selection; singleton chains; every selection allocation failure (ASan/UBSan)')

if __name__ == '__main__': main()

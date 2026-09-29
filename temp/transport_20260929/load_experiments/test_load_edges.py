from pathlib import Path
import os,subprocess,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'validation'))
from test_p2_native import SPI
SOURCE=SPI+r'''
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
#include "ctools_data_io.c"
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
    return 0;
}
'''
source_dir=Path(sys.argv[1]).resolve();out=Path(__file__).resolve().parent/source_dir.parent.name
cfile=out/'edges.c';cfile.write_text(SOURCE)
for omp in (False,True):
    binary=out/('edges_omp' if omp else 'edges_serial')
    flags=['-std=c11','-O2','-DSD_FASTMODE','-fno-fast-math','-fsanitize=undefined','-fno-omit-frame-pointer','-pthread','-I'+str(source_dir),'-DSYSTEM=APPLEMAC','-D_DARWIN_C_SOURCE','-Wl,-dead_strip','-Wl,-undefined,dynamic_lookup']
    if omp:flags+=['-Xpreprocessor','-fopenmp','-I/opt/homebrew/opt/libomp/include','/opt/homebrew/opt/libomp/lib/libomp.a']
    subprocess.run(['clang',*flags,str(cfile),*[str(source_dir/f) for f in ['ctools_types.c','ctools_arena.c','ctools_threads.c']],'-o',str(binary)],check=True)
    subprocess.run([str(binary)],check=True,timeout=120,env=dict(os.environ,UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1',KMP_BLOCKTIME='0',OMP_WAIT_POLICY='PASSIVE'))
    print(f'PASS load edges omp={omp} source={source_dir}',flush=True)

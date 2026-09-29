from pathlib import Path
import os,subprocess,sys
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'validation'))
from test_p2_native import SPI
SOURCE=SPI+r'''
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
'''
source_dir=Path(sys.argv[1]).resolve();out=Path(__file__).resolve().parent/source_dir.parent.name;cfile=out/'no_order_test.c';cfile.write_text(SOURCE)
for omp in (False,True):
 binary=out/('no_order_omp' if omp else 'no_order_serial')
 flags=['-std=c11','-O2','-DSD_FASTMODE','-fno-fast-math','-fsanitize=undefined','-fno-omit-frame-pointer','-pthread','-I'+str(source_dir),'-DSYSTEM=APPLEMAC','-D_DARWIN_C_SOURCE','-Wl,-dead_strip','-Wl,-undefined,dynamic_lookup']
 if omp:flags+=['-Xpreprocessor','-fopenmp','-I/opt/homebrew/opt/libomp/include','/opt/homebrew/opt/libomp/lib/libomp.a']
 subprocess.run(['clang',*flags,str(cfile),*[str(source_dir/f) for f in ['ctools_types.c','ctools_arena.c','ctools_threads.c']],'-o',str(binary)],check=True)
 subprocess.run([str(binary)],check=True,timeout=120,env=dict(os.environ,UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1',KMP_BLOCKTIME='0',OMP_WAIT_POLICY='PASSIVE'))
 print(f'PASS no sort order omp={omp} source={source_dir}',flush=True)

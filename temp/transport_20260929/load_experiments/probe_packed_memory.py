from pathlib import Path
import os,subprocess,sys,json
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'validation'))
from test_p2_native import SPI
SOURCE=SPI+r'''
#include <stdio.h>
#include <sys/resource.h>
#include <mach/mach.h>
#include <time.h>
#include "ctools_data_io.c"
static int ncols=4,content_len=8;
static ST_int wide_vars(void){return ncols;}
static ST_boolean text_col(ST_int v){(void)v;return 1;}
static ST_boolean no_strl(ST_int v){(void)v;return 0;}
static ST_retcode text(ST_int v,ST_int row,char *out){(void)v;memset(out,'a'+row%26,content_len);out[content_len]=0;return 0;}
static ST_retcode width_macro(char *name,char *out,ST_int size){(void)name;(void)size;strcpy(out,"2045,2045,2045,2045");return 0;}
static double now(void){struct timespec ts;clock_gettime(CLOCK_MONOTONIC,&ts);return ts.tv_sec+1e-9*ts.tv_nsec;}
int main(int argc,char **argv){
 setup();nrows=100000;content_len=argc>1?atoi(argv[1]):8;
 mock.nvars=mock.nvar=wide_vars;mock.isstr=text_col;mock.isstrl=no_strl;mock.sdata=text;mock.macuse=width_macro;
 ctools_set_max_threads(4);ctools_filtered_data fd;double start=now();
 assert(ctools_data_load(&fd,NULL,0,1,nrows,CTOOLS_LOAD_SKIP_IF)==0);
 double elapsed=now()-start;struct rusage usage;getrusage(RUSAGE_SELF,&usage);
 mach_task_basic_info_data_t mem;mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;
 assert(task_info(mach_task_self(),MACH_TASK_BASIC_INFO,(task_info_t)&mem,&count)==KERN_SUCCESS);
 for(size_t j=0;j<fd.data.nvars;j++)for(size_t i=0;i<fd.data.nobs;i++) {
    assert(strlen(fd.data.vars[j].data.str[i])==(size_t)content_len);
    if(content_len)assert(fd.data.vars[j].data.str[i][0]=='a'+(i+1)%26);
 }
 start=now();ctools_filtered_data_free(&fd);double cleanup=now()-start;
 printf("{\"length\":%d,\"load_s\":%.6f,\"cleanup_s\":%.6f,\"maxrss_bytes\":%ld,\"resident_bytes\":%llu,\"virtual_bytes\":%llu}\n",content_len,elapsed,cleanup,usage.ru_maxrss,(unsigned long long)mem.resident_size,(unsigned long long)mem.virtual_size);
 ctools_destroy_global_pool();return 0;
}
'''
out=Path(__file__).resolve().parent
results=[]
for name in ('baseline','str2045_packed','str2045'):
 source_dir=out.parent/'baseline/src' if name=='baseline' else out/name/'src'
 binary=out/('memory_'+name);cfile=out/'memory_probe.c';cfile.write_text(SOURCE)
 flags=['-std=c11','-O3','-DSD_FASTMODE','-fno-fast-math','-pthread','-I'+str(source_dir),'-DSYSTEM=APPLEMAC','-D_DARWIN_C_SOURCE','-Wl,-dead_strip','-Wl,-undefined,dynamic_lookup','-Xpreprocessor','-fopenmp','-I/opt/homebrew/opt/libomp/include','/opt/homebrew/opt/libomp/lib/libomp.a']
 subprocess.run(['clang',*flags,str(cfile),*[str(source_dir/f) for f in ['ctools_types.c','ctools_arena.c','ctools_threads.c']],'-o',str(binary)],check=True)
 for length in (0,8,64,2045):
  result=json.loads(subprocess.check_output([str(binary),str(length)],text=True,env=dict(os.environ,KMP_BLOCKTIME='0',OMP_WAIT_POLICY='PASSIVE')))
  result['variant']=name;results.append(result);print(json.dumps(result),flush=True)
(out/'memory_probe.json').write_text(json.dumps(results,indent=2)+'\n')

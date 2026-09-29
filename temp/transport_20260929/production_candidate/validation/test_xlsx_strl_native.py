"""ASan/UBSan checks for XLSX's dedicated length-aware strL loader."""
import os
import platform
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = r'''
#include <assert.h>
#include <string.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
static char payload[140000];
static int length, fail_read, fail_numeric, binary, cancel;
static ST_int first(void){return 1;}
static ST_int last(void){return 3;}
static ST_boolean selected(ST_int obs){return obs!=2;}
static ST_boolean isstring(ST_int var){return var!=2;}
static ST_boolean islong(ST_int var){return var==1;}
static ST_boolean isbinary(ST_int var,ST_int obs){(void)var;(void)obs;return binary;}
static ST_int bytes(ST_int var,ST_int obs){(void)var;return obs==1?length+7:8;}
static ST_retcode readlong(ST_int var,ST_int obs,char *out,ST_int cap){
    (void)var;int n=bytes(var,obs);assert(cap==n+1);
    if(fail_read)return -1;
    memset(out,0,(size_t)n);if(obs==1)memcpy(out,payload,(size_t)length);else memcpy(out,"short",5);
    return n;
}
static ST_retcode readfixed(ST_int var,ST_int obs,char *out){(void)var;(void)obs;strcpy(out,"fixed");return 0;}
static ST_retcode readnumber(ST_int var,ST_int obs,double *out){(void)var;if(fail_numeric)return 459;*out=obs+0.25;return 0;}
static ST_int poll(void){return cancel;}
#include "cexport/cexport_xlsx.c"
static void settext(const char *unit,int count){
    size_t width=strlen(unit);length=(int)width*count;
    for(int i=0;i<count;i++)memcpy(payload+(size_t)i*width,unit,width);
    payload[length]=0;
}
static void run(int expected){
    XLSXExportContext ctx={0};ctx.nvars=3;
    ctools_filtered_data fd;ctools_filtered_data_init(&fd);
    assert(xlsx_load_long_text(&ctx,&fd)==expected);
    if(!expected){
        assert(fd.data.nvars==3 && fd.data.nobs==2 && fd.was_filtered);
        assert(fd.obs_map[0]==1 && fd.obs_map[1]==3);
        assert(strlen(fd.data.vars[0].data.str[0])==(size_t)length);
        assert(!strcmp(fd.data.vars[0].data.str[0],payload));
        assert(!strcmp(fd.data.vars[0].data.str[1],"short"));
        assert(fd.data.vars[1].data.dbl[0]==1.25 && fd.data.vars[1].data.dbl[1]==3.25);
        assert(!strcmp(fd.data.vars[2].data.str[0],"fixed"));
    }
    ctools_filtered_data_free(&fd);ctools_filtered_data_free(&fd);
}
int main(void){
    mock.nobs1=first;mock.nobs2=last;mock.selobs=selected;mock.pollstd=poll;
    mock.isstr=isstring;mock.isstrl=islong;mock.isbinary=isbinary;
    mock.sdatalen=bytes;mock.strldata=readlong;mock.sdata=readfixed;mock.safevdata=readnumber;
    settext("a",5000);run(0);
    fail_read=1;run(609);fail_read=0;
    fail_numeric=1;run(459);fail_numeric=0;
    binary=1;run(109);binary=0;
    cancel=1;run(1);cancel=0;
    settext("\xc3\xa9",32767);run(0);
    settext("\xc3\xa9",32768);run(109);
    settext("\xf0\x9f\x98\x80",16000);run(0);
    settext("\xf0\x9f\x98\x80",17000);run(109);
    settext("a",100000);
    ctools_filtered_data fd;ctools_filtered_data_init(&fd);
    assert(!cexport_load_text_data(&fd,3,false));
    assert(strlen(fd.data.vars[0].data.str[0])==100000);
    ctools_filtered_data_free(&fd);
    return 0;
}
'''


def main():
    with tempfile.TemporaryDirectory(prefix='ctools-xlsx-strl-') as tmp:
        directory = Path(tmp)
        source = directory / 'test.c'
        source.write_text(SOURCE)
        binary = directory / 'test'
        flags = ['-std=c11', '-O1', '-g', '-DSD_FASTMODE', '-fno-fast-math',
                 '-ffp-contract=off', '-fsanitize=address,undefined',
                 '-fno-omit-frame-pointer', '-pthread', '-I', str(ROOT / 'src')]
        if platform.system() == 'Darwin':
            flags += ['-DSYSTEM=APPLEMAC', '-D_DARWIN_C_SOURCE', '-Wl,-dead_strip',
                      '-Wl,-undefined,dynamic_lookup']
        else:
            flags += ['-DSYSTEM=STUNIX', '-D_GNU_SOURCE', '-ffunction-sections',
                      '-fdata-sections', '-Wl,--gc-sections', '-lm']
        modules = ['ctools_data_io.c', 'ctools_types.c', 'ctools_arena.c', 'ctools_threads.c']
        subprocess.run([os.environ.get('CC', 'cc'), *flags, str(source),
                        *[str(ROOT / 'src' / name) for name in modules], '-o', str(binary)], check=True)
        subprocess.run([str(binary)], check=True, timeout=60,
                       env={**os.environ, 'UBSAN_OPTIONS': 'halt_on_error=1'})
        print('PASS XLSX long text, filtering, SPI failures, binary rejection, UTF-16 limit (ASan/UBSan)')


if __name__ == '__main__':
    main()

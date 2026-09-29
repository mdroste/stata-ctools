"""ASan/UBSan regression for BIFF continuation, formulas and OLE stream edits."""
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
ST_plugin *_stata_;
#include "io/cexport_xls.c"
static FILE *book(const char *text,unsigned row,unsigned col,int formula) {
    FILE *f=tmpfile();assert(f);assert(!bof(f,0));
    unsigned char b[32]={0};b[6]=4;b[7]=1;
    for(int i=0;i<4;i++)put16(b+8+2*i,"Data"[i]);
    long bound=ftell(f)+4;assert(!record(f,0x85,b,16));
    memset(b,0,sizeof(b));for(int i=0;i<3;i++)assert(!record(f,0xe0,b,20));
    xls_sst s={0};s.file=f;s.first=1;s.length=8;put32(s.data,1);put32(s.data+4,1);
    assert(!sst_string(&s,text));assert(!sst_flush(&s));assert(!record(f,0x0a,NULL,0));
    long offset=ftell(f);put32(b,(uint32_t)offset);assert(!fseek(f,bound,SEEK_SET));
    assert(fwrite(b,1,4,f)==4);assert(!fseek(f,offset,SEEK_SET));assert(!bof(f,1));
    memset(b,0,sizeof(b));put32(b+4,10);put16(b+10,8);assert(!record(f,0x200,b,14));
    assert(!string_cell(f,row,col,0));
    if(formula) {
        memset(b,0,sizeof(b));put16(b,8);put16(b+2,5);
        double value=99;uint64_t bits;memcpy(&bits,&value,8);put64(b+6,bits);
        put16(b+20,3);b[22]=0x1e;put16(b+23,99);assert(!record(f,6,b,25));
    }
    assert(!record(f,0x0a,NULL,0));return f;
}
int main(int argc,char **argv) {
    assert(argc==2);char path[4096],edited[4096];
    snprintf(path,sizeof(path),"%s/base.xls",argv[1]);snprintf(edited,sizeof(edited),"%s/edit.xls",argv[1]);
    char *a=malloc(24001),*b=malloc(40001);assert(a&&b);
    for(int i=0;i<12000;i++){a[2*i]=(char)0xce;a[2*i+1]=(char)0xa9;}a[24000]=0;
    memset(b,'q',40000);b[20000]=0;
    FILE *old=book(a,0,0,1),*fresh=book(b,1,1,0),*out=fopen(path,"wb");assert(out);
    size_t bytes=(size_t)ftell(old);assert(!compound(out,old,bytes));assert(!fclose(out));fclose(old);
    /* Install an unrelated OLE stream in a spare directory slot/sector. */
    out=fopen(path,"r+b");assert(out);assert(!fseek(out,0,SEEK_END));long length=ftell(out);
    unsigned char *raw=malloc((size_t)length+512);assert(raw);rewind(out);assert(fread(raw,1,length,out)==(size_t)length);
    unsigned sector=(unsigned)length/512-1, fatsector=xe32(raw+76),dirsector=xe32(raw+48);
    directory_entry(raw+512*(1+dirsector)+256,"Object",2,UINT32_MAX,sector,512);
    put32(raw+512*(1+dirsector)+128+72,2);put32(raw+512*(1+fatsector)+sector*4,0xfffffffe);
    memset(raw+length,0x5a,512);rewind(out);assert(fwrite(raw,1,length+512,out)==(size_t)length+512);fclose(out);
    FILE *result=NULL;xlsWorkBook *wb=NULL;assert(!xe_edit(path,"Data",1,0,1,1,1,1,fresh,&result,&wb));
    out=fopen(edited,"wb");assert(out);assert(!xe_compound(out,path,wb,result,(size_t)ftell(result)));fclose(out);
    fclose(result);fclose(fresh);xls_close_WB(wb);
    xls_error_t error;wb=xls_open_file(edited,"UTF-8",&error);assert(wb);
    xlsWorkSheet *ws=xls_getWorkSheet(wb,0);assert(ws && !xls_parseWorkSheet(ws));
    assert(!strcmp(xls_cell(ws,0,0)->str,a));assert(!strcmp(xls_cell(ws,1,1)->str,b));
    assert(xls_cell(ws,8,5)->d==99);xls_close_WS(ws);xls_close_WB(wb);
    out=fopen(edited,"rb");assert(out);assert(!fseek(out,length,SEEK_SET));
    unsigned char preserved[512];assert(fread(preserved,1,512,out)==512);assert(!memcmp(preserved,raw+length,512));fclose(out);
    for(size_t n=0;n<20;n++) {size_t pos=0;xe_record r;unsigned char bad[20]={0};put16(bad+2,8225);assert(!xe_next(bad,n,&pos,&r));}
    free(raw);free(a);free(b);return 0;
}
'''

def main():
    vendor = ROOT / 'src/io/vendor'
    with tempfile.TemporaryDirectory(prefix='ctools-xls-edit-') as tmp:
        directory = Path(tmp)
        source = directory / 'test.c'; source.write_text(SOURCE)
        binary = directory / 'test'
        flags = ['-std=c11', '-O1', '-g', '-DSD_FASTMODE', '-fno-fast-math',
                 '-ffp-contract=off', '-fsanitize=address,undefined', '-fno-omit-frame-pointer',
                 '-I', str(ROOT/'src'), '-I', str(vendor/'libxls/include'),
                 '-I', str(vendor/'libxls/include/libxls'), '-D_DARWIN_C_SOURCE', '-D_GNU_SOURCE']
        if platform.system() == 'Darwin':
            flags += ['-DSYSTEM=APPLEMAC', '-Wl,-dead_strip', '-Wl,-undefined,dynamic_lookup', '-liconv']
        else:
            flags += ['-DSYSTEM=STUNIX', '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections', '-lm']
        subprocess.run([os.environ.get('CC','cc'), *flags, str(source),
                        *map(str,sorted((vendor/'libxls/src').glob('*.c'))), '-o', str(binary)],check=True)
        subprocess.run([str(binary),str(directory)],check=True,timeout=60,
                       env={**os.environ,'UBSAN_OPTIONS':'halt_on_error=1'})
        print('PASS XLS long SST continuations, formulas, preserved OLE streams, record bounds (ASan/UBSan)')
if __name__ == '__main__':
    main()

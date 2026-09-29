"""ASan/UBSan coverage for binary codecs, corrupt files, and cache cleanup."""
from pathlib import Path
import os
import platform
import subprocess
import tempfile
from run_io_formats import ROOT, shapefiles, codec_fixtures
from io_excel_format_fixtures import fixtures

SOURCE = r'''
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdint.h>
#include <string.h>
#include "stplugin.h"
#include "ctools_threads.h"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
static char filename[4096];
static int excel_strings, raw_storage;
static ST_boolean missing(double x){return x>=mock.missval;}
static ST_int poll(void){return 0;}
static ST_retcode macro(char *name,char *out,ST_int capacity){
    if(!strcmp(name,"___cio_filename"))snprintf(out,(size_t)capacity+1,"%s",filename);
    else if(!strcmp(name,"___cio_firstrow"))strcpy(out,"1");
    else if(!strcmp(name,"___cio_allstring") && excel_strings)strcpy(out,"1");
    else if(!strcmp(name,"___cio_rawstorage") && raw_storage)strcpy(out,"1");
    else out[0]=0;
    return 0;
}
static ST_retcode save(char *name,char *value){(void)name;(void)value;return 0;}
void ctools_error(const char *command,const char *format,...){(void)command;(void)format;}
#include "io/cio.c"
double ctools_timer_seconds(void){return 0;}
static ST_retcode display(char *s){(void)s;return 0;}
static int store_error;
static ST_int vars(void){return cache.nvar;}
static ST_int observations(void){return (ST_int)cache.nobs;}
static ST_boolean isstring(ST_int v){return cache.cols[v-1].is_string;}
static ST_boolean islong(ST_int v){return cache.cols[v-1].binary || cache.cols[v-1].width>2045;}
static ST_retcode checked_store(ST_int v,ST_int i,double value){
    (void)value;assert(v>=1 && v<=cache.nvar && i>=1 && (size_t)i<=cache.nobs);
    return store_error?459:0;
}
static void valid(const char *directory,const char *file,const char *kind,int columns,size_t rows){
    snprintf(filename,sizeof(filename),"%s/%s",directory,file);
    assert(scan(kind)==0);assert(cache.nvar==columns && cache.nobs==rows);
    if(!strcmp(kind,"spss") || !strcmp(kind,"sas")){
        assert(strlen(cache.cols[2].strings[0])==2100);
        assert(cache.cols[2].width==2100);
        char type[40];assert(!strcmp(storage(&cache.cols[2],type,sizeof(type)),"strL"));
        assert(!strcmp(cache.cols[2].format,"%9s"));
        assert(blob_info(3,1,0)==0);assert(blob_info(3,1,2101)==198);
        assert(load()==0);store_error=1;assert(load()==459);store_error=0;
    }
    cio_cleanup();cio_cleanup();assert(cache.nvar==0 && cache.cols==NULL);
}
int main(int argc,char **argv){
    assert(argc==2);uint64_t bits=UINT64_C(0x7fe0000000000000);memcpy(&mock.missval,&bits,8);
    mock.ismissing=missing;mock.spoutsml=display;mock.spouterr=display;mock.pollstd=poll;mock.macuse=macro;mock.macresave=save;
    mock.nvars=vars;mock.nobs=observations;mock.isstr=isstring;mock.isstrl=islong;
    mock.safestore=checked_store;/* Raw store deliberately NULL. */
    valid(argv[1],"fixture.zsav","spss",3,3);
    valid(argv[1],"long.sas7bdat","sas",3,3);
    valid(argv[1],"shape1.shp","shp",4,1);
    valid(argv[1],"multipart.shp","shp",5,5);
    snprintf(filename,sizeof(filename),"%s/raw.v8xpt",argv[1]);
    assert(!scan("sasxport8"));assert(cache.nvar==2 && cache.nobs==3);
    char type[40];assert(!strcmp(storage(&cache.cols[0],type,sizeof(type)),"byte"));
    assert(cache.cols[1].width==3);cio_cleanup();raw_storage=1;
    assert(!scan("sasxport8"));assert(cache.cols[1].width==12);
    assert(!strcmp(storage(&cache.cols[0],type,sizeof(type)),"double"));
    assert(!strcmp(cache.cols[1].format,"%12s"));cio_cleanup();raw_storage=0;
    for(int i=0;i<2;i++) {
        snprintf(filename,sizeof(filename),"%s/display_formats.%s",argv[1],i ? "xlsx" : "xls");
        assert(!scan(i ? "xlsx" : "xls"));assert(cache.nvar==36 && cache.nobs==4);
        assert(!strcmp(cache.cols[27].format,"%10.0g"));
        assert(cache.cols[27].numbers[0]==44300.5123456);
        assert(cache.cols[35].numbers[0]==44300.5123456);
        assert(!strcmp(cache.cols[4].format,"%tchh:MM_AM"));
        cio_cleanup();excel_strings=1;
        assert(!scan(i ? "xlsx" : "xls"));
        assert(!strcmp(cache.cols[27].strings[0],"44300.5123456"));
        assert(!strcmp(cache.cols[35].strings[0],"44300.5123456"));
        assert(!strcmp(cache.cols[4].strings[2]," 3:00 AM"));
        cio_cleanup();excel_strings=0;
    }
    ctools_destroy_global_pool();
    const char *files[]={"corrupt.sav","corrupt.dbf","corrupt.shp","corrupt.xls"};
    const char *kinds[]={"spss","dbase","shp","xls"};
    for(int j=0;j<4;j++){
        snprintf(filename,sizeof(filename),"%s/%s",argv[1],files[j]);
        assert(scan(kinds[j])!=0);cio_cleanup();
    }
    /* Truncated SAS pages and ZSAV blocks must not leave partially allocated
       columns or dangling codec state across repeated import attempts. */
    for(int file=0;file<2;file++){
        char original[4096];snprintf(original,sizeof(original),"%s/%s",argv[1],file?"fixture.zsav":"long.sas7bdat");
        FILE *input=fopen(original,"rb");assert(input);unsigned char bytes[16384];size_t n=fread(bytes,1,sizeof(bytes),input);fclose(input);
        for(size_t length=0;length<n;length+=97){
            snprintf(filename,sizeof(filename),"%s/truncated.bin",argv[1]);FILE *out=fopen(filename,"wb");assert(out);assert(fwrite(bytes,1,length,out)==length);assert(!fclose(out));
            (void)scan(file?"spss":"sas");cio_cleanup();
        }
    }
    return 0;
}
'''


def main():
    vendor = ROOT / 'src/io/vendor'
    with tempfile.TemporaryDirectory(prefix='ctools-binary-asan-') as tmp:
        directory = Path(tmp)
        shapefiles(directory)
        codec_fixtures(directory)
        fixtures(directory)
        cfile = directory / 'test.c'
        cfile.write_text(SOURCE)
        sources = [*sorted((vendor / 'readstat').rglob('*.c')),
                   *sorted((vendor / 'zlib').glob('*.c')),
                   *sorted((vendor / 'libxls/src').glob('*.c')),
                   *sorted((ROOT/'src/cimport/miniz').glob('*.c')),
                   *sorted((ROOT/'src/cimport/libdeflate').rglob('*.c')),
                   ROOT / 'src/cexport/cexport_parse.c',
                   *[ROOT/'src'/name for name in ['cimport/cimport_xlsx.c','cimport/cimport_xlsx_xml.c','cimport/cimport_xlsx_zip.c','ctools_types.c','ctools_arena.c','ctools_threads.c']]]
        flags = ['-std=c11', '-O1', '-g', '-DHAVE_ZLIB=1', '-DSD_FASTMODE',
                 '-D_GNU_SOURCE', '-D_DARWIN_C_SOURCE', '-pthread', '-fno-fast-math', '-ffp-contract=off',
                 '-fsanitize=address,undefined', '-fno-omit-frame-pointer', '-I', str(ROOT / 'src'),
                 '-I', str(vendor / 'readstat'), '-I', str(vendor / 'zlib'),
                 '-I', str(vendor / 'libxls/include'), '-I', str(vendor / 'libxls/include/libxls')]
        if platform.system() == 'Darwin':
            flags += ['-DSYSTEM=APPLEMAC', '-Wl,-dead_strip', '-liconv']
        else:
            flags += ['-DSYSTEM=STUNIX', '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections', '-lm']
        binary = directory / 'test'
        subprocess.run([os.environ.get('CC', 'cc'), *flags, str(cfile), *map(str, sources), '-o', str(binary)], check=True)
        subprocess.run([str(binary), str(directory)], check=True, timeout=60,
                       env={**os.environ, 'UBSAN_OPTIONS': 'halt_on_error=1'})
        print('PASS binary codecs, truncation, and cache cleanup (ASan/UBSan)')


if __name__ == '__main__':
    main()

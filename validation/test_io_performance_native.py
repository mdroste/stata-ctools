"""Boundary and equivalence checks for the I/O performance paths (ASan/UBSan)."""
from native.support.compiler import compile_native
from pathlib import Path
import os
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]

CRC = r'''
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include "vendor/libdeflate/libdeflate.h"
#include "vendor/miniz/miniz.h"
int main(void) {
    unsigned char *data = malloc(100000 + 32);
    assert(data);
    uint32_t state = 95183;
    for (int i=0; i<100032; i++) {
        state ^= state << 13; state ^= state >> 17; state ^= state << 5;
        data[i] = (unsigned char)state;
    }
    assert(libdeflate_crc32(0,"123456789",9)==UINT32_C(0xcbf43926));
    for (int offset=0; offset<32; offset++) {
        for (size_t n=0; n<100000; n=n<128?n+1:n*2+1) {
            uint32_t expected=(uint32_t)mz_crc32(0,data+offset,n);
            assert(libdeflate_crc32(0,data+offset,n)==expected);
            size_t split=n/3;
            uint32_t part=libdeflate_crc32(0,data+offset,split);
            part=libdeflate_crc32(part,data+offset+split,n-split);
            assert(part==expected);
        }
    }
    free(data);
}
'''

XML = r'''
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
#include "cimport/cimport_xlsx.c"

static void check(const char *attrs, const char *body, int type,
                  double value, const char *text) {
    char xml[2048];
    int n=snprintf(xml,sizeof(xml),"<c r=\"A1\" %s>%s</c>",attrs,body);
    assert(n>0 && n<(int)sizeof(xml));
    /* Every prefix is separately allocated: ASan checks the bounded lookahead. */
    for (int length=0;length<=n;length++) {
        char *input=malloc(length ? (size_t)length : 1);
        assert(input);memcpy(input,xml,(size_t)length);
        XLSXContext x; xlsx_context_init(&x);
        assert(xlsx_cm_reserve(&x,1,1));
        xlsx_scan_local loc; xlsx_scan_local_init(&loc);
        loc.arena=&x.parse_arena;loc.grow=&x;
        const char *end=input+length;
        const char *gt=length ? xlsx_markup_end(input,end) : NULL;
        int prev=-1;
        const char *after=gt ? xlsx_scan_cell(&x,&loc,input+2,gt,end,1,&prev,false) : NULL;
        if (length==n) {
            assert(after==end && !loc.failed);
            assert(x.cm_types[0][0]==type);
            if (text) {
                const char *parsed=NULL;memcpy(&parsed,&x.cm_numeric[0][0],sizeof(parsed));
                assert(parsed && !strcmp(parsed,text));
            } else if (type!=XLSX_CELL_EMPTY) assert(x.cm_numeric[0][0]==value);
        }
        xlsx_context_free(&x);free(input);
    }
}
int main(void) {
    uint64_t bits=UINT64_C(0x7fe0000000000000);memcpy(&mock.missval,&bits,8);
    /* Each compact form also has a general-parser equivalent. */
    const char *values[]={"1.25", "-0.00025", "1e20", "1,5", " 2 ", "NaN", ""};
    double numbers[]={1.25,-0.00025,1e20,1,2,0,0};
    for (int i=0;i<7;i++) {
        char compact[128], general[128];
        snprintf(compact,sizeof(compact),"<v>%s</v>",values[i]);
        snprintf(general,sizeof(general),"<v>%s</v> ",values[i]);
        int type=i<5?XLSX_CELL_NUMBER:XLSX_CELL_EMPTY;
        check("",compact,type,numbers[i],NULL);check("",general,type,numbers[i],NULL);
    }
    check("t=\"inlineStr\"","<is><t>A&amp;B&lt;C</t></is>",XLSX_CELL_STRING,0,"A&B<C");
    check("t=\"inlineStr\"","<is><t xml:space=\"preserve\"> A </t></is>",XLSX_CELL_STRING,0," A ");
    check("t=\"inlineStr\"","<is><r><t>first</t></r><r><t>second</t></r></is>",XLSX_CELL_STRING,0,"first");
    check("t=\"s\"","<v>3</v>",XLSX_CELL_SHARED_STRING,3,NULL);
    check("t=\"b\"","<v>1</v>",XLSX_CELL_BOOLEAN,1,NULL);
    check("t=\"e\"","<v>#N/A</v>",XLSX_CELL_EMPTY,0,NULL);
    check("t=\"str\"","<v>A&amp;B</v>",XLSX_CELL_STRING,0,"A&B");
    check("","<f>SUM(A2:A3)</f><v>12</v>",XLSX_CELL_NUMBER,12,NULL);
    check("","<x:v>42</x:v>",XLSX_CELL_NUMBER,42,NULL);
    check("note='a>b'","<v>7</v><v>8</v>",XLSX_CELL_NUMBER,7,NULL);
    check("t=\"inlineStr\"","<is><t/></is>",XLSX_CELL_EMPTY,0,NULL);
}
'''


CSV = r'''
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
static char explicit_encoding[40];
static ST_retcode macro(char *key,char *out,ST_int cap) {
    snprintf(out,(size_t)cap+1,"%s",!strcmp(key,"___cimport_encoding_name")?explicit_encoding:"");
    return 0;
}
static ST_retcode display(char *s){(void)s;return 0;}
static ST_int poll(void){return 0;}
double ctools_timer_ms(void){return 0;}
#include "cimport/cimport_impl.c"
int main(void) {
    mock.macuse=macro;mock.spoutsml=display;mock.spouterr=display;mock.pollstd=poll;
    uint64_t bits=UINT64_C(0x7fe0000000000000);memcpy(&mock.missval,&bits,8);
    ctools_set_max_threads(1);
    char data[32000]="v1,v2\n";
    for(int i=0;i<1000;i++)strcat(data,"123.456789,-678.901234\n");
    size_t n=strlen(data);
    const char *charset=cimport_detect_charset(data,n);
    assert(!strcmp(charset,"ISO-8859-2"));
    for(int nonascii=0;nonascii<2;nonascii++) {
        if(nonascii)strcpy(data+n,"1,\xa3\n");
        size_t size=strlen(data);
        char path[]="/tmp/ctools-ascii-import-XXXXXX";
        int fd=mkstemp(path);assert(fd>=0);
        assert(write(fd,data,size)==(ssize_t)size);close(fd);
        CImportContext *ctx=cimport_parse_csv(path,',',true,0,false,CIMPORT_BINDQUOTES_LOOSE,
            CIMPORT_NUMTYPE_AUTO,'.',0,CIMPORT_EMPTYLINES_SKIP,20,CIMPORT_ENC_UNKNOWN,NULL,0,NULL,0);
        assert(ctx && !strcmp(ctx->reported_encoding,charset));
        char *expected=NULL;size_t expected_size;
        assert(!cimport_convert_named_to_utf8(data,size,charset,&expected,&expected_size));
        assert(ctx->file_size==expected_size && !memcmp(ctx->file_data,expected,expected_size));
        assert(ctx->total_rows==(size_t)(1000+nonascii));
        if(!nonascii)assert(ctx->converted_data==NULL);
        else assert(ctx->converted_data!=NULL);
        free(expected);cimport_free_context(ctx);
        strcpy(explicit_encoding,"bogus");
        ctx=cimport_parse_csv(path,',',true,0,false,CIMPORT_BINDQUOTES_LOOSE,
            CIMPORT_NUMTYPE_AUTO,'.',0,CIMPORT_EMPTYLINES_SKIP,20,CIMPORT_ENC_UNKNOWN,NULL,0,NULL,0);
        assert(!ctx && g_cimport_parse_status==198);explicit_encoding[0]=0;
        unlink(path);
    }
    ctools_destroy_global_pool();
}
'''


def main():
    mac = platform.system() == 'Darwin'
    vendor = ROOT / 'src/vendor'
    with tempfile.TemporaryDirectory(prefix='ctools-io-performance-') as tmp:
        directory = Path(tmp)
        flags = ['-std=c11', '-O1', '-g', '-DSD_FASTMODE', '-fno-fast-math',
                 '-ffp-contract=off', '-fno-strict-aliasing',
                 '-fsanitize=address,undefined', '-fno-omit-frame-pointer',
                 '-pthread', '-I', str(ROOT/'src')]
        flags += (['-DSYSTEM=APPLEMAC', '-D_DARWIN_C_SOURCE', '-Wl,-dead_strip',
                   '-Wl,-undefined,dynamic_lookup'] if mac else
                  ['-DSYSTEM=STUNIX', '-D_GNU_SOURCE', '-ffunction-sections',
                   '-fdata-sections', '-Wl,--gc-sections'])
        crc_sources = [vendor/'libdeflate/lib/crc32.c', vendor/'miniz/miniz.c',
                       vendor/'libdeflate/lib/arm/cpu_features.c',
                       vendor/'libdeflate/lib/x86/cpu_features.c']
        xml_sources = [ROOT/'src'/p for p in ('ctools_types.c', 'ctools_arena.c',
                       'ctools_threads.c', 'formats/xlsx_xml.c')]
        csv_sources = [ROOT/'src'/p for p in ('ctools_types.c', 'ctools_arena.c',
                       'ctools_threads.c', 'cimport/cimport_parse.c',
                       'cimport/cimport_mmap.c', 'cimport/cimport_encoding.c')]
        for name, source, modules in [('crc', CRC, crc_sources), ('xml', XML, xml_sources),
                                      ('csv', CSV, csv_sources)]:
            cfile=directory/(name+'.c'); cfile.write_text(source)
            binary=directory/name
            compile_native([os.environ.get('CC','cc'), *flags, str(cfile), *map(str,modules),
                            *(['-liconv'] if mac else ['-lm']), '-o', str(binary)],check=True)
            subprocess.run([str(binary)],check=True,
                           env={**os.environ,'UBSAN_OPTIONS':'halt_on_error=1'})
            print('PASS I/O performance '+name+' (ASan/UBSan)',flush=True)


if __name__ == '__main__':
    main()
